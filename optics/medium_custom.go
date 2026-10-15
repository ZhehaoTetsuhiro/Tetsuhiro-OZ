package optics

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"strings"
)

// ---------------------------------------------------------------------------
// User-defined inhomogeneous media (自定义不均匀介质, element type
// "custom_medium")
//
// The user gives the refractive index as an *expression* over the spatial
// coordinates and the wavelength, n(x, y, z, λ), optionally with an
// absorption/gain coefficient k(x, y, z, λ). The medium is a z-stack of
// 1..N layers, each with its own expression and its own thickness; every
// z-interval covered by a layer is that layer's medium, and the light travels
// through the stack in order. Outside the layers nothing is applied (vacuum).
//
// What the simulation does with it (the physics, in three lines):
//
//	per layer, per slice of thickness dz centered at depth zc:
//	    E <- ASM_{nref}(dz/2) [ exp(i k0 (n(x,y,zc) - nref) dz)
//	                            exp(-k0 k(x,y,zc) dz) ] ASM_{nref}(dz/2)
//
//	nref  the layer's reference index (below), used by the diffraction kernel.
//	k0    2*pi/lambda. Im(n) > 0 absorbs, < 0 amplifies.
//
// The reference index is what makes the difference between "a slab is a phase
// plate" (reference 1: the classic vacuum-reference BPM, no refraction at all)
// and "a slab is a slab". By default each layer takes its reference from its
// own index expression evaluated on the optical axis at the layer's mid-plane
// — for a homogeneous layer that value *is* the layer index, the screen term
// vanishes identically and the layer is propagated by the exact homogeneous
// transfer function of index n (Snell refraction at the faces included, axial
// phase k0 n d exact, focal shift d(1-1/n)). A layer whose index varies in
// space may set "reference" explicitly to pick the expansion point; the screen
// then carries the full deviation n - nref. See docs/PHYSICS.md §5.3 and
// medium.go's reference-index block comment.
//
// The expressions are parsed and compiled *once*, at element construction time
// (see medium_expr.go); the propagation loop evaluates closures over a reused
// environment struct and allocates nothing per point. Layers whose index and
// absorption do not vary in space skip the per-pixel evaluation entirely.
// ---------------------------------------------------------------------------

// Limits of a custom medium. They protect memory and CPU the same way the
// element/plane limits of the simulator do.
const (
	// MaxCustomMediumLayers bounds the number of z layers in one medium.
	MaxCustomMediumLayers = 64
	// MaxCustomMediumDistance bounds the total thickness of one medium (m).
	MaxCustomMediumDistance = 1.0
	// MaxCustomMediumSteps bounds the total number of split-step slices of one
	// medium (per-layer overrides included); each slice costs two FFTs.
	MaxCustomMediumSteps = 4096
	// DefaultCustomMediumSteps is the total slice budget when "steps" is not
	// given, matching the homogeneous "medium" element's default.
	DefaultCustomMediumSteps = 20
	// customMediumProbeWavelength is the wavelength at which the expressions
	// are probed (finiteness, n > 0, reference resolvability) when the element
	// is built, i.e. before a run wavelength exists.
	customMediumProbeWavelength = 632.8e-9
)

// CustomMediumLayerSpec is one z-interval of a custom medium.
//
//	{"thickness": 1e-3,                         // metres, > 0
//	 "index": "1.5 + 1e-3*z/1e-3",              // n(x,y,z,zl,λ), real part
//	 "absorption": "0.05*step(r-1e-3)",          // optional Im(n), >0 absorbs
//	 "reference": 1.5,                           // optional, default: n on axis
//	 "params": {"g": 100.0},                     // named constants for the exprs
//	 "steps": 40}                                // optional slice override
type CustomMediumLayerSpec struct {
	Thickness  float64            `json:"thickness"`
	Index      string             `json:"index"`
	Absorption string             `json:"absorption,omitempty"`
	Reference  *float64           `json:"reference,omitempty"`
	Params     map[string]float64 `json:"params,omitempty"`
	Steps      int                `json:"steps,omitempty"`
}

// CustomMediumSpec is the whole medium: either a z-stack of layers, or the flat
// single-layer form (distance + index + absorption) that mirrors the built-in
// homogeneous "medium" element — give one form or the other, not both.
type CustomMediumSpec struct {
	Layers []CustomMediumLayerSpec `json:"layers,omitempty"`
	// Flat single-layer form.
	Distance   float64            `json:"distance,omitempty"`
	Index      string             `json:"index,omitempty"`
	Absorption string             `json:"absorption,omitempty"`
	Params     map[string]float64 `json:"params,omitempty"`
	// Steps is the total split-step budget of the medium, shared by the layers
	// in proportion to their thickness (each layer gets at least one slice).
	Steps int `json:"steps,omitempty"`
}

// customLayer is one compiled layer.
type customLayer struct {
	spec      CustomMediumLayerSpec
	z0        float64 // depth of the layer's entrance face inside the medium
	thickness float64
	steps     int // per-layer override (0 = proportional share of the budget)
	n         *MediumExpr
	abs       *MediumExpr // nil when the layer is lossless
	refFixed  bool
	refValue  float64
	uniform   bool // index and absorption do not vary in space
}

// CustomMedium is a compiled medium: parsed expressions, resolved layer
// geometry and the slice plan, ready to propagate.
type CustomMedium struct {
	layers   []customLayer
	distance float64
	steps    int
}

// MediumLayerInfo describes one compiled layer for inspection, the catalog API
// and tests (no expressions are needed to interpret it).
type MediumLayerInfo struct {
	Index      int        `json:"index"`
	Depth      float64    `json:"depth"`      // m, front face inside the medium
	Thickness  float64    `json:"thickness"`  // m
	Steps      int        `json:"steps"`      // slices used by this layer
	IndexExpr  string     `json:"index_expr"` //
	AbsExpr    string     `json:"absorption_expr,omitempty"`
	Uniform    bool       `json:"uniform"`   // value independent of x,y,z
	Reference  float64    `json:"reference"` // at ResolveInfo's wavelength
	AutoRef    bool       `json:"reference_auto"`
	Names      []string   `json:"names,omitempty"` // variables/params referenced
	LayerEdges [2]float64 `json:"layer_edges,omitempty"`
}

// NewCustomMedium compiles a spec into a medium. Every expression is parsed and
// type/name-checked here, the layer geometry is validated, and the expressions
// are probed on a sample grid at the probe wavelength so a NaN, a non-positive
// index or an unresolvable reference index is reported at load time instead of
// spreading over the field.
func NewCustomMedium(spec CustomMediumSpec) (*CustomMedium, error) {
	layers := spec.Layers
	if len(layers) == 0 {
		if strings.TrimSpace(spec.Index) == "" {
			return nil, fmt.Errorf("custom_medium: no index expression (give \"layers\":[...] or the flat \"distance\"+\"index\" form)")
		}
		if !(spec.Distance > 0) || math.IsNaN(spec.Distance) {
			return nil, fmt.Errorf("custom_medium: thickness \"distance\" must be > 0 (got %g m)", spec.Distance)
		}
		layers = []CustomMediumLayerSpec{{
			Thickness: spec.Distance, Index: spec.Index, Absorption: spec.Absorption,
			Params: spec.Params, Steps: spec.Steps,
		}}
	} else if strings.TrimSpace(spec.Index) != "" {
		return nil, fmt.Errorf("custom_medium: give either \"layers\" or the flat \"distance\"+\"index\" form, not both")
	}
	if len(layers) > MaxCustomMediumLayers {
		return nil, fmt.Errorf("custom_medium: too many layers (%d, limit %d)", len(layers), MaxCustomMediumLayers)
	}
	m := &CustomMedium{steps: spec.Steps}
	z0 := 0.0
	for i := range layers {
		L, err := compileCustomLayer(layers[i], i, z0)
		if err != nil {
			return nil, err
		}
		m.layers = append(m.layers, L)
		z0 += L.thickness
	}
	m.distance = z0
	if m.distance > MaxCustomMediumDistance {
		return nil, fmt.Errorf("custom_medium: total thickness %g m exceeds the limit %g m", m.distance, MaxCustomMediumDistance)
	}
	if m.steps < 0 {
		return nil, fmt.Errorf("custom_medium: steps must be >= 0 (got %d)", m.steps)
	}
	if m.steps == 0 {
		m.steps = DefaultCustomMediumSteps
	}
	total := 0
	for i := range m.layers {
		total += m.stepsFor(i)
	}
	if total > MaxCustomMediumSteps {
		return nil, fmt.Errorf("custom_medium: %d split-step slices in total exceeds the limit %d — lower \"steps\" or a layer's \"steps\" (each slice costs two FFTs)", total, MaxCustomMediumSteps)
	}
	return m, nil
}

// compileCustomLayer validates one layer spec and compiles its expressions.
func compileCustomLayer(spec CustomMediumLayerSpec, i int, z0 float64) (customLayer, error) {
	where := fmt.Sprintf("custom_medium layer %d", i)
	L := customLayer{spec: spec, z0: z0, thickness: spec.Thickness, steps: spec.Steps}
	if !(spec.Thickness > 0) || math.IsNaN(spec.Thickness) || math.IsInf(spec.Thickness, 0) {
		return L, fmt.Errorf("%s: thickness must be > 0 (got %g m)", where, spec.Thickness)
	}
	if strings.TrimSpace(spec.Index) == "" {
		return L, fmt.Errorf("%s: index expression is required (n(x,y,z,zl,λ), real part)", where)
	}
	for k, v := range spec.Params {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return L, fmt.Errorf("%s: parameter %q must be finite (got %g)", where, k, v)
		}
	}
	n, err := CompileMediumExpr(spec.Index, spec.Params)
	if err != nil {
		return L, fmt.Errorf("%s index: %v", where, err)
	}
	L.n = n
	if s := strings.TrimSpace(spec.Absorption); s != "" {
		a, err := CompileMediumExpr(s, spec.Params)
		if err != nil {
			return L, fmt.Errorf("%s absorption: %v", where, err)
		}
		L.abs = a
	}
	L.uniform = !n.IsSpatial() && (L.abs == nil || !L.abs.IsSpatial())
	if spec.Reference != nil {
		if !(*spec.Reference > 0) || math.IsNaN(*spec.Reference) || math.IsInf(*spec.Reference, 0) {
			return L, fmt.Errorf("%s: reference index must be a finite number > 0 (got %g)", where, *spec.Reference)
		}
		L.refFixed, L.refValue = true, *spec.Reference
	}
	if L.steps < 0 {
		return L, fmt.Errorf("%s: steps must be >= 0 (got %d)", where, L.steps)
	}
	if L.steps > MaxCustomMediumSteps {
		return L, fmt.Errorf("%s: steps %d exceeds the limit %d", where, L.steps, MaxCustomMediumSteps)
	}
	if err := probeCustomLayer(&L, i, customMediumProbeWavelength); err != nil {
		return L, err
	}
	return L, nil
}

// probeCustomLayer evaluates a layer's expressions on a small sample grid at
// the given wavelength: every value must be finite and the index positive, and
// the reference index must be derivable. The probe catches a NaN or a negative
// index at load time — after one propagation step a NaN has already spread over
// the whole plane — and it is repeated with the real wavelength at run time
// (see PropagateCustomMedium), because a dispersion term can move the value
// with λ.
func probeCustomLayer(L *customLayer, li int, wl float64) error {
	samples := []float64{-4e-3, -1e-3, 0, 1e-3, 4e-3}
	for _, z := range []float64{0, L.thickness / 2, L.thickness} {
		for _, x := range samples {
			for _, y := range samples {
				if err := checkMediumIndex(L, li, x, y, z, z, wl); err != nil {
					return err
				}
			}
		}
	}
	if !L.refFixed {
		if _, err := L.autoReference(wl); err != nil {
			return err
		}
	}
	return nil
}

// checkMediumIndex validates one sample point of a layer.
func checkMediumIndex(L *customLayer, li int, x, y, z, zl, wl float64) error {
	where := layerName(li, L)
	v := L.n.Eval(x, y, z, zl, wl)
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return fmt.Errorf("%s: index expression is not finite (%g) at x=%g, y=%g, z=%g m — guard the singularity with if(...) or clamp(...)",
			where, v, x, y, z)
	}
	if v <= 0 {
		return fmt.Errorf("%s: index must be > 0 (got %g at x=%g, y=%g, z=%g m)", where, v, x, y, z)
	}
	if L.abs != nil {
		a := L.abs.Eval(x, y, z, zl, wl)
		if math.IsNaN(a) || math.IsInf(a, 0) {
			return fmt.Errorf("%s: absorption expression is not finite (%g) at x=%g, y=%g, z=%g m", where, a, x, y, z)
		}
	}
	return nil
}

// layerName renders a layer's ordinal and depth span for messages.
func layerName(li int, L *customLayer) string {
	return fmt.Sprintf("custom_medium layer %d (z = %g…%g m)", li, L.z0, L.z0+L.thickness)
}

// autoReference resolves the default reference index of a layer: its own index
// expression on the optical axis (x = y = 0) at the layer's mid-plane. For a
// homogeneous layer this is exactly the layer index, so the screen term of that
// layer vanishes and the propagation becomes the exact homogeneous-medium
// transfer function.
func (L *customLayer) autoReference(wl float64) (float64, error) {
	zm := L.thickness / 2
	v := L.n.Eval(0, 0, zm, zm, wl)
	if math.IsNaN(v) || math.IsInf(v, 0) || v <= 0 {
		return 0, fmt.Errorf("custom_medium layer at z = %g…%g m: cannot derive a reference index from the index expression on the optical axis (n(0,0,%g) = %g at λ = %g m); set \"reference\" explicitly",
			L.z0, L.z0+L.thickness, zm, v, wl)
	}
	return v, nil
}

// Distance returns the total thickness of the medium (m).
func (m *CustomMedium) Distance() float64 { return m.distance }

// SliceCount returns the total number of split-step slices the medium will run.
func (m *CustomMedium) SliceCount() int {
	total := 0
	for i := range m.layers {
		total += m.stepsFor(i)
	}
	return total
}

// frameSlots is the argument-arena size the layer's expressions need (the index
// and absorption expressions are evaluated one at a time).
func (L *customLayer) frameSlots() int {
	n := L.n.FrameSlots()
	if L.abs != nil && L.abs.FrameSlots() > n {
		n = L.abs.FrameSlots()
	}
	return n
}

// stepsFor returns the slice count of layer i: its own override, or its share
// of the medium's budget in proportion to its thickness (at least one slice).
func (m *CustomMedium) stepsFor(i int) int {
	L := &m.layers[i]
	if L.steps > 0 {
		return L.steps
	}
	if m.steps <= 0 {
		return 1
	}
	share := int(math.Round(float64(m.steps) * L.thickness / m.distance))
	if share < 1 {
		share = 1
	}
	if share > MaxCustomMediumSteps {
		share = MaxCustomMediumSteps
	}
	return share
}

// referenceIndices resolves every layer's reference index at the given
// wavelength (fixed values as given, automatic ones from the layer expression).
func (m *CustomMedium) referenceIndices(wl float64) ([]float64, error) {
	out := make([]float64, len(m.layers))
	for i := range m.layers {
		L := &m.layers[i]
		if L.refFixed {
			out[i] = L.refValue
			continue
		}
		v, err := L.autoReference(wl)
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// LayerInfo describes the compiled layers at the given wavelength: depths,
// thickness, slice counts, the reference index actually used and whether the
// layer is spatially uniform (and therefore propagated exactly).
func (m *CustomMedium) LayerInfo(wl float64) ([]MediumLayerInfo, error) {
	refs, err := m.referenceIndices(wl)
	if err != nil {
		return nil, err
	}
	out := make([]MediumLayerInfo, len(m.layers))
	for i := range m.layers {
		L := &m.layers[i]
		info := MediumLayerInfo{
			Index: i, Depth: L.z0, Thickness: L.thickness, Steps: m.stepsFor(i),
			IndexExpr: L.n.Source(), Uniform: L.uniform, Reference: refs[i],
			AutoRef: !L.refFixed, Names: L.n.Names(),
			LayerEdges: [2]float64{L.z0, L.z0 + L.thickness},
		}
		if L.abs != nil {
			info.AbsExpr = L.abs.Source()
		}
		out[i] = info
	}
	return out, nil
}

// PropagateCustomMedium advances f through the whole layer stack. It is the
// element's operator: each layer is propagated with the reference-index
// split-step (medium.go) at the wavelength of ctx.
func PropagateCustomMedium(f *Field, m *CustomMedium, ctx *Context) error {
	if m == nil || f == nil || len(m.layers) == 0 || m.distance <= 0 {
		return nil
	}
	if ctx == nil {
		ctx = &Context{Evanescent: "decay"}
	}
	wl := ctx.Wavelength
	if !(wl > 0) || math.IsNaN(wl) || math.IsInf(wl, 0) {
		return fmt.Errorf("custom_medium: wavelength must be > 0 (got %g m)", wl)
	}
	refs, err := m.referenceIndices(wl)
	if err != nil {
		return err
	}
	// One environment (with one argument arena sized for the widest layer
	// expression), one coordinate row and one phase row are reused across every
	// layer, slice and pixel: the propagation loop allocates nothing.
	frames := 0
	for i := range m.layers {
		if n := m.layers[i].frameSlots(); n > frames {
			frames = n
		}
	}
	env := &mediumEnv{wl: wl, k: 2 * math.Pi / wl}
	if frames > 0 {
		env.frames = make([]float64, frames)
	}
	xs := make([]float64, f.N)
	for i := range xs {
		xs[i] = f.X(i)
	}
	scratch := &mediumScratch{xs: xs}
	for i := range m.layers {
		if err := propagateCustomLayer(f, &m.layers[i], i, refs[i], env, scratch, m.stepsFor(i), ctx); err != nil {
			return err
		}
	}
	return nil
}

// mediumScratch holds the reusable buffers of one PropagateCustomMedium call.
type mediumScratch struct {
	xs      []float64 // transverse coordinates of one row
	prevRow []float64 // phase screen of the previous row (sampling check)
}

// propagateCustomLayer runs the symmetric split-step across one layer.
func propagateCustomLayer(f *Field, L *customLayer, li int, nref float64, env *mediumEnv, sc *mediumScratch, steps int, ctx *Context) error {
	if steps < 1 {
		steps = 1
	}
	dz := L.thickness / float64(steps)
	zero, reg, _, _ := ctx.evanescentPolicy(dz / 2)
	for s := 0; s < steps; s++ {
		zc := L.z0 + (float64(s)+0.5)*dz // depth from the medium's entrance
		zl := (float64(s) + 0.5) * dz    // depth from this layer's entrance
		pIn := f.Power()
		propagateMediumHalf(f, dz/2, nref, ctx)
		if err := applyCustomScreen(f, L, li, zc, zl, nref, dz, env, sc, ctx); err != nil {
			return err
		}
		propagateMediumHalf(f, dz/2, nref, ctx)
		asmWarn(f, dz/2, ctx, pIn, f.Power(), zero, reg)
	}
	return nil
}

// applyCustomScreen multiplies the field by the reference-relative phase screen
// exp(i k0 (n - nref) dz) * exp(-k0 Im(n) dz) of the slice centered at depth
// zc. Every value is evaluated once per pixel from the compiled closures.
//
// Two paths: a layer whose index and absorption do not vary in space gets a
// single constant factor per slice (a plain scalar multiply), which is both the
// common "glass plate" case and the fastest; a spatially varying layer is
// evaluated per pixel, and while doing so the transverse phase gradient of the
// screen is monitored: a screen faster than Pi per pixel aliases, and that is
// reported as a warning rather than silently producing a wrong field.
func applyCustomScreen(f *Field, L *customLayer, li int, zc, zl, nref, dz float64, env *mediumEnv, sc *mediumScratch, ctx *Context) error {
	k0 := env.k
	env.z, env.zl = zc, zl
	if L.uniform {
		env.x, env.y = 0, 0
		nv := L.n.evalEnv(env)
		if !mediumValueOK(nv) {
			return fmt.Errorf("custom_medium layer %d: index expression is not finite or not positive (%g) at z=%g m", li, nv, zc)
		}
		att := 0.0
		if L.abs != nil {
			att = L.abs.evalEnv(env)
			if math.IsNaN(att) || math.IsInf(att, 0) {
				return fmt.Errorf("custom_medium layer %d: absorption expression is not finite (%g) at z=%g m", li, att, zc)
			}
		}
		t := complex(math.Exp(-k0*att*dz), 0) * cexpI(k0*(nv-nref)*dz)
		if t == 1 {
			return nil
		}
		for i := range f.Ex {
			f.Ex[i] *= t
		}
		if f.Polarized {
			for i := range f.Ey {
				f.Ey[i] *= t
			}
		}
		return nil
	}

	nn := f.N
	xs := sc.xs
	if sc.prevRow == nil || len(sc.prevRow) != nn {
		sc.prevRow = make([]float64, nn)
	}
	prevRow := sc.prevRow
	maxStep := 0.0
	bad, badX, badY := 0, 0.0, 0.0
	for j := 0; j < nn; j++ {
		y := f.Y(j)
		env.y = y
		row := j * nn
		prevPh := 0.0
		for i := 0; i < nn; i++ {
			env.x = xs[i]
			nv := L.n.evalEnv(env)
			att := 0.0
			if L.abs != nil {
				att = L.abs.evalEnv(env)
			}
			if !mediumValueOK(nv) || math.IsNaN(att) || math.IsInf(att, 0) {
				if bad == 0 {
					badX, badY = xs[i], y
				}
				bad++
				f.Ex[row+i] = 0
				if f.Polarized {
					f.Ey[row+i] = 0
				}
				continue
			}
			ph := k0 * (nv - nref) * dz
			if i > 0 {
				if d := math.Abs(ph - prevPh); d > maxStep {
					maxStep = d
				}
			}
			if j > 0 {
				if d := math.Abs(ph - prevRow[i]); d > maxStep {
					maxStep = d
				}
			}
			prevPh = ph
			prevRow[i] = ph
			t := complex(math.Exp(-k0*att*dz), 0) * cexpI(ph)
			f.Ex[row+i] *= t
			if f.Polarized {
				f.Ey[row+i] *= t
			}
		}
	}
	if bad > 0 {
		return fmt.Errorf("custom_medium layer %d: index/absorption expression returned a non-finite or non-positive value at x=%g, y=%g, z=%g m (%d pixels of the slice) — guard it with if(...)/clamp(...) or fix the parameters",
			li, badX, badY, zc, bad)
	}
	if maxStep > math.Pi && ctx != nil {
		ctx.Warnings.Add("custom_medium_alias",
			fmt.Sprintf("自定义介质第 %d 层在 z≈%.4g m 处的相位屏横向梯度达到 %.2f rad/像素（>π）：该层采样不足，grid/steps 需要加密，否则结果混叠", li+1, zc, maxStep),
			maxStep/math.Pi)
	}
	return nil
}

// mediumValueOK reports whether a sampled index value is usable: finite and > 0.
func mediumValueOK(v float64) bool { return v > 0 && v < 1e9 }

// ---------------------------------------------------------------------------
// Element glue
// ---------------------------------------------------------------------------

// customMediumEl is the "custom_medium" element: a user-defined inhomogeneous
// medium (expression index, z layers, optional absorption/gain) as a simulator
// element, usable in element trains and as a positioned scene component.
type customMediumEl struct {
	m *CustomMedium
}

// newCustomMedium builds the element from its raw JSON params. Both the layered
// and the flat form are accepted, and "layers" may be given either as a JSON
// array or as a string holding that array (the GUI's text control).
func newCustomMedium(p map[string]any) (Element, error) {
	spec, err := decodeCustomMediumParams(p)
	if err != nil {
		return nil, err
	}
	m, err := NewCustomMedium(*spec)
	if err != nil {
		return nil, err
	}
	return &customMediumEl{m: m}, nil
}

// Apply propagates the field through the medium (the element's operator).
func (e *customMediumEl) Apply(f *Field, ctx *Context) error {
	return PropagateCustomMedium(f, e.m, ctx)
}

// decodeCustomMediumParams normalizes element params into a CustomMediumSpec.
// Unknown fields are rejected: a typo in a parameter name is a configuration
// error, not a silently ignored setting (the layers are therefore decoded
// strictly, see decodeStrictJSON).
func decodeCustomMediumParams(p map[string]any) (*CustomMediumSpec, error) {
	if p == nil {
		return nil, fmt.Errorf("custom_medium: no parameters (give \"layers\" or \"distance\"+\"index\")")
	}
	spec := &CustomMediumSpec{}
	rest := map[string]any{}
	for k, v := range p {
		switch k {
		case "layers":
			layers, err := decodeLayerList(v)
			if err != nil {
				return nil, fmt.Errorf("custom_medium: %v", err)
			}
			spec.Layers = layers
		case "distance", "thickness":
			f, err := asFloat(v)
			if err != nil {
				return nil, fmt.Errorf("custom_medium: %s must be a number", k)
			}
			spec.Distance = f
		case "index":
			s, ok := v.(string)
			if !ok {
				return nil, fmt.Errorf("custom_medium: index must be an expression string")
			}
			spec.Index = s
		case "absorption":
			s, ok := v.(string)
			if !ok {
				return nil, fmt.Errorf("custom_medium: absorption must be an expression string (Im(n), >0 absorbs)")
			}
			spec.Absorption = s
		case "steps":
			n, err := asFloat(v)
			if err != nil {
				return nil, fmt.Errorf("custom_medium: steps must be an integer")
			}
			spec.Steps = int(n)
		case "params":
			m, err := decodeParamMap(v)
			if err != nil {
				return nil, fmt.Errorf("custom_medium: params: %v", err)
			}
			spec.Params = m
		default:
			rest[k] = v
		}
	}
	if len(rest) > 0 {
		keys := make([]string, 0, len(rest))
		for k := range rest {
			keys = append(keys, k)
		}
		sortStrings(keys)
		return nil, fmt.Errorf("custom_medium: unknown parameter(s) %s (known: layers, distance, index, absorption, params, steps)", joinNames(keys))
	}
	return spec, nil
}

// decodeLayerList accepts a layer array as a JSON array or as a string holding
// one, and decodes it strictly (unknown fields are errors).
func decodeLayerList(v any) ([]CustomMediumLayerSpec, error) {
	var raw []byte
	switch t := v.(type) {
	case string:
		raw = []byte(t)
	case nil:
		return nil, fmt.Errorf("layers is empty")
	default:
		b, err := json.Marshal(v)
		if err != nil {
			return nil, fmt.Errorf("layers: %v", err)
		}
		raw = b
	}
	var layers []CustomMediumLayerSpec
	if err := decodeStrictJSON(raw, &layers); err != nil {
		return nil, fmt.Errorf("layers: %v", err)
	}
	if len(layers) == 0 {
		return nil, fmt.Errorf("layers is empty")
	}
	return layers, nil
}

// decodeParamMap accepts a parameter map as a JSON object or as a string
// holding one.
func decodeParamMap(v any) (map[string]float64, error) {
	var raw []byte
	switch t := v.(type) {
	case string:
		if strings.TrimSpace(t) == "" {
			return nil, nil
		}
		raw = []byte(t)
	case nil:
		return nil, nil
	default:
		b, err := json.Marshal(v)
		if err != nil {
			return nil, err
		}
		raw = b
	}
	var m map[string]float64
	if err := decodeStrictJSON(raw, &m); err != nil {
		return nil, err
	}
	return m, nil
}

// decodeStrictJSON decodes JSON with unknown fields rejected, so a misspelled
// layer key ("thicknes") is reported instead of silently defaulting to zero.
func decodeStrictJSON(raw []byte, out any) error {
	dec := json.NewDecoder(bytes.NewReader(raw))
	dec.DisallowUnknownFields()
	if err := dec.Decode(out); err != nil {
		return err
	}
	// A trailing token (two JSON values in one string) is a user mistake.
	if _, err := dec.Token(); err != io.EOF {
		return fmt.Errorf("trailing data after the JSON value")
	}
	return nil
}

func init() {
	// Registered from this file (elements.go holds the built-in registry but is
	// owned by another change in flight); the name is unique, so nothing is
	// replaced.
	RegisterElement("custom_medium", newCustomMedium)
}
