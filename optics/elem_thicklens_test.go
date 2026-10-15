package optics

import (
	"math"
	"testing"
)

// The thick lens is a two-surface lens collapsed onto one phase mask. The
// tests below check (1) the phase the element applies against the sag
// derivation of elem_thicklens.go, (2) the focus it produces under ASM
// propagation against the closed-form focal lengths, and (3) the geometry
// options (aperture, blocking, conics, dispersion).
//
// Reference numbers used throughout (N-BK7, λ = 632.8 nm, n = 1.515089):
//
//	R1 = +60 mm, R2 = −60 mm, d = 12 mm
//	f_mask = 1/((n−1)(1/R1 − 1/R2))                            = 58.24 mm
//	f_geo  = 1/((n−1)(1/R1 − 1/R2 + (n−1)d/(n·R1·R2)))          = 60.29 mm
//
// i.e. the power-stacking term of the classical formula is a 3.5% effect, big
// enough that the focus measurement can tell the two models apart.

const (
	thickTestWL  = 632.8e-9
	thickTestR1  = 0.06
	thickTestR2  = -0.06
	thickTestD   = 0.012
	thickTestApx = 0.8e-3
)

// thickTestElement builds the test lens through the public registry, so the
// constructor, the parameter names and the JSON-style map are exercised too.
func thickTestElement(t *testing.T, extra map[string]any) Element {
	t.Helper()
	p := map[string]any{
		"material": "N-BK7", "R1": thickTestR1, "R2": thickTestR2,
		"d": thickTestD, "aperture": thickTestApx,
	}
	for k, v := range extra {
		p[k] = v
	}
	el, err := NewElement(ElementSpec{Type: "thick_lens", Params: p})
	if err != nil {
		t.Fatalf("NewElement(thick_lens): %v", err)
	}
	return el
}

// TestThickLensPhaseExact checks the applied phase against an independent
// evaluation of the sag formula: t(r) = exp(i·k0·(n−1)·(d + sag2(r) − sag1(r)))
// with sag = R − sign(R)·sqrt(R² − r²) (the elementary sphere formula, written
// separately from ConicSag).
func TestThickLensPhaseExact(t *testing.T) {
	n := 128
	width := 4e-3
	wl := thickTestWL
	k0 := 2 * math.Pi / wl
	ng, _, _, err := ThickLensGlassIndex("N-BK7", wl)
	if err != nil {
		t.Fatal(err)
	}
	el := thickTestElement(t, nil)
	f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	if err := el.Apply(f, ctxFor(wl)); err != nil {
		t.Fatal(err)
	}
	sag := func(r2, R float64) float64 {
		if r2 >= R*R {
			t.Fatalf("test radius %g m outside the surface (|R| = %g)", math.Sqrt(r2), R)
		}
		return R - math.Copysign(math.Sqrt(R*R-r2), R) // R−sign(R)·sqrt(R²−r²)
	}
	checked := 0
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			x, y := f.X(i), f.Y(j)
			r2 := x*x + y*y
			if r2 > thickTestApx*thickTestApx {
				continue // outside the semi-aperture: covered separately
			}
			thick := thickTestD + sag(r2, thickTestR2) - sag(r2, thickTestR1)
			want := math.Mod(k0*(ng-1)*thick, 2*math.Pi)
			got := math.Atan2(imag(f.Ex[j*n+i]), real(f.Ex[j*n+i]))
			d := math.Remainder(got-want, 2*math.Pi)
			if math.Abs(d) > 1e-9 {
				t.Fatalf("phase at (%g, %g) m: got %g rad, want %g rad (Δ=%g)", x, y, got, want, d)
			}
			checked++
		}
	}
	if checked < 1000 {
		t.Fatalf("only %d pixels checked", checked)
	}
	// Paraxial check: for r ≪ R the mask is the thin lens
	// exp(−i·k0·r²/(2f_mask)) up to the piston exp(i·k0·(n−1)·d), and f_mask is
	// what ThickLensMaskFocalLength reports.
	fm := ThickLensMaskFocalLength(thickTestR1, thickTestR2, ng)
	if math.Abs(fm-0.058242) > 1e-5 {
		t.Fatalf("ThickLensMaskFocalLength = %g m, want ≈ 0.058242 m", fm)
	}
	fg := ThickLensGeometricFocalLength(thickTestR1, thickTestR2, thickTestD, ng)
	if math.Abs(fg-0.060292) > 1e-5 {
		t.Fatalf("ThickLensGeometricFocalLength = %g m, want ≈ 0.060292 m", fg)
	}
	if rel := (fg - fm) / fm; math.Abs(rel-0.0352) > 0.002 {
		t.Fatalf("power-stacking term is %g of f_mask, want ≈ 3.5%%", rel)
	}
}

// TestThickLensFocusScan propagates the masked field with ASM and locates the
// axial intensity maximum, comparing it with the two closed forms.
//
// Sampling: N = 512, width = 3 mm → dx = 5.86 µm, so the semi-aperture rim
// (0.8 mm) carries a phase step of 2π·r·dx/(λ·f_mask) = 0.80 rad/pixel, well
// below the π/pixel Nyquist limit of docs/PHYSICS.md §2; the focal spot
// λ·f/(2R) = 22.9 µm spans 3.9 pixels and the Rayleigh range is 3.2 mm, so the
// 0.6 mm z-scan steps resolve the axial peak. Tolerance: 1.5% of f (the scan
// plus the parabolic refinement is good to ≈0.3%; the remaining deviation is
// the focal shift of the truncated aperture, which is ≪1% at Fresnel number
// R²/(λf) = 17).
func TestThickLensFocusScan(t *testing.T) {
	n, width := 512, 3e-3
	wl := thickTestWL
	ng, _, _, err := ThickLensGlassIndex("N-BK7", wl)
	if err != nil {
		t.Fatal(err)
	}
	mask := ThickLensMaskFocalLength(thickTestR1, thickTestR2, ng)
	geom := ThickLensGeometricFocalLength(thickTestR1, thickTestR2, thickTestD, ng)
	measured := measureThickLensFocus(t, map[string]any{}, n, width, wl, mask)
	t.Logf("measured focus %.5f m: f_mask = %.5f m (Δ = %+.2f%%), f_geo = %.5f m (Δ = %+.2f%%)",
		measured, mask, 100*(measured-mask)/mask, geom, 100*(measured-geom)/geom)
	if rel := math.Abs(measured-mask) / mask; rel > 0.015 {
		t.Fatalf("focus at %g m differs from f_mask = %g m by %.2f%% (>1.5%%)", measured, mask, 100*rel)
	}
	if rel := math.Abs(measured-geom) / geom; rel < 0.02 {
		t.Fatalf("focus at %g m is within 2%% of the *geometric* focal length %g m; "+
			"the sag mask should sit at f_mask = %g m (the 3.5%% power-stacking difference must show)", measured, geom, mask)
	}
	// The diffraction-limited axial peak of a uniformly illuminated circular
	// pupil is P·πR²/(λ²f²); the measured value documents the mask quality.
	f, _ := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	el := thickTestElement(t, nil)
	if err := el.Apply(f, ctxFor(wl)); err != nil {
		t.Fatal(err)
	}
	pAp := f.Power()
	ideal := pAp * math.Pi * thickTestApx * thickTestApx / (wl * wl * mask * mask)
	peak := axialIntensityAt(t, f.Clone(), n, width, wl, measured)
	if rel := math.Abs(peak-ideal) / ideal; rel > 0.15 {
		t.Fatalf("axial peak at focus %g W/m², Airy ideal %g W/m² (rel %g)", peak, ideal, rel)
	}
	t.Logf("axial peak %.3g W/m² vs Airy ideal %.3g W/m² (%.1f%%)", peak, ideal, 100*(peak-ideal)/ideal)
}

// TestThickLensEquivalentModelFocusesAtGeometricFocalLength checks model =
// "equivalent": the element then carries the classical thick-lens power, so
// the focus sits at f_geo from the element plane (same scan, same tolerance).
func TestThickLensEquivalentModelFocusesAtGeometricFocalLength(t *testing.T) {
	n, width := 512, 3e-3
	wl := thickTestWL
	ng, _, _, err := ThickLensGlassIndex("N-BK7", wl)
	if err != nil {
		t.Fatal(err)
	}
	geom := ThickLensGeometricFocalLength(thickTestR1, thickTestR2, thickTestD, ng)
	measured := measureThickLensFocus(t, map[string]any{"model": "equivalent"}, n, width, wl, geom)
	t.Logf("equivalent model: measured %.5f m vs f_geo %.5f m (%+.2f%%)", measured, geom, 100*(measured-geom)/geom)
	if rel := math.Abs(measured-geom) / geom; rel > 0.015 {
		t.Fatalf("equivalent model focuses at %g m, want f_geo = %g m (rel %g)", measured, geom, rel)
	}
}

// measureThickLensFocus scans z around wantF with 0.6 mm steps and refines the
// axial maximum with a parabola. It returns the focus distance from the
// element plane in metres.
func measureThickLensFocus(t *testing.T, extra map[string]any, n int, width, wl, wantF float64) float64 {
	t.Helper()
	step := 0.6e-3
	nSteps := 21
	z0 := wantF - float64(nSteps/2)*step
	zs := make([]float64, 0, nSteps)
	vals := make([]float64, 0, nSteps)
	for k := 0; k < nSteps; k++ {
		z := z0 + float64(k)*step
		f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
		if err != nil {
			t.Fatal(err)
		}
		el := thickTestElement(t, extra)
		if err := el.Apply(f, ctxFor(wl)); err != nil {
			t.Fatal(err)
		}
		zs = append(zs, z)
		vals = append(vals, axialIntensityAt(t, f, n, width, wl, z))
	}
	best := 0
	for k := range vals {
		if vals[k] > vals[best] {
			best = k
		}
	}
	if best == 0 || best == len(vals)-1 {
		t.Fatalf("axial maximum on the scan boundary (z = %g m): the focus is outside [%g, %g] m",
			zs[best], zs[0], zs[len(zs)-1])
	}
	a, b, c := vals[best-1], vals[best], vals[best+1]
	den := a - 2*b + c
	off := 0.0
	if den != 0 {
		off = 0.5 * (a - c) / den
	}
	return zs[best] + off*step
}

// axialIntensityAt propagates f by z and returns the on-axis intensity |Ex|².
func axialIntensityAt(t *testing.T, f *Field, n int, width, wl, z float64) float64 {
	t.Helper()
	ctx := ctxFor(wl)
	if err := Propagate(f, z, MethodASM, ctx); err != nil {
		t.Fatal(err)
	}
	v := f.Ex[(n/2)*n+n/2]
	return real(v)*real(v) + imag(v)*imag(v)
}

// TestThickLensGeometryOptions covers the aperture/blocking behaviour, the
// conic constants, plane surfaces and determinism.
func TestThickLensGeometryOptions(t *testing.T) {
	n, width, wl := 256, 4e-3, thickTestWL
	apply := func(params map[string]any) *Field {
		f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
		if err != nil {
			t.Fatal(err)
		}
		el, err := NewElement(ElementSpec{Type: "thick_lens", Params: params})
		if err != nil {
			t.Fatal(err)
		}
		if err := el.Apply(f, ctxFor(wl)); err != nil {
			t.Fatal(err)
		}
		return f
	}
	base := func(extra map[string]any) map[string]any {
		p := map[string]any{"R1": 0.05, "R2": -0.05, "d": 0.004, "aperture": 1e-3}
		for k, v := range extra {
			p[k] = v
		}
		return p
	}
	// Default: blocked outside the semi-aperture.
	blk := apply(base(nil))
	at := func(f *Field, i, j int) complex128 { return f.Ex[j*n+i] }
	if at(blk, n/2+80, n/2) != 0 { // 80 px = 1.25 mm > 1 mm
		t.Fatal("default (block=true) transmits light outside the semi-aperture")
	}
	if at(blk, n/2+40, n/2) == 0 { // 40 px = 0.625 mm < 1 mm
		t.Fatal("block=true also blocks light inside the semi-aperture")
	}
	// block = false: light outside passes unchanged — the same value as an
	// element-free run (the source amplitude, not unit amplitude).
	ref, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	pass := apply(base(map[string]any{"block": false}))
	if v, w := at(pass, n/2+80, n/2), ref.Ex[(n/2)*n+n/2+80]; v != w {
		t.Fatalf("block=false does not pass light outside the aperture unchanged: got %v, want %v", v, w)
	}
	// Paraboloid (k = −1): sag is exactly r²/(2R); compare with the sphere it
	// differs at order r⁴, so the phases must differ.
	sph := apply(base(nil))
	par := apply(base(map[string]any{"k1": -1.0, "k2": -1.0}))
	k0 := 2 * math.Pi / wl
	ng, _, _, _ := ThickLensGlassIndex("N-BK7", wl)
	diff := false
	for i := 0; i < n; i++ {
		x := par.X(i)
		r2 := x * x
		if r2 == 0 || r2 > 1e-6 {
			continue
		}
		want := math.Mod(k0*(ng-1)*(0.004-r2/(2*0.05)-r2/(2*0.05)), 2*math.Pi)
		got := math.Atan2(imag(at(par, i, n/2)), real(at(par, i, n/2)))
		if d := math.Remainder(got-want, 2*math.Pi); math.Abs(d) > 1e-9 {
			t.Fatalf("parabolic conic phase at x = %g: got %g want %g", x, got, want)
		}
		if math.Abs(math.Remainder(got-math.Atan2(imag(at(sph, i, n/2)), real(at(sph, i, n/2))), 2*math.Pi)) > 1e-12 {
			diff = true
		}
	}
	if !diff {
		t.Fatal("sphere and paraboloid gave identical phases; the conic constant is ignored")
	}
	// Plane surfaces: R = 0, +Inf and −Inf all mean "no curvature" — the
	// element then applies only the piston exp(i·k0·(n−1)·d) (aperture 0:
	// no clear-aperture limit, so every pixel sees the same phase).
	f0 := apply(base(map[string]any{"R1": 0.0, "R2": math.Inf(1), "d": 1e-3, "aperture": 0.0}))
	wantPiston := cexpI(k0 * (ng - 1) * 1e-3)
	for i := range f0.Ex {
		want := ref.Ex[i] * wantPiston
		if math.Abs(real(f0.Ex[i]-want)) > 1e-12 || math.Abs(imag(f0.Ex[i]-want)) > 1e-12 {
			t.Fatalf("plane surfaces did not act as a pure piston (pixel %d = %v, want %v)", i, f0.Ex[i], want)
		}
	}
	fInf := apply(base(map[string]any{"R1": math.Inf(-1), "R2": 0.0, "d": 1e-3, "aperture": 0.0}))
	for i := range fInf.Ex {
		if fInf.Ex[i] != f0.Ex[i] {
			t.Fatalf("R = −Inf and R = +Inf/0 differ (pixel %d)", i)
		}
	}
	// Determinism: two identical applications agree bit for bit.
	a := apply(base(nil))
	b := apply(base(nil))
	for i := range a.Ex {
		if a.Ex[i] != b.Ex[i] {
			t.Fatalf("thick_lens is not deterministic at pixel %d", i)
		}
	}
	// Bad parameters are rejected in the constructor.
	for _, bad := range []map[string]any{
		{"material": "unobtainium"},
		{"d": -1e-3},
		{"model": "quantum"},
		{"aperture": -1.0},
		{"k1": math.NaN()},
		{"R1": 0.01, "R2": -0.05, "aperture": 0.02}, // aperture beyond the surface rim |R1|
		{"material": "constant", "n": 0},
	} {
		if _, err := NewElement(ElementSpec{Type: "thick_lens", Params: bad}); err == nil {
			t.Fatalf("thick_lens accepted %v", bad)
		}
	}
}

// TestThickLensGlassTable pins the built-in Sellmeier glasses to their
// published nd values (the fits are quoted in elem_thicklens.go).
func TestThickLensGlassTable(t *testing.T) {
	cases := []struct {
		name string
		nd   float64
	}{
		{"N-BK7", 1.51680},
		{"FusedSilica", 1.45846},
		{"N-SF5", 1.67270},
		{"N-SK16", 1.62041},
		{"F2", 1.62004},
		{"SiO2", 1.45846}, // alias
	}
	for _, c := range cases {
		n, lo, hi, err := ThickLensGlassIndex(c.name, 587.6e-9)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if math.Abs(n-c.nd) > 2e-4 {
			t.Fatalf("%s: nd = %.5f, want %.5f", c.name, n, c.nd)
		}
		if !(lo > 0 && hi > lo) {
			t.Fatalf("%s: validity range [%g, %g] m is not sane", c.name, lo, hi)
		}
	}
	if _, _, _, err := ThickLensGlassIndex("unobtainium", 500e-9); err == nil {
		t.Fatal("unknown glass must be rejected")
	}
	// Dispersion: the focal length shortens towards the blue for a positive
	// lens (n grows, so f = 1/((n−1)P) shrinks).
	f450 := ThickLensMaskFocalLength(0.06, -0.06, mustIndex(t, "N-BK7", 450e-9))
	f633 := ThickLensMaskFocalLength(0.06, -0.06, mustIndex(t, "N-BK7", 632.8e-9))
	if !(f450 < f633) {
		t.Fatalf("no chromatic shift: f(450 nm) = %g m, f(632.8 nm) = %g m", f450, f633)
	}
}

func mustIndex(t *testing.T, name string, wl float64) float64 {
	t.Helper()
	n, _, _, err := ThickLensGlassIndex(name, wl)
	if err != nil {
		t.Fatal(err)
	}
	return n
}

// TestThickLensWarnings covers the two diagnostics the element emits: a
// wavelength outside the glass' documented Sellmeier range, and light reaching
// beyond the rim of a surface (where the exact sag stops existing and the
// element clamps it).
func TestThickLensWarnings(t *testing.T) {
	// (1) Glass range: N-BK7 is quoted for 0.30–2.50 µm.
	el := thickTestElement(t, map[string]any{"material": "N-BK7"})
	ctx := ctxFor(4e-6)
	f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, 32, 1e-3, false, 4e-6)
	if err != nil {
		t.Fatal(err)
	}
	if err := el.Apply(f, ctx); err != nil {
		t.Fatal(err)
	}
	if !hasWarningCode(ctx, "thick_lens_glass_range") {
		t.Fatalf("no thick_lens_glass_range warning at 4 µm: %+v", ctx.Warnings.List())
	}
	for i := range f.Ex {
		if !finiteComplex(f.Ex[i]) {
			t.Fatalf("extrapolated glass produced a non-finite field at pixel %d", i)
		}
	}
	// (2) Sag clipping: a 2 cm-radius surface on a 6 cm-wide grid with no
	// aperture limit puts light past the rim.
	el2, err := NewElement(ElementSpec{Type: "thick_lens", Params: map[string]any{
		"material": "constant", "n": 1.5, "R1": 0.02, "R2": -0.02, "d": 2e-3}})
	if err != nil {
		t.Fatal(err)
	}
	ctx2 := ctxFor(thickTestWL)
	f2, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, 64, 0.06, false, thickTestWL)
	if err != nil {
		t.Fatal(err)
	}
	if err := el2.Apply(f2, ctx2); err != nil {
		t.Fatal(err)
	}
	if !hasWarningCode(ctx2, "thick_lens_sag_clip") {
		t.Fatalf("no thick_lens_sag_clip warning: %+v", ctx2.Warnings.List())
	}
	// ... and no warning when the aperture keeps the light inside the rim.
	el3, err := NewElement(ElementSpec{Type: "thick_lens", Params: map[string]any{
		"material": "constant", "n": 1.5, "R1": 0.02, "R2": -0.02, "d": 2e-3, "aperture": 0.015}})
	if err != nil {
		t.Fatal(err)
	}
	ctx3 := ctxFor(thickTestWL)
	f3, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, 64, 0.06, false, thickTestWL)
	if err != nil {
		t.Fatal(err)
	}
	if err := el3.Apply(f3, ctx3); err != nil {
		t.Fatal(err)
	}
	if hasWarningCode(ctx3, "thick_lens_sag_clip") {
		t.Fatal("sag-clip warning fired although the aperture keeps the light inside the rim")
	}
}

// hasWarningCode reports whether a run recorded a warning with the given code.
func hasWarningCode(ctx *Context, code string) bool {
	for _, w := range ctx.Warnings.List() {
		if w.Code == code {
			return true
		}
	}
	return false
}

// TestThickLensInSceneAndTrain wires the element into both execution paths:
// the legacy element train (Config.Elements) and a positioned scene, so the
// registry, the dispatcher and the scene JSON all see it.
func TestThickLensInSceneAndTrain(t *testing.T) {
	wl := thickTestWL
	ng := mustIndex(t, "N-BK7", wl)
	fm := ThickLensMaskFocalLength(thickTestR1, thickTestR2, ng)

	// (a) Legacy train.
	cfg := scalarCfg(512, 3e-3, wl)
	cfg.Source = planeSource()
	cfg.Elements = []ElementSpec{
		{Type: "thick_lens", Params: map[string]any{
			"material": "N-BK7", "R1": thickTestR1, "R2": thickTestR2,
			"d": thickTestD, "aperture": thickTestApx}},
		{Type: "propagate", Params: map[string]any{"distance": fm}},
		{Type: "sensor", Params: map[string]any{"label": "focus"}},
	}
	res := mustSim(t, cfg)
	if len(res.Planes) != 1 {
		t.Fatalf("train run recorded %d planes", len(res.Planes))
	}
	// The detector sees the aperture-truncated beam: power = P·πR²/W².
	pAp := 1e-3 * math.Pi * thickTestApx * thickTestApx / (3e-3 * 3e-3)
	if rel := math.Abs(res.Planes[0].Stats.Power-pAp) / pAp; rel > 1e-3 {
		t.Fatalf("focused power %g W, want %g W (rel %g)", res.Planes[0].Stats.Power, pAp, rel)
	}
	ideal := pAp * math.Pi * thickTestApx * thickTestApx / (wl * wl * fm * fm)
	if rel := math.Abs(res.Planes[0].Stats.Peak-ideal) / ideal; rel > 0.15 {
		t.Fatalf("scene peak %g W/m², Airy ideal %g W/m² (rel %g)", res.Planes[0].Stats.Peak, ideal, rel)
	}

	// (b) Positioned scene: same optics, geometry-driven distances. In a scene
	// the component outline *is* the clear aperture (buildElements zeroes the
	// element's own `aperture`, exactly as for `lens`), so the drawn shape uses
	// the same semi-aperture as the train case above.
	scene := &SceneSpec{Components: []ComponentSpec{
		{ID: "lens", Type: "thick_lens", Label: "厚透镜", Pos: v3(0, 0, 0),
			Shape: CircleOutline(thickTestApx),
			Params: map[string]any{
				"material": "N-BK7", "R1": thickTestR1, "R2": thickTestR2,
				"d": thickTestD, "aperture": thickTestApx}},
		{ID: "det", Type: "sensor", Label: "焦面", Pos: v3(0, 0, fm), Yaw: math.Pi,
			Shape: CircleOutline(3e-3)},
	}}
	scfg := scalarCfg(512, 3e-3, wl)
	scfg.Scene = scene
	srcPos, srcDir := v3(0, 0, -0.02), v3(0, 0, 1)
	scfg.Sources = []SourceSpec{{ID: "src", Type: "plane", Pos: &srcPos, Dir: &srcDir,
		Params: map[string]any{"power": 1e-3}}}
	sres := mustSim(t, scfg)
	if len(sres.Planes) != 1 {
		t.Fatalf("scene run recorded %d planes", len(sres.Planes))
	}
	if sres.Planes[0].Stats.Power <= 0 {
		t.Fatal("scene detector got no light through the thick lens")
	}
	if rel := math.Abs(sres.Planes[0].Stats.Peak-ideal) / ideal; rel > 0.15 {
		t.Fatalf("scene peak %g W/m², Airy ideal %g W/m² (rel %g)", sres.Planes[0].Stats.Peak, ideal, rel)
	}
}
