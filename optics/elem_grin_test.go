package optics

import (
	"encoding/json"
	"math"
	"testing"
)

// The GRIN plate applies t = exp(i·k0·d·(n0 − 1 + Δn(x,y))): the transverse
// phase k0·Δn·d plus the air-contrast piston of the slab. The tests check the
// pointwise phase against the polynomial (and the expression form against the
// polynomial), that a purely quadratic profile is indistinguishable from a thin
// lens of f = −1/(2·c2·d) after ASM propagation, and the aperture/error paths.

// grinTestElement builds the element through the registry.
func grinTestElement(t *testing.T, extra map[string]any) Element {
	t.Helper()
	el, err := NewElement(ElementSpec{Type: "grin_plate", Params: extra})
	if err != nil {
		t.Fatalf("NewElement(grin_plate): %v", err)
	}
	return el
}

func TestGRINPlatePhaseExact(t *testing.T) {
	n, width, wl := 128, 4e-3, 632.8e-9
	k0 := 2 * math.Pi / wl
	params := map[string]any{"d": 2e-3, "n0": 1.6, "c2": -60.0, "c4": 1e4, "c6": -2e7}
	el := grinTestElement(t, params)
	f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	if err := el.Apply(f, ctxFor(wl)); err != nil {
		t.Fatal(err)
	}
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			x, y := f.X(i), f.Y(j)
			r2 := x*x + y*y
			dn := r2 * (-60.0 + r2*(1e4+r2*(-2e7)))
			want := math.Mod(k0*2e-3*(1.6-1+dn), 2*math.Pi)
			got := math.Atan2(imag(f.Ex[j*n+i]), real(f.Ex[j*n+i]))
			if d := math.Remainder(got-want, 2*math.Pi); math.Abs(d) > 1e-9 {
				t.Fatalf("phase at (%g, %g) m: got %g, want %g (Δ=%g)", x, y, got, want, d)
			}
		}
	}
	// The piston of a quadratic-free plate is just k0·(n0−1)·d.
	flat, err := NewElement(ElementSpec{Type: "grin_plate", Params: map[string]any{
		"d": 1e-3, "n0": 1.5, "c2": 0.0, "c4": 0.0, "c6": 0.0}})
	if err != nil {
		t.Fatal(err)
	}
	g, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, 8, 1e-3, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	if err := flat.Apply(g, ctxFor(wl)); err != nil {
		t.Fatal(err)
	}
	wantPhase := k0 * 0.5 * 1e-3
	for i := range g.Ex {
		// The piston is ~5e3 rad, so the tolerance is limited by the absolute
		// accuracy of the phase computation, not by the element.
		if d := math.Remainder(math.Atan2(imag(g.Ex[i]), real(g.Ex[i]))-wantPhase, 2*math.Pi); math.Abs(d) > 1e-9 {
			t.Fatalf("uniform plate piston off by %g rad", d)
		}
	}
}

// TestGRINPlateEqualsThinLens: a purely quadratic Δn is a thin lens with
// f = −1/(2·c2·d); after the same propagation the two fields must agree
// (relative L2 error of the complex field), and the focal position must match
// the analytic focal length.
func TestGRINPlateEqualsThinLens(t *testing.T) {
	n, width, wl := 512, 3e-3, 632.8e-9
	f0 := 0.25
	dThick := 3e-3
	c2 := -1 / (2 * f0 * dThick) // f = −1/(2·c2·d)
	ap := 1.2e-3
	if got := GRINParaxialFocalLength(c2, dThick); math.Abs(got-f0) > 1e-9 {
		t.Fatalf("GRINParaxialFocalLength = %g, want %g", got, f0)
	}
	build := func(kind string) *Field {
		f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
		if err != nil {
			t.Fatal(err)
		}
		var el Element
		if kind == "grin" {
			el = grinTestElement(t, map[string]any{"d": dThick, "n0": 1.0, "c2": c2, "aperture": ap})
		} else {
			el, err = NewElement(ElementSpec{Type: "lens", Params: map[string]any{"f": f0, "aperture": ap}})
			if err != nil {
				t.Fatal(err)
			}
		}
		if err := el.Apply(f, ctxFor(wl)); err != nil {
			t.Fatal(err)
		}
		return f
	}
	grin, lens := build("grin"), build("lens")
	// The GRIN element is a *slab*: it also carries the piston of n0 = 1? With
	// n0 = 1 the piston vanishes, so the two fields carry exactly the same
	// phase profile inside the aperture.
	var num, den float64
	for i := range grin.Ex {
		num += complexAbs(grin.Ex[i] - lens.Ex[i])
		den += complexAbs(lens.Ex[i])
	}
	if rel := num / den; rel > 1e-12 {
		t.Fatalf("GRIN quadratic plate differs from the equivalent thin lens by %g (relative L2)", rel)
	}
	// Both must focus at the same place, and that place is f0.
	zg := measureFocus(t, func() *Field { return build("grin") }, n, width, wl, f0)
	zl := measureFocus(t, func() *Field { return build("lens") }, n, width, wl, f0)
	t.Logf("focus: GRIN plate %.5f m, thin lens %.5f m, analytic f = %.5f m", zg, zl, f0)
	if rel := math.Abs(zg-zl) / f0; rel > 0.01 {
		t.Fatalf("GRIN plate focuses at %g m, thin lens at %g m (rel %g)", zg, zl, rel)
	}
	if rel := math.Abs(zg-f0) / f0; rel > 0.02 {
		t.Fatalf("GRIN plate focus %g m differs from f = %g m by %g", zg, f0, rel)
	}
}

// measureFocus scans the axial intensity around wantF with 2 mm steps and
// refines with a parabola; build() must return a freshly masked field (the
// propagation mutates it).
func measureFocus(t *testing.T, build func() *Field, n int, width, wl, wantF float64) float64 {
	t.Helper()
	step := 2e-3
	steps := 21
	z0 := wantF - 0.5*step*float64(steps-1)
	best, bestV := -1, -1.0
	zs := make([]float64, steps)
	vals := make([]float64, steps)
	for k := 0; k < steps; k++ {
		z := z0 + float64(k)*step
		f := build()
		if err := Propagate(f, z, MethodASM, ctxFor(wl)); err != nil {
			t.Fatal(err)
		}
		v := f.Ex[(n/2)*n+n/2]
		val := real(v)*real(v) + imag(v)*imag(v)
		zs[k], vals[k] = z, val
		if val > bestV {
			best, bestV = k, val
		}
	}
	if best <= 0 || best >= steps-1 {
		t.Fatalf("axial maximum on the scan boundary (z = %g m)", zs[best])
	}
	a, b, c := vals[best-1], vals[best], vals[best+1]
	off := 0.0
	if den := a - 2*b + c; den != 0 {
		off = 0.5 * (a - c) / den
	}
	return zs[best] + off*step
}

// TestGRINPlateExpression: the expression form must reproduce the polynomial
// (same phase, deterministic) and be validated at construction.
func TestGRINPlateExpression(t *testing.T) {
	n, width, wl := 128, 4e-3, 632.8e-9
	poly := grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.0, "c2": -50, "c4": 2e4})
	expr := grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.0, "expr": "c2*r*r + c4*r*r*r*r"})
	// The expression element must be built with the same coefficients bound:
	// rebuild with explicit c2/c4 so both share the numbers.
	expr = grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.0, "c2": -50, "c4": 2e4, "expr": "c2*r*r + c4*r*r*r*r"})
	a, b := twoFields(t, n, width, wl, poly), twoFields(t, n, width, wl, expr)
	var num, den float64
	for i := range a.Ex {
		num += complexAbs(a.Ex[i] - b.Ex[i])
		den += complexAbs(a.Ex[i])
	}
	if rel := num / den; rel > 1e-12 {
		t.Fatalf("expression Δn differs from the polynomial by %g (relative L2)", rel)
	}
	// The expression may use the polar angle (th) and r; the element must feed
	// them per point (an unfilled env would silently read th = 0).
	azim := grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.0, "c2": -40.0,
		"expr": "c2*r*r*(1 + 0.5*cos(2*th))"})
	fa := twoFields(t, n, width, wl, azim)
	poly2 := twoFields(t, n, width, wl, grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.0, "c2": -40.0}))
	same := true
	for i := range fa.Ex {
		if fa.Ex[i] != poly2.Ex[i] {
			same = false
			break
		}
	}
	if same {
		t.Fatal("an azimuthal expression gave the same field as the plain polynomial: th/r are not fed per point")
	}
	// A custom expression that the polynomial cannot express (Gaussian dip).
	g := grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.0,
		"expr": "c2*(1 - exp(-r*r/(8e-7)))", "c2": -80.0})
	f := twoFields(t, n, width, wl, g)
	if f.Power() <= 0 {
		t.Fatal("expression element produced an empty field")
	}
	// Determinism.
	c := twoFields(t, n, width, wl, g)
	for i := range c.Ex {
		if c.Ex[i] != f.Ex[i] {
			t.Fatalf("grin_plate expression path is not deterministic at pixel %d", i)
		}
	}
	// Bad parameters.
	for _, bad := range []map[string]any{
		{"d": 0.0},
		{"d": -1e-3},
		{"n0": 0.0},
		{"c2": math.Inf(1)},
		{"aperture": -1.0},
		{"expr": "1/0"},          // non-finite at the probe
		{"expr": "unknownvar*2"}, // unknown identifier
		{"expr": "r*("},          // syntax error
	} {
		if _, err := NewElement(ElementSpec{Type: "grin_plate", Params: bad}); err == nil {
			t.Fatalf("grin_plate accepted %v", bad)
		}
	}
}

// twoFields builds a plane wave and applies el to it.
func twoFields(t *testing.T, n int, width, wl float64, el Element) *Field {
	t.Helper()
	f, err := BuildSource(SourceSpec{Type: "plane", Params: map[string]any{"power": 1e-3}}, n, width, false, wl)
	if err != nil {
		t.Fatal(err)
	}
	if err := el.Apply(f, ctxFor(wl)); err != nil {
		t.Fatal(err)
	}
	return f
}

// TestGRINPlateApertureAndPower: the aperture blocks (or not) and the element
// is lossless inside the aperture.
func TestGRINPlateApertureAndPower(t *testing.T) {
	n, width, wl := 256, 4e-3, 632.8e-9
	base := map[string]any{"d": 1e-3, "n0": 1.5, "c2": -30.0}
	ap := 1e-3
	withAp := map[string]any{"d": 1e-3, "n0": 1.5, "c2": -30.0, "aperture": ap}
	el := grinTestElement(t, withAp)
	f := twoFields(t, n, width, wl, el)
	at := func(i, j int) complex128 { return f.Ex[j*n+i] }
	if at(0, 0) != 0 {
		t.Fatal("light outside the semi-aperture should be blocked by default")
	}
	// Power: the transmitted power is the fraction of the grid inside the disk
	// (the blocked light is removed, there is no other loss).
	frac := math.Pi * ap * ap / (width * width)
	if rel := math.Abs(f.Power()/(1e-3*frac) - 1); rel > 0.01 {
		t.Fatalf("transmitted power %g W, want ≈ %g W (rel %g)", f.Power(), 1e-3*frac, rel)
	}
	// block = false: outside is untouched.
	pass := base
	pass["aperture"] = ap
	pass["block"] = false
	el2 := grinTestElement(t, pass)
	g := twoFields(t, n, width, wl, el2)
	ref := twoFields(t, n, width, wl, grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.0, "c2": 0.0}))
	if g.Ex[0] != ref.Ex[0] {
		t.Fatalf("block=false changed the field outside the aperture: %v vs %v", g.Ex[0], ref.Ex[0])
	}
	// A soft edge (edge_sigma) keeps the field finite everywhere and, as the
	// smoothing length goes to zero, converges back to the hard edge. (A wide
	// erf edge does change the transmitted power: the fade mixes the steered
	// and unsteered field, which is interference, not absorption — so the
	// convergence check uses a sub-pixel sigma.)
	soft := grinTestElement(t, map[string]any{"d": 1e-3, "n0": 1.5, "c2": -30.0,
		"aperture": ap, "edge_sigma": 2e-5})
	fs := twoFields(t, n, width, wl, soft)
	for i := range fs.Ex {
		if !finiteComplex(fs.Ex[i]) {
			t.Fatalf("soft-edged GRIN plate produced a non-finite field at pixel %d", i)
		}
	}
	// The residual difference is the pixelization of the rim (sigma ≈ 1.3 px
	// blends the partially covered boundary pixels), not a physics difference.
	if rel := math.Abs(fs.Power()-f.Power()) / f.Power(); rel > 0.05 {
		t.Fatalf("sub-pixel soft edge transmits %g W, hard edge %g W (rel %g)", fs.Power(), f.Power(), rel)
	}
}

// TestNewRefractiveElementsRegisteredAndDocumented checks the registration and
// the catalog entry of all three elements: they must be built by the registry,
// appear in /api/catalog (BuildCatalog) with their parameters, and be
// classified for the GUI.
func TestNewRefractiveElementsRegisteredAndDocumented(t *testing.T) {
	names := []string{"thick_lens", "coating", "grin_plate"}
	registered := map[string]bool{}
	for _, r := range RegisteredElements() {
		registered[r] = true
	}
	cat := BuildCatalog()
	byType := map[string]*ElementDoc{}
	for i := range cat.Elements {
		byType[cat.Elements[i].Type] = &cat.Elements[i]
	}
	for _, name := range names {
		if !registered[name] {
			t.Fatalf("%s is not in the element registry", name)
		}
		doc, ok := byType[name]
		if !ok {
			t.Fatalf("%s has no catalog entry (missing ElementDocs)", name)
		}
		if doc.Label == "" || doc.Help == "" {
			t.Fatalf("%s catalog entry needs a label and a help text", name)
		}
		if len(doc.Params) == 0 {
			t.Fatalf("%s catalog entry has no parameters", name)
		}
		for _, p := range doc.Params {
			if p.Key == "" || p.Kind == "" {
				t.Fatalf("%s: parameter %+v is incomplete", name, p)
			}
			if p.Default == nil && p.Kind != "bool" {
				t.Fatalf("%s: parameter %q has no default (the GUI fills inserted elements from it)", name, p.Key)
			}
		}
		if _, ok := ComponentClasses()[name]; !ok {
			t.Fatalf("%s is missing from ComponentClasses", name)
		}
		// The catalog must be complete enough to build the element from the
		// documented defaults (that is what the GUI's insert button does).
		params := map[string]any{}
		for _, p := range doc.Params {
			if p.Kind == "float" || p.Kind == "int" || p.Kind == "choice" {
				params[p.Key] = p.Default
			}
		}
		if _, err := NewElement(ElementSpec{Type: name, Params: params}); err != nil {
			t.Fatalf("%s cannot be built from its documented defaults: %v", name, err)
		}
	}
	// The catalog parameter keys must match what the elements read: a wrong
	// key would silently fall back to a default.
	wantKeys := map[string][]string{
		"thick_lens": {"material", "R1", "R2", "d", "aperture", "block", "k1", "k2", "model", "edge_sigma", "n"},
		"coating":    {"layers", "substrate", "incident", "theta0", "plane"},
		"grin_plate": {"d", "n0", "c2", "c4", "c6", "expr", "aperture", "block", "edge_sigma"},
	}
	for name, keys := range wantKeys {
		doc := byType[name]
		have := map[string]bool{}
		for _, p := range doc.Params {
			have[p.Key] = true
		}
		for _, k := range keys {
			if !have[k] {
				t.Fatalf("%s: catalog is missing parameter %q", name, k)
			}
		}
	}
}

// TestNewRefractiveElementsSceneJSON writes the three elements into a scene
// JSON document and routes it, so the scene path (BuildSceneGraph,
// TraceScene, Simulate) is exercised end to end.
func TestNewRefractiveElementsSceneJSON(t *testing.T) {
	doc := `{
	  "components": [
	    {"id":"l","type":"thick_lens","pos":{"x":0,"y":0,"z":0},"shape":{"kind":"circle","params":{"radius":1.5e-3}},
	     "params":{"material":"N-BK7","R1":2.0,"R2":-2.0,"d":0.012}},
	    {"id":"c","type":"coating","pos":{"x":0,"y":0,"z":0.02},"shape":{"kind":"circle","params":{"radius":1.5e-3}},
	     "params":{"substrate":"N-BK7","layers":[{"material":"MgF2","thickness_nm":99.7}]}},
	    {"id":"g","type":"grin_plate","pos":{"x":0,"y":0,"z":0.04},"shape":{"kind":"circle","params":{"radius":1.5e-3}},
	     "params":{"d":1e-3,"n0":1.0,"c2":-20.0}},
	    {"id":"det","type":"sensor","label":"out","pos":{"x":0,"y":0,"z":0.20},"yaw":3.14159265358979,
	     "shape":{"kind":"circle","params":{"radius":3e-3}}}
	  ]
	}`
	var scene SceneSpec
	if err := json.Unmarshal([]byte(doc), &scene); err != nil {
		t.Fatalf("scene JSON: %v", err)
	}
	tr, err := TraceScene(&scene, []SourceSpec{{Type: "plane"}}, 632.8e-9)
	if err != nil {
		t.Fatalf("TraceScene: %v", err)
	}
	if len(tr.Components) != 4 {
		t.Fatalf("trace has %d components", len(tr.Components))
	}
	cfg := scalarCfg(256, 4e-3, 632.8e-9)
	cfg.Scene = &scene
	srcPos, srcDir := v3(0, 0, -0.02), v3(0, 0, 1)
	cfg.Sources = []SourceSpec{{ID: "src", Type: "plane", Pos: &srcPos, Dir: &srcDir,
		Params: map[string]any{"power": 1e-3}}}
	res := mustSim(t, cfg)
	if len(res.Planes) != 1 {
		t.Fatalf("scene run recorded %d planes", len(res.Planes))
	}
	// Loss budget: the coating transmits T (QFR at the run wavelength ~ MgF2
	// QW designed at 550 nm, evaluated at 632.8 nm), the rest is lossless
	// (the GRIN plate has n0 = 1 and a quadratic term only).
	el, err := NewElement(ElementSpec{Type: "coating", Params: map[string]any{
		"substrate": "N-BK7", "layers": []any{map[string]any{"material": "MgF2", "thickness_nm": 99.7}}}})
	if err != nil {
		t.Fatal(err)
	}
	r, err := el.(*coatingEl).response(632.8e-9)
	if err != nil {
		t.Fatal(err)
	}
	// Expected: the source power times the fraction of the grid inside the
	// drawn 1.5 mm aperture (the scene outline *is* the clear aperture) times
	// the coating's T. Tolerance 3%: the three hard-edged masks diffract a
	// little light into evanescent components, which ASM does not carry.
	gridFrac := math.Pi * 1.5e-3 * 1.5e-3 / (4e-3 * 4e-3)
	want := 1e-3 * gridFrac * r.TsPow
	if rel := math.Abs(res.Planes[0].Stats.Power-want) / want; rel > 0.03 {
		t.Fatalf("scene power %g W, want %g W (T = %g, grid fraction %g)", res.Planes[0].Stats.Power, want, r.TsPow, gridFrac)
	}
	t.Logf("scene power %.6g W vs expected %.6g W (T = %.5f)", res.Planes[0].Stats.Power, want, r.TsPow)
}
