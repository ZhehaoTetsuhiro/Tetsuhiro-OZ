package optics

import (
	"fmt"
	"math"
	"strings"
	"testing"
)

// extendedScene builds a scene whose only source is an extended (incoherent)
// emitter at z = -0.30, illuminating a slit at -0.15 and a wide detector at
// +0.15. The default geometry is a genuinely far-field one (a 0.4 mm slit,
// detector 300 mm behind it), so the pattern is the textbook single-slit
// Fraunhofer shape.
func extendedScene(t *testing.T, params map[string]any, comps ...ComponentSpec) Config {
	t.Helper()
	pos := v3(0, 0, -0.30)
	dir := v3(0, 0, 1)
	if len(comps) == 0 {
		comps = []ComponentSpec{
			{Type: "aperture", Label: "单缝", Pos: v3(0, 0, -0.15),
				Shape: RectOutline(4e-4, 0.02)},
			bigDet("屏", v3(0, 0, 0.15), math.Pi, 0.03),
		}
	}
	return sceneCfg(256, 0.03, []SourceSpec{{
		ID: "src", Label: "扩展光源", Type: "extended",
		Pos: &pos, Dir: &dir, Params: params,
	}}, comps...)
}

// TestExtendedSourceExpandsToPointSources pins the expansion itself: an
// extended emitter becomes mutually incoherent point sources, one coherent
// unit per sample, with the total power split between them.
func TestExtendedSourceExpandsToPointSources(t *testing.T) {
	src := SourceSpec{ID: "src", Type: "extended", Params: map[string]any{
		"width": 1e-2, "height": 1e-3, "samples_x": 4, "samples_y": 2,
		"power": 2e-3, "standoff": 0.1,
	}}
	scene := &SceneSpec{Components: []ComponentSpec{
		{Type: "aperture", Label: "缝", Pos: v3(0, 0, 0.15)},
		{Type: "sensor", Label: "屏", Pos: v3(0, 0, 0.45)},
	}}
	out, warns := ExpandExtendedSources([]SourceSpec{src}, scene)
	if len(warns) != 0 {
		t.Fatalf("unexpected warnings: %v", warns)
	}
	if len(out) != 8 {
		t.Fatalf("4x2 samples: got %d sources, want 8", len(out))
	}
	var power float64
	seen := map[string]bool{}
	for i, s := range out {
		if s.Type != "spherical" {
			t.Errorf("source %d: type %q, want spherical (point emitter)", i, s.Type)
		}
		if s.Group != "" {
			t.Errorf("source %d: group %q, want empty (every point is its own incoherent unit)", i, s.Group)
		}
		if seen[s.GroupID(i)] {
			t.Errorf("source %d: duplicate unit key %q (points must be mutually incoherent)", i, s.GroupID(i))
		}
		seen[s.GroupID(i)] = true
		if r := pfd(s.Params, "radius", 0); r <= 0 {
			t.Errorf("source %d: radius %g must be > 0 (the wavefront is sampled R in front of the point)", i, r)
		}
		power += s.Power()
	}
	if math.Abs(power-2e-3)/2e-3 > 1e-12 {
		t.Errorf("total power after expansion = %.9g W, want 2e-3 W", power)
	}
}

// TestExtendedSourceExpansionGeometry pins where the points land: the sample
// centres are the cell centres of the emitter, and every virtual point sits
// on the emitter plane (the source is R in front of it with radius R).
func TestExtendedSourceExpansionGeometry(t *testing.T) {
	out, warns := ExpandExtendedSources([]SourceSpec{{
		Type: "extended", Params: map[string]any{
			"width": 1e-2, "height": 2e-3, "samples_x": 2, "samples_y": 2,
			"power": 1e-3, "standoff": 0.1,
		},
	}}, &SceneSpec{Components: []ComponentSpec{{Type: "sensor", Pos: v3(0, 0, 0.2)}}})
	if len(warns) != 0 {
		t.Fatalf("unexpected warnings: %v", warns)
	}
	if len(out) != 4 {
		t.Fatalf("got %d points, want 4", len(out))
	}
	want := [][2]float64{{-2.5e-3, -5e-4}, {2.5e-3, -5e-4}, {-2.5e-3, 5e-4}, {2.5e-3, 5e-4}}
	for i, s := range out {
		p := s.Position()
		if math.Abs(p.X-want[i][0]) > 1e-12 || math.Abs(p.Y-want[i][1]) > 1e-12 {
			t.Errorf("point %d at (%g, %g), want (%g, %g)", i, p.X, p.Y, want[i][0], want[i][1])
		}
		// The source sits R in front of the emitter plane along +z.
		if math.Abs(p.Z-0.1) > 1e-12 {
			t.Errorf("point %d at z=%g, want the standoff 0.1", i, p.Z)
		}
		if r := pfd(s.Params, "radius", 0); math.Abs(r-0.1) > 1e-12 {
			t.Errorf("point %d radius %g, want the standoff 0.1", i, r)
		}
	}
}

// TestExtendedSourceMatchesSinglePoint is the continuity check: a 1x1
// sampling must be exactly the plain spherical point source it stands for
// (same position, same radius, same power).
func TestExtendedSourceMatchesSinglePoint(t *testing.T) {
	const standoff = 0.14
	point := extendedScene(t, map[string]any{"power": 1e-3})
	// The extended source samples the wavefront R in front of the emitter
	// point, so its exact plain-source equivalent is a spherical source R in
	// front of the emitter plane with the same radius R.
	point.Sources = []SourceSpec{{
		ID: "src", Type: "spherical", Pos: &Vec3{Z: -0.30 + standoff}, Dir: &Vec3{Z: 1},
		Params: map[string]any{"power": 1e-3, "radius": standoff},
	}}
	wide := extendedScene(t, map[string]any{
		"width": 0, "height": 0, "samples_x": 1, "samples_y": 1, "power": 1e-3,
		"standoff": standoff,
	})
	a, err := Simulate(point)
	if err != nil {
		t.Fatal(err)
	}
	b, err := Simulate(wide)
	if err != nil {
		t.Fatal(err)
	}
	if len(a.Planes) != 1 || len(b.Planes) != 1 {
		t.Fatalf("want one plane each, got %v / %v", planeLabels(a), planeLabels(b))
	}
	pa, pb := a.Planes[0].Stats.Power, b.Planes[0].Stats.Power
	if math.Abs(pa-pb)/pa > 1e-9 {
		t.Errorf("1-sample extended source: power %.9g W, point source %.9g W", pb, pa)
	}
}

// TestExtendedSourceIsIncoherent is the physics pin: two emitter points that
// would interfere if they were coherent must simply add in intensity, so the
// detector power is the sum of the two single-point powers — not the coherent
// sum, which would show fringes and a different total.
func TestExtendedSourceIsIncoherent(t *testing.T) {
	one := extendedScene(t, map[string]any{
		"width": 0, "height": 0, "samples_x": 1, "samples_y": 1, "power": 5e-4,
		"standoff": 0.14,
	})
	two := extendedScene(t, map[string]any{
		"width": 2e-3, "height": 0, "samples_x": 2, "samples_y": 1, "power": 1e-3,
		"standoff": 0.14,
	})
	r1, err := Simulate(one)
	if err != nil {
		t.Fatal(err)
	}
	r2, err := Simulate(two)
	if err != nil {
		t.Fatal(err)
	}
	_ = r1
	// Each of the two off-axis points on its own: the two-point power must be
	// the sum of the two single-point powers (an incoherent pair adds
	// intensities, whatever each point delivers).
	left := extendedScene(t, map[string]any{"power": 1e-3})
	left.Sources = []SourceSpec{{ID: "src", Type: "spherical",
		Pos: &Vec3{X: -5e-4, Z: -0.16}, Dir: &Vec3{Z: 1},
		Params: map[string]any{"power": 5e-4, "radius": 0.14}}}
	right := extendedScene(t, map[string]any{"power": 1e-3})
	right.Sources = []SourceSpec{{ID: "src", Type: "spherical",
		Pos: &Vec3{X: 5e-4, Z: -0.16}, Dir: &Vec3{Z: 1},
		Params: map[string]any{"power": 5e-4, "radius": 0.14}}}
	rl, err := Simulate(left)
	if err != nil {
		t.Fatal(err)
	}
	rr, err := Simulate(right)
	if err != nil {
		t.Fatal(err)
	}
	pl, pr, p2 := rl.Planes[0].Stats.Power, rr.Planes[0].Stats.Power, r2.Planes[0].Stats.Power
	if math.Abs(p2-(pl+pr))/(pl+pr) > 0.02 {
		t.Errorf("two incoherent emitters: power %.6g W, want the sum %.6g W", p2, pl+pr)
	}
	if len(r2.Planes[0].Parts) != 2 {
		t.Errorf("two emitters: %d parts, want 2 (one incoherent unit per point)", len(r2.Planes[0].Parts))
	}
	// The same two points in one coherent group interfere instead: the total
	// power differs from the incoherent sum (fringes move light between the
	// bright and dark pixels, and the slit clips them differently).
	coherent := extendedScene(t, map[string]any{"power": 1e-3})
	coherent.Sources = []SourceSpec{
		{ID: "a", Type: "spherical", Pos: &Vec3{X: -5e-4, Z: -0.16}, Dir: &Vec3{Z: 1}, Group: "g",
			Params: map[string]any{"power": 5e-4, "radius": 0.14}},
		{ID: "b", Type: "spherical", Pos: &Vec3{X: 5e-4, Z: -0.16}, Dir: &Vec3{Z: 1}, Group: "g",
			Params: map[string]any{"power": 5e-4, "radius": 0.14}},
	}
	rc, err := Simulate(coherent)
	if err != nil {
		t.Fatal(err)
	}
	if len(rc.Planes[0].Parts) != 1 {
		t.Errorf("a coherent pair must be one unit, got %d parts", len(rc.Planes[0].Parts))
	}
	if math.Abs(rc.Planes[0].Stats.Power-p2)/p2 < 0.02 {
		t.Errorf("coherent pair power %.6g W equals the incoherent sum %.6g W: the two models are indistinguishable",
			rc.Planes[0].Stats.Power, p2)
	}
}

// TestExtendedSourcePowerSplitIsEven pins the sampling weights: a uniform
// profile gives every point the same power, a gaussian one weights the centre.
func TestExtendedSourcePowerSplitIsEven(t *testing.T) {
	mk := func(profile string) []SourceSpec {
		out, _ := ExpandExtendedSources([]SourceSpec{{
			Type: "extended", Params: map[string]any{
				"width": 1e-2, "height": 1e-2, "samples_x": 3, "samples_y": 3,
				"power": 9e-3, "profile": profile, "standoff": 0.1,
			},
		}}, &SceneSpec{Components: []ComponentSpec{{Type: "sensor", Pos: v3(0, 0, 0.2)}}})
		return out
	}
	uni := mk("uniform")
	for i, s := range uni {
		if math.Abs(s.Power()-1e-3)/1e-3 > 1e-12 {
			t.Errorf("uniform: point %d carries %.9g W, want 1e-3 W", i, s.Power())
		}
	}
	// A gaussian profile exp(-2 r^2/w^2) with w the 1/e^2 half-extent: the
	// centre sample (index 4 of a 3x3 grid) carries more than a corner
	// sample, and the weights still sum to the total power.
	gauss := mk("gaussian")
	centre := gauss[4].Power()
	corner := gauss[0].Power()
	if !(centre > corner) {
		t.Errorf("gaussian: centre %.6g W should dominate the corner %.6g W", centre, corner)
	}
	// The 3x3 cell centres sit at 0 and ±w/3 (the cells are w/3 wide), so with
	// the 1/e^2 half-extent convention (w/2) the raw weight is
	// exp(-2 r^2/(w/2)^2) = exp(-8 r^2/w^2): centre 1, edge exp(-8/9) = 0.411,
	// corner exp(-16/9) = 0.165. The normalisation divides all of them by the
	// same sum, so these are the exact ratios.
	if got, want := gauss[1].Power()/centre, math.Exp(-8.0/9); math.Abs(got-want) > 1e-12 {
		t.Errorf("gaussian edge weight ratio %.9g, want %.9g", got, want)
	}
	if got, want := corner/centre, math.Exp(-16.0/9); math.Abs(got-want) > 1e-12 {
		t.Errorf("gaussian corner weight ratio %.9g, want %.9g", got, want)
	}
	var sum float64
	for _, s := range gauss {
		sum += s.Power()
	}
	if math.Abs(sum-9e-3)/9e-3 > 1e-12 {
		t.Errorf("gaussian weights sum to %.9g W, want 9e-3 W", sum)
	}
}

// TestExtendedSourceLimits: absurd sample counts must be refused before the
// run starts (each point is a full simulation).
func TestExtendedSourceLimits(t *testing.T) {
	cfg := extendedScene(t, map[string]any{
		"width": 1e-2, "height": 1e-2, "samples_x": 500, "samples_y": 500,
		"power": 1e-3,
	})
	if issues := ValidateConfig(&cfg); len(issues) == 0 {
		t.Fatal("500x500 point emitters passed validation")
	} else if !strings.Contains(issues[0].Message, "samples_x") {
		t.Errorf("validation message %q does not mention the sample count", issues[0].Message)
	}
	if _, err := Simulate(cfg); err == nil {
		t.Error("Simulate accepted an over-large extended source")
	}
}

// TestExtendedSourceManyUnitsKeepsTotalIntensity guards the memory policy: a
// run with more incoherent units than the retained-part cap must still report
// the *summed* intensity (the plane field is an amplitude proxy), not the
// strongest single unit's field.
func TestExtendedSourceManyUnitsKeepsTotalIntensity(t *testing.T) {
	cfg := extendedScene(t, map[string]any{
		"width": 2e-3, "height": 0, "samples_x": maxRetainedParts + 4, "samples_y": 1,
		"power": 1e-3,
	})
	res, err := Simulate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	pl := res.Planes[0]
	if len(pl.Parts) != maxRetainedParts {
		t.Errorf("retained %d parts, want the cap %d", len(pl.Parts), maxRetainedParts)
	}
	if !pl.Merged {
		t.Error("a plane combining many incoherent units must be flagged merged")
	}
	tot := pl.TotalIntensity()
	var sum float64
	for _, v := range tot {
		sum += v
	}
	// The proxy field must carry the summed intensity, not one unit's.
	var proxy float64
	for i := range pl.Ex {
		proxy += real(pl.Ex[i])*real(pl.Ex[i]) + imag(pl.Ex[i])*imag(pl.Ex[i])
	}
	if math.Abs(proxy-sum)/sum > 1e-9 {
		t.Errorf("plane field carries %.6g of the summed intensity %.6g", proxy, sum)
	}
	warned := false
	for _, w := range res.Warnings {
		if strings.Contains(w.Code, "parts") {
			warned = true
		}
	}
	if !warned {
		t.Errorf("dropping per-source fields was silent: %v", res.Warnings)
	}
}

// TestExtendedSourceRejectsBadParams covers the parameter validation.
func TestExtendedSourceRejectsBadParams(t *testing.T) {
	for _, p := range []map[string]any{
		{"samples_x": 0},
		{"samples_x": -3},
		{"samples_x": 1.5},
		{"samples_y": 0},
		{"width": -1e-3},
		{"profile": "nope"},
		{"power": 0},
		{"standoff": -1e-3},
	} {
		cfg := extendedScene(t, p)
		if issues := ValidateConfig(&cfg); len(issues) == 0 {
			t.Errorf("params %v passed validation", p)
		}
	}
}

func TestExtendedSourceSampleCentres(t *testing.T) {
	// 4 samples over a 1 cm width: cell centres at ±1.25 mm and ±3.75 mm.
	pts := extendedSampleCentres(1e-2, 4)
	want := []float64{-3.75e-3, -1.25e-3, 1.25e-3, 3.75e-3}
	if len(pts) != len(want) {
		t.Fatalf("got %d centres, want %d", len(pts), len(want))
	}
	for i := range want {
		if math.Abs(pts[i]-want[i]) > 1e-15 {
			t.Errorf("centre %d = %.6g m, want %.6g m", i, pts[i], want[i])
		}
	}
	if got := fmt.Sprint(extendedSampleCentres(0, 3)); got != "[0]" {
		t.Errorf("zero width must give one on-axis sample, got %s", got)
	}
}

// TestExtendedSourceFarFieldMatchesAnalytic is the end-to-end physics check:
// a single emitter point illuminating a narrow slit in the far field must
// reproduce the textbook single-slit Fraunhofer pattern — zeros at
// k*lambda*z/b and a first side lobe at 4.7% of the peak. A 1x1 sampling is
// used so the comparison is against a pure Fraunhofer result.
func TestExtendedSourceFarFieldMatchesAnalytic(t *testing.T) {
	const (
		b  = 2e-4 // slit width (m)
		z  = 1.0  // slit -> detector distance (m)
		sl = 0.45 // emitter -> slit distance (m)
		wl = 632.8e-9
	)
	w := 0.02
	cfg := sceneCfg(1024, w, []SourceSpec{{
		ID: "src", Type: "extended", Pos: &Vec3{Z: -sl}, Dir: &Vec3{Z: 1},
		Params: map[string]any{"power": 1e-3, "samples_x": 1, "samples_y": 1},
	}},
		ComponentSpec{Type: "aperture", Label: "单缝", Pos: v3(0, 0, 0), Shape: RectOutline(b, 0.01)},
		bigDet("屏", v3(0, 0, z), math.Pi, w))
	cfg.Bandlimit = nil
	res, err := Simulate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	pl := res.Planes[0]
	prof, err := pl.ProfileOf("x", KindIntensity, nil)
	if err != nil {
		t.Fatal(err)
	}
	var peak, xPeak float64
	for i, v := range prof.V {
		if v > peak {
			peak, xPeak = v, prof.X[i]
		}
	}
	if math.Abs(xPeak) > 2*pl.DX {
		t.Errorf("pattern peak at x=%g m, want the optical axis", xPeak)
	}
	x1 := wl * z / b
	at := func(x float64) float64 {
		i := int(math.Round(x/pl.DX + float64(pl.Size)/2))
		if i < 0 || i >= pl.Size {
			return 0
		}
		return prof.V[i]
	}
	// Zeros at k*lambda*z/b: the minimum is narrow, so check the pixel next
	// to it as well (a 1024^2 grid over 20 mm has a 19.5 um pixel, the zero
	// is 2.6 mm from its neighbour).
	for _, k := range []float64{1, 2, 3} {
		xe := k * x1
		for _, d := range []float64{-1, 0, 1} {
			if got := at(xe + d*pl.DX); got > 0.06*peak {
				t.Errorf("intensity near the %g-th zero (x=%.4g mm) = %.4g, want < 6%% of the peak %.4g",
					k, (xe+d*pl.DX)*1e3, got, peak)
			}
		}
	}
	// First side lobe at 1.43*x1: the analytic value is 4.7% of I0.
	if lobe := at(1.43*x1) / peak; lobe < 0.02 || lobe > 0.10 {
		t.Errorf("first side lobe = %.4g of the peak, want ~0.047", lobe)
	}
}
