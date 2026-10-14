package coating

import (
	"math"
	"strings"
	"testing"
)

// This file holds the tests of the optimiser core: the damped least squares
// solver, the merit function and its target scales, the option bookkeeping and
// the layer-growth (GrowLayers) path. The physical guarantees the solver must
// preserve - energy conservation and the log-scale merit for low-reflectance
// goals - are checked here; the spectral promises of the templates are checked
// in templates_test.go.

// optimiseFixture returns a two-layer starting design and a single
// low-reflectance target that the optimiser has to work on.
func optimiseFixture() (Stack, []Target) {
	st := NewStack(Air, BK7,
		Layer{Material: MgF2, ThicknessNm: 40},
		Layer{Material: Ta2O5, ThicknessNm: 80},
	)
	return st, []Target{RTarget(550, 1e-4, 1)}
}

// TestOptimizeImprovesAndReports checks the report of a successful fit: the
// merit improves, the history starts at the initial merit and never rises, the
// RMS residual is sqrt(merit/m), Improvement() exceeds one, and the reported
// merit is the merit of the returned design.
func TestOptimizeImprovesAndReports(t *testing.T) {
	st, targets := optimiseFixture()
	res, err := Optimize(st, targets, Options{})
	if err != nil {
		t.Fatalf("Optimize: %v", err)
	}
	if res.Merit >= res.InitialMerit {
		t.Errorf("merit did not improve: initial %g, final %g", res.InitialMerit, res.Merit)
	}
	if len(res.History) < 2 {
		t.Fatalf("history has %d entries, want the initial merit plus at least one accepted step", len(res.History))
	}
	if res.History[0] != res.InitialMerit {
		t.Errorf("history[0] = %g, want the initial merit %g", res.History[0], res.InitialMerit)
	}
	for i := 1; i < len(res.History); i++ {
		if res.History[i] > res.History[i-1] {
			t.Errorf("history is not non-increasing at step %d: %g after %g", i, res.History[i], res.History[i-1])
			break
		}
	}
	if last := res.History[len(res.History)-1]; last != res.Merit {
		t.Errorf("history ends at %g, want the final merit %g", last, res.Merit)
	}
	closeFloat(t, "RmsResidual", res.RmsResidual, math.Sqrt(res.Merit/float64(len(targets))), 1e-12)
	if res.Improvement() <= 1 {
		t.Errorf("Improvement() = %g, want > 1", res.Improvement())
	}
	m, err := Merit(res.Stack, res.Targets)
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "reported merit", m, res.Merit, 1e-12)
	values, err := res.Achieved()
	if err != nil {
		t.Fatal(err)
	}
	if len(values) != len(targets) {
		t.Fatalf("Achieved returned %d values for %d targets", len(values), len(targets))
	}
	got, err := res.Stack.Reflectance(targets[0].LambdaNm, targets[0].AngleDeg, targets[0].Pol)
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "achieved reflectance", values[0], got, 1e-12)
	if !strings.Contains(res.Summary(), "MgF2") || !strings.Contains(res.Summary(), "merit") {
		t.Errorf("Summary does not describe the design:\n%s", res.Summary())
	}
}

// TestOptimizeBoundsAndFreeMask checks the box constraints and the freeze mask,
// and that the caller's stack is not mutated: Optimize works on its own copy.
func TestOptimizeBoundsAndFreeMask(t *testing.T) {
	st, targets := optimiseFixture()
	opts := Options{MinThicknessNm: 20, MaxThicknessNm: 120}
	res, err := Optimize(st, targets, opts)
	if err != nil {
		t.Fatalf("Optimize: %v", err)
	}
	for i, l := range res.Stack.Layers {
		if l.ThicknessNm < opts.MinThicknessNm || l.ThicknessNm > opts.MaxThicknessNm {
			t.Errorf("layer %d thickness %.4f nm outside the bounds [%g, %g]", i, l.ThicknessNm, opts.MinThicknessNm, opts.MaxThicknessNm)
		}
	}
	if st.Layers[0].ThicknessNm != 40 || st.Layers[1].ThicknessNm != 80 {
		t.Errorf("the caller's stack was mutated: %v", st.ThicknessesNm())
	}

	// Freezing the first layer must keep it bit-exactly and still move the
	// second one (the target is far from the starting design).
	opts.Free = []bool{false, true}
	res2, err := Optimize(st, targets, opts)
	if err != nil {
		t.Fatalf("Optimize with mask: %v", err)
	}
	if res2.Stack.Layers[0].ThicknessNm != 40 {
		t.Errorf("frozen layer moved from 40 nm to %g nm", res2.Stack.Layers[0].ThicknessNm)
	}
	if res2.Stack.Layers[1].ThicknessNm == 80 {
		t.Errorf("free layer stayed at its starting thickness %g nm", res2.Stack.Layers[1].ThicknessNm)
	}
	// A short mask leaves the layers it does not mention free (layer 1 is free
	// when only layer 0 is masked).
	opts.Free = []bool{false}
	res3, err := Optimize(st, targets, opts)
	if err != nil {
		t.Fatalf("Optimize with short mask: %v", err)
	}
	if res3.Stack.Layers[0].ThicknessNm != 40 {
		t.Errorf("masked layer moved with a short mask: %g nm", res3.Stack.Layers[0].ThicknessNm)
	}
	// Annihilation is allowed when the floor is explicitly zero: a redundant
	// layer inserted into an already-good design must be driven to zero. The
	// outer layers are frozen so that the only way to fit the reference design
	// is to remove the middle layer.
	ref := NewStack(Air, BK7,
		Layer{Material: SiO2, ThicknessNm: 120},
		Layer{Material: Ta2O5, ThicknessNm: 90},
	)
	var refTargets []Target
	for _, lambda := range []float64{450, 550, 650} {
		r, err := ref.Reflectance(lambda, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		refTargets = append(refTargets, RTarget(lambda, r, 1).LinearScale())
	}
	redundant := NewStack(Air, BK7,
		Layer{Material: SiO2, ThicknessNm: 120},
		Layer{Material: MgF2, ThicknessNm: 50},
		Layer{Material: Ta2O5, ThicknessNm: 90},
	)
	resR, err := Optimize(redundant, refTargets, Options{MinThicknessNm: 0, Free: []bool{false, true, false}, MaxIter: 100})
	if err != nil {
		t.Fatalf("Optimize with MinThicknessNm 0: %v", err)
	}
	if d := resR.Stack.Layers[1].ThicknessNm; d > 1e-6 {
		t.Errorf("redundant layer thickness %g nm, want it annihilated to zero", d)
	}
	refMerit, err := Merit(ref, refTargets)
	if err != nil {
		t.Fatal(err)
	}
	if resR.Merit > refMerit+1e-15 {
		t.Errorf("annihilated design merit %g, want the reference merit %g", resR.Merit, refMerit)
	}
	// With the template floor in place the same layer may not vanish.
	resF, err := Optimize(redundant, refTargets, Options{MinThicknessNm: 1, Free: []bool{false, true, false}, MaxIter: 100})
	if err != nil {
		t.Fatal(err)
	}
	if d := resF.Stack.Layers[1].ThicknessNm; d < 1 {
		t.Errorf("floored layer thickness %g nm, want >= 1 nm", d)
	}
}

// TestOptimizeDeterministicAndBudget checks that the search is reproducible for
// a given seed and that the iteration budget is honoured.
func TestOptimizeDeterministicAndBudget(t *testing.T) {
	st, targets := optimiseFixture()
	a, err := Optimize(st, targets, Options{Restarts: 2})
	if err != nil {
		t.Fatal(err)
	}
	b, err := Optimize(st, targets, Options{Restarts: 2})
	if err != nil {
		t.Fatal(err)
	}
	if a.Merit != b.Merit || a.Iterations != b.Iterations || a.Starts != b.Starts {
		t.Errorf("runs differ: merit %g/%g, iterations %d/%d, starts %d/%d", a.Merit, b.Merit, a.Iterations, b.Iterations, a.Starts, b.Starts)
	}
	da, db := a.Stack.ThicknessesNm(), b.Stack.ThicknessesNm()
	for i := range da {
		if da[i] != db[i] {
			t.Errorf("thickness %d differs between identical runs: %g vs %g", i, da[i], db[i])
		}
	}
	if a.Starts != 3 {
		t.Errorf("Starts = %d, want 1 + 2 restarts = 3", a.Starts)
	}

	budget := 5
	c, err := Optimize(st, targets, Options{MaxIter: budget, Restarts: 1})
	if err != nil {
		t.Fatal(err)
	}
	if c.Iterations > budget*c.Starts || c.Iterations < 1 {
		t.Errorf("Iterations = %d, want between 1 and MaxIter·Starts = %d", c.Iterations, budget*c.Starts)
	}
}

// TestOptimizeErrorPaths checks the input validation of the optimiser.
func TestOptimizeErrorPaths(t *testing.T) {
	st, targets := optimiseFixture()
	if _, err := Optimize(NewStack(Air, BK7), targets, Options{}); err == nil {
		t.Error("Optimize accepted a stack without layers")
	}
	if _, err := Optimize(st, nil, Options{}); err == nil {
		t.Error("Optimize accepted an empty target list")
	}
	if _, err := Optimize(st, []Target{{LambdaNm: 550, Pol: PolUnpolarized, Qty: Quantity(9), Value: 0.5, Weight: 1}}, Options{}); err == nil {
		t.Error("Optimize accepted an invalid target quantity")
	}
	bad := st.Clone()
	bad.Layers[1].ThicknessNm = -5
	if _, err := Optimize(bad, targets, Options{}); err == nil {
		t.Error("Optimize accepted a negative thickness")
	}
	if _, err := Optimize(st, targets, Options{Free: []bool{false, false}}); err == nil {
		t.Error("Optimize accepted an all-frozen mask")
	}
	if _, err := Merit(st, []Target{{LambdaNm: 550, Qty: Quantity(3)}}); err == nil {
		t.Error("Merit accepted an invalid quantity")
	}
	if _, err := TargetValues(NewStack(Air, BK7), []Target{{LambdaNm: 550, Qty: Quantity(3)}}); err == nil {
		t.Error("TargetValues accepted an invalid quantity")
	}
	// GrowLayers must refuse a material that carries no dispersion model.
	if _, err := GrowLayers(st, targets, []Material{{Name: "broken"}}, Options{GrowthRounds: 1}); err == nil {
		t.Error("GrowLayers accepted an invalid growth material")
	}
}

// TestMeritScales checks the two residual scales of the merit function.
//
//   - A logarithmic target measures the relative error: changing the target
//     value by two decades changes the residual by exactly 2, and a vanishing
//     computed quantity stays finite thanks to the log floor.
//   - A linear target measures the absolute error.
//   - The weight enters through its square root, so the merit is the weighted
//     sum of squared residuals.
func TestMeritScales(t *testing.T) {
	bare := NewStack(Air, BK7)
	const lambda = 550.0
	R, err := bare.Reflectance(lambda, 0, PolUnpolarized)
	if err != nil {
		t.Fatal(err)
	}
	if R <= 0 || R >= 1 {
		t.Fatalf("bare BK7 reflectance %g is not a usable reference", R)
	}

	logTarget := RTarget(lambda, 1e-6, 1)
	if !logTarget.Log {
		t.Fatal("RTarget(1e-6) must select the logarithmic scale")
	}
	mLog, err := Merit(bare, []Target{logTarget})
	if err != nil {
		t.Fatal(err)
	}
	residual := math.Log10(R) - math.Log10(1e-6)
	closeFloat(t, "log merit", mLog, residual*residual, 1e-12)

	// Two decades of target value are two units of residual.
	mLog2, err := Merit(bare, []Target{RTarget(lambda, 1e-4, 1)})
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "decade spacing", math.Sqrt(mLog)-math.Sqrt(mLog2), 2, 1e-12)

	// Linear scale: absolute error.
	mLin, err := Merit(bare, []Target{logTarget.LinearScale()})
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "linear merit", mLin, (R-1e-6)*(R-1e-6), 1e-12)

	// The weight enters as its square root: four times the weight, four times
	// the merit.
	heavy := logTarget
	heavy.Weight = 4
	mHeavy, err := Merit(bare, []Target{heavy})
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "weighted merit", mHeavy, 4*mLog, 1e-12)

	// A zero target value is legal on the log scale: the floor keeps the
	// residual finite instead of producing -Inf.
	zero := Target{LambdaNm: lambda, Pol: PolUnpolarized, Qty: Reflectance, Value: 0, Weight: 1, Log: true}
	mZero, err := Merit(bare, []Target{zero})
	if err != nil {
		t.Fatal(err)
	}
	if math.IsInf(mZero, 0) || math.IsNaN(mZero) {
		t.Fatalf("log merit of a zero target = %g, want a finite floor value", mZero)
	}
	floored := math.Log10(R) - math.Log10(logFloor)
	closeFloat(t, "log floor merit", mZero, floored*floored, 1e-12)

	// A high-reflectance goal is linear: the log scale is flat as R -> 1.
	hr := RTarget(lambda, 0.9999, 1)
	if hr.Log {
		t.Error("RTarget(0.9999) must stay on the linear scale")
	}
}

// TestTargetHelpers checks the small conveniences of Target: the logarithmic
// scale selection of RTarget/TTarget, the fluent At/LogScale/LinearScale
// copies, and the report string.
func TestTargetHelpers(t *testing.T) {
	if r := RTarget(550, 0.01, 2); !r.Log || r.Qty != Reflectance || r.Pol != PolUnpolarized || r.AngleDeg != 0 {
		t.Errorf("RTarget(0.01) = %+v, want a logarithmic normal-incidence reflectance target", r)
	}
	if r := RTarget(550, 0.2, 2); r.Log {
		t.Errorf("RTarget(0.2) = %+v, want the linear scale", r)
	}
	if tr := TTarget(550, 1e-3, 1); !tr.Log || tr.Qty != Transmittance {
		t.Errorf("TTarget(1e-3) = %+v, want a logarithmic transmittance target", tr)
	}
	if a := ATarget(550, 0.01, 1); a.Qty != Absorptance {
		t.Errorf("ATarget = %+v, want the absorptance quantity", a)
	}
	angled := RTarget(550, 0.2, 1).At(45, PolP)
	if angled.AngleDeg != 45 || angled.Pol != PolP {
		t.Errorf("At did not copy the operating point: %+v", angled)
	}
	if s := angled.String(); !strings.Contains(s, "45") || !strings.Contains(s, "p") {
		t.Errorf("String = %q, want the angle and polarization", s)
	}
	if !angled.LogScale().Log || angled.LinearScale().Log {
		t.Error("LogScale/LinearScale did not switch the scale")
	}
	if len(HRCandidates()) == 0 || len(ARCandidates()) == 0 {
		t.Error("candidate lists must not be empty")
	}
}

// TestBandTargets checks the band-target constructor: the wavelength grid is
// inclusive, every target carries the requested operating point, and the usual
// range errors surface.
func TestBandTargets(t *testing.T) {
	ts, err := BandTargets(400, 500, 25, 30, PolS, Transmittance, 0.99, 2, false)
	if err != nil {
		t.Fatal(err)
	}
	if len(ts) != 5 {
		t.Fatalf("got %d targets, want 5 (400, 425, ..., 500)", len(ts))
	}
	for i, tg := range ts {
		want := 400 + 25*float64(i)
		if tg.LambdaNm != want {
			t.Errorf("target %d wavelength %g, want %g", i, tg.LambdaNm, want)
		}
		if tg.AngleDeg != 30 || tg.Pol != PolS || tg.Qty != Transmittance || tg.Value != 0.99 || tg.Weight != 2 || tg.Log {
			t.Errorf("target %d = %+v, want the requested fields", i, tg)
		}
	}
	if _, err := BandTargets(500, 400, 10, 0, PolUnpolarized, Reflectance, 0.01, 1, true); err == nil {
		t.Error("BandTargets accepted a descending range")
	}
	if _, err := BandTargets(400, 500, 0, 0, PolUnpolarized, Reflectance, 0.01, 1, true); err == nil {
		t.Error("BandTargets accepted a zero step")
	}
}

// TestOptionsDefaults documents the mapping of unset option fields. The zero
// value is usable; notably MinThicknessNm == 0 is not "unset" but really means
// "layers may be annihilated" (DefaultOptions selects the 1 nm floor instead).
func TestOptionsDefaults(t *testing.T) {
	d := DefaultOptions()
	if d.MaxIter != 200 || d.Tolerance != 1e-9 || d.Restarts != 0 || d.Seed != 1 ||
		d.RestartJitter != 0.6 || d.MinThicknessNm != 1 || d.MaxThicknessNm != 2000 ||
		d.GrowthRounds != 0 || d.GrowthThicknessNm != 10 {
		t.Errorf("DefaultOptions = %+v", d)
	}
	z := Options{}.withDefaults()
	if z.MaxIter != 200 || z.Tolerance != 1e-9 || z.Restarts != 0 || z.Seed != 1 ||
		z.RestartJitter != 0.6 || z.MaxThicknessNm != 2000 || z.GrowthThicknessNm != 10 {
		t.Errorf("zero Options filled as %+v", z)
	}
	if z.MinThicknessNm != 0 {
		t.Errorf("zero Options: MinThicknessNm = %g, want 0 (annihilation allowed)", z.MinThicknessNm)
	}
	if got := (Options{MinThicknessNm: -1}).withDefaults().MinThicknessNm; got != 1 {
		t.Errorf("negative MinThicknessNm resolved to %g, want the 1 nm default", got)
	}
	if got := (Options{MinThicknessNm: 50, MaxThicknessNm: 10}).withDefaults().MaxThicknessNm; got != 50 {
		t.Errorf("MaxThicknessNm below MinThicknessNm resolved to %g, want it raised to the minimum 50", got)
	}
}

// TestGrowLayersDesign exercises the layer-growth path end to end: a two-layer
// stack that is far from a broadband antireflection design is grown with four
// needle-like rounds. The measured outcome is a four-layer design whose merit
// improves by a factor 2.1 (107.3 -> 51.0) and whose mean reflectance over
// 420-680 nm drops to 0.89% - as good as the six-layer BroadbandAR template on
// the same band, which is the point of the growth path.
func TestGrowLayersDesign(t *testing.T) {
	t.Parallel()
	targets, err := BandTargets(420, 680, 20, 0, PolUnpolarized, Reflectance, 1e-4, 1, true)
	if err != nil {
		t.Fatal(err)
	}
	start := NewStack(Air, BK7,
		Layer{Material: MgF2, ThicknessNm: 60},
		Layer{Material: Ta2O5, ThicknessNm: 80},
	)
	startMerit, err := Merit(start, targets)
	if err != nil {
		t.Fatal(err)
	}
	res, err := GrowLayers(start, targets, nil, Options{GrowthRounds: 4, Restarts: 2})
	if err != nil {
		t.Fatalf("GrowLayers: %v", err)
	}
	if len(res.Stack.Layers) <= len(start.Layers) {
		t.Errorf("growth did not add a layer: %d layers at the start, %d after", len(start.Layers), len(res.Stack.Layers))
	}
	if res.Merit >= startMerit {
		t.Errorf("merit did not improve: start %g, grown %g", startMerit, res.Merit)
	}
	if improvement := startMerit / res.Merit; improvement < 1.5 {
		t.Errorf("merit improved only %.2fx (measured 2.1x), want at least 1.5x", improvement)
	}
	mean, _, maxR := spectrumStats(t, res.Stack, 420, 680, 2, Reflectance)
	if mean > 0.015 {
		t.Errorf("grown design: mean R over 420-680 nm = %.5f, want < 1.5%% (measured 0.89%%)", mean)
	}
	if maxR > 0.05 {
		t.Errorf("grown design: peak R = %.5f, want < 5%% (measured 1.86%%)", maxR)
	}
	for i, l := range res.Stack.Layers {
		if !l.Material.Valid() {
			t.Errorf("layer %d has no material", i)
		}
		if l.ThicknessNm < 1 || l.ThicknessNm > 2000 {
			t.Errorf("layer %d thickness %.4f nm outside the optimiser bounds", i, l.ThicknessNm)
		}
	}
	// The growth report must account for the work done.
	if res.Starts < 1 || res.Iterations < 1 {
		t.Errorf("report does not describe the search: %d starts, %d iterations", res.Starts, res.Iterations)
	}
	if m, err := Merit(res.Stack, res.Targets); err != nil {
		t.Fatal(err)
	} else {
		closeFloat(t, "reported merit", m, res.Merit, 1e-12)
	}
	// Growth never degrades a design: continuing from the grown stack can only
	// keep or improve the merit (every round that would end worse is rejected).
	again, err := GrowLayers(res.Stack, targets, nil, Options{GrowthRounds: 2, Restarts: 0})
	if err != nil {
		t.Fatalf("GrowLayers from the grown stack: %v", err)
	}
	if again.Merit > res.Merit {
		t.Errorf("a further growth run worsened the merit: %g after %g", again.Merit, res.Merit)
	}
}

// TestGrowLayersZeroRounds checks that a growth run without rounds is exactly
// an Optimize call: the same start and options must give the same design.
func TestGrowLayersZeroRounds(t *testing.T) {
	st, targets := optimiseFixture()
	opts := Options{Restarts: 1}
	plain, err := Optimize(st, targets, opts)
	if err != nil {
		t.Fatal(err)
	}
	grown, err := GrowLayers(st, targets, nil, opts)
	if err != nil {
		t.Fatal(err)
	}
	if grown.Merit != plain.Merit || grown.Starts != plain.Starts || grown.Iterations != plain.Iterations {
		t.Errorf("GrowLayers without rounds differs from Optimize: merit %g/%g, starts %d/%d, iterations %d/%d",
			grown.Merit, plain.Merit, grown.Starts, plain.Starts, grown.Iterations, plain.Iterations)
	}
	ga, pa := grown.Stack.ThicknessesNm(), plain.Stack.ThicknessesNm()
	for i := range ga {
		if ga[i] != pa[i] {
			t.Errorf("thickness %d differs: %g vs %g", i, ga[i], pa[i])
		}
	}
}
