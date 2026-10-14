package lens

import (
	"math"
	"testing"
)

func toleranceAchromat(t *testing.T) *Prescription {
	t.Helper()
	p, _, err := DesignAchromat(100, "N-BK7", "F2", 10, nil)
	if err != nil {
		t.Fatal(err)
	}
	return p
}

// TestTolerancesAchromat runs the full Monte-Carlo analysis on a designed
// achromat with realistic workshop tolerances and checks the report shape.
func TestTolerancesAchromat(t *testing.T) {
	p := toleranceAchromat(t)
	spec := ToleranceSpec{RadiusMM: 0.05, ThicknessMM: 0.05, IndexDelta: 0.0005}
	crit := ToleranceCriteria{
		TargetEFLMM:     100,
		MaxFocalErrMM:   0.1,
		MaxColorShiftMM: 0.05,
		MaxRMSOnAxisMM:  0.02,
		MaxRMSOffAxisMM: 0.03,
	}
	rep, err := RunTolerances(p, spec, crit, &ToleranceOptions{Runs: 400, Seed: 7})
	if err != nil {
		t.Fatal(err)
	}
	if rep.Runs != 400 {
		t.Fatalf("runs = %d", rep.Runs)
	}
	if rep.Yield < 0 || rep.Yield > 1 {
		t.Fatalf("yield out of range: %g", rep.Yield)
	}
	for name, d := range map[string]Distribution{
		"on": rep.OnAxisRMS, "off": rep.OffAxisRMS,
		"color": rep.ColorShift, "focal": rep.FocalErr,
	} {
		if !(d.Min <= d.Mean && d.Mean <= d.Max) {
			t.Fatalf("%s distribution unordered: %+v", name, d)
		}
		if !(d.P50 <= d.P90 && d.P90 <= d.P99 && d.P99 <= d.Max) {
			t.Fatalf("%s percentiles unordered: %+v", name, d)
		}
		if d.StdDev < 0 {
			t.Fatalf("%s negative stddev", name)
		}
	}
	if rep.FailedTraces > 0 && len(rep.FailReasons) == 0 {
		t.Fatalf("failed traces without reasons")
	}
	t.Logf("report:\n%s", rep.String())
}

// TestTolerancesDeterministic pins bit-identical repeatability for a seed.
func TestTolerancesDeterministic(t *testing.T) {
	p := toleranceAchromat(t)
	spec := ToleranceSpec{RadiusMM: 0.05, ThicknessMM: 0.05, IndexDelta: 0.0005}
	crit := ToleranceCriteria{TargetEFLMM: 100, MaxFocalErrMM: 0.1, MaxColorShiftMM: 0.05, MaxRMSOnAxisMM: 0.02}
	a, err := RunTolerances(p, spec, crit, &ToleranceOptions{Runs: 150, Seed: 3})
	if err != nil {
		t.Fatal(err)
	}
	b, err := RunTolerances(p, spec, crit, &ToleranceOptions{Runs: 150, Seed: 3})
	if err != nil {
		t.Fatal(err)
	}
	if a.Yield != b.Yield || a.OnAxisRMS.Mean != b.OnAxisRMS.Mean ||
		a.ColorShift.P99 != b.ColorShift.P99 || a.FailedTraces != b.FailedTraces {
		t.Fatalf("not deterministic:\n%+v\nvs\n%+v", a, b)
	}
}

// TestTolerancesTrend checks the physics of the analysis: looser tolerances
// must widen the distributions and must not improve the yield.
func TestTolerancesTrend(t *testing.T) {
	p := toleranceAchromat(t)
	tight := ToleranceSpec{RadiusMM: 0.01, ThicknessMM: 0.01, IndexDelta: 0.0001}
	loose := ToleranceSpec{RadiusMM: 0.1, ThicknessMM: 0.1, IndexDelta: 0.001}
	crit := ToleranceCriteria{TargetEFLMM: 100, MaxFocalErrMM: 0.15, MaxColorShiftMM: 0.08, MaxRMSOnAxisMM: 0.02}
	rt, err := RunTolerances(p, tight, crit, &ToleranceOptions{Runs: 250, Seed: 11})
	if err != nil {
		t.Fatal(err)
	}
	rl, err := RunTolerances(p, loose, crit, &ToleranceOptions{Runs: 250, Seed: 11})
	if err != nil {
		t.Fatal(err)
	}
	if rl.OnAxisRMS.StdDev <= rt.OnAxisRMS.StdDev {
		t.Fatalf("looser tolerances did not widen the RMS spread: tight σ=%g loose σ=%g",
			rt.OnAxisRMS.StdDev, rl.OnAxisRMS.StdDev)
	}
	if rl.FocalErr.StdDev <= rt.FocalErr.StdDev {
		t.Fatalf("looser tolerances did not widen the focal spread: tight σ=%g loose σ=%g",
			rt.FocalErr.StdDev, rl.FocalErr.StdDev)
	}
	if rl.Yield > rt.Yield+1e-12 {
		t.Fatalf("looser tolerances gave better yield: tight=%g loose=%g", rt.Yield, rl.Yield)
	}
	t.Logf("tight: σ_on=%.5f σ_focal=%.5f yield=%.3f | loose: σ_on=%.5f σ_focal=%.5f yield=%.3f",
		rt.OnAxisRMS.StdDev, rt.FocalErr.StdDev, rt.Yield,
		rl.OnAxisRMS.StdDev, rl.FocalErr.StdDev, rl.Yield)
}

// TestTolerancesValidation checks the input guards.
func TestTolerancesValidation(t *testing.T) {
	p := toleranceAchromat(t)
	if _, err := RunTolerances(p, ToleranceSpec{}, ToleranceCriteria{MaxRMSOnAxisMM: 0.02}, nil); err == nil {
		t.Fatal("empty spec accepted")
	}
	if _, err := RunTolerances(p, ToleranceSpec{RadiusMM: 0.05}, ToleranceCriteria{}, nil); err == nil {
		t.Fatal("empty criteria accepted")
	}
	_ = math.Abs
}
