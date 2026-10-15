package lens

import (
	"math"
	"testing"
)

// tolerancePerturbation is the fixed tolerance draw used by the compensator
// tests: a curvature change of 2e-5 1/mm (radius change ≈ 0.05 mm at R ≈ 50 mm),
// a thickness change of 0.03 mm and a glass index change of 3e-4.
func tolerancePerturbation() Perturbation {
	return Perturbation{
		Params: []Param{
			{Kind: ParamCurvature, Card: 0},
			{Kind: ParamCurvature, Card: 1},
			{Kind: ParamThickness, Card: 0},
			{Kind: ParamIndex, Card: 0},
		},
		Offsets: []float64{2e-5, -3e-5, 0.03, 3e-4},
	}
}

// TestCompensateImageShiftMatchesBestFocus checks the compensator against an
// independent reference. With the image plane as the only compensator, the
// compensated merit must equal the best-focus merit of the perturbed system,
// found by a coarse scan and a golden-section search (not by the DLS).
func TestCompensateImageShiftMatchesBestFocus(t *testing.T) {
	p := solvedAchromat(t)
	m, err := MeritSpec{Kind: MeritSpotRMS}.Fixed(p)
	if err != nil {
		t.Fatal(err)
	}
	pert := tolerancePerturbation()
	res, err := Compensate(p, pert, m, CompensatorSpec{ImageShift: true, ShiftRange: 1})
	if err != nil {
		t.Fatal(err)
	}
	pq, err := applyParams(p, pert.Params, pert.Offsets)
	if err != nil {
		t.Fatal(err)
	}
	zStar, mStar := bestFocus(t, m, pq)
	t.Logf("perturbed merit %.6e mm, nominal %.6e mm; compensated (image shift %+.6f mm) %.6e mm; best focus %+.6f mm, merit %.6e mm; recovery %.4f",
		res.Perturbed, res.Nominal, res.ImageShift, res.Compensated, zStar, mStar, res.Recovery)
	if res.Compensated > res.Perturbed {
		t.Errorf("compensated merit %.9g exceeds perturbed merit %.9g", res.Compensated, res.Perturbed)
	}
	if math.Abs(res.ImageShift-zStar) > 1e-5 {
		t.Errorf("compensator image shift %+.6f mm differs from the best focus %+.6f mm", res.ImageShift, zStar)
	}
	if relDiff(res.Compensated, mStar) > 1e-6 {
		t.Errorf("compensated merit %.9g differs from the best-focus merit %.9g", res.Compensated, mStar)
	}
	if !(res.Perturbed > res.Nominal) {
		t.Errorf("the tolerance draw did not degrade the merit (%.6g vs %.6g)", res.Perturbed, res.Nominal)
	}
}

// TestCompensateRecovery runs a compensator set of an image shift and the front
// radius on a DLS-polished achromat, after the fixed tolerance draw. The
// acceptance is Compensated ≤ Perturbed (and CompensatedNominal ≤ Nominal); the
// recovery share is reported and its observed value is frozen as a regression
// threshold.
func TestCompensateRecovery(t *testing.T) {
	p, _, err := DesignAchromat(100, "N-BK7", "F2", 10, nil)
	if err != nil {
		t.Fatal(err)
	}
	m, err := MeritSpec{Kind: MeritSpotRMS}.Fixed(p)
	if err != nil {
		t.Fatal(err)
	}
	res, err := Compensate(p, tolerancePerturbation(), m, CompensatorSpec{
		Vars:       []VarSpec{{Card: 0, Kind: VarRadius}},
		ImageShift: true,
	})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("nominal %.6e mm (compensated %.6e), perturbed %.6e mm, compensated %.6e mm, image shift %+.6f mm, s0.R -> %.6f mm, recovery %.4f (DLS iterations %d, converged %v)",
		res.Nominal, res.CompensatedNominal, res.Perturbed, res.Compensated, res.ImageShift,
		res.Prescription.Cards[0].RadiusMM, res.Recovery, res.Opt.Iter, res.Opt.Converged)
	if res.Compensated > res.Perturbed {
		t.Errorf("compensated merit %.9g exceeds perturbed merit %.9g", res.Compensated, res.Perturbed)
	}
	if res.CompensatedNominal > res.Nominal {
		t.Errorf("compensating the nominal worsened it: %.9g > %.9g", res.CompensatedNominal, res.Nominal)
	}
	if !(res.Recovery > 0) {
		t.Errorf("no recovery (%.4f)", res.Recovery)
	}
	// Frozen regression thresholds (observed in the log above: merit reduction
	// 0.820, recovery 1.0006; thresholds set below the observed values).
	if red := (res.Perturbed - res.Compensated) / res.Perturbed; red < 0.75 {
		t.Errorf("merit reduction after compensation %.4f < frozen threshold 0.75", red)
	}
	if res.Recovery < 0.9 {
		t.Errorf("recovery %.4f < frozen threshold 0.9", res.Recovery)
	}
}

// TestCompensateRejectsBadSpec pins the input guards: no compensator, a floating
// focus with an image shift, and a range list of the wrong length.
func TestCompensateRejectsBadSpec(t *testing.T) {
	p := solvedAchromat(t)
	m, err := MeritSpec{Kind: MeritSpotRMS}.Fixed(p)
	if err != nil {
		t.Fatal(err)
	}
	pert := tolerancePerturbation()
	if _, err := Compensate(p, pert, m, CompensatorSpec{}); err == nil {
		t.Error("no compensator accepted")
	}
	ff := m
	ff.FloatingFocus = true
	if _, err := Compensate(p, pert, ff, CompensatorSpec{ImageShift: true}); err == nil {
		t.Error("image shift with a floating focus accepted")
	}
	if _, err := Compensate(p, pert, m, CompensatorSpec{Vars: []VarSpec{{Card: 2, Kind: VarRadius}}, VarRange: []float64{1, 2}}); err == nil {
		t.Error("range list of the wrong length accepted")
	}
}

// TestCompensateAirSpacing uses an air spacing as the compensator: two
// biconvex N-BK7 singlets with a 20 mm air gap (card 1 is the air gap, its
// thickness is the spacing). A tolerance draw on the four radii moves the
// focus; the spacing alone is re-optimised. The acceptance is the same
// monotonicity as before, plus a non-trivial spacing change.
func TestCompensateAirSpacing(t *testing.T) {
	p := &Prescription{
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{RadiusMM: 60, ThicknessMM: 4, Glass: "N-BK7", SemiAperMM: 4},
			{RadiusMM: -60, ThicknessMM: 20, Glass: "", SemiAperMM: 4},
			{RadiusMM: 60, ThicknessMM: 4, Glass: "N-BK7", SemiAperMM: 4},
			{RadiusMM: -60, Glass: "", SemiAperMM: 4},
		},
	}
	m, err := MeritSpec{Kind: MeritSpotRMS}.Fixed(p)
	if err != nil {
		t.Fatal(err)
	}
	pert := Perturbation{
		Params:  []Param{{Kind: ParamCurvature, Card: 0}, {Kind: ParamCurvature, Card: 2}},
		Offsets: []float64{3e-5, -3e-5},
	}
	res, err := Compensate(p, pert, m, CompensatorSpec{Vars: []VarSpec{{Card: 1, Kind: VarThickness}}})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("air spacing: nominal %.6e mm (compensated %.6e), perturbed %.6e mm, compensated %.6e mm, spacing %.4f -> %.4f mm, recovery %.4f (converged %v)",
		res.Nominal, res.CompensatedNominal, res.Perturbed, res.Compensated,
		p.Cards[1].ThicknessMM, res.Prescription.Cards[1].ThicknessMM, res.Recovery, res.Opt.Converged)
	if res.Compensated > res.Perturbed {
		t.Errorf("compensated merit %.9g exceeds perturbed merit %.9g", res.Compensated, res.Perturbed)
	}
	if res.Prescription.Cards[1].ThicknessMM == p.Cards[1].ThicknessMM {
		t.Error("the air spacing did not move")
	}
	if res.Prescription.Cards[1].Glass != "" {
		t.Error("the spacing compensator altered the glass of card 1")
	}
}
