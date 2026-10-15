package lens

import (
	"math"
	"testing"
)

// TestAllocateEqualContributionLinear checks the inverse-sensitivity allocation
// in the worst-case (linear) budget. Acceptance: the linear-model increment
// predicted from the allocated tolerances, Σ|s_i|·t_i, equals the budget to
// ≤ 1e-9 (relative). Non-capped parameters contribute equal shares |s_i|·t_i.
func TestAllocateEqualContributionLinear(t *testing.T) {
	params := []Param{
		{Kind: ParamCurvature, Card: 0},
		{Kind: ParamCurvature, Card: 1},
		{Kind: ParamThickness, Card: 0},
		{Kind: ParamIndex, Card: 0},
		{Kind: ParamCurvature, Card: 2},
	}
	sens := []float64{2.0, -0.5, 0.1, 0, -4.0}
	const budget = 0.01
	inf := math.Inf(1)

	// No binding limit: four sensitive parameters, zero-sensitivity one held at its cap.
	a, err := AllocateTolerances(params, sens, budget, AllocLinear, []float64{inf, inf, inf, 0.02, inf})
	if err != nil {
		t.Fatal(err)
	}
	checkAllocation(t, a, sens, budget, 1e-9)
	share := budget / 4
	for i := range params {
		if i == 3 {
			if !a.Capped[i] || a.Tol[i] != 0.02 || a.Share[i] != 0 {
				t.Errorf("zero-sensitivity parameter: tol %g capped %v share %g", a.Tol[i], a.Capped[i], a.Share[i])
			}
			continue
		}
		if a.Capped[i] {
			t.Errorf("%s capped without a binding limit", params[i])
		}
		if relDiff(math.Abs(sens[i])*a.Tol[i], share) > 1e-12 {
			t.Errorf("%s: share %.12g, want equal share %.12g", params[i], math.Abs(sens[i])*a.Tol[i], share)
		}
	}
	t.Logf("linear, no binding limit: tol = %v, predicted %.15g vs budget %.15g (rel %.2e)",
		a.Tol, a.Predicted, a.Budget, a.RelErr)

	// A binding limit on the strongest sensitivity: the rest share the remaining budget.
	a2, err := AllocateTolerances(params, sens, budget, AllocLinear, []float64{1e-5, inf, inf, 0.02, inf})
	if err != nil {
		t.Fatal(err)
	}
	checkAllocation(t, a2, sens, budget, 1e-9)
	if !a2.Capped[0] || a2.Tol[0] != 1e-5 {
		t.Fatalf("binding limit not applied: tol %g capped %v", a2.Tol[0], a2.Capped[0])
	}
	rest := (budget - math.Abs(sens[0])*1e-5) / 3
	for _, i := range []int{1, 2, 4} {
		if a2.Capped[i] {
			t.Errorf("%s capped unexpectedly", params[i])
		}
		if relDiff(math.Abs(sens[i])*a2.Tol[i], rest) > 1e-12 {
			t.Errorf("%s: share %.12g, want water-filled share %.12g", params[i], math.Abs(sens[i])*a2.Tol[i], rest)
		}
	}
	t.Logf("linear, binding limit on s0.c: tol = %v, predicted rel err %.2e", a2.Tol, a2.RelErr)
}

// TestAllocateEqualContributionRSS checks the statistical (root-sum-square)
// budget: sqrt(Σ(s_i·t_i)²) equals the budget to ≤ 1e-9, with equal squared
// contributions among the non-capped parameters.
func TestAllocateEqualContributionRSS(t *testing.T) {
	params := []Param{
		{Kind: ParamCurvature, Card: 0},
		{Kind: ParamThickness, Card: 0},
		{Kind: ParamIndex, Card: 0},
	}
	sens := []float64{3.0, -0.2, 1.5}
	const budget = 0.004
	a, err := AllocateTolerances(params, sens, budget, AllocRSS, nil)
	if err != nil {
		t.Fatal(err)
	}
	checkAllocation(t, a, sens, budget, 1e-9)
	for i := range params {
		if relDiff((sens[i]*a.Tol[i])*(sens[i]*a.Tol[i]), budget*budget/3) > 1e-12 {
			t.Errorf("%s: squared share %.12g, want %.12g", params[i], (sens[i]*a.Tol[i])*(sens[i]*a.Tol[i]), budget*budget/3)
		}
	}
	t.Logf("rss: tol = %v, predicted %.15g vs budget %.15g (rel %.2e)", a.Tol, a.Predicted, a.Budget, a.RelErr)
}

// checkAllocation asserts the acceptance of the linear model: the predicted
// increment recomputed from the allocated tolerances equals the budget.
func checkAllocation(t *testing.T, a *Allocation, sens []float64, budget, tol float64) {
	t.Helper()
	var lin, sq float64
	for i := range sens {
		c := math.Abs(sens[i]) * a.Tol[i]
		lin += c
		sq += c * c
	}
	want := lin
	if a.Mode == AllocRSS {
		want = math.Sqrt(sq)
	}
	if relDiff(want, budget) > tol {
		t.Errorf("%v: predicted increment %.15g vs budget %.15g: relative difference %.3e > %.0e",
			a.Mode, want, budget, relDiff(want, budget), tol)
	}
	if a.RelErr > tol || relDiff(a.Predicted, budget) > tol {
		t.Errorf("%v: reported prediction %.15g (RelErr %.3e) off the budget %.15g", a.Mode, a.Predicted, a.RelErr, budget)
	}
	for i := range sens {
		if !(a.Tol[i] > 0) {
			t.Errorf("parameter %d: non-positive tolerance %g", i, a.Tol[i])
		}
	}
}

// TestAllocateErrors pins the rejected inputs: a zero-sensitivity parameter
// without a limit, a budget that cannot be placed, and inconsistent lengths.
func TestAllocateErrors(t *testing.T) {
	params := []Param{{Kind: ParamCurvature, Card: 0}, {Kind: ParamIndex, Card: 0}}
	inf := math.Inf(1)
	if _, err := AllocateTolerances(params, []float64{1, 0}, 0.01, AllocLinear, nil); err == nil {
		t.Error("zero sensitivity without a limit accepted")
	}
	if _, err := AllocateTolerances(params, []float64{1, 1}, 0.01, AllocLinear, []float64{1e-6, 1e-6}); err == nil {
		t.Error("budget that cannot be placed accepted")
	}
	if _, err := AllocateTolerances(params, []float64{1}, 0.01, AllocLinear, nil); err == nil {
		t.Error("length mismatch accepted")
	}
	if _, err := AllocateTolerances(params, []float64{1, 1}, -1, AllocRSS, []float64{inf, inf}); err == nil {
		t.Error("negative budget accepted")
	}
}

// TestAllocateRealLens drives the inverse sensitivity from the Richardson
// sensitivities of the fixed-plane spot merit on the solved achromat. The
// acceptance (predicted increment within 1e-9 of the budget) applies to both
// budget modes. The nonlinear response to the allocated tolerances, applied at
// the worst-case corner (each parameter moved by sign(∂M/∂x)·t), is logged as a
// diagnostic of the linear model's validity; it is not part of the acceptance.
func TestAllocateRealLens(t *testing.T) {
	p := solvedAchromat(t)
	m, err := MeritSpec{Kind: MeritSpotRMS}.Fixed(p)
	if err != nil {
		t.Fatal(err)
	}
	params := []Param{
		{Kind: ParamCurvature, Card: 0},
		{Kind: ParamCurvature, Card: 2},
		{Kind: ParamThickness, Card: 0},
		{Kind: ParamIndex, Card: 0},
	}
	sr, err := Sensitivity(p, params, m.Func(), SensOptions{})
	if err != nil {
		t.Fatal(err)
	}
	m0, err := m.Value(p)
	if err != nil {
		t.Fatal(err)
	}
	for _, budget := range []float64{2e-4, 2e-5} {
		a, err := AllocateTolerances(params, sr.Deriv, budget, AllocLinear, nil)
		if err != nil {
			t.Fatal(err)
		}
		checkAllocation(t, a, sr.Deriv, budget, 1e-9)
		x := make([]float64, len(params))
		for i := range params {
			x[i] = math.Copysign(a.Tol[i], sr.Deriv[i])
		}
		q, err := applyParams(p, params, x)
		if err != nil {
			t.Fatal(err)
		}
		m1, err := m.Value(q)
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("linear budget %.1e mm: tol %v; worst-corner ΔM actual %.4e mm (ratio to budget %.4f, nominal M0 %.4e)",
			budget, a.Tol, m1-m0, (m1-m0)/budget, m0)
		aR, err := AllocateTolerances(params, sr.Deriv, budget, AllocRSS, nil)
		if err != nil {
			t.Fatal(err)
		}
		checkAllocation(t, aR, sr.Deriv, budget, 1e-9)
		t.Logf("rss budget %.1e mm: tol %v, predicted rel err %.2e", budget, aR.Tol, aR.RelErr)
	}
}
