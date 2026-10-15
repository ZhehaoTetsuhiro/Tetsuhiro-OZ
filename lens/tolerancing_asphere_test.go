package lens

import (
	"math"
	"math/rand"
	"testing"
)

// ---------------------------------------------------------------------------
// CoefAbs tolerancing of polynomial aspheres and Forbes Q aspheres
// ---------------------------------------------------------------------------

// TestTolerancingCoefAbsStreamLocked locks the random-draw order of one
// Monte-Carlo run, which is the part of the tolerance engine that must not
// change: per card the radius (finite non-zero), the thickness (all but the
// last card), the conic constant, then the aspheric coefficients in stored
// order — Coeffs (A4, A6, …), then PolyTerms, then QAsphere.Coeffs — and
// finally the single glass-melt Δn of the run. The expected values below are
// replayed here with the same generator and the documented order, so any
// reshuffle (for instance drawing the poly terms before the even-asphere
// coefficients, or perturbing the coefficients of a card in a different
// order) fails this test even though the report would still look plausible.
func TestTolerancingCoefAbsStreamLocked(t *testing.T) {
	p := &Prescription{
		Name:      "stream lock",
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{Name: "even asphere", RadiusMM: 100, ConicK: -0.5, ThicknessMM: 6, Glass: "N-BK7",
				Coeffs: []float64{-3e-7, 1e-10}},
			{Name: "poly asphere", RadiusMM: -100, ThicknessMM: 40, Glass: "",
				PolyTerms: []PolyTermCard{{Power: 4, Coeff: -2e-7}, {Power: 6, Coeff: 2e-11}}},
			{Name: "q asphere", RadiusMM: math.Inf(1), ConicK: 0, ThicknessMM: 20, Glass: "F2",
				QAsphere: &QCard{Form: "qcon", NormRadiusMM: 10, Coeffs: []float64{-0.01, 0.002, -0.0005}}},
		},
		Stop: 0,
	}
	spec := ToleranceSpec{RadiusMM: 0.05, ThicknessMM: 0.03, ConicK: 0.01, CoefAbs: 0.02, IndexDelta: 0.0005}
	const seed = 42
	q, err := perturbPrescription(p, spec, rand.New(rand.NewSource(seed)))
	if err != nil {
		t.Fatal(err)
	}
	// Replay the documented draw order by hand into a second copy.
	want := p.Clone()
	rng := rand.New(rand.NewSource(seed))
	u := func(scale float64) float64 { return (rng.Float64()*2 - 1) * scale }
	for i := range want.Cards {
		c := &want.Cards[i]
		src := p.Cards[i]
		if spec.RadiusMM > 0 && src.RadiusMM != 0 && !math.IsInf(src.RadiusMM, 0) {
			c.RadiusMM += u(spec.RadiusMM)
		}
		if spec.ThicknessMM > 0 && i < len(want.Cards)-1 {
			c.ThicknessMM += u(spec.ThicknessMM)
		}
		if spec.ConicK > 0 {
			c.ConicK += u(spec.ConicK)
		}
		if spec.CoefAbs > 0 {
			for j := range c.Coeffs {
				c.Coeffs[j] += u(spec.CoefAbs)
			}
			for j := range c.PolyTerms {
				c.PolyTerms[j].Coeff += u(spec.CoefAbs)
			}
			if c.QAsphere != nil {
				for j := range c.QAsphere.Coeffs {
					c.QAsphere.Coeffs[j] += u(spec.CoefAbs)
				}
			}
		}
	}
	// The melt Δn is drawn once, after every card.
	_ = u(spec.IndexDelta)
	for i := range q.Cards {
		a, b := q.Cards[i], want.Cards[i]
		if a.RadiusMM != b.RadiusMM || a.ThicknessMM != b.ThicknessMM || a.ConicK != b.ConicK {
			t.Fatalf("card %d geometry mismatch: R %g/%g T %g/%g K %g/%g",
				i, a.RadiusMM, b.RadiusMM, a.ThicknessMM, b.ThicknessMM, a.ConicK, b.ConicK)
		}
		if len(a.Coeffs) != len(b.Coeffs) || len(a.PolyTerms) != len(b.PolyTerms) {
			t.Fatalf("card %d coefficient list lengths changed", i)
		}
		for j := range a.Coeffs {
			if a.Coeffs[j] != b.Coeffs[j] {
				t.Fatalf("card %d Coeffs[%d] = %g, want %g", i, j, a.Coeffs[j], b.Coeffs[j])
			}
		}
		for j := range a.PolyTerms {
			if a.PolyTerms[j] != b.PolyTerms[j] {
				t.Fatalf("card %d PolyTerms[%d] = %+v, want %+v", i, j, a.PolyTerms[j], b.PolyTerms[j])
			}
		}
		switch {
		case a.QAsphere == nil && b.QAsphere == nil:
		case a.QAsphere == nil || b.QAsphere == nil:
			t.Fatalf("card %d q_asphere presence changed", i)
		default:
			if a.QAsphere.Form != b.QAsphere.Form || a.QAsphere.NormRadiusMM != b.QAsphere.NormRadiusMM {
				t.Fatalf("card %d q_asphere header changed: %+v vs %+v", i, a.QAsphere, b.QAsphere)
			}
			if len(a.QAsphere.Coeffs) != len(b.QAsphere.Coeffs) {
				t.Fatalf("card %d q_asphere coefficient list length changed", i)
			}
			for j := range b.QAsphere.Coeffs {
				if a.QAsphere.Coeffs[j] != b.QAsphere.Coeffs[j] {
					t.Fatalf("card %d QAsphere.Coeffs[%d] = %g, want %g", i, j, a.QAsphere.Coeffs[j], b.QAsphere.Coeffs[j])
				}
			}
		}
	}
	// A plain even-asphere prescription (the pre-existing case) must draw its
	// coefficients exactly as it always did: radius, thickness, conic, Coeffs,
	// per card in order.
	old := &Prescription{
		Name: "legacy", ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{Name: "a", RadiusMM: 120, ThicknessMM: 5, Glass: "N-BK7", Coeffs: []float64{-1e-7, 2e-11}},
			{Name: "b", RadiusMM: -80, ThicknessMM: 30, Glass: ""},
		},
	}
	legacySpec := ToleranceSpec{RadiusMM: 0.02, ThicknessMM: 0.01, CoefAbs: 0.05}
	q2, err := perturbPrescription(old, legacySpec, rand.New(rand.NewSource(7)))
	if err != nil {
		t.Fatal(err)
	}
	rng2 := rand.New(rand.NewSource(7))
	u2 := func(scale float64) float64 { return (rng2.Float64()*2 - 1) * scale }
	_ = u2(legacySpec.RadiusMM)    // card 0 radius
	_ = u2(legacySpec.ThicknessMM) // card 0 thickness
	wantA4 := old.Cards[0].Coeffs[0] + u2(legacySpec.CoefAbs)
	wantA6 := old.Cards[0].Coeffs[1] + u2(legacySpec.CoefAbs)
	_ = u2(legacySpec.RadiusMM) // card 1 radius
	if q2.Cards[0].Coeffs[0] != wantA4 || q2.Cards[0].Coeffs[1] != wantA6 {
		t.Fatalf("legacy even-asphere draws changed: got %v, want [%g %g]", q2.Cards[0].Coeffs, wantA4, wantA6)
	}
}

// TestTolerancingCoefAbsPolyAsphere is the acceptance test of the CoefAbs
// extension: a single lens corrected by a polynomial-asphere r⁴ term is
// toleranced with and without a coefficient drift. With the coefficient held
// (CoefAbs = 0) every run of the radius tolerance passes the 20 µm RMS
// criteria; adding the coefficient drift visibly widens the RMS distribution
// and drops the yield, and a repeated run is bit-identical.
func TestTolerancingCoefAbsPolyAsphere(t *testing.T) {
	const a4 = -8.4e-7
	crit := ToleranceCriteria{TargetEFLMM: 100, MaxRMSOnAxisMM: 0.02, MaxRMSOffAxisMM: 0.02}
	opts := &ToleranceOptions{Runs: 200, Seed: 5, Fields: asphereFields, Samples: asphereSamples}
	without, err := RunTolerances(polySinglet(a4), ToleranceSpec{RadiusMM: 0.05}, crit, opts)
	if err != nil {
		t.Fatal(err)
	}
	spec := ToleranceSpec{RadiusMM: 0.05, CoefAbs: 1e-7}
	with, err := RunTolerances(polySinglet(a4), spec, crit, opts)
	if err != nil {
		t.Fatal(err)
	}
	if without.Yield != 1 {
		t.Fatalf("the radius-only reference must pass every run, got yield %g", without.Yield)
	}
	if !(with.Yield < without.Yield) {
		t.Fatalf("the coefficient tolerance did not change the yield: %g vs %g", with.Yield, without.Yield)
	}
	if with.Yield <= 0 {
		t.Fatalf("the coefficient tolerance wiped out the yield: %g", with.Yield)
	}
	if with.OnAxisRMS.StdDev < 5*without.OnAxisRMS.StdDev {
		t.Fatalf("the coefficient tolerance did not widen the RMS spread: σ = %g vs %g",
			with.OnAxisRMS.StdDev, without.OnAxisRMS.StdDev)
	}
	if with.FailedTraces != 0 {
		t.Fatalf("perturbed prescriptions failed to trace: %d", with.FailedTraces)
	}
	// The perturbation lands on the polynomial coefficient, not somewhere
	// else: the on-axis spot of a perturbed design is degraded by an r⁴
	// wavefront term, so the RMS distribution grows monotonically with
	// CoefAbs.
	bigger, err := RunTolerances(polySinglet(a4), ToleranceSpec{RadiusMM: 0.05, CoefAbs: 2e-7}, crit, opts)
	if err != nil {
		t.Fatal(err)
	}
	if !(bigger.OnAxisRMS.Mean > with.OnAxisRMS.Mean) {
		t.Fatalf("a larger CoefAbs did not degrade the mean RMS: %g vs %g", bigger.OnAxisRMS.Mean, with.OnAxisRMS.Mean)
	}
	// Deterministic repeats.
	again, err := RunTolerances(polySinglet(a4), spec, crit, opts)
	if err != nil {
		t.Fatal(err)
	}
	if again.Yield != with.Yield || again.OnAxisRMS.Mean != with.OnAxisRMS.Mean ||
		again.OnAxisRMS.StdDev != with.OnAxisRMS.StdDev || again.OnAxisRMS.Max != with.OnAxisRMS.Max ||
		again.OffAxisRMS.Mean != with.OffAxisRMS.Mean || again.Nominal.OnAxisRMSMM != with.Nominal.OnAxisRMSMM {
		t.Fatalf("repeated runs differ:\n%+v\nvs\n%+v", with, again)
	}
}

// TestTolerancingCoefAbsQAsphere runs the same exercise for a Forbes qcon
// surface: the coefficient drift of the q_asphere coefficients must change
// the yield, and the report must be reproducible. The equivalent coefficient
// amplitude is s₀ = a₄·ρ⁴, so a drift of 1e-3 is ≈ 6 % of the −0.0172 that
// corrects the lens.
func TestTolerancingCoefAbsQAsphere(t *testing.T) {
	s0 := -8.4e-7 * math.Pow(12, 4)
	crit := ToleranceCriteria{TargetEFLMM: 100, MaxRMSOnAxisMM: 0.02, MaxRMSOffAxisMM: 0.02}
	opts := &ToleranceOptions{Runs: 200, Seed: 5, Fields: asphereFields, Samples: asphereSamples}
	without, err := RunTolerances(qSinglet(s0, 0), ToleranceSpec{RadiusMM: 0.05}, crit, opts)
	if err != nil {
		t.Fatal(err)
	}
	spec := ToleranceSpec{RadiusMM: 0.05, CoefAbs: 1e-3}
	with, err := RunTolerances(qSinglet(s0, 0), spec, crit, opts)
	if err != nil {
		t.Fatal(err)
	}
	if without.Yield != 1 {
		t.Fatalf("the radius-only reference must pass every run, got yield %g", without.Yield)
	}
	if !(with.Yield < without.Yield) || with.Yield <= 0 {
		t.Fatalf("the q coefficient tolerance did not change the yield properly: %g vs %g", with.Yield, without.Yield)
	}
	if with.OnAxisRMS.StdDev < 5*without.OnAxisRMS.StdDev {
		t.Fatalf("the q coefficient tolerance did not widen the RMS spread: σ = %g vs %g",
			with.OnAxisRMS.StdDev, without.OnAxisRMS.StdDev)
	}
	again, err := RunTolerances(qSinglet(s0, 0), spec, crit, opts)
	if err != nil {
		t.Fatal(err)
	}
	if again.Yield != with.Yield || again.OnAxisRMS.Mean != with.OnAxisRMS.Mean ||
		again.OnAxisRMS.P99 != with.OnAxisRMS.P99 || again.FailReasons["spot_on_axis"] != with.FailReasons["spot_on_axis"] {
		t.Fatalf("repeated runs differ:\n%+v\nvs\n%+v", with, again)
	}
}
