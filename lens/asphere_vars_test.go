package lens

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
)

// ---------------------------------------------------------------------------
// design variables for polynomial aspheres and Forbes Q aspheres
// ---------------------------------------------------------------------------

// polySinglet is the reference single-lens prescription of the asphere
// variable tests: a BK7 singlet whose front surface carries a polynomial
// asphere term (r⁴) and whose stop is the front surface. With a₄ = 0 the
// surface is the plain sphere of the card radius; a₄ ≈ −8.4e-7 mm⁻³ nulls the
// third-order spherical aberration of the F/3.8 lens.
func polySinglet(a4 float64) *Prescription {
	return &Prescription{
		Name:      "asphere singlet (poly)",
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{Name: "front asphere", RadiusMM: 100, ThicknessMM: 6, Glass: "N-BK7", SemiAperMM: 12.7,
				PolyTerms: []PolyTermCard{{Power: 4, Coeff: a4}}},
			{Name: "rear", RadiusMM: -100, ThicknessMM: 95, Glass: "", SemiAperMM: 12.7},
		},
		Stop: 0,
	}
}

// qSinglet is the same lens with the front surface described as a Forbes
// qcon asphere. The r⁴ coefficient of the departure is s₀/ρ⁴, so the
// coefficient that nulls the lens's spherical aberration is s₀ ≈ −0.0172 for
// ρ = 12 mm.
func qSinglet(s0, s1 float64) *Prescription {
	return &Prescription{
		Name:      "asphere singlet (qcon)",
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{Name: "front asphere", RadiusMM: 100, ThicknessMM: 6, Glass: "N-BK7", SemiAperMM: 12.7,
				QAsphere: &QCard{Form: "qcon", NormRadiusMM: 12, Coeffs: []float64{s0, s1}}},
			{Name: "rear", RadiusMM: -100, ThicknessMM: 95, Glass: "", SemiAperMM: 12.7},
		},
		Stop: 0,
	}
}

// asphereFields and asphereSamples are the on/off-axis merit sample set of the
// asphere tests (the same as the tolerance defaults: on axis and 0.5°).
var (
	asphereFields  = []Field{{0, 0}, {0, 0.0087}}
	asphereSamples = UnitDiskSamples(2, 6)
)

// spotMerit is the merit function of the asphere tests: the spot residuals of
// the two fields against a 10 µm target.
func spotMerit(q *Prescription) ([]float64, error) {
	return SpotResiduals(q, geo.Wavelengthd, asphereFields, asphereSamples, 0.01)
}

// meritOf sums the squared residuals (the Problem's S(x) convention).
func meritOf(t *testing.T, q *Prescription) float64 {
	t.Helper()
	res, err := spotMerit(q)
	if err != nil {
		t.Fatal(err)
	}
	s := 0.0
	for _, r := range res {
		s += r * r
	}
	return s
}

// TestVarSpecAsphereKinds pins the plumbing of the two new variable kinds:
// labels, reads/writes by index and the out-of-range errors. The index
// convention is fixed (stored order of the coefficient lists).
func TestVarSpecAsphereKinds(t *testing.T) {
	p := polySinglet(-3e-7)
	p.Cards[0].PolyTerms = []PolyTermCard{{Power: 4, Coeff: -3e-7}, {Power: 6, Coeff: 1e-10}}
	specs := []VarSpec{
		{Card: 0, Kind: VarPolyTerm, Order: 0},
		{Card: 0, Kind: VarPolyTerm, Order: 1},
	}
	if got := specs[0].String(); got != "s0.poly0" {
		t.Errorf("VarPolyTerm label = %q, want s0.poly0", got)
	}
	x, err := p.Vars(specs)
	if err != nil {
		t.Fatal(err)
	}
	if x[0] != -3e-7 || x[1] != 1e-10 {
		t.Fatalf("Vars = %v, want [-3e-07 1e-10]", x)
	}
	if err := p.SetVars(specs, []float64{-4e-7, 2e-10}); err != nil {
		t.Fatal(err)
	}
	if p.Cards[0].PolyTerms[0].Coeff != -4e-7 || p.Cards[0].PolyTerms[1].Coeff != 2e-10 {
		t.Fatalf("SetVars did not land: %+v", p.Cards[0].PolyTerms)
	}
	// The power is untouched: the index addresses the term, not the monomial.
	if p.Cards[0].PolyTerms[1].Power != 6 {
		t.Fatalf("the term power changed: %d", p.Cards[0].PolyTerms[1].Power)
	}
	if _, err := p.Vars([]VarSpec{{Card: 0, Kind: VarPolyTerm, Order: 2}}); err == nil {
		t.Error("out-of-range polynomial term accepted")
	}
	if _, err := p.Vars([]VarSpec{{Card: 1, Kind: VarPolyTerm, Order: 0}}); err == nil {
		t.Error("polynomial term of a card without PolyTerms accepted")
	}

	q := qSinglet(-0.017, 0)
	qSpecs := []VarSpec{{Card: 0, Kind: VarQCoef, Order: 0}, {Card: 0, Kind: VarQCoef, Order: 1}}
	if got := qSpecs[1].String(); got != "s0.q1" {
		t.Errorf("VarQCoef label = %q, want s0.q1", got)
	}
	if v, err := q.GetVar(qSpecs[0]); err != nil || v != -0.017 {
		t.Fatalf("GetVar q0 = %v, %v", v, err)
	}
	if err := q.SetVar(qSpecs[1], 0.001); err != nil {
		t.Fatal(err)
	}
	if q.Cards[0].QAsphere.Coeffs[1] != 0.001 {
		t.Fatalf("SetVar q1 did not land: %v", q.Cards[0].QAsphere.Coeffs)
	}
	if _, err := q.Vars([]VarSpec{{Card: 0, Kind: VarQCoef, Order: 2}}); err == nil {
		t.Error("out-of-range q coefficient accepted")
	}
	if _, err := p.Vars([]VarSpec{{Card: 0, Kind: VarQCoef, Order: 0}}); err == nil {
		t.Error("q coefficient of a card without QAsphere accepted")
	}
	// Clone detaches the coefficient storage (the optimizer clones per run).
	c := q.Clone()
	c.Cards[0].QAsphere.Coeffs[0] = 1
	c.Cards[0].PolyTerms = append(c.Cards[0].PolyTerms, PolyTermCard{Power: 2, Coeff: 1})
	if q.Cards[0].QAsphere.Coeffs[0] == 1 || len(q.Cards[0].PolyTerms) != 0 {
		t.Error("Clone shares the coefficient storage")
	}
}

// TestOptimizePolyAsphereVariable optimizes the r⁴ coefficient of the
// front-surface polynomial asphere with the DLS solver: the merit must fall
// monotonically, the run must converge, the spherical aberration must be
// corrected (the start is a plain sphere, the end is the aspheric singlet),
// and a repeated run must be bit-for-bit identical.
func TestOptimizePolyAsphereVariable(t *testing.T) {
	p := polySinglet(0)
	specs := []VarSpec{{Card: 0, Kind: VarPolyTerm, Order: 0}}
	bounds := func(s VarSpec, v float64) (float64, float64) {
		if s.Kind == VarPolyTerm {
			return -1e-5, 1e-5
		}
		return v - 1, v + 1
	}
	run := func() *Result {
		t.Helper()
		prob, err := PrescriptionProblem(p, specs, bounds, spotMerit)
		if err != nil {
			t.Fatal(err)
		}
		if prob.Variables[0].Name != "s0.poly0" {
			t.Fatalf("variable label = %q", prob.Variables[0].Name)
		}
		res, err := Minimize(prob)
		if err != nil {
			t.Fatal(err)
		}
		return res
	}
	start := meritOf(t, p)
	a := run()
	if !a.Converged {
		t.Fatalf("not converged after %d iterations (merit %g)", a.Iter, a.Merit)
	}
	// Merit monotone: every accepted iteration of Minimize strictly reduces S.
	for i := 1; i < len(a.History); i++ {
		if !(a.History[i] < a.History[i-1]) {
			t.Fatalf("merit not decreasing at iteration %d: %g -> %g", i, a.History[i-1], a.History[i])
		}
	}
	if a.History[0] != start {
		t.Fatalf("history starts at %g, want the initial merit %g", a.History[0], start)
	}
	if a.Merit > start/100 {
		t.Fatalf("the aspheric term did not correct the lens: merit %g -> %g", start, a.Merit)
	}
	// The optimum corrects the third-order spherical aberration of the F/3.8
	// singlet (≈ −8.4e-7 mm⁻³); the coefficient must have moved there.
	if x := a.X[0]; x > -5e-7 || x < -1.2e-6 {
		t.Fatalf("r⁴ optimum %g is outside the SA-correcting range [-1.2e-6, -5e-7]", x)
	}
	rms0, _, _, err := SpotRMS(p, geo.Wavelengthd, asphereFields[0], asphereSamples)
	if err != nil {
		t.Fatal(err)
	}
	opt := p.Clone()
	if err := opt.SetVars(specs, a.X); err != nil {
		t.Fatal(err)
	}
	rms1, _, _, err := SpotRMS(opt, geo.Wavelengthd, asphereFields[0], asphereSamples)
	if err != nil {
		t.Fatal(err)
	}
	if !(rms1 < 0.2*rms0) {
		t.Fatalf("on-axis RMS %g -> %g, want a large reduction", rms0, rms1)
	}
	// The caller's prescription is untouched (the residual works on clones).
	if p.Cards[0].PolyTerms[0].Coeff != 0 {
		t.Fatalf("the optimization mutated the input prescription: %g", p.Cards[0].PolyTerms[0].Coeff)
	}
	// Repeat: bit-for-bit identical.
	b := run()
	if a.X[0] != b.X[0] || a.Merit != b.Merit || a.RMS != b.RMS || a.Iter != b.Iter ||
		a.Evals != b.Evals || a.Converged != b.Converged || len(a.History) != len(b.History) {
		t.Fatalf("runs differ:\n%+v\nvs\n%+v", a, b)
	}
	for i := range a.History {
		if a.History[i] != b.History[i] {
			t.Fatalf("history[%d] differs: %g vs %g", i, a.History[i], b.History[i])
		}
	}
}

// TestOptimizeQAsphereVariable repeats the DLS exercise with the Forbes qcon
// coefficients s₀ and s₁ of the front surface as the design variables: two
// coupled asphere coefficients must converge to the same corrected lens, with
// a monotone merit and bit-identical repeats.
func TestOptimizeQAsphereVariable(t *testing.T) {
	p := qSinglet(0, 0)
	specs := []VarSpec{{Card: 0, Kind: VarQCoef, Order: 0}, {Card: 0, Kind: VarQCoef, Order: 1}}
	bounds := func(s VarSpec, v float64) (float64, float64) {
		if s.Order == 0 {
			return -0.2, 0.2
		}
		return -0.05, 0.05
	}
	run := func() *Result {
		prob, err := PrescriptionProblem(p, specs, bounds, spotMerit)
		if err != nil {
			t.Fatal(err)
		}
		if prob.Variables[0].Name != "s0.q0" || prob.Variables[1].Name != "s0.q1" {
			t.Fatalf("variable labels = %q, %q", prob.Variables[0].Name, prob.Variables[1].Name)
		}
		res, err := Minimize(prob)
		if err != nil {
			t.Fatal(err)
		}
		return res
	}
	start := meritOf(t, p)
	a := run()
	if !a.Converged {
		t.Fatalf("not converged after %d iterations (merit %g)", a.Iter, a.Merit)
	}
	for i := 1; i < len(a.History); i++ {
		if !(a.History[i] < a.History[i-1]) {
			t.Fatalf("merit not decreasing at iteration %d: %g -> %g", i, a.History[i-1], a.History[i])
		}
	}
	if a.Merit > start/100 {
		t.Fatalf("the q departure did not correct the lens: merit %g -> %g", start, a.Merit)
	}
	// s₀ ≈ −0.0172 corrects the spherical aberration; s₁ adjusts the
	// higher-order residual and must stay small against it.
	if a.X[0] > -0.012 || a.X[0] < -0.03 {
		t.Fatalf("qcon s₀ optimum %g is outside the SA-correcting range", a.X[0])
	}
	if math.Abs(a.X[1]) > 0.01 {
		t.Fatalf("qcon s₁ optimum %g is implausibly large", a.X[1])
	}
	b := run()
	if a.X[0] != b.X[0] || a.X[1] != b.X[1] || a.Merit != b.Merit || a.RMS != b.RMS ||
		a.Iter != b.Iter || a.Evals != b.Evals || a.Converged != b.Converged {
		t.Fatalf("runs differ:\n%+v\nvs\n%+v", a, b)
	}
	for i := range a.History {
		if a.History[i] != b.History[i] {
			t.Fatalf("history[%d] differs: %g vs %g", i, a.History[i], b.History[i])
		}
	}
	// The converged qcon prescription must at least match the single-r⁴ poly
	// asphere (its basis with s₁ = 0 reproduces that term exactly, and s₁ only
	// adds freedom), and its on-axis spot must sit at the merit target.
	best := p.Clone()
	if err := best.SetVars(specs, a.X); err != nil {
		t.Fatal(err)
	}
	rmsQ, _, _, err := SpotRMS(best, geo.Wavelengthd, asphereFields[0], asphereSamples)
	if err != nil {
		t.Fatal(err)
	}
	polyOpt := polySinglet(-8.42e-7)
	rmsP, _, _, err := SpotRMS(polyOpt, geo.Wavelengthd, asphereFields[0], asphereSamples)
	if err != nil {
		t.Fatal(err)
	}
	if rmsQ > rmsP {
		t.Fatalf("qcon optimum RMS %g is worse than the equivalent single-r⁴ poly asphere %g", rmsQ, rmsP)
	}
	if rmsQ > 0.01 {
		t.Fatalf("qcon optimum on-axis RMS %g exceeds the 10 µm merit target", rmsQ)
	}
}
