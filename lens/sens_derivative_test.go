package lens

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
)

// singletAnalytic is the paraxial thick-lens model in air used as the
// analytic reference of the sensitivity tests:
//
//	φ1 = (n−1)/R1,  φ2 = (1−n)/R2,  φ = φ1 + φ2 − (t/n)·φ1·φ2,  f = 1/φ.
//
// For t → 0 it is the thin-lens formula f = 1/((n−1)(1/R1 − 1/R2)).
type singletAnalytic struct{ R1, R2, t, n float64 }

func (s singletAnalytic) phis() (p1, p2 float64) {
	return (s.n - 1) / s.R1, (1 - s.n) / s.R2
}

func (s singletAnalytic) efl() float64 {
	p1, p2 := s.phis()
	return 1 / (p1 + p2 - (s.t/s.n)*p1*p2)
}

// dR1 = ∂f/∂R1 = f²(n−1)(1 − t·φ2/n)/R1².
func (s singletAnalytic) dR1() float64 {
	f := s.efl()
	_, p2 := s.phis()
	return f * f * (s.n - 1) * (1 - s.t*p2/s.n) / (s.R1 * s.R1)
}

// dR2 = ∂f/∂R2 = −f²(n−1)(1 − t·φ1/n)/R2².
func (s singletAnalytic) dR2() float64 {
	f := s.efl()
	p1, _ := s.phis()
	return -f * f * (s.n - 1) * (1 - s.t*p1/s.n) / (s.R2 * s.R2)
}

// dt = ∂f/∂t = f²·φ1·φ2/n.
func (s singletAnalytic) dt() float64 {
	f := s.efl()
	p1, p2 := s.phis()
	return f * f * p1 * p2 / s.n
}

// dn = ∂f/∂n = −f²·∂φ/∂n with
// ∂φ/∂n = 1/R1 − 1/R2 − t·[(φ2/R1 − φ1/R2)/n − φ1φ2/n²].
func (s singletAnalytic) dn() float64 {
	f := s.efl()
	p1, p2 := s.phis()
	dphi := 1/s.R1 - 1/s.R2 - s.t*((p2/s.R1-p1/s.R2)/s.n-p1*p2/(s.n*s.n))
	return -f * f * dphi
}

// singletPrescription builds the singlet of the analytic model (N-BK7 front
// element, air behind, stop on the front surface).
func singletPrescription(R1, R2, t float64) *Prescription {
	return &Prescription{
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{Name: "front", RadiusMM: R1, ThicknessMM: t, Glass: "N-BK7", SemiAperMM: 12},
			{Name: "rear", RadiusMM: R2, Glass: "", SemiAperMM: 12},
		},
	}
}

// solvedAchromat is the analytic achromat solve (no DLS polish), the fixed
// test system of the merit tests.
func solvedAchromat(t *testing.T) *Prescription {
	t.Helper()
	p, _, err := DesignAchromat(100, "N-BK7", "F2", 10, &DesignOptions{SkipPolish: true})
	if err != nil {
		t.Fatal(err)
	}
	return p
}

func relDiff(a, b float64) float64 {
	return math.Abs(a-b) / math.Max(math.Abs(b), 1e-300)
}

// TestSensThinLensEFLAnalytic checks the Richardson sensitivities of the
// paraxial EFL against the closed-form thick-lens derivatives. A thin-lens
// element (t = 0.5 mm) and a thicker one (t = 3 mm) are both covered. The
// acceptance is a relative error ≤ 1e-6 for every derivative with respect to
// curvature, radius, thickness and index.
func TestSensThinLensEFLAnalytic(t *testing.T) {
	g, ok := geo.Glass("N-BK7")
	if !ok {
		t.Fatal("N-BK7 missing")
	}
	n := g.N(geo.Wavelengthd)
	for _, thick := range []float64{0.5, 3} {
		R1, R2 := 80.0, -80.0
		p := singletPrescription(R1, R2, thick)
		an := singletAnalytic{R1: R1, R2: R2, t: thick, n: n}

		// The analytic EFL must match the engine before its derivatives are
		// compared (guards against a wrong reference model).
		eflNum, err := EFLScalar(geo.Wavelengthd)(p)
		if err != nil {
			t.Fatal(err)
		}
		if relDiff(eflNum, an.efl()) > 1e-12 {
			t.Fatalf("t=%g: engine EFL %.15g vs analytic %.15g", thick, eflNum, an.efl())
		}

		// Curvature, thickness and index in one call; radius separately, since
		// radius and curvature of one card are the same storage (see checkParams).
		curv := []Param{
			{Kind: ParamCurvature, Card: 0},
			{Kind: ParamCurvature, Card: 1},
			{Kind: ParamThickness, Card: 0},
			{Kind: ParamIndex, Card: 0},
		}
		wantCurv := []float64{
			-R1 * R1 * an.dR1(), // ∂f/∂c1 = −R1²·∂f/∂R1
			-R2 * R2 * an.dR2(),
			an.dt(),
			an.dn(),
		}
		sc, err := Sensitivity(p, curv, EFLScalar(geo.Wavelengthd), SensOptions{})
		if err != nil {
			t.Fatal(err)
		}
		sr, err := Sensitivity(p, []Param{{Kind: ParamRadius, Card: 0}}, EFLScalar(geo.Wavelengthd), SensOptions{})
		if err != nil {
			t.Fatal(err)
		}
		worst := 0.0
		check := func(name string, num, ana float64, s *SensResult, i int) {
			rel := relDiff(num, ana)
			worst = math.Max(worst, rel)
			t.Logf("t=%.1f %-6s num=%+.12e ana=%+.12e rel=%.2e (RelStep %.0e, uncertainty %.1e, converged %v)",
				thick, name, num, ana, rel, s.RelSteps[i], s.Err[i], s.Converged[i])
			if rel > 1e-6 {
				t.Errorf("t=%g %s: relative error %.3e > 1e-6 (num %g, analytic %g)", thick, name, rel, num, ana)
			}
		}
		for i, pr := range curv {
			check(pr.String(), sc.Deriv[i], wantCurv[i], sc, i)
		}
		check("s0.R", sr.Deriv[0], an.dR1(), sr, 0)
		t.Logf("t=%.1f: worst relative error = %.3e", thick, worst)
	}
}

// TestSensAliasingRejected pins the guard against parameters that write the
// same prescription field, which would otherwise silently overwrite each other.
func TestSensAliasingRejected(t *testing.T) {
	p := singletPrescription(80, -80, 3)
	f := EFLScalar(geo.Wavelengthd)
	for _, pair := range [][]Param{
		{{Kind: ParamRadius, Card: 0}, {Kind: ParamCurvature, Card: 0}},
		{{Kind: ParamCurvature, Card: 0}, {Kind: ParamCurvature, Card: 0}},
	} {
		if _, err := Sensitivity(p, pair, f, SensOptions{}); err == nil {
			t.Fatalf("aliased parameters %v accepted", pair)
		}
	}
	// Distinct fields of one card are fine.
	if _, err := Sensitivity(p, []Param{{Kind: ParamRadius, Card: 0}, {Kind: ParamThickness, Card: 0}, {Kind: ParamIndex, Card: 0}}, f, SensOptions{}); err != nil {
		t.Fatalf("distinct fields rejected: %v", err)
	}
}

// TestSensStepConsistency checks that the sensitivities of the fixed-plane spot
// and wavefront merits do not depend on the step size, and that the adaptive
// step selection returns a value consistent with the plateau.
//
// The step scan is logged for transparency. Measured on the solved achromat
// (fixed image plane, best focus 0.074 mm from the paraxial focus): the spot
// derivatives with respect to curvature and index are unreliable at RelStep
// 1e-3 (truncation, 0.3–18 % in the logged scan) and converged from 1e-4 down; the wavefront
// thickness derivative is roundoff-limited below 1e-5. The pair 1e-4 / 1e-5 is
// inside the converged region for all four parameters of both merits, so the
// acceptance ≤ 1e-4 is applied to that pair and to the adaptive result, and the
// forward difference (O(h)) is checked at the plateau step.
func TestSensStepConsistency(t *testing.T) {
	p := solvedAchromat(t)
	params := []Param{
		{Kind: ParamCurvature, Card: 0},
		{Kind: ParamCurvature, Card: 2},
		{Kind: ParamThickness, Card: 0},
		{Kind: ParamIndex, Card: 0},
	}
	scanSteps := []float64{1e-3, 1e-4, 1e-5, 1e-6, 1e-7}
	for _, kind := range []MeritKind{MeritSpotRMS, MeritWavefrontRMS} {
		m, err := MeritSpec{Kind: kind}.Fixed(p)
		if err != nil {
			t.Fatal(err)
		}
		f := m.Func()
		scan := make([]*SensResult, len(scanSteps))
		for s, rs := range scanSteps {
			scan[s], err = Sensitivity(p, params, f, SensOptions{Fixed: true, RelStep: rs})
			if err != nil {
				t.Fatal(err)
			}
		}
		adapt, err := Sensitivity(p, params, f, SensOptions{})
		if err != nil {
			t.Fatal(err)
		}
		plateauA, plateauB := scan[1], scan[2] // RelStep 1e-4 and 1e-5
		for i := range params {
			t.Logf("%s %-7s scan dM/dx: 1e-3 %+.8e | 1e-4 %+.8e | 1e-5 %+.8e | 1e-6 %+.8e | 1e-7 %+.8e",
				kind, params[i], scan[0].Deriv[i], scan[1].Deriv[i], scan[2].Deriv[i], scan[3].Deriv[i], scan[4].Deriv[i])
			step := relDiff(plateauA.Deriv[i], plateauB.Deriv[i])
			ad := relDiff(adapt.Deriv[i], plateauB.Deriv[i])
			// The forward difference has O(h) truncation error, so it is run at a
			// step ten times smaller than the plateau step (its error is then ~1e-3).
			fwd, err := forwardDiff(p, params[i], f, plateauB.Steps[i]/10)
			if err != nil {
				t.Fatal(err)
			}
			t.Logf("%s %-7s RelStep 1e-4 vs 1e-5: %.2e | adaptive (RelStep %.0e, u %.1e, converged %v) vs 1e-5: %.2e | forward (RelStep 1e-6) vs 1e-5: %.2e",
				kind, params[i], step, adapt.RelSteps[i], adapt.Err[i], adapt.Converged[i], ad,
				relDiff(fwd, plateauB.Deriv[i]))
			if step > 1e-4 {
				t.Errorf("%s %s: step dependence %.3e > 1e-4 inside the plateau", kind, params[i], step)
			}
			if !adapt.Converged[i] || ad > 1e-4 {
				t.Errorf("%s %s: adaptive result %.3e from the plateau value (converged %v)",
					kind, params[i], ad, adapt.Converged[i])
			}
			if relDiff(fwd, plateauB.Deriv[i]) > 1e-2 {
				t.Errorf("%s %s: forward difference disagrees by %.3e > 1%%", kind, params[i], relDiff(fwd, plateauB.Deriv[i]))
			}
		}
	}
}

// forwardDiff is the plain one-sided difference (f(x+h) − f(x))/h of a single
// parameter, used only as an independent cross-check.
func forwardDiff(p *Prescription, pr Param, f ScalarFunc, h float64) (float64, error) {
	f0, err := f(p)
	if err != nil {
		return 0, err
	}
	q, err := applyParams(p, []Param{pr}, []float64{h})
	if err != nil {
		return 0, err
	}
	f1, err := f(q)
	if err != nil {
		return 0, err
	}
	return (f1 - f0) / h, nil
}

// TestSensDeterministic pins bit-identical repeatability of the sensitivity
// engine and the agreement of the extrapolated value with the plain base-step
// central difference.
func TestSensDeterministic(t *testing.T) {
	p := singletPrescription(80, -80, 3)
	params := []Param{{Kind: ParamCurvature, Card: 0}, {Kind: ParamThickness, Card: 0}}
	a, err := Sensitivity(p, params, EFLScalar(geo.Wavelengthd), SensOptions{})
	if err != nil {
		t.Fatal(err)
	}
	b, err := Sensitivity(p, params, EFLScalar(geo.Wavelengthd), SensOptions{})
	if err != nil {
		t.Fatal(err)
	}
	for i := range params {
		if a.Deriv[i] != b.Deriv[i] {
			t.Fatalf("%s: runs differ: %v vs %v", params[i], a.Deriv[i], b.Deriv[i])
		}
		if !a.Converged[i] {
			t.Errorf("%s: step selection did not converge (u %.2e)", params[i], a.Err[i])
		}
		if relDiff(a.Naive[i], a.Deriv[i]) > 1e-3 {
			t.Errorf("%s: naive %.6g far from extrapolated %.6g", params[i], a.Naive[i], a.Deriv[i])
		}
	}
}

// minimizeScalar is a golden-section search for the minimum of f on [a, b].
// The fixed-plane merit is unimodal on the bracket used below.
func minimizeScalar(f func(float64) float64, a, b, tol float64) (float64, float64) {
	const g = 0.6180339887498949
	c := b - g*(b-a)
	d := a + g*(b-a)
	fc, fd := f(c), f(d)
	for b-a > tol {
		if fc < fd {
			b, d, fd = d, c, fc
			c = b - g*(b-a)
			fc = f(c)
		} else {
			a, c, fc = c, d, fd
			d = a + g*(b-a)
			fd = f(d)
		}
	}
	x := 0.5 * (a + b)
	return x, f(x)
}

// bestFocus returns the image-plane offset (mm, from the paraxial focus) that
// minimises the fixed-plane merit: a coarse scan over ±0.3 mm followed by a
// golden-section refinement. It also returns the merit at the minimum.
func bestFocus(t *testing.T, m MeritSpec, p *Prescription) (float64, float64) {
	t.Helper()
	at := func(dz float64) float64 {
		mm := m
		mm.ImageZ = m.ImageZ + dz
		v, err := mm.Value(p)
		if err != nil {
			t.Fatal(err)
		}
		return v
	}
	best, bestZ := math.Inf(1), 0.0
	for dz := -0.3; dz <= 0.3+1e-12; dz += 0.01 {
		if v := at(dz); v < best {
			best, bestZ = v, dz
		}
	}
	return minimizeScalar(at, bestZ-0.01, bestZ+0.01, 1e-8)
}

// TestMeritFocusMinimum checks the best-focus physics of the fixed-plane merit.
// The minimum of the spot RMS over the image plane is the best focus, which
// spherical aberration moves away from the paraxial focus; the paraxial
// (floating) merit is the value at the paraxial focus. For small pupils the
// third-order theory predicts that the best-focus shift scales as r² (the
// pupil radius), so doubling the pupil radius should multiply the shift by 4.
func TestMeritFocusMinimum(t *testing.T) {
	p := solvedAchromat(t)
	axis := []Field{{0, 0}}

	// Full aperture, on-axis: the best focus is displaced from the paraxial focus.
	m, err := MeritSpec{Kind: MeritSpotRMS, Fields: axis}.Fixed(p)
	if err != nil {
		t.Fatal(err)
	}
	zStar, mStar := bestFocus(t, m, p)
	m0, err := m.Value(p)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("full aperture on-axis: best focus dz* = %+.6f mm, M* = %.6e mm; paraxial focus M(0) = %.6e mm",
		zStar, mStar, m0)
	if math.Abs(zStar) < 0.01 {
		t.Errorf("best focus %+.6f mm not displaced from the paraxial focus by ≥ 0.01 mm", zStar)
	}
	if !(mStar < m0/2) {
		t.Errorf("best-focus merit %.6g not well below the paraxial-focus merit %.6g", mStar, m0)
	}
	// Floating focus = fixed merit at the paraxial focus (same plane, same rays).
	ff := m
	ff.FloatingFocus = true
	vf, err := ff.Value(p)
	if err != nil {
		t.Fatal(err)
	}
	if relDiff(vf, m0) > 1e-9 {
		t.Errorf("floating-focus merit %.12g differs from the paraxial-focus value %.12g", vf, m0)
	}

	// Third-order scaling of the best-focus shift with pupil radius.
	shift := func(r float64) float64 {
		s := []Sample{{0, 0}, {r, 0}, {0, r}, {-r, 0}, {0, -r}}
		mr, err := MeritSpec{Kind: MeritSpotRMS, Fields: axis, Samples: s}.Fixed(p)
		if err != nil {
			t.Fatal(err)
		}
		z, _ := bestFocus(t, mr, p)
		return z
	}
	z05, z10 := shift(0.05), shift(0.10)
	ratio := z10 / z05
	t.Logf("small pupils on-axis: dz*(0.05) = %+.6e mm, dz*(0.10) = %+.6e mm, ratio = %.4f (r² law: 4)",
		z05, z10, ratio)
	if math.Abs(ratio-4) > 0.08 {
		t.Errorf("best-focus shift ratio %.4f not within 2%% of the third-order value 4", ratio)
	}
	if math.Abs(z05) > 1e-3 {
		t.Errorf("small-pupil best-focus shift %.3e mm exceeds 1e-3 mm", z05)
	}
}

// TestMeritLengthAndUnsetPlane checks the residual layout and the guard
// against an unset image plane.
func TestMeritLengthAndUnsetPlane(t *testing.T) {
	p := singletPrescription(80, -80, 3)
	spot := MeritSpec{Kind: MeritSpotRMS}
	if spot.Len() != 2*2*len(UnitDiskSamples(3, 6)) {
		t.Fatalf("spot residual length %d", spot.Len())
	}
	wf := MeritSpec{Kind: MeritWavefrontRMS}
	if wf.Len() != 2*len(UnitDiskSamples(3, 6)) {
		t.Fatalf("wavefront residual length %d", wf.Len())
	}
	if _, err := spot.Value(p); err == nil {
		t.Fatal("unset image plane accepted")
	}
}
