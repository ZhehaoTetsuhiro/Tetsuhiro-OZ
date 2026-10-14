package coating

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
)

// TestSingleLayerARZeroReflectance checks the analytic optimum of a single
// quarter-wave layer: with n1 = sqrt(n0·ns) the amplitude reflection at the
// design wavelength vanishes exactly, so R has to be at round-off level.
func TestSingleLayerARZeroReflectance(t *testing.T) {
	const lambda = 550.0
	for _, tc := range []struct {
		name      string
		n0, ns    float64
	}{
		{"air/BK7", 1.0, 1.5168},
		{"air/fused silica", 1.0, 1.4585},
		{"air/SF11", 1.0, 1.7847},
		{"cemented glass", 1.52, 1.62},
	} {
		t.Run(tc.name, func(t *testing.T) {
			inc := NewMaterial("incident", ConstIndex{N: tc.n0})
			sub := NewMaterial("substrate", ConstIndex{N: tc.ns})
			st, err := SingleLayerARExact(inc, sub, lambda)
			if err != nil {
				t.Fatalf("SingleLayerARExact: %v", err)
			}
			smp, err := st.Evaluate(lambda, 0, PolUnpolarized)
			if err != nil {
				t.Fatalf("Evaluate: %v", err)
			}
			if smp.R > 1e-24 {
				t.Errorf("R at the design wavelength = %g, want < 1e-24", smp.R)
			}
			closeFloat(t, "T", smp.T, 1, 1e-12)
			// Away from the design wavelength the layer stops being a quarter
			// wave and the reflectance rises: the single-layer AR is narrow.
			far, err := st.Reflectance(lambda*0.5, 0, PolUnpolarized)
			if err != nil {
				t.Fatal(err)
			}
			if far < 1e-3 {
				t.Errorf("R at 0.5·lambda0 = %g, want a clearly non-zero value (quarter-wave AR is narrow)", far)
			}
		})
	}
}

// TestAiryFormulaComparison verifies the characteristic-matrix result against
// the closed-form Airy summation of one layer, for random real and complex
// indices, thicknesses, angles and both polarizations.
func TestAiryFormulaComparison(t *testing.T) {
	rng := rand.New(rand.NewSource(20240))
	const cases = 300
	for i := 0; i < cases; i++ {
		lambda := 300 + rng.Float64()*1200
		n0 := 1 + rng.Float64()*0.6
		n1 := 1.2 + rng.Float64()*2.0
		ns := 1.3 + rng.Float64()*1.0
		k1 := 0.0
		ks := 0.0
		if i%3 == 0 {
			k1 = rng.Float64() * 0.05
		}
		if i%4 == 0 {
			ks = rng.Float64() * 0.02
		}
		d := rng.Float64() * 1500
		theta := rng.Float64() * 80
		pol := PolS
		if i%2 == 1 {
			pol = PolP
		}
		c0 := complex(n0, 0)
		c1 := complex(n1, k1)
		cs := complex(ns, ks)
		st := NewStack(
			NewMaterial("incident", ConstIndex{N: n0}),
			NewMaterial("substrate", ConstIndex{N: ns, K: ks}),
			Layer{Material: NewMaterial("film", ConstIndex{N: n1, K: k1}), ThicknessNm: d},
		)
		rc, tc, R, T, err := st.polCoefficients(lambda, theta, pol)
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		wantR, wantT := airySingleLayer(lambda, c0, c1, cs, d, theta, pol)
		closeComplex(t, "r", rc, wantR, 1e-12)
		closeComplex(t, "t", tc, wantT, 1e-12)
		closeFloat(t, "R", R, real(wantR*cmplx.Conj(wantR)), 1e-12)
		wantTr := real(refAdmittance(cs, c0*complex(math.Sin(theta*math.Pi/180), 0), pol)) /
			real(refAdmittance(c0, c0*complex(math.Sin(theta*math.Pi/180), 0), pol))
		closeFloat(t, "T", T, wantTr*real(wantT*cmplx.Conj(wantT)), 1e-12)
	}
}

// TestAdmittanceRecursionComparison cross-checks the matrix product against the
// Macleod admittance recurrence for random multilayer stacks, including
// absorbing layers.
func TestAdmittanceRecursionComparison(t *testing.T) {
	rng := rand.New(rand.NewSource(777))
	for i := 0; i < 150; i++ {
		nLayers := 1 + rng.Intn(12)
		st := randomStack(rng, nLayers, i%2 == 0, 1.2, 2.6, 5, 900)
		lambda := 350 + rng.Float64()*1000
		theta := rng.Float64() * 85
		pol := PolS
		if i%3 == 0 {
			pol = PolP
		}
		rc, _, _, _, err := st.polCoefficients(lambda, theta, pol)
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		want := recursionReflectance(lambda, st.Incident.IndexNm(lambda), st.Layers, st.Substrate.IndexNm(lambda), theta, pol)
		closeComplex(t, "r", rc, want, 1e-11)
	}
}

// TestQuarterWaveStackPeakReflectance checks the analytic peak reflectance of a
// quarter-wave stack, R = ((n0-Y)/(n0+Y))^2 with Y = (nH/nL)^(2p)·ns, which the
// classical construction predicts at the design wavelength.
func TestQuarterWaveStackPeakReflectance(t *testing.T) {
	const lambda = 1064.0
	inc, sub := Air, BK7
	high, low := Ta2O5, SiO2
	for _, pairs := range []int{1, 2, 5, 10, 18} {
		st, _, err := HighReflectorPair(inc, sub, lambda, pairs, high, low, Options{MaxIter: 1})
		if err != nil {
			t.Fatalf("pairs=%d: %v", pairs, err)
		}
		// options with MaxIter: 1 still runs one iteration; rebuild the purely
		// analytic stack for the reference comparison.
		dh, _ := QWOTNm(high, lambda)
		dl, _ := QWOTNm(low, lambda)
		var layers []Layer
		for i := 0; i < pairs; i++ {
			layers = append(layers, Layer{Material: high, ThicknessNm: dh}, Layer{Material: low, ThicknessNm: dl})
		}
		analytic := NewStack(inc, sub, layers...)
		got, err := analytic.Reflectance(lambda, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		nH, nL, ns := high.N(lambda), low.N(lambda), sub.N(lambda)
		Y := math.Pow(nH/nL, float64(2*pairs)) * ns
		want := math.Pow((1-Y)/(1+Y), 2)
		closeFloat(t, "peak R", got, want, 1e-12)
		// The real template, with its (mild) refinement against the three
		// reflectance goals, must stay close to the analytic peak and be a high
		// reflector for a realistic pair count.
		refined, err := st.Reflectance(lambda, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if pairs >= 10 {
			if refined < 0.99 {
				t.Errorf("pairs=%d: refined R = %g, want > 0.99", pairs, refined)
			}
			if math.Abs(refined-want) > 1e-3 {
				t.Errorf("pairs=%d: refined R = %g drifted from the analytic peak %g", pairs, refined, want)
			}
		}
		if refined > 1 {
			t.Errorf("pairs=%d: reflectance %g above unity", pairs, refined)
		}
	}
}

// TestEnergyConservation verifies R + T + A = 1 and, for lossless stacks,
// R + T = 1 to 1e-12 over random multilayer designs, angles and polarizations.
func TestEnergyConservation(t *testing.T) {
	rng := rand.New(rand.NewSource(4242))
	for i := 0; i < 200; i++ {
		nLayers := rng.Intn(8)
		lossy := i%2 == 0
		st := randomStack(rng, nLayers, lossy, 1.1, 2.5, 10, 1200)
		lambda := 300 + rng.Float64()*1300
		theta := rng.Float64()*89.5 * sign(rng)
		for _, pol := range []Polarization{PolS, PolP, PolUnpolarized} {
			smp, err := st.Evaluate(lambda, theta, pol)
			if err != nil {
				t.Fatalf("case %d: %v", i, err)
			}
			if math.Abs(smp.R+smp.T+smp.A-1) > 1e-12 {
				t.Errorf("case %d %s: R+T+A = %.17g, want 1", i, pol, smp.R+smp.T+smp.A)
			}
			if smp.A < -1e-12 || smp.A > 1+1e-12 {
				t.Errorf("case %d %s: A = %g out of range", i, pol, smp.A)
			}
			if !lossy && math.Abs(smp.R+smp.T-1) > 1e-12 {
				t.Errorf("case %d %s: lossless R+T = %.17g, want 1", i, pol, smp.R+smp.T)
			}
			if lossy && smp.R+smp.T > 1+1e-12 {
				t.Errorf("case %d %s: R+T = %.17g > 1 for an absorbing stack", i, pol, smp.R+smp.T)
			}
		}
	}
}

// sign returns +1 or -1 at random.
func sign(rng *rand.Rand) float64 {
	if rng.Intn(2) == 0 {
		return 1
	}
	return -1
}

// TestAbsorptionFromInternalField is the physical check of the absorptance: the
// power absorbed by a lossy layer is computed by integrating the ohmic loss
// density (omega·eps0·n'·k·|E|^2) over the standing-wave field profile and must
// match A = 1 - R - T obtained from the complex amplitudes.
func TestAbsorptionFromInternalField(t *testing.T) {
	const lambda = 532.0
	cases := []struct {
		name    string
		n1, k1  float64
		d       float64
		ns      float64
		theta   float64
	}{
		{"normal incidence", 1.6, 0.02, 400, 1.52, 0},
		{"oblique", 2.1, 0.005, 250, 1.52, 45},
		{"strong absorption", 1.9, 0.15, 700, 1.52, 20},
		{"absorbing substrate", 1.5, 0.01, 300, 1.52, 0}, // substrate kept lossless below
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			st := NewStack(Air, NewMaterial("sub", ConstIndex{N: tc.ns}),
				Layer{Material: NewMaterial("absorbing film", ConstIndex{N: tc.n1, K: tc.k1}), ThicknessNm: tc.d})
			smp, err := st.Evaluate(lambda, tc.theta, PolS)
			if err != nil {
				t.Fatal(err)
			}
			// Numerical integration of the loss over the layer, s-polarization.
			const steps = 20001
			k0 := 2 * math.Pi / lambda
			depths := make([]float64, steps)
			for i := range depths {
				depths[i] = tc.d * float64(i) / float64(steps-1)
			}
			field, err := st.FieldProfile(lambda, tc.theta, PolS, depths)
			if err != nil {
				t.Fatal(err)
			}
			integral := 0.0
			for i, e := range field {
				w := 1.0
				if i == 0 || i == steps-1 {
					w = 0.5
				}
				integral += w * real(e*cmplx.Conj(e))
			}
			integral *= tc.d / float64(steps-1) // trapezoid
			n0 := math.Cos(tc.theta * math.Pi / 180)
			wantA := 2 * k0 * tc.n1 * tc.k1 * integral / n0
			// The tolerance covers the trapezoidal integration error only; the
			// two routes to the absorptance are otherwise independent.
			closeFloat(t, "A from field", wantA, smp.A, 1e-4)
			if smp.A <= 0 {
				t.Errorf("absorptance %g is not positive for an absorbing layer", smp.A)
			}
		})
	}
}

// TestBareSubstrateFresnel compares a stack without layers against the Fresnel
// equations, and checks that the p-polarized reflectance vanishes at the
// Brewster angle while the s reflectance does not.
func TestBareSubstrateFresnel(t *testing.T) {
	inc, sub := Air, BK7
	bare := NewStack(inc, sub)
	for _, lambda := range []float64{400, 550, 633, 1064} {
		n0, ns := inc.N(lambda), sub.N(lambda)
		brewster := math.Atan(ns/n0) * 180 / math.Pi
		for _, theta := range []float64{0, 15, 30, 45, 60, brewster, 75} {
			th := theta * math.Pi / 180
			ct0, ct1 := math.Cos(th), math.Sqrt(1-math.Pow(n0*math.Sin(th)/ns, 2))
			wantRs := math.Pow((n0*ct0-ns*ct1)/(n0*ct0+ns*ct1), 2)
			wantRp := math.Pow((ns*ct0-n0*ct1)/(ns*ct0+n0*ct1), 2)
			got, err := bare.Evaluate(lambda, theta, PolS)
			if err != nil {
				t.Fatal(err)
			}
			closeFloat(t, "Rs", got.R, wantRs, 1e-12)
			closeFloat(t, "Ts", got.T, 1-wantRs, 1e-12)
			gotP, err := bare.Evaluate(lambda, theta, PolP)
			if err != nil {
				t.Fatal(err)
			}
			closeFloat(t, "Rp", gotP.R, wantRp, 1e-12)
		}
		// Brewster: Rp must vanish to round-off, Rs must not.
		gotP, err := bare.Reflectance(lambda, brewster, PolP)
		if err != nil {
			t.Fatal(err)
		}
		if gotP > 1e-24 {
			t.Errorf("lambda=%g: Rp at Brewster = %g, want < 1e-24", lambda, gotP)
		}
		gotS, err := bare.Reflectance(lambda, brewster, PolS)
		if err != nil {
			t.Fatal(err)
		}
		if gotS < 0.05 {
			t.Errorf("lambda=%g: Rs at Brewster = %g, want a clearly non-zero value", lambda, gotS)
		}
		// Brewster angle of an air-glass interface is about 56.6 deg for BK7.
		if math.Abs(brewster-56.6) > 0.5 {
			t.Errorf("lambda=%g: Brewster angle %g deg is far from 56.6 deg", lambda, brewster)
		}
	}
}

// TestUnpolarizedIsMean checks that the unpolarized result is exactly the mean
// of the two polarizations.
func TestUnpolarizedIsMean(t *testing.T) {
	rng := rand.New(rand.NewSource(99))
	st := randomStack(rng, 5, false, 1.3, 2.4, 20, 500)
	for i := 0; i < 50; i++ {
		lambda := 400 + rng.Float64()*600
		theta := rng.Float64() * 80
		s, err := st.Evaluate(lambda, theta, PolS)
		if err != nil {
			t.Fatal(err)
		}
		p, err := st.Evaluate(lambda, theta, PolP)
		if err != nil {
			t.Fatal(err)
		}
		u, err := st.Evaluate(lambda, theta, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if u.R != 0.5*(s.R+p.R) || u.T != 0.5*(s.T+p.T) || u.A != 0.5*(s.A+p.A) {
			// Go may fuse or reassociate the averaging expression, so the
			// identity is checked to one ulp of the values involved rather
			// than bit-exactly.
			closeFloat(t, "unpolarized R", u.R, 0.5*(s.R+p.R), 1e-15)
			closeFloat(t, "unpolarized T", u.T, 0.5*(s.T+p.T), 1e-15)
			closeFloat(t, "unpolarized A", u.A, 0.5*(s.A+p.A), 1e-15)
		}
		if u.Rc != 0 || u.Tc != 0 {
			t.Errorf("unpolarized complex amplitudes must be zero, got %v, %v", u.Rc, u.Tc)
		}
	}
}

// TestNormalIncidencePolarizationsAgree checks that s and p coincide at normal
// incidence and that a scan agrees point-by-point with single evaluations.
func TestNormalIncidencePolarizationsAgree(t *testing.T) {
	st := NewStack(Air, BK7,
		Layer{Material: MgF2, ThicknessNm: 100},
		Layer{Material: Ta2O5, ThicknessNm: 60},
		Layer{Material: SiO2, ThicknessNm: 30},
	)
	spec, err := st.SpectrumRange(420, 680, 20, 0, PolS)
	if err != nil {
		t.Fatal(err)
	}
	specP, err := st.SpectrumRange(420, 680, 20, 0, PolP)
	if err != nil {
		t.Fatal(err)
	}
	if len(spec) != len(specP) || len(spec) == 0 {
		t.Fatalf("unexpected scan lengths %d/%d", len(spec), len(specP))
	}
	for i := range spec {
		// At normal incidence s and p are the same physical quantity; the two
		// code paths differ only by round-off.
		closeFloat(t, "s vs p R", spec[i].R, specP[i].R, 1e-14)
		closeFloat(t, "s vs p T", spec[i].T, specP[i].T, 1e-14)
		if spec[i].R != specP[i].R && math.Abs(spec[i].R-specP[i].R) > 1e-15 {
			t.Errorf("normal incidence: s and p differ at %g nm: %g vs %g", spec[i].LambdaNm, spec[i].R, specP[i].R)
		}
		want, err := st.Evaluate(spec[i].LambdaNm, 0, PolS)
		if err != nil {
			t.Fatal(err)
		}
		if want != spec[i] {
			t.Errorf("scan sample at %g nm is not the individual evaluation", spec[i].LambdaNm)
		}
	}
	// The endpoint must be included even when the step does not divide the range.
	last := spec[len(spec)-1].LambdaNm
	if last != 680 {
		t.Errorf("scan ended at %g nm, want 680", last)
	}
}

// TestAngleScanSymmetryAndGrid checks the angular scan: the grid is inclusive,
// the reflectance is symmetric about the normal and the scan matches direct
// evaluations.
func TestAngleScanSymmetryAndGrid(t *testing.T) {
	st := NewStack(Air, BK7, Layer{Material: TiO2, ThicknessNm: 55}, Layer{Material: SiO2, ThicknessNm: 90})
	scan, err := st.AngleScanRange(550, -60, 60, 15, PolS)
	if err != nil {
		t.Fatal(err)
	}
	if len(scan) != 9 {
		t.Fatalf("angle scan has %d points, want 9", len(scan))
	}
	for i := range scan {
		j := len(scan) - 1 - i
		if math.Abs(scan[i].R-scan[j].R) > 1e-14 {
			t.Errorf("asymmetric reflectance: %g deg -> %g, %g deg -> %g", scan[i].AngleDeg, scan[i].R, scan[j].AngleDeg, scan[j].R)
		}
		want, err := st.Evaluate(550, scan[i].AngleDeg, PolS)
		if err != nil {
			t.Fatal(err)
		}
		closeFloat(t, "angle scan R", scan[i].R, want.R, 1e-15)
	}
}

// TestPhaseAndThinFilmInterference checks the phase and interference
// bookkeeping on the classical phase-thickness cases, with dispersion switched
// off so that the phase thickness is the only variable: a quarter-wave layer
// gives the analytic reflectance minimum, a half-wave layer is an absentee
// layer for reflectance (but inverts the transmitted amplitude, the well-known
// pi phase shift of a half-wave layer in transmission), and three quarter waves
// reproduce the quarter-wave result exactly.
func TestPhaseAndThinFilmInterference(t *testing.T) {
	const lambda = 500.0
	const n1 = 1.38
	const ns = 1.52
	inc := NewMaterial("incident", ConstIndex{N: 1})
	sub := NewMaterial("substrate", ConstIndex{N: ns})
	film := NewMaterial("film n=1.38", ConstIndex{N: n1})
	bare := NewStack(inc, sub)
	qw := NewStack(inc, sub, Layer{Material: film, ThicknessNm: lambda / (4 * n1)})
	hw := NewStack(inc, sub, Layer{Material: film, ThicknessNm: lambda / (2 * n1)})

	bareSmp, err := bare.Evaluate(lambda, 0, PolS)
	if err != nil {
		t.Fatal(err)
	}
	qwSmp, err := qw.Evaluate(lambda, 0, PolS)
	if err != nil {
		t.Fatal(err)
	}
	// A quarter-wave layer of index n1 on a substrate ns has
	// R = ((n0·ns - n1^2)/(n0·ns + n1^2))^2 at the design wavelength.
	wantQW := math.Pow((n1*n1-ns)/(n1*n1+ns), 2)
	closeFloat(t, "quarter-wave R", qwSmp.R, wantQW, 1e-12)
	if qwSmp.R >= bareSmp.R {
		t.Errorf("quarter-wave R = %g should be below the bare substrate %g", qwSmp.R, bareSmp.R)
	}

	// Half-wave layer: identical reflectance, transmitted amplitude negated.
	hwSmp, err := hw.Evaluate(lambda, 0, PolS)
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "half-wave R", hwSmp.R, bareSmp.R, 1e-12)
	closeComplex(t, "half-wave r", hwSmp.Rc, bareSmp.Rc, 1e-12)
	closeComplex(t, "half-wave t", hwSmp.Tc, -bareSmp.Tc, 1e-12)
	// The transmitted phase differs from the bare substrate by exactly 180 deg.
	d := math.Mod(hwSmp.PhaseTdeg()-bareSmp.PhaseTdeg(), 360)
	if d < -180 {
		d += 360
	}
	if d > 180 {
		d -= 360
	}
	if math.Abs(math.Abs(d)-180) > 1e-9 {
		t.Errorf("half-wave transmitted phase differs by %g deg, want 180", math.Abs(d))
	}

	// The same layer evaluated at lambda/3 is three quarter waves thick, where
	// exp(2i·delta) = -1 again, so the reflectance must repeat exactly.
	threeQ, err := qw.Evaluate(lambda/3, 0, PolS)
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "three quarter waves R", threeQ.R, qwSmp.R, 1e-12)
}

// TestZeroThicknessAndReversal checks two structural identities: a layer of
// zero thickness is the identity, and reversing a stack (layers and media
// swapped) preserves |r| for the reversed incidence.
func TestZeroThicknessAndReversal(t *testing.T) {
	rng := rand.New(rand.NewSource(31337))
	st := randomStack(rng, 4, false, 1.3, 2.4, 30, 400)
	withZero := st.Clone()
	withZero.Layers = append(withZero.Layers, Layer{Material: SiO2, ThicknessNm: 0})
	withZero.Layers = append([]Layer{{Material: MgF2, ThicknessNm: 0}}, withZero.Layers...)
	for _, theta := range []float64{0, 30, 65} {
		a, err := st.Evaluate(550, theta, PolS)
		if err != nil {
			t.Fatal(err)
		}
		b, err := withZero.Evaluate(550, theta, PolS)
		if err != nil {
			t.Fatal(err)
		}
		closeComplex(t, "zero-thickness r", b.Rc, a.Rc, 1e-14)
		closeComplex(t, "zero-thickness t", b.Tc, a.Tc, 1e-14)
	}
}

// TestFieldProfileContinuity checks the field profile: tangential E is
// continuous across every interface, the substrate field is the transmitted
// amplitude and the incident-medium value at depth 0 is 1 + r.
func TestFieldProfileContinuity(t *testing.T) {
	const lambda = 632.8
	st := NewStack(Air, BK7,
		Layer{Material: TiO2, ThicknessNm: 60},
		Layer{Material: SiO2, ThicknessNm: 110},
		Layer{Material: Ta2O5, ThicknessNm: 25},
	)
	for _, pol := range []Polarization{PolS, PolP} {
		for _, theta := range []float64{0, 37, 70} {
			total := st.TotalThicknessNm()
			offset := 1e-7
			bounds := []float64{0.5 * st.Layers[0].ThicknessNm, st.Layers[0].ThicknessNm, total - 0.5*st.Layers[2].ThicknessNm, total}
			field, err := st.FieldProfile(lambda, theta, pol, bounds)
			if err != nil {
				t.Fatal(err)
			}
			before, err := st.FieldProfile(lambda, theta, pol, []float64{bounds[0] - offset})
			if err != nil {
				t.Fatal(err)
			}
			after, err := st.FieldProfile(lambda, theta, pol, []float64{bounds[0] + offset})
			if err != nil {
				t.Fatal(err)
			}
			closeComplex(t, "continuity", after[0], before[0], 1e-6)

			// At the front surface the total field is 1 + r.
			smp, err := st.Evaluate(lambda, theta, pol)
			if err != nil {
				t.Fatal(err)
			}
			zero, err := st.FieldProfile(lambda, theta, pol, []float64{0})
			if err != nil {
				t.Fatal(err)
			}
			closeComplex(t, "front field", zero[0], 1+smp.Rc, 1e-12)

			// In the substrate the field keeps the transmitted amplitude: deep
			// inside a lossless substrate only the phase advances, so compare
			// magnitudes there, and the exact complex value at the interface.
			deep, err := st.FieldProfile(lambda, theta, pol, []float64{total + 1e6})
			if err != nil {
				t.Fatal(err)
			}
			closeFloat(t, "substrate field magnitude", cmplx.Abs(deep[0]), cmplx.Abs(smp.Tc), 1e-12)
			// ... and exactly at the exit interface.
			exit, err := st.FieldProfile(lambda, theta, pol, []float64{total})
			if err != nil {
				t.Fatal(err)
			}
			closeComplex(t, "exit field", exit[0], smp.Tc, 1e-12)
			_ = field
		}
	}
}

// TestErrorPaths checks the input validation of the optics entry points.
func TestErrorPaths(t *testing.T) {
	good := NewStack(Air, BK7, Layer{Material: SiO2, ThicknessNm: 100})
	if _, err := good.Evaluate(0, 0, PolS); err == nil {
		t.Error("zero wavelength accepted")
	}
	if _, err := good.Evaluate(550, 90, PolS); err == nil {
		t.Error("90 degree incidence accepted")
	}
	if _, err := good.Evaluate(550, 0, Polarization(99)); err == nil {
		t.Error("invalid polarization accepted")
	}
	bad := good.Clone()
	bad.Substrate = Material{Name: "unset"}
	if _, err := bad.Evaluate(550, 0, PolS); err == nil {
		t.Error("stack without a dispersion model accepted")
	}
	neg := good.Clone()
	neg.Layers[0].ThicknessNm = -1
	if err := neg.Validate(); err == nil {
		t.Error("negative thickness accepted by Validate")
	}
	if _, err := neg.Evaluate(550, 0, PolS); err == nil {
		t.Error("negative thickness accepted by Evaluate")
	}
	if _, err := good.SpectrumRange(700, 400, 10, 0, PolS); err == nil {
		t.Error("descending spectrum range accepted")
	}
	if _, err := good.SpectrumRange(400, 700, 0, 0, PolS); err == nil {
		t.Error("zero step accepted")
	}
	if _, err := good.FieldProfile(550, 0, PolUnpolarized, []float64{0}); err == nil {
		t.Error("unpolarized field profile accepted")
	}
	if _, err := good.FieldProfile(550, 0, PolS, []float64{math.NaN()}); err == nil {
		t.Error("NaN depth accepted")
	}
	for _, label := range []string{"s", "p", "unpolarized", "avg", ""} {
		if _, err := ParsePolarization(label); err != nil {
			t.Errorf("ParsePolarization(%q): unexpected error %v", label, err)
		}
	}
	if _, err := ParsePolarization("circular"); err == nil {
		t.Error("unknown polarization label accepted")
	}
	// A Sellmeier pole must surface as an error, not as a NaN result.
	pole := NewStack(Air, BK7, Layer{Material: NewMaterial("pole", Sellmeier{B: [3]float64{1, 0, 0}, C: [3]float64{0.3025, 0, 0}}), ThicknessNm: 100})
	if _, err := pole.Evaluate(550, 0, PolS); err == nil {
		t.Error("Sellmeier pole accepted (expected a non-finite index error)")
	}
	if _, err := good.BandAverage(400, 700, 10, 0, PolS, Quantity(42)); err == nil {
		t.Error("invalid quantity accepted by BandAverage")
	}
}
