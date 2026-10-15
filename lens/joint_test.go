package lens

import (
	"math"
	"testing"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/geo"
)

// fresnelPower returns the Fresnel power reflectances (Rs, Rp) of a plane
// interface from medium n1 into n2 at incidence angle thetaI (radians), from
// the textbook amplitude coefficients. Total internal reflection gives 1.
func fresnelPower(n1, n2, thetaI float64) (rs, rp float64) {
	sinT := n1 * math.Sin(thetaI) / n2
	if math.Abs(sinT) >= 1 {
		return 1, 1
	}
	cosI := math.Cos(thetaI)
	cosT := math.Sqrt(1 - sinT*sinT)
	a := (n1*cosI - n2*cosT) / (n1*cosI + n2*cosT)
	b := (n2*cosI - n1*cosT) / (n2*cosI + n1*cosT)
	return a * a, b * b
}

// normalIncidenceR returns the normal-incidence Fresnel reflectance ((n1−n2)/(n1+n2))².
func normalIncidenceR(n1, n2 float64) float64 {
	d := (n1 - n2) / (n1 + n2)
	return d * d
}

// planeInterface is a single plane surface between air and N-BK7 (the test
// system of the single-interface checks). A ray launched at polar angle θ from
// z = −5 mm through the origin meets the surface at incidence θ.
func planeInterface() *Prescription {
	return &Prescription{
		ObjectZMM: math.Inf(-1),
		Cards:     []Card{{Name: "plane", RadiusMM: 0, Glass: "N-BK7", SemiAperMM: 0}},
	}
}

// rayAtAngle launches the ray of planeInterface at incidence thetaDeg.
func rayAtAngle(thetaDeg, lambdaNm float64) geo.Ray {
	th := thetaDeg * math.Pi / 180
	origin := geo.V(-5*math.Tan(th), 0, -5)
	return geo.NewRay(origin, geo.V(math.Sin(th), 0, math.Cos(th)), lambdaNm)
}

// TestJointUncoatedNormalIncidence checks acceptance (a): for an uncoated
// system, the on-axis chief ray at normal incidence transmits Π(1 − R_i), where
// R_i = ((n_before − n_after)/(n_before + n_after))² is the Fresnel reflectance
// of each interface from the geometric indices. Acceptance: |T − Π(1 − R_i)| ≤ 1e-12.
func TestJointUncoatedNormalIncidence(t *testing.T) {
	p, _, err := DesignAchromat(100, "N-BK7", "F2", 10, &DesignOptions{SkipPolish: true})
	if err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	for _, lam := range []float64{geo.WavelengthF, geo.Wavelengthd, geo.WavelengthC} {
		par, err := sys.Paraxial(lam, p.ObjectZMM)
		if err != nil {
			t.Fatal(err)
		}
		chief, err := launchRay(par, p.ObjectZMM, Field{0, 0}, Sample{0, 0})
		if err != nil {
			t.Fatal(err)
		}
		chief.WavelengthNm = lam
		rt, err := TraceCoated(sys, chief, nil)
		if err != nil {
			t.Fatal(err)
		}
		if !rt.Passed || len(rt.Hits) != 3 {
			t.Fatalf("%g nm: chief ray did not traverse the three interfaces (passed %v, hits %d)", lam, rt.Passed, len(rt.Hits))
		}
		want := 1.0
		for i, h := range rt.Hits {
			if h.IncDeg != 0 {
				t.Errorf("%g nm surface %d: chief incidence %g°, want exactly 0", lam, i, h.IncDeg)
			}
			ri := normalIncidenceR(sys.Mediums[h.Surface].N(lam), sys.Mediums[h.Surface+1].N(lam))
			if math.Abs(h.Rs-ri) > 1e-12 || math.Abs(h.Rp-ri) > 1e-12 {
				t.Errorf("%g nm surface %d: TMM R (%.15g, %.15g) vs Fresnel %.15g", lam, i, h.Rs, h.Rp, ri)
			}
			want *= 1 - ri
		}
		t.Logf("%g nm: T = %.15f, Π(1−R_i) = %.15f, |diff| = %.2e", lam, rt.TS, want, math.Abs(rt.TS-want))
		if math.Abs(rt.TS-want) > 1e-12 || math.Abs(rt.TP-want) > 1e-12 {
			t.Errorf("%g nm: T (s %.15g, p %.15g) vs Π(1−R_i) %.15g: difference > 1e-12", lam, rt.TS, rt.TP, want)
		}
	}
	// End to end through Transmittance (launch, trace, aggregation).
	tab, err := p.Transmittance(TransmittanceSpec{
		Wavelengths: []float64{geo.Wavelengthd},
		Fields:      []Field{{0, 0}},
		Samples:     []Sample{{0, 0}},
	})
	if err != nil {
		t.Fatal(err)
	}
	want := 1.0
	for i := 0; i < 3; i++ {
		want *= 1 - normalIncidenceR(sys.Mediums[i].N(geo.Wavelengthd), sys.Mediums[i+1].N(geo.Wavelengthd))
	}
	for k, pol := range tab.Pols {
		if math.Abs(tab.T[0][0][k]-want) > 1e-12 {
			t.Errorf("Transmittance %v: %.15g vs Π(1−R_i) %.15g", pol, tab.T[0][0][k], want)
		}
	}
}

// TestJointIdealAR checks acceptance (b): with an ideal quarter-wave
// antireflection layer on every interface (index √(n_before·n_after), optical
// thickness λ0/4 at the design wavelength λ0), the normal-incidence chief ray
// has R = 0 on every interface and T = 1. Acceptance: |T − 1| ≤ 1e-12.
func TestJointIdealAR(t *testing.T) {
	p, _, err := DesignAchromat(100, "N-BK7", "F2", 10, &DesignOptions{SkipPolish: true})
	if err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	lam0 := geo.Wavelengthd
	coatings := make([]Coating, len(sys.Surfaces))
	for i := range coatings {
		nb := sys.Mediums[i].N(lam0)
		na := sys.Mediums[i+1].N(lam0)
		nAR := math.Sqrt(nb * na)
		coatings[i] = Coating{
			Name:   "ideal quarter-wave AR",
			Layers: []coating.Layer{{Material: coating.NewMaterial("AR", coating.ConstIndex{N: nAR}), ThicknessNm: lam0 / (4 * nAR)}},
		}
	}
	par, err := sys.Paraxial(lam0, p.ObjectZMM)
	if err != nil {
		t.Fatal(err)
	}
	chief, err := launchRay(par, p.ObjectZMM, Field{0, 0}, Sample{0, 0})
	if err != nil {
		t.Fatal(err)
	}
	chief.WavelengthNm = lam0
	rt, err := TraceCoated(sys, chief, coatings)
	if err != nil {
		t.Fatal(err)
	}
	for i, h := range rt.Hits {
		t.Logf("surface %d: R_s = %.3e, R_p = %.3e (ideal AR)", i, h.Rs, h.Rp)
		if h.Rs > 1e-12 || h.Rp > 1e-12 {
			t.Errorf("surface %d: residual reflectance %.3e / %.3e > 1e-12", i, h.Rs, h.Rp)
		}
	}
	t.Logf("ideal AR system: T_s = %.15f, T_p = %.15f", rt.TS, rt.TP)
	if math.Abs(rt.TS-1) > 1e-12 || math.Abs(rt.TP-1) > 1e-12 {
		t.Errorf("ideal AR: T_s %.15g, T_p %.15g, want 1 within 1e-12", rt.TS, rt.TP)
	}
}

// TestJointSingleInterfaceVsCoating checks acceptance (c): for a single coated
// interface, the per-ray R and T that the ray tracing produces at each
// wavelength, incidence angle and polarization agree with a direct evaluation
// of the coating package on the same stack, point by point, to 1e-12. The
// traced incidence angle is checked against the launch angle as well.
func TestJointSingleInterfaceVsCoating(t *testing.T) {
	mgf2, ok1 := coating.LookupMaterial("MgF2")
	ta2o5, ok2 := coating.LookupMaterial("Ta2O5")
	if !ok1 || !ok2 {
		t.Fatal("library materials missing")
	}
	layers := []coating.Layer{{Material: mgf2, ThicknessNm: 100}, {Material: ta2o5, ThicknessNm: 60}}
	c := Coating{Name: "two-layer test stack", Layers: layers}
	coatings := []Coating{c}
	sys, err := planeInterface().System()
	if err != nil {
		t.Fatal(err)
	}
	g, ok := geo.Glass("N-BK7")
	if !ok {
		t.Fatal("N-BK7 missing")
	}
	points := 0
	worst := 0.0
	for _, lam := range []float64{450, 550, 650} {
		direct := coating.Stack{
			Incident:  coating.Air,
			Substrate: coating.NewMaterial("N-BK7", coating.ConstIndex{N: g.N(lam)}),
			Layers:    layers,
		}
		for theta := 0.0; theta <= 80; theta += 5 {
			rt, err := TraceCoated(sys, rayAtAngle(theta, lam), coatings)
			if err != nil {
				t.Fatal(err)
			}
			if !rt.Passed || len(rt.Hits) != 1 {
				t.Fatalf("θ=%g° at %g nm: ray did not reach the interface", theta, lam)
			}
			h := rt.Hits[0]
			if d := math.Abs(h.IncDeg - theta); d > 1e-9 {
				t.Fatalf("traced incidence %.12g° differs from the launch angle %g° by %.2e", h.IncDeg, theta, d)
			}
			for _, pol := range []coating.Polarization{coating.PolS, coating.PolP} {
				smp, err := direct.Evaluate(lam, h.IncDeg, pol)
				if err != nil {
					t.Fatal(err)
				}
				gotR, gotT := h.Rs, h.Ts
				if pol == coating.PolP {
					gotR, gotT = h.Rp, h.Tp
				}
				dR := math.Abs(gotR - smp.R)
				dT := math.Abs(gotT - smp.T)
				worst = math.Max(worst, math.Max(dR, dT))
				points++
				if dR > 1e-12 || dT > 1e-12 {
					t.Errorf("%g nm θ=%g° %v: joint (R %.15g, T %.15g) vs coating (R %.15g, T %.15g)",
						lam, theta, pol, gotR, gotT, smp.R, smp.T)
				}
			}
		}
	}
	t.Logf("single interface: %d point comparisons, worst |difference| = %.2e", points, worst)
}

// TestJointFresnelOblique checks acceptance (d): for an uncoated interface at
// oblique incidence, the s and p reflectances (and transmittances, 1 − R for a
// lossless glass) of the traced rays agree with the Fresnel formulas at the
// traced angle, to 1e-9. The F, d and C lines and angles 0–85° are covered.
func TestJointFresnelOblique(t *testing.T) {
	sys, err := planeInterface().System()
	if err != nil {
		t.Fatal(err)
	}
	g, _ := geo.Glass("N-BK7")
	worst := 0.0
	for _, lam := range []float64{geo.WavelengthF, geo.Wavelengthd, geo.WavelengthC} {
		n2 := g.N(lam)
		for theta := 0.0; theta <= 85; theta += 5 {
			rt, err := TraceCoated(sys, rayAtAngle(theta, lam), nil)
			if err != nil {
				t.Fatal(err)
			}
			if !rt.Passed || len(rt.Hits) != 1 {
				t.Fatalf("θ=%g° at %g nm: ray did not reach the interface", theta, lam)
			}
			h := rt.Hits[0]
			rs, rp := fresnelPower(1, n2, h.IncDeg*math.Pi/180)
			worst = math.Max(worst, math.Max(math.Abs(h.Rs-rs), math.Abs(h.Rp-rp)))
			if math.Abs(h.Rs-rs) > 1e-9 || math.Abs(h.Rp-rp) > 1e-9 {
				t.Errorf("%g nm θ=%g°: TMM (Rs %.12g, Rp %.12g) vs Fresnel (%.12g, %.12g)",
					lam, theta, h.Rs, h.Rp, rs, rp)
			}
			if math.Abs(h.Ts-(1-rs)) > 1e-9 || math.Abs(h.Tp-(1-rp)) > 1e-9 {
				t.Errorf("%g nm θ=%g°: lossless T (%.12g, %.12g) vs 1−R (%.12g, %.12g)",
					lam, theta, h.Ts, h.Tp, 1-rs, 1-rp)
			}
		}
	}
	t.Logf("oblique Fresnel: worst |difference| over 3 lines × 18 angles = %.2e", worst)
}

// TestJointSpectrumAndPolarization checks the bookkeeping of the transmittance
// table on an AR-coated doublet: unpolarized T equals the mean of the s and p
// values, Throughput equals T × Pass, a field with no coating has coating-free
// Fresnel values, and the coated system transmits more than the uncoated one
// over the visible band (spectrum logged for the report).
func TestJointSpectrumAndPolarization(t *testing.T) {
	p, _, err := DesignAchromat(100, "N-BK7", "F2", 10, &DesignOptions{SkipPolish: true})
	if err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	var wls []float64
	for lam := 400.0; lam <= 700; lam += 50 {
		wls = append(wls, lam)
	}
	coatings := make([]Coating, len(sys.Surfaces))
	for i := range coatings {
		lam0 := geo.Wavelengthd
		nAR := math.Sqrt(sys.Mediums[i].N(lam0) * sys.Mediums[i+1].N(lam0))
		coatings[i] = Coating{Layers: []coating.Layer{{Material: coating.NewMaterial("AR", coating.ConstIndex{N: nAR}), ThicknessNm: lam0 / (4 * nAR)}}}
	}
	fields := []Field{{0, 0}, {0, 0.0087}}
	bare, err := p.Transmittance(TransmittanceSpec{Wavelengths: wls, Fields: fields})
	if err != nil {
		t.Fatal(err)
	}
	coated, err := p.Transmittance(TransmittanceSpec{Wavelengths: wls, Fields: fields, Coatings: coatings})
	if err != nil {
		t.Fatal(err)
	}
	unpol := 2
	for w, lam := range wls {
		for f := range fields {
			uncT := bare.T[w][f][unpol]
			cT := coated.T[w][f][unpol]
			t.Logf("%4.0f nm field %d: uncoated T %.5f, AR-coated T %.5f, pass %.2f",
				lam, f, uncT, cT, coated.Pass[w][f])
			if d := math.Abs(cT - 0.5*(coated.T[w][f][0]+coated.T[w][f][1])); d > 1e-12 {
				t.Errorf("%g nm field %d: unpolarized T differs from (s+p)/2 by %.2e", lam, f, d)
			}
			if d := math.Abs(coated.Throughput[w][f][unpol] - cT*coated.Pass[w][f]); d > 1e-15 {
				t.Errorf("%g nm field %d: throughput != T·pass (%.2e)", lam, f, d)
			}
		}
	}
	// Mean over the visible band, on-axis: the AR coating must not lose light.
	var meanBare, meanCoated float64
	for w := range wls {
		meanBare += bare.T[w][0][unpol]
		meanCoated += coated.T[w][0][unpol]
	}
	meanBare /= float64(len(wls))
	meanCoated /= float64(len(wls))
	t.Logf("on-axis band mean (400–700 nm): uncoated %.5f, AR-coated %.5f", meanBare, meanCoated)
	if !(meanCoated > meanBare) {
		t.Errorf("AR-coated band mean %.5f not above uncoated %.5f", meanCoated, meanBare)
	}
}

// TestJointRejections pins the guards: a mirror surface is outside the
// transmittance model, and more coatings than surfaces is an error.
func TestJointRejections(t *testing.T) {
	mirror := &Prescription{ObjectZMM: math.Inf(-1), Cards: []Card{
		{RadiusMM: -100, Glass: "", SemiAperMM: 10, Mirror: true},
		{RadiusMM: 0, Glass: "", SemiAperMM: 10},
	}}
	sys, err := mirror.System()
	if err != nil {
		t.Fatal(err)
	}
	r := rayAtAngle(0, geo.Wavelengthd)
	if _, err := TraceCoated(sys, r, nil); err == nil {
		t.Error("mirror surface accepted")
	}
	p := planeInterface()
	if _, err := p.Transmittance(TransmittanceSpec{Coatings: make([]Coating, 2)}); err == nil {
		t.Error("two coatings for a one-surface system accepted")
	}
}
