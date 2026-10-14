package geo

import (
	"math"
	"testing"
)

// testSinglet returns a biconvex N-BK7 singlet of about 76 mm focal length at
// the d line (two equal radii, 6 mm centre thickness), with the stop at its
// front surface — so the classical "no lateral colour when the stop touches
// the lens" applies to it.
func testSinglet(t *testing.T) *System {
	t.Helper()
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{6}, []Medium{Air{}, N_BK7, Air{}})
}

// TestAxialChromaticSinglet checks the axial (longitudinal) chromatic
// aberration of a singlet against the classical thin-lens estimate
// Δz = f_F − f_C ≈ −f_d/V_d, including sign, and cross-checks the value
// against real ray tracing of the axial fan at the two wavelengths.
func TestAxialChromaticSinglet(t *testing.T) {
	sys := testSinglet(t)
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	v := abbeOf(N_BK7)
	checkRel(t, "Abbe number of N-BK7", v, 64.17, 5e-3)

	dz, err := sys.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	want := -par.EFL / v
	checkRel(t, "F−C focus shift vs −f/V", dz, want, 0.05)
	if dz >= 0 {
		t.Errorf("F−C shift = %g, want negative (blue focuses closer)", dz)
	}

	// Real rays: the axial focus of the F and C lines (crossing of the
	// marginal rays) must differ by the same amount.
	focusAt := func(lambda float64) float64 {
		res := sys.Trace(NewRay(V(0, 1e-3, -100), V(0, 0, 1), lambda))
		if !res.OK() {
			t.Fatalf("trace failed at λ=%g: %v", lambda, res.Error())
		}
		return axisCrossingZ(res)
	}
	real := focusAt(WavelengthF) - focusAt(WavelengthC)
	checkRel(t, "real F−C focus shift", real, dz, 0.05)

	// FocusAt and RefocusZ are consistent definitions.
	zF, err := sys.FocusAt(math.Inf(-1), WavelengthF)
	if err != nil {
		t.Fatal(err)
	}
	zC, err := sys.FocusAt(math.Inf(-1), WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "FocusAt difference", zF-zC, dz, 1e-12)
	ref, err := sys.RefocusZ(math.Inf(-1), Wavelengthd, WavelengthF)
	if err != nil {
		t.Fatal(err)
	}
	zd, err := sys.FocusAt(math.Inf(-1), Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "RefocusZ", ref, zF-zd, 1e-12)

	// The chromatic focal-shift curve is monotonic in wavelength for a
	// simple positive lens (red focuses farther away).
	lambdas := []float64{WavelengthF, Wavelengthe, Wavelengthd, WavelengthC}
	shifts, err := sys.ChromaticFocalShift(math.Inf(-1), lambdas)
	if err != nil {
		t.Fatal(err)
	}
	if len(shifts) != len(lambdas) || shifts[0] != 0 {
		t.Fatalf("ChromaticFocalShift = %v, want the first entry zero", shifts)
	}
	for i := 1; i < len(shifts); i++ {
		if shifts[i] <= shifts[i-1] {
			t.Errorf("focus shifts not increasing with wavelength: %v", shifts)
		}
	}
	checkRel(t, "secondary point of the curve", shifts[len(shifts)-1], -dz, 0.02)
}

// TestAchromatDoublet designs and checks a cemented crown/flint achromat. The
// classic thin-lens rule φ₁/V₁ + φ₂/V₂ = 0 fixes only the ratio of the powers;
// the thickness of the two elements adds a small chromatic power of its own,
// so the design solves the real conditions — focal length at the d line equal
// to the target and exactly zero F−C power difference — for both powers with
// nested deterministic bisections. The resulting doublet must then show only a
// vanishing F−C focus shift, far below a singlet of comparable focal length.
func TestAchromatDoublet(t *testing.T) {
	nd1 := N_BK7.N(Wavelengthd)
	nd2 := N_SF5.N(Wavelengthd)
	v1, v2 := abbeOf(N_BK7), abbeOf(N_SF5)
	const target = 200.0

	build := func(phi1, phi2 float64) *System {
		R1 := (nd1 - 1) / phi1  // crown: plane-convex, its plane is the cement
		R3 := -(nd2 - 1) / phi2 // flint: plane-concave, its plane is the cement
		return mustSystem(t, Air{},
			[]Surface{
				{Profile: Sphere{R: R1}, Aperture: CircleAperture(12)},
				{Profile: Plane{}, Aperture: CircleAperture(12)},
				{Profile: Sphere{R: R3}, Aperture: CircleAperture(12)},
			},
			[]float64{4, 3},
			[]Medium{Air{}, N_BK7, N_SF5, Air{}})
	}
	efl := func(sys *System, lambda float64) float64 {
		p, err := sys.Paraxial(lambda, math.Inf(-1))
		if err != nil {
			t.Fatal(err)
		}
		return p.EFL
	}
	chromPower := func(phi1, phi2 float64) float64 {
		sys := build(phi1, phi2)
		return 1/efl(sys, WavelengthF) - 1/efl(sys, WavelengthC)
	}
	// solveEqual finds phi2 making f(phi1, phi2) = 0 by scanning a bracket
	// and bisecting it.
	solveEqual := func(f func(float64) float64, lo, hi float64) float64 {
		flo, fhi := f(lo), f(hi)
		for i := 0; i < 60 && flo*fhi > 0; i++ {
			lo, hi = lo+(hi-lo)*0.25, lo+(hi-lo)*0.75
			flo, fhi = f(lo), f(hi)
		}
		for i := 0; i < 200 && hi-lo > 1e-15*(1+math.Abs(hi)); i++ {
			mid := 0.5 * (lo + hi)
			fm := f(mid)
			if flo*fm <= 0 {
				hi, fhi = mid, fm
			} else {
				lo, flo = mid, fm
			}
		}
		_ = fhi
		return 0.5 * (lo + hi)
	}
	// Outer: the d-line focal length; inner: the F−C power match.
	phi1Nom := (1 / target) / (1 - v2/v1)
	phi2Of := func(phi1 float64) float64 {
		return solveEqual(func(p2 float64) float64 { return chromPower(phi1, p2) },
			-2.5*phi1*v2/v1, 0.5*phi1*v2/v1)
	}
	powerOf := func(phi1 float64) float64 {
		return 1/efl(build(phi1, phi2Of(phi1)), Wavelengthd) - 1/target
	}
	phi1 := solveEqual(powerOf, 0.7*phi1Nom, 1.4*phi1Nom)
	phi2 := phi2Of(phi1)
	doublet := build(phi1, phi2)
	par, err := doublet.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "achromat design focal length", par.EFL, target, 1e-4)

	dz, err := doublet.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	single := testSinglet(t)
	dSingle, err := single.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	// The design makes the focal LENGTH achromatic; the image-plane shift
	// still carries the (much smaller) chromatic movement of the principal
	// plane with the glass thickness, so it is bounded rather than zero.
	checkRel(t, "achromatic focal length", efl(doublet, WavelengthF), efl(doublet, WavelengthC), 1e-6)
	if math.Abs(dz) > 0.1 {
		t.Errorf("achromat F−C focus shift = %g mm, want a small residual", dz)
	}
	if math.Abs(dz) > 0.1*math.Abs(dSingle) {
		t.Errorf("achromat F−C shift %g is not strongly suppressed vs the singlet's %g", dz, dSingle)
	}
	if math.Abs(dSingle) < 1e-3 {
		t.Fatalf("the singlet reference shows no chromatic shift (%g)", dSingle)
	}
	// The system really is dispersive (otherwise the test would be vacuous).
	if !doublet.DispersiveSystem() {
		t.Error("the doublet must be reported as dispersive")
	}
	mirror := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: -200}, Mirror: true, Aperture: CircleAperture(50)}},
		nil, nil)
	if mirror.DispersiveSystem() {
		t.Error("a mirror in air must not be dispersive")
	}
}

// TestLateralChromatic checks the lateral (magnification) chromatic
// aberration. A lens with the stop in contact has essentially none (the
// classical "no lateral colour when the stop is at the lens"), while moving
// the stop away from the lens produces a measurable height difference between
// the F and C chief rays; the paraxial prediction and the real (aimed) chief
// rays must agree.
func TestLateralChromatic(t *testing.T) {
	f := FieldTan(0, 0.05)
	single := testSinglet(t)
	par, err := single.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	dh, err := single.LateralChromatic(f, par.ImageZ, WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	// Stop at the lens: the classic result is a vanishing lateral colour.
	checkLess(t, "lateral colour with the stop at the lens", math.Abs(dh), 0.01)

	// The same lens with the stop 30 mm in front: the chief rays of the two
	// colours land at clearly different heights.
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	displaced := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(6)},
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{30, 6}, []Medium{Air{}, Air{}, N_BK7, Air{}})
	displaced.StopIndex = 0
	par2, err := displaced.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	dh2, err := displaced.LateralChromatic(f, par2.ImageZ, WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(dh2) < 5*math.Abs(dh) {
		t.Errorf("displaced-stop lateral colour %g is not clearly larger than the in-contact %g", dh2, dh)
	}

	// Real chief rays at the two wavelengths, at the same image plane.
	hitY := func(sys *System, lambda float64) float64 {
		chief := sys.ChiefRay(f, lambda)
		res := sys.Trace(chief)
		if !res.OK() {
			t.Fatalf("chief trace failed at λ=%g: %v", lambda, res.Error())
		}
		d := res.Ray.Dir
		tt := (par2.ImageZ - res.Ray.Origin.Z) / d.Z
		return res.Ray.Origin.Y + tt*d.Y
	}
	real := hitY(displaced, WavelengthF) - hitY(displaced, WavelengthC)
	if math.Abs(dh2-real) > 0.05*math.Abs(real) {
		t.Errorf("lateral colour: paraxial %g, real chief rays %g", dh2, real)
	}
	// On axis there is no lateral colour to speak of.
	dh0, err := displaced.LateralChromatic(FieldTan(0, 0), par2.ImageZ, WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(dh0) > 1e-9 {
		t.Errorf("on-axis lateral colour = %g, want 0", dh0)
	}
}

// TestChromaticErrors checks the error paths of the chromatic helpers.
func TestChromaticErrors(t *testing.T) {
	mirror := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: -200}, Mirror: true, Aperture: CircleAperture(50)}},
		nil, nil)
	if _, err := mirror.FocusAt(math.Inf(-1), Wavelengthd); err == nil {
		t.Error("FocusAt must reject mirror systems (no paraxial engine)")
	}
	if _, err := mirror.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC); err == nil {
		t.Error("AxialChromatic must reject mirror systems")
	}
	if _, err := mirror.LateralChromatic(FieldTan(0, 0.05), -100, WavelengthF, WavelengthC); err == nil {
		t.Error("LateralChromatic must reject mirror systems")
	}
	sys := testSinglet(t)
	if _, err := sys.ChromaticFocalShift(math.Inf(-1), nil); err == nil {
		t.Error("ChromaticFocalShift must reject an empty wavelength list")
	}
}
