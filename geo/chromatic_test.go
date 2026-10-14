package geo

import (
	"math"
	"testing"
)

// singlet151 returns a biconvex N-BK7 singlet of about 151 mm focal length at
// the d line (two equal radii, 6 mm centre thickness).
func singlet151(t *testing.T) *System {
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
	sys := singlet151(t)
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

// TestAchromatDoublet checks the classical achromat condition: a cemented
// crown/flint doublet whose powers satisfy φ₁/V₁ + φ₂/V₂ = 0 has a vanishing
// F−C focus shift while a single crown lens of similar focal length shows
// millimetres of it.
func TestAchromatDoublet(t *testing.T) {
	nd1 := N_BK7.N(Wavelengthd)
	nd2 := N_SF5.N(Wavelengthd)
	v1, v2 := abbeOf(N_BK7), abbeOf(N_SF5)
	const f = 200.0 // target focal length of the doublet at the d line
	// Split the power according to the achromat condition.
	phi1 := (1 / f) / (1 - v2/v1)
	phi2 := -phi1 * v2 / v1
	// Crown: plane-convex, its plane surface is the cement.
	R1 := (nd1 - 1) / phi1
	// Flint: plane-concave, its plane surface is the cement.
	R3 := -(nd2 - 1) / phi2
	doublet := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R1}, Aperture: CircleAperture(12)},
			{Profile: Plane{}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: R3}, Aperture: CircleAperture(12)},
		},
		[]float64{4, 3},
		[]Medium{Air{}, N_BK7, N_SF5, Air{}})
	par, err := doublet.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "doublet EFL", par.EFL, f, 0.02)

	dz, err := doublet.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	// The remaining shift is the secondary spectrum; the achromat must beat
	// the equivalent singlet by a large factor.
	single := singlet151(t)
	dSingle, err := single.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(dz) > 0.05*math.Abs(dSingle) {
		t.Errorf("achromat F−C shift %g is not strongly suppressed vs the singlet's %g", dz, dSingle)
	}
	if math.Abs(dz) > 0.05 {
		t.Errorf("achromat F−C shift = %g mm, want a small secondary spectrum", dz)
	}
	// The system really is dispersive (otherwise the test would be vacuous).
	if !doublet.DispersiveSystem() {
		t.Error("the doublet must be reported as dispersive")
	}
	// A mirror system is not dispersive.
	mirror := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: -200}, Mirror: true, Aperture: CircleAperture(50)}},
		nil, nil)
	if mirror.DispersiveSystem() {
		t.Error("a mirror in air must not be dispersive")
	}
}

// TestLateralChromatic checks the lateral (magnification) chromatic
// aberration: the paraxial chief-ray image height changes with wavelength, and
// the real chief rays reproduce the change.
func TestLateralChromatic(t *testing.T) {
	sys := singlet151(t)
	f := FieldTan(0, 0.05)
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	dh, err := sys.LateralChromatic(f, par.ImageZ, WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	// Paraxial estimate: Δh = (f_F − f_C)·tanθ = Δz·tanθ.
	dz, err := sys.AxialChromatic(math.Inf(-1), WavelengthF, WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "lateral vs axial·tanθ", dh, dz*f.TanY, 0.05)

	// Real chief rays at the two wavelengths.
	hitY := func(lambda float64) float64 {
		chief := sys.ChiefRay(f, lambda)
		res := sys.Trace(chief)
		if !res.OK() {
			t.Fatalf("chief trace failed at λ=%g: %v", lambda, res.Error())
		}
		// The chief exits into air; propagate to the image plane.
		d := res.Ray.Dir
		tt := (par.ImageZ - res.Ray.Origin.Z) / d.Z
		return res.Ray.Origin.Y + tt*d.Y
	}
	real := hitY(WavelengthF) - hitY(WavelengthC)
	if math.Abs(dh-real) > 0.05*math.Abs(real) {
		t.Errorf("lateral colour: paraxial %g, real chief rays %g", dh, real)
	}
	// On axis there is no lateral colour to speak of.
	dh0, err := sys.LateralChromatic(FieldTan(0, 0), par.ImageZ, WavelengthF, WavelengthC)
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
	sys := singlet151(t)
	if _, err := sys.ChromaticFocalShift(math.Inf(-1), nil); err == nil {
		t.Error("ChromaticFocalShift must reject an empty wavelength list")
	}
}
