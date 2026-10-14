package geo

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// helpers: least-squares fit of an OPD map onto a third-order basis
// ---------------------------------------------------------------------------

// basisW is the fit basis of the OPD checks: the monomials of the third-order
// wavefront polynomial in a meridional frame (field along +y), up to the
// orders that the sampled aberrations can populate:
//
//	1, px, py, ρ², ρ²cos2θ, ρ²sin2θ, ρ³cosθ, ρ³sinθ, ρ⁴, ρ⁴cos2θ,
//	ρ⁵cosθ, ρ⁵sinθ, ρ⁶
func basisW(p PupilPoint) []float64 {
	px, py := p.PX, p.PY
	r2 := px*px + py*py
	r := math.Sqrt(r2)
	ct, st := 1.0, 0.0
	if r > 0 {
		ct, st = px/r, py/r
	}
	cos2t := ct*ct - st*st
	sin2t := 2 * ct * st
	return []float64{1, px, py, r2, r2 * cos2t, r2 * sin2t,
		r2 * r * ct, r2 * r * st, r2 * r2, r2 * r2 * cos2t,
		r2 * r2 * r * ct, r2 * r2 * r * st, r2 * r2 * r2}
}

// fitBasis fits the coefficients of the given basis onto the samples of an
// OPD grid (NaN samples excluded) with normal equations.
func fitBasis(pts []PupilPoint, v []float64, basis func(PupilPoint) []float64) []float64 {
	nf := len(basis(PupilPoint{}))
	a := make([][]float64, nf)
	for i := range a {
		a[i] = make([]float64, nf)
	}
	b := make([]float64, nf)
	for i, p := range pts {
		if i >= len(v) || math.IsNaN(v[i]) {
			continue
		}
		f := basis(p)
		for r := 0; r < nf; r++ {
			for c := 0; c < nf; c++ {
				a[r][c] += f[r] * f[c]
			}
			b[r] += f[r] * v[i]
		}
	}
	for col := 0; col < nf; col++ {
		piv := col
		for r := col + 1; r < nf; r++ {
			if math.Abs(a[r][col]) > math.Abs(a[piv][col]) {
				piv = r
			}
		}
		a[col], a[piv] = a[piv], a[col]
		b[col], b[piv] = b[piv], b[col]
		for r := 0; r < nf; r++ {
			if r == col {
				continue
			}
			f := a[r][col] / a[col][col]
			for c := col; c < nf; c++ {
				a[r][c] -= f * a[col][c]
			}
			b[r] -= f * b[col]
		}
	}
	for i := range b {
		b[i] /= a[i][i]
	}
	return b
}

// ---------------------------------------------------------------------------
// structural checks of the sums
// ---------------------------------------------------------------------------

// TestSeidelEllipsoidFociZero checks that every third-order sum vanishes for
// an ellipsoidal mirror with the object at one of its foci: the ellipse images
// one focus exactly onto the other at any aperture, so there is no spherical
// aberration at all, and the axial field makes the field-dependent sums
// vanish identically. This exercises the conic (aspheric) term of S1 for a
// conic constant other than the parabola's −1.
func TestSeidelEllipsoidFociZero(t *testing.T) {
	const (
		Rp = -200.0
		K  = -0.75
	)
	e := math.Sqrt(-K)
	a := -Rp / (1 + K)
	far := -a - a*e
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: Rp, K: K}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	sd, err := sys.SeidelCoefficients(Wavelengthd, FieldHeight(far, 0, 0), 55)
	if err != nil {
		t.Fatal(err)
	}
	scale := math.Abs(sd.MarginalHeight) // y_k sets the natural scale of the sums
	for name, v := range map[string]float64{"S1": sd.S1, "S2": sd.S2, "S3": sd.S3, "S4": sd.S4, "S5": sd.S5} {
		if math.Abs(v) > 1e-12*scale {
			t.Errorf("%s = %g for a stigmatic ellipsoid, want 0 (scale %g)", name, v, scale)
		}
	}
	// Mirror Petzval: c·Δ(1/n) with Δ(1/n) = −2/n = −2 for n = 1.
	checkRel(t, "ellipsoid Petzval", sd.Petzval, -2/Rp, 1e-12)
}

// TestSeidelParabolaZeroSA checks that the parabola has no third-order
// spherical aberration, and that the same surface built as an EvenAsphere
// with the equivalent r⁴ coefficient reproduces it.
func TestSeidelParabolaZeroSA(t *testing.T) {
	const Rp = -200.0
	parab := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: Rp, K: -1}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	sphere := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: Rp}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	f := FieldTan(0, 0)
	sdP, err := parab.SeidelCoefficients(Wavelengthd, f, 40)
	if err != nil {
		t.Fatal(err)
	}
	sdS, err := sphere.SeidelCoefficients(Wavelengthd, f, 40)
	if err != nil {
		t.Fatal(err)
	}
	checkClose(t, "parabola S1", sdP.S1, 0, 1e-12*math.Abs(sdS.S1))
	if sdS.S1 == 0 {
		t.Fatal("the spherical mirror must have spherical aberration")
	}
	// The parabola's r⁴ sag term is −c³/8 (K = −1 on a spherical base).
	c := 1 / Rp
	asphere := mustSystem(t, Air{},
		[]Surface{{Profile: EvenAsphere{R: Rp, K: 0, Coeffs: []float64{-c * c * c / 8}}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	sdA, err := asphere.SeidelCoefficients(Wavelengthd, f, 40)
	if err != nil {
		t.Fatal(err)
	}
	checkClose(t, "asphere S1 = parabola S1", sdA.S1, sdP.S1, 1e-9*math.Abs(sdS.S1))
}

// TestSeidelThinLensNoDistortion checks the classical result that a thin lens
// with the stop in contact has no distortion (S5 = 0 and, by the validated
// conversion, no image-height error of the chief ray).
func TestSeidelThinLensNoDistortion(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{0}, []Medium{Air{}, N_BK7, Air{}})
	sys.StopIndex = 0
	f := FieldTan(0, 0.05)
	sd, err := sys.SeidelCoefficients(Wavelengthd, f, 10)
	if err != nil {
		t.Fatal(err)
	}
	scale := math.Abs(sd.S2) + math.Abs(sd.S3) + 1e-30
	if math.Abs(sd.S5) > 1e-10*scale {
		t.Errorf("S5 = %g for a thin lens with the stop at the lens, want 0", sd.S5)
	}
	if math.Abs(sd.RelativeDistortion()) > 1e-10 {
		t.Errorf("relative distortion = %g, want 0", sd.RelativeDistortion())
	}
	// The real chief ray agrees: the paraxial image height is exact.
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	d, err := sys.Distortion(f, par.ImageZ, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	checkLess(t, "real thin-lens distortion", math.Abs(d), 1e-10)
}

// TestSeidelSymmetricSystem checks the classical symmetry result: a system
// that mirrors itself about the stop plane, used at its symmetric (unit
// magnification) conjugates, has no coma and no distortion — the aberration
// function must be antisymmetric in the field about the symmetry plane, which
// kills the even-in-field terms (S2, S5). The check compares the symmetric
// stop against the same lens with the stop displaced to its front surface: the
// symmetric configuration must be smaller by orders of magnitude, and the real
// (aimed) chief ray must land on the paraxial image height.
func TestSeidelSymmetricSystem(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	build := func(stop int) *System {
		sys := mustSystem(t, Air{},
			[]Surface{
				{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
				{Profile: Plane{}, Aperture: CircleAperture(4)}, // middle (or stop)
				{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
			},
			[]float64{3, 3},
			[]Medium{Air{}, N_BK7, N_BK7, Air{}})
		sys.StopIndex = stop
		return sys
	}
	sym := build(1)
	objZ, err := magnificationMinusOne(sym, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	par, err := sym.Paraxial(Wavelengthd, objZ)
	if err != nil {
		t.Fatal(err)
	}
	checkClose(t, "unit magnification", par.Magnification, -1, 1e-8)

	h := 0.5
	f := FieldHeight(objZ, 0, h)
	mid, err := sym.SeidelCoefficients(Wavelengthd, f, 4)
	if err != nil {
		t.Fatal(err)
	}
	// The same lens with the stop at its front surface: no longer symmetric
	// about the stop, so the even-in-field sums become large.
	front := build(0)
	asym, err := front.SeidelCoefficients(Wavelengthd, f, 4)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(asym.S2) < 1e-6 || math.Abs(asym.S5) < 1e-9 {
		t.Fatalf("the displaced-stop reference must have coma/distortion: S2=%g S5=%g", asym.S2, asym.S5)
	}
	if r := math.Abs(mid.S2) / math.Abs(asym.S2); r > 1e-2 {
		t.Errorf("symmetric S2 is %.3g of the displaced-stop S2 (%.3g), want a strong suppression", r, math.Abs(mid.S2))
	}
	if r := math.Abs(mid.S5) / math.Abs(asym.S5); r > 1e-4 {
		t.Errorf("symmetric S5 is %.3g of the displaced-stop S5 (%.3g), want a strong suppression", r, math.Abs(mid.S5))
	}
	// The real chief ray lands on the paraxial image height.
	d, err := sym.Distortion(f, par.ImageZ, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	checkLess(t, "symmetric system real distortion", math.Abs(d), 1e-6)
}

// magnificationMinusOne returns the object plane of transverse magnification
// −1 (the symmetric conjugate of a symmetric system), by a deterministic scan
// plus bisection on the paraxial magnification.
func magnificationMinusOne(sys *System, lambdaNm float64) (float64, error) {
	mag := func(z float64) (float64, error) {
		p, err := sys.Paraxial(lambdaNm, z)
		if err != nil {
			return 0, err
		}
		return p.Magnification, nil
	}
	// The unit-magnification object plane of a symmetric system sits near
	// −2·EFL (in front of the first surface).
	pInf, err := sys.Paraxial(lambdaNm, math.Inf(-1))
	if err != nil {
		return 0, err
	}
	lo, hi := -4*math.Abs(pInf.EFL), -0.5*math.Abs(pInf.EFL)
	flo, err := mag(lo)
	if err != nil {
		return 0, err
	}
	fhi, err := mag(hi)
	if err != nil {
		return 0, err
	}
	for i := 0; i < 200 && (flo+1)*(fhi+1) > 0; i++ {
		mid := 0.5 * (lo + hi)
		fm, err := mag(mid)
		if err != nil {
			return 0, err
		}
		if (flo+1)*(fm+1) <= 0 {
			hi, fhi = mid, fm
		} else {
			lo, flo = mid, fm
		}
	}
	// Secant refinement.
	x0, x1 := lo, hi
	f0, f1 := flo+1, fhi+1
	for i := 0; i < 60 && math.Abs(f1) > 1e-12; i++ {
		if f1 == f0 {
			break
		}
		x2 := x1 - f1*(x1-x0)/(f1-f0)
		f2, err := mag(x2)
		if err != nil {
			return 0, err
		}
		x0, f0, x1, f1 = x1, f1, x2, f2+1
	}
	return x1, nil
}

// TestSeidelPetzval checks the Petzval sum against its closed forms: a thin
// lens in air has Σ c·Δ(1/n) = −φ/n, and a mirror has Δ(1/n) = −2/n.
func TestSeidelPetzval(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{6}, []Medium{Air{}, N_BK7, Air{}})
	sd, err := sys.SeidelCoefficients(Wavelengthd, FieldTan(0, 0), 10)
	if err != nil {
		t.Fatal(err)
	}
	phi := (n - 1) * (1/R - (-1 / R)) // total thin-lens power
	checkRel(t, "singlet Petzval", sd.Petzval, -phi/n, 1e-12)
	checkRel(t, "singlet Petzval radius", sd.PetzvalRadius(), -n/phi, 1e-12)

	// Mirror: P = Σ c·(−2/n).
	const Rm = -200.0
	mirror := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: Rm}, Mirror: true, Aperture: CircleAperture(50)}},
		nil, nil)
	sdM, err := mirror.SeidelCoefficients(Wavelengthd, FieldTan(0, 0), 40)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "mirror Petzval", sdM.Petzval, (1/Rm)*(-2.0), 1e-12)
}

// TestSeidelFieldMagnitude checks the rotational reduction of the field: a
// field of the same magnitude but a different azimuth gives the same sums
// (the meridional magnitude is what enters).
func TestSeidelFieldMagnitude(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
			{Profile: Plane{}, Aperture: CircleAperture(3)},
		},
		[]float64{6, 30}, []Medium{Air{}, N_BK7, Air{}, Air{}})
	sys.StopIndex = 2
	a, err := sys.SeidelCoefficients(Wavelengthd, FieldTan(0, 0.05), 5)
	if err != nil {
		t.Fatal(err)
	}
	b, err := sys.SeidelCoefficients(Wavelengthd, FieldTan(0.03, 0.04), 5) // same magnitude
	if err != nil {
		t.Fatal(err)
	}
	for name, pair := range map[string][2]float64{
		"S1": {a.S1, b.S1}, "S2": {a.S2, b.S2}, "S3": {a.S3, b.S3},
		"S4": {a.S4, b.S4}, "S5": {a.S5, b.S5}, "P": {a.Petzval, b.Petzval},
	} {
		checkRel(t, name+" under field rotation", pair[1], pair[0], 1e-12)
	}
	checkClose(t, "H under field rotation", b.H, a.H, 1e-12)
}

// ---------------------------------------------------------------------------
// the sums against real ray tracing (small-aperture limit)
// ---------------------------------------------------------------------------

// TestSeidelWavefrontVersusOPD is the central validation of the aberration
// layer: the third-order wavefront coefficients derived from the Seidel sums
// are compared with the exact OPD maps computed by real ray tracing, for
// three configurations (a mirror on axis, an off-axis mirror and an off-axis
// lens, all sampled in the chief-aligned pupil frame):
//
//	W040 = −S1/8           ρ⁴        (spherical aberration)
//	W131 = −S2/2           ρ³sinθ    (coma of a y-field)
//	W222 = −S3/2           ρ²sin²θ   (astigmatism)
//	W220 = −(S3+S4)/4      ρ²        (field curvature)
//	W311 = −S5/2           ρ sinθ    (distortion)
//
// Each fitted coefficient must agree with the prediction of the sums in the
// small-aperture limit.
func TestSeidelWavefrontVersusOPD(t *testing.T) {
	const (
		Rp  = -200.0
		rPu = 40.0
		zPu = -50.0
	)
	type config struct {
		name  string
		sys   *System
		field Field
	}
	mirror := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: Rp}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	parab := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: Rp, K: -1}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	Rl := 155.04
	lens := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(6)},
			{Profile: Sphere{R: Rl}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -Rl}, Aperture: CircleAperture(12)},
		},
		[]float64{30, 6}, []Medium{Air{}, Air{}, N_BK7, Air{}})
	lens.StopIndex = 0
	lensPar, err := lens.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	lensPup := Pupil{Z: lensPar.EntrancePupilZ, Radius: lensPar.EntrancePupilRadius}

	for _, tc := range []config{
		{"mirror on axis", mirror, FieldTan(0, 0)},
		{"mirror off axis", mirror, FieldTan(0, 0.03)},
		{"parabola off axis", parab, FieldTan(0, 0.03)},
		{"lens off axis", lens, FieldTan(0, 0.05)},
	} {
		// Pupil frame: chief-aligned, and the image plane is the paraxial
		// focus of the configuration.
		var pup Pupil
		var zImg float64
		var pupilRadius float64
		if tc.sys == lens {
			pup, zImg = lensPup, lensPar.ImageZ
			pupilRadius = lensPup.Radius
			hPar, err := lens.ParaxialImageHeight(tc.field, zImg, Wavelengthd)
			if err != nil {
				t.Fatal(err)
			}
			_ = hPar
		} else {
			zImg = Rp / 2
			pup = Pupil{Z: zPu, Radius: rPu}
			if !tc.field.IsInfinite() || tc.field.TanY != 0 || tc.field.TanX != 0 {
				off, err := tc.sys.ChiefOffset(tc.field, pup.Z)
				if err != nil {
					t.Fatal(err)
				}
				pup.Offset = off
			}
			pupilRadius = rPu
		}
		sd, err := tc.sys.SeidelCoefficients(Wavelengthd, tc.field, pupilRadius)
		if err != nil {
			t.Fatalf("%s: Seidel: %v", tc.name, err)
		}
		w := sd.Wavefront()
		g, err := tc.sys.OPDGrid(tc.field, pup, Rings(6, 12), zImg, Wavelengthd, OPDOptions{})
		if err != nil {
			t.Fatalf("%s: OPD: %v", tc.name, err)
		}
		c := fitBasis(g.Pupils, g.OPD, basisW)
		c40, c31s, c22, c20 := c[8], c[7], c[4], c[3]

		if math.Abs(w.W040) > 1e-9 {
			checkRel(t, tc.name+" W040 (ρ⁴)", c40, w.W040, 0.02)
		} else {
			// A parabola has no spherical aberration at any field: the fitted
			// ρ⁴ term must be negligible against the other aberrations of the
			// same map (what remains is higher-order contamination).
			scale := math.Max(math.Abs(w.W131), math.Abs(w.W222))
			if math.Abs(c40) > 0.02*scale {
				t.Errorf("%s: W040 = 0 but the fitted ρ⁴ coefficient is %g (scale %g)", tc.name, c40, scale)
			}
		}
		checkRel(t, tc.name+" W131 (ρ³sinθ)", c31s, w.W131, 0.03)
		// The astigmatism term of a y-field is W222·ρ²sin²θ, i.e. the fitted
		// ρ²cos2θ coefficient is −W222/2.
		checkRel(t, tc.name+" W222 (via c22)", c22, -w.W222/2, 0.03)
		// The rotationally symmetric ρ² term is W222/2 + W220.
		checkRel(t, tc.name+" W220 (via c20)", c20, w.W222/2+w.W220, 0.05)
	}
}

// TestSeidelDistortionVersusRealRays compares the third-order distortion of
// the sums with the distortion of the real (aimed) chief ray for a lens whose
// stop is displaced in front of it — the configuration that produces
// distortion.
func TestSeidelDistortionVersusRealRays(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	_ = n
	R := 150 * (n - 1)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(6)},
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{30, 6}, []Medium{Air{}, Air{}, N_BK7, Air{}})
	sys.StopIndex = 0
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	for _, tan := range []float64{0.03, 0.05, 0.07} {
		f := FieldTan(0, tan)
		sd, err := sys.SeidelCoefficients(Wavelengthd, f, par.EntrancePupilRadius)
		if err != nil {
			t.Fatal(err)
		}
		real, err := sys.Distortion(f, par.ImageZ, Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		pred := sd.RelativeDistortion()
		if math.Abs(real-pred) > 0.10*math.Abs(real) {
			t.Errorf("tan=%.2f: real distortion %.6g, third-order prediction %.6g (%.1f%% off)",
				tan, real, pred, 100*math.Abs(real-pred)/math.Abs(real))
		}
	}
}

// TestSeidelApertureScaling checks the aperture scaling of the sums for a
// fixed system: S1 ∝ y₀⁴, S2 ∝ y₀³, S3 ∝ y₀², S4 ∝ y₀² and S5 ∝ y₀, so the
// wavefront coefficients of the sums in the normalized ρ basis
// (WavefrontCoeffs) are aperture-independent.
func TestSeidelApertureScaling(t *testing.T) {
	const Rp = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: Rp}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	f := FieldTan(0, 0.03)
	s1, err := sys.SeidelCoefficients(Wavelengthd, f, 10)
	if err != nil {
		t.Fatal(err)
	}
	s2, err := sys.SeidelCoefficients(Wavelengthd, f, 20)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "S1 ∝ y⁴", s2.S1, s1.S1*16, 1e-9)
	checkRel(t, "S2 ∝ y³", s2.S2, s1.S2*8, 1e-9)
	checkRel(t, "S3 ∝ y²", s2.S3, s1.S3*4, 1e-9)
	checkRel(t, "S4 ∝ y²", s2.S4, s1.S4*4, 1e-9)
	checkRel(t, "S5 ∝ y", s2.S5, s1.S5*2, 1e-9)
	// The wavefront coefficients are values of the polynomial at ρ = 1, so
	// they carry the same aperture powers as the sums (W040 ∝ y₀⁴, W131 ∝ y₀³,
	// W222 ∝ y₀², W311 ∝ y₀); the aperture-independent quantities are the
	// physical coefficients per unit r, i.e. these divided by the power of y₀.
	w1, w2 := s1.Wavefront(), s2.Wavefront()
	checkRel(t, "W040 ∝ y⁴", w2.W040, w1.W040*16, 1e-9)
	checkRel(t, "W131 ∝ y³", w2.W131, w1.W131*8, 1e-9)
	checkRel(t, "W222 ∝ y²", w2.W222, w1.W222*4, 1e-9)
	checkRel(t, "W311 ∝ y", w2.W311, w1.W311*2, 1e-9)
}

// TestSeidelErrors exercises the validation paths of SeidelCoefficients.
func TestSeidelErrors(t *testing.T) {
	const Rp = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: Rp}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	if _, err := sys.SeidelCoefficients(Wavelengthd, FieldTan(0, 0), 0); err == nil {
		t.Error("zero pupil radius must be rejected")
	}
	if _, err := sys.SeidelCoefficients(Wavelengthd, FieldTan(0, 0), -5); err == nil {
		t.Error("negative pupil radius must be rejected")
	}
	if _, err := sys.SeidelCoefficients(Wavelengthd, FieldHeight(10, 0, 0), 40); err == nil {
		t.Error("an object behind the first surface must be rejected")
	}
	// A finite object in front is fine.
	if _, err := sys.SeidelCoefficients(Wavelengthd, FieldHeight(-100, 0, 5), 40); err != nil {
		t.Errorf("finite object rejected: %v", err)
	}
}
