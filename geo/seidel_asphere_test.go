package geo

import (
	"fmt"
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// helpers for the asphere Seidel tests
// ---------------------------------------------------------------------------

// checkRelStrict fails when |got − want| > tol·scale. It is a strictly
// relative comparison: the shared checkRel of the package floors the scale at
// 1, which would turn the equality tests below (whose values are ~mm·1e-3 and
// must agree to 1e-12 relative) into absolute tests.
func checkRelStrict(t *testing.T, name string, got, want, scale, tol float64) {
	t.Helper()
	if math.IsNaN(got) || math.IsNaN(want) {
		t.Errorf("%s: NaN (got %v, want %v)", name, got, want)
		return
	}
	if scale <= 0 {
		scale = math.Max(math.Abs(got), math.Abs(want))
	}
	if scale <= 0 {
		return // both exactly zero
	}
	if d := math.Abs(got - want); d > tol*scale {
		t.Errorf("%s: got %.17g, want %.17g (|diff| = %.3g > %.3g·%.3g)", name, got, want, d, tol, scale)
	}
}

// seidelTestRig builds the reference system of the asphere tests: a BK7
// singlet whose front surface (surface 0, also the stop) carries the profile
// under test and whose rear surface is spherical. The lens is used off axis,
// so every Seidel sum responds to the profile: S1 through the conic/r⁴ term
// of the front surface, the field sums through its curvature.
func seidelTestRig(t *testing.T, prof Profile) *System {
	t.Helper()
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: prof, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{6}, []Medium{Air{}, N_BK7, Air{}})
}

const aspherePupil = 6.0

// asphereSums returns the Seidel solution of a rig for the shared field and
// pupil radius of the asphere tests.
func asphereSums(t *testing.T, sys *System) *Seidel {
	t.Helper()
	f := FieldTan(0, 0.05)
	sd, err := sys.SeidelCoefficients(Wavelengthd, f, aspherePupil)
	if err != nil {
		t.Fatalf("SeidelCoefficients: %v", err)
	}
	return sd
}

// checkSeidelEqual asserts that all five sums and the Petzval sum of two
// solutions agree to a relative tolerance; scale is the magnitude of the
// largest sum involved, so the comparison is genuinely relative. It also
// refuses a vacuous comparison of two identically zero solutions.
func checkSeidelEqual(t *testing.T, name string, a, b *Seidel, tol float64) {
	t.Helper()
	scale := 0.0
	for _, v := range []float64{a.S1, a.S2, a.S3, a.S4, a.S5, b.S1, b.S2, b.S3, b.S4, b.S5} {
		scale = math.Max(scale, math.Abs(v))
	}
	if scale == 0 {
		t.Fatalf("%s: both solutions are identically zero — vacuous comparison", name)
	}
	sums := [5][2]float64{{a.S1, b.S1}, {a.S2, b.S2}, {a.S3, b.S3}, {a.S4, b.S4}, {a.S5, b.S5}}
	for i, pair := range sums {
		checkRelStrict(t, fmt.Sprintf("%s: S%d", name, i+1), pair[0], pair[1], scale, tol)
	}
	if a.Petzval != 0 || b.Petzval != 0 {
		checkRelStrict(t, name+": Petzval", a.Petzval, b.Petzval, math.Abs(a.Petzval), tol)
	}
}

// fitSag4 interpolates the sag of p at r, r/2, r/4 and r/8 with
// c2·r² + c4·r⁴ + c6·r⁶ + c8·r⁸ (Vandermonde system in t = r², solved with
// partial pivoting) and returns c2 and c4. The interpolation is exact for the
// terms through r⁸, so the extracted c2/c4 are the Taylor coefficients of the
// sag up to the contamination of the next neglected order (r¹⁰), which is
// far below the tolerances used here for small r.
func fitSag4(p Profile, r float64) (c2, c4 float64) {
	var a [4][5]float64
	for i, rr := range [4]float64{r, r / 2, r / 4, r / 8} {
		t2 := rr * rr
		a[i] = [5]float64{t2, t2 * t2, t2 * t2 * t2, t2 * t2 * t2 * t2, p.Sag(rr)}
	}
	for col := 0; col < 4; col++ {
		piv := col
		for row := col + 1; row < 4; row++ {
			if math.Abs(a[row][col]) > math.Abs(a[piv][col]) {
				piv = row
			}
		}
		a[col], a[piv] = a[piv], a[col]
		for row := col + 1; row < 4; row++ {
			f := a[row][col] / a[col][col]
			for c := col; c < 5; c++ {
				a[row][c] -= f * a[col][c]
			}
		}
	}
	var x [4]float64
	for row := 3; row >= 0; row-- {
		s := a[row][4]
		for c := row + 1; c < 4; c++ {
			s -= a[row][c] * x[c]
		}
		x[row] = s / a[row][row]
	}
	return x[0], x[1]
}

// ---------------------------------------------------------------------------
// asphereParams mapping
// ---------------------------------------------------------------------------

// TestSeidelAsphereParamsMapping pins the (K, A4) mapping of the two polynomial
// asphere types directly: the base conic constant, the r⁴ coefficient (terms
// of equal power add), and the fact that neither the higher even powers nor
// the odd powers enter the third-order pair. An r² term is the one exception:
// it shifts the vertex curvature, and the returned A4 then carries the exact
// base-shift term of polyAsphereParams (see TestSeidelAsphereParamsFromSag
// and TestSeidelPolyAsphereCurvatureShift).
func TestSeidelAsphereParamsMapping(t *testing.T) {
	const R = -200.0
	if K, A4 := asphereParams(PolyAsphere{R: R, K: -0.6, Terms: []PolyTerm{{Power: 4, Coeff: 1e-7}}}); K != -0.6 || A4 != 1e-7 {
		t.Errorf("PolyAsphere single r⁴: got (%g %g), want (-0.6 1e-07)", K, A4)
	}
	// Equal powers add; r⁶, r⁸ and the odd powers stay out of the pair.
	poly := PolyAsphere{R: R, K: 0.5, Terms: []PolyTerm{
		{Power: 8, Coeff: 3e-9}, {Power: 4, Coeff: 1e-7},
		{Power: 3, Coeff: 1e-9}, {Power: 4, Coeff: 2e-7}, {Power: 6, Coeff: 4e-9},
	}}
	if K, A4 := asphereParams(poly); K != 0.5 || A4 != 1e-7+2e-7 {
		t.Errorf("PolyAsphere power-4 sum: got (%g %g), want (0.5 %g)", K, A4, 1e-7+2e-7)
	}
	// No r⁴ term and no r² term: A4 = 0 exactly, whatever the other powers do.
	for _, terms := range [][]PolyTerm{
		{{Power: 6, Coeff: 1e-9}, {Power: 10, Coeff: 1e-12}},
		{{Power: 1, Coeff: 1e-3}, {Power: 3, Coeff: 1e-5}},
		nil,
	} {
		if K, A4 := asphereParams(PolyAsphere{R: R, K: -1.2, Terms: terms}); K != -1.2 || A4 != 0 {
			t.Errorf("PolyAsphere %v: got (%g %g), want (-1.2 0)", terms, K, A4)
		}
	}
	// An r² term shifts the curvature by 2a₂; A4 carries the base-shift term
	// (1+K)(c₀³ − c³)/8 exactly.
	kk := -1.2
	only2 := PolyAsphere{R: R, K: kk, Terms: []PolyTerm{{Power: 2, Coeff: 1e-4}}}
	c0 := baseCurvature(R)
	c := only2.Curvature()
	wantShift := (1 + kk) * (c0*c0*c0 - c*c*c) / 8
	if K, A4 := asphereParams(only2); K != kk || A4 != wantShift {
		t.Errorf("PolyAsphere r² only: got (%g %g), want (%g %g)", K, A4, kk, wantShift)
	}
	// ... which vanishes exactly for K = −1 (the parabolic base has no r⁴
	// term to re-express) and for a parabolic shift-free surface.
	if _, A4 := asphereParams(PolyAsphere{R: R, K: -1, Terms: []PolyTerm{{Power: 2, Coeff: 1e-4}}}); A4 != 0 {
		t.Errorf("PolyAsphere r² only, K = −1: A4 = %g, want 0", A4)
	}
	// Pointer forms behave like the values.
	pp := PolyAsphere{R: R, K: 0.25, Terms: []PolyTerm{{Power: 4, Coeff: -5e-8}}}
	if K, A4 := asphereParams(&pp); K != 0.25 || A4 != -5e-8 {
		t.Errorf("*PolyAsphere: got (%g %g), want (0.25 -5e-08)", K, A4)
	}
	var nilPoly *PolyAsphere
	if K, A4 := asphereParams(nilPoly); K != 0 || A4 != 0 {
		t.Errorf("nil *PolyAsphere: got (%g %g), want (0 0)", K, A4)
	}
	// QAsphere without coefficients: the base conic constant, no r⁴ term.
	qa := QAsphere{Form: QCon, R: R, K: -0.7, NormR: 20}
	if K, A4 := asphereParams(qa); K != -0.7 || A4 != 0 {
		t.Errorf("QAsphere zero coeffs: got (%g %g), want (-0.7 0)", K, A4)
	}
	qb := QAsphere{Form: QBfs, R: R, K: 0, NormR: 20, Coeffs: []float64{0, 0, 0}}
	if K, A4 := asphereParams(&qb); K != 0 || A4 != 0 {
		t.Errorf("*QAsphere zero coeffs: got (%g %g), want (0 0)", K, A4)
	}
	var nilQ *QAsphere
	if K, A4 := asphereParams(nilQ); K != 0 || A4 != 0 {
		t.Errorf("nil *QAsphere: got (%g %g), want (0 0)", K, A4)
	}
}

// TestQconBasisAxisValues checks the identity behind the Qcon r⁴ coefficient,
// Q_m^con(0) = P_m^(0,4)(−1) = (−1)^m·C(m+4, 4) (the standard value of a Jacobi
// polynomial at the lower end point), against the recurrence that evaluates
// the basis, up to an order well beyond any practical asphere.
func TestQconBasisAxisValues(t *testing.T) {
	const m = 14
	z := make([]float64, m)
	dz := make([]float64, m)
	qconBasis(0, z, dz)
	for i := range z {
		fm := float64(i)
		want := math.Pow(-1, fm) * (fm + 1) * (fm + 2) * (fm + 3) * (fm + 4) / 24
		checkRelStrict(t, fmt.Sprintf("Q_%d^con(0)", i), z[i], want, math.Abs(want), 1e-12)
	}
}

// TestQbfsBasisAxisValues pins the axis values of the Qbfs basis against the
// exact values of the recurrence constants of Forbes (2010) Eq. (A.14)–(A.16),
// which are Q_m(0) = 1, 13/√19, 29√190/95, 207√5090/2545 and
// Q_m′(0) = 0, −16/√19, −20√190/19, −252√5090/509 for m = 0…3. They were
// derived symbolically from the recurrence (the derivation is reproduced in
// the task notes) and are what qAsphereParams folds into its closed form.
func TestQbfsBasisAxisValues(t *testing.T) {
	coeffs := []float64{1}
	sum, dsum := qbfsSum(coeffs, 0)
	checkClose(t, "Q_0(0)", sum, 1, 1e-15)
	checkClose(t, "Q_0′(0)", dsum, 0, 1e-15)
	coeffs = []float64{0, 1}
	sum, dsum = qbfsSum(coeffs, 0)
	checkRelStrict(t, "Q_1(0)", sum, 13/math.Sqrt(19), 13/math.Sqrt(19), 1e-14)
	checkRelStrict(t, "Q_1′(0)", dsum, -16/math.Sqrt(19), 16/math.Sqrt(19), 1e-14)
	coeffs = []float64{0, 0, 1}
	sum, dsum = qbfsSum(coeffs, 0)
	w2 := 29 * math.Sqrt(190) / 95
	checkRelStrict(t, "Q_2(0)", sum, w2, w2, 1e-14)
	d2 := -20 * math.Sqrt(190) / 19
	checkRelStrict(t, "Q_2′(0)", dsum, d2, math.Abs(d2), 1e-14)
	coeffs = []float64{0, 0, 0, 1}
	sum, dsum = qbfsSum(coeffs, 0)
	w3 := 207 * math.Sqrt(5090) / 2545
	checkRelStrict(t, "Q_3(0)", sum, w3, w3, 1e-14)
	d3 := -252 * math.Sqrt(5090) / 509
	checkRelStrict(t, "Q_3′(0)", dsum, d3, math.Abs(d3), 1e-14)
}

// TestSeidelAsphereParamsFromSag validates the whole (K, A4) mapping
// numerically and representation-independently: the sag of every profile is
// interpolated at four small radii, which yields its curvature c and its r⁴
// coefficient E, and the pair returned by asphereParams must satisfy the
// invariant of the S1 formula, K·c³ + 8·A4 = 8E − c³, for the conic-based
// profiles (PolyAsphere, QAsphere, EvenAsphere, Conic).
func TestSeidelAsphereParamsFromSag(t *testing.T) {
	const (
		rFit = 0.8 // mm, small against every aperture below
		R    = 77.5
	)
	cases := []struct {
		name string
		p    Profile
	}{
		{"Conic", Conic{R: R, K: -0.6}},
		{"EvenAsphere", EvenAsphere{R: R, K: -0.6, Coeffs: []float64{-3e-7, 2e-9}}},
		{"PolyAsphere r⁴/r⁶", PolyAsphere{R: R, K: -0.6, Terms: []PolyTerm{{Power: 4, Coeff: -3e-7}, {Power: 6, Coeff: 2e-9}}}},
		{"PolyAsphere r²/r⁴", PolyAsphere{R: R, K: -0.6, Terms: []PolyTerm{{Power: 2, Coeff: 1e-4}, {Power: 4, Coeff: -3e-7}}}},
		{"PolyAsphere r² only", PolyAsphere{R: R, K: -0.6, Terms: []PolyTerm{{Power: 2, Coeff: 1e-4}}}},
		{"QCon", QAsphere{Form: QCon, R: R, K: -0.4, NormR: 12, Coeffs: []float64{1e-4, -2e-5, 5e-6, -1e-6}}},
		{"QCon flat base", QAsphere{Form: QCon, R: math.Inf(1), K: 0, NormR: 12, Coeffs: []float64{1e-4, -2e-5}}},
		{"QBfs", QAsphere{Form: QBfs, R: R, K: 0, NormR: 12, Coeffs: []float64{1e-4, -2e-5, 5e-6, -1e-6}}},
	}
	for _, tc := range cases {
		K, A4 := asphereParams(tc.p)
		c2, c4 := fitSag4(tc.p, rFit)
		// The interpolated r² coefficient is the vertex curvature of the sag.
		sc := math.Abs(tc.p.Curvature())
		if sc == 0 {
			sc = 1
		}
		checkRelStrict(t, tc.name+": fitted curvature", 2*c2, tc.p.Curvature(), sc, 1e-8)
		cFit := 2 * c2
		lhs := K*cFit*cFit*cFit + 8*A4
		rhs := 8*c4 - cFit*cFit*cFit
		scale := math.Abs(K*cFit*cFit*cFit) + math.Abs(8*A4) + math.Abs(8*c4) + math.Abs(cFit*cFit*cFit)
		checkRelStrict(t, tc.name+": r⁴ invariant", lhs, rhs, scale, 1e-6)
	}
}

// ---------------------------------------------------------------------------
// acceptance (a): PolyAsphere with a single r⁴ term vs. an equivalent
// EvenAsphere
// ---------------------------------------------------------------------------

// TestSeidelPolyAsphereMatchesEvenAsphere checks that a PolyAsphere whose
// terms are the r⁴/r⁶ coefficients of an EvenAsphere gives the identical
// third-order solution. The contrast systems (the bare conic and an asphere
// with only the r⁶ term) show that the comparison has teeth: dropping the r⁴
// term moves W040 far beyond the 1e-12 tolerance, while dropping the r⁶ term
// changes nothing at all.
func TestSeidelPolyAsphereMatchesEvenAsphere(t *testing.T) {
	const (
		K  = -0.7
		a4 = -1.5e-7
		a6 = 2e-10
	)
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	even := asphereSums(t, seidelTestRig(t, EvenAsphere{R: R, K: K, Coeffs: []float64{a4, a6}}))
	poly := asphereSums(t, seidelTestRig(t, PolyAsphere{R: R, K: K, Terms: []PolyTerm{{Power: 4, Coeff: a4}, {Power: 6, Coeff: a6}}}))
	checkSeidelEqual(t, "PolyAsphere vs EvenAsphere", even, poly, 1e-12)
	// The wavefront coefficients of the task's acceptance criterion.
	wE, wP := even.Wavefront(), poly.Wavefront()
	checkRelStrict(t, "W040", wP.W040, wE.W040, math.Abs(wE.W040), 1e-12)
	checkRelStrict(t, "W131", wP.W131, wE.W131, math.Abs(wE.W131), 1e-12)
	if wE.W040 == 0 || wE.W131 == 0 {
		t.Fatalf("the rig must have spherical aberration and coma (W040=%g W131=%g)", wE.W040, wE.W131)
	}
	// Teeth: without the r⁴ term the sums differ far beyond the tolerance.
	cone := asphereSums(t, seidelTestRig(t, Conic{R: R, K: K}))
	dW := math.Abs(wE.W040 - cone.Wavefront().W040)
	if dW < 1e-3*math.Abs(wE.W040) {
		t.Errorf("the r⁴ term of the asphere is not observable in W040 (Δ = %g, W040 = %g)", dW, wE.W040)
	}
	// The r⁶ term alone must not move the sums at all (it is a fifth-order
	// contribution): an asphere with Coeffs = [0, a6] is the conic for
	// third-order purposes.
	r6only := asphereSums(t, seidelTestRig(t, EvenAsphere{R: R, K: K, Coeffs: []float64{0, a6}}))
	checkSeidelEqual(t, "r⁶ term is outside the third order", cone, r6only, 1e-15)
}

// TestSeidelPolyAsphereCurvatureShift pins the handling of a power-2 term,
// which shifts the vertex curvature of a PolyAsphere (Curvature = c₀ + 2a₂)
// and therefore the osculating sphere of the third-order sums. In third order
// the surface is identical to the conic of the shifted curvature c and the
// equivalent conic constant K' = 8E/c³ − 1 (E = (1+K)c₀³/8 + a₄); the A4 of
// polyAsphereParams carries exactly the base-shift correction, so the two
// solutions agree. The difference from the literal reading of A4 — an
// EvenAsphere on the shifted radius with the bare polynomial coefficient — is
// the predicted term (n'−n)(1+K)(c₀³ − c³)y₀⁴.
func TestSeidelPolyAsphereCurvatureShift(t *testing.T) {
	const (
		K  = -0.4
		a2 = 1e-3 // r² term: shifts the curvature by 2e-3 mm⁻¹
		a4 = -1e-6
	)
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	c0 := 1 / R
	poly := PolyAsphere{R: R, K: K, Terms: []PolyTerm{{Power: 2, Coeff: a2}, {Power: 4, Coeff: a4}}}
	c := poly.Curvature()
	if c == c0 {
		t.Fatal("the r² term must shift the curvature")
	}
	E := (1+K)*c0*c0*c0/8 + a4
	keq := 8*E/(c*c*c) - 1
	equiv := Conic{R: 1 / c, K: keq}

	sdPoly := asphereSums(t, seidelTestRig(t, poly))
	sdConic := asphereSums(t, seidelTestRig(t, equiv))
	checkSeidelEqual(t, "PolyAsphere(r²) vs equivalent conic", sdPoly, sdConic, 1e-12)

	// The equivalence itself: the two sags have the same Taylor coefficients
	// through r⁴ (they are different surfaces at r⁶ and beyond).
	const rFit = 0.5
	c2P, c4P := fitSag4(poly, rFit)
	c2C, c4C := fitSag4(equiv, rFit)
	checkRelStrict(t, "curvature through r⁴", c2P, c2C, math.Abs(c2C), 1e-8)
	checkRelStrict(t, "r⁴ coefficient", 8*c4P, 8*c4C, 8*math.Abs(c4C), 1e-7)

	// The correction matters: the literal A4 = a₄ on the shifted radius is a
	// different third-order surface by exactly (n'−n)(1+K)(c₀³ − c³)y₀⁴.
	literal := asphereSums(t, seidelTestRig(t, EvenAsphere{R: 1 / c, K: K, Coeffs: []float64{a4}}))
	nA := N_BK7.N(Wavelengthd)
	nB := Air{}.N(Wavelengthd)
	want := -(nA - nB) * (1 + K) * (c0*c0*c0 - c*c*c) * math.Pow(aspherePupil, 4)
	got := literal.S1 - sdPoly.S1
	checkRelStrict(t, "base-shift term of S1", got, want, math.Abs(want), 1e-9)
	if math.Abs(got) < 1e-3*math.Abs(sdPoly.S1) {
		t.Errorf("the base-shift correction is not observable: ΔS1 = %g vs S1 = %g", got, sdPoly.S1)
	}
}

// ---------------------------------------------------------------------------
// acceptance (b): a QAsphere without coefficients is its base surface
// ---------------------------------------------------------------------------

// TestSeidelQAsphereZeroCoeffsEqualBase checks that a QAsphere with an empty,
// nil or all-zero coefficient list gives the Seidel solution of its base
// surface: the Conic(R, K) for QCon, the Sphere(R) for QBfs (K must be 0
// there). K, c and A4 all have to reduce correctly for this to hold.
func TestSeidelQAsphereZeroCoeffsEqualBase(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	for _, K := range []float64{0, -0.5, 0.8} {
		for _, coeffs := range [][]float64{nil, {}, {0, 0, 0, 0}} {
			qa := asphereSums(t, seidelTestRig(t, QAsphere{Form: QCon, R: R, K: K, NormR: 15, Coeffs: coeffs}))
			base := asphereSums(t, seidelTestRig(t, Conic{R: R, K: K}))
			checkSeidelEqual(t, fmt.Sprintf("QCon K=%.1f coeffs=%v vs Conic", K, coeffs), qa, base, 1e-12)
		}
	}
	for _, coeffs := range [][]float64{nil, {}, {0, 0, 0}} {
		qa := asphereSums(t, seidelTestRig(t, QAsphere{Form: QBfs, R: R, K: 0, NormR: 15, Coeffs: coeffs}))
		base := asphereSums(t, seidelTestRig(t, Sphere{R: R}))
		checkSeidelEqual(t, fmt.Sprintf("QBfs coeffs=%v vs Sphere", coeffs), qa, base, 1e-12)
	}
	// A positive normalisation radius is only required while coefficients are
	// active: the bare base profile is still exact.
	bare := asphereSums(t, seidelTestRig(t, QAsphere{Form: QCon, R: R, K: -0.3, NormR: 0}))
	base := asphereSums(t, seidelTestRig(t, Conic{R: R, K: -0.3}))
	checkSeidelEqual(t, "QCon NormR=0 without coefficients", bare, base, 1e-12)
}

// ---------------------------------------------------------------------------
// acceptance (c): self-consistency of the Qcon/QBfs r⁴ coefficients
// ---------------------------------------------------------------------------

// TestSeidelQAsphereQconEquivalentPoly verifies the Qcon r⁴ coefficient the
// way the acceptance criterion prescribes: the QAsphere is compared with the
// PolyAsphere that has the same base conic (R, K) and a single r⁴ term whose
// coefficient is the A4 returned by asphereParams.
//
//   - With only the m = 0 coefficient the departure is exactly the single term
//     u⁴·s₀·Q₀(x) = s₀u⁴, so the two profiles must agree pointwise at every
//     radius and their Seidel solutions must be identical.
//   - With several coefficients the single-r⁴ PolyAsphere can only match the
//     sag through r⁴ order: the signed pointwise difference qa.Sag − poly.Sag
//     must equal the leading neglected term of the departure, u⁶·S′(0), and
//     must fall off like r⁶ (a factor 2⁶ = 64 when the radius is halved),
//     while the Seidel solutions stay identical because both profiles are
//     described by the same (K, A4).
func TestSeidelQAsphereQconEquivalentPoly(t *testing.T) {
	const (
		K   = -0.5
		rho = 12.0
	)
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	s0 := 3e-4

	one := QAsphere{Form: QCon, R: R, K: K, NormR: rho, Coeffs: []float64{s0}}
	polyOne := PolyAsphere{R: R, K: K, Terms: []PolyTerm{{Power: 4, Coeff: s0 / (rho * rho * rho * rho)}}}
	if Kgot, A4 := asphereParams(one); Kgot != K || A4 != s0/(rho*rho*rho*rho) {
		t.Fatalf("asphereParams(QCon s₀): got (%g %g), want (%g %g)", Kgot, A4, K, s0/(rho*rho*rho*rho))
	}
	for _, r := range []float64{0.3, 1, 3, rho / 2, 0.9 * rho} {
		got, want := one.Sag(r), polyOne.Sag(r)
		checkRelStrict(t, fmt.Sprintf("pointwise sag at r=%.2f", r), got, want, math.Abs(want), 1e-13)
	}
	checkSeidelEqual(t, "QCon single coefficient vs PolyAsphere",
		asphereSums(t, seidelTestRig(t, one)), asphereSums(t, seidelTestRig(t, polyOne)), 1e-12)

	// Several coefficients: A4 must be the r⁴ limit of the whole departure.
	coeffs := []float64{1e-4, -2e-6, 5e-7, -1e-7}
	multi := QAsphere{Form: QCon, R: R, K: K, NormR: rho, Coeffs: coeffs}
	Kgot, A4 := asphereParams(multi)
	// Independent closed form: S(0) = Σ (−1)^m C(m+4, 4) s_m.
	wantA4 := 0.0
	for m, s := range coeffs {
		fm := float64(m)
		wantA4 += math.Pow(-1, fm) * (fm + 1) * (fm + 2) * (fm + 3) * (fm + 4) / 24 * s
	}
	wantA4 /= rho * rho * rho * rho
	checkRelStrict(t, "QCon A4 closed form", A4, wantA4, math.Abs(wantA4), 1e-12)
	if Kgot != K {
		t.Errorf("QCon K = %g, want %g", Kgot, K)
	}
	polyMulti := PolyAsphere{R: R, K: K, Terms: []PolyTerm{{Power: 4, Coeff: A4}}}
	// The signed difference is the leading neglected term of the departure.
	_, dS := qconSum(coeffs, 0) // S′(0)
	for _, r := range []float64{rho / 8, rho / 16} {
		got := multi.Sag(r) - polyMulti.Sag(r)
		want := math.Pow(r/rho, 6) * dS
		checkRelStrict(t, fmt.Sprintf("residual is the r⁶ term at r=ρ/%d", int(rho/r)), got, want, math.Abs(want), 0.15)
	}
	d1 := multi.Sag(rho/8) - polyMulti.Sag(rho/8)
	d2 := multi.Sag(rho/16) - polyMulti.Sag(rho/16)
	if ratio := d1 / d2; ratio < 32 || ratio > 128 {
		t.Errorf("sag gap falls by %g×(want ≈ 2⁶ = 64 for an r⁶ residual)", ratio)
	}
	checkSeidelEqual(t, "QCon multi-coefficient vs PolyAsphere",
		asphereSums(t, seidelTestRig(t, multi)), asphereSums(t, seidelTestRig(t, polyMulti)), 1e-12)
}

// TestSeidelQAsphereQbfsEquivalentPoly is the self-consistency check of the
// QBfs r⁴ coefficient. The departure u²(1−u²)S(x)/ψ has an r² term,
// 2·S(0)/ρ², which QAsphere.Curvature adds to the base curvature; the
// equivalent single-r⁴ PolyAsphere therefore has to be built on the shifted
// radius R′ = 1/c (same K = 0). Its r⁴ coefficient is the A4 that
// asphereParams returns, and the profile pair must then agree pointwise
// through r⁴ order (the residual falls like r⁶, a factor 64 per halving) and
// in every Seidel sum.
func TestSeidelQAsphereQbfsEquivalentPoly(t *testing.T) {
	const (
		R   = 120.0
		rho = 10.0
	)
	coeffs := []float64{2e-4, -3e-5, 4e-6, -5e-7}
	qa := QAsphere{Form: QBfs, R: R, K: 0, NormR: rho, Coeffs: coeffs}
	K, A4 := asphereParams(qa)
	if K != 0 {
		t.Errorf("QBfs K = %g, want 0", K)
	}
	c := qa.Curvature()
	c0 := 1 / R
	if c == c0 {
		t.Fatal("the QBfs departure must shift the curvature")
	}
	// Independent closed form of A4, from the axis values of the basis (see
	// TestQbfsBasisAxisValues for their exact values):
	//   S(0) = a₀ + 13/√19·a₁ + 29√190/95·a₂ + 207√5090/2545·a₃
	//   S′(0) = −16/√19·a₁ − 20√190/19·a₂ − 252√5090/509·a₃
	//   c = c₀ + 2S(0)/ρ²,  E = c₀³/8 + (S′(0)−S(0))/ρ⁴ + c₀²S(0)/(2ρ²),
	//   A4 = E − c³/8.
	s0 := coeffs[0] + 13/math.Sqrt(19)*coeffs[1] + 29*math.Sqrt(190)/95*coeffs[2] + 207*math.Sqrt(5090)/2545*coeffs[3]
	ds0 := -16/math.Sqrt(19)*coeffs[1] - 20*math.Sqrt(190)/19*coeffs[2] - 252*math.Sqrt(5090)/509*coeffs[3]
	cExp := c0 + 2*s0/(rho*rho)
	eExp := c0*c0*c0/8 + (ds0-s0)/(rho*rho*rho*rho) + c0*c0*s0/(2*rho*rho)
	a4Exp := eExp - cExp*cExp*cExp/8
	checkRelStrict(t, "QBfs A4 closed form", A4, a4Exp, math.Abs(a4Exp), 1e-12)
	checkRelStrict(t, "QBfs curvature", c, cExp, math.Abs(cExp), 1e-13)

	poly := PolyAsphere{R: 1 / c, K: 0, Terms: []PolyTerm{{Power: 4, Coeff: A4}}}
	// Pointwise agreement through r⁴: the residual falls like r⁶ and stays
	// small against the r⁴ term of the departure.
	ref := math.Abs(A4) * math.Pow(rho/8, 4)
	d1 := qa.Sag(rho/8) - poly.Sag(rho/8)
	d2 := qa.Sag(rho/16) - poly.Sag(rho/16)
	if math.Abs(d1) > 0.05*ref {
		t.Errorf("QBfs vs single-r⁴ PolyAsphere: sag gap %g at r=ρ/8 is large against the r⁴ term %g", d1, ref)
	}
	if ratio := d1 / d2; ratio < 32 || ratio > 128 {
		t.Errorf("sag gap falls by %g×(want ≈ 2⁶ = 64 for an r⁶ residual)", ratio)
	}
	checkSeidelEqual(t, "QBfs vs single-r⁴ PolyAsphere",
		asphereSums(t, seidelTestRig(t, qa)), asphereSums(t, seidelTestRig(t, poly)), 1e-12)
}

// ---------------------------------------------------------------------------
// the new profiles against real ray tracing
// ---------------------------------------------------------------------------

// TestSeidelQAsphereVersusOPD extends the central validation of the aberration
// layer (TestSeidelWavefrontVersusOPD) to the Q-type aspheres: a concave
// mirror carries a Qcon departure chosen to change its ρ⁴ wavefront term by
// a known amount, and the third-order prediction of the sums — which now
// includes A4 = S(0)/ρ⁴ — must match the exact OPD map of the real rays. The
// sphere-only prediction is off by tens of percent, so the check is
// meaningful.
func TestSeidelQAsphereVersusOPD(t *testing.T) {
	const (
		Rp   = -200.0
		rPu  = 40.0
		zPu  = -50.0
		norm = 40.0
	)
	// Reference: the bare sphere mirror.
	sphere := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: Rp}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	sdSphere, err := sphere.SeidelCoefficients(Wavelengthd, FieldTan(0, 0), rPu)
	if err != nil {
		t.Fatal(err)
	}
	w0 := sdSphere.Wavefront().W040
	if w0 == 0 {
		t.Fatal("the mirror must have spherical aberration")
	}
	// Add a departure that changes W040 by +35 %: ΔW040 = −(n′−n)·A4·y⁴ with
	// n′−n = −2 for a mirror, so A4 = ΔW040/(2y⁴).
	const frac = 0.35
	a4 := frac * w0 / (2 * math.Pow(rPu, 4))
	s0 := a4 * math.Pow(norm, 4)
	qa := QAsphere{Form: QCon, R: Rp, K: 0, NormR: norm, Coeffs: []float64{s0, 0.2 * s0 / norm}}
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: qa, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	for _, tc := range []struct {
		name  string
		field Field
	}{{"on axis", FieldTan(0, 0)}, {"off axis", FieldTan(0, 0.03)}} {
		sd, err := sys.SeidelCoefficients(Wavelengthd, tc.field, rPu)
		if err != nil {
			t.Fatalf("%s: %v", tc.name, err)
		}
		w := sd.Wavefront()
		pup := Pupil{Z: zPu, Radius: rPu}
		if !tc.field.IsInfinite() || tc.field.TanY != 0 || tc.field.TanX != 0 {
			off, err := sys.ChiefOffset(tc.field, pup.Z)
			if err != nil {
				t.Fatal(err)
			}
			pup.Offset = off
		}
		g, err := sys.OPDGrid(tc.field, pup, Rings(6, 12), Rp/2, Wavelengthd, OPDOptions{})
		if err != nil {
			t.Fatalf("%s: OPD: %v", tc.name, err)
		}
		c := fitBasis(g.Pupils, g.OPD, basisW)
		c40 := c[8]
		// The sums must have moved the ρ⁴ term by exactly the intended
		// fraction (a departure is a pure ρ⁴ term, so this holds at every
		// field) ...
		checkRelStrict(t, tc.name+" sums W040", w.W040, (1+frac)*w0, math.Abs(w0), 0.02)
		// ... the OPD ρ⁴ term must be clearly away from the sphere-only
		// prediction ...
		if math.Abs(c40-w0) < 0.1*math.Abs(c40) {
			t.Errorf("%s: the OPD ρ⁴ term %g is too close to the sphere-only value %g — the departure is not observable", tc.name, c40, w0)
		}
		// ... and the OPD fit must agree with the sums.
		checkRelStrict(t, tc.name+" W040 versus OPD ρ⁴", c40, w.W040, math.Abs(w.W040), 0.03)
	}
}

// TestSeidelPolyAsphereVersusOPD repeats the real-ray validation for a
// PolyAsphere mirror with an r² term (the curvature-shift case): the third-
// order sums of the shifted-curvature surface must match the OPD fit of the
// real rays, which is the physical statement behind the base-shift correction
// of polyAsphereParams.
func TestSeidelPolyAsphereVersusOPD(t *testing.T) {
	const (
		Rp  = -200.0
		rPu = 40.0
	)
	a2 := 2e-4
	a4 := -1e-8
	prof := PolyAsphere{R: Rp, K: 0, Terms: []PolyTerm{{Power: 2, Coeff: a2}, {Power: 4, Coeff: a4}}}
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: prof, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	field := FieldTan(0, 0)
	sd, err := sys.SeidelCoefficients(Wavelengthd, field, rPu)
	if err != nil {
		t.Fatal(err)
	}
	w := sd.Wavefront()
	if w.W040 == 0 {
		t.Fatal("the mirror must have spherical aberration")
	}
	// The paraxial focus of the shifted-curvature mirror: the marginal ray of
	// the Seidel solution crosses the axis there (the mirror vertex sits at
	// z = 0, the folded axis runs towards negative z).
	zImg := 0.0 - sd.MarginalHeight/sd.MarginalSlope
	pup := Pupil{Z: -50, Radius: rPu}
	g, err := sys.OPDGrid(field, pup, Rings(6, 12), zImg, Wavelengthd, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	c := fitBasis(g.Pupils, g.OPD, basisW)
	checkRelStrict(t, "PolyAsphere(r²) W040 versus OPD ρ⁴", c[8], w.W040, math.Abs(w.W040), 0.03)
	// Teeth: the sphere-only surface has a very different ρ⁴ term.
	sphere := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: Rp}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	sdS, err := sphere.SeidelCoefficients(Wavelengthd, field, rPu)
	if err != nil {
		t.Fatal(err)
	}
	ws := sdS.Wavefront()
	if math.Abs(ws.W040-w.W040) < 0.1*math.Abs(w.W040) {
		t.Errorf("the r² term of the profile is not observable: W040 = %g vs sphere %g", w.W040, ws.W040)
	}
}
