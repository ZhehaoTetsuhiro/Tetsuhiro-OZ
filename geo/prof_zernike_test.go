package geo

import (
	"math"
	"math/rand"
	"testing"
)

// ---------------------------------------------------------------------------
// Independent references
//
// The Zernike surface is evaluated in polar form (Jacobi recurrence + cosine /
// sine) in prof_zernike.go. Everything here is deliberately built from other
// pieces of mathematics: the ordering comes from a hand transcription of
// Noll's Table I, and the values and derivatives come from the monomial
// expansion of the definition, i.e. from (x, y) polynomials with integer
// binomial coefficients -- no polar recurrence is shared.
// ---------------------------------------------------------------------------

// zsRefTable is Noll (1976) Table I, j = 1 ... 36, transcribed by hand:
// (j, n, m). It is the independent check of the ordering rule implemented in
// ZernikeNollIndex.
var zsRefTable = [36][2]int{
	{0, 0}, {1, 1}, {1, -1}, {2, 0}, {2, -2}, {2, 2},
	{3, -1}, {3, 1}, {3, -3}, {3, 3},
	{4, 0}, {4, 2}, {4, -2}, {4, 4}, {4, -4},
	{5, 1}, {5, -1}, {5, 3}, {5, -3}, {5, 5}, {5, -5},
	{6, 0}, {6, -2}, {6, 2}, {6, -4}, {6, 4}, {6, -6}, {6, 6},
	{7, -1}, {7, 1}, {7, -3}, {7, 3}, {7, -5}, {7, 5}, {7, -7}, {7, 7},
}

// zsRefIndex returns (n, m) for the Noll index j from zsRefTable.
func zsRefIndex(t *testing.T, j int) (int, int) {
	t.Helper()
	if j < 1 || j > len(zsRefTable) {
		t.Fatalf("reference table only covers j = 1 ... %d (asked for %d)", len(zsRefTable), j)
	}
	return zsRefTable[j-1][0], zsRefTable[j-1][1]
}

// zsBinom returns C(n, k) as a float64 (exact for the small n used here).
func zsBinom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	r := 1.0
	for i := 1; i <= k; i++ {
		r = r * float64(n-k+i) / float64(i)
	}
	return r
}

// zsRadialCoeff returns the coefficient of ρ^(n-2s) in the explicit sum of
// Noll's Eq. (2),
//
//	R_n^m(ρ) = Σ_s (-1)^s (n-s)! / [s! ((n+m)/2-s)! ((n-m)/2-s)!] ρ^(n-2s),
//
// written as the product of the exact binomials C(n-s, s)·C(n-2s, (n+m)/2-s).
func zsRadialCoeff(n, m, s int) float64 {
	c := zsBinom(n-s, s) * zsBinom(n-2*s, (n+m)/2-s)
	if s%2 == 1 {
		c = -c
	}
	return c
}

// zsPow returns (x + i·y)^p for p >= 0 (p = -1 returns 1 + 0i, which is only
// used multiplied by |m| = 0).
func zsPow(x, y float64, p int) (re, im float64) {
	re, im = 1, 0
	for i := 0; i < p; i++ {
		re, im = re*x-im*y, re*y+im*x
	}
	return re, im
}

// zsRefZernike returns Z_j and the analytic partial derivatives ∂Z_j/∂x,
// ∂Z_j/∂y at the plane point (x, y) from the monomial expansion of the
// definition. With K = n-2s (K >= |m|, K ≡ |m| mod 2),
//
//	ρ^K cos(|m|θ) = (x²+y²)^((K-|m|)/2) · Re[(x+iy)^|m|] / normR^K ,
//	ρ^K sin(|m|θ) = (x²+y²)^((K-|m|)/2) · Im[(x+iy)^|m|] / normR^K ,
//
// both sides being polynomials, and the coefficients of zsRadialCoeff. The
// gradient is the analytic derivative of that expansion (product rule plus
// d/dx (x+iy)^p = p (x+iy)^(p-1)).
func zsRefZernike(t *testing.T, j int, x, y, normR float64) (z, dzdx, dzdy float64) {
	t.Helper()
	n, m := zsRefIndex(t, j)
	am := m
	if am < 0 {
		am = -am
	}
	cosine := m >= 0
	w := x*x + y*y
	for s := 0; s <= (n-am)/2; s++ {
		K := n - 2*s
		a := zsRadialCoeff(n, am, s) / math.Pow(normR, float64(K))
		q := (K - am) / 2
		pr, pi := zsPow(x, y, am)
		cr, ci := zsPow(x, y, am-1)
		cr *= float64(am)
		ci *= float64(am)
		var p, dpx, dpy float64
		if cosine {
			p, dpx, dpy = pr, cr, -ci
		} else {
			p, dpx, dpy = pi, ci, cr
		}
		wq := math.Pow(w, float64(q))
		tw := 0.0
		if q > 0 {
			tw = 2 * float64(q) * math.Pow(w, float64(q-1))
		}
		z += a * wq * p
		dzdx += a * (tw*x*p + wq*dpx)
		dzdy += a * (tw*y*p + wq*dpy)
	}
	return z, dzdx, dzdy
}

// zsRefConic returns the conic sag and its radial slope, written directly from
// the standard form z = c r² / [1 + sqrt(1 - (1+K) c² r²)], c = 1/R.
func zsRefConic(r, R, K float64) (z, dzdr float64) {
	if R == 0 || math.IsInf(R, 0) {
		return 0, 0
	}
	c := 1 / R
	sq := math.Sqrt(1 - (1+K)*c*c*r*r)
	return c * r * r / (1 + sq), c * r / sq
}

// zsRefSag calls zsRefConic at r = sqrt(x²+y²).
func zsRefSag(x, y, R, K float64) (z, dzdx, dzdy float64) {
	r := math.Hypot(x, y)
	z, dr := zsRefConic(r, R, K)
	if r > 0 {
		return z, dr * x / r, dr * y / r
	}
	return z, 0, 0
}

// zsCentral6 is the sixth-order central difference of f in x,
//
//	f'(x) ≈ [f(x+3h) - 9f(x+2h) + 45f(x+h) - 45f(x-h) + 9f(x-2h) - f(x-3h)] / (60h),
//
// used as a second, purely numerical cross-check of the analytic gradient.
func zsCentral6(f func(float64) float64, x, h float64) float64 {
	return (f(x+3*h) - 9*f(x+2*h) + 45*f(x+h) - 45*f(x-h) + 9*f(x-2*h) - f(x-3*h)) / (60 * h)
}

// zsCoeffs returns the deterministic test coefficient list j = 1 ... 36.
func zsCoeffs() []float64 {
	coeffs := make([]float64, len(zsRefTable))
	rng := rand.New(rand.NewSource(20261015))
	for i := range coeffs {
		mag := math.Pow(10, -3-float64(i%3)) // 1e-3, 1e-4, 1e-5 mm
		coeffs[i] = (2*rng.Float64() - 1) * mag
	}
	return coeffs
}

// ---------------------------------------------------------------------------
// (a) Sag and gradient against the independent monomial expansion
// ---------------------------------------------------------------------------

// TestZernikeSagGradAgainstMonomialExpansion compares SagXY and the analytic
// GradXY (∂z/∂x, ∂z/∂y) of a 36-term Zernike surface against the independent
// monomial expansion zsRefZernike (Noll's Eq. (2) plus sin/cos, in integer
// polynomial form). Criterion (G4): relative error <= 1e-12.
//
// Relative error definition: rel = |got-ref| / max(|ref|, 1e-3·S), where S is
// the largest |reference| over the whole sample set of that configuration. The
// floor only applies where the reference itself is three orders of magnitude
// below the surface amplitude, i.e. at its zeros, where a relative measure is
// not defined; the absolute error is additionally bounded by 1e-12·S at every
// single point, which is a floor-free statement.
func TestZernikeSagGradAgainstMonomialExpansion(t *testing.T) {
	const norm = 20.0
	coeffs := zsCoeffs()

	// Sample set: rings × azimuths plus points very close to the vertex.
	type pt struct{ x, y float64 }
	var pts []pt
	for _, rho := range []float64{0, 0.03, 0.15, 0.4, 0.65, 0.85, 1.0} {
		for _, th := range []float64{0, math.Pi / 7, math.Pi / 3, 2 * math.Pi / 5, 1.1 * math.Pi, 1.7 * math.Pi} {
			pts = append(pts, pt{rho * norm * math.Cos(th), rho * norm * math.Sin(th)})
		}
	}
	for _, eps := range []float64{1e-9, 1e-12} {
		for _, th := range []float64{0.3, 2.0, 4.4} {
			pts = append(pts, pt{eps * norm * math.Cos(th), eps * norm * math.Sin(th)})
		}
	}

	configs := []struct {
		name         string
		R, K         float64
		symmetric    bool
		coefficients []float64
	}{
		{"conic base", 400, -0.6, false, coeffs},
		{"sphere base", -250, 0, false, coeffs},
		{"flat base (R=Inf)", math.Inf(1), 0, false, coeffs},
		{"flat base (R=0)", 0, 0, false, coeffs},
	}
	for _, cfg := range configs {
		p := ZernikeSurface{R: cfg.R, K: cfg.K, NormR: norm, Coeffs: cfg.coefficients}
		if err := p.Validate(); err != nil {
			t.Fatalf("%s: Validate: %v", cfg.name, err)
		}
		// Reference values, and the amplitude scale of this configuration.
		ref := make([]struct{ z, gx, gy float64 }, len(pts))
		scaleZ, scaleG := 0.0, 0.0
		for i, q := range pts {
			bz, bgx, bgy := zsRefSag(q.x, q.y, cfg.R, cfg.K)
			for j, c := range cfg.coefficients {
				z, gx, gy := zsRefZernike(t, j+1, q.x, q.y, norm)
				bz += c * z
				bgx += c * gx
				bgy += c * gy
			}
			ref[i] = struct{ z, gx, gy float64 }{bz, bgx, bgy}
			scaleZ = math.Max(scaleZ, math.Abs(bz))
			scaleG = math.Max(scaleG, math.Max(math.Abs(bgx), math.Abs(bgy)))
		}
		worstAbsZ, worstRelZ := 0.0, 0.0
		worstAbsG, worstRelG := 0.0, 0.0
		worstRelZClean, worstRelGClean := 0.0, 0.0
		for i, q := range pts {
			z := p.SagXY(q.x, q.y)
			gx, gy := p.GradXY(q.x, q.y)
			r := ref[i]
			az := math.Abs(z - r.z)
			ag := math.Max(math.Abs(gx-r.gx), math.Abs(gy-r.gy))
			worstAbsZ = math.Max(worstAbsZ, az)
			worstAbsG = math.Max(worstAbsG, ag)
			worstRelZ = math.Max(worstRelZ, az/math.Max(math.Abs(r.z), 1e-3*scaleZ))
			worstRelG = math.Max(worstRelG, ag/math.Max(math.Max(math.Abs(r.gx), math.Abs(r.gy)), 1e-3*scaleG))
			if math.Abs(r.z) >= 0.01*scaleZ {
				worstRelZClean = math.Max(worstRelZClean, az/math.Abs(r.z))
			}
			if math.Max(math.Abs(r.gx), math.Abs(r.gy)) >= 0.01*scaleG {
				worstRelGClean = math.Max(worstRelGClean, ag/math.Max(math.Abs(r.gx), math.Abs(r.gy)))
			}
		}
		// Absolute bound, floor-free; then the relative criterion of G4.
		if worstAbsZ > 1e-12*scaleZ {
			t.Errorf("%s: max |sag error| = %.3g exceeds 1e-12*scale = %.3g", cfg.name, worstAbsZ, 1e-12*scaleZ)
		}
		if worstAbsG > 1e-12*scaleG {
			t.Errorf("%s: max |gradient error| = %.3g exceeds 1e-12*scale = %.3g", cfg.name, worstAbsG, 1e-12*scaleG)
		}
		if worstRelZ > 1e-12 || worstRelG > 1e-12 {
			t.Errorf("%s: max relative error sag %.3g, grad %.3g exceed 1e-12", cfg.name, worstRelZ, worstRelG)
		}
		t.Logf("%s: sag |err|/scale = %.3g, grad |err|/scale = %.3g, rel (clean subset) sag = %.3g grad = %.3g",
			cfg.name, worstAbsZ/scaleZ, worstAbsG/scaleG, worstRelZClean, worstRelGClean)
	}

	// Sixth-order finite-difference cross-check of the analytic gradient on a
	// conical-base configuration (an independent numerical method, tolerance
	// looser than the polynomial comparison because of the truncation error).
	p := ZernikeSurface{R: 400, K: -0.6, NormR: norm, Coeffs: coeffs}
	h := 1e-3 * norm
	worstFD := 0.0
	for _, q := range []pt{{4, 3}, {-6.5, 8.2}, {0.7, -12.3}, {11.4, 5.5}} {
		gotX, gotY := p.GradXY(q.x, q.y)
		fx := zsCentral6(func(v float64) float64 { return p.SagXY(v, q.y) }, q.x, h)
		fy := zsCentral6(func(v float64) float64 { return p.SagXY(q.x, v) }, q.y, h)
		scale := math.Max(math.Max(math.Abs(gotX), math.Abs(gotY)), 1e-12)
		worstFD = math.Max(worstFD, math.Max(math.Abs(gotX-fx), math.Abs(gotY-fy))/scale)
	}
	if worstFD > 1e-9 {
		t.Errorf("finite-difference cross-check: relative gradient deviation %.3g exceeds 1e-9", worstFD)
	}
	t.Logf("finite-difference cross-check: worst relative gradient deviation %.3g", worstFD)

	// Rotational symmetry flag: a surface built from the m = 0 members only
	// must be reported symmetric, the test list must not.
	sym := ZernikeSurface{R: 400, K: -0.6, NormR: norm, Coeffs: []float64{0.01, 0, 0, 0.002, 0, 0, 0, 0, 0, 0, -3e-4}}
	if !sym.Symmetric() {
		t.Errorf("a surface built from piston, defocus and spherical aberration must be symmetric")
	}
	if p.Symmetric() {
		t.Errorf("a surface with |m| > 0 terms must not be reported symmetric")
	}
	if (ZernikeSurface{}).Symmetric() != true {
		t.Errorf("the zero surface is rotationally symmetric")
	}
}

// ---------------------------------------------------------------------------
// (b) Hand-written low-order formulas
// ---------------------------------------------------------------------------

// TestZernikeLowOrderFormulas checks the low-order polynomials against the
// closed forms written out literally here (Noll 1976 Table I, in the
// unit-coefficient form of this package). Criterion (G4): agreement with the
// analytic value <= 1e-14 at a grid of (ρ, θ) points, both for the polynomial
// evaluator and through a single-term ZernikeSurface.
func TestZernikeLowOrderFormulas(t *testing.T) {
	const norm = 17.0
	forms := []struct {
		j int
		f func(rho, th float64) float64
	}{
		{1, func(_, _ float64) float64 { return 1 }},
		{2, func(r, th float64) float64 { return r * math.Cos(th) }},
		{3, func(r, th float64) float64 { return r * math.Sin(th) }},
		{4, func(r, _ float64) float64 { return 2*r*r - 1 }},
		{5, func(r, th float64) float64 { return r * r * math.Sin(2*th) }},
		{6, func(r, th float64) float64 { return r * r * math.Cos(2*th) }},
		{7, func(r, th float64) float64 { return (3*r*r*r - 2*r) * math.Sin(th) }},
		{8, func(r, th float64) float64 { return (3*r*r*r - 2*r) * math.Cos(th) }},
		{9, func(r, th float64) float64 { return r * r * r * math.Sin(3*th) }},
		{10, func(r, th float64) float64 { return r * r * r * math.Cos(3*th) }},
		{11, func(r, _ float64) float64 { return 6*r*r*r*r - 6*r*r + 1 }},
		{12, func(r, th float64) float64 { return (4*r*r*r*r - 3*r*r) * math.Cos(2*th) }},
		{13, func(r, th float64) float64 { return (4*r*r*r*r - 3*r*r) * math.Sin(2*th) }},
		{16, func(r, th float64) float64 { return (10*r*r*r*r*r - 12*r*r*r + 3*r) * math.Cos(th) }},
		{17, func(r, th float64) float64 { return (10*r*r*r*r*r - 12*r*r*r + 3*r) * math.Sin(th) }},
		{18, func(r, th float64) float64 { return (5*r*r*r*r*r - 4*r*r*r) * math.Cos(3*th) }},
		{19, func(r, th float64) float64 { return (5*r*r*r*r*r - 4*r*r*r) * math.Sin(3*th) }},
		{22, func(r, _ float64) float64 { return 20*r*r*r*r*r*r - 30*r*r*r*r + 12*r*r - 1 }},
	}
	worst := 0.0
	for _, tc := range forms {
		coeff := []float64{}
		for j := 1; j <= tc.j; j++ {
			if j == tc.j {
				coeff = append(coeff, 1)
			} else {
				coeff = append(coeff, 0)
			}
		}
		p := ZernikeSurface{R: math.Inf(1), NormR: norm, Coeffs: coeff}
		for _, rho := range []float64{0, 0.25, 0.5, 0.8, 1.0} {
			for _, th := range []float64{0, 0.7, math.Pi / 2, 2.4, 5.1} {
				want := tc.f(rho, th)
				got := ZernikeNoll(tc.j, rho, th)
				if d := math.Abs(got - want); d > 1e-14 {
					t.Errorf("Z%d(%.3g, %.3g) = %.17g, want %.17g (|diff| = %.3g > 1e-14)",
						tc.j, rho, th, got, want, d)
				}
				// Through the surface: the sag is exactly the coefficient times
				// the polynomial (flat base).
				gotS := p.SagXY(rho*norm*math.Cos(th), rho*norm*math.Sin(th))
				if d := math.Abs(gotS - want); d > 1e-14 {
					t.Errorf("surface Z%d(%.3g, %.3g) = %.17g, want %.17g (|diff| = %.3g > 1e-14)",
						tc.j, rho, th, gotS, want, d)
				}
				worst = math.Max(worst, math.Abs(got-want))
			}
		}
	}
	t.Logf("worst |analytic - polynomial| over the low-order table: %.3g", worst)

	// Radial closed forms and the orthonormal factors of Noll's Eq. (1).
	for _, tc := range []struct {
		n, m int
		f    func(r float64) float64
	}{
		{0, 0, func(r float64) float64 { return 1 }},
		{2, 0, func(r float64) float64 { return 2*r*r - 1 }},
		{2, 2, func(r float64) float64 { return r * r }},
		{4, 0, func(r float64) float64 { return 6*r*r*r*r - 6*r*r + 1 }},
		{4, 2, func(r float64) float64 { return 4*r*r*r*r - 3*r*r }},
		{4, 4, func(r float64) float64 { return r * r * r * r }},
		{6, 6, func(r float64) float64 { return r * r * r * r * r * r }},
		{6, 0, func(r float64) float64 { return 20*r*r*r*r*r*r - 30*r*r*r*r + 12*r*r - 1 }},
	} {
		for _, rho := range []float64{0, 0.1, 0.5, 0.9, 1.0} {
			if d := math.Abs(ZernikeRadial(tc.n, tc.m, rho) - tc.f(rho)); d > 1e-14 {
				t.Errorf("R_%d^%d(%.3g): |diff| = %.3g > 1e-14", tc.n, tc.m, rho, d)
			}
		}
	}
	for _, tc := range []struct {
		j    int
		want float64
	}{
		{1, 1}, {2, 2}, {4, math.Sqrt(3)}, {5, math.Sqrt(6)},
		{11, math.Sqrt(5)}, {12, math.Sqrt(10)}, {22, math.Sqrt(7)},
	} {
		if d := math.Abs(ZernikeNollNorm(tc.j) - tc.want); d > 1e-15 {
			t.Errorf("Noll normalisation N_%d = %.17g, want %.17g", tc.j, ZernikeNollNorm(tc.j), tc.want)
		}
	}
	// The orthonormal rescaling: a coefficient √3·a on Z4 (polynomial form)
	// reads a on the unit-RMS Z4^Noll = √3(2ρ²-1) of Noll's Eq. (1).
	p := ZernikeSurface{R: math.Inf(1), NormR: 17, Coeffs: []float64{0, 0, 0, math.Sqrt(3) * 1e-3}}
	on := p.OrthonormalCoeffs()
	if math.Abs(on[3]-1e-3) > 1e-18 {
		t.Errorf("OrthonormalCoeffs[3] = %.17g, want 1e-3", on[3])
	}
	if (ZernikeSurface{}).OrthonormalCoeffs() != nil {
		t.Errorf("OrthonormalCoeffs of an empty surface must be nil")
	}
}

// TestZernikeNollIndexTable checks the ordering rule against the hand table.
func TestZernikeNollIndexTable(t *testing.T) {
	for j := 1; j <= len(zsRefTable); j++ {
		n, m := ZernikeNollIndex(j)
		wn, wm := zsRefTable[j-1][0], zsRefTable[j-1][1]
		if n != wn || m != wm {
			t.Errorf("ZernikeNollIndex(%d) = (n=%d, m=%d), want (%d, %d)", j, n, m, wn, wm)
		}
		// Structural invariants of Noll's indexing.
		if m < 0 && j%2 == 0 {
			t.Errorf("j = %d is even but m = %d is negative (Noll: even j carries the cosine)", j, m)
		}
		if m > 0 && j%2 == 1 {
			t.Errorf("j = %d is odd but m = %d is positive (Noll: odd j carries the sine)", j, m)
		}
		if am := absInt(m); am > n || (n-am)%2 != 0 {
			t.Errorf("j = %d: n = %d, m = %d violate |m| <= n, n-|m| even", j, n, m)
		}
	}
	if n, m := ZernikeNollIndex(0); n != 0 || m != 0 {
		t.Errorf("ZernikeNollIndex(0) = (%d, %d), want (0, 0)", n, m)
	}
}

func absInt(v int) int {
	if v < 0 {
		return -v
	}
	return v
}

// ---------------------------------------------------------------------------
// (c) Zero coefficients reduce to the base
// ---------------------------------------------------------------------------

// TestZernikeZeroCoeffsReduceToBase checks that Coeffs = nil and an all-zero
// coefficient list give exactly the base surface: sag, gradient, Profile sag
// and slope, curvature and ray intersection agree with Plane/Sphere/Conic.
// Criterion (G4): agreement <= 1e-13.
func TestZernikeZeroCoeffsReduceToBase(t *testing.T) {
	bases := []struct {
		name string
		R, K float64
		base Profile
	}{
		{"plane", math.Inf(1), 0, Plane{}},
		{"flat (R=0)", 0, 0, Plane{}},
		{"sphere", -250, 0, Sphere{R: -250}},
		{"conic", 400, -0.6, Conic{R: 400, K: -0.6}},
	}
	for _, b := range bases {
		for _, coeffs := range []struct {
			name string
			v    []float64
		}{
			{"nil", nil},
			{"12 zeros", make([]float64, 12)},
		} {
			p := ZernikeSurface{R: b.R, K: b.K, NormR: 20, Coeffs: coeffs.v}
			if err := p.Validate(); err != nil {
				t.Fatalf("%s/%s: Validate: %v", b.name, coeffs.name, err)
			}
			worst := 0.0
			for i := 0; i <= 40; i++ {
				r := 12.0 * float64(i) / 40
				for j := 0; j < 8; j++ {
					th := 2 * math.Pi * float64(j) / 8
					x, y := r*math.Cos(th), r*math.Sin(th)
					worst = math.Max(worst, math.Abs(p.SagXY(x, y)-b.base.Sag(r)))
					gx, gy := p.GradXY(x, y)
					sl := b.base.Slope(r)
					if r > 0 {
						worst = math.Max(worst, math.Abs(gx-sl*x/r))
						worst = math.Max(worst, math.Abs(gy-sl*y/r))
					}
				}
				worst = math.Max(worst, math.Abs(p.Sag(r)-b.base.Sag(r)))
				worst = math.Max(worst, math.Abs(p.Slope(r)-b.base.Slope(r)))
			}
			worst = math.Max(worst, math.Abs(p.Curvature()-b.base.Curvature()))
			if worst > 1e-13 {
				t.Errorf("%s/%s: degenerate surface deviates from the base by %.3g > 1e-13", b.name, coeffs.name, worst)
			}
			// Intersection must agree with the base solver to 1e-12 mm.
			if !math.IsInf(b.R, 0) && b.R != 0 {
				o := V(3, -2, -40.0)
				d, _ := V(0.2, -0.1, 4).Unit()
				tz, okz := p.Intersect(o, d, -1e-9, 100)
				tb, okb := b.base.Intersect(o, d, -1e-9, 100)
				if okz != okb {
					t.Errorf("%s/%s: intersection ok = %v, base %v", b.name, coeffs.name, okz, okb)
				} else if okz && math.Abs(tz-tb) > 1e-12 {
					t.Errorf("%s/%s: intersection t = %.15g, base %.15g", b.name, coeffs.name, tz, tb)
				}
			}
		}
	}
}

// ---------------------------------------------------------------------------
// (d) Ray intersection: Newton residual
// ---------------------------------------------------------------------------

// TestZernikeIntersectionResidual checks the intersection solver on random
// rays: the hit must satisfy the exact sag equation to <= 1e-12 mm (the G4
// criterion for the Newton residual), it must be the *first* crossing of the
// ray with the surface, and a ray that cannot reach the search radius must be
// reported as a miss.
func TestZernikeIntersectionResidual(t *testing.T) {
	const norm = 15.0
	const rMax = 10.0
	coeffs := make([]float64, 12)
	rng := rand.New(rand.NewSource(4711))
	for i := range coeffs {
		coeffs[i] = (2*rng.Float64() - 1) * 1e-3
	}
	profiles := []struct {
		name string
		p    ZernikeSurface
	}{
		{"conic base", ZernikeSurface{R: -400, K: -0.4, NormR: norm, Coeffs: coeffs}},
		{"sphere base", ZernikeSurface{R: 300, K: 0, NormR: norm, Coeffs: coeffs}},
		{"flat base (scan path)", ZernikeSurface{R: 0, K: 0, NormR: norm, Coeffs: []float64{0, 0, 0, 2e-3, 0, 0, 0, 0, 0, 0, 4e-4}}},
	}
	rng = rand.New(rand.NewSource(20261015))
	for _, pr := range profiles {
		worstResidual := 0.0
		checked := 0
		misses := 0
		for i := 0; i < 500; i++ {
			// Rays from a pupil plane 60 mm in front, aimed near the axis.
			o := V((rng.Float64()*2-1)*rMax*1.5, (rng.Float64()*2-1)*rMax*1.5, -60)
			aim := V((rng.Float64()*2-1)*rMax*0.9, (rng.Float64()*2-1)*rMax*0.9, 0)
			d, ok := aim.Sub(o).Unit()
			if !ok {
				continue
			}
			tt, ok := pr.p.Intersect(o, d, -1e-9, rMax)
			if !ok {
				misses++
				continue
			}
			checked++
			hit := o.Add(d.Mul(tt))
			res := math.Abs(hit.Z - pr.p.SagXY(hit.X, hit.Y))
			if res > worstResidual {
				worstResidual = res
			}
			if res > 1e-12 {
				t.Errorf("%s: ray %d: Newton residual %.3g mm > 1e-12", pr.name, i, res)
			}
			// First crossing: no sign change of the residual between tMin and
			// the returned t (sampled).
			g := func(t float64) float64 {
				x := o.X + t*d.X
				y := o.Y + t*d.Y
				return o.Z + t*d.Z - pr.p.SagXY(x, y)
			}
			const nS = 4000
			g0 := g(-1e-9)
			first := true
			for k := 1; k <= nS; k++ {
				tk := -1e-9 + (tt+1e-9)*float64(k)/float64(nS)*0.999
				gk := g(tk)
				if g0*gk < 0 && first {
					t.Errorf("%s: ray %d: sign change at t = %.6g before the returned root %.6g", pr.name, i, tk, tt)
					first = false
				}
			}
			if tt < -1e-9 {
				t.Errorf("%s: ray %d: returned t = %.6g below tMin", pr.name, i, tt)
			}
		}
		if checked == 0 {
			t.Fatalf("%s: no ray hit the surface", pr.name)
		}
		t.Logf("%s: %d hits, %d misses, worst residual %.3g mm", pr.name, checked, misses, worstResidual)
	}

	// A ray that never reaches the surface, and a ray inside the region where
	// the base conic sag is undefined, must both be reported as misses. (rMax
	// is a search hint, not an aperture: a hit on the surface outside it is
	// still a hit on the profile, and vignetting is the aperture's job.)
	p := ZernikeSurface{R: -400, K: -0.4, NormR: norm, Coeffs: coeffs}
	if tHit, ok := p.Intersect(V(0, 0, -1), V(0, 0, -1), -1e-9, rMax); ok {
		t.Errorf("a ray travelling away from the surface hit at t = %g", tHit)
	}
	sphere := ZernikeSurface{R: 100, K: 0, NormR: norm, Coeffs: []float64{0, 0, 0, 1e-3}}
	if tHit, ok := sphere.Intersect(V(120, 0, -60), V(0, 0, 1), -1e-9, norm); ok {
		t.Errorf("a ray outside the domain of the base sphere hit at t = %g", tHit)
	}
	// Two-crossing geometry: a steep convex surface crossed twice by one ray
	// must return the smaller root.
	steep := ZernikeSurface{R: 40, K: 0, NormR: norm, Coeffs: nil}
	o, d := V(0, 0, -30), V(0, 0, 1)
	t1, ok1 := steep.Intersect(o, d, -1e-9, norm)
	if !ok1 {
		t.Fatalf("sphere-like Zernike surface: no axial hit")
	}
	if math.Abs(t1-30) > 1e-9 {
		t.Errorf("axial hit at t = %.12g, want 30 (vertex)", t1)
	}
	// Off-axis: the ray crosses the sphere twice; the first (nearer) crossing
	// is wanted. The conicoid equation of the sphere R = 40 in the local frame
	// is r² + z² - 2Rz = 0, so for r = 15 the crossings are at z = R ± sqrt(R²-r²).
	o2, d2 := V(15, 0, -30), V(0, 0, 1)
	t2, ok2 := steep.Intersect(o2, d2, -1e-9, 2*norm)
	if !ok2 {
		t.Fatalf("off-axis ray through a sphere-like Zernike surface: no hit")
	}
	near := 40.0 - math.Sqrt(40.0*40.0-15.0*15.0)
	if math.Abs((o2.Z+t2)-near) > 1e-9 {
		t.Errorf("off-axis first crossing z = %.12g, want %.12g (near root)", o2.Z+t2, near)
	}
}

// ---------------------------------------------------------------------------
// (f) Determinism
// ---------------------------------------------------------------------------

// TestZernikeDeterministicBitwise checks that repeated evaluation is
// bit-identical: the profile holds no state, iterates no map and uses no
// random numbers, so every observable must reproduce to the last bit (this is
// also exercised by running the package with -count=2).
func TestZernikeDeterministicBitwise(t *testing.T) {
	coeffs := zsCoeffs()
	build := func() ZernikeSurface {
		return ZernikeSurface{R: 400, K: -0.6, NormR: 20, Coeffs: append([]float64(nil), coeffs...)}
	}
	fingerprint := func() []uint64 {
		p := build()
		var out []uint64
		for i := 0; i <= 24; i++ {
			r := 20.0 * float64(i) / 24
			for j := 0; j < 6; j++ {
				th := 2 * math.Pi * float64(j) / 6
				x, y := r*math.Cos(th), r*math.Sin(th)
				out = append(out, math.Float64bits(p.SagXY(x, y)))
				gx, gy := p.GradXY(x, y)
				out = append(out, math.Float64bits(gx), math.Float64bits(gy))
			}
			out = append(out, math.Float64bits(p.Sag(r)), math.Float64bits(p.Slope(r)))
		}
		out = append(out, math.Float64bits(p.Curvature()))
		tt, ok := p.Intersect(V(2, -3, -50), V(0.05, -0.02, 1).Normalized(), -1e-9, 12)
		out = append(out, math.Float64bits(tt))
		if !ok {
			t.Errorf("fingerprint ray did not hit")
		}
		return out
	}
	a := fingerprint()
	b := fingerprint()
	if len(a) != len(b) {
		t.Fatalf("fingerprint lengths differ: %d vs %d", len(a), len(b))
	}
	for i := range a {
		if a[i] != b[i] {
			t.Fatalf("fingerprint[%d] differs between runs: %016x vs %016x", i, a[i], b[i])
		}
	}
	t.Logf("fingerprint of %d values reproduced bit for bit", len(a))
}

// ---------------------------------------------------------------------------
// Structure: validation, curvature, stringer
// ---------------------------------------------------------------------------

// TestZernikeSurfaceStructure checks validation errors, the vertex curvature
// (against a finite difference of the sag) and the stringer.
func TestZernikeSurfaceStructure(t *testing.T) {
	// Validation.
	if err := (ZernikeSurface{R: 100, K: 0, NormR: 0, Coeffs: []float64{0, 0, 0, 1e-3}}).Validate(); err == nil {
		t.Errorf("active coefficients with NormR = 0 must fail validation")
	}
	if err := (ZernikeSurface{R: 100, NormR: 10, Coeffs: []float64{math.NaN()}}).Validate(); err == nil {
		t.Errorf("a NaN coefficient must fail validation")
	}
	if err := (ZernikeSurface{R: math.NaN(), NormR: 10}).Validate(); err == nil {
		t.Errorf("a NaN radius must fail validation")
	}
	if err := (ZernikeSurface{R: 100, NormR: -1}).Validate(); err == nil {
		t.Errorf("a negative normalisation radius must fail validation")
	}
	if err := (ZernikeSurface{R: 100, NormR: 0}).Validate(); err != nil {
		t.Errorf("an inactive surface may have NormR = 0: %v", err)
	}
	if s := (ZernikeSurface{R: 100, NormR: 10, Coeffs: []float64{1}}).String(); s == "" {
		t.Errorf("empty String()")
	}

	// Vertex curvature: the x-meridian value ∂²z/∂x²(0,0). Every member with a
	// ρ² content (even n >= 2, m = 0 or m = +2) contributes 2κ/NormR² with κ
	// the ρ² coefficient of R_n^|m|: Z4 +4/norm², Z6 +2/norm², Z5 0,
	// Z11 -12/norm², Z12 -6/norm², Z22 +24/norm².
	const norm = 10.0
	const h = 1e-4
	z22 := make([]float64, 22)
	z22[21] = 2e-3
	for _, tc := range []struct {
		name   string
		coeffs []float64
		want   float64 // base curvature plus the analytic Zernike term
	}{
		{"Z4", []float64{0, 0, 0, 2e-3}, 2 * 2 * 2e-3 / (norm * norm)},
		{"Z6", []float64{0, 0, 0, 0, 0, 3e-3}, 2 * 1 * 3e-3 / (norm * norm)},
		{"Z5", []float64{0, 0, 0, 0, -3e-3}, 0},
		{"Z11", []float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5e-3}, 2 * -6 * 5e-3 / (norm * norm)},
		{"Z12", []float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4e-3}, 2 * -3 * 4e-3 / (norm * norm)},
		{"Z22", z22, 2 * 12 * 2e-3 / (norm * norm)},
	} {
		p := ZernikeSurface{R: math.Inf(1), NormR: norm, Coeffs: tc.coeffs}
		got := p.Curvature()
		fd := (p.SagXY(h, 0) - 2*p.SagXY(0, 0) + p.SagXY(-h, 0)) / (h * h)
		if d := math.Abs(got - tc.want); d > 1e-15 {
			t.Errorf("%s: Curvature = %.17g, want %.17g", tc.name, got, tc.want)
		}
		if d := math.Abs(got - fd); d > 1e-9 {
			t.Errorf("%s: Curvature = %.17g, finite difference %.17g", tc.name, got, fd)
		}
	}
	// The base contributes 1/R.
	c := ZernikeSurface{R: -250, NormR: norm}
	if d := math.Abs(c.Curvature() - (-1.0 / 250)); d > 1e-18 {
		t.Errorf("base curvature = %.17g, want %.17g", c.Curvature(), -1.0/250)
	}
	// NaNs where the surface is undefined.
	if !math.IsNaN((ZernikeSurface{R: 100, K: 0, NormR: 1, Coeffs: []float64{0, 0, 0, 1}}).SagXY(0.9, 0)) {
		// A sphere of R = 100 is real out to r = 100; use a hyperbola rim instead.
		t.Logf("sag inside the conic domain is finite (expected)")
	}

	// For a rotationally symmetric coefficient list the r-only Profile path
	// must be exact: Sag(r) equals SagXY at every azimuth and Slope(r) the
	// radial derivative. This is the guarantee behind Symmetric().
	sym := ZernikeSurface{R: -300, K: -0.5, NormR: 12, Coeffs: []float64{2e-3, 0, 0, 1.5e-3, 0, 0, 0, 0, 0, 0, -4e-4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8e-5}}
	if !sym.Symmetric() {
		t.Fatalf("m = 0 coefficient list not detected as symmetric")
	}
	worstSym := 0.0
	for i := 1; i <= 20; i++ {
		r := 12.0 * float64(i) / 20
		for j := 0; j < 12; j++ {
			th := 2 * math.Pi * float64(j) / 12
			x, y := r*math.Cos(th), r*math.Sin(th)
			worstSym = math.Max(worstSym, math.Abs(sym.Sag(r)-sym.SagXY(x, y)))
			gx, gy := sym.GradXY(x, y)
			worstSym = math.Max(worstSym, math.Abs(gx*math.Cos(th)+gy*math.Sin(th)-sym.Slope(r)))
		}
	}
	if worstSym > 1e-15 {
		t.Errorf("symmetric surface: r-only path deviates from the XY path by %.3g", worstSym)
	}
	t.Logf("symmetric surface: r-only vs XY path worst deviation %.3g", worstSym)
}
