package geo

import (
	"fmt"
	"math"
)

// This file adds the Zernike polynomial surface profile. The formulas and the
// Noll index table below follow the primary source,
//
//	R. J. Noll, "Zernike polynomials and atmospheric turbulence",
//	J. Opt. Soc. Am. 66(3), 207-211 (1976),
//
// Eq. (1)-(2) and Table I; the same mode order is OEIS A176988.

// ===========================================================================
// Zernike polynomials, Noll's sequential index
// ===========================================================================
//
// # Definition
//
// With the normalised radius ρ = r/R_norm (the unit disk is ρ <= 1) and the
// azimuth θ measured from +x towards +y (right-handed, θ = atan2(y, x)),
//
//	Z_j(ρ, θ) = R_n^|m|(ρ) · cos(m θ)    for m >= 0 (Noll's "even j"),
//	Z_j(ρ, θ) = R_n^|m|(ρ) · sin(|m| θ)  for m <  0 (Noll's "odd j"),
//
// where j -- Noll's sequential index, j = 1, 2, 3, ... -- is the function of
// the radial degree n and the signed azimuthal frequency m given by the
// ordering rule below, and R_n^m is the radial polynomial of Noll's Eq. (2),
//
//	R_n^m(ρ) = Σ_{s=0}^{(n-m)/2} (-1)^s
//	           · (n-s)! / [ s! · ((n+m)/2 - s)! · ((n-m)/2 - s)! ] · ρ^(n-2s) ,
//
// with n >= 0, m >= 0, m <= n and n - m even. In particular
//
//	R_n^n(ρ) = ρ^n ,     R_n^0(ρ) = P_n(2ρ²-1)  (Legendre),
//	R_2^0 = 2ρ²-1 ,      R_4^0 = 6ρ⁴-6ρ²+1 .
//
// The values in this package are the *polynomial* (unit-coefficient) forms:
// the coefficient c_j of a term multiplies exactly R_n^|m| times the sine or
// cosine above. Noll's Eq. (1) additionally divides them by the orthonormal
// factor
//
//	N_j = sqrt(2(n+1))   (m != 0),      N_j = sqrt(n+1)   (m = 0),
//
// which makes ∫∫ Z_j Z_j' d²r / π = δ_jj' on the unit disk (Noll's Eq. (3));
// ZernikeNollNorm returns it, so coefficients between the two conventions
// convert as c_orthonormal = c_polynomial / N_j. The polynomial form is the
// one in which the textbook low-order expressions read 2ρ²-1 (defocus) and
// 6ρ⁴-6ρ²+1 (primary spherical) without prefactors.
//
// # Ordering rule (Noll 1976, Table I caption)
//
// Modes are ordered by increasing radial degree n; within one n the modes
// with the lower |m| come first, and even j carries the cosine (m >= 0) while
// odd j carries the sine (m < 0) member of each |m| pair. Concretely, n is
// the smallest integer with j <= (n+1)(n+2)/2, k = j - n(n+1)/2 is the
// 1-based position inside the row (k = 1 ... n+1), and
//
//	|m| = 0                                        (n even, k = 1),
//	|m| = 2·floor(k/2)                             (n even, k > 1),
//	|m| = 2·floor((k+1)/2) - 1                     (n odd),
//	m   = +|m|  (cosine)  if j is even,  -|m| (sine) if j is odd.
//
// The resulting table (j, n, m, polynomial) begins
//
//	j  n  m    Z_j(ρ, θ)
//	1  0  0    1
//	2  1  1    ρ cos θ
//	3  1 -1    ρ sin θ
//	4  2  0    2ρ² - 1                                     (defocus)
//	5  2 -2    ρ² sin 2θ
//	6  2  2    ρ² cos 2θ                                   (astigmatism)
//	7  3 -1    (3ρ³ - 2ρ) sin θ
//	8  3  1    (3ρ³ - 2ρ) cos θ                            (coma)
//	9  3 -3    ρ³ sin 3θ
//	10 3  3    ρ³ cos 3θ                                   (trefoil)
//	11 4  0    6ρ⁴ - 6ρ² + 1                               (primary spherical)
//	12 4  2    (4ρ⁴ - 3ρ²) cos 2θ
//	13 4 -2    (4ρ⁴ - 3ρ²) sin 2θ
//	14 4  4    ρ⁴ cos 4θ
//	15 4 -4    ρ⁴ sin 4θ
//	16 5  1    (10ρ⁵ - 12ρ³ + 3ρ) cos θ
//	17 5 -1    (10ρ⁵ - 12ρ³ + 3ρ) sin θ
//	18 5  3    (5ρ⁵ - 4ρ³) cos 3θ
//	19 5 -3    (5ρ⁵ - 4ρ³) sin 3θ
//	20 5  5    ρ⁵ cos 5θ
//	21 5 -5    ρ⁵ sin 5θ
//	22 6  0    20ρ⁶ - 30ρ⁴ + 12ρ² - 1                      (5th order spherical)
//
// (These are the unit-coefficient polynomials; Noll's Table I prints the same
// modes with the factor of Eq. (1) included, e.g. Z2 = 2ρ cos θ there and
// Z4 = √3(2ρ²-1).)
//
// # Derivatives
//
// The code evaluates the radial polynomial as the Jacobi polynomial
//
//	R_n^m(ρ) = ρ^m · P_k^{(0,m)}(2ρ² - 1) ,    k = (n-m)/2 ,
//
// (equivalent to the sum above; both evaluate the same polynomial, the sum is
// used as the independent reference in prof_zernike_test.go) and forms the
// partial derivatives of Z_j = R(ρ)·Θ(θ) by the chain rule,
//
//	∂Z_j/∂x = (dR/dρ)·Θ·(x/r) / R_norm + R·(dΘ/dθ)·(-y/r²) ,
//	∂Z_j/∂y = (dR/dρ)·Θ·(y/r) / R_norm + R·(dΘ/dθ)·(+x/r²) ,
//
//	R' = m ρ^(m-1) P_k(2ρ²-1) + 4ρ·ρ^m·P_k'(2ρ²-1) ,
//	dΘ/dθ = -|m| sin(|m|θ)  (cosine members),  +|m| cos(|m|θ)  (sine members).
//
// At the vertex r = 0 the polar form is singular; the gradient there is the
// exact limit of the polynomial, which vanishes for m = 0 and |m| >= 2 and
// equals P_k^{(0,|m|)}(-1)/R_norm along +x (cosine members) or +y (sine
// members) for |m| = 1 -- this is the linear (tilt) term. No finite
// differences are used anywhere.

// ZernikeNollIndex returns the radial degree n and the signed azimuthal
// frequency m of Noll's sequential index j >= 1 (see the file comment for the
// ordering rule). j < 1 returns (0, 0).
func ZernikeNollIndex(j int) (n, m int) {
	if j < 1 {
		return 0, 0
	}
	n = 0
	for (n+1)*(n+2)/2 < j {
		n++
	}
	k := j - n*(n+1)/2 // 1 ... n+1
	am := 0
	if n%2 == 0 {
		am = 2 * (k / 2) // |m| = 0, 2, 2, 4, 4, ...
	} else {
		am = 2*((k+1)/2) - 1 // |m| = 1, 1, 3, 3, ...
	}
	if j%2 == 0 {
		m = am // even j: cosine (m >= 0)
	} else {
		m = -am // odd j: sine (m < 0)
	}
	return n, m
}

// ZernikeNollNorm returns Noll's orthonormal factor N_j of Eq. (1) that turns
// the polynomial form of this package into the unit-RMS form,
//
//	N_j = sqrt(2(n+1))  (m != 0),   N_j = sqrt(n+1)  (m = 0).
//
// Coefficients convert with c_orthonormal = c_polynomial / N_j.
func ZernikeNollNorm(j int) float64 {
	n, m := ZernikeNollIndex(j)
	if m != 0 {
		return math.Sqrt(2 * float64(n+1))
	}
	return math.Sqrt(float64(n + 1))
}

// ZernikeNoll returns the value Y_j(ρ, θ) of the j-th Zernike polynomial in
// Noll's ordering (the polynomial form, without Noll's orthonormal factor; see
// the file comment). j < 1 returns 0.
func ZernikeNoll(j int, rho, theta float64) float64 {
	if j < 1 {
		return 0
	}
	n, m := ZernikeNollIndex(j)
	am := m
	if am < 0 {
		am = -am
	}
	r := zernikeRadial(n, am, rho)
	if m >= 0 {
		return r * math.Cos(float64(am)*theta)
	}
	return r * math.Sin(float64(am)*theta)
}

// ZernikeRadial returns the radial polynomial R_n^m(ρ) for n >= m >= 0 with
// n - m even (Noll 1976 Eq. (2)); it returns 0 for m > n or odd n - m.
func ZernikeRadial(n, m int, rho float64) float64 {
	if m < 0 || n < 0 || m > n || (n-m)%2 != 0 {
		return 0
	}
	return zernikeRadial(n, m, rho)
}

// zernikeRadial evaluates R_n^m(ρ), n >= m >= 0, n - m even.
func zernikeRadial(n, m int, rho float64) float64 {
	k := (n - m) / 2
	p, _ := zernikeJacobi(k, m, 2*rho*rho-1)
	return powInt(rho, m) * p
}

// zernikeJacobi returns P_k^{(0,m)}(x) and dP/dx by the three-term Jacobi
// recurrence of, e.g., Abramowitz & Stegun 22.7.1 with α = 0, β = m,
//
//	P_0 = 1,  P_1 = [(α-β) + (α+β+2) x]/2 = [(m+2) x - m]/2,
//	P_{i+1} = (A_i x + B_i) P_i - C_i P_{i-1},
//	A_i = (2i+m+1)(2i+m+2) / (2(i+1)(i+m+1)),
//	B_i = -m² (2i+m+1) / (2(i+1)(i+m+1)(2i+m)),
//	C_i = i(i+m)(2i+m+2) / ((i+1)(i+m+1)(2i+m))        (C_0 = 0).
//
// (The classical reduction R_n^m(ρ) = ρ^m P_{(n-m)/2}^{(0,m)}(2ρ²-1) is
// checked against the explicit sum of Noll's Eq. (2) in prof_zernike_test.go.)
func zernikeJacobi(k, m int, x float64) (p, dp float64) {
	fm := float64(m)
	p0, dp0 := 1.0, 0.0
	if k == 0 {
		return p0, dp0
	}
	p1 := 0.5 * ((fm+2)*x - fm)
	dp1 := 0.5 * (fm + 2)
	if k == 1 {
		return p1, dp1
	}
	pm2, dpm2 := p0, dp0
	pm1, dpm1 := p1, dp1
	for i := 1; i < k; i++ {
		fi := float64(i)
		a := (2*fi + fm + 1) * (2*fi + fm + 2) / (2 * (fi + 1) * (fi + fm + 1))
		b := 0.0
		if m != 0 {
			b = -fm * fm * (2*fi + fm + 1) / (2 * (fi + 1) * (fi + fm + 1) * (2*fi + fm))
		}
		c := 0.0
		if i > 0 {
			c = fi * (fi + fm) * (2*fi + fm + 2) / ((fi + 1) * (fi + fm + 1) * (2*fi + fm))
		}
		lin := a*x + b
		pm := lin*pm1 - c*pm2
		dpm := a*pm1 + lin*dpm1 - c*dpm2
		pm2, dpm2 = pm1, dpm1
		pm1, dpm1 = pm, dpm
	}
	return pm1, dpm1
}

// zernikeTerm returns Z_j and its partial derivatives ∂Z_j/∂x, ∂Z_j/∂y at the
// plane point (x, y), for the unit-disk normalisation radius normR > 0.
func zernikeTerm(j int, x, y, normR float64) (z, dzdx, dzdy float64) {
	n, m := ZernikeNollIndex(j)
	am := m
	if am < 0 {
		am = -am
	}
	r := math.Hypot(x, y)
	rho := r / normR
	k := (n - am) / 2
	p, dp := zernikeJacobi(k, am, 2*rho*rho-1)
	rm := powInt(rho, am)
	theta := math.Atan2(y, x)
	var th, dth float64 // Θ(θ) and dΘ/dθ
	if m >= 0 {
		th, dth = math.Cos(float64(am)*theta), -float64(am)*math.Sin(float64(am)*theta)
	} else {
		th, dth = math.Sin(float64(am)*theta), float64(am)*math.Cos(float64(am)*theta)
	}
	z = rm * p * th
	if r > 0 {
		// dR/dρ = am·ρ^(am-1)·P + 4ρ·ρ^am·P' ; ∂ρ/∂x = x/(r·normR), ∂θ/∂x = -y/r².
		dR := float64(am)*powInt(rho, am-1)*p + 4*rho*rm*dp
		dzdx = dR*th*x/(r*normR) - rm*p*dth*y/(r*r)
		dzdy = dR*th*y/(r*normR) + rm*p*dth*x/(r*r)
		return z, dzdx, dzdy
	}
	// Vertex r = 0, where the polar basis is singular: take the exact limit of
	// the polynomial. The gradient vanishes for m = 0 and |m| >= 2; for
	// |m| = 1 the leading term is P_k^{(0,1)}(-1)·(x or y)/normR.
	switch am {
	case 0:
		return z, 0, 0
	case 1:
		p, _ := zernikeJacobi(k, 1, -1)
		if m > 0 {
			return z, p / normR, 0
		}
		return z, 0, p / normR
	}
	return z, 0, 0
}

// ===========================================================================
// ZernikeSurface
// ===========================================================================

// ZernikeSurface is a surface whose sag is a base conic plus a sum of Zernike
// polynomials over a normalisation circle of radius NormR,
//
//	z(x, y) = conic(R, K; r) + Σ_j Coeffs[j-1] · Z_j(ρ, θ) ,
//	r = sqrt(x² + y²) ,  ρ = r / NormR ,  θ = atan2(y, x) ,
//
// where Z_j is the j-th Zernike polynomial in Noll's sequential index (see the
// file comment for the definition, the ordering rule and the table; Noll 1976
// Eq. (1)-(2), Table I) and Coeffs[j-1] is the coefficient of Noll index j,
// so the coefficient list is a deterministic, positional enumeration starting
// at piston: Coeffs[3] multiplies Z4 (defocus), Coeffs[10] multiplies Z11
// (primary spherical). Trailing zero coefficients are simply omitted; nil or
// all-zero Coeffs gives the bare base, which is Conic{R, K} (Plane for
// R = 0 or ±Inf).
//
// The base conic has the meaning of Conic: R is the vertex radius of curvature
// (mm) with the package sign convention and K the conic constant; R = 0 or
// ±Inf is a flat base. The Zernike sum is evaluated on the whole plane (as a
// polynomial it extends beyond ρ = 1), but it describes the surface only
// inside the normalisation circle, so the surface should carry a clear
// aperture of at most NormR: NormR is normally the semi-aperture of the part,
// exactly like the ρ_max of the Q-type aspheres.
//
// # Rotationally symmetric and general terms
//
// The Profile interface of this package describes z = Sag(r), a surface of
// revolution: it has no θ. ZernikeSurface implements it with the θ = 0 cut,
// Sag(r) = SagXY(r, 0) and Slope(r) = ∂z/∂x(r, 0), which is *exact* for
// surfaces built from the rotationally symmetric members only (m = 0: piston,
// defocus, spherical aberration, ...). Such a surface traces, paraxially
// analyses and focuses correctly through the whole existing machinery; use
// Symmetric to test for it. A surface with |m| > 0 terms (tilt, astigmatism,
// coma, trefoil, ...) is *not* a surface of revolution and must be evaluated
// through SagXY/GradXY/NormalXY, which are exact for every term; feeding it to
// the r-only machinery (Surface.Normal, a trace) is a user error and is
// documented as such in the report that introduced this file.
//
// # Verification status
//
// All formulas above are the ones in the cited literature; the sag, the
// analytic partial derivatives, the low-order polynomials, the degeneracy to
// the base and the ray intersection are verified numerically in
// prof_zernike_test.go (sag and derivatives against an independent monomial
// expansion, <= 1e-12 relative) and prof_zernike_wave_test.go (Z4 on a mirror
// against the wave bridge).
type ZernikeSurface struct {
	R      float64   // base vertex radius of curvature (mm); 0 or ±Inf = flat base
	K      float64   // base conic constant
	NormR  float64   // normalisation radius R_norm (mm), ρ = r/NormR; must be > 0 when any coefficient is non-zero
	Coeffs []float64 // Coeffs[j-1] multiplies Z_j, Noll index j = 1, 2, 3, ...
}

// Validate reports structural errors: NaN base or normalisation parameters,
// non-finite coefficients and a non-positive normalisation radius while the
// Zernike part is active.
func (p ZernikeSurface) Validate() error {
	if math.IsNaN(p.R) || math.IsNaN(p.K) || math.IsNaN(p.NormR) {
		return fmt.Errorf("geo: zernike surface: NaN parameter (R=%v K=%v norm=%v)", p.R, p.K, p.NormR)
	}
	if p.NormR < 0 {
		return fmt.Errorf("geo: zernike surface: normalisation radius must not be negative (got %g)", p.NormR)
	}
	for i, c := range p.Coeffs {
		if math.IsNaN(c) || math.IsInf(c, 0) {
			return fmt.Errorf("geo: zernike surface: coefficient %d (Noll j=%d) is not finite (%v)", i, i+1, c)
		}
	}
	if p.active() && !(p.NormR > 0) {
		return fmt.Errorf("geo: zernike surface: normalisation radius must be > 0 (got %g)", p.NormR)
	}
	return nil
}

// active reports whether any Zernike coefficient is non-zero.
func (p ZernikeSurface) active() bool {
	for _, c := range p.Coeffs {
		if c != 0 {
			return true
		}
	}
	return false
}

// Symmetric reports whether every non-zero coefficient belongs to the
// rotationally symmetric members (m = 0, i.e. Noll j = 1, 4, 11, 22, ...). For
// such a surface the whole r-only Profile path (Sag, Slope, a trace through
// Surface) is exact; otherwise those routines describe the θ = 0 cut only and
// the surface must be evaluated through SagXY/GradXY/NormalXY.
func (p ZernikeSurface) Symmetric() bool {
	for i, c := range p.Coeffs {
		if c == 0 {
			continue
		}
		if _, m := ZernikeNollIndex(i + 1); m != 0 {
			return false
		}
	}
	return true
}

// SagXY returns the sag z(x, y) at the local point (x, y) — the base conic
// plus the Zernike sum, see ZernikeSurface. It is the exact sag of every term
// and the entry point for non-rotationally symmetric surfaces.
func (p ZernikeSurface) SagXY(x, y float64) float64 {
	z := conicBaseSag(math.Hypot(x, y), p.R, p.K)
	if math.IsNaN(z) || !p.active() {
		return z
	}
	if !(p.NormR > 0) {
		return math.NaN()
	}
	for i, c := range p.Coeffs {
		if c == 0 {
			continue
		}
		zk, _, _ := zernikeTerm(i+1, x, y, p.NormR)
		z += c * zk
	}
	return z
}

// GradXY returns the analytic partial derivatives ∂z/∂x and ∂z/∂y of SagXY:
// the base conic slope and the Zernike derivatives of the file comment. No
// finite differences are involved.
func (p ZernikeSurface) GradXY(x, y float64) (dzdx, dzdy float64) {
	r := math.Hypot(x, y)
	sl := conicBaseSlope(r, p.R, p.K)
	if math.IsNaN(sl) {
		return math.NaN(), math.NaN()
	}
	if r > 0 {
		dzdx, dzdy = sl*x/r, sl*y/r
	}
	if !p.active() {
		return dzdx, dzdy
	}
	if !(p.NormR > 0) {
		return math.NaN(), math.NaN()
	}
	for i, c := range p.Coeffs {
		if c == 0 {
			continue
		}
		_, ax, ay := zernikeTerm(i+1, x, y, p.NormR)
		dzdx += c * ax
		dzdy += c * ay
	}
	return dzdx, dzdy
}

// NormalXY returns the unit surface normal at the local point (x, y),
// oriented towards +z exactly like Surface.Normal,
//
//	n ∝ (-∂z/∂x, -∂z/∂y, 1) ,
//
// computed from the analytic gradient GradXY. It is exact for every Zernike
// term, including the ones with |m| > 0. A zero Vec3 is returned for a
// degenerate gradient (never the case for a finite surface).
func (p ZernikeSurface) NormalXY(x, y float64) Vec3 {
	dzdx, dzdy := p.GradXY(x, y)
	n, ok := Vec3{-dzdx, -dzdy, 1}.Unit()
	if !ok {
		return Vec3{}
	}
	return n
}

// Sag implements Profile with the θ = 0 cut, Sag(r) = SagXY(r, 0). It is the
// exact sag of a surface built from the rotationally symmetric members alone
// (see Symmetric) and the θ = 0 meridian otherwise.
func (p ZernikeSurface) Sag(r float64) float64 { return p.SagXY(r, 0) }

// Slope implements Profile with the θ = 0 cut, Slope(r) = ∂z/∂x at (r, 0),
// i.e. d/dr Sag(r). See Sag.
func (p ZernikeSurface) Slope(r float64) float64 {
	dzdx, _ := p.GradXY(r, 0)
	return dzdx
}

// Curvature implements Profile: the meridional vertex curvature
// ∂²z/∂x²(0, 0), which is the vertex curvature 1/R for the base and the only
// shape information first-order optics uses. Every Zernike member whose radial
// polynomial contains a ρ² term contributes, i.e. every even radial degree
// n >= 2 with |m| <= 2 (m = 0 and the cosine member m = 2; the sine member
// m = -2 has no x² content). With κ_n^m the coefficient of ρ² in R_n^|m|,
//
//	Z4   (2ρ²-1)          contributes 2κ/NormR² = 4/NormR² ,
//	Z6   (ρ² cos 2θ)      contributes 2κ/NormR² = 2/NormR² ,
//	Z11  (6ρ⁴-6ρ²+1)      contributes 2κ/NormR² = -12/NormR² ,
//	Z22  (20ρ⁶-30ρ⁴+12ρ²-1) contributes 2κ/NormR² = +24/NormR² ,
//
// because ρ² = (x²+y²)/NormR². For a surface with |m| > 2 or odd m terms first
// order optics is not defined; this method returns the x-meridian value.
func (p ZernikeSurface) Curvature() float64 {
	c := baseCurvature(p.R)
	if !p.active() || !(p.NormR > 0) {
		return c
	}
	n2 := p.NormR * p.NormR
	for i, v := range p.Coeffs {
		if v == 0 {
			continue
		}
		n, m := ZernikeNollIndex(i + 1)
		if n < 2 || n%2 != 0 || (m != 0 && m != 2) {
			continue // no ρ² content (m = -2 is the sine member, also none)
		}
		c += 2 * zernikeRho2(n, m) * v / n2
	}
	return c
}

// zernikeRho2 returns the coefficient κ of ρ² in R_n^m(ρ) for even n >= 2 and
// m in {0, 2} (the members with a quadratic part), and 0 otherwise:
//
//	κ = (-1)^s · C(n-s, s) · C(2, (n+m)/2 - s) ,   s = (n-2)/2 ,
//
// i.e. the s-th coefficient of Noll's Eq. (2) sum with n - 2s = 2:
// κ = 2 for (n,m) = (2,0), κ = -6 for (4,0), κ = +12 for (6,0),
// κ = 1 for (2,2), κ = -3 for (4,2), κ = +6 for (6,2), as in
// R_2^0 = 2ρ²-1, R_4^0 = 6ρ⁴-6ρ²+1, R_2^2 = ρ², R_4^2 = 4ρ⁴-3ρ².
func zernikeRho2(n, m int) float64 {
	if n < 2 || n%2 != 0 || (m != 0 && m != 2) {
		return 0
	}
	s := (n - 2) / 2
	κ := zernikeBinom(n-s, s) * zernikeBinom(2, (n+m)/2-s)
	if s%2 == 1 {
		κ = -κ
	}
	return κ
}

// zernikeBinom returns C(n, k) as a float64, exactly for the small arguments
// used by zernikeRho2.
func zernikeBinom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	r := 1.0
	for i := 1; i <= k; i++ {
		r = r * float64(n-k+i) / float64(i)
	}
	return r
}

// OrthonormalCoeffs returns the coefficient list rescaled to Noll's unit-RMS
// normalisation of Eq. (1): out[j-1] = Coeffs[j-1] / ZernikeNollNorm(j), so
// that the same surface reads Σ_j out[j-1] · Z_j^Noll with the orthonormal
// polynomials Z_j^Noll (∫∫ Z_j Z_k d²r/π = δ_jk). The result is nil when no
// coefficient is set.
func (p ZernikeSurface) OrthonormalCoeffs() []float64 {
	if len(p.Coeffs) == 0 {
		return nil
	}
	out := make([]float64, len(p.Coeffs))
	for i, c := range p.Coeffs {
		out[i] = c / ZernikeNollNorm(i+1)
	}
	return out
}

// residual evaluates the intersection residual g(t) = z(t) - SagXY(x(t), y(t))
// and its derivative dg/dt = d.z - ∇z·d_xy for the ray o + t·d (d unit).
func (p ZernikeSurface) residual(o, d Vec3, t float64) (g, gd float64) {
	x := o.X + t*d.X
	y := o.Y + t*d.Y
	z := o.Z + t*d.Z
	s := p.SagXY(x, y)
	if math.IsNaN(s) {
		return math.NaN(), math.NaN()
	}
	dzdx, dzdy := p.GradXY(x, y)
	if math.IsNaN(dzdx) {
		return math.NaN(), math.NaN()
	}
	return z - s, d.Z - dzdx*d.X - dzdy*d.Y
}

// newton polishes t towards a root of the sag residual with guarded Newton
// steps on the analytic derivative, reporting ok only at convergence level for
// the local surface scale.
func (p ZernikeSurface) newton(o, d Vec3, t, tMin float64) (float64, bool) {
	for i := 0; i < 60; i++ {
		g, gd := p.residual(o, d, t)
		if math.IsNaN(g) {
			return 0, false
		}
		z := o.Z + t*d.Z
		r := math.Hypot(o.X+t*d.X, o.Y+t*d.Y)
		if math.Abs(g) <= 1e-14*(1+math.Abs(z)+r) {
			if t < tMin {
				return 0, false
			}
			return t, true
		}
		if gd == 0 || math.IsNaN(gd) || math.IsInf(gd, 0) {
			return 0, false
		}
		step := -g / gd
		if math.IsNaN(step) || math.Abs(step) > 1e12 {
			return 0, false
		}
		t += step
	}
	g, _ := p.residual(o, d, t)
	z := o.Z + t*d.Z
	r := math.Hypot(o.X+t*d.X, o.Y+t*d.Y)
	return t, t >= tMin && math.Abs(g) <= 1e-11*(1+math.Abs(z)+r)
}

// polish squeezes the last ulps out of an intersection parameter: Newton steps
// are accepted only while they keep t >= tMin and strictly reduce |g|, so the
// polish can never jump to another root or behind the ray origin.
func (p ZernikeSurface) polish(o, d Vec3, t, tMin float64) float64 {
	g, _ := p.residual(o, d, t)
	for i := 0; i < 8 && !math.IsNaN(g) && g != 0; i++ {
		_, gd := p.residual(o, d, t)
		if gd == 0 || math.IsNaN(gd) || math.IsInf(gd, 0) {
			break
		}
		tn := t - g/gd
		gn, _ := p.residual(o, d, tn)
		if math.IsNaN(gn) || tn < tMin || math.Abs(gn) >= math.Abs(g) {
			break
		}
		t, g = tn, gn
	}
	return t
}

// zernikeWindow derives a conservative [tA, tB] window that contains the
// intersections of the ray with the surface inside the radial search radius.
// Unlike the rotationally symmetric searchWindow it samples SagXY over a polar
// grid of the disk, which bounds the sag of a non-symmetric surface.
func (p ZernikeSurface) zernikeWindow(o, d Vec3, tMin, rMax float64) (float64, float64, bool) {
	rSearch := rMax
	if !(rSearch > 0) || math.IsInf(rSearch, 1) {
		rSearch = o.PerpNorm() + 1
	}
	lo := math.Inf(1)
	hi := math.Inf(-1)
	add := func(t float64) {
		if math.IsNaN(t) || math.IsInf(t, 0) {
			return
		}
		if t < lo {
			lo = t
		}
		if t > hi {
			hi = t
		}
	}
	// z extent of the sag over the search disk.
	zmin, zmax := math.Inf(1), math.Inf(-1)
	const nR, nT = 24, 16
	for i := 0; i <= nR; i++ {
		r := rSearch * float64(i) / float64(nR)
		for j := 0; j < nT; j++ {
			th := 2 * math.Pi * float64(j) / float64(nT)
			z := p.SagXY(r*math.Cos(th), r*math.Sin(th))
			if math.IsNaN(z) || math.IsInf(z, 0) {
				continue
			}
			if z < zmin {
				zmin = z
			}
			if z > zmax {
				zmax = z
			}
		}
	}
	if zmin <= zmax {
		pad := 0.5*(zmax-zmin) + 1e-6
		zmin -= pad
		zmax += pad
		if math.Abs(d.Z) > 1e-12 {
			add((zmin - o.Z) / d.Z)
			add((zmax - o.Z) / d.Z)
		}
	}
	// radial crossings |o_perp + t d_perp| = rSearch.
	a := d.X*d.X + d.Y*d.Y
	b := 2 * (o.X*d.X + o.Y*d.Y)
	c := o.X*o.X + o.Y*o.Y - rSearch*rSearch
	if a > 1e-300 {
		disc := b*b - 4*a*c
		if disc >= 0 {
			sq := math.Sqrt(disc)
			add((-b - sq) / (2 * a))
			add((-b + sq) / (2 * a))
		}
	}
	if math.IsInf(lo, 1) || math.IsInf(hi, -1) || hi <= lo {
		return 0, 0, false
	}
	if lo < tMin {
		lo = tMin
	}
	if hi <= lo {
		return 0, 0, false
	}
	return lo, hi, true
}

// scanIntersect brackets the first sign change of the residual inside the
// search window, bisects it and polishes it with Newton. It is the fallback
// for intersections the base-conic seeds do not reach.
func (p ZernikeSurface) scanIntersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	tA, tB, ok := p.zernikeWindow(o, d, tMin, rMax)
	if !ok {
		return 0, false
	}
	const n = 400
	prevT := tA
	prevG, _ := p.residual(o, d, prevT)
	for i := 1; i <= n; i++ {
		t := tA + (tB-tA)*float64(i)/float64(n)
		g, _ := p.residual(o, d, t)
		if !math.IsNaN(g) && !math.IsNaN(prevG) && prevG*g <= 0 {
			lo, hi := prevT, t
			glo := prevG
			for j := 0; j < 120 && hi-lo > 1e-15*(1+math.Abs(lo)); j++ {
				mid := 0.5 * (lo + hi)
				gm, _ := p.residual(o, d, mid)
				if math.IsNaN(gm) {
					break
				}
				if (glo > 0) == (gm > 0) {
					lo, glo = mid, gm
				} else {
					hi = mid
				}
			}
			root := 0.5 * (lo + hi)
			if rt, ok := p.newton(o, d, root, tMin); ok && rt >= tMin {
				return p.polish(o, d, rt, tMin), true
			}
			if root >= tMin {
				return p.polish(o, d, root, tMin), true
			}
		}
		prevT, prevG = t, g
	}
	return 0, false
}

// Intersect implements Profile: the intersection with the exact sag equation
// z = SagXY(x, y). The analytic roots of the base conicoid (when there is one)
// seed a Newton iteration on the true residual; when those miss, a bracketed
// scan over a sag-bounded window finds the first crossing. The winner is the
// smallest valid t >= tMin, polished on the analytic derivative.
func (p ZernikeSurface) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	if o.IsNaN() || d.IsNaN() {
		return 0, false
	}
	if p.active() && !(p.NormR > 0) {
		return 0, false
	}
	hasBase := !math.IsInf(p.R, 0) && p.R != 0
	if hasBase {
		best := math.Inf(1)
		found := false
		for _, t0 := range conicRoots(p.R, p.K, o, d) {
			if t0 < tMin {
				continue
			}
			if t, ok := p.newton(o, d, t0, tMin); ok && t >= tMin && t < best {
				best, found = t, true
			}
		}
		if found {
			return p.polish(o, d, best, tMin), true
		}
	}
	return p.scanIntersect(o, d, tMin, rMax)
}

// String implements Profile.
func (p ZernikeSurface) String() string {
	return fmt.Sprintf("zernike R=%g K=%g norm=%g (%d coeffs)", p.R, p.K, p.NormR, len(p.Coeffs))
}
