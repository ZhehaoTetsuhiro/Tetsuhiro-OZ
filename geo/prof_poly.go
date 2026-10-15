package geo

import (
	"fmt"
	"math"
)

// PolyTerm is one monomial a·r^Power of a polynomial asphere. Power must be a
// positive integer; odd powers (r, r³, r⁵, …) are allowed as well as even
// ones. Coeff has units mm^(1−Power).
type PolyTerm struct {
	Power int
	Coeff float64
}

// PolyAsphere is a rotationally symmetric general-polynomial asphere on a conic
// base:
//
//	z(r) = c r² / (1 + sqrt(1 − (1+K) c² r²)) + Σ_i a_i r^i ,   c = 1/R ,
//
// where the sum runs over the Terms, each with an arbitrary positive integer
// power i. R and K have the meaning of Conic; R = 0 or ±Inf gives a flat base,
// i.e. a pure polynomial surface. Terms with the same power add.
//
// Notes on odd powers:
//   - a₁ r makes the vertex a cone (|r| is not differentiable on the axis). The
//     sag, slope and ray intersection are well defined away from the axis, but
//     the paraxial curvature (Curvature) ignores it, because a cone has no
//     vertex curvature. Use real-ray analysis for such surfaces.
//   - a₃, a₅, … do not change the vertex curvature; first-order optics is
//     insensitive to them.
//
// The even-only case reproduces EvenAsphere pointwise (see prof_poly_test.go).
type PolyAsphere struct {
	R     float64
	K     float64
	Terms []PolyTerm
}

// Validate reports structural errors: NaN base parameters, powers below 1 and
// non-finite coefficients.
func (p PolyAsphere) Validate() error {
	if math.IsNaN(p.R) || math.IsNaN(p.K) {
		return fmt.Errorf("geo: poly asphere: NaN base parameter (R=%v K=%v)", p.R, p.K)
	}
	for i, t := range p.Terms {
		if t.Power < 1 {
			return fmt.Errorf("geo: poly asphere term %d: power %d is not a positive integer", i, t.Power)
		}
		if math.IsNaN(t.Coeff) || math.IsInf(t.Coeff, 0) {
			return fmt.Errorf("geo: poly asphere term %d (r^%d): coefficient %v is not finite", i, t.Power, t.Coeff)
		}
	}
	return nil
}

// Sag implements Profile: base conic plus the polynomial sum, evaluated term by
// term with powInt so that sparse and odd powers need no special casing.
func (p PolyAsphere) Sag(r float64) float64 {
	z := conicBaseSag(r, p.R, p.K)
	if math.IsNaN(z) {
		return z
	}
	s := 0.0
	for _, t := range p.Terms {
		s += t.Coeff * powInt(r, t.Power)
	}
	return z + s
}

// Slope implements Profile: d/dr of the base plus Σ i·a_i·r^(i−1).
func (p PolyAsphere) Slope(r float64) float64 {
	s := conicBaseSlope(r, p.R, p.K)
	if math.IsNaN(s) {
		return s
	}
	for _, t := range p.Terms {
		s += float64(t.Power) * t.Coeff * powInt(r, t.Power-1)
	}
	return s
}

// Curvature implements Profile: the paraxial vertex curvature
//
//	c_eff = c + 2·a₂ ,
//
// i.e. twice the r² coefficient of the sag, which is what first-order optics
// sees. The base contributes c = 1/R (0 for a flat base); the r⁴ and higher
// terms, and the r¹ cone, do not change the vertex curvature.
func (p PolyAsphere) Curvature() float64 {
	c := baseCurvature(p.R)
	for _, t := range p.Terms {
		if t.Power == 2 {
			c += 2 * t.Coeff
		}
	}
	return c
}

// Intersect implements Profile: the conic base gives analytic Newton seeds (see
// solveIntersection); the result is then polished to machine precision on the
// exact sag equation.
func (p PolyAsphere) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	hasBase := !math.IsInf(p.R, 0) && p.R != 0
	t, ok := solveIntersection(p, p.R, p.K, hasBase, o, d, tMin, rMax)
	if !ok {
		return 0, false
	}
	return polishIntersection(p, o, d, t, tMin), true
}

// String implements Profile.
func (p PolyAsphere) String() string {
	s := fmt.Sprintf("poly asphere R=%g K=%g", p.R, p.K)
	for _, t := range p.Terms {
		s += fmt.Sprintf(" a%d=%g", t.Power, t.Coeff)
	}
	return s
}

// baseCurvature returns 1/R for a finite non-zero radius and 0 for a flat base
// (R = 0 or ±Inf). Unlike conicCurvature it never returns ±Inf for R = 0.
func baseCurvature(R float64) float64 {
	if R == 0 || math.IsInf(R, 0) {
		return 0
	}
	return 1 / R
}

// powInt returns x^n for an integer n >= 0 by binary exponentiation (a few
// ulp accurate; x^0 = 1 for every x, including 0).
func powInt(x float64, n int) float64 {
	result := 1.0
	base := x
	for n > 0 {
		if n&1 == 1 {
			result *= base
		}
		base *= base
		n >>= 1
	}
	return result
}

// polishIntersection refines the ray parameter t with guarded Newton steps on
// the sag residual g(t) = z(t) − Sag(r(t)). A step is accepted only when it
// keeps t >= tMin and strictly reduces |g|, so the polish can never move the
// hit to a different root or behind the ray origin.
func polishIntersection(p Profile, o, d Vec3, t, tMin float64) float64 {
	g, _ := sagResidual(p, o, d, t)
	for i := 0; i < 4 && !math.IsNaN(g) && g != 0; i++ {
		_, gd := sagResidual(p, o, d, t)
		if gd == 0 || math.IsNaN(gd) || math.IsInf(gd, 0) {
			break
		}
		tn := t - g/gd
		gn, _ := sagResidual(p, o, d, tn)
		if math.IsNaN(gn) || tn < tMin || math.Abs(gn) >= math.Abs(g) {
			break
		}
		t, g = tn, gn
	}
	return t
}
