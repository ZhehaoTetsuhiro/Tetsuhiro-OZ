package geo

import (
	"fmt"
	"math"
	"sort"
)

// ---------------------------------------------------------------------------
// Freeform XY polynomial surface
// ---------------------------------------------------------------------------
//
// The formulas and conventions of the freeform surface, written once here and
// repeated in abridged form on the profiles themselves.
//
// Sag. An XYPolynomial is a freeform surface whose sag in its own local frame
// (vertex at the origin, local optical axis = +z) is
//
//	z(x, y) = conicBase(r) + Σ_k a_k · (x/N)^i_k · (y/N)^j_k ,   r = √(x²+y²),
//
// where
//
//   - conicBase is the optional conic base of revolution of Conic (R, K),
//     evaluated exactly as Conic does: R is the vertex radius of curvature
//     (mm, same sign convention as Conic: R > 0 puts the centre of curvature
//     on the +z side), K the conic constant, and R = 0 or ±Inf selects a flat
//     base (the pure freeform case). The base is rotationally symmetric and
//     contributes only to the sag through r;
//   - the departure is a polynomial in the *local* transverse coordinates with
//     integer exponents i_k, j_k ≥ 0 and i_k + j_k ≥ 2 (no piston and no
//     tilt: the vertex is at the origin with a horizontal tangent plane, and
//     the surface slope at the vertex is exactly zero);
//   - N = NormR is the optional normalisation radius (mm). With NormR > 0 the
//     coefficients a_k carry the same meaning as the sag in mm at the point
//     (i_k = 2, j_k = 0: z = a·x²/N², so a is a sag in mm at x = N); with
//     NormR ≤ 0 the raw coordinates are used and a_k has units mm^(1−i−j).
//     Normalised coordinates keep the coefficients O(sag) even for high
//     orders, which is what vendors' "normalised XY polynomial" surfaces do.
//
// Terms are kept in a canonical order — ascending total degree i+j, then
// ascending i, then j — which makes the floating-point evaluation order a
// deterministic function of the term set (see Canonical). Terms with equal
// (i, j) are summed. Evaluation is bit-for-bit repeatable.
//
// Derivative and normal. The partial derivatives are analytic, term by term:
//
//	∂z/∂x = c x / √(1−(1+K)c²r²) + (1/N) Σ i_k a_k (x/N)^(i_k−1) (y/N)^j_k ,
//	∂z/∂y = c y / √(1−(1+K)c²r²) + (1/N) Σ j_k a_k (x/N)^i_k (y/N)^(j_k−1) ,
//
// with c = 1/R for the base (zero for a flat base), and the unit surface normal
// follows the convention of the package (Normal):
//
//	n ∝ (−∂z/∂x, −∂z/∂y, 1) ,  oriented towards +z.
//
// At the vertex r = 0 the gradient vanishes identically (every term has
// i+j ≥ 2), so the normal there is exactly (0, 0, 1).
//
// Intersection. There is no closed-form ray intersection for a general XY
// polynomial. IntersectXY solves the residual
//
//	g(t) = z(t) − z(x(t), y(t)) = 0 ,  (x(t), y(t), z(t)) = o + t·d ,
//
// by Newton iteration (dt = −g/g'), seeded with the analytic roots of the base
// conicoid when the base is active and with the vertex tangent plane otherwise,
// followed by a bracketed bisection scan when Newton fails. The scan window is
// derived from the sag range of the surface over the aperture disc and from the
// crossings of the ray with that disc, so it is valid off the axis as well as
// on it: nothing in the solver assumes rotational symmetry (in particular the
// intersection is *not* reduced to the r-meridian). A failure to converge is
// reported as ok = false / an error (IntersectRay); no NaN or ±Inf ever leaves
// the solver. The returned t is the smallest root with t ≥ tMin, up to the
// resolution of the fallback scan (see IntersectXY).
//
// Applicability. Freeform surfaces are non-coaxial by nature. The paraxial
// (first-order) engine and the third-order (Seidel) sums of this package are
// defined for rotationally symmetric coaxial systems only; a freeform surface
// whose departure is not a polynomial in r² is rejected by them exactly as a
// decentred or tilted surface is (see System.requireCoaxial), and
// System.Reversed refuses it as well. A rotationally symmetric XYPolynomial
// (see RotationallySymmetric) degenerates to an ordinary rotationally
// symmetric profile and stays fully usable with them; its Curvature is the
// exact vertex curvature. Genuinely asymmetric freeform systems must be
// analysed with real rays (Trace, SpotDiagram, RayFan, OPDGrid and the
// geometric MTF built on them), placed with the coordinate-break machinery
// (Surface.Pose / Surface.Break).

// XYTerm is one monomial a·(x/N)^i·(y/N)^j of a freeform XY polynomial
// surface. I and J are the exponents of x and y (both ≥ 0, I+J ≥ 2) and Coeff
// is the coefficient a in mm (with a normalisation radius NormR > 0) or
// mm^(1−I−J) (without one). See XYPolynomial for the exact definition.
type XYTerm struct {
	I     int
	J     int
	Coeff float64
}

// String formats the monomial compactly.
func (t XYTerm) String() string {
	return fmt.Sprintf("a_%d%d=%g", t.I, t.J, t.Coeff)
}

// XYPolynomial is a freeform surface: an optional conic base of revolution plus
// a two-dimensional polynomial departure in the local transverse coordinates
// (see the package section "Freeform XY polynomial surface" above for the
// formulas and conventions). It implements Profile through the +y meridian, and
// the additional interface Freeform for the true two-dimensional shape.
//
// Field order for a literal: optional base (R, K), optional normalisation
// radius NormR, then the terms. A zero-value XYPolynomial is the plane z ≡ 0
// (no base, no terms), which is the identity surface of the family: with R = 0
// (or ±Inf) and all coefficients zero it reproduces Plane exactly.
type XYPolynomial struct {
	R     float64 // vertex radius of the conic base (mm); 0 or ±Inf: flat base
	K     float64 // conic constant of the base
	NormR float64 // normalisation radius N (mm); > 0: normalised terms, 0: raw x, y
	Terms []XYTerm
}

// Validate reports structural errors: non-finite base parameters, a negative
// normalisation radius, non-finite coefficients and monomials below the
// quadratic order (I+J < 2 would add a piston or a tilt, which the vertex
// convention of the package excludes).
func (p XYPolynomial) Validate() error {
	if math.IsNaN(p.R) || math.IsNaN(p.K) || math.IsNaN(p.NormR) {
		return fmt.Errorf("geo: xy polynomial: NaN parameter (R=%v K=%v norm=%v)", p.R, p.K, p.NormR)
	}
	if p.NormR < 0 {
		return fmt.Errorf("geo: xy polynomial: normalisation radius %g must be positive (or zero for raw coordinates)", p.NormR)
	}
	for k, t := range p.Terms {
		if t.I < 0 || t.J < 0 {
			return fmt.Errorf("geo: xy polynomial term %d: exponents (%d, %d) must not be negative", k, t.I, t.J)
		}
		if t.I+t.J < 2 {
			return fmt.Errorf("geo: xy polynomial term %d: x^%d y^%d is below the quadratic order (i+j = %d < 2)",
				k, t.I, t.J, t.I+t.J)
		}
		if math.IsNaN(t.Coeff) || math.IsInf(t.Coeff, 0) {
			return fmt.Errorf("geo: xy polynomial term %d (x^%d y^%d): coefficient %v is not finite", k, t.I, t.J, t.Coeff)
		}
	}
	return nil
}

// Canonical returns a copy of the surface whose terms are in the canonical
// order: ascending total degree i+j, then ascending i, then ascending j, with
// terms of equal (i, j) summed. The order is a deterministic function of the
// term set, so two surfaces describing the same shape evaluate in the same
// order and give bit-identical sags wherever they agree mathematically.
func (p XYPolynomial) Canonical() XYPolynomial {
	out := p
	out.Terms = append([]XYTerm(nil), p.Terms...)
	sort.Slice(out.Terms, func(a, b int) bool {
		ta, tb := out.Terms[a], out.Terms[b]
		if ta.I+ta.J != tb.I+tb.J {
			return ta.I+ta.J < tb.I+tb.J
		}
		if ta.I != tb.I {
			return ta.I < tb.I
		}
		return ta.J < tb.J
	})
	// Merge duplicate monomials; the canonical order makes the loop below a
	// fixed schedule, so the merged coefficients are deterministic too.
	w := 0
	for _, t := range out.Terms {
		if w > 0 {
			prev := &out.Terms[w-1]
			if prev.I == t.I && prev.J == t.J {
				prev.Coeff += t.Coeff
				continue
			}
		}
		out.Terms[w] = t
		w++
	}
	out.Terms = out.Terms[:w]
	return out
}

// NewXYPolynomial builds a freeform surface from a base (R, K), a
// normalisation radius (0 for raw coordinates) and a term list. The terms are
// canonicalised (see Canonical) and validated; the error is nil for a valid
// surface.
func NewXYPolynomial(R, K, normR float64, terms []XYTerm) (XYPolynomial, error) {
	p := XYPolynomial{R: R, K: K, NormR: normR, Terms: terms}.Canonical()
	if err := p.Validate(); err != nil {
		return XYPolynomial{}, err
	}
	return p, nil
}

// xyCoords maps (x, y) to the coordinates in which the terms are evaluated:
// (x/N, y/N) when a normalisation radius is set, the raw coordinates
// otherwise.
func (p XYPolynomial) xyCoords(x, y float64) (u, v float64) {
	if p.NormR > 0 {
		return x / p.NormR, y / p.NormR
	}
	return x, y
}

// departure evaluates Σ a_k u^i v^j in the order of the stored terms.
func (p XYPolynomial) departure(u, v float64) float64 {
	var s float64
	for _, t := range p.Terms {
		if t.Coeff == 0 {
			continue
		}
		s += t.Coeff * powInt(u, t.I) * powInt(v, t.J)
	}
	return s
}

// departureGrad evaluates (∂/∂u, ∂/∂v) of the departure.
func (p XYPolynomial) departureGrad(u, v float64) (du, dv float64) {
	for _, t := range p.Terms {
		if t.Coeff == 0 {
			continue
		}
		if t.I > 0 {
			du += float64(t.I) * t.Coeff * powInt(u, t.I-1) * powInt(v, t.J)
		}
		if t.J > 0 {
			dv += float64(t.J) * t.Coeff * powInt(u, t.I) * powInt(v, t.J-1)
		}
	}
	return du, dv
}

// SagXY returns the sag z(x, y) of the freeform surface in its local frame.
func (p XYPolynomial) SagXY(x, y float64) float64 {
	z := conicBaseSag(math.Hypot(x, y), p.R, p.K)
	if math.IsNaN(z) {
		return z
	}
	u, v := p.xyCoords(x, y)
	return z + p.departure(u, v)
}

// GradXY returns the analytic partial derivatives ∂z/∂x and ∂z/∂y at (x, y).
// They are NaN wherever the base conic is outside its domain (the same points
// where SagXY is NaN).
func (p XYPolynomial) GradXY(x, y float64) (zx, zy float64) {
	r := math.Hypot(x, y)
	slope := conicBaseSlope(r, p.R, p.K)
	if math.IsNaN(slope) {
		return math.NaN(), math.NaN()
	}
	if r > 0 {
		zx = slope * x / r
		zy = slope * y / r
	}
	u, v := p.xyCoords(x, y)
	du, dv := p.departureGrad(u, v)
	if p.NormR > 0 {
		inv := 1 / p.NormR
		du, dv = du*inv, dv*inv
	}
	return zx + du, zy + dv
}

// NormalXY returns the unit surface normal at the local point (x, y), oriented
// towards +z: n ∝ (−∂z/∂x, −∂z/∂y, 1). At the vertex it is exactly (0, 0, 1).
// The zero vector is returned where the surface is undefined (NaN gradient) or
// the normal cannot be normalised.
func (p XYPolynomial) NormalXY(x, y float64) Vec3 {
	zx, zy := p.GradXY(x, y)
	if math.IsNaN(zx) || math.IsNaN(zy) {
		return Vec3{}
	}
	n, ok := Vec3{-zx, -zy, 1}.Unit()
	if !ok {
		return Vec3{}
	}
	return n
}

// Sag implements Profile as the sag along the +y meridian, z(0, r): the
// rotationally symmetric surrogate of the freeform surface (exact when the
// surface itself is rotationally symmetric). The two-dimensional shape is
// SagXY.
func (p XYPolynomial) Sag(r float64) float64 { return p.SagXY(0, r) }

// Slope implements Profile as dz/dy along the +y meridian at (0, r).
func (p XYPolynomial) Slope(r float64) float64 {
	_, zy := p.GradXY(0, r)
	return zy
}

// Curvature implements Profile: the mean vertex curvature of the surface,
//
//	c_mean = c + (a_20 + a_02)/N² ,   c = 1/R (0 for a flat base),
//
// i.e. half the sum of the two vertex curvatures 2a_20/N² (along x) and
// 2a_02/N² (along y). For a rotationally symmetric freeform (a_20 = a_02) this
// is exactly the vertex curvature 1/R + 2a_20/N² of the equivalent
// rotationally symmetric profile, which is the only shape information the
// first-order engine uses. Terms of order i+j ≠ 2 and the cross term x·y do
// not change any vertex curvature. For an asymmetric surface the single number
// is a surrogate: first-order optics has no unique curvature there, and such
// systems must be analysed with real rays (see the package section on freeform
// surfaces).
func (p XYPolynomial) Curvature() float64 {
	c := baseCurvature(p.R)
	scale := 1.0
	if p.NormR > 0 {
		scale = 1 / (p.NormR * p.NormR)
	}
	for _, t := range p.Terms {
		if (t.I == 2 && t.J == 0) || (t.I == 0 && t.J == 2) {
			c += t.Coeff * scale
		}
	}
	return c
}

// String implements Profile.
func (p XYPolynomial) String() string {
	base := fmt.Sprintf("R=%g K=%g", p.R, p.K)
	if math.IsInf(p.R, 0) || p.R == 0 {
		base = "flat"
	}
	norm := ""
	if p.NormR > 0 {
		norm = fmt.Sprintf(" N=%g", p.NormR)
	}
	return fmt.Sprintf("xy polynomial %s%s (%d terms)", base, norm, len(p.Terms))
}

// RotationallySymmetric reports whether the surface (base plus departure) is
// invariant under rotations about the local z axis, i.e. whether the freeform
// polynomial is a polynomial in r² = x²+y². The test is exact on the
// coefficients: every monomial must have even exponents (2k, 2m−2k) and, with
// b_m the coefficient of x^(2m) in the departure, every monomial of total
// degree 2m must satisfy a_(2k, 2m−2k) = C(m, k)·b_m. (The base conic is
// rotationally symmetric by construction.)
//
// The property matters because the paraxial and third-order engines of the
// package are defined for rotationally symmetric surfaces only: a symmetric
// XYPolynomial degenerates to an ordinary rotationally symmetric profile and
// can be used with them, an asymmetric one is rejected (see
// System.requireCoaxial) and must be analysed with real rays.
func (p XYPolynomial) RotationallySymmetric() bool {
	q := p.Canonical()
	maxDeg := 0
	for _, t := range q.Terms {
		if t.I%2 != 0 || t.J%2 != 0 {
			return false
		}
		if d := t.I + t.J; d > maxDeg {
			maxDeg = d
		}
	}
	// b[m] is the coefficient of x^(2m), i.e. of the (2m, 0) monomial.
	b := make([]float64, maxDeg/2+1)
	for _, t := range q.Terms {
		if t.J == 0 {
			b[t.I/2] = t.Coeff
		}
	}
	for _, t := range q.Terms {
		m := (t.I + t.J) / 2
		if t.Coeff != b[m]*binom(m, t.I/2) {
			return false
		}
	}
	return true
}

// binom returns the binomial coefficient C(n, k) as a float64 (exact for the
// small n of freeform term sets).
func binom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	c := 1.0
	for i := 1; i <= k; i++ {
		c = c * float64(n-k+i) / float64(i)
	}
	return c
}

// ---------------------------------------------------------------------------
// FreeformProfile interface and the intersection solver
// ---------------------------------------------------------------------------

// FreeformProfile is implemented by profiles whose sag is a function of the two
// local transverse coordinates (x, y) and not of r alone. Surface and System
// dispatch to it wherever the two-dimensional shape matters (Sag, Normal,
// Intersect, ZRange); the rotationally symmetric Profile methods
// (Sag(r), Slope(r), Curvature) keep their meaning as the meridian surrogate
// described on each implementation, so a freeform surface can also be dropped
// into code that only understands the radial interface.
//
// Implementations must be pure functions of (x, y) and safe for concurrent
// use, and IntersectXY must follow the Profile.Intersect contract: the smallest
// parameter t ≥ tMin with o + t·d on the surface, d unit, ok = false when
// there is none — never a NaN t.
type FreeformProfile interface {
	Profile
	// SagXY returns the local sag z(x, y).
	SagXY(x, y float64) float64
	// NormalXY returns the unit surface normal at (x, y), oriented towards +z.
	NormalXY(x, y float64) Vec3
	// IntersectXY returns the smallest t ≥ tMin with o + t·d on the surface
	// (d unit), ok = false when the ray misses it or the solver fails.
	IntersectXY(o, d Vec3, tMin, rMax float64) (float64, bool)
	// ZRangeXY returns the smallest and largest sag over the aperture disc of
	// radius rMax centered on the vertex (an unbounded interval for a
	// non-positive or +Inf radius). It is a conservative bracketing bound for
	// the intersection search.
	ZRangeXY(rMax float64) (float64, float64)
	// RotationallySymmetric reports whether the surface is invariant under
	// rotations about its local z axis. A symmetric freeform surface is an
	// ordinary rotationally symmetric profile (the radial Profile methods are
	// exact for it and the first-order engine is valid); an asymmetric one is
	// rejected by the coaxial-only engines and by System.Reversed.
	RotationallySymmetric() bool
}

// Intersection tolerances of the freeform solver. xyResidualTol is the
// relative residual accepted by the Newton iteration (relative to
// 1 + |z| + r, the local scale of the geometry), xyScanTol the relaxed level
// accepted for a root produced by the bisection fallback. Both are well below
// the 1e-12 mm acceptance level of the freeform test suite; the residual that
// a returned intersection leaves behind is measured in
// prof_freeform_test.go.
const (
	xyResidualTol = 1e-14
	xyScanTol     = 1e-13
	// xyScanCoarse and xyScanFine are the sample counts of the pre-scan (a
	// cheap check for a root before the Newton solution) and of the full
	// fallback scan.
	xyScanCoarse = 64
	xyScanFine   = 512
)

// residualXY evaluates the intersection residual g(t) = z(t) − z(x(t), y(t))
// and its derivative dg/dt for the ray o + t·d (d unit).
func (p XYPolynomial) residualXY(o, d Vec3, t float64) (g, gd float64) {
	x := o.X + t*d.X
	y := o.Y + t*d.Y
	z := o.Z + t*d.Z
	sag := p.SagXY(x, y)
	if math.IsNaN(sag) {
		return math.NaN(), math.NaN()
	}
	zx, zy := p.GradXY(x, y)
	if math.IsNaN(zx) || math.IsNaN(zy) {
		return math.NaN(), math.NaN()
	}
	return z - sag, d.Z - (zx*d.X + zy*d.Y)
}

// xyTol returns the residual acceptance level at the ray parameter t.
func xyTol(o, d Vec3, t, rel float64) float64 {
	z := o.Z + t*d.Z
	r := math.Hypot(o.X+t*d.X, o.Y+t*d.Y)
	return rel * (1 + math.Abs(z) + r)
}

// newtonXY refines t towards a root of the freeform residual with guarded
// Newton steps. A step is rejected (and the iteration abandoned) when it is
// non-finite, leaves the ray (t < tMin), overshoots by more than a factor
// 1e9 · (1 + |t|) — which is how a polynomial at a large distance from the
// axis blows up — or drives the point off the finite range of float64. ok is
// true only when the residual reaches xyResidualTol at the local scale, so a
// returned t is always a genuine hit.
func (p XYPolynomial) newtonXY(o, d Vec3, t, tMin float64) (float64, bool) {
	for i := 0; i < 100; i++ {
		g, gd := p.residualXY(o, d, t)
		if math.IsNaN(g) {
			return 0, false
		}
		if math.Abs(g) <= xyTol(o, d, t, xyResidualTol) {
			return t, true
		}
		if gd == 0 || math.IsNaN(gd) || math.IsInf(gd, 0) {
			return 0, false
		}
		step := -g / gd
		if math.IsNaN(step) || math.IsInf(step, 0) || math.Abs(step) > 1e9*(1+math.Abs(t)) {
			return 0, false
		}
		t += step
		if t < tMin || math.Abs(t) > 1e12 {
			return 0, false
		}
	}
	g, _ := p.residualXY(o, d, t)
	return t, !math.IsNaN(g) && math.Abs(g) <= xyTol(o, d, t, xyScanTol)
}

// zRangeXYDisc returns the smallest and largest finite sag of the surface over
// the disc of radius r (a deterministic polar grid; see ZRangeXY).
func (p XYPolynomial) zRangeXYDisc(r float64) (float64, float64) {
	const nr, na = 16, 24
	zmin, zmax := math.Inf(1), math.Inf(-1)
	upd := func(x, y float64) {
		z := p.SagXY(x, y)
		if math.IsNaN(z) || math.IsInf(z, 0) {
			return
		}
		if z < zmin {
			zmin = z
		}
		if z > zmax {
			zmax = z
		}
	}
	upd(0, 0)
	for k := 1; k <= nr; k++ {
		rr := r * float64(k) / nr
		for j := 0; j < na; j++ {
			a := 2 * math.Pi * float64(j) / na
			upd(rr*math.Cos(a), rr*math.Sin(a))
		}
	}
	if zmin > zmax {
		return math.Inf(-1), math.Inf(1)
	}
	return zmin, zmax
}

// ZRangeXY implements FreeformProfile: the sag range over the aperture disc of
// radius rMax (a polar grid with 385 samples, deterministic). A non-positive
// or +Inf radius means "unbounded" and yields an infinite interval.
func (p XYPolynomial) ZRangeXY(rMax float64) (float64, float64) {
	if !(rMax > 0) || math.IsInf(rMax, 1) {
		return math.Inf(-1), math.Inf(1)
	}
	return p.zRangeXYDisc(rMax)
}

// searchWindowXY derives a conservative [tA, tB] window around the first
// intersection of the ray with the surface, using the sag bounding box over
// the search disc (radius rMax, or the ray's own transverse extent when rMax
// does not bound the surface) and the crossings of the ray with the disc.
// It is the freeform analogue of searchWindow and makes no symmetry assumption.
func (p XYPolynomial) searchWindowXY(o, d Vec3, tMin, rMax float64) (float64, float64, bool) {
	rSearch := rMax
	if !(rSearch > 0) || math.IsInf(rSearch, 1) {
		rSearch = o.PerpNorm() + 1
	}
	lo, hi := math.Inf(1), math.Inf(-1)
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
	zmin, zmax := p.ZRangeXY(rSearch)
	if !math.IsInf(zmin, -1) && !math.IsInf(zmax, 1) {
		pad := 0.25*(zmax-zmin) + 1e-3*(1+math.Abs(zmin)+math.Abs(zmax))
		zmin -= pad
		zmax += pad
		if math.Abs(d.Z) > 1e-12 {
			add((zmin - o.Z) / d.Z)
			add((zmax - o.Z) / d.Z)
		}
	}
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

// bisectXY bisects a sign-change bracket [lo, hi] of the residual to machine
// precision, polishes the root with guarded Newton steps and returns it when
// the residual reaches xyScanTol.
func (p XYPolynomial) bisectXY(o, d Vec3, lo, hi, tMin float64) (float64, bool) {
	glo, _ := p.residualXY(o, d, lo)
	if math.IsNaN(glo) {
		return 0, false
	}
	for i := 0; i < 200 && hi-lo > 1e-16*(1+math.Abs(lo)+math.Abs(hi)); i++ {
		mid := 0.5 * (lo + hi)
		gm, _ := p.residualXY(o, d, mid)
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
	if !finite(root) || root < tMin {
		return 0, false
	}
	if r, ok := p.newtonXY(o, d, root, tMin); ok {
		return r, true
	}
	// Newton failed (a nearly tangential hit leaves a tiny derivative): accept
	// the bisection root when its residual is at the scan acceptance level.
	g, _ := p.residualXY(o, d, root)
	if !math.IsNaN(g) && math.Abs(g) <= xyTol(o, d, root, xyScanTol) {
		return root, true
	}
	return 0, false
}

// scanXY walks the search window in n uniform steps and bisects the first sign
// change found. tStop caps the scan (< +Inf to look only for a root before
// tStop). ok is false when no sign change is bracketed.
func (p XYPolynomial) scanXY(o, d Vec3, tMin, tStop, rMax float64, n int) (float64, bool) {
	tA, tB, ok := p.searchWindowXY(o, d, tMin, rMax)
	if !ok {
		return 0, false
	}
	if tStop < tB {
		tB = tStop
	}
	if !(tB > tA) {
		return 0, false
	}
	prevT := tA
	prevG, _ := p.residualXY(o, d, prevT)
	for i := 1; i <= n; i++ {
		t := tA + (tB-tA)*float64(i)/float64(n)
		g, _ := p.residualXY(o, d, t)
		if !math.IsNaN(g) && !math.IsNaN(prevG) && (prevG == 0 || prevG*g <= 0) {
			if r, ok := p.bisectXY(o, d, prevT, t, tMin); ok {
				return r, true
			}
		}
		prevT, prevG = t, g
	}
	return 0, false
}

// IntersectXY implements FreeformProfile and the Profile.Intersect contract:
// the smallest t ≥ tMin with o + t·d on the surface (d unit). The candidates
// are the analytic roots of the base conicoid (Newton-polished on the exact
// freeform residual) and the vertex-plane crossing; when Newton finds a root, a
// coarse scan of the window before it checks for an earlier crossing, and when
// no seed converges the full bracketed scan runs. ok = false when there is no
// intersection: the solver never returns NaN or ±Inf.
func (p XYPolynomial) IntersectXY(o, d Vec3, tMin, rMax float64) (float64, bool) {
	if o.IsNaN() || d.IsNaN() || !d.UnitValid() || math.IsNaN(tMin) {
		return 0, false
	}
	best := math.Inf(1)
	found := false
	consider := func(t0 float64) {
		if math.IsNaN(t0) || math.IsInf(t0, 0) || t0 < tMin {
			return
		}
		if t, ok := p.newtonXY(o, d, t0, tMin); ok && t < best {
			best, found = t, true
		}
	}
	if !math.IsInf(p.R, 0) && p.R != 0 {
		// The base conicoid gives up to two analytic seeds; they bracket the
		// near and far branches of the surface.
		for _, t0 := range conicRoots(p.R, p.K, o, d) {
			consider(t0)
		}
	}
	if math.Abs(d.Z) > 1e-300 {
		consider(-o.Z / d.Z) // crossing of the vertex tangent plane
	}
	if found {
		// Newton converges to *a* root; make sure no earlier crossing of the
		// surface lies in front of it (a freeform can be crossed several
		// times). The check is coarse and only replaces the root when the
		// earlier crossing is unambiguous.
		if t, ok := p.scanXY(o, d, tMin, best, rMax, xyScanCoarse); ok &&
			t < best-1e-9*(1+math.Abs(best)) {
			return t, true
		}
		return best, true
	}
	return p.scanXY(o, d, tMin, math.Inf(1), rMax, xyScanFine)
}

// IntersectRay is IntersectXY with the failure reported as an error, for
// callers that need a reason instead of an ok flag. A successful call returns
// a finite t ≥ tMin and a nil error; no error message ever carries a NaN.
func (p XYPolynomial) IntersectRay(o, d Vec3, tMin, rMax float64) (float64, error) {
	if o.IsNaN() || d.IsNaN() {
		return 0, fmt.Errorf("geo: xy polynomial: ray origin or direction is not finite (o=%s d=%s)", o, d)
	}
	if !d.UnitValid() {
		return 0, fmt.Errorf("geo: xy polynomial: ray direction %s is not a unit vector", d)
	}
	t, ok := p.IntersectXY(o, d, tMin, rMax)
	if !ok {
		return 0, fmt.Errorf("geo: xy polynomial: ray from %s along %s does not intersect the surface for t >= %g", o, d, tMin)
	}
	return t, nil
}
