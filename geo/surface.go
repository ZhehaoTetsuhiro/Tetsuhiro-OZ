package geo

import (
	"fmt"
	"math"
)

// Profile is the rotationally symmetric sag description of a surface in its own
// local frame (vertex at the origin, optical axis = +z):
//
//	z = Sag(r),  r = sqrt(x² + y²) >= 0 .
//
// Slope returns dz/dr at r, which is used for the analytic surface normal
// (see Surface.Normal). Curvature returns the vertex curvature 1/R (0 for a
// plane), which is the only shape information first-order (paraxial) optics
// uses. Implementations are pure functions of r and must be safe for
// concurrent use.
//
// Intersect returns the smallest parameter t >= tMin with o + t·d on the
// profile, where d is a unit vector and rMax is a radial search hint (mm, <= 0
// or +Inf means "unknown/unbounded"); ok is false when there is no such
// intersection, when the profile is undefined there (NaN sag) or when the
// numerical solver fails to converge. tMin = 0 accepts an intersection at the
// ray's own origin, which is what zero-thickness (thin lens) systems need; a
// slightly negative tMin tolerates round-off for rays that start on a surface
// (see tMinDefault in system.go).
type Profile interface {
	Sag(r float64) float64
	Slope(r float64) float64
	Curvature() float64
	Intersect(o, d Vec3, tMin, rMax float64) (float64, bool)
	String() string
}

// ---------------------------------------------------------------------------
// Plane
// ---------------------------------------------------------------------------

// Plane is a flat surface, z ≡ 0 in the local frame (radius of curvature
// ±Infinity).
type Plane struct{}

// Sag returns 0.
func (Plane) Sag(r float64) float64 { return 0 }

// Slope returns 0.
func (Plane) Slope(r float64) float64 { return 0 }

// Curvature returns 0.
func (Plane) Curvature() float64 { return 0 }

// Intersect solves the linear equation o.z + t d.z = 0. A ray travelling
// parallel to the plane (d.z = 0) does not intersect it.
func (Plane) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	if d.Z == 0 || math.IsNaN(o.Z) {
		return 0, false
	}
	t := -o.Z / d.Z
	if t < tMin || math.IsInf(t, 0) || math.IsNaN(t) {
		return 0, false
	}
	return t, true
}

// String implements Profile.
func (Plane) String() string { return "plane" }

// ---------------------------------------------------------------------------
// Conic section: sphere (K = 0), paraboloid (K = -1), ellipsoid, hyperboloid
// ---------------------------------------------------------------------------

// Conic is a conicoid of revolution with vertex radius R (mm) and conic
// constant K. The sag is the standard, numerically stable form
//
//	z(r) = c r² / (1 + sqrt(1 - (1+K) c² r²)),  c = 1/R ,
//
// which equals r²/(2R) to leading order and reproduces
//
//	K = 0    sphere of radius R,
//	K = -1   paraboloid (paraxial radius R, exactly z = r²/(2R)),
//	-1 < K < 0   prolate ellipsoid,      K > 0  oblate ellipsoid,
//	K < -1   hyperboloid.
//
// The sag is real only where 1 - (1+K) c² r² >= 0, i.e. r <= |R|/sqrt(1+K) for
// K > -1 (the full aperture of a sphere is |R|); outside that domain Sag
// returns NaN and the surface cannot be hit there. For K <= -1 the domain is
// unbounded.
//
// R = ±Inf yields a plane.
type Conic struct {
	R float64 // vertex radius of curvature (mm); sign convention: R > 0 puts
	// the centre of curvature on the +z side of the vertex.
	K float64 // conic constant (dimensionless)
}

// conicSag evaluates the conic sag z(r) for radius R and conic constant K.
func conicSag(r, R, K float64) float64 {
	if math.IsInf(R, 0) {
		return 0
	}
	if R == 0 || math.IsNaN(R) || math.IsNaN(K) {
		return math.NaN()
	}
	if r == 0 {
		return 0
	}
	c := 1 / R
	a := 1 - (1+K)*c*c*r*r
	if a < 0 {
		return math.NaN()
	}
	// c r² / (1 + sqrt(a)): no catastrophic cancellation, since 1 + sqrt(a) ∈ [1, 2].
	return c * r * r / (1 + math.Sqrt(a))
}

// conicBaseSag evaluates the conic sag used as the base of an asphere. R = 0
// and R = ±Inf both mean "flat base" (pure polynomial asphere).
func conicBaseSag(r, R, K float64) float64 {
	if R == 0 || math.IsInf(R, 0) {
		return 0
	}
	return conicSag(r, R, K)
}

// conicBaseSlope is the slope of conicBaseSag.
func conicBaseSlope(r, R, K float64) float64 {
	if R == 0 || math.IsInf(R, 0) {
		return 0
	}
	return conicSlope(r, R, K)
}

// conicCurvature returns 1/R, and 0 for a plane (|R| = Inf).
func conicCurvature(R float64) float64 {
	if math.IsInf(R, 0) {
		return 0
	}
	return 1 / R
}

// conicSlope evaluates dz/dr = c r / sqrt(1 - (1+K) c² r²).
func conicSlope(r, R, K float64) float64 {
	if math.IsInf(R, 0) {
		return 0
	}
	if R == 0 || math.IsNaN(R) || math.IsNaN(K) {
		return math.NaN()
	}
	c := 1 / R
	a := 1 - (1+K)*c*c*r*r
	if a <= 0 {
		// At (or beyond) the rim of a sphere-like surface the slope is
		// vertical; report an infinity signs-compatible with c.
		return math.Copysign(math.Inf(1), c*r)
	}
	return c * r / math.Sqrt(a)
}

// Sag implements Profile.
func (p Conic) Sag(r float64) float64 { return conicSag(r, p.R, p.K) }

// Slope implements Profile.
func (p Conic) Slope(r float64) float64 { return conicSlope(r, p.R, p.K) }

// Curvature returns c = 1/R (0 for a plane).
func (p Conic) Curvature() float64 { return conicCurvature(p.R) }

// String implements Profile.
func (p Conic) String() string { return fmt.Sprintf("conic R=%g K=%g", p.R, p.K) }

// Intersect implements Profile: analytic quadratic solution followed by a
// Newton polish on the exact sag equation (which also enforces the correct
// branch of the conicoid). tMin = 0 is accepted, see Profile.
func (p Conic) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	return solveIntersection(p, p.R, p.K, true, o, d, tMin, rMax)
}

// Sphere is a spherical surface of radius R (mm); equivalently Conic{R, 0}. See
// Conic for the sign convention and the domain of the sag.
type Sphere struct {
	R float64
}

// Sag implements Profile.
func (p Sphere) Sag(r float64) float64 { return conicSag(r, p.R, 0) }

// Slope implements Profile.
func (p Sphere) Slope(r float64) float64 { return conicSlope(r, p.R, 0) }

// Curvature implements Profile.
func (p Sphere) Curvature() float64 { return conicCurvature(p.R) }

// Intersect implements Profile.
func (p Sphere) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	return solveIntersection(p, p.R, 0, true, o, d, tMin, rMax)
}

// String implements Profile.
func (p Sphere) String() string { return fmt.Sprintf("sphere R=%g", p.R) }

// ---------------------------------------------------------------------------
// Even-order polynomial asphere
// ---------------------------------------------------------------------------

// EvenAsphere is a rotationally symmetric even-order polynomial asphere on a
// conic base:
//
//	z(r) = c r² / (1 + sqrt(1 - (1+K) c² r²)) + Σ_i Coeffs[i] · r^(2i+4) ,
//	c = 1/R ,
//
// i.e. Coeffs[0] = A4, Coeffs[1] = A6, Coeffs[2] = A8 and so on (the standard
// A4/A6/A8... notation). R and K are interpreted exactly as in Conic; K = 0
// gives a spherical base, K = -1 a parabolic base. R = ±Inf gives a pure
// polynomial surface on a flat base.
//
// The aspheric series is not truncated by the type: give as many coefficients
// as needed, unused trailing coefficients are simply zero. The intersection
// uses the analytic conic solution as an initial guess and falls back to a
// bracketed scan, so surfaces with several ray crossings (steep, non-monotonic
// aspheres) are handled correctly: the smallest valid t > tMin wins.
type EvenAsphere struct {
	R      float64
	K      float64
	Coeffs []float64 // A4, A6, A8, ... (mm^(1-2i))
}

// Sag implements Profile.
func (p EvenAsphere) Sag(r float64) float64 {
	z := conicBaseSag(r, p.R, p.K)
	if math.IsNaN(z) {
		return z
	}
	r2 := r * r
	// Horner evaluation in r², starting at the lowest coefficient.
	var s float64
	rp := r2 * r2 // r⁴
	for i := len(p.Coeffs) - 1; i >= 0; i-- {
		s += p.Coeffs[i]
		if i > 0 {
			s *= r2
		}
	}
	return z + s*rp
}

// Slope implements Profile.
func (p EvenAsphere) Slope(r float64) float64 {
	s := conicBaseSlope(r, p.R, p.K)
	if math.IsNaN(s) {
		return s
	}
	r2 := r * r
	// d/dr Σ A_i r^(2i+4) = Σ (2i+4) A_i r^(2i+3) = r³ Σ (2i+4) A_i r^(2i)
	r3 := r2 * r
	var acc float64
	for i := len(p.Coeffs) - 1; i >= 0; i-- {
		acc += float64(2*i+4) * p.Coeffs[i]
		if i > 0 {
			acc *= r2
		}
	}
	return s + r3*acc
}

// Intersect implements Profile.
func (p EvenAsphere) Intersect(o, d Vec3, tMin, rMax float64) (float64, bool) {
	hasBase := !math.IsInf(p.R, 0) && p.R != 0
	return solveIntersection(p, p.R, p.K, hasBase, o, d, tMin, rMax)
}

// Curvature implements Profile: the paraxial vertex curvature 1/R.
func (p EvenAsphere) Curvature() float64 { return conicCurvature(p.R) }

// String implements Profile.
func (p EvenAsphere) String() string {
	return fmt.Sprintf("even asphere R=%g K=%g A4..A%d (%d coeffs)",
		p.R, p.K, 2*len(p.Coeffs)+2, len(p.Coeffs))
}

// ---------------------------------------------------------------------------
// Intersection machinery
// ---------------------------------------------------------------------------

// sagResidual evaluates g(t) = z(t) - Sag(r(t)) and its derivative dg/dt for
// the ray o + t·d (d unit).
func sagResidual(p Profile, o, d Vec3, t float64) (g, gd float64) {
	x := o.X + t*d.X
	y := o.Y + t*d.Y
	z := o.Z + t*d.Z
	r := math.Hypot(x, y)
	sag := p.Sag(r)
	if math.IsNaN(sag) {
		return math.NaN(), math.NaN()
	}
	var drdt float64
	if r > 0 {
		drdt = (x*d.X + y*d.Y) / r
	}
	return z - sag, d.Z - p.Slope(r)*drdt
}

// newtonRefine polishes t towards a root of the sag residual. It reports ok
// only when the residual is at convergence level for the local surface scale.
func newtonRefine(p Profile, o, d Vec3, t float64) (float64, bool) {
	for i := 0; i < 80; i++ {
		g, gd := sagResidual(p, o, d, t)
		if math.IsNaN(g) {
			return 0, false
		}
		z := o.Z + t*d.Z
		r := math.Hypot(o.X+t*d.X, o.Y+t*d.Y)
		tol := 1e-13 * (1 + math.Abs(z) + r)
		if math.Abs(g) <= tol {
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
	g, _ := sagResidual(p, o, d, t)
	z := o.Z + t*d.Z
	r := math.Hypot(o.X+t*d.X, o.Y+t*d.Y)
	return t, math.Abs(g) <= 1e-10*(1+math.Abs(z)+r)
}

// conicRoots returns the real roots of the base conicoid equation
//
//	r² + (1+K) z² - 2 R z = 0 ,  x = o.x + t d.x , ... ,
//
// sorted ascending. It is used as an analytic initial guess for Conic and
// EvenAsphere intersections.
func conicRoots(R, K float64, o, d Vec3) []float64 {
	if math.IsInf(R, 0) || R == 0 {
		return nil
	}
	a := d.X*d.X + d.Y*d.Y + (1+K)*d.Z*d.Z
	b := 2 * (o.X*d.X + o.Y*d.Y + (1+K)*o.Z*d.Z - R*d.Z)
	c := o.X*o.X + o.Y*o.Y + (1+K)*o.Z*o.Z - 2*R*o.Z
	if math.Abs(a) < 1e-300 {
		// Degenerate: the ray is a straight line in the vertex tangent plane.
		if b == 0 {
			return nil
		}
		return []float64{-c / b}
	}
	disc := b*b - 4*a*c
	if disc < 0 {
		return nil
	}
	sq := math.Sqrt(disc)
	t1 := (-b - sq) / (2 * a)
	t2 := (-b + sq) / (2 * a)
	if t1 > t2 {
		t1, t2 = t2, t1
	}
	if t1 == t2 {
		return []float64{t1}
	}
	return []float64{t1, t2}
}

// sagRange samples the sag of a profile between r = 0 and rMax and returns the
// smallest and largest finite values. Unbounded (rMax <= 0 or +Inf) requests
// return an infinite interval.
func sagRange(p Profile, rMax float64) (float64, float64) {
	if rMax <= 0 || math.IsInf(rMax, 1) {
		return math.Inf(-1), math.Inf(1)
	}
	zmin, zmax := math.Inf(1), math.Inf(-1)
	const n = 128
	for i := 0; i <= n; i++ {
		r := rMax * float64(i) / float64(n)
		z := p.Sag(r)
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
	if zmin > zmax {
		return math.Inf(-1), math.Inf(1)
	}
	return zmin, zmax
}

// searchWindow derives a conservative [tA, tB] window that contains every
// intersection of the ray with the profile that lies inside the radial search
// radius. It returns ok = false when no useful window exists.
func searchWindow(p Profile, o, d Vec3, tMin, rMax float64) (float64, float64, bool) {
	rSearch := rMax
	if !(rSearch > 0) || math.IsInf(rSearch, 1) {
		// Unbounded: search around the ray's own radial extent.
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
	// z crossings of the sag bounding box (padded).
	zmin, zmax := sagRange(p, rSearch)
	if !math.IsInf(zmin, -1) && !math.IsInf(zmax, 1) {
		pad := 0.5*(zmax-zmin) + 1e-6
		zmin -= pad
		zmax += pad
		if math.Abs(d.Z) > 1e-12 {
			add((zmin - o.Z) / d.Z)
			add((zmax - o.Z) / d.Z)
		}
	}
	// radial crossings |o_perp + t d_perp| = rSearch
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

// scanIntersect brackets the first sign change of the sag residual inside the
// search window, bisects it and polishes with Newton. It is the fallback for
// profiles without an analytic solution.
func scanIntersect(p Profile, o, d Vec3, tMin, rMax float64) (float64, bool) {
	tA, tB, ok := searchWindow(p, o, d, tMin, rMax)
	if !ok {
		return 0, false
	}
	const n = 400
	prevT := tA
	prevG, _ := sagResidual(p, o, d, prevT)
	for i := 1; i <= n; i++ {
		t := tA + (tB-tA)*float64(i)/float64(n)
		g, _ := sagResidual(p, o, d, t)
		if !math.IsNaN(g) && !math.IsNaN(prevG) && prevG*g <= 0 {
			lo, hi := prevT, t
			glo := prevG
			for j := 0; j < 120 && hi-lo > 1e-15*(1+math.Abs(lo)); j++ {
				mid := 0.5 * (lo + hi)
				gm, _ := sagResidual(p, o, d, mid)
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
			if rt, ok := newtonRefine(p, o, d, root); ok && rt >= tMin {
				return rt, true
			}
			if root >= tMin {
				return root, true
			}
		}
		prevT, prevG = t, g
	}
	return 0, false
}

// solveIntersection finds the smallest t > tMin with o + t·d on the profile p.
// When hasBase is true and baseR/baseK describe a conicoid, its analytic roots
// are tried first as Newton seeds.
func solveIntersection(p Profile, baseR, baseK float64, hasBase bool, o, d Vec3, tMin, rMax float64) (float64, bool) {
	if o.IsNaN() || d.IsNaN() {
		return 0, false
	}
	if hasBase {
		roots := conicRoots(baseR, baseK, o, d)
		best := math.Inf(1)
		found := false
		for _, t0 := range roots {
			if t0 < tMin {
				continue
			}
			if t, ok := newtonRefine(p, o, d, t0); ok && t >= tMin && t < best {
				best, found = t, true
			}
		}
		if found {
			return best, true
		}
	}
	return scanIntersect(p, o, d, tMin, rMax)
}

// ---------------------------------------------------------------------------
// Surface: profile + aperture + interaction type
// ---------------------------------------------------------------------------

// Surface is one optical surface of a sequential system.
//
// The surface lives in its own local frame: the vertex sits at the local
// origin and the local optical axis is the surface's z axis. A System places
// each surface on the chain of frames (see the package documentation and
// CoordBreak): with no poses and no break elements the chain is the global
// axis and every surface sits at the sum of the preceding thicknesses.
//
// Pose decentres and tilts this surface alone (the chain is restored for the
// following surfaces: auto-return). Break marks a pure coordinate-break element
// of the sequence: it has no profile and no optical action, its Pose is the
// break, and it changes the chain for every following element until a return
// break restores it.
//
// Mirror selects the interaction: when false the surface refracts between the
// medium before and after it (Snell's law, with total internal reflection
// detection); when true it reflects and the medium does not change.
type Surface struct {
	Profile  Profile  // sag description, must not be nil unless Break is set
	Aperture Aperture // clear aperture, nil = unlimited
	Mirror   bool     // true: reflect instead of refract
	Name     string   // optional label used in traces and reports
	// Pose decentres and tilts the surface (a coordinate break applied to this
	// surface only). The zero value leaves it on the chain axis.
	Pose CoordBreak
	// Break turns the element into a pure coordinate break: no profile, no
	// interaction with rays, no Step in traces; Pose is its break. The media on
	// both sides of a break element must be the same medium.
	Break bool
}

// NewSurface returns a surface with the given profile and a circular clear
// aperture of radius semiAperture (mm, <= 0 for unlimited).
func NewSurface(p Profile, semiAperture float64) Surface {
	s := Surface{Profile: p}
	if semiAperture > 0 {
		s.Aperture = CircleAperture(semiAperture)
	}
	return s
}

// AsMirror returns a copy of the surface that reflects instead of refracting.
func (s Surface) AsMirror() Surface {
	s.Mirror = true
	return s
}

// WithAperture returns a copy of the surface with the given clear aperture.
func (s Surface) WithAperture(a Aperture) Surface { s.Aperture = a; return s }

// WithName returns a copy of the surface carrying the given label.
func (s Surface) WithName(n string) Surface { s.Name = n; return s }

// Label returns the surface name, falling back to the profile description.
func (s Surface) Label() string {
	if s.Name != "" {
		return s.Name
	}
	if s.Break {
		return "coordinate break"
	}
	if s.Profile == nil {
		return "unnamed surface"
	}
	return s.Profile.String()
}

// Sag returns the local sag z at (x, y).
func (s Surface) Sag(x, y float64) float64 {
	if s.Profile == nil {
		return math.NaN()
	}
	return s.Profile.Sag(math.Hypot(x, y))
}

// SagR returns the local sag z at radius r.
func (s Surface) SagR(r float64) float64 {
	if s.Profile == nil {
		return math.NaN()
	}
	return s.Profile.Sag(r)
}

// SlopeR returns dz/dr at radius r.
func (s Surface) SlopeR(r float64) float64 {
	if s.Profile == nil {
		return math.NaN()
	}
	return s.Profile.Slope(r)
}

// Normal returns the unit surface normal at the local point (x, y), oriented
// towards +z (positive z component). With f'(r) = dz/dr and r = sqrt(x²+y²),
//
//	∂z/∂x = f'(r) x/r ,   ∂z/∂y = f'(r) y/r ,
//	n ∝ (-∂z/∂x, -∂z/∂y, 1) .
//
// At the vertex (r = 0) the normal is (0, 0, 1). The normal is analytic: no
// finite differences are involved.
func (s Surface) Normal(x, y float64) Vec3 {
	if s.Profile == nil {
		return Vec3{}
	}
	r := math.Hypot(x, y)
	if r == 0 {
		return Vec3{0, 0, 1}
	}
	sl := s.Profile.Slope(r)
	n, ok := Vec3{-sl * x / r, -sl * y / r, 1}.Unit()
	if !ok {
		return Vec3{}
	}
	return n
}

// NormalR returns the unit normal at radius r in the meridional plane
// (a positive r lies along +y): the unit normal has components
// (-f'(r)/sqrt(1+f'²), 0, 1/sqrt(1+f'²)).
func (s Surface) NormalR(r float64) Vec3 { return s.Normal(0, r) }

// MaxRadius returns the aperture's maximum radial extent (mm), +Inf when the
// aperture is unlimited.
func (s Surface) MaxRadius() float64 { return apertureOf(s.Aperture).MaxRadius() }

// Vignetted reports whether the local point (x, y) lies outside the clear
// aperture.
func (s Surface) Vignetted(x, y float64) bool {
	return !apertureOf(s.Aperture).Contains(x, y)
}

// Intersect returns the first intersection of the ray (o, d) — d unit, o in
// the surface's local frame — with the surface, using the profile's solver
// with the aperture radius as the search hint. The returned point is local.
func (s Surface) Intersect(o, d Vec3, tMin float64) (Vec3, float64, bool) {
	if s.Profile == nil {
		return Vec3{}, 0, false
	}
	t, ok := s.Profile.Intersect(o, d, tMin, s.MaxRadius())
	if !ok {
		return Vec3{}, 0, false
	}
	return o.Add(d.Mul(t)), t, true
}

// ZRange returns the smallest and largest z the surface reaches over the clear
// aperture (or over radius rMax when the aperture is unlimited and rMax > 0).
// It is used to place launch planes and search windows.
func (s Surface) ZRange(rMax float64) (float64, float64) {
	if s.Profile == nil {
		return 0, 0
	}
	r := s.MaxRadius()
	if math.IsInf(r, 1) {
		r = rMax
	}
	return sagRange(s.Profile, r)
}

// String implements fmt.Stringer.
func (s Surface) String() string {
	kind := "refract"
	if s.Mirror {
		kind = "mirror"
	}
	return fmt.Sprintf("%s [%s, %s]", s.Label(), s.Profile, kind)
}
