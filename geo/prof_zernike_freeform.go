package geo

import "math"

// FreeformProfile wiring for ZernikeSurface. The two-dimensional dispatch of
// Surface/System (Surface.Sag, Surface.Normal, Surface.Intersect,
// Surface.ZRange) keys on the FreeformProfile interface, which the polynomial
// freeform implements directly. A Zernike surface is also a genuine function
// of (x, y) — its sag and gradient are evaluated from the two-dimensional
// Zernike sum, not from a radial surrogate — so it satisfies the same
// interface. Without these methods an asymmetric Zernike surface (|m| > 0)
// would silently fall back to the θ = 0 cut in Surface.Normal/Sag and evade
// the freeform rejection of the first-order engines.
var _ FreeformProfile = ZernikeSurface{}

// IntersectXY implements FreeformProfile. ZernikeSurface.Intersect already
// solves the full two-dimensional residual (x, y), so this is the same solver.
func (p ZernikeSurface) IntersectXY(o, d Vec3, tMin, rMax float64) (float64, bool) {
	return p.Intersect(o, d, tMin, rMax)
}

// ZRangeXY implements FreeformProfile: the sag range over the aperture disc of
// radius rMax, sampled on a deterministic polar grid (the same scheme as the
// polynomial freeform's ZRangeXY; Zernike polynomials of the orders used for
// optical surfaces are smooth, and the interval is only a conservative
// bracketing bound for the intersection search). A non-positive or +Inf radius
// means "unbounded" and yields an infinite interval.
func (p ZernikeSurface) ZRangeXY(rMax float64) (float64, float64) {
	if !(rMax > 0) || math.IsInf(rMax, 1) {
		return math.Inf(-1), math.Inf(1)
	}
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
	const nr, na = 16, 24
	upd(0, 0)
	for k := 1; k <= nr; k++ {
		rr := rMax * float64(k) / float64(nr)
		for j := 0; j < na; j++ {
			a := 2 * math.Pi * float64(j) / float64(na)
			upd(rr*math.Cos(a), rr*math.Sin(a))
		}
	}
	if zmin > zmax {
		return math.Inf(-1), math.Inf(1)
	}
	return zmin, zmax
}

// RotationallySymmetric implements FreeformProfile: a Zernike surface is
// rotationally symmetric exactly when only m = 0 members are active (see
// Symmetric). Asymmetric Zernike surfaces are rejected by the first-order and
// third-order engines (System.requireCoaxial) like any other freeform.
func (p ZernikeSurface) RotationallySymmetric() bool { return p.Symmetric() }
