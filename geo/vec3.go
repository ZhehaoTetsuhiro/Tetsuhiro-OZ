package geo

import (
	"fmt"
	"math"
)

// Vec3 is a three-component vector in millimetres in the coordinate frame of
// package geo (x right, y up, z along the optical axis, right-handed).
//
// The zero value is the origin. Vec3 values are immutable in practice; every
// operation returns a new value.
type Vec3 struct {
	X, Y, Z float64
}

// V is shorthand for Vec3{X: x, Y: y, Z: z}.
func V(x, y, z float64) Vec3 { return Vec3{X: x, Y: y, Z: z} }

// Zero is the origin / null vector.
var Zero = Vec3{}

// AxisZ is the unit vector along the optical axis (+z).
var AxisZ = Vec3{0, 0, 1}

// AxisX and AxisY are the unit vectors along +x and +y.
var (
	AxisX = Vec3{1, 0, 0}
	AxisY = Vec3{0, 1, 0}
)

// Add returns v + w.
func (v Vec3) Add(w Vec3) Vec3 { return Vec3{v.X + w.X, v.Y + w.Y, v.Z + w.Z} }

// Sub returns v - w.
func (v Vec3) Sub(w Vec3) Vec3 { return Vec3{v.X - w.X, v.Y - w.Y, v.Z - w.Z} }

// Mul returns the scalar multiple s*v.
func (v Vec3) Mul(s float64) Vec3 { return Vec3{v.X * s, v.Y * s, v.Z * s} }

// Div returns v/s. Division by zero yields ±Inf/NaN as in IEEE arithmetic.
func (v Vec3) Div(s float64) Vec3 { return Vec3{v.X / s, v.Y / s, v.Z / s} }

// Neg returns -v.
func (v Vec3) Neg() Vec3 { return Vec3{-v.X, -v.Y, -v.Z} }

// Dot returns the scalar (inner) product v·w.
func (v Vec3) Dot(w Vec3) float64 { return v.X*w.X + v.Y*w.Y + v.Z*w.Z }

// Cross returns the vector (outer) product v×w. In the right-handed frame of
// this package the components are
//
//	(v×w).x = v.y*w.z - v.z*w.y
//	(v×w).y = v.z*w.x - v.x*w.z
//	(v×w).z = v.x*w.y - v.y*w.x .
func (v Vec3) Cross(w Vec3) Vec3 {
	return Vec3{
		v.Y*w.Z - v.Z*w.Y,
		v.Z*w.X - v.X*w.Z,
		v.X*w.Y - v.Y*w.X,
	}
}

// Norm2 returns the squared Euclidean norm v·v.
func (v Vec3) Norm2() float64 { return v.Dot(v) }

// Norm returns the Euclidean norm |v|.
func (v Vec3) Norm() float64 { return math.Sqrt(v.Dot(v)) }

// Normalized returns v/|v|. If |v| is zero (or not finite) the zero vector is
// returned, so callers that need to detect the degenerate case should use
// Unit instead.
func (v Vec3) Normalized() Vec3 {
	n := v.Norm()
	if n == 0 || math.IsNaN(n) || math.IsInf(n, 0) {
		return Vec3{}
	}
	return v.Div(n)
}

// Unit returns v/|v| and reports whether the result is meaningful: ok is false
// when |v| is not a strictly positive finite number.
func (v Vec3) Unit() (Vec3, bool) {
	n := v.Norm()
	if !(n > 0) || math.IsInf(n, 0) || math.IsNaN(n) {
		return Vec3{}, false
	}
	return v.Div(n), true
}

// Perp returns the component of v perpendicular to the optical axis, (x, y, 0).
func (v Vec3) Perp() Vec3 { return Vec3{v.X, v.Y, 0} }

// PerpNorm returns the distance from the optical axis, sqrt(x²+y²).
func (v Vec3) PerpNorm() float64 { return math.Hypot(v.X, v.Y) }

// Distance returns |v - w|.
func (v Vec3) Distance(w Vec3) float64 { return v.Sub(w).Norm() }

// Lerp returns the linear interpolation (1-t)*v + t*w.
func Lerp(v, w Vec3, t float64) Vec3 {
	return Vec3{
		v.X + (w.X-v.X)*t,
		v.Y + (w.Y-v.Y)*t,
		v.Z + (w.Z-v.Z)*t,
	}
}

// IsNaN reports whether any component is NaN.
func (v Vec3) IsNaN() bool {
	return math.IsNaN(v.X) || math.IsNaN(v.Y) || math.IsNaN(v.Z)
}

// Equal reports whether v and w agree component-wise within the absolute
// tolerance tol.
func (v Vec3) Equal(w Vec3, tol float64) bool {
	return math.Abs(v.X-w.X) <= tol && math.Abs(v.Y-w.Y) <= tol && math.Abs(v.Z-w.Z) <= tol
}

// Reflect returns the direction of the mirror reflection of the unit direction
// d off a surface with unit normal n (n need not be oriented towards d; the
// result is independent of sign(n)):
//
//	d' = d - 2 (d·n) n .
//
// The result is a unit vector when d and n are unit vectors.
func Reflect(d, n Vec3) Vec3 {
	return d.Sub(n.Mul(2 * d.Dot(n)))
}

// Refract applies Snell's law in vector form to the unit incident direction d
// and the unit surface normal n, for a ray travelling from a medium of index
// n1 into a medium of index n2. n may point either way with respect to d (the
// sign is detected from d·n), so callers may pass the +z oriented normal of a
// Surface directly.
//
// With i the incidence angle (cos i = -d·n for the correctly oriented normal)
// the transmitted direction is
//
//	d' = (n1/n2) d + ( (n1/n2) cos i - cos t ) n ,
//	cos t = sqrt(1 - (n1/n2)² (1 - cos²i)) .
//
// ok is false when 1 - (n1/n2)²(1-cos²i) < 0, i.e. for total internal
// reflection; the returned vector is then the zero vector. Grazing
// transmission (cos t = 0) is returned as ok with a direction in the tangent
// plane.
func Refract(d, n Vec3, n1, n2 float64) (Vec3, bool) {
	// Orient the normal against the incoming ray so that cos i >= 0.
	nn := n
	cosi := -d.Dot(nn)
	if cosi < 0 {
		nn = nn.Neg()
		cosi = -cosi
	}
	if n2 == 0 {
		return Vec3{}, false
	}
	eta := n1 / n2
	k := 1 - eta*eta*(1-cosi*cosi)
	if k < 0 {
		return Vec3{}, false
	}
	cost := math.Sqrt(k)
	out := d.Mul(eta).Add(nn.Mul(eta*cosi - cost))
	u, ok := out.Unit()
	if !ok {
		return Vec3{}, false
	}
	return u, true
}

// IncidenceAngle returns the angle of incidence (in radians, 0 <= θ <= π/2) of
// the unit direction d on the surface with unit normal n, independent of the
// orientation of n.
func IncidenceAngle(d, n Vec3) float64 {
	c := math.Abs(d.Dot(n))
	if c > 1 {
		c = 1
	}
	return math.Acos(c)
}

// String formats a vector as "(x, y, z)" with %g precision.
func (v Vec3) String() string {
	return fmt.Sprintf("(%g, %g, %g)", v.X, v.Y, v.Z)
}
