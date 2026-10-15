package geo

import (
	"math"
	"math/cmplx"
)

// Ellipse is the polarization ellipse of a Jones vector. Its conventions are
// those of optics.StokesAt in optics/polar.go, and a test compares the two
// field by field.
type Ellipse struct {
	S0, S1, S2, S3 float64 // Stokes parameters: S0 = |Ex|²+|Ey|², S2 = 2Re(Ex*Ey), S3 = 2Im(Ex*Ey)
	Azimuth        float64 // rad, angle of the major axis against ex, in [−π/2, π/2]
	Ellipticity    float64 // rad, ½·asin(S3/S0) in [−π/4, π/4]; > 0 is left-handed
	AxisRatio      float64 // minor/major axis ratio |tan χ| in [0, 1]
	Handedness     string  // "左旋" (left), "右旋" (right) or "线偏振" (linear), as in optics
	Degree         float64 // degree of polarization sqrt(S1²+S2²+S3²)/S0; 1 for a coherent state
}

// EllipseOf returns the polarization ellipse of the Jones vector (ex, ey). S3 > 0
// means that the electric vector rotates counterclockwise as seen by a viewer
// facing the beam (the left-handed case of optics/polar.go). The rotation sense
// is checked against the time-domain field in polar_convention_test.go.
func EllipseOf(ex, ey complex128) Ellipse {
	var el Ellipse
	ax2 := real(ex)*real(ex) + imag(ex)*imag(ex)
	ay2 := real(ey)*real(ey) + imag(ey)*imag(ey)
	cr := cmplx.Conj(ex) * ey
	el.S0 = ax2 + ay2
	el.S1 = ax2 - ay2
	el.S2 = 2 * real(cr)
	el.S3 = 2 * imag(cr)
	if el.S0 <= 0 {
		return el
	}
	el.Azimuth = polWrapHalfPi(0.5 * math.Atan2(el.S2, el.S1))
	el.Ellipticity = 0.5 * math.Asin(polClamp(el.S3/el.S0, -1, 1))
	el.AxisRatio = math.Abs(math.Tan(el.Ellipticity))
	el.Degree = math.Hypot(math.Hypot(el.S1, el.S2), el.S3) / el.S0
	switch {
	case el.AxisRatio < 1e-3:
		el.Handedness = "线偏振"
	case el.S3 > 0:
		el.Handedness = "左旋"
	default:
		el.Handedness = "右旋"
	}
	return el
}

// polWrapHalfPi maps an angle into [−π/2, π/2].
func polWrapHalfPi(a float64) float64 {
	for a > math.Pi/2 {
		a -= math.Pi
	}
	for a < -math.Pi/2 {
		a += math.Pi
	}
	return a
}

// polClamp limits v to [lo, hi].
func polClamp(v, lo, hi float64) float64 {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}
