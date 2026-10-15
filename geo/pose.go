package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Chief rays of decentred and tilted systems
// ---------------------------------------------------------------------------

// chiefLaunch returns the chief ray of the sampling of field f on the pupil
// pup at wavelength lambdaNm:
//
//   - for a coaxial system (no coordinate break, no pose) the pupil-centre ray,
//     LaunchRay with the zero PupilPoint, as before;
//   - for a system with coordinate breaks or posed surfaces the real chief ray
//     (ChiefRay): the pupil centre is not the chief ray of such a layout, so
//     the reference of every analysis routine is taken from the real chief.
func (s *System) chiefLaunch(f Field, pup Pupil, lambdaNm float64) Ray {
	if s.hasPose() {
		return s.ChiefRay(f, lambdaNm)
	}
	return LaunchRay(f, pup, PupilPoint{}, lambdaNm)
}

// chiefRayPosed finds the chief ray of a system with coordinate breaks or posed
// surfaces: the ray of the field whose trace crosses the centre of the
// aperture stop (the vertex of the stop element, local transverse coordinates
// (0, 0)). The unknowns are the transverse launch point (x0, y0) of the ray on a
// plane in front of the system (for an object at infinity the direction is the
// field direction; for a finite object the ray starts at the object point). The
// residual is the local transverse position of the trace on the stop surface;
// it is solved by Newton's method with a forward-difference Jacobian, a fixed
// iteration count and the straight line through the stop centre as the
// starting point, so the result is deterministic.
//
// When the trace cannot reach the stop or the iteration does not converge to
// 1e-12 mm the returned ray has Status StatusNotConverged (or StatusDegenerate
// for a degenerate field) and must not be used as a reference.
func (s *System) chiefRayPosed(f Field, lambdaNm float64) Ray {
	stop := s.StopIndex
	if stop < 0 || stop >= len(s.Surfaces) {
		stop = 0
	}
	stopFrame := s.elementFrame(stop)
	centre := stopFrame.o
	zmin := s.launchPlaneZ()

	var launch func(x, y float64) Ray
	var x0, y0 float64
	if !f.Finite {
		d, ok := f.Direction()
		if !ok || d.Z == 0 {
			return Ray{WavelengthNm: lambdaNm, Status: StatusDegenerate}
		}
		zL := zmin - 1
		launch = func(x, y float64) Ray { return NewRay(V(x, y, zL), d, lambdaNm) }
		t := (zL - centre.Z) / d.Z
		x0, y0 = centre.X+d.X*t, centre.Y+d.Y*t
	} else {
		o := f.ObjectPoint()
		if o.IsNaN() || !(o.Z < zmin) || centre.Z == o.Z {
			return Ray{WavelengthNm: lambdaNm, Status: StatusDegenerate}
		}
		// The launch plane lies between the object and the first surface.
		zL := 0.5 * (o.Z + zmin)
		launch = func(x, y float64) Ray { return RayThrough(o, V(x, y, zL), lambdaNm) }
		t := (zL - o.Z) / (centre.Z - o.Z)
		x0, y0 = o.X+(centre.X-o.X)*t, o.Y+(centre.Y-o.Y)*t
	}

	// residual returns the local transverse position of the trace on the stop.
	residual := func(x, y float64) (float64, float64, bool) {
		res := s.Trace(launch(x, y))
		for _, st := range res.Steps {
			if st.Surface == stop {
				h := stopFrame.toLocal(st.Point)
				return h.X, h.Y, true
			}
		}
		return 0, 0, false
	}

	const (
		maxIter = 60
		tol     = 1e-12
		step    = 1e-7
	)
	converged := false
	for it := 0; it < maxIter; it++ {
		gx, gy, ok := residual(x0, y0)
		if !ok {
			break
		}
		if math.Hypot(gx, gy) <= tol {
			converged = true
			break
		}
		ax, ay, ok1 := residual(x0+step, y0)
		bx, by, ok2 := residual(x0, y0+step)
		if !ok1 || !ok2 {
			break
		}
		j00, j10 := (ax-gx)/step, (ay-gy)/step
		j01, j11 := (bx-gx)/step, (by-gy)/step
		det := j00*j11 - j01*j10
		if det == 0 || math.IsNaN(det) || math.IsInf(det, 0) {
			break
		}
		x0 -= (j11*gx - j01*gy) / det
		y0 -= (-j10*gx + j00*gy) / det
	}
	if !converged {
		if gx, gy, ok := residual(x0, y0); ok && math.Hypot(gx, gy) <= tol {
			converged = true
		}
	}
	ray := launch(x0, y0)
	if !converged {
		ray.Status = StatusNotConverged
	}
	return ray
}

// launchPlaneZ returns the axial position of the plane in front of every
// surface that the chief-ray search launches from: the lowest z over the clear
// apertures, or over a 10 mm radius for surfaces with unlimited apertures.
func (s *System) launchPlaneZ() float64 {
	zmin, _ := s.BoundsZ(0)
	if math.IsInf(zmin, 0) || math.IsNaN(zmin) {
		zmin, _ = s.BoundsZ(10)
	}
	return zmin
}

// ChiefCrossing returns the transverse position (x, y, zPlane) at which the real
// chief ray of the field crosses the plane z = zPlane, which must lie in the
// object space in front of the system. It is the pupil offset that centres a
// sampling of the pupil plane zPlane on the real chief ray: set
// Pupil{Offset: ChiefCrossing(...)}. It is the replacement of ChiefOffset for
// systems with coordinate breaks or posed surfaces, which need a wavelength.
func (s *System) ChiefCrossing(f Field, zPlane, lambdaNm float64) (Vec3, error) {
	if err := s.Validate(); err != nil {
		return Zero, err
	}
	if !finite(zPlane) {
		return Zero, fmt.Errorf("geo: chief crossing: plane z is not finite")
	}
	zmin := s.launchPlaneZ()
	if !(zPlane < zmin) {
		return Zero, fmt.Errorf("geo: chief crossing: plane z=%g is not in front of the system (z=%g)", zPlane, zmin)
	}
	ray := s.ChiefRay(f, lambdaNm)
	if ray.Status != StatusOK {
		return Zero, fmt.Errorf("geo: chief crossing: the real chief ray did not converge (%s)", ray.Status)
	}
	if ray.Dir.Z == 0 {
		return Zero, fmt.Errorf("geo: chief crossing: the chief ray is parallel to the plane")
	}
	t := (zPlane - ray.Origin.Z) / ray.Dir.Z
	p := ray.At(t)
	return V(p.X, p.Y, zPlane), nil
}
