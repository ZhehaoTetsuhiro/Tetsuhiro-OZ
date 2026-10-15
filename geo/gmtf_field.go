package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Field curvature from real rays (G11)
// ---------------------------------------------------------------------------
//
// The tangential (meridional) and sagittal focal positions of a field are
// measured with real ray tracing in the paraxial limit of the entrance pupil:
// the tangential focus is where the two rays through (0, ±ε) of the pupil
// cross each other in the meridional plane, the sagittal focus is where the
// two rays through (±ε, 0) cross the meridional plane (by the symmetry of a
// rotationally symmetric system the two sagittal rays meet on that plane).
// For ε → 0 both points converge to the focal points of the chief ray. The
// result is the set of focal points (y, z) of the field, from which the
// focal-surface sags are fitted.
//
// Third-order structure (classical: the tangential image lies three times as
// far from the Petzval surface as the sagittal one, z_T − z_P = 3(z_S − z_P);
// cf. Conrady 1929 as quoted in MNRAS 342, 33 (2003), "Exact optics II";
// textbook page citation 待核对). The Petzval sag is δz_P = −S_IV·h²/(2n′u′²).
// In terms of the Seidel sums of this package, the sagittal and tangential
// focus shifts are q + a and q + 3a with q = −S4/(2NA²) (Petzval) and
// a = −S3/(2NA²) (astigmatism). The Petzval surface is therefore
// P = (3·S − T)/2, which needs only the two measured surfaces. The test in
// gmtf_field_test.go checks this against the analytic thin-lens Petzval curvature.

// FocalPoint is the tangential and sagittal focus of one field of an
// object-at-infinity system, in the meridional (y–z) plane.
type FocalPoint struct {
	ThetaDeg float64 // field angle in the meridional plane (degrees)
	YT, ZT   float64 // tangential focus (mm)
	YS, ZS   float64 // sagittal focus (mm)
}

// parallelRay launches the parallel-bundle ray of an infinite field through
// the pupil point p (EntrancePupil coordinates), starting on the plane zStart
// in front of the system so that surfaces at any z are reached. The ray is the
// one a parallel bundle carries through p, not a ray from a finite point.
func parallelRay(f Field, pup Pupil, p PupilPoint, zStart, lambdaNm float64) (Ray, error) {
	d, ok := f.Direction()
	if !ok {
		return Ray{}, fmt.Errorf("geo: field focus: the field must be at infinity")
	}
	q := V(pup.Offset.X+p.PX*pup.Radius, pup.Offset.Y+p.PY*pup.Radius, pup.Z)
	s := (pup.Z - zStart) / d.Z
	r := NewRay(q.Sub(d.Mul(s)), d, lambdaNm)
	if !r.Valid() {
		return r, fmt.Errorf("geo: field focus: degenerate launch ray")
	}
	return r, nil
}

// FieldFoci returns the tangential and sagittal foci of each field angle
// (degrees, in the meridional y–z plane) of an object-at-infinity system, using
// pupil pairs at ±rimFrac of the entrance-pupil radius (0 < rimFrac ≪ 1). The
// entrance pupil is the paraxial one (EntrancePupil); each ray is a real ray
// traced through the surfaces. It is an error if a ray is not traversed or if
// the tangential pair leaves the meridional plane, which would mean the system
// or the field is not the configuration assumed.
func (s *System) FieldFoci(thetaDeg []float64, lambdaNm, rimFrac float64) ([]FocalPoint, error) {
	if !(rimFrac > 0 && rimFrac < 0.1) {
		return nil, fmt.Errorf("geo: field foci: rimFrac %g outside (0, 0.1)", rimFrac)
	}
	pup, err := s.EntrancePupil(Field{}, lambdaNm)
	if err != nil {
		return nil, err
	}
	zmin, _ := s.BoundsZ(pup.Radius)
	zStart := zmin - 1
	out := make([]FocalPoint, len(thetaDeg))
	for i, th := range thetaDeg {
		f := FieldAngle(0, th)
		trace := func(p PupilPoint) (Vec3, Vec3, error) {
			r, err := parallelRay(f, pup, p, zStart, lambdaNm)
			if err != nil {
				return Vec3{}, Vec3{}, err
			}
			res := s.Trace(r)
			if !res.OK() {
				return Vec3{}, Vec3{}, fmt.Errorf("geo: field foci: θ=%g° ray at pupil %v: %v", th, p, res.Error())
			}
			return res.Ray.Origin, res.Ray.Dir, nil
		}
		// Tangential pair: (0, ±ε), crossing in the meridional plane x = 0.
		p1, d1, err := trace(PupilPoint{PX: 0, PY: rimFrac})
		if err != nil {
			return nil, err
		}
		p2, d2, err := trace(PupilPoint{PX: 0, PY: -rimFrac})
		if err != nil {
			return nil, err
		}
		if math.Abs(p1.X)+math.Abs(d1.X)+math.Abs(p2.X)+math.Abs(d2.X) > 1e-9 {
			return nil, fmt.Errorf("geo: field foci: tangential rays leave the meridional plane at θ=%g°", th)
		}
		// Solve p1 + t·d1 = p2 + u·d2 in the (y, z) plane.
		det := -d1.Y*d2.Z + d2.Y*d1.Z
		if det == 0 {
			return nil, fmt.Errorf("geo: field foci: parallel tangential rays at θ=%g°", th)
		}
		dy, dz := p2.Y-p1.Y, p2.Z-p1.Z
		t := (dy*(-d2.Z) - (-d2.Y)*dz) / det
		tY := p1.Y + t*d1.Y
		tZ := p1.Z + t*d1.Z

		// Sagittal pair: (±ε, 0). Each ray meets x = 0; by symmetry the two
		// crossings coincide, and their mean is the sagittal focus.
		q1, e1, err := trace(PupilPoint{PX: rimFrac, PY: 0})
		if err != nil {
			return nil, err
		}
		q2, e2, err := trace(PupilPoint{PX: -rimFrac, PY: 0})
		if err != nil {
			return nil, err
		}
		if e1.X == 0 || e2.X == 0 {
			return nil, fmt.Errorf("geo: field foci: sagittal ray parallel to the meridional plane at θ=%g°", th)
		}
		ta := -q1.X / e1.X
		tb := -q2.X / e2.X
		sY := 0.5 * ((q1.Y + ta*e1.Y) + (q2.Y + tb*e2.Y))
		sZ := 0.5 * ((q1.Z + ta*e1.Z) + (q2.Z + tb*e2.Z))
		out[i] = FocalPoint{ThetaDeg: th, YT: tY, ZT: tZ, YS: sY, ZS: sZ}
	}
	return out, nil
}
