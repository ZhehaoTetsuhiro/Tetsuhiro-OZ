package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Fields: the object points a ray bundle belongs to
// ---------------------------------------------------------------------------

// Field identifies one object point: one "field" in the classical
// object–stop–pupil description of an optical system. Following the package
// conventions two kinds of object are supported:
//
//   - object at infinity (Finite = false, the zero value): TanX and TanY are
//     the tangents of
//     the field angles, (tan θx, tan θy). A source above the axis (TanY > 0)
//     sends light travelling along d = (−TanX, −TanY, 1)/|…|, i.e. downwards
//     and towards the image (see the package documentation).
//   - finite object (Finite = true, see FieldHeight): the object point sits at
//     (HX, HY, ObjectZ) with ObjectZ in front of the first surface (whose
//     vertex is at z = 0); HX and HY are the transverse object coordinates in
//     mm.
//
// The zero value is the on-axis field of an object at infinity, so the
// tangents can be set directly in a literal: Field{TanY: 0.1} is a tilted
// infinite field.
type Field struct {
	TanX, TanY float64 // tangents of the field angles (object at infinity)
	HX, HY     float64 // object height in mm (finite object)
	ObjectZ    float64 // object plane in mm (finite object)
	Finite     bool    // true for a finite object (use FieldHeight)
}

// FieldAngle returns the field of an object at infinity built from field
// angles in degrees (the public constructors of this package take angles in
// degrees, see the package documentation).
func FieldAngle(thetaXDeg, thetaYDeg float64) Field {
	return Field{
		TanX: math.Tan(thetaXDeg * math.Pi / 180),
		TanY: math.Tan(thetaYDeg * math.Pi / 180),
	}
}

// FieldTan returns the field of an object at infinity from the tangents of
// the field angles directly (TanX = tan θx, TanY = tan θy).
func FieldTan(tanX, tanY float64) Field {
	return Field{TanX: tanX, TanY: tanY}
}

// FieldHeight returns the field of a finite object at axial position objectZ
// (mm, in front of the first surface) with transverse coordinates (hx, hy) in
// mm.
func FieldHeight(objectZ, hx, hy float64) Field {
	return Field{HX: hx, HY: hy, ObjectZ: objectZ, Finite: true}
}

// IsInfinite reports whether the field belongs to an object at infinity.
func (f Field) IsInfinite() bool { return !f.Finite }

// ObjectPlane returns the axial object position in mm, math.Inf(-1) for an
// object at infinity (the convention of System.Paraxial and
// System.EntrancePupil).
func (f Field) ObjectPlane() float64 {
	if !f.Finite {
		return math.Inf(-1)
	}
	return f.ObjectZ
}

// ObjectPoint returns the object point of a finite field; the zero vector for
// an infinite field.
func (f Field) ObjectPoint() Vec3 {
	if !f.Finite {
		return Zero
	}
	return V(f.HX, f.HY, f.ObjectZ)
}

// Direction returns the unit propagation direction of the parallel beam of an
// infinite field,
//
//	d = (−TanX, −TanY, 1) / |(−TanX, −TanY, 1)| .
//
// ok is false for a finite field (whose rays are not parallel).
func (f Field) Direction() (Vec3, bool) {
	if f.Finite {
		return Zero, false
	}
	return V(-f.TanX, -f.TanY, 1).Unit()
}

// Transverse returns the transverse field coordinates: (TanX, TanY) for an
// infinite field, (HX, HY) for a finite one.
func (f Field) Transverse() (x, y float64) {
	if f.Finite {
		return f.HX, f.HY
	}
	return f.TanX, f.TanY
}

// Label formats the field compactly.
func (f Field) Label() string {
	if f.Finite {
		return fmt.Sprintf("object h=(%.6g, %.6g) mm at z=%.6g", f.HX, f.HY, f.ObjectZ)
	}
	return fmt.Sprintf("∞ field tan=(%.6g, %.6g)", f.TanX, f.TanY)
}

// String implements fmt.Stringer.
func (f Field) String() string { return "field[" + f.Label() + "]" }

// ---------------------------------------------------------------------------
// Pupils
// ---------------------------------------------------------------------------

// Pupil identifies where a ray bundle is launched: the rays of a sampling
// cross the plane z = Z at
//
//	Offset + (PX·Radius, PY·Radius)   with |(PX, PY)| <= 1
//
// in the normalized pupil coordinates of the unit disc. It is normally the
// paraxial entrance pupil (see System.EntrancePupil), whose centre lies on the
// optical axis (Offset = 0), but any reference plane works — mirror systems,
// for which the paraxial engine is not available, simply pass a plane in front
// of the first surface with the clear semi-diameter as Radius and the chief
// offset of System.ChiefOffset as Offset, which puts the pupil centre on the
// real chief ray (the classical pupil convention that the third-order sums
// assume; for an entrance pupil the offset is zero because the chief crosses
// the stop centre, whose image is the entrance pupil centre on the axis).
type Pupil struct {
	Z      float64 // axial position of the pupil plane (mm)
	Radius float64 // pupil semi-diameter (mm)
	Offset Vec3    // transverse offset of the pupil centre in the plane (mm)
}

// String implements fmt.Stringer.
func (p Pupil) String() string {
	return fmt.Sprintf("pupil[z=%g mm r=%g mm offset=(%g, %g)]", p.Z, p.Radius, p.Offset.X, p.Offset.Y)
}

// PupilPoint is a point of the unit pupil disc used to sample a ray bundle.
// The zero value is the pupil centre (the chief ray of the sampling).
type PupilPoint struct {
	PX, PY float64
}

// PupilPointAt returns the normalized pupil point (px, py).
func PupilPointAt(px, py float64) PupilPoint { return PupilPoint{PX: px, PY: py} }

// R returns the distance of the point from the pupil centre.
func (p PupilPoint) R() float64 { return math.Hypot(p.PX, p.PY) }

// Inside reports whether the point lies in the closed unit disc. The boundary
// is accepted with a round-off tolerance so that points generated on the rim
// by the samplings are never rejected by the last bit of a cosine.
func (p PupilPoint) Inside() bool { return p.PX*p.PX+p.PY*p.PY <= 1+1e-12 }

// Angle returns the azimuth of the point in radians in [0, 2π).
func (p PupilPoint) Angle() float64 {
	a := math.Atan2(p.PY, p.PX)
	if a < 0 {
		a += 2 * math.Pi
	}
	return a
}

// String implements fmt.Stringer.
func (p PupilPoint) String() string {
	return fmt.Sprintf("p=(%.6g, %.6g)", p.PX, p.PY)
}

// ---------------------------------------------------------------------------
// Pupil samplings
// ---------------------------------------------------------------------------

// Hexapolar returns the points of a hexapolar sampling of the pupil: the
// centre plus nRings concentric rings, ring k (1..nRings) carrying 6k points
// at radius k/nRings, the first of them on the +x axis. Total count
// 1 + 3·nRings·(nRings+1). This is the standard mixed radial/azimuthal
// sampling used for spot diagrams: every ring has the same coverage in angle
// and the spacing in radius is uniform.
func Hexapolar(nRings int) []PupilPoint {
	if nRings <= 0 {
		return []PupilPoint{{}}
	}
	out := make([]PupilPoint, 0, 1+3*nRings*(nRings+1))
	out = append(out, PupilPoint{})
	for k := 1; k <= nRings; k++ {
		r := float64(k) / float64(nRings)
		n := 6 * k
		for j := 0; j < n; j++ {
			a := 2 * math.Pi * float64(j) / float64(n)
			out = append(out, PupilPoint{r * math.Cos(a), r * math.Sin(a)})
		}
	}
	return out
}

// Rings returns nRings concentric rings at radii i/nRings (i = 1..nRings),
// each carrying perRing points uniformly spaced in angle with the first point
// on the +x axis, plus the pupil centre. Total count
// 1 + nRings·perRing. Successive rings are staggered by half a step in angle
// to avoid purely radial alignments.
func Rings(nRings, perRing int) []PupilPoint {
	out := []PupilPoint{{}}
	if nRings <= 0 || perRing <= 0 {
		return out
	}
	out = make([]PupilPoint, 0, 1+nRings*perRing)
	out = append(out, PupilPoint{})
	for k := 1; k <= nRings; k++ {
		r := float64(k) / float64(nRings)
		off := math.Pi * float64(k%2) / float64(perRing)
		for j := 0; j < perRing; j++ {
			a := 2*math.Pi*float64(j)/float64(perRing) + off
			out = append(out, PupilPoint{r * math.Cos(a), r * math.Sin(a)})
		}
	}
	return out
}

// Grid returns the nx × ny points of an evenly spaced grid over the square
// [−1, 1]² (boundary included), keeping only the points of the unit disc, so
// the effective count is somewhat smaller than nx·ny. The centre is included
// when both nx and ny are odd.
func Grid(nx, ny int) []PupilPoint {
	if nx <= 0 || ny <= 0 {
		return nil
	}
	out := make([]PupilPoint, 0, nx*ny)
	for j := 0; j < ny; j++ {
		for i := 0; i < nx; i++ {
			px, py := gridCoord(i, nx), gridCoord(j, ny)
			p := PupilPoint{px, py}
			if p.Inside() {
				out = append(out, p)
			}
		}
	}
	return out
}

// gridCoord maps i in [0, n) onto [−1, 1].
func gridCoord(i, n int) float64 {
	if n == 1 {
		return 0
	}
	return 2*float64(i)/float64(n-1) - 1
}

// Fan returns n points along the pupil diameter in the direction
// azimuthDeg (degrees, 0 = the +x/sagittal direction, 90 = the +y/meridional
// direction), from −1 to +1, endpoints included. For n = 1 the pupil centre is
// returned. This is the classical "ray fan" sampling of transverse aberration
// plots.
func Fan(n int, azimuthDeg float64) []PupilPoint {
	if n <= 0 {
		return nil
	}
	if n == 1 {
		return []PupilPoint{{}}
	}
	a := azimuthDeg * math.Pi / 180
	ca, sa := math.Cos(a), math.Sin(a)
	// Snap the components of the cardinal directions to exact zeros so that
	// the meridional and sagittal fans are exactly on-axis (as their
	// definition requires), not ±6e-17 off.
	if math.Abs(ca) < 1e-15 {
		ca = 0
	}
	if math.Abs(sa) < 1e-15 {
		sa = 0
	}
	out := make([]PupilPoint, n)
	for j := 0; j < n; j++ {
		t := 2*float64(j)/float64(n-1) - 1
		out[j] = PupilPoint{t * ca, t * sa}
	}
	return out
}

// MeridionalFan returns n points along the meridional (y) diameter of the
// pupil, from (0, −1) to (0, +1): the fan used for tangential transverse
// aberration.
func MeridionalFan(n int) []PupilPoint { return Fan(n, 90) }

// SagittalFan returns n points along the sagittal (x) diameter of the pupil,
// from (−1, 0) to (+1, 0): the fan used for sagittal transverse aberration.
func SagittalFan(n int) []PupilPoint { return Fan(n, 0) }

// CrossFan returns the union of the meridional and the sagittal fan. The two
// fans share only the pupil centre, so the count is 2n−1 for n ≥ 2.
func CrossFan(n int) []PupilPoint {
	out := MeridionalFan(n)
	for _, p := range SagittalFan(n) {
		if p == (PupilPoint{}) {
			continue // the centre is already present
		}
		out = append(out, p)
	}
	return out
}

// ---------------------------------------------------------------------------
// Launching rays
// ---------------------------------------------------------------------------

// LaunchRay returns the ray of the given field that crosses the pupil plane
// pup at the normalized pupil point p:
//
//   - object at infinity: all rays of the field are parallel, so the ray
//     starts at the pupil point and travels along the field direction
//     (−TanX, −TanY, 1)/|…|;
//   - finite object: the ray starts at the object point
//     (HX, HY, ObjectZ) and travels towards the pupil point.
//
// The pupil point p = (0, 0) is the chief ray of the sampling. A degenerate
// configuration (zero direction, non-finite data) is reported through the
// status of the returned ray (StatusDegenerate), following the status model of
// the package.
func LaunchRay(f Field, pup Pupil, p PupilPoint, lambdaNm float64) Ray {
	if !finite(pup.Z) || !finite(pup.Radius) || !finite(p.PX) || !finite(p.PY) ||
		!finite(pup.Offset.X) || !finite(pup.Offset.Y) {
		return Ray{WavelengthNm: lambdaNm, Status: StatusDegenerate}
	}
	q := V(pup.Offset.X+p.PX*pup.Radius, pup.Offset.Y+p.PY*pup.Radius, pup.Z)
	if !f.Finite {
		d, ok := f.Direction()
		if !ok {
			return Ray{WavelengthNm: lambdaNm, Status: StatusDegenerate}
		}
		return NewRay(q, d, lambdaNm)
	}
	o := f.ObjectPoint()
	if !finite(o.X) || !finite(o.Y) || !finite(o.Z) {
		return Ray{WavelengthNm: lambdaNm, Status: StatusDegenerate}
	}
	return RayThrough(o, q, lambdaNm)
}

// LaunchSet launches the whole sampling (field, pupil, points) as rays, in the
// same order as the points.
func LaunchSet(f Field, pup Pupil, pts []PupilPoint, lambdaNm float64) []Ray {
	out := make([]Ray, len(pts))
	for i, p := range pts {
		out[i] = LaunchRay(f, pup, p, lambdaNm)
	}
	return out
}

// finite reports whether x is a finite number.
func finite(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }

// ---------------------------------------------------------------------------
// System-level pupil helpers
// ---------------------------------------------------------------------------

// EntrancePupil returns the paraxial entrance pupil of the system for the
// object position implied by the field f (the object position does not affect
// the entrance pupil; it is accepted so that the call reads like the rest of
// the sampling API), in the form of a Pupil plane: its axial position
// (EntrancePupilZ) and its semi-diameter (EntrancePupilRadius).
//
// The entrance pupil is the image of the aperture stop in the object space; it
// is derived from the first-order model, so the system must be refractive (the
// first-order engine does not handle mirrors; use an explicit Pupil for mirror
// systems). Systems whose entrance pupil lies at infinity (telecentric object
// space) or whose stop has no finite aperture are reported as errors.
func (s *System) EntrancePupil(f Field, lambdaNm float64) (Pupil, error) {
	p, err := s.Paraxial(lambdaNm, f.ObjectPlane())
	if err != nil {
		return Pupil{}, err
	}
	if math.IsInf(p.EntrancePupilZ, 0) {
		return Pupil{}, fmt.Errorf("geo: entrance pupil at infinity (telecentric object space)")
	}
	if !finitePositive(p.EntrancePupilRadius) {
		return Pupil{}, fmt.Errorf("geo: aperture stop has no finite semi-diameter; cannot size the entrance pupil")
	}
	return Pupil{Z: p.EntrancePupilZ, Radius: p.EntrancePupilRadius}, nil
}

// ExitPupil returns the paraxial exit pupil of the system as a Pupil plane
// (axial position and semi-diameter). See EntrancePupil for the restrictions.
func (s *System) ExitPupil(f Field, lambdaNm float64) (Pupil, error) {
	p, err := s.Paraxial(lambdaNm, f.ObjectPlane())
	if err != nil {
		return Pupil{}, err
	}
	if math.IsInf(p.ExitPupilZ, 0) {
		return Pupil{}, fmt.Errorf("geo: exit pupil at infinity (telecentric image space)")
	}
	if !finitePositive(p.ExitPupilRadius) {
		return Pupil{}, fmt.Errorf("geo: aperture stop has no finite semi-diameter; cannot size the exit pupil")
	}
	return Pupil{Z: p.ExitPupilZ, Radius: p.ExitPupilRadius}, nil
}

// ChiefRay returns the real chief ray of the field: the ray that crosses the
// centre of the aperture stop (the stop vertex at
// (0, 0, VertexZ(StopIndex))).
//
// The ray is found by deterministic ray aiming: the meridional equivalent of
// the field (its meridional magnitude) is used and the launch height in front
// of the system is refined with a fixed-count secant iteration on the
// condition "the traced ray hits the stop vertex"; the solution is then
// rotated about the axis to the azimuth of the actual field. When the trial
// rays cannot be evaluated (they do not reach the stop) the paraxial chief
// ray — the straight line from the object (or along the field direction)
// through the stop vertex — is returned instead, so a result is always
// available. The returned ray starts in front of the system (1 mm before the
// frontmost sag) and can be traced directly.
//
// Note the difference from the chief ray of a pupil sampling (the ray through
// the centre of the paraxial entrance pupil, LaunchRay with the zero
// PupilPoint): the two agree to first order but the real chief ray defined
// here is the one whose image height defines the classical distortion.
func (s *System) ChiefRay(f Field, lambdaNm float64) Ray {
	if s.hasPose() {
		return s.chiefRayPosed(f, lambdaNm)
	}
	stop := s.StopIndex
	if stop < 0 || stop >= len(s.Surfaces) {
		stop = 0
	}
	zStop := s.VertexZ(stop)
	zmin, _ := s.BoundsZ(0)
	zL := zmin - 1 // launch plane, in front of every sag

	// Meridional equivalent of the field: magnitude h along +y, azimuth φ.
	tx, ty := f.Transverse()
	hMag := math.Hypot(tx, ty)
	azimuth := math.Atan2(ty, tx)
	meridional := func(y0 float64) Ray {
		if f.Finite {
			return RayThrough(V(0, hMag, f.ObjectZ), V(0, y0, zL), lambdaNm)
		}
		d := V(0, -hMag, 1)
		return NewRay(V(0, y0, zL), d, lambdaNm)
	}
	// Paraxial starting height: the straight line through the stop vertex.
	y0p := hMag * (zStop - zL)
	if f.Finite {
		// Line from (0, hMag, ObjectZ) through (0, 0, zStop), at z = zL.
		y0p = hMag + (0-hMag)*(zL-f.ObjectZ)/(zStop-f.ObjectZ)
	}
	if hMag == 0 {
		return meridional(0) // on-axis: the axis itself
	}
	// Merit: the height of the hit on the stop surface (0 when the ray passes
	// through the stop centre).
	merit := func(y0 float64) (float64, bool) {
		res := s.Trace(meridional(y0))
		if !res.OK() || len(res.Steps) <= stop {
			return 0, false
		}
		return res.Steps[stop].Point.Y, true
	}
	bestY, bestM, haveBest := y0p, math.Inf(1), false
	if m, ok := merit(y0p); ok {
		bestM, haveBest = math.Abs(m), true
	}
	// Secant iteration from two nearby starting points (deterministic).
	a, b := y0p, y0p+0.05*hMag+1e-3
	fa, oka := merit(a)
	fb, okb := merit(b)
	for i := 0; i < 20 && oka && okb; i++ {
		if fb == fa {
			break
		}
		c := b - fb*(b-a)/(fb-fa)
		if !finite(c) {
			break
		}
		fc, okc := merit(c)
		if okc && math.Abs(fc) < bestM {
			bestM, bestY, haveBest = math.Abs(fc), c, true
		}
		a, fa, oka = b, fb, okb
		b, fb, okb = c, fc, okc
		if okc && math.Abs(fc) < 1e-12*(1+math.Abs(c)) {
			break
		}
	}
	ray := meridional(bestY)
	if !haveBest {
		ray = meridional(y0p)
	}
	// Rotate from the meridional solution (which lives at azimuth +90°) to the
	// azimuth of the actual field. The system is rotationally symmetric, so
	// the rotated ray is exactly the chief ray of the rotated field.
	rotAngle := azimuth - math.Pi/2
	if math.Abs(rotAngle) > 1e-15 {
		ca, sa := math.Cos(rotAngle), math.Sin(rotAngle)
		rot := func(v Vec3) Vec3 { return V(v.X*ca-v.Y*sa, v.X*sa+v.Y*ca, v.Z) }
		ray = Ray{Origin: rot(ray.Origin), Dir: rot(ray.Dir), WavelengthNm: ray.WavelengthNm, Status: ray.Status}
	}
	return ray
}

// ChiefOffset returns the transverse position (its x and y components in mm)
// of the line through the stop vertex along the field direction (or from the
// object point through the stop vertex) at the plane z = zPlane: the offset
// that turns a bare plane into the classical pupil frame centred on the chief
// ray. Set it as Pupil.Offset.
//
// It is the launch-point correction that aligns a pupil sampling with the
// third-order model: the Seidel sums assume the pupil centre to be the chief
// ray, which crosses the stop centre and therefore crosses a general launch
// plane off-axis whenever the stop is not the entrance pupil itself (a single
// mirror, for instance, has the mirror surface as entrance pupil). When
// refracting surfaces precede the stop the line is the chief ray's path in the
// object space to first order (the exact real chief is System.ChiefRay).
//
// The plane must lie in the object space (in front of the first surface).
func (s *System) ChiefOffset(f Field, zPlane float64) (Vec3, error) {
	if s.hasPose() {
		return Zero, fmt.Errorf("geo: chief offset: the system has coordinate breaks or posed surfaces; " +
			"use ChiefCrossing(field, zPlane, wavelength), which needs the wavelength of the trace")
	}
	stop := s.StopIndex
	if stop < 0 || stop >= len(s.Surfaces) {
		stop = 0
	}
	zc := s.VertexZ(stop)
	if !f.Finite {
		d, ok := f.Direction()
		if !ok || d.Z == 0 {
			return Zero, fmt.Errorf("geo: chief offset: degenerate field direction")
		}
		t := (zPlane - zc) / d.Z
		p := V(0, 0, zc).Add(d.Mul(t))
		return V(p.X, p.Y, 0), nil
	}
	o := f.ObjectPoint()
	dir := V(0, 0, zc).Sub(o)
	if dir.Z == 0 {
		return Zero, fmt.Errorf("geo: chief offset: stop plane is parallel to the chief ray")
	}
	t := (zPlane - o.Z) / dir.Z
	p := o.Add(dir.Mul(t))
	return V(p.X, p.Y, 0), nil
}

// ParaxialImageHeight returns the paraxial image height (mm) of the field at
// the image plane z = imageZ, evaluated in the meridional (y) direction. For
// an object at infinity this is EFL·tan θy; for a finite object the transverse
// magnification times the object height. The system must be refractive (the
// first-order engine has no mirror support).
func (s *System) ParaxialImageHeight(f Field, imageZ, lambdaNm float64) (float64, error) {
	p, err := s.Paraxial(lambdaNm, f.ObjectPlane())
	if err != nil {
		return 0, err
	}
	var t float64
	if f.Finite {
		t = f.HY
	} else {
		t = f.TanY
	}
	c := p.ChiefRay(t)
	k := len(c.Y) - 1
	if k < 0 {
		return 0, fmt.Errorf("geo: empty system")
	}
	return c.Y[k] + (imageZ-c.VertexZ[k])*c.UAfter[k], nil
}
