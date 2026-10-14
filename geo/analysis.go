package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Reference points
// ---------------------------------------------------------------------------

// ReferenceMode selects the transverse reference point about which spot and
// OPD data are evaluated.
type ReferenceMode int

const (
	// RefChief: the image-plane hit of the chief ray of the sampling (the ray
	// through the centre of the pupil). This is the zero value and the
	// classical choice.
	RefChief ReferenceMode = iota
	// RefCentroid: the centroid of the successfully traced image-plane hits.
	RefCentroid
	// RefVertex: the point where the optical axis pierces the image plane,
	// (0, 0, imageZ).
	RefVertex
	// RefExplicit: the point Reference.Point.
	RefExplicit
)

// Reference identifies the transverse reference point used by the analysis
// routines. The zero value is the chief ray reference.
type Reference struct {
	Mode  ReferenceMode
	Point Vec3 // used when Mode == RefExplicit
}

// ChiefReference returns the chief-ray reference (the default).
func ChiefReference() Reference { return Reference{} }

// CentroidReference returns the centroid reference.
func CentroidReference() Reference { return Reference{Mode: RefCentroid} }

// VertexReference returns the on-axis (optical axis) reference.
func VertexReference() Reference { return Reference{Mode: RefVertex} }

// ExplicitReference returns a reference at the given image-plane point.
func ExplicitReference(p Vec3) Reference { return Reference{Mode: RefExplicit, Point: p} }

// String implements fmt.Stringer.
func (r Reference) String() string {
	switch r.Mode {
	case RefCentroid:
		return "centroid"
	case RefVertex:
		return "vertex"
	case RefExplicit:
		return fmt.Sprintf("explicit%s", r.Point)
	}
	return "chief"
}

// ---------------------------------------------------------------------------
// Tracing a ray to an image plane
// ---------------------------------------------------------------------------

// imageHit holds what the analysis routines need from tracing one ray to a
// plane z = imageZ: the image-plane intersection and the optical path length
// to it, plus the exit data at the last surface (used by the OPD computation,
// whose reference sphere is anchored at the exit point of the chief ray).
type imageHit struct {
	Status  Status
	Hit     Vec3    // image-plane intersection (mm, valid when Status is OK)
	OPL     float64 // optical path length from the ray origin to the plane (mm)
	Exit    Vec3    // last surface hit point (mm)
	ExitOPL float64
	Dir     Vec3    // propagation direction after the last surface
	Index   float64 // refractive index after the last surface
}

// traceToPlane traces r through the system and propagates the exit ray to the
// plane z = imageZ, returning the intersection point, the optical path length
// from the launch point to the plane and the exit data at the last surface.
//
// Analysis routines require a system that describes the optical surfaces only;
// an image plane appended with WithImagePlane would be treated as the "last
// surface" and would flatten the OPD reference geometry (the plane z = imageZ
// is implied by the imageZ argument instead).
func (s *System) traceToPlane(r Ray, imageZ float64) imageHit {
	out := imageHit{}
	res := s.Trace(r)
	if len(res.Steps) == 0 {
		out.Status = res.Status
		return out
	}
	last := res.Steps[len(res.Steps)-1]
	out.Exit = last.Point
	out.ExitOPL = last.OPL
	out.Dir = res.Ray.Dir
	out.Index = s.Mediums[len(s.Mediums)-1].N(r.WavelengthNm)
	if !res.OK() {
		out.Status = res.Status
		return out
	}
	d := res.Ray.Dir
	if d.Z == 0 {
		out.Status = StatusEscape
		return out
	}
	t := (imageZ - res.Ray.Origin.Z) / d.Z
	if !finite(t) || t < 0 {
		// The ray does not reach the image plane travelling forwards (a
		// backwards ray, or a plane behind the exit direction).
		out.Status = StatusEscape
		return out
	}
	out.Status = StatusOK
	out.Hit = res.Ray.Origin.Add(d.Mul(t))
	out.OPL = res.Ray.OPL + out.Index*t
	return out
}

// ---------------------------------------------------------------------------
// Spot diagrams
// ---------------------------------------------------------------------------

// SpotPoint is the outcome of tracing one sampled ray to the image plane.
type SpotPoint struct {
	Pupil  PupilPoint // normalized pupil coordinates of the ray
	Hit    Vec3       // image-plane intersection (valid when Status is OK)
	OPL    float64    // optical path length from the launch point to the image plane
	Status Status     // StatusOK when the ray reached the image plane
}

// SpotDiagram is the image-plane spot of one field at one wavelength: one
// entry per requested pupil point, in the order of the request.
type SpotDiagram struct {
	Field  Field
	Lambda float64
	Pupil  Pupil
	ImageZ float64
	Points []SpotPoint // one per requested pupil point
	// Chief is the image-plane hit of the chief ray (the ray through the pupil
	// centre), traced in addition to the requested points.
	Chief       Vec3
	ChiefStatus Status
	// Centroid is the centroid of the successful image-plane hits.
	Centroid Vec3
}

// SpotDiagram traces the sampling (field, pupil, points) to the image plane
// z = imageZ and returns the spot. Rays that fail (vignetting, total internal
// reflection, escape) are reported with their status and excluded from the
// statistics; everything is deterministic, so equal inputs give equal outputs.
func (s *System) SpotDiagram(f Field, pup Pupil, pts []PupilPoint, imageZ, lambdaNm float64) (*SpotDiagram, error) {
	if err := s.Validate(); err != nil {
		return nil, err
	}
	if !finite(imageZ) {
		return nil, fmt.Errorf("geo: spot diagram: image plane z is not finite")
	}
	if !finitePositive(pup.Radius) {
		return nil, fmt.Errorf("geo: spot diagram: pupil radius %g is not positive", pup.Radius)
	}
	sd := &SpotDiagram{Field: f, Lambda: lambdaNm, Pupil: pup, ImageZ: imageZ}
	if len(pts) == 0 {
		return nil, fmt.Errorf("geo: spot diagram: no pupil points")
	}
	sd.Points = make([]SpotPoint, len(pts))
	sum := Zero
	nOK := 0
	for i, p := range pts {
		h := s.traceToPlane(LaunchRay(f, pup, p, lambdaNm), imageZ)
		sd.Points[i] = SpotPoint{Pupil: p, Hit: h.Hit, OPL: h.OPL, Status: h.Status}
		if h.Status.OK() {
			sum = sum.Add(h.Hit)
			nOK++
		}
	}
	if nOK > 0 {
		sd.Centroid = sum.Div(float64(nOK))
	}
	h := s.traceToPlane(LaunchRay(f, pup, PupilPoint{}, lambdaNm), imageZ)
	sd.Chief, sd.ChiefStatus = h.Hit, h.Status
	return sd, nil
}

// SpotStats summarizes a spot diagram about a reference point: the RMS radius
// (root mean square distance of the hits from the reference), the GEO radius
// (largest distance from the reference) and the centroid of the hits. All
// values are in mm and are computed over the successfully traced rays only.
type SpotStats struct {
	Reference Vec3 // the resolved reference point
	Centroid  Vec3 // centroid of the valid hits
	N         int  // number of valid hits
	RMS       float64
	GEO       float64 // largest distance from Reference
}

// Stats computes the spot statistics about the given reference point. RefChief
// uses the chief-ray hit (an error of the chief trace makes the reference
// undefined and Stats falls back to the centroid); RefCentroid uses the
// centroid; RefVertex the optical axis of the image plane; RefExplicit the
// given point.
func (sp *SpotDiagram) Stats(ref Reference) SpotStats {
	st := SpotStats{Centroid: sp.Centroid}
	origin := V(0, 0, sp.ImageZ)
	switch ref.Mode {
	case RefCentroid:
		var sum Vec3
		n := 0
		for _, p := range sp.Points {
			if p.Status.OK() {
				sum = sum.Add(p.Hit)
				n++
			}
		}
		if n == 0 {
			return st
		}
		st.Centroid = sum.Div(float64(n))
		origin = st.Centroid
	case RefVertex:
		origin = V(0, 0, sp.ImageZ)
	case RefExplicit:
		origin = ref.Point
	default: // RefChief
		if sp.ChiefStatus.OK() {
			origin = sp.Chief
		} else {
			var sum Vec3
			n := 0
			for _, p := range sp.Points {
				if p.Status.OK() {
					sum = sum.Add(p.Hit)
					n++
				}
			}
			if n == 0 {
				return st
			}
			st.Centroid = sum.Div(float64(n))
			origin = st.Centroid
		}
	}
	st.Reference = origin
	sum2 := 0.0
	max2 := 0.0
	n := 0
	for _, p := range sp.Points {
		if !p.Status.OK() {
			continue
		}
		d2 := p.Hit.Sub(origin).Norm2()
		sum2 += d2
		if d2 > max2 {
			max2 = d2
		}
		n++
	}
	st.N = n
	if n > 0 {
		st.RMS = math.Sqrt(sum2 / float64(n))
		st.GEO = math.Sqrt(max2)
	}
	return st
}

// StatsAt is Stats with an explicit reference point.
func (sp *SpotDiagram) StatsAt(p Vec3) SpotStats { return sp.Stats(ExplicitReference(p)) }

// ---------------------------------------------------------------------------
// Transverse aberration (ray fans)
// ---------------------------------------------------------------------------

// FanPoint is one sample of a transverse-aberration curve: the pupil
// coordinate P along the fan (−1 to +1), the image-plane hit and the two
// transverse aberration components, i.e. the hit minus the reference point
// (mm). DX is the aberration along x (sagittal), DY along y (meridional).
type FanPoint struct {
	P      float64    // fan parameter, −1 … +1
	Pupil  PupilPoint // normalized pupil coordinates
	Hit    Vec3
	DX, DY float64
	OPL    float64
	Status Status
}

// RayFan is the transverse-aberration curve of a fan of rays: one entry per
// sampled pupil point.
type RayFan struct {
	Field      Field
	Lambda     float64
	Pupil      Pupil
	ImageZ     float64
	AzimuthDeg float64 // fan direction in the pupil
	Reference  Reference
	RefPoint   Vec3 // resolved reference point (mm)
	Points     []FanPoint
}

// RayFan traces n rays across the pupil along the direction azimuthDeg
// (degrees; 0 = sagittal/x, 90 = meridional/y) and returns the transverse
// aberration curve relative to the reference point. Use MeridionalRayFan and
// SagittalRayFan for the two classical fans.
func (s *System) RayFan(f Field, pup Pupil, n int, azimuthDeg, imageZ, lambdaNm float64, ref Reference) (*RayFan, error) {
	pts := Fan(n, azimuthDeg)
	if len(pts) == 0 {
		return nil, fmt.Errorf("geo: ray fan: n must be at least 1")
	}
	rf := &RayFan{Field: f, Lambda: lambdaNm, Pupil: pup, ImageZ: imageZ, AzimuthDeg: azimuthDeg, Reference: ref}
	// Resolve the reference point (chief hit, centroid of the fan, axis or
	// explicit) from the traced rays.
	rf.Points = make([]FanPoint, len(pts))
	sum := Zero
	nOK := 0
	for i, p := range pts {
		h := s.traceToPlane(LaunchRay(f, pup, p, lambdaNm), imageZ)
		// Fan parameter: the projection of the pupil point onto the fan
		// direction, t = p·(cosφ, sinφ). The cardinal fans (φ = 0/90°) make
		// one term vanish, so this is exact for them, and it is also exact
		// for oblique azimuths where a single-component projection would be
		// only approximate.
		t := 0.0
		if n > 1 {
			ca := math.Cos(azimuthDeg * math.Pi / 180)
			sa := math.Sin(azimuthDeg * math.Pi / 180)
			t = p.PX*ca + p.PY*sa
		}
		rf.Points[i] = FanPoint{P: t, Pupil: p, Hit: h.Hit, OPL: h.OPL, Status: h.Status}
		if h.Status.OK() {
			sum = sum.Add(h.Hit)
			nOK++
		}
	}
	origin := V(0, 0, imageZ)
	switch ref.Mode {
	case RefCentroid:
		if nOK == 0 {
			return nil, fmt.Errorf("geo: ray fan: no valid rays to compute the centroid")
		}
		origin = sum.Div(float64(nOK))
	case RefExplicit:
		origin = ref.Point
	default: // RefChief
		h := s.traceToPlane(LaunchRay(f, pup, PupilPoint{}, lambdaNm), imageZ)
		if !h.Status.OK() {
			return nil, fmt.Errorf("geo: ray fan: the chief ray (%s) did not reach the image plane", h.Status)
		}
		origin = h.Hit
	}
	rf.RefPoint = origin
	for i := range rf.Points {
		p := &rf.Points[i]
		if p.Status.OK() {
			p.DX = p.Hit.X - origin.X
			p.DY = p.Hit.Y - origin.Y
		}
	}
	return rf, nil
}

// MeridionalRayFan is RayFan with the fan along the meridional (y) diameter
// (azimuth 90°): the classical tangential transverse-aberration curve, whose
// DY carries the aberration of a field in the y direction.
func (s *System) MeridionalRayFan(f Field, pup Pupil, n int, imageZ, lambdaNm float64, ref Reference) (*RayFan, error) {
	return s.RayFan(f, pup, n, 90, imageZ, lambdaNm, ref)
}

// SagittalRayFan is RayFan with the fan along the sagittal (x) diameter
// (azimuth 0°).
func (s *System) SagittalRayFan(f Field, pup Pupil, n int, imageZ, lambdaNm float64, ref Reference) (*RayFan, error) {
	return s.RayFan(f, pup, n, 0, imageZ, lambdaNm, ref)
}

// RMS returns the root mean square of sqrt(DX²+DY²) over the valid fan points
// (mm), a single-number measure of the fan's aberration.
func (rf *RayFan) RMS() float64 {
	sum, n := 0.0, 0
	for _, p := range rf.Points {
		if !p.Status.OK() {
			continue
		}
		sum += p.DX*p.DX + p.DY*p.DY
		n++
	}
	if n == 0 {
		return math.NaN()
	}
	return math.Sqrt(sum / float64(n))
}

// ---------------------------------------------------------------------------
// OPD: optical path difference over the pupil
// ---------------------------------------------------------------------------

// OPDOptions controls the OPD computation.
//
// Reference selects the image point C the reference sphere is centred at
// (RefChief, the default, uses the image-plane hit of the chief ray;
// RefCentroid the centroid of the traced hits; RefVertex the optical axis;
// RefExplicit the given point).
//
// RemovePiston subtracts the mean OPD over the valid samples and RemoveTilt
// the least-squares best-fit linear term c1·px + c2·py (the sampled pupil
// coordinates). The raw OPD is already referenced to the chief ray, so its
// value at the chief sample is zero; the two options are useful when an
// off-centre reference point or a curved reference plane leaves a piston or
// tilt in the map.
type OPDOptions struct {
	Reference    Reference
	RemovePiston bool
	RemoveTilt   bool
}

// OPDGrid is the optical path difference (wavefront aberration) of one field
// over a set of pupil points, in mm (divide by the wavelength in mm for
// waves; see RMSWaves).
//
// Definition. Let C be the reference image point (see OPDOptions) and Q the
// point where the chief ray leaves the last surface of the system; the
// reference sphere is the sphere centred at C through Q, of radius
// R = |Q − C|. For every ray, O_i is the optical path length from the launch
// point to its own exit point Q_i at the last surface and t_i the distance
// along the exit ray from Q_i to its first crossing of the reference sphere;
// then
//
//	OPD_i = O_i + n'·t_i − O_chief ,
//
// with n' the index after the last surface. This is the OPL of the ray to the
// reference sphere relative to the chief ray (the wavefront error measured on
// the reference sphere); for a perfect system all rays cross the sphere with
// equal optical path and the OPD is zero. The sphere is anchored on the exit
// surface of the system rather than on the paraxial exit pupil, so a pure
// refocusing of the reference plane appears as a ρ² (defocus) term, which the
// piston/tilt options and a focus scan around the best-focus plane remove in
// practice.
type OPDGrid struct {
	Field   Field
	Lambda  float64
	Pupil   Pupil
	ImageZ  float64
	Options OPDOptions
	// ReferencePoint is the resolved reference image point C (mm) and
	// SphereRadius the radius of the reference sphere |Q − C|.
	ReferencePoint Vec3
	SphereRadius   float64
	Pupils         []PupilPoint
	Raw            []float64 // OPD before piston/tilt removal (mm, NaN for failed rays)
	OPD            []float64 // OPD after piston/tilt removal (mm)
	RMS            float64   // RMS of OPD over the valid samples (mm)
	PV             float64   // peak-to-valley of OPD (mm)
	N              int       // number of valid samples
}

// OPDGrid traces the sampling to the image plane z = imageZ and computes the
// optical path difference over the pupil (see the OPDGrid type for the exact
// definition). It needs at least one ray that reaches the plane and a chief
// ray that reaches it too.
func (s *System) OPDGrid(f Field, pup Pupil, pts []PupilPoint, imageZ, lambdaNm float64, opts OPDOptions) (*OPDGrid, error) {
	if err := s.Validate(); err != nil {
		return nil, err
	}
	if !finite(imageZ) {
		return nil, fmt.Errorf("geo: OPD: image plane z is not finite")
	}
	if !finitePositive(pup.Radius) {
		return nil, fmt.Errorf("geo: OPD: pupil radius %g is not positive", pup.Radius)
	}
	if len(pts) == 0 {
		return nil, fmt.Errorf("geo: OPD: no pupil points")
	}
	type rayData struct {
		hit  imageHit
		ok   bool
		tSph float64
		opl  float64 // OPL to the reference sphere
	}
	data := make([]rayData, len(pts))
	var sumHit Vec3
	nHit := 0
	for i, p := range pts {
		h := s.traceToPlane(LaunchRay(f, pup, p, lambdaNm), imageZ)
		data[i].hit = h
		data[i].ok = h.Status.OK()
		if data[i].ok {
			sumHit = sumHit.Add(h.Hit)
			nHit++
		}
	}
	chief := s.traceToPlane(LaunchRay(f, pup, PupilPoint{}, lambdaNm), imageZ)
	if !chief.Status.OK() {
		return nil, fmt.Errorf("geo: OPD: the chief ray (%s) did not reach the image plane", chief.Status)
	}
	var c Vec3
	switch opts.Reference.Mode {
	case RefCentroid:
		if nHit == 0 {
			return nil, fmt.Errorf("geo: OPD: no valid rays to compute the centroid")
		}
		c = sumHit.Div(float64(nHit))
	case RefVertex:
		c = V(0, 0, imageZ)
	case RefExplicit:
		c = opts.Reference.Point
	default:
		c = chief.Hit
	}
	rS := chief.Exit.Distance(c)
	if !finitePositive(rS) {
		return nil, fmt.Errorf("geo: OPD: degenerate reference sphere (radius %g)", rS)
	}
	g := &OPDGrid{
		Field: f, Lambda: lambdaNm, Pupil: pup, ImageZ: imageZ, Options: opts,
		ReferencePoint: c, SphereRadius: rS,
		Pupils: append([]PupilPoint(nil), pts...),
		Raw:    make([]float64, len(pts)),
		OPD:    make([]float64, len(pts)),
	}
	// The OPL of an infinite-field ray is accumulated from its own launch
	// point, and for a tilted beam those points are not all on one wavefront
	// of the incident plane wave. Projecting the launch points onto the plane
	// through the pupil centre perpendicular to the beam direction — i.e.
	// adding (P_i − P_chief)·d — puts every ray on the same incident
	// wavefront before the OPLs are compared. (For a finite object all rays
	// start at the same object point and the correction vanishes.)
	incident := func(p PupilPoint) float64 {
		if f.Finite {
			return 0
		}
		d, ok := f.Direction()
		if !ok {
			return 0
		}
		return V(p.PX*pup.Radius, p.PY*pup.Radius, 0).Dot(d)
	}
	oplChief := chief.ExitOPL
	for i := range pts {
		g.Raw[i] = math.NaN()
		if !data[i].ok {
			continue
		}
		h := data[i].hit
		t, ok := sphereCrossing(h.Exit, h.Dir, c, rS)
		if !ok {
			continue
		}
		opl := h.ExitOPL + h.Index*t + incident(pts[i])
		g.Raw[i] = opl - oplChief
	}
	// Piston/tilt removal, least-squares over the valid samples.
	copy(g.OPD, g.Raw)
	if opts.RemovePiston || opts.RemoveTilt {
		if c0, c1, c2, ok := fitPistonTilt(g.Pupils, g.Raw); ok {
			for i := range g.OPD {
				if math.IsNaN(g.OPD[i]) {
					continue
				}
				if opts.RemovePiston {
					g.OPD[i] -= c0
				}
				if opts.RemoveTilt {
					g.OPD[i] -= c1*g.Pupils[i].PX + c2*g.Pupils[i].PY
				}
			}
		}
	}
	// Statistics over the valid samples.
	sum, sum2, vmin, vmax, n := 0.0, 0.0, math.Inf(1), math.Inf(-1), 0
	for _, v := range g.OPD {
		if math.IsNaN(v) {
			continue
		}
		sum += v
		sum2 += v * v
		if v < vmin {
			vmin = v
		}
		if v > vmax {
			vmax = v
		}
		n++
	}
	g.N = n
	if n > 0 {
		mean := sum / float64(n)
		var2 := sum2/float64(n) - mean*mean
		if var2 < 0 {
			var2 = 0
		}
		g.RMS = math.Sqrt(var2)
		g.PV = vmax - vmin
	}
	return g, nil
}

// sphereCrossing returns the parameter t at which the line q + t·d (d unit)
// crosses the sphere of radius r centred at c, choosing the crossing nearest
// to the ray's starting point q (t = 0). t may be slightly negative when q
// lies just inside the sphere; the crossing is then the backward extension of
// the ray, which keeps the OPD continuous as the exit surface of a system
// passes through the reference sphere. ok is false when the line does not meet
// the sphere (it passes farther than r from the centre).
func sphereCrossing(q, d, c Vec3, r float64) (float64, bool) {
	v := q.Sub(c)
	u := d.Dot(v)
	w2 := v.Norm2()
	disc := u*u - w2 + r*r
	if disc < 0 || math.IsNaN(disc) {
		return 0, false
	}
	sq := math.Sqrt(disc)
	t1, t2 := -u-sq, -u+sq
	if math.Abs(t1) <= math.Abs(t2) {
		return t1, true
	}
	return t2, true
}

// fitPistonTilt fits c0 + c1·px + c2·py to the valid samples by least squares
// and returns the coefficients. ok is false when fewer than three valid
// samples are available (the tilt part is then under-determined).
func fitPistonTilt(pts []PupilPoint, v []float64) (c0, c1, c2 float64, ok bool) {
	// Normal equations of the 3-parameter linear fit; the basis is small and
	// well-conditioned for any sampling symmetric about the pupil centre.
	var s [3][3]float64
	var b [3]float64
	basis := func(p PupilPoint) [3]float64 { return [3]float64{1, p.PX, p.PY} }
	n := 0
	for i, p := range pts {
		if i >= len(v) || math.IsNaN(v[i]) {
			continue
		}
		bf := basis(p)
		for r := 0; r < 3; r++ {
			for c := 0; c < 3; c++ {
				s[r][c] += bf[r] * bf[c]
			}
			b[r] += bf[r] * v[i]
		}
		n++
	}
	if n < 3 {
		return 0, 0, 0, false
	}
	x, ok := solve3(s, b)
	if !ok {
		return 0, 0, 0, false
	}
	return x[0], x[1], x[2], true
}

// solve3 solves the 3×3 linear system a·x = b by Gauss–Jordan elimination with
// partial pivoting. ok is false for a singular system.
func solve3(a [3][3]float64, b [3]float64) ([3]float64, bool) {
	m := a
	v := b
	for col := 0; col < 3; col++ {
		p := col
		for r := col + 1; r < 3; r++ {
			if math.Abs(m[r][col]) > math.Abs(m[p][col]) {
				p = r
			}
		}
		if math.Abs(m[p][col]) < 1e-300 {
			return [3]float64{}, false
		}
		m[col], m[p] = m[p], m[col]
		v[col], v[p] = v[p], v[col]
		inv := 1 / m[col][col]
		for c := col; c < 3; c++ {
			m[col][c] *= inv
		}
		v[col] *= inv
		for r := 0; r < 3; r++ {
			if r == col || m[r][col] == 0 {
				continue
			}
			f := m[r][col]
			for c := col; c < 3; c++ {
				m[r][c] -= f * m[col][c]
			}
			v[r] -= f * v[col]
		}
	}
	return v, true
}

// Waves converts an OPD value in mm to waves at the grid's wavelength.
func (g *OPDGrid) Waves(mm float64) float64 { return mm / (g.Lambda * 1e-6) }

// RMSWaves returns the RMS of the (removed) OPD in waves.
func (g *OPDGrid) RMSWaves() float64 { return g.Waves(g.RMS) }

// PVWaves returns the peak-to-valley OPD in waves.
func (g *OPDGrid) PVWaves() float64 { return g.Waves(g.PV) }

// MaxWaves returns the largest |OPD| in waves.
func (g *OPDGrid) MaxWaves() float64 {
	m := 0.0
	for _, v := range g.OPD {
		if math.IsNaN(v) {
			continue
		}
		if a := math.Abs(v); a > m {
			m = a
		}
	}
	return g.Waves(m)
}

// ---------------------------------------------------------------------------
// Distortion (real rays)
// ---------------------------------------------------------------------------

// Distortion returns the relative distortion of the field at the image plane
// z = imageZ: the difference between the image height of the real chief ray
// (the ray through the exact centre of the aperture stop, see System.ChiefRay)
// and the paraxial image height, divided by the paraxial image height (so the
// result is a fraction: 0.01 = 1% pincushion for positive values). The
// meridional component is returned (the transverse y direction).
func (s *System) Distortion(f Field, imageZ, lambdaNm float64) (float64, error) {
	if err := s.Validate(); err != nil {
		return 0, err
	}
	chief := s.ChiefRay(f, lambdaNm)
	h := s.traceToPlane(chief, imageZ)
	if !h.Status.OK() {
		return 0, fmt.Errorf("geo: distortion: chief ray (%s) did not reach the image plane", h.Status)
	}
	par, err := s.ParaxialImageHeight(f, imageZ, lambdaNm)
	if err != nil {
		return 0, err
	}
	if !finitePositive(math.Abs(par)) {
		return 0, fmt.Errorf("geo: distortion: paraxial image height is zero (on-axis field)")
	}
	return (h.Hit.Y - par) / par, nil
}
