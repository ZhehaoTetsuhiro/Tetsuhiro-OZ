package geo

import (
	"fmt"
	"math"
)

// tMinDefault is the smallest ray parameter accepted at a surface (mm). It is
// slightly negative so that a ray which starts exactly on a surface (reversed
// traces, double passes, zero-thickness systems) is not rejected by the
// round-off of the previous intersection, which is of the order of 1e-16 mm.
// The path length contributed by such a segment is clamped to zero.
const tMinDefault = -1e-9

// System is a sequential optical system: an ordered list of surfaces separated
// by axial thicknesses and filled with media.
//
// Layout: the vertex of surface 0 sits at z = 0; the vertex of surface i sits at
// the sum of Thickness[0..i-1]. Thickness has one entry fewer than Surfaces
// (Thickness[i] is the axial distance from surface i to surface i+1, mm) and
// may be negative when a fold or an unusual layout requires it. Mediums has one
// entry more than Surfaces: Mediums[i] is the medium between surface i-1 and
// surface i, so Mediums[0] is the object-side medium and Mediums[len(Surfaces)]
// the image-side medium. Every medium is evaluated at the ray wavelength
// through the Medium interface, so dispersive systems are described simply by
// placing Sellmeier glasses in the gaps.
type System struct {
	Surfaces  []Surface
	Thickness []float64
	Mediums   []Medium
	StopIndex int // index of the aperture stop surface
	// OriginZ is the axial position (mm) of the vertex of surface 0 in the
	// global frame. It is zero for a freshly built system and non-zero for
	// sub-systems and time-reversed systems so that they occupy the same
	// physical space as their parent.
	OriginZ float64
	Name    string // optional label
}

// NewSystem builds a system from its parts and validates it. Pass nil for
// mediums to fill every gap with air, and nil for thicknesses only when the
// system has a single surface.
func NewSystem(medium0 Medium, surfaces []Surface, thicknesses []float64, mediums []Medium) (*System, error) {
	s := &System{Surfaces: surfaces, Thickness: thicknesses, Mediums: mediums}
	if mediums == nil {
		s.Mediums = make([]Medium, len(surfaces)+1)
		for i := range s.Mediums {
			s.Mediums[i] = Air{}
		}
	}
	if medium0 != nil && len(s.Mediums) > 0 {
		s.Mediums[0] = medium0
	}
	if err := s.Validate(); err != nil {
		return nil, err
	}
	return s, nil
}

// SimpleSystem builds a system whose gaps are all filled with medium m (air,
// typically — use NewSystem for cemented or immersed systems).
func SimpleSystem(m Medium, surfaces []Surface, thicknesses []float64) (*System, error) {
	meds := make([]Medium, len(surfaces)+1)
	if m == nil {
		m = Air{}
	}
	for i := range meds {
		meds[i] = m
	}
	return NewSystem(m, surfaces, thicknesses, meds)
}

// Validate checks the internal consistency of the system.
func (s *System) Validate() error {
	n := len(s.Surfaces)
	if n == 0 {
		return fmt.Errorf("geo: system has no surfaces")
	}
	for i := range s.Surfaces {
		sf := &s.Surfaces[i]
		if sf.Break {
			if sf.Mirror {
				return fmt.Errorf("geo: element %d is a coordinate break and cannot be a mirror", i)
			}
			if sf.Pose != nil && !sf.Pose.finite() {
				return fmt.Errorf("geo: coordinate break %d has non-finite parameters (%+v)", i, *sf.Pose)
			}
			continue
		}
		if sf.Profile == nil {
			return fmt.Errorf("geo: surface %d has a nil profile", i)
		}
		if sf.Pose != nil && !sf.Pose.finite() {
			return fmt.Errorf("geo: surface %d has a non-finite pose (%+v)", i, *sf.Pose)
		}
	}
	if len(s.Mediums) != n+1 {
		return fmt.Errorf("geo: got %d medium entries for %d surfaces, want %d",
			len(s.Mediums), n, n+1)
	}
	if len(s.Thickness) != n-1 {
		return fmt.Errorf("geo: got %d thicknesses for %d surfaces, want %d",
			len(s.Thickness), n, n-1)
	}
	for i, m := range s.Mediums {
		if m == nil {
			return fmt.Errorf("geo: medium %d is nil", i)
		}
	}
	for i, t := range s.Thickness {
		if math.IsNaN(t) || math.IsInf(t, 0) {
			return fmt.Errorf("geo: thickness %d is not finite (%v)", i, t)
		}
	}
	for i := range s.Surfaces {
		if s.Surfaces[i].Break && s.Mediums[i].Name() != s.Mediums[i+1].Name() {
			return fmt.Errorf("geo: coordinate break %d separates different media (%s, %s)",
				i, s.Mediums[i].Name(), s.Mediums[i+1].Name())
		}
	}
	if s.StopIndex < 0 || s.StopIndex >= n {
		return fmt.Errorf("geo: stop index %d out of range [0,%d]", s.StopIndex, n-1)
	}
	if s.Surfaces[s.StopIndex].Break {
		return fmt.Errorf("geo: stop index %d is a coordinate break, not a surface", s.StopIndex)
	}
	return nil
}

// NumSurfaces returns the number of surfaces.
func (s *System) NumSurfaces() int { return len(s.Surfaces) }

// VertexZ returns the axial position (mm) of the vertex of surface i in the
// global frame: for a system without poses this is OriginZ plus the thicknesses
// in front of it; for a system with coordinate breaks or posed surfaces it is
// the z coordinate of the vertex on the chain of frames (see CoordBreak).
func (s *System) VertexZ(i int) float64 {
	c := s.chainStart()
	n := len(s.Surfaces)
	for k := 0; k < i && k < n && k < len(s.Thickness); k++ {
		_, ref := s.place(k, c)
		c = s.advance(ref, k)
	}
	if i >= 0 && i < n {
		e, _ := s.place(i, c)
		return e.o.Z
	}
	return c.o.Z
}

// IndexBefore returns the refractive index in front of surface i at wavelength
// wavelengthNm.
func (s *System) IndexBefore(i int, wavelengthNm float64) float64 {
	return s.Mediums[i].N(wavelengthNm)
}

// IndexAfter returns the refractive index behind surface i.
func (s *System) IndexAfter(i int, wavelengthNm float64) float64 {
	return s.Mediums[i+1].N(wavelengthNm)
}

// BoundsZ returns the smallest and largest z (mm) reached by any surface over
// the clear aperture (or over radius rMax for surfaces with unlimited
// apertures). Rays may be launched safely in front of zmin.
func (s *System) BoundsZ(rMax float64) (float64, float64) {
	zmin, zmax := math.Inf(1), math.Inf(-1)
	c := s.chainStart()
	for i, sf := range s.Surfaces {
		e, ref := s.place(i, c)
		c = s.advance(ref, i)
		if sf.Break {
			continue
		}
		var lo, hi float64
		if e.plain {
			z := e.o.Z
			lo, hi = sf.ZRange(rMax)
			lo, hi = z+lo, z+hi
		} else {
			lo, hi = frameZRange(sf, e, rMax)
		}
		if lo < zmin {
			zmin = lo
		}
		if hi > zmax {
			zmax = hi
		}
	}
	return zmin, zmax
}

// frameZRange returns the smallest and largest global z reached by the clear
// aperture of a surface placed in the frame e (the aperture radius is
// MaxRadius, or rMax for an unlimited aperture). The surface is sampled on a
// polar grid because the axial sag range does not apply to a tilted frame.
func frameZRange(sf Surface, e frame, rMax float64) (float64, float64) {
	rr := sf.MaxRadius()
	if math.IsInf(rr, 1) {
		rr = rMax
	}
	if !(rr > 0) || math.IsInf(rr, 1) {
		return math.Inf(-1), math.Inf(1)
	}
	const nr, na = 48, 96
	zmin, zmax := math.Inf(1), math.Inf(-1)
	for k := 0; k <= nr; k++ {
		r := rr * float64(k) / nr
		for j := 0; j < na; j++ {
			a := 2 * math.Pi * float64(j) / na
			x, y := r*math.Cos(a), r*math.Sin(a)
			zl := sf.Sag(x, y)
			if math.IsNaN(zl) || math.IsInf(zl, 0) {
				continue
			}
			z := e.toGlobal(Vec3{x, y, zl}).Z
			if z < zmin {
				zmin = z
			}
			if z > zmax {
				zmax = z
			}
		}
	}
	if zmin > zmax {
		return math.Inf(-1), math.Inf(1)
	}
	return zmin, zmax
}

// Step records one ray–surface interaction of a trace.
type Step struct {
	Surface int    // surface index
	Label   string // surface label
	Mirror  bool   // reflective interaction
	Point   Vec3   // intersection point, global frame (mm)
	Normal  Vec3   // unit normal at the hit point, +z oriented (+z may be flipped
	// for the refraction formula; the stored vector is the surface normal)
	DirIn    Vec3    // unit direction before the interaction
	DirOut   Vec3    // unit direction after (equal to DirIn when the ray stopped)
	T        float64 // geometrical path length of this segment (mm)
	Path     float64 // cumulative geometrical path length (mm)
	OPL      float64 // cumulative optical path length (mm)
	NBefore  float64 // refractive index in front of the surface
	NAfter   float64 // refractive index behind the surface
	CosInc   float64 // |cos| of the incidence angle
	IncAngle float64 // incidence angle (rad) in [0, π/2]
	Status   Status  // outcome of this interaction
}

// TraceResult is the outcome of tracing one ray through a System.
type TraceResult struct {
	Ray      Ray    // final ray state (origin = last interaction point)
	Steps    []Step // one entry per surface reached
	Status   Status // final status: StatusOK when the ray traversed every surface
	FailedAt int    // index of the surface that stopped the trace, -1 when OK
}

// OK reports whether the ray traversed the whole system.
func (t TraceResult) OK() bool { return t.Status == StatusOK }

// Final returns the last point reached by the ray (mm, global frame).
func (t TraceResult) Final() Vec3 { return t.Ray.Origin }

// Direction returns the final propagation direction.
func (t TraceResult) Direction() Vec3 { return t.Ray.Dir }

// OPL returns the accumulated optical path length (mm).
func (t TraceResult) OPL() float64 { return t.Ray.OPL }

// Error returns a descriptive error for a failed trace, nil when the ray
// traversed the system.
func (t TraceResult) Error() error {
	if t.Status.OK() {
		return nil
	}
	where := "before the first surface"
	if t.FailedAt >= 0 {
		where = fmt.Sprintf("at surface %d", t.FailedAt)
	}
	return fmt.Errorf("geo: ray %s %s", t.Status, where)
}

// String implements fmt.Stringer.
func (t TraceResult) String() string {
	return fmt.Sprintf("trace[status=%s steps=%d final=%s opl=%.9gmm]",
		t.Status, len(t.Steps), t.Ray.Origin, t.Ray.OPL)
}

// Trace propagates ray r through every surface of the system, applying
// refraction (Snell's law in vector form), reflection at mirrors and
// recording the accumulated optical path length. It stops at the first
// non-OK interaction; the returned TraceResult always carries the partial path
// travelled so far (its Steps and final Ray origin).
//
// The ray's Status field is overwritten; use a fresh ray or NewRay to start a
// new trace.
func (s *System) Trace(r Ray) TraceResult {
	res := TraceResult{FailedAt: -1}
	if !r.Valid() {
		r.Status = StatusDegenerate
		res.Ray, res.Status = r, StatusDegenerate
		return res
	}
	if err := s.Validate(); err != nil {
		r.Status = StatusDegenerate
		res.Ray, res.Status = r, StatusDegenerate
		return res
	}
	r.Status = StatusOK
	lambda := r.WavelengthNm
	c := s.chainStart()
	dir := r.Dir
	for i := range s.Surfaces {
		surf := s.Surfaces[i]
		e, ref := s.place(i, c)
		c = s.advance(ref, i)
		if surf.Break {
			// A pure coordinate break has no interaction: it only changed the chain.
			continue
		}
		nB := s.Mediums[i].N(lambda)
		nA := s.Mediums[i+1].N(lambda)
		if !finitePositive(nB) || !finitePositive(nA) {
			r.Status = StatusDegenerate
			res.Ray, res.Status, res.FailedAt = r, StatusDegenerate, i
			return res
		}
		local := e.toLocal(r.Origin)
		point, t, hit := surf.Intersect(local, e.dirToLocal(dir), tMinDefault)
		step := Step{
			Surface: i, Label: surf.Label(), Mirror: surf.Mirror,
			DirIn: dir, DirOut: dir, NBefore: nB, NAfter: nA,
		}
		if !hit {
			step.Status = StatusEscape
			res.Steps = append(res.Steps, step)
			r.Status = StatusEscape
			res.Ray, res.Status, res.FailedAt = r, StatusEscape, i
			return res
		}
		if t < 0 {
			t = 0
		}
		r.OPL += nB * t
		r.Path += t
		r.Origin = e.toGlobal(point)
		nrm := e.dirToGlobal(surf.Normal(point.X, point.Y))
		if nrm.Norm2() == 0 {
			step.Status = StatusDegenerate
			res.Steps = append(res.Steps, step)
			r.Status = StatusDegenerate
			res.Ray, res.Status, res.FailedAt = r, StatusDegenerate, i
			return res
		}
		cos := math.Abs(dir.Dot(nrm))
		if cos > 1 {
			cos = 1
		}
		step.Point = r.Origin
		step.Normal = nrm
		step.T, step.Path, step.OPL = t, r.Path, r.OPL
		step.CosInc = cos
		step.IncAngle = math.Acos(cos)
		if surf.Vignetted(point.X, point.Y) {
			step.Status = StatusVignetted
			res.Steps = append(res.Steps, step)
			r.Status = StatusVignetted
			res.Ray, res.Status, res.FailedAt = r, StatusVignetted, i
			return res
		}
		newDir := dir
		if surf.Mirror {
			newDir = Reflect(dir, nrm)
		} else if nB != nA {
			nd, ok := Refract(dir, nrm, nB, nA)
			if !ok {
				step.Status = StatusTIR
				res.Steps = append(res.Steps, step)
				r.Status = StatusTIR
				res.Ray, res.Status, res.FailedAt = r, StatusTIR, i
				return res
			}
			newDir = nd
		}
		step.Status = StatusOK
		step.DirOut = newDir
		res.Steps = append(res.Steps, step)
		dir = newDir
		r.Dir = newDir
	}
	r.Status = StatusOK
	res.Ray, res.Status = r, StatusOK
	return res
}

// Sub returns the sub-system made of surfaces [from, to) (from included, to
// excluded). The object-side medium of the sub-system is the medium in front of
// surface 'from'.
func (s *System) Sub(from, to int) (*System, error) {
	if s.hasPose() {
		return nil, fmt.Errorf("geo: sub-systems of systems with coordinate breaks or posed surfaces are not supported")
	}
	if from < 0 || to > len(s.Surfaces) || from >= to {
		return nil, fmt.Errorf("geo: invalid sub-system range [%d,%d)", from, to)
	}
	sub := &System{
		Surfaces:  append([]Surface(nil), s.Surfaces[from:to]...),
		Thickness: append([]float64(nil), s.Thickness[from:min(to, len(s.Thickness))]...),
		Mediums:   append([]Medium(nil), s.Mediums[from:to+1]...),
		OriginZ:   s.VertexZ(from),
		Name:      s.Name,
	}
	// Thickness has to describe the gaps between the selected surfaces only.
	if n := len(sub.Surfaces); len(sub.Thickness) > n-1 {
		sub.Thickness = sub.Thickness[:max(n-1, 0)]
	}
	stop := s.StopIndex - from
	if stop < 0 {
		stop = 0
	}
	if stop > len(sub.Surfaces)-1 {
		stop = len(sub.Surfaces) - 1
	}
	sub.StopIndex = stop
	if err := sub.Validate(); err != nil {
		return nil, err
	}
	return sub, nil
}

// WithImagePlane returns a copy of the system with a flat image plane appended
// at axial position z (mm). The ray can end on this plane, which is what the
// analysis routines of this package (spot diagram, transverse aberration, OPD)
// expect. The plane is placed along the axis of the last element's reference
// frame, so that it lies at global z for systems with coordinate breaks too.
func (s *System) WithImagePlane(z float64) *System {
	out := &System{
		Surfaces:  append([]Surface(nil), s.Surfaces...),
		Thickness: append([]float64(nil), s.Thickness...),
		Mediums:   append([]Medium(nil), s.Mediums...),
		StopIndex: s.StopIndex,
		OriginZ:   s.OriginZ,
		Name:      s.Name,
	}
	c := s.chainStart()
	var ref frame
	for i := range s.Surfaces {
		_, ref = s.place(i, c)
		c = s.advance(ref, i)
	}
	out.Surfaces = append(out.Surfaces, Surface{Profile: Plane{}, Name: "image plane"})
	out.Thickness = append(out.Thickness, (z-ref.o.Z)/ref.r.col(2).Z)
	out.Mediums = append(out.Mediums, s.Mediums[len(s.Mediums)-1])
	return out
}

// Reversed returns the time-reversed system, or nil when the system contains
// coordinate breaks or posed surfaces (a time-reversed chain of tilted frames is
// not provided). For a coaxial system it describes exactly the same
// physical surfaces (the profiles are unchanged, the sag still measured from
// each vertex in +z), traversed in the opposite direction: the surface list is
// reversed, the gaps are listed in reverse order and therefore become negative,
// the media are reversed, the stop index is remapped and OriginZ is set to the
// forward system's last vertex. Tracing the reversed ray through it undoes the
// forward trace, which is the standard way to solve conjugates and double-pass
// systems.
func (s *System) Reversed() *System {
	if s.hasPose() {
		return nil
	}
	n := len(s.Surfaces)
	out := &System{Name: s.Name, StopIndex: n - 1 - s.StopIndex}
	out.OriginZ = s.VertexZ(n - 1)
	out.Surfaces = make([]Surface, n)
	for i := 0; i < n; i++ {
		out.Surfaces[i] = s.Surfaces[n-1-i]
	}
	out.Thickness = make([]float64, len(s.Thickness))
	for j := 0; j < n-1; j++ {
		out.Thickness[j] = -s.Thickness[n-2-j]
	}
	out.Mediums = make([]Medium, len(s.Mediums))
	for i := range s.Mediums {
		out.Mediums[i] = s.Mediums[len(s.Mediums)-1-i]
	}
	return out
}

// finitePositive reports whether x is a positive finite number.
func finitePositive(x float64) bool {
	return x > 0 && !math.IsInf(x, 0) && !math.IsNaN(x)
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
