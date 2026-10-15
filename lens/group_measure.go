package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// ---------------------------------------------------------------------------
// Real-ray measurements of a group in 3D
// ---------------------------------------------------------------------------
//
// All the measurements here launch real rays through the geo.System of the
// group (Group.System), trace them with geo.System.Trace and work on the traced
// exit rays, so they are valid for mirrored, decentred and folded systems for
// which the paraxial engine of package geo is not defined. Every search is
// deterministic: fixed scans, golden-section refinement, no randomness.

// Bundle describes the ray bundle of a real-ray measurement.
type Bundle struct {
	// Object is the object: geo.FieldTan / geo.FieldAngle for an object at
	// infinity, geo.FieldHeight for a finite one.
	Object geo.Field
	// Pupil is the launch plane of the sampling, in the sense of
	// geo.LaunchRay. Leave Radius unset (<= 0) to use the automatic entrance
	// pupil of the group: the clear semi-aperture of the stop element on a
	// plane 1 mm in front of the frontmost surface.
	Pupil geo.Pupil
	// Samples is the pupil sampling in the unit disc; nil selects
	// geo.Hexapolar(3). Sample [0] (the pupil centre) is the chief ray.
	Samples []geo.PupilPoint
	// LambdaNm is the wavelength in nm; 0 selects the d line
	// (geo.Wavelengthd).
	LambdaNm float64
}

// ImageOptions configures the best-focus search of ImagePosition.
type ImageOptions struct {
	// RangeMM is the half-width (mm) of the axial search window, centred on
	// the automatic estimate of the focus; 0 selects the automatic window
	// (half the estimated focus distance plus 10 mm). The window must bracket
	// the minimum of the spot RMS.
	RangeMM float64
	// Scan is the number of coarse samples of the window (default 65, at least
	// 3). The best bracket is then refined by a golden-section search.
	Scan int
	// TolMM is the axial tolerance of the refinement (default 1e-9 mm).
	TolMM float64
}

// defaults fills the zero fields.
func (o ImageOptions) defaults() ImageOptions {
	if o.Scan < 3 {
		o.Scan = 65
	}
	if !(o.TolMM > 0) {
		o.TolMM = 1e-9
	}
	if o.RangeMM < 0 {
		o.RangeMM = 0
	}
	return o
}

// ImageResult is the outcome of ImagePosition.
type ImageResult struct {
	// Point is the image point: the centroid of the traced rays on the best
	// focus plane, in the global frame (mm).
	Point geo.Vec3
	// ChiefHit is the point where the chief ray (the centre sample of the
	// sampling) crosses the best focus plane (mm, global frame) and ChiefDir
	// its exit direction (unit, global frame): together they are the traced
	// chief ray, so a caller can evaluate the image height on any plane.
	ChiefHit geo.Vec3
	ChiefDir geo.Vec3
	// RMS is the RMS radius of the traced rays about Point in the plane (mm);
	// GeoRadius is the largest radius.
	RMS, GeoRadius float64
	// Axis is the unit normal of the best focus plane: the exit axis of the
	// group (the folded axis for a folding mirror). PlaneO is a point of the
	// plane, the vertex of the last surface.
	Axis, PlaneO geo.Vec3
	// FocusMM is the signed distance from PlaneO to the best focus plane along
	// Axis (mm): the focus distance of the group. It is measured along the
	// beam, so it grows away from the object.
	FocusMM float64
	// Rays is the number of rays that reached the last surface and entered the
	// statistics, Failed the number of the sampling that did not.
	Rays, Failed int
}

// exitRay is one traced ray leaving the last surface of the group.
type exitRay struct{ o, d geo.Vec3 }

// prepared is a bundle together with the system it is traced through.
type prepared struct {
	pl      *plan
	sys     *geo.System
	pupil   geo.Pupil
	samples []geo.PupilPoint
	lambda  float64
}

// prepare builds the system of the group and resolves the bundle.
func (g *Group) prepare(b Bundle) (*prepared, error) {
	pl, err := g.plan()
	if err != nil {
		return nil, err
	}
	sys, err := pl.system(g)
	if err != nil {
		return nil, err
	}
	lam := b.LambdaNm
	if lam == 0 {
		lam = geo.Wavelengthd
	}
	if !(lam > 0) || !lensFinite(lam) {
		return nil, fmt.Errorf("lens: group %q: wavelength %v nm is not usable", g.Name, b.LambdaNm)
	}
	samples := b.Samples
	if samples == nil {
		samples = geo.Hexapolar(3)
	}
	if len(samples) == 0 {
		return nil, fmt.Errorf("lens: group %q: empty pupil sampling", g.Name)
	}
	pupil := b.Pupil
	if !(pupil.Radius > 0) || !lensFinite(pupil.Radius) || !lensFinite(pupil.Z) {
		// The automatic entrance pupil: the stop aperture in front of the
		// system.
		e, err := g.element(g.Stop)
		if err != nil {
			return nil, err
		}
		if !(e.Front.SemiAperMM > 0) {
			return nil, fmt.Errorf("lens: group %q: the stop element %d has no clear aperture to size the pupil from; set Bundle.Pupil", g.Name, g.Stop)
		}
		zmin, _ := sys.BoundsZ(10)
		if !lensFinite(zmin) {
			zmin = pl.vertex[0].Z
		}
		pupil = geo.Pupil{Z: zmin - 1, Radius: e.Front.SemiAperMM}
		// The outermost ring of the sampling would graze the rim of the stop
		// surface and be lost to round-off (the pupilInset convention of this
		// package), so the automatic pupil insets the sampling by 0.1%.
		samples = insetSamples(samples, 1-pupilInset)
	}
	if !lensFinite(pupil.Z) || !(pupil.Radius > 0) || !lensFinite(pupil.Radius) {
		return nil, fmt.Errorf("lens: group %q: unusable pupil (z=%v, radius=%v)", g.Name, pupil.Z, pupil.Radius)
	}
	if pupil.Offset == (geo.Vec3{}) {
		// Centre the sampling on the real chief ray. For a coaxial system this
		// is the (on-axis) chief anyway; for a decentred, tilted or folded one
		// it is what makes the centre sample the chief ray. When the chief
		// search fails the sampling stays on the axis (the result is then
		// referred to the chief sample of the sampling, as geo does).
		if off, err := sys.ChiefCrossing(b.Object, pupil.Z, lam); err == nil {
			pupil.Offset = off
		}
	}
	return &prepared{pl: pl, sys: sys, pupil: pupil, samples: samples, lambda: lam}, nil
}

// insetSamples scales a sampling towards the pupil centre.
func insetSamples(pts []geo.PupilPoint, f float64) []geo.PupilPoint {
	out := make([]geo.PupilPoint, len(pts))
	for i, p := range pts {
		out[i] = geo.PupilPoint{PX: p.PX * f, PY: p.PY * f}
	}
	return out
}

// exitRays traces the whole sampling of a field and returns the exit rays and
// the number of samples that did not reach the last surface.
func (p *prepared) exitRays(obj geo.Field) ([]exitRay, int) {
	out := make([]exitRay, 0, len(p.samples))
	failed := 0
	for _, s := range p.samples {
		r := geo.LaunchRay(obj, p.pupil, s, p.lambda)
		if r.Status != geo.StatusOK {
			failed++
			continue
		}
		res := p.sys.Trace(r)
		if !res.OK() {
			failed++
			continue
		}
		out = append(out, exitRay{o: res.Ray.Origin, d: res.Ray.Dir})
	}
	return out, failed
}

// traceOne traces one pupil point and returns its exit ray.
func (p *prepared) traceOne(obj geo.Field, pt geo.PupilPoint) (exitRay, error) {
	r := geo.LaunchRay(obj, p.pupil, pt, p.lambda)
	if r.Status != geo.StatusOK {
		return exitRay{}, fmt.Errorf("lens: degenerate launch ray")
	}
	res := p.sys.Trace(r)
	if !res.OK() {
		return exitRay{}, fmt.Errorf("lens: the ray does not traverse the group (%s at surface %d)", res.Status, res.FailedAt)
	}
	return exitRay{o: res.Ray.Origin, d: res.Ray.Dir}, nil
}

// solve3 solves the 3×3 system a·x = b by Gaussian elimination with partial
// pivoting; ok is false when the system is singular to working precision.
func solve3(a [3][3]float64, b [3]float64) ([3]float64, bool) {
	scale := 0.0
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			scale = math.Max(scale, math.Abs(a[i][j]))
		}
	}
	if scale == 0 {
		return [3]float64{}, false
	}
	m := a
	v := b
	for c := 0; c < 3; c++ {
		p := c
		for r := c + 1; r < 3; r++ {
			if math.Abs(m[r][c]) > math.Abs(m[p][c]) {
				p = r
			}
		}
		if math.Abs(m[p][c]) < 1e-12*scale {
			return [3]float64{}, false
		}
		if p != c {
			m[p], m[c] = m[c], m[p]
			v[p], v[c] = v[c], v[p]
		}
		for r := c + 1; r < 3; r++ {
			f := m[r][c] / m[c][c]
			for k := c; k < 3; k++ {
				m[r][k] -= f * m[c][k]
			}
			v[r] -= f * v[c]
		}
	}
	var x [3]float64
	for r := 2; r >= 0; r-- {
		s := v[r]
		for k := r + 1; k < 3; k++ {
			s -= m[r][k] * x[k]
		}
		x[r] = s / m[r][r]
	}
	return x, true
}

// focusEstimate returns the axial position, measured from o along axis, of the
// point that is closest to every exit ray in the least-squares sense, i.e. the
// point minimizing Σ d(i)² with d(i) the distance from the point to ray i. ok is
// false when the rays are parallel, which is the case of a bundle without a
// focus (an afocal or degenerate system).
func focusEstimate(rays []exitRay, o, axis geo.Vec3) (float64, bool) {
	var m [3][3]float64
	var rhs [3]float64
	for _, r := range rays {
		d := [3]float64{r.d.X, r.d.Y, r.d.Z}
		po := [3]float64{r.o.X, r.o.Y, r.o.Z}
		var pr [3][3]float64
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				pr[i][j] = -d[i] * d[j]
			}
			pr[i][i] += 1 // the projector I − d dᵀ
		}
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				m[i][j] += pr[i][j]
			}
			s := 0.0
			for j := 0; j < 3; j++ {
				s += pr[i][j] * po[j]
			}
			rhs[i] += s
		}
	}
	q, ok := solve3(m, rhs)
	if !ok {
		return 0, false
	}
	return geo.V(q[0]-o.X, q[1]-o.Y, q[2]-o.Z).Dot(axis), true
}

// spotStats intersects every exit ray with the plane through p with unit normal
// axis and returns the centroid of the hits, their RMS radius about the
// centroid and the largest radius (mm). ok is false for a degenerate bundle
// (an empty one, or a ray parallel to the plane).
func spotStats(rays []exitRay, p, axis geo.Vec3) (centroid geo.Vec3, rms, geoR float64, ok bool) {
	if len(rays) == 0 {
		return geo.Vec3{}, 0, 0, false
	}
	hits := make([]geo.Vec3, 0, len(rays))
	var sum geo.Vec3
	for _, r := range rays {
		den := r.d.Dot(axis)
		if den == 0 {
			return geo.Vec3{}, 0, 0, false
		}
		tau := p.Sub(r.o).Dot(axis) / den
		if !lensFinite(tau) {
			return geo.Vec3{}, 0, 0, false
		}
		h := r.o.Add(r.d.Mul(tau))
		hits = append(hits, h)
		sum = sum.Add(h)
	}
	centroid = sum.Div(float64(len(hits)))
	var ss, worst float64
	for _, h := range hits {
		d2 := h.Sub(centroid).Norm2()
		ss += d2
		if d2 > worst {
			worst = d2
		}
	}
	return centroid, math.Sqrt(ss / float64(len(hits))), math.Sqrt(worst), true
}

// goldenRatio is the golden-section step (1/φ).
const goldenRatio = 0.6180339887498949

// searchFocus minimizes the spot RMS over the plane family through o + t·axis.
// The window [lo, hi] is sampled with scan points; the best bracket is then
// refined by a golden-section search to the axial tolerance tol. The returned
// point, RMS and radius belong to the best plane found.
func searchFocus(rays []exitRay, o, axis geo.Vec3, lo, hi float64, scan int, tol float64) (t float64, centroid geo.Vec3, rms, geoR float64, ok bool) {
	rms = math.Inf(1)
	eval := func(tt float64) float64 {
		c, r, gr, good := spotStats(rays, o.Add(axis.Mul(tt)), axis)
		if !good {
			return math.Inf(1)
		}
		if r < rms {
			t, centroid, rms, geoR, ok = tt, c, r, gr, true
		}
		return r
	}
	for i := 0; i < scan; i++ {
		tt := lo
		if scan > 1 {
			tt = lo + (hi-lo)*float64(i)/float64(scan-1)
		}
		eval(tt)
	}
	if !ok {
		return 0, geo.Vec3{}, 0, 0, false
	}
	// Refine the bracket around the best sample.
	step := (hi - lo) / float64(scan-1)
	a := math.Max(lo, t-step)
	b := math.Min(hi, t+step)
	if !(b > a) {
		return t, centroid, rms, geoR, true
	}
	c := b - goldenRatio*(b-a)
	d := a + goldenRatio*(b-a)
	fc, fd := eval(c), eval(d)
	for i := 0; i < 200 && b-a > tol; i++ {
		if fc < fd {
			b, d, fd = d, c, fc
			c = b - goldenRatio*(b-a)
			fc = eval(c)
		} else {
			a, c, fc = c, d, fd
			d = a + goldenRatio*(b-a)
			fd = eval(d)
		}
	}
	return t, centroid, rms, geoR, true
}

// ImagePosition finds the image of a bundle: it traces the sampling of the
// object through the group, then searches the plane perpendicular to the exit
// axis that minimizes the RMS radius of the traced rays (a coarse scan of the
// axial window followed by a golden-section refinement) and returns the image
// point in 3D, the residual RMS and the axis of the focal plane.
//
// The search window is centred on the point that is closest to all the exit
// rays in the least-squares sense (their common focus) and is ImageOptions
// wide. The traced rays are independent of the plane, so only the intersection
// moves during the search.
func (g *Group) ImagePosition(b Bundle, opt ImageOptions) (*ImageResult, error) {
	o := opt.defaults()
	p, err := g.prepare(b)
	if err != nil {
		return nil, err
	}
	rays, failed := p.exitRays(b.Object)
	if len(rays) == 0 {
		return nil, fmt.Errorf("lens: group %q: no ray of the sampling reached the last surface (%d failed)", g.Name, failed)
	}
	anchor := p.pl.exit.o
	axis := p.pl.exit.dirToGlobal(p.pl.beam)
	t0, ok := focusEstimate(rays, anchor, axis)
	if !ok || !lensFinite(t0) || math.Abs(t0) > 1e8 {
		return nil, fmt.Errorf("lens: group %q: the exit bundle has no focus along the exit axis (an afocal or divergent system)", g.Name)
	}
	r := o.RangeMM
	if r == 0 {
		r = 0.5*math.Abs(t0) + 10
	}
	t, centroid, rms, geoR, ok := searchFocus(rays, anchor, axis, t0-r, t0+r, o.Scan, o.TolMM)
	if !ok {
		return nil, fmt.Errorf("lens: group %q: the focus search failed", g.Name)
	}
	res := &ImageResult{
		Point: centroid, RMS: rms, GeoRadius: geoR,
		Axis: axis, PlaneO: anchor, FocusMM: t, Rays: len(rays), Failed: failed,
	}
	if chief, err := p.traceOne(b.Object, p.samples[0]); err == nil {
		plane := anchor.Add(axis.Mul(t))
		den := chief.d.Dot(axis)
		res.ChiefDir = chief.d
		if den != 0 {
			res.ChiefHit = chief.o.Add(chief.d.Mul(plane.Sub(chief.o).Dot(axis) / den))
		} else {
			res.ChiefHit = plane
		}
	} else {
		res.ChiefHit, res.ChiefDir = res.Point, axis
	}
	return res, nil
}

// AngularResult is the outcome of AngularMagnification.
type AngularResult struct {
	// Gamma is the signed angular magnification Γ (see AngularMagnification).
	Gamma float64
	// BeamRatio is D_in/D_out, the ratio of the entrance beam semi-diameter to
	// the exit one, measured on the on-axis bundle: a real-ray cross-check of
	// |Γ| for an afocal group.
	BeamRatio float64
	// EntryDir and ExitDir are the directions of the probe bundle tilted by
	// +delta in the sagittal direction, before and after the group.
	EntryDir, ExitDir geo.Vec3
	// AxisIn and AxisOut are the object-space axis (the global +z) and the exit
	// axis of the group.
	AxisIn, AxisOut geo.Vec3
	// SpreadRad is the largest angle between an exit ray of the on-axis bundle
	// and the exit axis (rad): the collimation diagnostic of the afocal check.
	SpreadRad float64
	// Rays is the number of rays of the on-axis bundle that were used.
	Rays int
}

// angularDelta is the field tangent of the two probe bundles of
// AngularMagnification. It is small enough that the second-order term of the
// measured slope (O(δ²)) is below 1e-12 and large enough that the intersection
// with the exit axis is computed to full precision.
const angularDelta = 1e-6

// AngularMagnification measures the angular magnification Γ of an afocal group
// with real rays.
//
// Convention. The probe bundles are collimated beams of an object at infinity,
// tilted in the sagittal direction x of the object space (geo's field
// convention d = (−TanX, −TanY, 1)/|…|, so a bundle with TanX = −δ travels with
// a positive x component). With s = d·êx / d·êz the slope of a direction in the
// exit frame of the group (êx, êz of Group.ExitFrame, which the folds rotate
// with the beam), Γ is the central difference
//
//	Γ = ( s(+δ) − s(−δ) ) / 2δ ,  δ = 1e-6 .
//
// Γ < 0 means the field is inverted (a Kepler telescope, Γ = −1 for two equal
// lenses); Γ > 0 is erect (a Galileo telescope). The x direction is the
// transverse direction that a fold in the y–z plane (a tilt about x, the usual
// folding plane) leaves unchanged, so Γ is invariant under folding: the same
// telescope with a plane fold mirror in its collimated space has exactly the
// same Γ.
//
// The group must be afocal: the exit bundle of the on-axis field must be
// collimated to better than 1e-2 rad, otherwise an error is returned (a focal
// system has no angular magnification).
func (g *Group) AngularMagnification(b Bundle) (*AngularResult, error) {
	if b.Object.Finite {
		return nil, fmt.Errorf("lens: angular magnification: the object must be at infinity (a collimated bundle)")
	}
	p, err := g.prepare(b)
	if err != nil {
		return nil, err
	}
	axisOut := p.pl.exit.dirToGlobal(p.pl.beam)
	ex := p.pl.exit.ex
	slope := func(d geo.Vec3) (float64, error) {
		den := d.Dot(axisOut)
		if math.Abs(den) < 0.5 {
			return 0, fmt.Errorf("lens: angular magnification: the exit ray is not near the exit axis (d·axis = %g)", den)
		}
		return d.Dot(ex) / den, nil
	}
	// The probe bundles: TanX = −δ has a positive x component.
	probe := func(tanX float64) (exitRay, float64, error) {
		r, err := p.traceOne(geo.FieldTan(tanX, 0), p.samples[0])
		if err != nil {
			return exitRay{}, 0, err
		}
		s, err := slope(r.d)
		if err != nil {
			return exitRay{}, 0, err
		}
		return r, s, nil
	}
	uP, sP, err := probe(-angularDelta)
	if err != nil {
		return nil, err
	}
	_, sM, err := probe(+angularDelta)
	if err != nil {
		return nil, err
	}
	res := &AngularResult{
		Gamma:    (sP - sM) / (2 * angularDelta),
		EntryDir: geo.V(angularDelta, 0, 1).Normalized(),
		ExitDir:  uP.d,
		AxisIn:   geo.V(0, 0, 1),
		AxisOut:  axisOut,
	}
	// Collimation diagnostic and the beam ratio, on the on-axis bundle.
	rays, _ := p.exitRays(geo.FieldTan(0, 0))
	res.Rays = len(rays)
	maxP := 0.0
	for _, s := range p.samples {
		maxP = math.Max(maxP, math.Hypot(s.PX, s.PY))
	}
	dIn := p.pupil.Radius * maxP
	dOut := 0.0
	for _, r := range rays {
		v := r.o.Sub(p.pl.exit.o)
		perp := v.Sub(axisOut.Mul(v.Dot(axisOut))).Norm()
		dOut = math.Max(dOut, perp)
		if r.d.Dot(axisOut) > 0 {
			ang := math.Atan2(r.d.Cross(axisOut).Norm(), r.d.Dot(axisOut))
			res.SpreadRad = math.Max(res.SpreadRad, ang)
		}
	}
	if dIn > 0 && dOut > 0 {
		res.BeamRatio = dIn / dOut
	}
	if res.SpreadRad > 1e-2 {
		return nil, fmt.Errorf("lens: group %q: the exit bundle is not collimated (spread %.3g rad): the group is not afocal, it has no angular magnification", g.Name, res.SpreadRad)
	}
	return res, nil
}

// DispersionResult is the outcome of Dispersion.
type DispersionResult struct {
	// ShiftMM is the axial F−C style focus shift of the best focus planes:
	// FocusMM(lambdaA) − FocusMM(lambdaB) measured along the exit axis from the
	// vertex of the last surface (mm). It is negative when the lambdaA (blue)
	// focus lies in front of the lambdaB (red) one, the sign convention of
	// geo.System.AxialChromatic.
	ShiftMM float64
	// A and B are the two image searches.
	A, B *ImageResult
}

// Dispersion measures the axial chromatic shift of the best focus of the group
// with real rays: the focus planes of the two wavelengths are searched with
// ImagePosition and their axial distance is returned (mm, along the exit axis).
// A wavelength of 0 selects the d line.
//
// Mirrors do not disperse: a group whose media are all constant (a pure mirror
// system, or a system in air) has ShiftMM exactly 0, because every traced ray
// is independent of the wavelength. A plane fold mirror does not change the
// dispersion of the system it folds, because it neither disperses nor aberrates
// a plane wave.
func (g *Group) Dispersion(b Bundle, lambdaA, lambdaB float64, opt ImageOptions) (*DispersionResult, error) {
	ba, bb := b, b
	ba.LambdaNm, bb.LambdaNm = lambdaA, lambdaB
	ra, err := g.ImagePosition(ba, opt)
	if err != nil {
		return nil, err
	}
	rb, err := g.ImagePosition(bb, opt)
	if err != nil {
		return nil, err
	}
	return &DispersionResult{ShiftMM: ra.FocusMM - rb.FocusMM, A: ra, B: rb}, nil
}
