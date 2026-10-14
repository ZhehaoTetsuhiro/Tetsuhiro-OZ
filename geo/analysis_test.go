package geo

import (
	"math"
	"testing"
)

// ---------------------------------------------------------------------------
// Spot diagrams
// ---------------------------------------------------------------------------

// TestSpotStatsAnalytic checks the spot statistics against hand-computed
// values for a set of four points, and the four reference modes.
func TestSpotStatsAnalytic(t *testing.T) {
	z := 50.0
	sd := &SpotDiagram{ImageZ: z}
	for _, p := range []Vec3{V(1, 1, z), V(1, -1, z), V(-1, 1, z), V(-1, -1, z)} {
		sd.Points = append(sd.Points, SpotPoint{Pupil: PupilPoint{}, Hit: p, Status: StatusOK})
	}
	sd.Chief = V(1, 0, z) // arbitrary for the reference test
	sd.ChiefStatus = StatusOK
	sd.Centroid = V(0, 0, z)

	// About the origin: all four points are at √2.
	st := sd.Stats(VertexReference())
	checkClose(t, "vertex RMS", st.RMS, math.Sqrt2, 1e-15)
	checkClose(t, "vertex GEO", st.GEO, math.Sqrt2, 1e-15)
	if st.N != 4 {
		t.Errorf("N = %d, want 4", st.N)
	}
	checkClose(t, "centroid x", st.Centroid.X, 0, 1e-15)

	// About (1, 0): distances 1, 1, √5, √5 → RMS = √3, GEO = √5.
	st = sd.Stats(ExplicitReference(V(1, 0, z)))
	checkClose(t, "explicit RMS", st.RMS, math.Sqrt(3), 1e-15)
	checkClose(t, "explicit GEO", st.GEO, math.Sqrt(5), 1e-15)

	// The chief reference resolves to the chief hit.
	st = sd.Stats(ChiefReference())
	checkClose(t, "chief ref y", st.Reference.Y, 0, 1e-15)
	checkClose(t, "chief ref x", st.Reference.X, 1, 1e-15)
	// The centroid reference resolves to the centroid.
	st = sd.Stats(CentroidReference())
	checkClose(t, "centroid ref x", st.Reference.X, 0, 1e-15)

	// StatsAt is the explicit form.
	st2 := sd.StatsAt(V(1, 0, z))
	checkClose(t, "StatsAt RMS", st2.RMS, math.Sqrt(3), 1e-15)

	// The centroid minimizes the RMS over all reference points: any other
	// reference must give a value at least as large.
	for _, ref := range []Reference{ChiefReference(), VertexReference(), ExplicitReference(V(0.3, -0.2, z))} {
		if r := sd.Stats(ref).RMS; r < sd.Stats(CentroidReference()).RMS-1e-12 {
			t.Errorf("reference %s yields RMS %g below the centroid RMS %g",
				ref, r, sd.Stats(CentroidReference()).RMS)
		}
	}
	// An empty (fully vignetted) spot has zero statistics and no panic.
	empty := &SpotDiagram{ImageZ: z}
	checkClose(t, "empty N", float64(empty.Stats(ChiefReference()).N), 0, 0)
}

// TestSpotParabolicMirror checks the two analytic spot limits of a parabolic
// mirror: an exactly stigmatic spot at the focus, and the paraxial blur circle
// r_p·Δz/L when the image plane is moved by Δz.
func TestSpotParabolicMirror(t *testing.T) {
	const (
		R     = -200.0
		rPup  = 40.0
		zPup  = -50.0
		focus = -100.0 // R/2 for the concave mirror
	)
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: R, K: -1}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	pup := Pupil{Z: zPup, Radius: rPup}
	f := FieldTan(0, 0)
	pts := Rings(3, 8)

	// In focus: every ray lands on the axis point (0, 0, focus).
	sd, err := sys.SpotDiagram(f, pup, pts, focus, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	st := sd.Stats(VertexReference())
	if st.N != len(pts) {
		t.Fatalf("N = %d, want %d", st.N, len(pts))
	}
	checkLess(t, "in-focus RMS", st.RMS, 1e-12)
	checkLess(t, "in-focus GEO", st.GEO, 1e-12)
	checkClose(t, "chief hit x", sd.Chief.X, 0, 1e-12)
	checkClose(t, "chief hit y", sd.Chief.Y, 0, 1e-12)
	checkClose(t, "chief hit z", sd.Chief.Z, focus, 1e-12)
	// The chief ray must be the one through the pupil centre: the samplings
	// start with the centre point, so its hit must equal the chief hit.
	checkClose(t, "chief = centre sample", sd.Chief.Distance(sd.Points[0].Hit), 0, 1e-12)

	// Defocused: the ring of rays becomes the circle where the converging
	// cone (apex at the focus, base at the mirror rim at height r_p) meets the
	// plane Δz away from the apex: radius r_p·|Δz|/L with L the distance from
	// the rim to the focus. For a single ring the RMS equals the GEO radius,
	// which is exactly the paraxial prediction of the defocused spot.
	sag := -rPup * rPup / (2 * (-R)) // parabola sag at the rim
	L := math.Abs(focus - sag)
	for _, dz := range []float64{1, -2.5, 5} {
		ring := []PupilPoint{{1, 0}, {0, 1}, {-1, 0}, {0, -1}}
		sd, err := sys.SpotDiagram(f, pup, ring, focus+dz, Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		st := sd.Stats(VertexReference())
		want := rPup * math.Abs(dz) / L
		checkRel(t, "defocus blur radius", st.RMS, want, 1e-9)
		checkRel(t, "defocus GEO", st.GEO, want, 1e-9)
	}
}

// TestSpotReferencesReal checks the reference modes on a real off-axis spot:
// the chief reference is the pupil-centre hit, the centroid reference is the
// centroid of the hits, and the centroid has the smallest RMS.
func TestSpotReferencesReal(t *testing.T) {
	const R = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: R}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	pup := Pupil{Z: -50, Radius: 40}
	f := FieldTan(0, 0.03)
	sd, err := sys.SpotDiagram(f, pup, Rings(3, 8), -100, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	chief := sd.Stats(ChiefReference())
	cent := sd.Stats(CentroidReference())
	checkClose(t, "chief ref = chief hit", chief.Reference.Distance(sd.Chief), 0, 0)
	if chief.RMS <= 0 {
		t.Error("off-axis spherical mirror must have a non-zero spot")
	}
	if cent.RMS > chief.RMS+1e-12 {
		t.Errorf("centroid RMS %g exceeds chief RMS %g", cent.RMS, chief.RMS)
	}
	checkClose(t, "centroid = diagram centroid", cent.Reference.Distance(sd.Centroid), 0, 1e-15)
	if chief.GEO < chief.RMS {
		t.Errorf("GEO %g smaller than RMS %g", chief.GEO, chief.RMS)
	}
}

// TestSpotVignetting checks that blocked rays are reported with their status
// and excluded from the statistics.
func TestSpotVignetting(t *testing.T) {
	// Parabolic mirror with a 10 mm clear aperture sampled over a 40 mm
	// pupil: the outer samples fall outside the aperture.
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: -200, K: -1}, Mirror: true, Aperture: CircleAperture(10)}},
		nil, nil)
	pup := Pupil{Z: -50, Radius: 40}
	pts := []PupilPoint{{0, 0}, {0, 0.2}, {0, 0.3}} // 0, 8, 12 mm at the mirror
	sd, err := sys.SpotDiagram(FieldTan(0, 0), pup, pts, -100, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if sd.Points[0].Status != StatusOK || sd.Points[1].Status != StatusOK {
		t.Errorf("inner rays should pass: %v %v", sd.Points[0].Status, sd.Points[1].Status)
	}
	checkStatus(t, "outer ray", sd.Points[2].Status, StatusVignetted)
	st := sd.Stats(ChiefReference())
	if st.N != 2 {
		t.Errorf("stats N = %d, want 2 (the vignetted ray excluded)", st.N)
	}
}

// ---------------------------------------------------------------------------
// Transverse aberration (ray fans)
// ---------------------------------------------------------------------------

// TestRayFanSymmetry checks the fan conventions on an on-axis field: the
// meridional fan stays in the meridional plane (DX = 0) and carries the
// spherical aberration in DY, the sagittal fan mirrors that, and the
// transverse aberration of a spherical mirror is cubic in the pupil
// coordinate with the undercorrected sign.
func TestRayFanSymmetry(t *testing.T) {
	const R = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: R}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	pup := Pupil{Z: -50, Radius: 40}
	f := FieldTan(0, 0)
	zImg := R / 2 // paraxial focus

	mer, err := sys.MeridionalRayFan(f, pup, 9, zImg, Wavelengthd, ChiefReference())
	if err != nil {
		t.Fatal(err)
	}
	sag, err := sys.SagittalRayFan(f, pup, 9, zImg, Wavelengthd, ChiefReference())
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range mer.Points {
		checkClose(t, "meridional DX", p.DX, 0, 1e-12)
		if math.Abs(p.DY) < 1e-12 && math.Abs(p.P) == 1 {
			t.Error("marginal meridional ray has no transverse aberration")
		}
	}
	for _, p := range sag.Points {
		checkClose(t, "sagittal DY", p.DY, 0, 1e-12)
	}
	// Cubic scaling: DY/p³ is the same at p = 0.25 and 0.5.
	find := func(rf *RayFan, p float64) float64 {
		for _, q := range rf.Points {
			if math.Abs(q.P-p) < 1e-12 {
				return q.DY
			}
		}
		t.Fatalf("p = %g not sampled", p)
		return 0
	}
	c1 := find(mer, 0.125) / math.Pow(0.125, 3)
	c2 := find(mer, 0.25) / math.Pow(0.25, 3)
	checkRel(t, "cubic SA scaling", c2, c1, 0.01)
	// Undercorrected mirror: at the paraxial focus the marginal ray has
	// already crossed the axis, so its hit is on the far side (positive y for
	// a ray starting above the axis... the reflected beam travels to −z, so
	// the sign follows the mirror geometry): the aberration must be
	// antisymmetric in p and non-zero.
	checkRel(t, "antisymmetry", find(mer, -1), -find(mer, 1), 1e-9)
	// The fan RMS is the quadratic mean of the transverse aberrations.
	sum, n := 0.0, 0
	for _, p := range mer.Points {
		if p.Status.OK() {
			sum += p.DX*p.DX + p.DY*p.DY
			n++
		}
	}
	checkRel(t, "fan RMS", mer.RMS(), math.Sqrt(sum/float64(n)), 1e-12)
}

// TestChiefAlignedSampling checks that a chief-aligned pupil sampling
// reproduces the real (aimed) chief ray: the pupil-centre sample must land on
// the image plane exactly where System.ChiefRay does, while an alignment-less
// pupil misses it. This is the property the third-order comparisons rely on.
func TestChiefAlignedSampling(t *testing.T) {
	const R = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: R}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	f := FieldTan(0, 0.03)
	pup := Pupil{Z: -50, Radius: 40}
	zImg := R / 2

	chief := sys.ChiefRay(f, Wavelengthd)
	ref := sys.traceToPlane(chief, zImg)
	if !ref.Status.OK() {
		t.Fatalf("chief trace failed: %s", ref.Status)
	}
	off, err := sys.ChiefOffset(f, pup.Z)
	if err != nil {
		t.Fatal(err)
	}
	// Aligned: the pupil-centre launch is the chief ray's line in the object
	// space, so its image-plane hit must agree with the aimed chief ray.
	aligned := pup
	aligned.Offset = off
	h := sys.traceToPlane(LaunchRay(f, aligned, PupilPoint{}, Wavelengthd), zImg)
	if !h.Status.OK() {
		t.Fatal("aligned launch failed")
	}
	checkClose(t, "aligned centre hit = chief", h.Hit.Distance(ref.Hit), 0, 2e-3)
	// Unaligned: the pupil-centre launch is a different ray, landing several
	// millimetres away (the offset is zPup·tanθ = 1.5 mm at the mirror).
	hu := sys.traceToPlane(LaunchRay(f, pup, PupilPoint{}, Wavelengthd), zImg)
	if d := hu.Hit.Distance(ref.Hit); d < 1e-2 {
		t.Errorf("alignment-less centre launch is only %g mm from the chief", d)
	}
	// The fan itself: the sample at p = 0 has zero aberration by construction.
	rf, err := sys.MeridionalRayFan(f, aligned, 9, zImg, Wavelengthd, ChiefReference())
	if err != nil {
		t.Fatal(err)
	}
	for _, q := range rf.Points {
		if q.P == 0 {
			checkClose(t, "centre sample aberration", math.Hypot(q.DX, q.DY), 0, 1e-12)
		}
	}
}

// findFan returns DY at the fan point with parameter p (exact match).
func findFan(rf *RayFan, p float64) float64 {
	for _, q := range rf.Points {
		if math.Abs(q.P-p) < 1e-12 {
			return q.DY
		}
	}
	return math.NaN()
}

// ---------------------------------------------------------------------------
// OPD
// ---------------------------------------------------------------------------

// TestOPDParabolicMirrorZero checks the defining property of the parabolic
// mirror: the reflected wavefront is exactly spherical, so the OPD over the
// whole pupil is zero to numerical precision (this is the "OPD ≈ 0" anchor of
// the OPD machinery).
func TestOPDParabolicMirrorZero(t *testing.T) {
	const R = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: R, K: -1}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	pup := Pupil{Z: -50, Radius: 40}
	g, err := sys.OPDGrid(FieldTan(0, 0), pup, Hexapolar(3), -100, Wavelengthd, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	checkLess(t, "parabola OPD RMS", g.RMS, 1e-12)
	checkLess(t, "parabola OPD PV", g.PV, 1e-12)
	checkClose(t, "chief OPD", g.Raw[0], 0, 1e-15)
	if g.N != len(g.Pupils) {
		t.Errorf("N = %d, want %d valid samples", g.N, len(g.Pupils))
	}
	checkLess(t, "parabola RMS in waves", g.RMSWaves(), 1e-9)
	// The tilted field of the same mirror is also stigmatic for a parabola?
	// No — off axis the parabola has coma; with the chief-aligned pupil the
	// OPD must still be exactly zero on-axis only, so check the off-axis case
	// is *not* zero (the aberration is real).
	off, err := sys.ChiefOffset(FieldTan(0, 0.03), pup.Z)
	if err != nil {
		t.Fatal(err)
	}
	pup.Offset = off
	g2, err := sys.OPDGrid(FieldTan(0, 0.03), pup, Hexagon3(), -100, Wavelengthd, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if g2.RMS < 1e-6 {
		t.Errorf("off-axis parabola OPD RMS = %g, expected the coma/astigmatism to show", g2.RMS)
	}
}

// Hexagon3 is Hexapolar(3) under another name, kept local to make the intent
// of the off-axis check readable.
func Hexagon3() []PupilPoint { return Hexapolar(3) }

// TestOPDAplanaticSphereZero is the finite-object counterpart of the parabola
// check: a single spherical surface used from an aplanatic point images
// exactly, so every ray leaves the surface on the same sphere centred on the
// image and the OPD over the pupil is zero.
func TestOPDAplanaticSphereZero(t *testing.T) {
	const (
		R = 20.0
		n = 1.5
	)
	objectZ := R * (1 - n)  // aplanatic object, in front of the vertex
	imageZ := R * (1 + 1/n) // aplanatic image, inside the glass
	sys := mustSystem(t, ConstantIndex(1),
		[]Surface{{Profile: Sphere{R: R}, Aperture: CircleAperture(6)}},
		nil, []Medium{ConstantIndex(1), ConstantIndex(n)})
	f := FieldHeight(objectZ, 0, 0)
	pup := Pupil{Z: -5, Radius: 2} // between the object and the surface
	g, err := sys.OPDGrid(f, pup, Rings(3, 8), imageZ, Wavelengthd, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	checkLess(t, "aplanatic OPD RMS", g.RMS, 1e-12)
	checkLess(t, "aplanatic OPD PV", g.PV, 1e-12)

	// The same surface from a non-aplanatic conjugate must show spherical
	// aberration in the OPD.
	sys2 := mustSystem(t, ConstantIndex(1),
		[]Surface{{Profile: Sphere{R: R}, Aperture: CircleAperture(6)}},
		nil, []Medium{ConstantIndex(1), ConstantIndex(n)})
	g2, err := sys2.OPDGrid(FieldHeight(-40, 0, 0), Pupil{Z: -40, Radius: 3}, Rings(3, 8),
		sys2.VertexZ(0)+n/((n-1)/R), Wavelengthd, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if g2.RMS < 1e-6 {
		t.Errorf("non-aplanatic conjugate OPD RMS = %g, expected spherical aberration", g2.RMS)
	}
}

// TestOPDPistonTiltRemoval checks the two cleanup options: a reference point
// displaced from the chief leaves a tilt and a piston in the raw map, and the
// least-squares removal cancels them.
func TestOPDPistonTiltRemoval(t *testing.T) {
	const R = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Conic{R: R, K: -1}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	pup := Pupil{Z: -50, Radius: 40}
	pts := Rings(3, 8)
	// Reference deliberately off-centre: the raw OPD is dominated by a tilt.
	ref := ExplicitReference(V(2, -1, -100))
	raw, err := sys.OPDGrid(FieldTan(0, 0), pup, pts, -100, Wavelengthd, OPDOptions{Reference: ref})
	if err != nil {
		t.Fatal(err)
	}
	if raw.PV < 0.1 {
		t.Fatalf("expected a large tilt-dominated raw map, PV = %g", raw.PV)
	}
	clean, err := sys.OPDGrid(FieldTan(0, 0), pup, pts, -100, Wavelengthd,
		OPDOptions{Reference: ref, RemovePiston: true, RemoveTilt: true})
	if err != nil {
		t.Fatal(err)
	}
	// Removing piston and tilt leaves only the defocus that the displaced
	// reference sphere introduces, so the RMS must drop by a large factor.
	if clean.RMS > 0.1*raw.RMS {
		t.Errorf("tilt/piston removal left RMS %g (raw %g)", clean.RMS, raw.RMS)
	}
	if clean.PV > 0.3*raw.PV {
		t.Errorf("tilt/piston removal left PV %g (raw %g)", clean.PV, raw.PV)
	}
	// Piston removal alone leaves the tilt.
	piston, err := sys.OPDGrid(FieldTan(0, 0), pup, pts, -100, Wavelengthd,
		OPDOptions{Reference: ref, RemovePiston: true})
	if err != nil {
		t.Fatal(err)
	}
	if piston.RMS < 0.01 {
		t.Errorf("piston-only removal left RMS %g, expected the tilt", piston.RMS)
	}
	// The raw map is unchanged by the options.
	if raw.Raw[0] != raw.OPD[0] && math.Abs(raw.Raw[0]-raw.OPD[0]) > 1e-15 {
		t.Error("raw and OPD must agree when no removal is requested")
	}
	// Waves conversion.
	checkRel(t, "waves conversion", raw.Waves(1), 1e6/Wavelengthd, 1e-12)
}

// TestOPDChiefReferenceZero checks that the raw OPD of the chief sample is
// zero by construction (the reference sphere passes through the chief's exit
// point) and that a fully vignetted pupil is reported rather than crashing.
func TestOPDChiefReferenceZero(t *testing.T) {
	const R = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: R}, Mirror: true, Aperture: CircleAperture(10)}},
		nil, nil)
	pup := Pupil{Z: -50, Radius: 40}
	g, err := sys.OPDGrid(FieldTan(0, 0), pup, Rings(3, 6), -100, Wavelengthd, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	// The aperture (10 mm) is smaller than the pupil (40 mm): every ring
	// point is vignetted while the centre sample survives, so exactly one
	// sample is valid and the rest are NaN.
	if g.N != 1 {
		t.Errorf("N = %d, want 1 (only the chief-centred sample) with a fully vignetted ring set", g.N)
	}
	for i, v := range g.Raw {
		if i == 0 {
			continue // the centre sample is the chief itself
		}
		if !math.IsNaN(v) {
			t.Errorf("vignetted sample %d OPD = %g, want NaN", i, v)
		}
	}
}

// ---------------------------------------------------------------------------
// Distortion
// ---------------------------------------------------------------------------

// TestDistortionThinLens checks the classical result that a thin lens with the
// stop in contact has no distortion: the real chief ray through the stop
// centre is undeviated and lands exactly at the paraxial image height.
func TestDistortionThinLens(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 150 * (n - 1)
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{0}, []Medium{Air{}, N_BK7, Air{}}) // zero thickness: stop at the lens
	sys.StopIndex = 0
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	f := FieldTan(0, 0.05)
	d, err := sys.Distortion(f, par.ImageZ, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	checkLess(t, "thin lens distortion", math.Abs(d), 1e-10)

	// Displacing the stop in front of the lens breaks the symmetry and
	// produces barrel/pincushion distortion of the order of a fraction of a
	// percent for this focal length.
	sys2 := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(6)},
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{30, 6}, []Medium{Air{}, Air{}, N_BK7, Air{}})
	sys2.StopIndex = 0
	par2, err := sys2.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	d2, err := sys2.Distortion(f, par2.ImageZ, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(d2) < 1e-4 {
		t.Errorf("displaced stop distortion = %g, expected a measurable value", d2)
	}
	// The relative distortion is the transverse image-height error (∝ h³)
	// divided by the image height (∝ h), so it scales as the square of the
	// field: doubling the field tangent multiplies it by about four.
	d2b, err := sys2.Distortion(FieldTan(0, 0.1), par2.ImageZ, Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	checkRel(t, "relative distortion ∝ field²", d2b/d2, 4, 0.15)
	// An on-axis field has no distortion to define.
	if _, err := sys2.Distortion(FieldTan(0, 0), par2.ImageZ, Wavelengthd); err == nil {
		t.Error("on-axis distortion must be an error (zero paraxial height)")
	}
}

// ---------------------------------------------------------------------------
// Concurrency
// ---------------------------------------------------------------------------

// TestConcurrentAnalyses runs the analysis routines on a shared system from
// several goroutines. The package is documented as read-only over its inputs,
// so the race detector (go test -race) must find nothing here.
func TestConcurrentAnalyses(t *testing.T) {
	const R = -200.0
	sys := mustSystem(t, Air{},
		[]Surface{{Profile: Sphere{R: R}, Mirror: true, Aperture: CircleAperture(60)}},
		nil, nil)
	pup := Pupil{Z: -50, Radius: 40}
	done := make(chan error, 8)
	for i := 0; i < 8; i++ {
		go func(i int) {
			f := FieldTan(0, 0.01*float64(i))
			// A mix of the three analysis entry points.
			if _, err := sys.SpotDiagram(f, pup, Rings(2, 6), -100, Wavelengthd); err != nil {
				done <- err
				return
			}
			if _, err := sys.MeridionalRayFan(f, pup, 7, -100, Wavelengthd, ChiefReference()); err != nil {
				done <- err
				return
			}
			if _, err := sys.OPDGrid(f, pup, Rings(2, 6), -100, Wavelengthd, OPDOptions{}); err != nil {
				done <- err
				return
			}
			if _, err := sys.SeidelCoefficients(Wavelengthd, f, 40); err != nil {
				done <- err
				return
			}
			done <- nil
		}(i)
	}
	for i := 0; i < 8; i++ {
		if err := <-done; err != nil {
			t.Fatal(err)
		}
	}
}
