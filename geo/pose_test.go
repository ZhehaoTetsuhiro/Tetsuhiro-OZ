package geo

import (
	"math"
	"math/rand"
	"testing"
)

// posePtr returns a pointer to a copy of cb (the Surface.Pose convention).
func posePtr(cb CoordBreak) *CoordBreak { return &cb }

// Acceptance tests of the decentre/tilt machinery (see the report of the
// geo-frame task). Every numerical threshold is justified in the test text.

const testLambda = 587.5618 // d-line (nm)

// decentredLens is a biconvex thin lens (R = ±100 mm, n = 1.5, f = 100 mm)
// whose two surfaces, both at zero thickness, are decentred by dx and tilted by
// tiltY (zero for a pure decentre). The first surface is the stop.
func decentredLens(t *testing.T, dx, tiltY float64) *System {
	t.Helper()
	pose := CoordBreak{Dx: dx, TiltY: tiltY}
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 100}, Aperture: CircleAperture(15), Pose: posePtr(pose), Name: "L front"},
			{Profile: Sphere{R: -100}, Aperture: CircleAperture(15), Pose: posePtr(pose), Name: "L rear"},
		},
		[]float64{0},
		[]Medium{Air{}, ConstantIndex(1.5), Air{}})
}

// nominalLens is decentredLens without any pose.
func nominalLens(t *testing.T) *System { return decentredLens(t, 0, 0) }

// TestTiltedPlaneMirrorDeflectsByTwoTheta is acceptance (a): an axial ray
// reflected by a plane mirror tilted by θ about x (or y) comes back deflected by
// 2θ from the reversed incident direction. The analytic reflected direction is
// d' = d − 2(d·n)n with n = Q·ẑ = (0, −sinθ, cosθ) (tilt about x) and
// (−sin2θ, 0, −cos2θ) for the tilt about y. Threshold 1e-12 rad on the angle
// (computed with atan2 of cross and dot products, so it is not limited by
// arccos conditioning near 0).
func TestTiltedPlaneMirrorDeflectsByTwoTheta(t *testing.T) {
	for _, about := range []string{"x", "y"} {
		for _, th := range []float64{0.05, 0.3, 0.7} {
			pose := CoordBreak{}
			if about == "x" {
				pose.TiltX = th
			} else {
				pose.TiltY = th
			}
			sys := mustSystem(t, Air{}, []Surface{{Profile: Plane{}, Mirror: true, Pose: posePtr(pose)}}, nil, nil)
			in := NewRay(V(0, 0, -10), V(0, 0, 1), testLambda)
			res := sys.Trace(in)
			if !res.OK() {
				t.Fatalf("tilt %s %.2f: trace failed: %v", about, th, res.Status)
			}
			out := res.Ray.Dir
			back := in.Dir.Neg()
			angle := math.Atan2(out.Cross(back).Norm(), out.Dot(back))
			t.Logf("tilt %s θ=%.2f: |deflection − 2θ| = %.3g rad", about, th, math.Abs(angle-2*th))
			checkClose(t, "deflection 2θ ("+about+")", angle, 2*th, 1e-12)
			var want Vec3
			if about == "x" {
				want = V(0, math.Sin(2*th), -math.Cos(2*th))
			} else {
				want = V(-math.Sin(2*th), 0, -math.Cos(2*th))
			}
			if d := out.Sub(want).Norm(); d > 1e-12 {
				t.Errorf("tilt %s %.2f: reflected direction %v, want %v (|diff| %.3g)", about, th, out, want, d)
			}
			// The mirror passes through the decentred vertex: the hit is the origin.
			if h := res.Steps[0].Point; h.Norm() > 1e-12 {
				t.Errorf("tilt %s %.2f: hit point %v, want origin", about, th, h)
			}
		}
	}
}

// TestDecentredThinLensImageShiftMatchesAnalytic is acceptance (b), lens part.
// A thin lens decentred by dx (both surfaces, zero thickness) has the chief ray
// through its vertex, which is undeviated, so the chief image point is exact:
//
//	object at infinity:   Δx = dx          (image point moves with the lens),
//	finite object:        Δx = dx·(1 − m)  (m = nominal transverse magnification).
//
// The image distance and m come from the nominal paraxial trace. Threshold 1e-9 mm.
func TestDecentredThinLensImageShiftMatchesAnalytic(t *testing.T) {
	const dx = 0.3
	// Object at infinity: EFL 100 mm, image plane at the focus.
	f := nominalLens(t)
	parInf, err := f.Paraxial(testLambda, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	dec := decentredLens(t, dx, 0)
	chief := dec.ChiefRay(FieldTan(0, 0), testLambda)
	if chief.Status != StatusOK {
		t.Fatalf("chief ray of the decentred lens: status %v", chief.Status)
	}
	h := dec.traceToPlane(chief, parInf.ImageZ)
	if !h.Status.OK() {
		t.Fatalf("chief ray does not reach the focal plane: %v", h.Status)
	}
	t.Logf("infinite object: |Δx − dx| = %.3g mm", math.Abs(h.Hit.X-dx))
	checkClose(t, "infinite object: image shift x", h.Hit.X, dx, 1e-9)
	checkClose(t, "infinite object: image shift y", h.Hit.Y, 0, 1e-9)

	// Finite object at 200 mm: image at 200 mm, m = -1, so the shift is 2·dx.
	parFin, err := f.Paraxial(testLambda, -200)
	if err != nil {
		t.Fatal(err)
	}
	obj := FieldHeight(-200, 0, 0)
	chiefF := dec.ChiefRay(obj, testLambda)
	if chiefF.Status != StatusOK {
		t.Fatalf("chief ray of the finite-object field: status %v", chiefF.Status)
	}
	hf := dec.traceToPlane(chiefF, parFin.ImageZ)
	if !hf.Status.OK() {
		t.Fatalf("finite-object chief does not reach the image plane: %v", hf.Status)
	}
	want := dx * (1 - parFin.Magnification)
	t.Logf("finite object (m=%.6g): |Δx − dx(1−m)| = %.3g mm", parFin.Magnification, math.Abs(hf.Hit.X-want))
	checkClose(t, "finite object: magnification", parFin.Magnification, -1, 1e-9)
	checkClose(t, "finite object: image shift x = dx(1-m)", hf.Hit.X, want, 1e-9)
}

// TestDecentredPlateLateralShiftIsInvariant is acceptance (b), plate part. A
// plane-parallel plate is translation invariant: decentring it leaves the ray
// geometry unchanged, and an oblique ray leaves it displaced by the classical
// lateral shift D = T·sin(θ − θ')/cos θ' with sin θ' = sin θ / n. Both the
// invariance (between dx = 0 and dx = 3 mm) and the closed form are checked to
// 1e-9 mm.
func TestDecentredPlateLateralShiftIsInvariant(t *testing.T) {
	const (
		T     = 10.0
		n     = 1.5
		theta = 0.3
	)
	plate := func(dx float64) *System {
		return mustSystem(t, Air{},
			[]Surface{
				{Profile: Plane{}, Aperture: CircleAperture(50), Pose: posePtr(CoordBreak{Dx: dx})},
				{Profile: Plane{}, Aperture: CircleAperture(50), Pose: posePtr(CoordBreak{Dx: dx})},
			},
			[]float64{T},
			[]Medium{Air{}, ConstantIndex(n), Air{}})
	}
	d := V(math.Sin(theta), 0, math.Cos(theta))
	in := NewRay(V(0, 0, -10), d, testLambda)
	var shifts [2]float64
	for k, dx := range []float64{0, 3} {
		res := plate(dx).Trace(in)
		if !res.OK() || len(res.Steps) != 2 {
			t.Fatalf("decentre %.1f: trace failed (%v, %d steps)", dx, res.Status, len(res.Steps))
		}
		pin := res.Steps[0].Point
		pout := res.Steps[1].Point
		if e := res.Ray.Dir.Sub(d).Norm(); e > 1e-12 {
			t.Errorf("decentre %.1f: emergent direction differs from incident by %.3g", dx, e)
		}
		shifts[k] = pout.Sub(pin).Cross(d).Norm()
	}
	checkClose(t, "lateral shift invariant under decentre", shifts[1], shifts[0], 1e-9)
	sinT := math.Sin(theta) / n
	tp := math.Asin(sinT)
	want := T * math.Sin(theta-tp) / math.Cos(tp)
	t.Logf("plate: D(dx=0)=%.12g, D(dx=3)=%.12g, analytic %.12g, |D−analytic|=%.3g mm",
		shifts[0], shifts[1], want, math.Abs(shifts[0]-want))
	checkClose(t, "plate lateral shift D = T sin(θ−θ')/cosθ'", shifts[0], want, 1e-9)
}

// TestTranslatedSystemIsTranslatedImage checks the analysis references on a
// posed system: when the whole singlet (all surfaces) and the pupil are
// translated by (d, e), the spot about the chief reference is unchanged, its
// centroid moves by (d, e), and the chief sample of the OPD is zero. This
// exercises the posed chief ray, the pupil offset, the chief OPL referred to
// the pupil plane (a missing n·t correction would make Raw[0] = n·t ≠ 0) and
// the OPD reference. Thresholds 1e-10 mm (sums of ~50 traced rays).
func TestTranslatedSystemIsTranslatedImage(t *testing.T) {
	const d, e = 0.8, -0.45
	bk7, _ := Glass("N-BK7")
	base := func(pose CoordBreak) *System {
		return mustSystem(t, Air{},
			[]Surface{
				{Profile: Sphere{R: 50}, Aperture: CircleAperture(12.5), Pose: posePtr(pose)},
				{Profile: Sphere{R: -50}, Aperture: CircleAperture(12.5), Pose: posePtr(pose)},
			},
			[]float64{5},
			[]Medium{Air{}, bk7, Air{}})
	}
	plain := base(CoordBreak{})
	moved := base(CoordBreak{Dx: d, Dy: e})
	field := FieldAngle(0, 2)
	pup := Pupil{Z: -20, Radius: 8}
	par, err := plain.Paraxial(testLambda, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	imageZ := par.ImageZ
	// Centre both samplings on the real chief ray (for a coaxial system the
	// pupil offset of the chief is ChiefOffset; for the posed system it is
	// ChiefCrossing). The chief crosses the stop (the first surface) at its
	// centre, so the two crossings differ by exactly (d, e).
	offP, err := plain.ChiefOffset(field, pup.Z)
	if err != nil {
		t.Fatal(err)
	}
	off, err := moved.ChiefCrossing(field, pup.Z, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	checkClose(t, "chief crossing x shifted by d", off.X-offP.X, d, 1e-9)
	checkClose(t, "chief crossing y shifted by e", off.Y-offP.Y, e, 1e-9)
	pupP := pup
	pupP.Offset = offP
	pupM := pup
	pupM.Offset = off

	pts := Hexapolar(4)
	sdA, err := plain.SpotDiagram(field, pupP, pts, imageZ, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	sdB, err := moved.SpotDiagram(field, pupM, pts, imageZ, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	stA := sdA.Stats(ChiefReference())
	stB := sdB.Stats(ChiefReference())
	worstPt := 0.0
	for i := range pts {
		shift := sdB.Points[i].Hit.Sub(sdA.Points[i].Hit)
		worstPt = math.Max(worstPt, math.Hypot(shift.X-d, shift.Y-e))
	}
	t.Logf("translation: worst image-point deviation from the (d,e) shift = %.3g mm", worstPt)
	checkClose(t, "image points shifted by (d,e)", worstPt, 0, 1e-10)
	if stA.N != len(pts) || stB.N != len(pts) {
		t.Fatalf("expected all %d rays to reach the image plane, got %d and %d", len(pts), stA.N, stB.N)
	}
	t.Logf("translation: |RMS_B−RMS_A|=%.3g mm, centroid shift−(d,e)=(%.3g,%.3g) mm",
		math.Abs(stB.RMS-stA.RMS), stB.Centroid.X-stA.Centroid.X-d, stB.Centroid.Y-stA.Centroid.Y-e)
	checkClose(t, "spot RMS invariant under translation", stB.RMS, stA.RMS, 1e-10)
	checkClose(t, "spot centroid x shifted by d", stB.Centroid.X-stA.Centroid.X, d, 1e-10)
	checkClose(t, "spot centroid y shifted by e", stB.Centroid.Y-stA.Centroid.Y, e, 1e-10)
	if stA.RMS <= 0 {
		t.Fatalf("degenerate test: zero spot RMS")
	}

	gA, err := plain.OPDGrid(field, pupP, pts, imageZ, testLambda, OPDOptions{Reference: ChiefReference()})
	if err != nil {
		t.Fatal(err)
	}
	gB, err := moved.OPDGrid(field, pupM, pts, imageZ, testLambda, OPDOptions{Reference: ChiefReference()})
	if err != nil {
		t.Fatal(err)
	}
	checkClose(t, "OPD RMS invariant under translation", gB.RMS, gA.RMS, 1e-10)
	checkClose(t, "OPD of the chief sample is zero (posed)", gB.Raw[0], 0, 1e-12)
	checkClose(t, "OPD of the chief sample is zero (coaxial)", gA.Raw[0], 0, 1e-12)
}

// TestRollInvarianceOfSpotAndOPD rolls every surface about the optical axis by
// γ and rotates the field direction and the pupil sampling by the same γ. A
// singlet is symmetric under rotation about its axis, so the spot RMS and the
// OPD RMS must be unchanged (1e-10 mm). The roll enters every stage of the
// posed machinery (frames, chief ray, sampling in global coordinates).
func TestRollInvarianceOfSpotAndOPD(t *testing.T) {
	bk7, _ := Glass("N-BK7")
	const gamma = 0.37
	mk := func(roll float64) *System {
		p := CoordBreak{TiltZ: roll}
		return mustSystem(t, Air{},
			[]Surface{
				{Profile: Sphere{R: 60}, Aperture: CircleAperture(12), Pose: posePtr(p)},
				{Profile: EvenAsphere{R: -60, K: 0, Coeffs: []float64{1e-6}}, Aperture: CircleAperture(12), Pose: posePtr(p)},
			},
			[]float64{4},
			[]Medium{Air{}, bk7, Air{}})
	}
	tx, ty := 0.03, 0.02
	c, s := math.Cos(gamma), math.Sin(gamma)
	fieldR := FieldTan(c*tx-s*ty, s*tx+c*ty)
	pts := Hexapolar(3)
	ptsR := make([]PupilPoint, len(pts))
	for i, p := range pts {
		ptsR[i] = PupilPoint{PX: c*p.PX - s*p.PY, PY: s*p.PX + c*p.PY}
	}
	pup := Pupil{Z: -15, Radius: 7}
	sys0, sys1 := mk(0), mk(gamma)
	par, err := sys0.Paraxial(testLambda, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	imageZ := par.ImageZ
	// Each sampling is centred on its own real chief ray.
	off0, err := sys0.ChiefOffset(FieldTan(tx, ty), pup.Z)
	if err != nil {
		t.Fatal(err)
	}
	off1, err := sys1.ChiefCrossing(fieldR, pup.Z, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	pup0, pup1 := pup, pup
	pup0.Offset, pup1.Offset = off0, off1
	sd0, err := sys0.SpotDiagram(FieldTan(tx, ty), pup0, pts, imageZ, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	sd1, err := sys1.SpotDiagram(fieldR, pup1, ptsR, imageZ, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	st0, st1 := sd0.Stats(ChiefReference()), sd1.Stats(ChiefReference())
	if st1.N != len(pts) {
		t.Fatalf("rolled system: %d of %d rays reached the image plane", st1.N, len(pts))
	}
	// RMS alone is blind to the sign of the roll (the image may rotate as a
	// whole), so compare the image points one by one: the rolled spot must be
	// the original spot rotated by γ about the axis (sign included).
	worstPt := 0.0
	for i := range pts {
		h0, h1 := sd0.Points[i].Hit, sd1.Points[i].Hit
		rx, ry := c*h0.X-s*h0.Y, s*h0.X+c*h0.Y
		worstPt = math.Max(worstPt, math.Hypot(h1.X-rx, h1.Y-ry))
	}
	t.Logf("roll %.2f rad: worst image-point deviation from the rotated original = %.3g mm", gamma, worstPt)
	checkClose(t, "image points rotate with the roll", worstPt, 0, 1e-10)
	t.Logf("roll %.2f rad: |RMS_roll−RMS_0|=%.3g mm", gamma, math.Abs(st1.RMS-st0.RMS))
	checkClose(t, "spot RMS invariant under roll", st1.RMS, st0.RMS, 1e-10)
	g0, err := sys0.OPDGrid(FieldTan(tx, ty), pup0, pts, imageZ, testLambda, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	g1, err := sys1.OPDGrid(fieldR, pup1, ptsR, imageZ, testLambda, OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("roll %.2f rad: |OPD RMS_roll−OPD RMS_0|=%.3g mm", gamma, math.Abs(g1.RMS-g0.RMS))
	checkClose(t, "OPD RMS invariant under roll", g1.RMS, g0.RMS, 1e-10)
}

// TestChiefRayOfPosedSystemHitsStopCentre checks the real chief ray of a
// system with a decentred and tilted stop and a tilted first surface: its trace
// meets the stop surface at the stop's centre (local transverse coordinates
// within 1e-10 mm) and it reaches the image plane.
func TestChiefRayOfPosedSystemHitsStopCentre(t *testing.T) {
	bk7, _ := Glass("N-BK7")
	f2, _ := Glass("F2")
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 80}, Aperture: CircleAperture(12), Pose: posePtr(CoordBreak{TiltX: -0.03})},
			{Profile: Sphere{R: -40}, Aperture: CircleAperture(10), Pose: posePtr(CoordBreak{Dx: 0.5, TiltY: 0.04})},
			{Profile: Sphere{R: -200}, Aperture: CircleAperture(12)},
		},
		[]float64{5, 2},
		[]Medium{Air{}, bk7, f2, Air{}})
	sys.StopIndex = 1
	for _, f := range []Field{FieldAngle(0, 3), FieldAngle(-2, 1), FieldTan(0, 0)} {
		ray := sys.ChiefRay(f, testLambda)
		if ray.Status != StatusOK {
			t.Fatalf("field %v: chief ray status %v", f, ray.Status)
		}
		res := sys.Trace(ray)
		var stop *Step
		for i := range res.Steps {
			if res.Steps[i].Surface == sys.StopIndex {
				stop = &res.Steps[i]
			}
		}
		if stop == nil {
			t.Fatalf("field %v: chief ray does not reach the stop", f)
		}
		h := sys.elementFrame(sys.StopIndex).toLocal(stop.Point)
		t.Logf("field %v: chief miss of the stop centre %.3g mm", f, math.Hypot(h.X, h.Y))
		if math.Hypot(h.X, h.Y) > 1e-10 {
			t.Errorf("field %v: chief ray misses the stop centre by (%.3g, %.3g) mm", f, h.X, h.Y)
		}
		if !res.OK() {
			t.Errorf("field %v: chief ray trace status %v", f, res.Status)
		}
	}
}

// TestCoordBreakAndReturnIsIdentity is acceptance (c): a coordinate break
// followed by its return break (Inverse) at zero spacing, placed inside a
// glass gap of a singlet, leaves every traced ray unchanged. The composite
// system has 4 elements (the two breaks are not surfaces, so the physical steps
// are the two singlet surfaces). Random rays are compared with the plain
// singlet: positions and OPL agree to 1e-12 mm, directions to 1e-12. A negative
// control shows that the break alone is NOT the identity.
func TestCoordBreakAndReturnIsIdentity(t *testing.T) {
	g := ConstantIndex(1.5)
	cb := CoordBreak{Dx: 0.4, Dy: -0.25, TiltX: 0.3, TiltY: -0.2, TiltZ: 0.15}
	surfaces := func() []Surface {
		return []Surface{
			{Profile: Sphere{R: 60}, Aperture: CircleAperture(15)},
			{Profile: Sphere{R: -60}, Aperture: CircleAperture(15)},
		}
	}
	plain := mustSystem(t, Air{}, surfaces(), []float64{4}, []Medium{Air{}, g, Air{}})

	const a = 1.7 // S0 -> break, the rest of the 4 mm gap follows the return
	withPair := mustSystem(t, Air{},
		[]Surface{
			surfaces()[0],
			{Break: true, Pose: posePtr(cb)},
			{Break: true, Pose: posePtr(cb.Inverse())},
			surfaces()[1],
		},
		[]float64{a, 0, 4 - a},
		[]Medium{Air{}, g, g, g, Air{}})

	withBreakOnly := mustSystem(t, Air{},
		[]Surface{
			surfaces()[0],
			{Break: true, Pose: posePtr(cb)},
			surfaces()[1],
		},
		[]float64{a, 4 - a},
		[]Medium{Air{}, g, g, Air{}})

	rng := rand.New(rand.NewSource(20261015))
	worstPos, worstDir, worstOPL, worstBreakOnly := 0.0, 0.0, 0.0, 0.0
	for k := 0; k < 200; k++ {
		origin := V(rng.Float64()*4-2, rng.Float64()*4-2, -20)
		dir := V(rng.Float64()*0.2-0.1, rng.Float64()*0.2-0.1, 1).Normalized()
		r := NewRay(origin, dir, testLambda)
		ra := plain.Trace(r)
		rb := withPair.Trace(r)
		if !ra.OK() || !rb.OK() || len(ra.Steps) != 2 || len(rb.Steps) != 2 {
			t.Fatalf("ray %d: traces must reach both surfaces (plain %v/%d, pair %v/%d)",
				k, ra.Status, len(ra.Steps), rb.Status, len(rb.Steps))
		}
		for i := range ra.Steps {
			worstPos = math.Max(worstPos, ra.Steps[i].Point.Sub(rb.Steps[i].Point).Norm())
			worstDir = math.Max(worstDir, ra.Steps[i].DirOut.Sub(rb.Steps[i].DirOut).Norm())
		}
		worstPos = math.Max(worstPos, ra.Ray.Origin.Sub(rb.Ray.Origin).Norm())
		worstDir = math.Max(worstDir, ra.Ray.Dir.Sub(rb.Ray.Dir).Norm())
		worstOPL = math.Max(worstOPL, math.Abs(ra.Ray.OPL-rb.Ray.OPL))
		if rc := withBreakOnly.Trace(r); rc.OK() {
			worstBreakOnly = math.Max(worstBreakOnly, rc.Ray.Origin.Sub(ra.Ray.Origin).Norm())
		}
	}
	checkClose(t, "break + return: worst position error (mm)", worstPos, 0, 1e-12)
	checkClose(t, "break + return: worst direction error", worstDir, 0, 1e-12)
	checkClose(t, "break + return: worst OPL error (mm)", worstOPL, 0, 1e-12)
	if worstBreakOnly < 1e-3 {
		t.Errorf("negative control: the break alone changed the trace by only %.3g mm", worstBreakOnly)
	}
	t.Logf("break+return vs plain: worst position %.3g mm, direction %.3g, OPL %.3g mm; break alone %.3g mm",
		worstPos, worstDir, worstOPL, worstBreakOnly)
}

// TestPosedSystemRejectsFirstOrderEngines is acceptance (e): every first-order
// and third-order routine, and the sub-system and reversal helpers, return an
// error (or nil) for a system with a decentred or tilted surface, instead of a
// silently wrong number.
func TestPosedSystemRejectsFirstOrderEngines(t *testing.T) {
	sys := decentredLens(t, 0.2, 0.01)
	inf := math.Inf(-1)
	field := FieldAngle(0, 1)
	checks := map[string]error{}
	_, checks["Paraxial"] = sys.Paraxial(testLambda, inf)
	_, checks["EntrancePupil"] = sys.EntrancePupil(field, testLambda)
	_, checks["ExitPupil"] = sys.ExitPupil(field, testLambda)
	_, checks["ParaxialImageHeight"] = sys.ParaxialImageHeight(field, 100, testLambda)
	_, checks["Distortion"] = sys.Distortion(field, 100, testLambda)
	_, checks["SeidelCoefficients"] = sys.SeidelCoefficients(testLambda, field, 5)
	_, checks["FocusAt"] = sys.FocusAt(-200, testLambda)
	_, checks["AxialChromatic"] = sys.AxialChromatic(-200, WavelengthF, WavelengthC)
	_, checks["LateralChromatic"] = sys.LateralChromatic(field, 100, WavelengthF, WavelengthC)
	_, checks["ChiefOffset"] = sys.ChiefOffset(field, -20)
	_, checks["Sub"] = sys.Sub(0, 1)
	for name, err := range checks {
		if err == nil {
			t.Errorf("%s on a decentred/tilted system: want an error, got nil", name)
		}
	}
	if sys.Reversed() != nil {
		t.Errorf("Reversed of a posed system: want nil")
	}
}

// TestRollAndAxialBreaksKeepFirstOrderOptics checks the cases that the
// first-order engine does accept: a roll about z and an axial break with zero
// decentre and tilt leave the axis unchanged, so the paraxial numbers equal the
// plain system's to 1e-12.
func TestRollAndAxialBreaksKeepFirstOrderOptics(t *testing.T) {
	plain := nominalLens(t)
	rolled := decentredLens(t, 0, 0)
	rolled.Surfaces[0].Pose = posePtr(CoordBreak{TiltZ: 0.4})
	rolled.Surfaces[1].Pose = posePtr(CoordBreak{TiltZ: 0.4})
	withAxialBreak := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 100}, Aperture: CircleAperture(15)},
			{Break: true, Pose: posePtr(CoordBreak{})},
			{Profile: Sphere{R: -100}, Aperture: CircleAperture(15)},
		},
		[]float64{0, 0},
		[]Medium{Air{}, ConstantIndex(1.5), ConstantIndex(1.5), Air{}})
	p0, err := plain.Paraxial(testLambda, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	for name, sys := range map[string]*System{"roll": rolled, "axial break": withAxialBreak} {
		p, err := sys.Paraxial(testLambda, math.Inf(-1))
		if err != nil {
			t.Fatalf("%s: paraxial must be accepted: %v", name, err)
		}
		checkClose(t, name+": EFL", p.EFL, p0.EFL, 1e-12)
		checkClose(t, name+": BFL", p.BFL, p0.BFL, 1e-12)
	}
}

// TestBreakAndPoseValidation checks the new validation rules.
func TestBreakAndPoseValidation(t *testing.T) {
	g := ConstantIndex(1.5)
	mk := func(surf []Surface, th []float64, med []Medium, stop int) error {
		s := &System{Surfaces: surf, Thickness: th, Mediums: med, StopIndex: stop}
		return s.Validate()
	}
	plane := Surface{Profile: Plane{}}
	if err := mk([]Surface{plane, {Break: true, Pose: posePtr(CoordBreak{Dx: 1})}, plane},
		[]float64{1, 1}, []Medium{Air{}, g, ConstantIndex(1.7), Air{}}, 0); err == nil {
		t.Error("a break between different media: want an error")
	}
	if err := mk([]Surface{plane, {Break: true, Mirror: true}, plane},
		[]float64{1, 1}, []Medium{Air{}, g, g, Air{}}, 0); err == nil {
		t.Error("a mirror break: want an error")
	}
	if err := mk([]Surface{plane, {Break: true, Pose: posePtr(CoordBreak{Dx: 1})}, plane},
		[]float64{1, 1}, []Medium{Air{}, g, g, Air{}}, 1); err == nil {
		t.Error("a break as the stop: want an error")
	}
	if err := mk([]Surface{{Profile: Plane{}, Pose: posePtr(CoordBreak{Dx: math.NaN()})}},
		nil, []Medium{Air{}, Air{}}, 0); err == nil {
		t.Error("a NaN pose: want an error")
	}
	if err := mk([]Surface{{Break: true, Pose: posePtr(CoordBreak{Dx: 1})}, plane},
		[]float64{1}, []Medium{Air{}, Air{}, Air{}}, 1); err != nil {
		t.Errorf("a valid break before a surface rejected: %v", err)
	}
}

// TestWithImagePlaneLiesOnGlobalZ checks that the image plane appended to a
// posed system lies at the requested global z (its frame origin).
func TestWithImagePlaneLiesOnGlobalZ(t *testing.T) {
	sys := decentredLens(t, 0.3, 0.05)
	withPlane := sys.WithImagePlane(87.25)
	e := withPlane.elementFrame(len(withPlane.Surfaces) - 1)
	checkClose(t, "image plane global z", e.o.Z, 87.25, 1e-12)
}

// TestAutoReturnAndPoseFrame checks the frame of a posed surface and the
// auto-return of the chain: the decentre is applied in the chain frame, the
// tilt about the decentred vertex, and the next surface stays on the axis at
// the nominal thickness.
func TestAutoReturnAndPoseFrame(t *testing.T) {
	const th = 0.2
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(20), Pose: posePtr(CoordBreak{Dx: 0.4, TiltX: th})},
			{Profile: Plane{}, Aperture: CircleAperture(20)},
		},
		[]float64{3},
		[]Medium{Air{}, ConstantIndex(1.5), Air{}})
	e0 := sys.elementFrame(0)
	checkClose(t, "posed vertex x (decentre)", e0.o.X, 0.4, 0)
	checkClose(t, "posed vertex y", e0.o.Y, 0, 0)
	checkClose(t, "posed vertex z", e0.o.Z, 0, 0)
	ax := e0.r.col(2)
	checkClose(t, "posed axis y = -sin(tilt)", ax.Y, -math.Sin(th), 1e-15)
	checkClose(t, "posed axis z = cos(tilt)", ax.Z, math.Cos(th), 1e-15)
	e1 := sys.elementFrame(1)
	if e1.o != V(0, 0, 3) || e1.r != ident3() {
		t.Errorf("following surface not restored to the axis: origin %v, axes %v", e1.o, e1.r)
	}
	if z := sys.VertexZ(1); z != 3 {
		t.Errorf("VertexZ(1) = %v, want 3 (thickness measured on the chain)", z)
	}
}

// TestChiefRayUnlimitedApertures checks the chief-ray search of a posed system
// whose surfaces have no aperture: the launch plane must stay finite and the
// chief ray must still hit the stop centre.
func TestChiefRayUnlimitedApertures(t *testing.T) {
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 100}, Pose: posePtr(CoordBreak{Dx: 0.3, TiltY: 0.02})},
			{Profile: Sphere{R: -100}, Pose: posePtr(CoordBreak{Dx: 0.3, TiltY: 0.02})},
		},
		[]float64{0},
		[]Medium{Air{}, ConstantIndex(1.5), Air{}})
	if z := sys.launchPlaneZ(); math.IsInf(z, 0) || math.IsNaN(z) {
		t.Fatalf("launch plane not finite: %v", z)
	}
	ray := sys.ChiefRay(FieldAngle(0, 1), testLambda)
	if ray.Status != StatusOK {
		t.Fatalf("chief ray status %v", ray.Status)
	}
	res := sys.Trace(ray)
	h := sys.elementFrame(0).toLocal(res.Steps[0].Point)
	if math.Hypot(h.X, h.Y) > 1e-10 {
		t.Errorf("chief ray misses the stop centre by (%.3g, %.3g)", h.X, h.Y)
	}
}

// TestAnalysisReferenceModesOnPosedSystem checks the three transverse
// references on a decentred singlet (all surfaces moved by (d, e), sampling
// centred on the real chief): chief = the real chief hit, centroid = the mean of
// the valid hits, vertex = the global axis point (0, 0, imageZ) (the image of
// the nominal axis; it is not the image of the decentred system's own axis).
// The spot RMS about any point P obeys the parallel-axis identity
// RMS(P)² = RMS(centroid)² + |P − centroid|², checked to 1e-12 mm.
func TestAnalysisReferenceModesOnPosedSystem(t *testing.T) {
	bk7, _ := Glass("N-BK7")
	const d, e = 0.8, -0.45
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 50}, Aperture: CircleAperture(12.5), Pose: posePtr(CoordBreak{Dx: d, Dy: e})},
			{Profile: Sphere{R: -50}, Aperture: CircleAperture(12.5), Pose: posePtr(CoordBreak{Dx: d, Dy: e})},
		},
		[]float64{5},
		[]Medium{Air{}, bk7, Air{}})
	field := FieldAngle(0, 2)
	pup := Pupil{Z: -20, Radius: 8}
	off, err := sys.ChiefCrossing(field, pup.Z, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	pup.Offset = off
	imageZ := 100.0 // any plane near focus; the references are defined on it
	pts := Hexapolar(4)
	sd, err := sys.SpotDiagram(field, pup, pts, imageZ, testLambda)
	if err != nil {
		t.Fatal(err)
	}
	chief := sd.Stats(ChiefReference())
	cen := sd.Stats(CentroidReference())
	vert := sd.Stats(VertexReference())
	if chief.Reference != sd.Chief {
		t.Errorf("chief reference %v, want the chief hit %v", chief.Reference, sd.Chief)
	}
	if cen.Reference.Sub(sd.Centroid).Norm() > 1e-12 {
		t.Errorf("centroid reference %v, want %v", cen.Reference, sd.Centroid)
	}
	if vert.Reference != V(0, 0, imageZ) {
		t.Errorf("vertex reference %v, want (0, 0, %g)", vert.Reference, imageZ)
	}
	for name, st := range map[string]SpotStats{"chief": chief, "vertex": vert} {
		lhs := st.RMS * st.RMS
		rhs := cen.RMS*cen.RMS + st.Reference.Sub(cen.Centroid).Norm2()
		if math.Abs(lhs-rhs) > 1e-12 {
			t.Errorf("parallel-axis identity violated for the %s reference: %.15g vs %.15g", name, lhs, rhs)
		}
	}
	// The OPD reference point follows the same definitions.
	for _, c := range []struct {
		ref  Reference
		want Vec3
	}{
		{ChiefReference(), sd.Chief},
		{CentroidReference(), sd.Centroid},
		{VertexReference(), V(0, 0, imageZ)},
	} {
		g, err := sys.OPDGrid(field, pup, pts, imageZ, testLambda, OPDOptions{Reference: c.ref})
		if err != nil {
			t.Fatalf("OPD with %v reference: %v", c.ref, err)
		}
		if g.ReferencePoint.Sub(c.want).Norm() > 1e-10 {
			t.Errorf("OPD %v reference %v, want %v", c.ref, g.ReferencePoint, c.want)
		}
	}
	t.Logf("references: chief %v, centroid %v, vertex %v; RMS chief %.6g, centroid %.6g, vertex %.6g mm",
		sd.Chief, sd.Centroid, vert.Reference, chief.RMS, cen.RMS, vert.RMS)
}
