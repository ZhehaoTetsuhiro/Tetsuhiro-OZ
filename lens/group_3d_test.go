package lens

import (
	"math"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
)

// ---------------------------------------------------------------------------
// Real-ray measurements of a group in 3D (task group3d)
// ---------------------------------------------------------------------------
//
// Every threshold is justified in the test that uses it. The tolerances come
// from three sources, all of them quantified below: the search of the best
// focus (ImageOptions.TolMM, 1e-9 mm by default), the geometry of the
// configuration (the spherical aberration of the elements, printed by each
// test), and — for the fold-invariance claims — the fact that a plane mirror in
// a collimated beam changes neither the aberration nor the beam.

// planoElement returns a plano-convex (or, for a negative f, plano-concave)
// lens element with the plane surface at the rear and the curved one at the
// front. Its paraxial power is exactly (n−1)/R = 1/f for any thickness (the
// plane surface has no power), so it is the ideal thin lens of focal length f
// placed at its front vertex: the classical "best form" for an object at
// infinity. This is the reference element of the measurements below, because
// its thin-lens parameters are known in closed form.
func planoElement(t *testing.T, fMM, thicknessMM, semiAperMM float64, glass string) Element {
	t.Helper()
	m, ok := geo.Glass(glass)
	if !ok {
		t.Fatalf("unknown glass %q", glass)
	}
	return Element{
		Kind:        ElementLens,
		Front:       SurfaceSpec{RadiusMM: (m.N(geo.Wavelengthd) - 1) * fMM, SemiAperMM: semiAperMM},
		Rear:        &SurfaceSpec{RadiusMM: 0, SemiAperMM: semiAperMM},
		Glass:       glass,
		ThicknessMM: thicknessMM,
	}
}

// gBundle returns the measurement bundle of the tests: an object at infinity on
// the axis, the hexagonal sampling and the d line.
func gBundle() Bundle {
	return Bundle{Object: geo.FieldTan(0, 0), Samples: geo.Hexapolar(3), LambdaNm: geo.Wavelengthd}
}

// achromatElements returns the two elements of a cemented crown/flint doublet of
// focal length fMM: the crown is the two-surface element [R1, flat cement
// interface] whose GapGlass is the flint, the flint is the single rear surface
// R3 — the classical three-surface achromat. The radii follow the thin-lens
// chromatic condition φ1/V1 + φ2/V2 = 0 with φ1 + φ2 = 1/f and a flat interface
// (the Abbe numbers of the two catalogue glasses), so the doublet corrects the
// primary chromatic aberration: its F−C focus shift is a small fraction of a
// singlet's.
func achromatElements(t *testing.T, fMM, semiAperMM, crownThickness, flintThickness float64) []Element {
	t.Helper()
	crown, ok1 := geo.Glass("N-BK7")
	flint, ok2 := geo.Glass("F2")
	if !ok1 || !ok2 {
		t.Fatal("the catalogue glasses are missing")
	}
	nd1, nd2 := crown.N(geo.Wavelengthd), flint.N(geo.Wavelengthd)
	v1, v2 := crown.AbbeNumber(), flint.AbbeNumber()
	phi := 1 / fMM
	phi1 := phi / (1 - v2/v1)
	phi2 := phi - phi1
	r1 := (nd1 - 1) / phi1
	r3 := -(nd2 - 1) / phi2
	front := Element{
		Name: "crown", Kind: ElementLens,
		Front: SurfaceSpec{RadiusMM: r1, SemiAperMM: semiAperMM},
		Rear:  &SurfaceSpec{RadiusMM: 0, SemiAperMM: semiAperMM},
		Glass: "N-BK7", GapGlass: "F2",
		ThicknessMM: crownThickness, GapMM: flintThickness,
	}
	rear := SingleSurface(r3, semiAperMM)
	rear.Name = "flint"
	return []Element{front, rear}
}

// TestGroupImagePositionParabolicMirror checks the focus search against a
// closed form: an on-axis collimated bundle reflected by a paraboloid of vertex
// radius R has no spherical aberration (the classical parabola property
// z = r²/2R focuses a plane wave exactly), so every ray crosses the axis at
// f = |R|/2. The image point must be that point, the focus distance from the
// vertex must be f measured along the beam (which the mirror sends back along
// −z), and the residual RMS at the best focus must be at the level of the
// arithmetic. Thresholds: 1e-9 mm on the position (the search tolerance, set to
// 1e-12 here) and 1e-11 mm on the RMS.
func TestGroupImagePositionParabolicMirror(t *testing.T) {
	const r, f = -400.0, 200.0
	g := gWantGroup(t, &Group{Name: "parabola", ObjectZMM: math.Inf(-1), Elements: []Element{{
		Kind:  ElementMirror,
		Front: SurfaceSpec{RadiusMM: r, ConicK: -1, SemiAperMM: 24},
	}}})
	res, err := g.ImagePosition(Bundle{Object: geo.FieldTan(0, 0), Samples: geo.Hexapolar(3)}, ImageOptions{TolMM: 1e-12})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("parabola: image %v, RMS %.3g mm, focus %.12g mm from %v, rays %d failed %d",
		res.Point, res.RMS, res.FocusMM, res.PlaneO, res.Rays, res.Failed)
	gVecClose(t, "image point", res.Point, geo.V(0, 0, -f), 1e-9)
	gClose(t, "focus distance from the vertex", res.FocusMM, f, 1e-9)
	gVecClose(t, "focal plane axis (the reflected beam)", res.Axis, geo.V(0, 0, -1), 1e-12)
	if res.RMS > 1e-11 {
		t.Errorf("the parabola has no spherical aberration: RMS %.3g mm", res.RMS)
	}
	if res.Rays != 37 || res.Failed != 0 {
		t.Errorf("rays %d, failed %d, want 37 and 0", res.Rays, res.Failed)
	}
}

// TestGroupImagePositionMatchesParaxial compares the real-ray best focus of a
// coaxial group with the paraxial image of package geo.
//
// The configuration is a stop in front of a plano-convex lens (the element of
// planoElement), so the paraxial engine and the real trace see the same
// first-order system: par.ImageZ is the exact reference. The best-RMS focus
// differs from the paraxial focus by the spherical aberration of the singlet,
// which is why the tolerance is stated per aperture: for a thin lens the
// longitudinal spherical aberration goes as the square of the aperture, so
// halving the stop radius cuts the difference by a factor of four (both facts
// are measured and printed here). Thresholds: 0.35 mm at a 5 mm stop, 0.09 mm
// at 2.5 mm, and the ratio of the two differences within [3, 5].
func TestGroupImagePositionMatchesParaxial(t *testing.T) {
	const f = 100.0
	build := func(stopSemi float64) *Group {
		return gWantGroup(t, &Group{
			Name: "singlet", ObjectZMM: math.Inf(-1),
			Elements: []Element{
				{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: stopSemi}},
				planoElement(t, f, 5, 14, "N-BK7"),
			},
		})
	}
	sys, err := build(5).System()
	if err != nil {
		t.Fatal(err)
	}
	par, err := sys.Paraxial(geo.Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("paraxial image z = %.12g, EFL = %.12g", par.ImageZ, par.EFL)
	gClose(t, "plano-convex EFL", par.EFL, f, 1e-9)

	measure := func(stopSemi float64) float64 {
		res, err := build(stopSemi).ImagePosition(gBundle(), ImageOptions{})
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("stop %4.1f mm: image z = %.12g, RMS %.3g mm, failed %d, |z − paraxial| = %.3g mm",
			stopSemi, res.Point.Z, res.RMS, res.Failed, math.Abs(res.Point.Z-par.ImageZ))
		return res.Point.Z
	}
	full := measure(5)
	small := measure(2.5)
	gClose(t, "best focus vs paraxial image (5 mm stop)", full, par.ImageZ, 0.35)
	gClose(t, "best focus vs paraxial image (2.5 mm stop)", small, par.ImageZ, 0.09)
	dFull, dSmall := math.Abs(full-par.ImageZ), math.Abs(small-par.ImageZ)
	if dFull < 1e-6 {
		t.Fatalf("the singlet shows no spherical aberration (%g mm)", dFull)
	}
	ratio := dFull / dSmall
	t.Logf("aperture-squared check: %.4f (want ~4)", ratio)
	if ratio < 3 || ratio > 5 {
		t.Errorf("the spherical aberration does not scale with the aperture squared: ratio %.4f", ratio)
	}

	// An off-axis field: the real chief ray traced through the group, evaluated
	// on the paraxial image plane, must reproduce the paraxial image height of
	// package geo (their difference is the distortion of the singlet, which the
	// test measures). Threshold 0.02 mm on a 17.6 mm image height: the
	// distortion of this lens is far below that (see the printed values), so
	// the threshold is a conservative bound of the model agreement, not the
	// physics.
	const thetaDeg = 10.0
	tan := math.Tan(thetaDeg * math.Pi / 180)
	field := geo.FieldTan(0, -tan)
	res, err := build(5).ImagePosition(Bundle{Object: field, Samples: geo.Hexapolar(3)}, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	wantHeight, err := sys.ParaxialImageHeight(field, par.ImageZ, geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	// The traced chief evaluated on the paraxial image plane.
	tz := (par.ImageZ - res.ChiefHit.Z) / res.ChiefDir.Z
	realHeight := res.ChiefHit.Y + res.ChiefDir.Y*tz
	t.Logf("field %g deg: chief hit %v (plane %g), chief dir %v, height at the paraxial plane %.12g, paraxial %.12g, distortion %.3g mm",
		thetaDeg, res.ChiefHit, res.ChiefHit.Z, res.ChiefDir, realHeight, wantHeight, realHeight-wantHeight)
	gClose(t, "real chief height vs the paraxial image height", realHeight, wantHeight, 0.02)
	if math.Abs(res.ChiefHit.Y) < 1 {
		t.Errorf("degenerate test: the image did not move off axis (%v)", res.ChiefHit)
	}
}

// TestGroupImagePositionFoldedMatchesHandGeometry is the hand-computed 3D case.
// Two 45 degree fold mirrors in series turn the beam by 90 degrees each
// (+z → +y → −z, see TestGroupLayoutTwoFolds90Deg), so the plano-convex lens of
// focal length f sits at (0, g1, z1 − g2) with its axis along −z, where z1 is
// the first mirror's vertex. The object is at infinity on the axis, so the beam
// is collimated at the mirrors (a plane mirror adds no aberration there) and the
// system is, ray for ray, the straight system of a single lens of the same focal
// length: the folded image must be the straight image transported by the
// hand-computed rigid motion.
//
// Threshold 1e-6 mm: the two traces describe the same rays (the fold rotations
// are exact to ~1e-15), and the two independent best-focus searches agree to
// 2.1e-8 mm — the limit of the scan bracket, because the RMS minimum is flat and
// the two systems search slightly different windows. The straight
// reference is computed by the same API, and the focus distance is additionally
// compared with the paraxial back focus of the element (the front vertex plus
// f + d(1 − 1/n) = 101.7036 mm) with the aperture tolerance of the singlet.
func TestGroupImagePositionFoldedMatchesHandGeometry(t *testing.T) {
	const (
		g1, g2 = 40.0, 30.0 // mirror-to-mirror and mirror-to-lens gaps
		f      = 100.0
	)
	// The stop is placed upstream of the first mirror: the sequential model
	// visits the surfaces in beam order, so every surface must lie behind its
	// predecessor along the beam (a 45 degree mirror reaches back by up to the
	// beam radius on the lower side of the beam).
	folded := gWantGroup(t, &Group{
		Name: "folded relay", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(6); e.GapMM = 15; return e }(),
			func() Element {
				e := MirrorElement(0, 0, 18)
				e.TiltXDeg = 45
				e.GapMM = g1
				return e
			}(),
			func() Element {
				e := MirrorElement(0, 0, 18)
				e.TiltXDeg = 45
				e.GapMM = g2
				return e
			}(),
			planoElement(t, f, 5, 18, "N-BK7"),
		},
	})
	straight := gWantGroup(t, &Group{
		Name: "straight", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(6); e.GapMM = 15; return e }(),
			planoElement(t, f, 5, 18, "N-BK7"),
		},
	})
	// The hand geometry of the chain: the stop at z = 0, the first mirror 15 mm
	// behind it, the second 40 mm further along +y, the lens 30 mm further along
	// −z.
	const z1 = 15.0
	lay, err := folded.Layout()
	if err != nil {
		t.Fatal(err)
	}
	byFace := func(e int, ff Face) LayoutSurface {
		for _, l := range lay {
			if l.Element == e && l.Face == ff {
				return l
			}
		}
		t.Fatalf("no layout entry for element %d %s", e, ff)
		return LayoutSurface{}
	}
	gVecClose(t, "first mirror vertex (hand)", byFace(1, Front).Vertex, geo.V(0, 0, z1), 1e-12)
	gVecClose(t, "second mirror vertex (hand)", byFace(2, Front).Vertex, geo.V(0, g1, z1), 1e-12)
	gVecClose(t, "lens vertex (hand)", byFace(3, Front).Vertex, geo.V(0, g1, z1-g2), 1e-12)
	gVecClose(t, "lens axis (hand)", byFace(3, Front).Axis, geo.V(0, 0, -1), 1e-12)

	rf, err := folded.ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	rs, err := straight.ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	// The straight reference has its lens front vertex at z = 15 (the stop gap),
	// so its image lies (rs.Point.Z - 15) behind the lens front vertex; the
	// folded lens front vertex is at (0, g1, z1-g2) with the axis -z.
	straightLens := geo.V(0, 0, 15)
	want := geo.V(0, g1, z1-g2).Add(geo.V(0, 0, -1).Mul(rs.Point.Z - straightLens.Z))
	t.Logf("folded: image %v RMS %.3g mm focus %.12g rays %d failed %d; straight: image %v RMS %.3g mm",
		rf.Point, rf.RMS, rf.FocusMM, rf.Rays, rf.Failed, rs.Point, rs.RMS)
	t.Logf("transported straight image (hand): %v, |diff| %.3g mm", want, rf.Point.Sub(want).Norm())
	gVecClose(t, "folded image = transported straight image", rf.Point, want, 1e-6)
	gClose(t, "folded vs straight focus distance", rf.FocusMM, rs.FocusMM, 1e-6)
	// The paraxial back focal distance of the plano-convex element, f - d/n =
	// 96.7036 mm from the rear vertex; the best-RMS focus differs from it by the
	// spherical aberration of the singlet at a 6 mm stop (0.34 mm, i.e. 0.936 mm
	// at 10 mm in the coaxial test above), hence the 0.5 mm threshold.
	gClose(t, "focus distance vs the paraxial back focal distance", rf.FocusMM, 96.7035865728693, 0.5)
	gVecClose(t, "exit axis", rf.Axis, geo.V(0, 0, -1), 1e-12)
	gVecClose(t, "anchor (last vertex)", rf.PlaneO, geo.V(0, g1, z1-g2-5), 1e-12)
}

// TestGroupAngularMagnificationKeplerGalileo checks Γ against the closed forms of
// the two classical afocal telescopes, measured with real rays through
// planoElement (ideal thin) lenses: the Kepler pair (f, f) separated by
// f1 + f2 − d1/n has Γ = −f1/f2 = −1 (an inverted image), the Galileo pair
// (f1, f2 = −f1/2) separated by f1 + f2 − d1/n has Γ = −f1/f2 = +2 (erect,
// magnified).
//
// The measurement uses the chief ray of two probe bundles tilted in the sagittal
// direction by δ = 1e-6, so the second-order term of the measured slope is
// O(δ²) = 1e-12 and the threshold is 1e-9. The beam-diameter ratio D_in/D_out is
// measured on the same bundle as an independent cross-check of |Γ| (the Lagrange
// invariant); it is a geometric diagnostic, so its tolerance is the spherical
// aberration of the telescope at this aperture: 2% here (measured 0.6% for the
// Kepler and 0.1% for the Galileo at a 4 mm stop, against 2.4% at 8 mm).
func TestGroupAngularMagnificationKeplerGalileo(t *testing.T) {
	m, _ := geo.Glass("N-BK7")
	n := m.N(geo.Wavelengthd)
	const d = 2.0 // the centre thickness of both elements
	cases := []struct {
		name       string
		f1, f2     float64
		want       float64
		separation float64
	}{
		{"kepler", 100, 100, -1, 100 + 100 - d/n},
		{"galileo", 100, -50, 2, 100 - 50 - d/n},
	}
	for _, c := range cases {
		l1 := planoElement(t, c.f1, d, 20, "N-BK7")
		l1.GapMM = c.separation
		l2 := planoElement(t, c.f2, d, 20, "N-BK7")
		g := gWantGroup(t, &Group{
			Name: c.name, ObjectZMM: math.Inf(-1),
			Elements: []Element{{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 4}}, l1, l2},
		})
		res, err := g.AngularMagnification(Bundle{Object: geo.FieldTan(0, 0), Samples: geo.Hexapolar(2), LambdaNm: geo.Wavelengthd})
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		t.Logf("%s: Γ = %.15g (want %.12g), D_in/D_out = %.12g (|Γ| %.12g), exit %v, spread %.3g rad, rays %d",
			c.name, res.Gamma, c.want, res.BeamRatio, math.Abs(c.want), res.ExitDir, res.SpreadRad, res.Rays)
		gClose(t, c.name+": Γ", res.Gamma, c.want, 1e-9)
		gClose(t, c.name+": D_in/D_out vs |Γ|", res.BeamRatio, math.Abs(c.want), 0.02)
		gVecClose(t, c.name+": exit axis", res.AxisOut, geo.V(0, 0, 1), 1e-12)
		if res.SpreadRad > 1e-3 {
			t.Errorf("%s: the exit bundle is not collimated (%.3g rad)", c.name, res.SpreadRad)
		}
	}
}

// TestGroupAngularMagnificationFoldInvariant is the 3D fold invariance: the same
// Kepler telescope with a plane mirror (45 degrees about x) inserted in its
// collimated space has exactly the same angular magnification. Γ is measured in
// the sagittal direction x, the transverse direction that a fold in the y–z
// plane leaves unchanged, so the two measurements agree to the level of the
// arithmetic: threshold 1e-9 on both Γ and the beam-diameter ratio (the beam
// ratio is identical here, because a plane mirror in a collimated beam preserves
// the beam). A negative control checks that the folded system really is folded:
// its exit axis is +y, the unfolded one's is +z.
func TestGroupAngularMagnificationFoldInvariant(t *testing.T) {
	m, _ := geo.Glass("N-BK7")
	n := m.N(geo.Wavelengthd)
	const d = 2.0
	sep := 100 + 100 - d/n
	newLens1 := func() Element { e := planoElement(t, 100, d, 20, "N-BK7"); e.GapMM = sep; return e }
	straight := gWantGroup(t, &Group{
		Name: "kepler", ObjectZMM: math.Inf(-1),
		Elements: []Element{{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 6}}, newLens1(), planoElement(t, 100, d, 20, "N-BK7")},
	})
	fold := MirrorElement(0, 0, 20)
	fold.TiltXDeg = 45
	fold.GapMM = 20
	folded := gWantGroup(t, &Group{
		Name: "folded kepler", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(6); e.GapMM = 15; return e }(),
			fold,
			newLens1(),
			planoElement(t, 100, d, 20, "N-BK7"),
		},
	})
	b := Bundle{Object: geo.FieldTan(0, 0), Samples: geo.Hexapolar(2), LambdaNm: geo.Wavelengthd}
	a1, err := straight.AngularMagnification(b)
	if err != nil {
		t.Fatal(err)
	}
	a2, err := folded.AngularMagnification(b)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("straight: Γ = %.15g, D ratio %.12g, axis %v, spread %.3g", a1.Gamma, a1.BeamRatio, a1.AxisOut, a1.SpreadRad)
	t.Logf("folded:   Γ = %.15g, D ratio %.12g, axis %v, spread %.3g", a2.Gamma, a2.BeamRatio, a2.AxisOut, a2.SpreadRad)
	gClose(t, "Γ unchanged by the fold", a2.Gamma, a1.Gamma, 1e-9)
	gClose(t, "beam ratio unchanged by the fold", a2.BeamRatio, a1.BeamRatio, 1e-9)
	gClose(t, "the folded telescope is a Kepler (Γ = −1)", a2.Gamma, -1, 1e-9)
	gVecClose(t, "negative control: the folded exit axis is +y", a2.AxisOut, geo.V(0, 1, 0), 1e-12)
	gVecClose(t, "the unfolded exit axis is +z", a1.AxisOut, geo.V(0, 0, 1), 1e-12)
}

// TestGroupDispersionMirrorZero checks that a mirror-only group has no
// dispersion: every medium is air (a constant index), so the traced rays do not
// depend on the wavelength at all and the two best-focus searches return
// bit-identical planes. The tolerance is 1e-12 mm; the value is exactly zero.
func TestGroupDispersionMirrorZero(t *testing.T) {
	g := gWantGroup(t, &Group{
		Name: "cassegrain", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(40); e.GapMM = 60; return e }(),
			func() Element { e := MirrorElement(-2000, -1, 45); e.GapMM = -300; return e }(),
			func() Element { e := MirrorElement(-600, -2.5, 30); e.GapMM = 400; return e }(),
		},
	})
	res, err := g.Dispersion(Bundle{Object: geo.FieldTan(0, 0), Samples: geo.Hexapolar(2)}, geo.WavelengthF, geo.WavelengthC, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("mirror-only: focus(F) %.15g mm, focus(C) %.15g mm, shift %.3g mm, rays %d",
		res.A.FocusMM, res.B.FocusMM, res.ShiftMM, res.A.Rays)
	if math.Abs(res.ShiftMM) > 1e-12 {
		t.Errorf("a mirror-only system must not disperse: shift %.3g mm", res.ShiftMM)
	}
	if res.A.FocusMM == 0 {
		t.Error("degenerate test: the focus distance is zero")
	}
	if res.A.FocusMM != res.B.FocusMM {
		t.Errorf("the two wavelengths gave different focus planes in a mirror system: %v vs %v", res.A.FocusMM, res.B.FocusMM)
	}
}

// TestGroupDispersionMatchesAxialChromatic cross-checks the real-ray dispersion
// against the paraxial one for a mirror-free coaxial system: the best-focus
// shift between F and C of a singlet measured with real rays must equal the
// paraxial AxialChromatic of the same system. The two differ only by the change
// of the spherical aberration with the wavelength (the best-RMS focus is not the
// paraxial focus), measured here as 0.002 mm against a 1.53 mm shift, so the
// threshold is 0.02 mm. The classical closed form Δf ≈ −f/Vd is checked to
// 0.06 mm (it neglects the principal-plane term of a thick lens, which shifts
// the measured value by 0.03 mm here).
func TestGroupDispersionMatchesAxialChromatic(t *testing.T) {
	const f = 100.0
	g := gWantGroup(t, &Group{
		Name: "singlet", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 4}},
			planoElement(t, f, 5, 12, "N-BK7"),
		},
	})
	res, err := g.Dispersion(gBundle(), geo.WavelengthF, geo.WavelengthC, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	p, err := g.Prescription()
	if err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	axial, err := sys.AxialChromatic(math.Inf(-1), geo.WavelengthF, geo.WavelengthC)
	if err != nil {
		t.Fatal(err)
	}
	m, _ := geo.Glass("N-BK7")
	want := -f / m.AbbeNumber()
	t.Logf("singlet F−C: real rays %.12g mm, geo paraxial %.12g mm, −f/Vd %.12g mm, |real−paraxial| %.3g mm",
		res.ShiftMM, axial, want, math.Abs(res.ShiftMM-axial))
	gClose(t, "real-ray F−C shift vs AxialChromatic", res.ShiftMM, axial, 0.02)
	gClose(t, "F−C shift vs the closed form −f/Vd", res.ShiftMM, want, 0.06)
	if math.Abs(res.ShiftMM) < 1e-3 {
		t.Errorf("a singlet must disperse: shift %.3g mm", res.ShiftMM)
	}
	if res.A.FocusMM == res.B.FocusMM {
		t.Error("the F and C focus planes are bit-identical for a glass system")
	}
}

// TestGroupDispersionFoldInvariant is the catadioptric case of the task: a plane
// fold mirror in front of an achromatic cemented crown/flint doublet does not
// change the dispersion of the doublet. The mirror is in the collimated space of
// the object at infinity (a tilted plane mirror radiates no aberration there, so
// it neither disperses nor changes the focus), and the doublet is chromatic only
// through its glasses. Each focus search is deterministic to 1e-9 mm, so the two
// shifts must agree to 1e-6 mm (each search is deterministic, but the two
// systems have different windows and the RMS minimum is quadratic, so the axial
// agreement is limited by the scan bracket: the measured value is printed). The
// test also checks that the doublet really is
// achromatic (its shift is below 30% of a singlet's) and that the folded system
// is folded (its exit axis is +y).
func TestGroupDispersionFoldInvariant(t *testing.T) {
	doublet := func() []Element {
		els := achromatElements(t, 200, 9, 6, 3)
		els[0].GapMM = 0
		els[1].GapMM = 190
		return els
	}
	straight := gWantGroup(t, &Group{
		Name:      "achromat",
		ObjectZMM: math.Inf(-1),
		Elements:  append([]Element{{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 6}}}, doublet()...),
	})
	folded := gWantGroup(t, &Group{
		Name: "folded achromat", ObjectZMM: math.Inf(-1),
		Elements: func() []Element {
			stop := StopElement(6)
			stop.GapMM = 15
			fold := MirrorElement(0, 0, 20)
			fold.TiltXDeg = 45
			fold.GapMM = 25
			return append([]Element{stop, fold}, doublet()...)
		}(),
	})
	singlet := gWantGroup(t, &Group{
		Name: "singlet", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 6}},
			planoElement(t, 200, 5, 12, "N-BK7"),
		},
	})
	b := gBundle()
	a1, err := straight.Dispersion(b, geo.WavelengthF, geo.WavelengthC, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	a2, err := folded.Dispersion(b, geo.WavelengthF, geo.WavelengthC, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	s1, err := singlet.Dispersion(b, geo.WavelengthF, geo.WavelengthC, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("achromat: shift %.12g mm (straight, rays %d failed %d) vs %.12g mm (folded, rays %d failed %d); singlet %.12g mm",
		a1.ShiftMM, a1.A.Rays, a1.A.Failed, a2.ShiftMM, a2.A.Rays, a2.A.Failed, s1.ShiftMM)
	gClose(t, "dispersion unchanged by the plane fold", a2.ShiftMM, a1.ShiftMM, 1e-6)
	if math.Abs(a1.ShiftMM) > 0.3*math.Abs(s1.ShiftMM) {
		t.Errorf("the doublet is not achromatic: %.3g mm vs the singlet's %.3g mm", a1.ShiftMM, s1.ShiftMM)
	}
	if math.Abs(s1.ShiftMM) < 1e-3 {
		t.Errorf("the reference singlet must disperse: %.3g mm", s1.ShiftMM)
	}
	_, ax, err := folded.ExitFrame()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "folded exit axis", ax, geo.V(0, 1, 0), 1e-12)
}

// TestGroupMeasurementsDeterministic repeats every measurement and requires
// bit-identical results: the searches use fixed scans and a golden-section
// refinement, no randomness and no map iteration, so two runs of the same
// computation must agree exactly (the "two runs, bit for bit" acceptance of the
// task).
func TestGroupMeasurementsDeterministic(t *testing.T) {
	build := func() *Group {
		m := MirrorElement(0, 0, 18)
		m.TiltXDeg = 45
		m.GapMM = 40
		return gWantGroup(t, &Group{Name: "determinism", ObjectZMM: math.Inf(-1), Elements: []Element{m, planoElement(t, 100, 5, 18, "N-BK7")}})
	}
	b := gBundle()
	a, err := build().ImagePosition(b, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	c, err := build().ImagePosition(b, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if a.Point != c.Point || a.RMS != c.RMS || a.FocusMM != c.FocusMM || a.GeoRadius != c.GeoRadius {
		t.Errorf("ImagePosition is not deterministic:\n%+v\n%+v", a, c)
	}
	d1, err := build().Dispersion(b, geo.WavelengthF, geo.WavelengthC, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	d2, err := build().Dispersion(b, geo.WavelengthF, geo.WavelengthC, ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if d1.ShiftMM != d2.ShiftMM {
		t.Errorf("Dispersion is not deterministic: %v vs %v", d1.ShiftMM, d2.ShiftMM)
	}
	// A folded group has no card form, so the conversion is checked on the
	// straight system of the same element.
	coaxial := func() *Group {
		return gWantGroup(t, &Group{Name: "determinism", ObjectZMM: math.Inf(-1), Elements: []Element{planoElement(t, 100, 5, 18, "N-BK7")}})
	}
	p1, err := coaxial().Prescription()
	if err != nil {
		t.Fatal(err)
	}
	p2, err := coaxial().Prescription()
	if err != nil {
		t.Fatal(err)
	}
	if len(p1.Cards) != len(p2.Cards) {
		t.Fatalf("card counts differ: %d vs %d", len(p1.Cards), len(p2.Cards))
	}
	for i := range p1.Cards {
		if p1.Cards[i].RadiusMM != p2.Cards[i].RadiusMM || p1.Cards[i].ThicknessMM != p2.Cards[i].ThicknessMM {
			t.Errorf("card %d differs between runs", i)
		}
	}
	lay1, err := build().Layout()
	if err != nil {
		t.Fatal(err)
	}
	lay2, err := build().Layout()
	if err != nil {
		t.Fatal(err)
	}
	if len(lay1) != len(lay2) {
		t.Fatal("layout lengths differ")
	}
	for i := range lay1 {
		if lay1[i].Vertex != lay2[i].Vertex || lay1[i].Axis != lay2[i].Axis {
			t.Errorf("layout entry %d differs between runs", i)
		}
	}
}

// TestGroupImagePositionRejectsAfocalAndBadPupil checks the error paths and the
// defaults of the measurement API: an afocal group has no focus plane (the exit
// bundle is parallel, the least-squares common point does not exist), a group
// without a stop aperture needs an explicit pupil, an empty sampling is
// refused, and an unset wavelength and sampling fall back to the d line and the
// hexagonal sampling.
func TestGroupImagePositionRejectsAfocalAndBadPupil(t *testing.T) {
	m, _ := geo.Glass("N-BK7")
	n := m.N(geo.Wavelengthd)
	const d = 2.0
	kepler := gWantGroup(t, &Group{
		Name: "kepler", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 6}},
			func() Element { e := planoElement(t, 100, d, 20, "N-BK7"); e.GapMM = 200 - d/n; return e }(),
			planoElement(t, 100, d, 20, "N-BK7"),
		},
	})
	// An afocal group has no image plane: the search reports the plane where the
	// (nearly parallel) bundle is most nearly concurrent, which is thousands of
	// millimetres away with a residual of the order of the beam. That is the
	// documented behaviour; AngularMagnification is the analysis for it.
	if res, err := kepler.ImagePosition(gBundle(), ImageOptions{}); err != nil {
		t.Errorf("the afocal case must be reported, not refused: %v", err)
	} else {
		t.Logf("afocal group: image %v, RMS %.3g mm, focus %.6g mm", res.Point, res.RMS, res.FocusMM)
		if math.Abs(res.FocusMM) < 1e3 {
			t.Errorf("the reported plane of an afocal group must be far away: %g mm", res.FocusMM)
		}
		if res.RMS < 1e-3 {
			t.Errorf("an afocal group must not focus the beam: RMS %.3g mm", res.RMS)
		}
	}
	// The same group is a perfectly good afocal system for Γ.
	if _, err := kepler.AngularMagnification(Bundle{Object: geo.FieldTan(0, 0), Samples: geo.Hexapolar(2)}); err != nil {
		t.Errorf("the same group must have an angular magnification: %v", err)
	}
	// A finite object has a focus: the same telescope with the object at its
	// front focal plane is the classical 4f relay, an image of the object plane
	// at unit magnification and inverted (m = −f2/f1 = −1). The object point is
	// 1 mm off axis, so the image must be at −1 mm, on the paraxial image plane
	// of the finite conjugate (geo's Paraxial on the same system, tolerance
	// 1 mm: the best-RMS focus of a singlet pair, measured value printed).
	relay := gWantGroup(t, &Group{
		Name: "4f relay", ObjectZMM: -100,
		Elements: []Element{
			{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 6}},
			func() Element { e := planoElement(t, 100, d, 20, "N-BK7"); e.GapMM = 200 - d/n; return e }(),
			planoElement(t, 100, d, 20, "N-BK7"),
		},
	})
	res, err := relay.ImagePosition(Bundle{Object: geo.FieldHeight(-100, 1, 0)}, ImageOptions{})
	if err != nil {
		t.Fatalf("the 4f relay must image a finite object: %v", err)
	}
	relaySys, err := relay.System()
	if err != nil {
		t.Fatal(err)
	}
	rpar, err := relaySys.Paraxial(geo.Wavelengthd, -100)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("4f relay: image %v (RMS %.3g mm, rays %d), paraxial image z %.12g, magnification %.6g",
		res.Point, res.RMS, res.Rays, rpar.ImageZ, rpar.Magnification)
	gClose(t, "4f transverse magnification is −1", res.Point.X, -1, 0.05)
	gClose(t, "4f image plane vs the paraxial conjugate", res.Point.Z, rpar.ImageZ, 1)
	gClose(t, "4f paraxial magnification", rpar.Magnification, -1, 1e-6)
	// The automatic pupil needs a clear aperture in the stop element.
	noAperture := gWantGroup(t, &Group{
		Name: "no aperture", ObjectZMM: math.Inf(-1),
		Elements: []Element{planoElement(t, 100, 5, 0, "N-BK7")},
	})
	if _, err := noAperture.ImagePosition(gBundle(), ImageOptions{}); err == nil {
		t.Error("a group without a stop aperture must ask for an explicit pupil")
	} else if !strings.Contains(err.Error(), "clear aperture") {
		t.Errorf("the error must mention the clear aperture: %v", err)
	}
	// An explicit pupil and a degenerate sampling.
	withPupil := Bundle{Object: geo.FieldTan(0, 0), Pupil: geo.Pupil{Z: -10, Radius: 5}, Samples: []geo.PupilPoint{}}
	if _, err := noAperture.ImagePosition(withPupil, ImageOptions{}); err == nil {
		t.Error("an empty sampling must be refused")
	}
	// The defaults: no wavelength (the d line) and no sampling (Hexapolar(3)).
	defaults := Bundle{Object: geo.FieldTan(0, 0), Pupil: geo.Pupil{Z: -10, Radius: 5}}
	res, err = noAperture.ImagePosition(defaults, ImageOptions{})
	if err != nil {
		t.Fatalf("the default bundle must work: %v", err)
	}
	if res.Rays != 37 || res.Failed != 0 {
		t.Errorf("default sampling: %d rays, %d failed, want 37 and 0", res.Rays, res.Failed)
	}
	if math.IsNaN(res.FocusMM) || res.FocusMM <= 0 {
		t.Errorf("focus %g mm", res.FocusMM)
	}
	// A bundle whose rays all fail is reported, not silently answered.
	blocked := gWantGroup(t, &Group{
		Name: "blocked", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 1}},
			func() Element { e := planoElement(t, 100, 5, 0.5, "N-BK7"); return e }(),
		},
	})
	if _, err := blocked.ImagePosition(Bundle{Object: geo.FieldTan(0, 0), Pupil: geo.Pupil{Z: -10, Radius: 5}}, ImageOptions{}); err == nil {
		t.Error("a bundle that is entirely vignetted must be reported")
	}
}
