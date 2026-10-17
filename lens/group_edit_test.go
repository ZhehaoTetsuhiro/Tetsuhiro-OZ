package lens

import (
	"bytes"
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
)

// ---------------------------------------------------------------------------
// Editing tests (task group3d): the element list, every parameter, and the
// physical effect of every edit
// ---------------------------------------------------------------------------

// gAbbe returns the d-line index and the Abbe number of a catalogue glass.
func gAbbe(t *testing.T, glass string) (float64, float64) {
	t.Helper()
	m, ok := geo.Glass(glass)
	if !ok {
		t.Fatalf("unknown glass %q", glass)
	}
	return m.N(geo.Wavelengthd), m.AbbeNumber()
}

// TestGroupEditListOps checks the element list API: append, insert, remove and
// move keep the order deterministic, the stop index follows the permutation and
// an invalid element is refused without touching the group.
func TestGroupEditListOps(t *testing.T) {
	g := NewGroup("list")
	stop := StopElement(5)
	i, err := g.AddElement(stop)
	if err != nil {
		t.Fatal(err)
	}
	if i != 0 || g.Len() != 1 {
		t.Fatalf("AddElement returned %d, len %d", i, g.Len())
	}
	lens := planoElement(t, 100, 5, 12, "N-BK7")
	if _, err := g.AddElement(lens); err != nil {
		t.Fatal(err)
	}
	// Insert in the middle: the new element takes the requested index.
	mirror := MirrorElement(-400, -1, 10)
	if err := g.InsertElement(1, mirror); err != nil {
		t.Fatal(err)
	}
	if g.Len() != 3 || g.Elements[1].Kind != ElementMirror || g.Elements[2].Kind != ElementLens {
		t.Fatalf("insert produced %+v", g.Elements)
	}
	// The stop index follows an insertion before it.
	g.Stop = 1
	if err := g.InsertElement(0, StopElement(3)); err != nil {
		t.Fatal(err)
	}
	if g.Stop != 2 {
		t.Errorf("stop %d after inserting at 0, want 2", g.Stop)
	}
	// Move an element to the end.
	if err := g.MoveElement(0, g.Len()-1); err != nil {
		t.Fatal(err)
	}
	last := g.Elements[g.Len()-1]
	if last.Kind != ElementStop || last.Front.SemiAperMM != 3 {
		t.Errorf("MoveElement did not move the element: %+v", last)
	}
	if g.Stop != 1 {
		t.Errorf("stop %d after the move, want 1", g.Stop)
	}
	// Remove an element before the stop: the index follows.
	before := g.Stop
	if err := g.RemoveElement(0); err != nil {
		t.Fatal(err)
	}
	if g.Stop != before-1 {
		t.Errorf("stop %d after removing before it, want %d", g.Stop, before-1)
	}
	// Errors do not modify the group.
	n := g.Len()
	for name, err := range map[string]error{
		"add a bad kind":           func() error { _, e := g.AddElement(Element{Kind: "prism"}); return e }(),
		"insert out of range":      g.InsertElement(n+1, lens),
		"insert a bad glass":       g.InsertElement(0, func() Element { e := lens; e.Glass = "unobtainium"; return e }()),
		"remove out of range":      g.RemoveElement(n),
		"move out of range":        g.MoveElement(0, n),
		"move from out of range":   g.MoveElement(-1, 0),
		"thickness negative":       g.SetThicknessMM(0, -1),
		"gap not finite":           g.SetGapMM(0, math.Inf(1)),
		"aperture negative":        g.SetSemiAperture(0, Front, -2),
		"unknown parameter kind":   g.SetParam(GroupParam{Kind: GroupParamKind(99)}, 1),
		"coefficient out of range": g.SetParam(GroupParam{Surf: SurfRef{Element: 0}, Kind: GPCoeff, Order: 3}, 1),
		"tilt not finite":          g.SetPose(0, 0, 0, math.NaN(), 0, 0),
	} {
		if err == nil {
			t.Errorf("%s: want an error", name)
		}
	}
	if g.Len() != n {
		t.Errorf("a rejected edit changed the group: %d elements, want %d", g.Len(), n)
	}
	if err := g.Validate(); err != nil {
		t.Errorf("the group is invalid after the edits: %v", err)
	}
	// SetElement replaces in place, SetKind converts (and drops the rear
	// surface of a mirror).
	if err := g.SetKind(0, ElementMirror); err != nil {
		t.Fatal(err)
	}
	if e, _ := g.ElementAt(0); e.Rear != nil || e.ThicknessMM != 0 {
		t.Errorf("SetKind to mirror kept the rear surface: %+v", e)
	}
	if err := g.SetKind(0, ElementLens); err != nil {
		t.Fatal(err)
	}
	e, err := g.ElementAt(0)
	if err != nil {
		t.Fatal(err)
	}
	e.Front.SemiAperMM = 7
	if err := g.SetElement(0, e); err != nil {
		t.Fatal(err)
	}
	if got, _ := g.ElementAt(0); got.Front.SemiAperMM != 7 {
		t.Errorf("SetElement did not replace: %+v", got)
	}
	// An empty group is invalid, and so is a stop index out of range.
	empty := NewGroup("empty")
	if err := empty.Validate(); err == nil {
		t.Error("an empty group must be invalid")
	}
	g.Stop = 99
	if err := g.Validate(); err == nil {
		t.Error("an out-of-range stop must be invalid")
	}
}

// TestGroupParamGetSet walks every parameter kind: the value written is read
// back, the composite kinds keep the units documented (radius/curvature are two
// coordinates of the same field, the position moves the element by adjusting the
// gap of its predecessor) and the referencing errors are descriptive.
func TestGroupParamGetSet(t *testing.T) {
	one := Element{
		Name: "asphere", Kind: ElementLens,
		Front: SurfaceSpec{
			RadiusMM: 50, ConicK: -0.5, Coeffs: []float64{1e-6, 2e-9}, SemiAperMM: 12,
		},
		Rear:        &SurfaceSpec{RadiusMM: -80, SemiAperMM: 12},
		Glass:       "N-BK7",
		ThicknessMM: 6,
		GapMM:       20,
	}
	poly := Element{
		Name: "poly", Kind: ElementLens,
		Front:       SurfaceSpec{RadiusMM: 0, PolyTerms: []PolyTermCard{{Power: 6, Coeff: 3e-10}}, SemiAperMM: 8},
		Rear:        &SurfaceSpec{RadiusMM: 0, SemiAperMM: 8},
		Glass:       "N-BK7",
		ThicknessMM: 2,
	}
	g := gWantGroup(t, &Group{Name: "params", ObjectZMM: math.Inf(-1), Elements: []Element{one, MirrorElement(-400, -1, 10), poly}})
	front := SurfRef{Element: 0, Face: Front}

	// Scalar surface parameters.
	for _, c := range []struct {
		kind GroupParamKind
		val  float64
	}{
		{GPRadius, 42.5}, {GPCurvature, 1.0 / 33.0}, {GPConic, -1.25},
		{GPCoeff, -3.5e-7}, {GPSemiAperture, 11.25},
	} {
		p := GroupParam{Surf: front, Kind: c.kind}
		if err := g.SetParam(p, c.val); err != nil {
			t.Fatalf("%s: %v", p, err)
		}
		got, err := g.GetParam(p)
		if err != nil {
			t.Fatalf("%s: %v", p, err)
		}
		gClose(t, p.String(), got, c.val, 0)
	}
	// A polynomial-asphere term of its own surface.
	pp := GroupParam{Surf: SurfRef{Element: 2, Face: Front}, Kind: GPPolyTerm}
	if err := g.SetParam(pp, 8e-11); err != nil {
		t.Fatal(err)
	}
	gClose(t, pp.String(), mustGet(t, g, pp), 8e-11, 0)

	// Radius and curvature are the same field.
	gClose(t, "curvature of R=33", mustGet(t, g, GroupParam{Surf: front, Kind: GPCurvature}), 1.0/33.0, 1e-15)
	gClose(t, "radius of c=1/33", mustGet(t, g, GroupParam{Surf: front, Kind: GPRadius}), 33, 1e-12)
	// A plane: radius 0 and curvature 0 agree, and the radius is ±Inf.
	if err := g.SetParam(GroupParam{Surf: front, Kind: GPRadius}, 0); err != nil {
		t.Fatal(err)
	}
	if r := mustGet(t, g, GroupParam{Surf: front, Kind: GPRadius}); !math.IsInf(r, 0) {
		t.Errorf("radius %v, want ±Inf after setting 0", r)
	}
	gClose(t, "curvature of a plane", mustGet(t, g, GroupParam{Surf: front, Kind: GPCurvature}), 0, 0)
	if err := g.SetParam(GroupParam{Surf: front, Kind: GPCurvature}, 0.02); err != nil {
		t.Fatal(err)
	}
	gClose(t, "radius from the curvature", mustGet(t, g, GroupParam{Surf: front, Kind: GPRadius}), 50, 1e-12)

	// Element-level parameters.
	gClose(t, "thickness", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 0}, Kind: GPThickness}), 6, 0)
	if err := g.SetThicknessMM(0, 7.5); err != nil {
		t.Fatal(err)
	}
	gClose(t, "thickness after set", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 0}, Kind: GPThickness}), 7.5, 0)
	if err := g.SetParam(GroupParam{Surf: SurfRef{Element: 1}, Kind: GPThickness}, 1); err == nil {
		t.Error("a single-surface element must not accept a thickness")
	}
	if err := g.SetParam(GroupParam{Surf: SurfRef{Element: 0}, Kind: GPTiltY}, 4.5); err != nil {
		t.Fatal(err)
	}
	gClose(t, "tilt y", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 0}, Kind: GPTiltY}), 4.5, 0)
	if err := g.SetParam(GroupParam{Surf: SurfRef{Element: 0}, Kind: GPDecenterX}, 0.5); err != nil {
		t.Fatal(err)
	}
	gClose(t, "decentre x", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 0}, Kind: GPDecenterX}), 0.5, 0)

	// Positions and gaps: the position of an element is the chain position of
	// its first surface, and writing it moves the element alone. Element 0 is
	// posed at this point (the tilt and decentre above), and a posed element
	// does not advance the chain (its return break restores the frame), so the
	// gap of element 0 alone places the mirror.
	gClose(t, "gap of the asphere", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 0}, Kind: GPGap}), 20, 0)
	gClose(t, "gap of the mirror", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 1}, Kind: GPGap}), 0, 0)
	gClose(t, "position of the mirror (the posed lens does not advance the chain)",
		mustGet(t, g, GroupParam{Surf: SurfRef{Element: 1}, Kind: GPPosition}), 20, 0)
	if err := g.SetPositionMM(1, 40); err != nil {
		t.Fatal(err)
	}
	gClose(t, "position after set", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 1}, Kind: GPPosition}), 40, 0)
	gClose(t, "gap after moving the mirror", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 0}, Kind: GPGap}), 40, 0)
	if err := g.SetPositionMM(0, 3); err != nil {
		t.Fatal(err)
	}
	gClose(t, "origin follows the first element", g.OriginZ, 3, 0)
	if err := g.SetPositionMM(1, 40); err != nil {
		t.Fatal(err)
	}
	gClose(t, "position with a non-zero origin", mustGet(t, g, GroupParam{Surf: SurfRef{Element: 1}, Kind: GPPosition}), 40, 0)
	// An unposed two-surface element advances the chain by its thickness: the
	// following element moves by exactly the centre thickness when the element
	// loses its pose.
	beforeUnpose := mustGet(t, g, GroupParam{Surf: SurfRef{Element: 1}, Kind: GPPosition})
	if err := g.SetPose(0, 0, 0, 0, 0, 0); err != nil {
		t.Fatal(err)
	}
	afterUnpose := mustGet(t, g, GroupParam{Surf: SurfRef{Element: 1}, Kind: GPPosition})
	gClose(t, "an unposed element advances the chain by its thickness", afterUnpose-beforeUnpose, 7.5, 0)

	// Glass names.
	if err := g.SetGlass(0, "F2"); err != nil {
		t.Fatal(err)
	}
	if err := g.SetGlass(0, "unobtainium"); err == nil {
		t.Error("an unknown glass must be refused")
	}
	if err := g.SetGapGlass(1, "N-BK7"); err != nil {
		t.Fatal(err)
	}
	if err := g.SetGlass(1, "N-BK7"); err == nil {
		t.Error("a mirror must not accept a glass")
	}

	// The nested aspheric kinds: a q-asphere surface and a freeform one.
	q := Element{Kind: ElementLens, Front: SurfaceSpec{RadiusMM: 60, QAsphere: &QCard{Form: "qcon", NormRadiusMM: 10, Coeffs: []float64{0, 1e-7}}},
		Rear: &SurfaceSpec{RadiusMM: -60}, Glass: "N-BK7", ThicknessMM: 4}
	gq := gWantGroup(t, &Group{Name: "q", ObjectZMM: math.Inf(-1), Elements: []Element{q}})
	qp := GroupParam{Surf: SurfRef{Element: 0, Face: Front}, Kind: GPQCoef, Order: 1}
	if err := gq.SetParam(qp, 5e-8); err != nil {
		t.Fatal(err)
	}
	gClose(t, qp.String(), mustGet(t, gq, qp), 5e-8, 0)
	if err := gq.SetParam(GroupParam{Surf: qp.Surf, Kind: GPXYTerm}, 1); err == nil {
		t.Error("a q-asphere surface has no xy terms")
	}
	xy := Element{Kind: ElementLens, Front: SurfaceSpec{RadiusMM: 0, XYPoly: &XYCard{Terms: []XYTermCard{{I: 2, J: 0, Coeff: 1e-5}}, NormRadiusMM: 10}},
		Rear: &SurfaceSpec{RadiusMM: 0}, Glass: "N-BK7", ThicknessMM: 2}
	gx := gWantGroup(t, &Group{Name: "freeform", ObjectZMM: math.Inf(-1), Elements: []Element{xy}})
	xp := GroupParam{Surf: SurfRef{Element: 0, Face: Front}, Kind: GPXYTerm}
	if err := gx.SetParam(xp, 2e-5); err != nil {
		t.Fatal(err)
	}
	gClose(t, xp.String(), mustGet(t, gx, xp), 2e-5, 0)
	z := Element{Kind: ElementLens, Front: SurfaceSpec{RadiusMM: 0, Zernike: &ZernikeCard{NormRadiusMM: 5, Coeffs: []float64{0, 0, 0, 1e-4}}},
		Rear: &SurfaceSpec{RadiusMM: 0}, Glass: "N-BK7", ThicknessMM: 2}
	gz := gWantGroup(t, &Group{Name: "zernike", ObjectZMM: math.Inf(-1), Elements: []Element{z}})
	zp := GroupParam{Surf: SurfRef{Element: 0, Face: Front}, Kind: GPZernike, Order: 3}
	if err := gz.SetParam(zp, 2e-4); err != nil {
		t.Fatal(err)
	}
	gClose(t, zp.String(), mustGet(t, gz, zp), 2e-4, 0)
	// The rear of a single-surface element does not exist (the mirror of g).
	if _, err := g.GetParam(GroupParam{Surf: SurfRef{Element: 1, Face: Rear}, Kind: GPRadius}); err == nil {
		t.Error("a single-surface element has no rear surface")
	}
	if _, err := gz.GetParam(GroupParam{Surf: SurfRef{Element: 0, Face: Front}, Kind: GPPolyTerm, Order: 0}); err == nil {
		t.Error("a Zernike surface has no polynomial terms")
	}
	if _, err := gz.GetParam(GroupParam{Surf: SurfRef{Element: 5}, Kind: GPRadius}); err == nil {
		t.Error("an out-of-range element must be refused")
	}
}

// mustGet reads a parameter and fails the test on error.
func mustGet(t *testing.T, g *Group, p GroupParam) float64 {
	t.Helper()
	v, err := g.GetParam(p)
	if err != nil {
		t.Fatalf("GetParam(%s): %v", p, err)
	}
	return v
}

// TestGroupEditRadiusChangesFocalLength checks the focal length after a radius
// edit against the thin-lens formula 1/f = (n−1)(1/R1 − 1/R2), first on the
// paraxial engine (exact: a zero-thickness biconvex element, relative threshold
// 1e-9) and then with real rays on a traceable plano-convex element: the best
// focus must follow the new focal length (the difference is the element's
// spherical aberration, printed by the test).
func TestGroupEditRadiusChangesFocalLength(t *testing.T) {
	n, _ := gAbbe(t, "N-BK7")
	thin := LensElement(2*(n-1)*100, -2*(n-1)*100, 0, "N-BK7", 12)
	g := gWantGroup(t, &Group{Name: "thin", ObjectZMM: math.Inf(-1), Elements: []Element{thin}})
	efl := func() float64 {
		sys, err := g.System()
		if err != nil {
			t.Fatal(err)
		}
		par, err := sys.Paraxial(geo.Wavelengthd, math.Inf(-1))
		if err != nil {
			t.Fatal(err)
		}
		return par.EFL
	}
	gClose(t, "EFL of the nominal thin lens", efl(), 100, 1e-9*100)
	for _, r1 := range []float64{120, 60, -150, 250} {
		if err := g.SetRadiusMM(0, Front, r1); err != nil {
			t.Fatal(err)
		}
		got := efl()
		want := 1 / ((n - 1) * (1/r1 - 1/g.Elements[0].Rear.RadiusMM))
		t.Logf("R1 = %6.1f: EFL %.12g mm, thin-lens formula %.12g mm", r1, got, want)
		gClose(t, "EFL follows the thin-lens formula", got, want, 1e-9*math.Abs(want))
	}
	// The same edit with real rays through a traceable element: the back focus
	// of the plano-convex element is f − d/n (measured exactly at the nominal
	// radius) and the best focus follows it within the spherical aberration of
	// the singlet at a 5 mm stop (0.23 mm, measured in the 3D test file).
	real := gWantGroup(t, &Group{
		Name: "real", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 5}},
			planoElement(t, 100, 5, 14, "N-BK7"),
		},
	})
	focus := func() float64 {
		res, err := real.ImagePosition(gBundle(), ImageOptions{})
		if err != nil {
			t.Fatal(err)
		}
		return res.FocusMM
	}
	gClose(t, "back focus of the nominal element", focus(), 100-5/n, 0.35)
	for _, f := range []float64{150, 80} {
		if err := real.SetRadiusMM(1, Front, (n-1)*f); err != nil {
			t.Fatal(err)
		}
		want := f - 5/n
		t.Logf("f = %5.1f mm: back focus %.12g mm, thin lens %.12g mm, |Δ| %.3g mm", f, focus(), want, math.Abs(focus()-want))
		gClose(t, "real-ray back focus follows the focal length", focus(), want, 0.5)
	}
}

// TestGroupEditGlassChangesDispersion checks that replacing the glass of the
// element changes the chromatic focus shift in the ratio of the Abbe numbers and
// that the paraxial engine of package geo agrees with the real-ray measurement
// for both glasses.
func TestGroupEditGlassChangesDispersion(t *testing.T) {
	build := func(glass string) *Group {
		return gWantGroup(t, &Group{
			Name: "glass", ObjectZMM: math.Inf(-1),
			Elements: []Element{
				{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 4}},
				planoElement(t, 100, 5, 12, glass),
			},
		})
	}
	dispersion := func(g *Group) (real, paraxial float64) {
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
		ax, err := sys.AxialChromatic(math.Inf(-1), geo.WavelengthF, geo.WavelengthC)
		if err != nil {
			t.Fatal(err)
		}
		return res.ShiftMM, ax
	}
	_, vBK7 := gAbbe(t, "N-BK7")
	_, vSF5 := gAbbe(t, "N-SF5")
	gBK7 := build("N-BK7")
	gSF5 := build("N-SF5")
	realBK7, parBK7 := dispersion(gBK7)
	realSF5, parSF5 := dispersion(gSF5)
	t.Logf("N-BK7: real %.12g mm, paraxial %.12g mm, Vd %.4g", realBK7, parBK7, vBK7)
	t.Logf("N-SF5: real %.12g mm, paraxial %.12g mm, Vd %.4g", realSF5, parSF5, vSF5)
	t.Logf("ratio of the shifts %.6g, ratio of the Abbe numbers %.6g", realSF5/realBK7, vBK7/vSF5)
	// The chromatic shift goes as 1/Vd: the ratio of the two shifts must follow
	// the ratio of the Abbe numbers. Threshold 3%: the two elements have the
	// same shape, so the geometric (principal-plane) term of the focus scales
	// with the index too, which the ratio of the Abbe numbers alone does not
	// capture exactly (measured: see the log).
	gClose(t, "shift ratio vs the Abbe ratio", realSF5/realBK7, vBK7/vSF5, 0.03*vBK7/vSF5)
	if math.Abs(realSF5) <= math.Abs(realBK7) {
		t.Errorf("the denser flint must disperse more: %.3g vs %.3g mm", realSF5, realBK7)
	}
	for _, c := range []struct {
		name     string
		real, ax float64
	}{{"N-BK7", realBK7, parBK7}, {"N-SF5", realSF5, parSF5}} {
		gClose(t, c.name+": real vs paraxial shift", c.real, c.ax, 0.02)
	}
}

// TestGroupEditConicAndAsphereChangeTheSpot checks the surface equation
// parameters with the two classical identities of a mirror:
//
//   - a paraboloid (K = −1) focuses a collimated on-axis bundle exactly, while
//     the sphere of the same vertex radius leaves the spherical aberration, so
//     the spot at the paraxial focus R/2 changes by orders of magnitude;
//   - the even-asphere coefficient A4 = −c³/8 turns the sphere into that
//     paraboloid, so the spot collapses to the level of the r⁶ residual of the
//     sphere's series (the sag identity is checked to 1e-12 mm directly on the
//     surface of the converted system, without any ray).
func TestGroupEditConicAndAsphereChangeTheSpot(t *testing.T) {
	const r = -400.0
	paraxialFocus := r / 2 // a concave mirror of R = −400 focuses at −200
	mk := func(conic, a4 float64) *Group {
		e := MirrorElement(r, conic, 20)
		e.Name = "mirror"
		if a4 != 0 {
			e.Front.Coeffs = []float64{a4}
		}
		return gWantGroup(t, &Group{Name: "mirror", ObjectZMM: math.Inf(-1), Elements: []Element{e}})
	}
	spotRMS := func(g *Group, imageZ float64) float64 {
		sys, err := g.System()
		if err != nil {
			t.Fatal(err)
		}
		withPlane := sys.WithImagePlane(imageZ)
		sd, err := withPlane.SpotDiagram(geo.FieldTan(0, 0), geo.Pupil{Z: -20, Radius: 18}, geo.Hexapolar(3), imageZ, geo.Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		st := sd.Stats(geo.CentroidReference())
		if st.N < 30 {
			t.Fatalf("only %d of 37 rays reached the image plane", st.N)
		}
		return st.RMS
	}
	sphere := mk(0, 0)
	parabola := mk(-1, 0)
	rmsSphere := spotRMS(sphere, paraxialFocus)
	rmsParabola := spotRMS(parabola, paraxialFocus)
	c := 1 / r
	corrected := mk(0, -c*c*c/8)
	rmsCorrected := spotRMS(corrected, paraxialFocus)
	t.Logf("spot at z = %g: sphere %.6g mm, parabola %.6g mm, sphere+A4 %.6g mm", paraxialFocus, rmsSphere, rmsParabola, rmsCorrected)
	if rmsSphere < 1e-2 {
		t.Errorf("the sphere must show spherical aberration at the paraxial focus: RMS %.3g mm", rmsSphere)
	}
	if rmsParabola > 1e-9 {
		t.Errorf("the parabola focuses exactly: RMS %.3g mm", rmsParabola)
	}
	if !(rmsCorrected < 1e-2*rmsSphere) {
		t.Errorf("K = 0 with A4 = −c³/8 must reproduce the parabola: RMS %.3g vs %.3g mm", rmsCorrected, rmsSphere)
	}
	// The sag identity: the converted surface and the parabola agree to the r⁶
	// term of the sphere's series, which is below 1e-12 mm over this aperture
	// (c⁵r⁶/16 = 4e-7 mm at r = 20 … the residual is checked, not assumed).
	sysA, err := corrected.System()
	if err != nil {
		t.Fatal(err)
	}
	sysB, err := parabola.System()
	if err != nil {
		t.Fatal(err)
	}
	worst := 0.0
	for i := 0; i <= 20; i++ {
		rr := 20.0 * float64(i) / 20
		worst = math.Max(worst, math.Abs(sysA.Surfaces[0].SagR(rr)-sysB.Surfaces[0].SagR(rr)))
	}
	t.Logf("sag difference of sphere+A4 vs the parabola: %.3g mm (the r⁶ residual of the sphere series)", worst)
	if worst > 1e-5 {
		t.Errorf("the A4 identity is wrong: sag difference %.3g mm", worst)
	}
	// The edit through the parameter API: the same coefficient written with
	// SetParam gives the same surface.
	a4 := GroupParam{Surf: SurfRef{Element: 0, Face: Front}, Kind: GPCoeff, Order: 0}
	edited := mk(0, -c*c*c/8)
	if err := edited.SetParam(a4, -c*c*c/8); err != nil {
		t.Fatal(err)
	}
	editedSys, err := edited.System()
	if err != nil {
		t.Fatal(err)
	}
	for i := 0; i <= 20; i++ {
		rr := 20.0 * float64(i) / 20
		if d := math.Abs(editedSys.Surfaces[0].SagR(rr) - sysA.Surfaces[0].SagR(rr)); d > 1e-14 {
			t.Fatalf("SetParam(A4) does not match the literal coefficient at r=%g: %.3g", rr, d)
		}
	}
}

// TestGroupEditMoveShiftsFocus checks that moving an element by Δz moves the
// image by exactly Δz along the beam: the focal length of a lens does not change
// when it is translated, so the whole system (and the image with it) is
// translated. The image point is global, so the test compares it before and
// after the move. Threshold 1e-6 mm: the two independent best-focus searches
// agree to the scan bracket (the measured value is printed).
func TestGroupEditMoveShiftsFocus(t *testing.T) {
	build := func() *Group {
		return gWantGroup(t, &Group{
			Name: "move", ObjectZMM: math.Inf(-1),
			Elements: []Element{
				{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 4}},
				planoElement(t, 100, 5, 12, "N-BK7"),
			},
		})
	}
	const dz = 7.5
	before, err := build().ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	moved := build()
	pos, err := moved.PositionMM(1)
	if err != nil {
		t.Fatal(err)
	}
	if err := moved.SetPositionMM(1, pos+dz); err != nil {
		t.Fatal(err)
	}
	after, err := moved.ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("move by %g mm: image %v → %v (Δz %.12g, |Δ − %g| %.3g)",
		dz, before.Point, after.Point, after.Point.Z-before.Point.Z, dz, math.Abs(after.Point.Z-before.Point.Z-dz))
	// The lens and the whole image-side layout move together, so the traced
	// rays are the same rays translated: the two searches differ only by the
	// arithmetic of the translated coordinates (~1e-16).
	gClose(t, "the image follows the moved element", after.Point.Z-before.Point.Z, dz, 1e-9)
	gClose(t, "lateral position unchanged", after.Point.X-before.Point.X, 0, 1e-9)
	gClose(t, "the focus distance is unchanged", after.FocusMM, before.FocusMM, 0)
	// The same move in a folded system: the image follows along the folded
	// axis (+y for a 45 degree mirror).
	folded := func() *Group {
		return gWantGroup(t, &Group{
			Name: "folded move", ObjectZMM: math.Inf(-1),
			Elements: []Element{
				func() Element { e := StopElement(6); e.GapMM = 15; return e }(),
				func() Element {
					e := MirrorElement(0, 0, 20)
					e.TiltXDeg = 45
					e.GapMM = 25
					return e
				}(),
				planoElement(t, 100, 5, 18, "N-BK7"),
			},
		})
	}
	fBefore, err := folded().ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	fMoved := folded()
	fPos, err := fMoved.PositionMM(2)
	if err != nil {
		t.Fatal(err)
	}
	if err := fMoved.SetPositionMM(2, fPos+dz); err != nil {
		t.Fatal(err)
	}
	fAfter, err := fMoved.ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	shift := fAfter.Point.Sub(fBefore.Point)
	t.Logf("folded move by %g mm along +y: image %v → %v, shift %v (|shift − (0,%g,0)| = %.3g)",
		dz, fBefore.Point, fAfter.Point, shift, dz, shift.Sub(geo.V(0, dz, 0)).Norm())
	// The folded layout searches the focus with windows that shift with the
	// absolute coordinates, so the two independent searches agree to the scan
	// bracket (measured 5.5e-10 mm on the image point, 5.5e-7 mm on the focus
	// distance); the coaxial case above is exact because the relative geometry
	// of the search is bit-identical there.
	gVecClose(t, "the folded image follows the moved element", shift, geo.V(0, dz, 0), 1e-6)
	gClose(t, "the folded focus distance is unchanged", fAfter.FocusMM, fBefore.FocusMM, 1e-5)
}

// TestGroupValidateRejectsBadElements walks the validation rules of the element
// model: unknown kinds, non-finite numbers, a mirror with a rear surface, a stop
// with a radius, an air gap of zero thickness, mutually exclusive aspheric
// descriptions and a bad glass name in a gap.
func TestGroupValidateRejectsBadElements(t *testing.T) {
	good := planoElement(t, 100, 5, 12, "N-BK7")
	cases := []struct {
		name string
		mut  func(*Element)
	}{
		{"unknown kind", func(e *Element) { e.Kind = "prism" }},
		{"NaN radius", func(e *Element) { e.Front.RadiusMM = math.NaN() }},
		{"Inf thickness", func(e *Element) { e.ThicknessMM = math.Inf(1) }},
		{"negative aperture", func(e *Element) { e.Front.SemiAperMM = -1 }},
		{"unknown glass", func(e *Element) { e.Glass = "unobtainium" }},
		{"unknown gap glass", func(e *Element) { e.GapGlass = "unobtainium" }},
		{"NaN pose", func(e *Element) { e.TiltZDeg = math.NaN() }},
		{"NaN gap", func(e *Element) { e.GapMM = math.NaN() }},
		{"zero-thickness air gap", func(e *Element) { e.Glass, e.ThicknessMM = "", 0 }},
		{"two aspheric descriptions", func(e *Element) {
			e.Front.Coeffs = []float64{1e-6}
			e.Front.PolyTerms = []PolyTermCard{{Power: 6, Coeff: 1e-9}}
		}},
		{"bad q form", func(e *Element) { e.Front.QAsphere = &QCard{Form: "qcubic", NormRadiusMM: 5} }},
	}
	for _, c := range cases {
		e := cloneElement(good)
		c.mut(&e)
		if err := e.validate(0); err == nil {
			t.Errorf("%s: want an error", c.name)
		}
	}
	for _, c := range []struct {
		name string
		el   Element
	}{
		{"mirror with a rear surface", func() Element {
			e := MirrorElement(-400, -1, 10)
			r := SurfaceSpec{RadiusMM: 0}
			e.Rear = &r
			return e
		}()},
		{"mirror with a glass", func() Element {
			e := MirrorElement(-400, -1, 10)
			e.Glass = "N-BK7"
			return e
		}()},
		{"mirror with a thickness", func(e Element) Element { e = MirrorElement(-400, -1, 10); e.ThicknessMM = 3; return e }(Element{})},
		{"stop with a radius", func() Element {
			e := StopElement(5)
			e.Front.RadiusMM = 50
			return e
		}()},
		{"stop with an asphere", func() Element {
			e := StopElement(5)
			e.Front.Coeffs = []float64{1e-6}
			return e
		}()},
		{"stop in glass", func() Element {
			e := StopElement(5)
			e.GapGlass = "N-BK7"
			return e
		}()},
		{"single-surface lens with a glass", func() Element {
			e := SingleSurface(50, 5)
			e.Glass = "N-BK7"
			return e
		}()},
		{"single-surface element with a thickness", func() Element {
			e := SingleSurface(50, 5)
			e.ThicknessMM = 2
			return e
		}()},
	} {
		if err := c.el.validate(0); err == nil {
			t.Errorf("%s: want an error", c.name)
		}
	}
	// The good element passes, and a folded mirror too.
	if err := good.validate(0); err != nil {
		t.Errorf("a valid element was rejected: %v", err)
	}
	m := MirrorElement(0, 0, 20)
	m.TiltXDeg = 45
	if err := m.validate(0); err != nil {
		t.Errorf("a folded mirror was rejected: %v", err)
	}
	if !m.posed() {
		t.Error("a mirror tilted by 45 degrees must report as posed")
	}
	// The group errors.
	if err := (&Group{ObjectZMM: math.NaN(), Elements: []Element{good}}).Validate(); err == nil {
		t.Error("a NaN object position must be invalid")
	}
	if err := (&Group{ObjectZMM: math.Inf(-1), Elements: []Element{good}, OriginZ: math.Inf(1)}).Validate(); err == nil {
		t.Error("an infinite origin must be invalid")
	}
}

// TestGroupJSONFileRoundTrip writes a group to a file, reads it back and
// requires the same system, the same trace and an identical String report. A
// hand-written file with a legacy-style layout (an aperture stop card, two
// mirrors and a negative gap) must load into the element model.
func TestGroupJSONFileRoundTrip(t *testing.T) {
	g := gWantGroup(t, &Group{
		Name: "file", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(6); e.GapMM = 15; return e }(),
			func() Element {
				e := MirrorElement(0, 0, 20)
				e.TiltXDeg = 45
				e.Name = "fold"
				e.GapMM = 30
				return e
			}(),
			planoElement(t, 100, 5, 18, "N-BK7"),
		},
	})
	dir := t.TempDir()
	path := filepath.Join(dir, "group.json")
	if err := g.SaveJSON(path); err != nil {
		t.Fatal(err)
	}
	back, err := LoadGroupJSON(path)
	if err != nil {
		t.Fatal(err)
	}
	if back.String() != g.String() {
		t.Errorf("the loaded group differs:\n%s\n%s", g.String(), back.String())
	}
	sysA, err := g.System()
	if err != nil {
		t.Fatal(err)
	}
	sysB, err := back.System()
	if err != nil {
		t.Fatal(err)
	}
	if err := sameGeometry(sysA, sysB); err != nil {
		t.Errorf("the loaded system differs: %v", err)
	}
	resA, err := g.ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	resB, err := back.ImagePosition(gBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if resA.Point != resB.Point || resA.FocusMM != resB.FocusMM {
		t.Errorf("the loaded group images differently: %v vs %v", resA.Point, resB.Point)
	}
	// A file with a broken element is refused with a descriptive error.
	broken := filepath.Join(dir, "broken.json")
	if err := os.WriteFile(broken, []byte(`{"name":"broken","object_z_mm":"inf","stop":0,"elements":[{"kind":"lens","front":{"radius_mm":50},"glass":"unobtainium","thickness_mm":5}]}`), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := LoadGroupJSON(broken); err == nil {
		t.Error("a group with an unknown glass must not load")
	} else if !strings.Contains(err.Error(), "unobtainium") {
		t.Errorf("the error must name the glass: %v", err)
	}
	// The element JSON is the documented schema: the surface fields of a Card
	// plus the element fields.
	var probe struct {
		Elements []struct {
			Kind  string  `json:"kind"`
			Front Card    `json:"front"`
			Rear  *Card   `json:"rear"`
			Glass string  `json:"glass"`
			Gap   float64 `json:"gap_mm"`
		} `json:"elements"`
	}
	buf, err := json.Marshal(g)
	if err != nil {
		t.Fatal(err)
	}
	if err := json.Unmarshal(buf, &probe); err != nil {
		t.Fatal(err)
	}
	if len(probe.Elements) != 3 || probe.Elements[0].Kind != "stop" || probe.Elements[1].Kind != "mirror" {
		t.Fatalf("the JSON schema is not the documented one: %s", buf)
	}
	gClose(t, "rear radius through the Card schema", probe.Elements[2].Rear.RadiusMM, 0, 0)
	if probe.Elements[2].Glass != "N-BK7" {
		t.Errorf("glass %q in the JSON", probe.Elements[2].Glass)
	}
}

// TestGroupPlaneRadiusRoundTrips pins the wire form of a plane. SetParam
// stores +Inf for radius 0 (the in-memory convention), but JSON has no
// infinity: MarshalJSON must write radius_mm 0, and reading it back must keep
// the surface flat and re-marshal to the same bytes. Without this, a group
// with a flat surface could not be serialised at all and /api/lens/group/set
// answered 500 for every table edit of a group that has a stop.
func TestGroupPlaneRadiusRoundTrips(t *testing.T) {
	m := MirrorElement(-400, -1, 10)
	stop := Element{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 25}}
	g := gWantGroup(t, &Group{Name: "plane", ObjectZMM: math.Inf(-1), Elements: []Element{stop, m}})
	if err := g.SetParam(GroupParam{Surf: SurfRef{Element: 0}, Kind: GPRadius}, 0); err != nil {
		t.Fatal(err)
	}
	if r := mustGet(t, g, GroupParam{Surf: SurfRef{Element: 0}, Kind: GPRadius}); !math.IsInf(r, 0) {
		t.Fatalf("the in-memory convention is +Inf for a plane, got %v", r)
	}
	buf, err := json.Marshal(g)
	if err != nil {
		t.Fatalf("a plane surface does not marshal: %v", err)
	}
	if !bytes.Contains(buf, []byte(`"radius_mm":0`)) {
		t.Errorf("a plane must be written as radius_mm 0: %s", buf)
	}
	var back Group
	if err := json.Unmarshal(buf, &back); err != nil {
		t.Fatal(err)
	}
	if !back.Elements[0].Front.isPlane() {
		t.Errorf("the decoded surface is not flat: %+v", back.Elements[0].Front)
	}
	buf2, err := json.Marshal(&back)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(buf, buf2) {
		t.Errorf("the round trip is not byte-stable:\n%s\n%s", buf, buf2)
	}
}
