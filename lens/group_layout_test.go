package lens

import (
	"encoding/json"
	"fmt"
	"math"
	"reflect"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
)

// ---------------------------------------------------------------------------
// Layout, chain and conversion tests (task group3d)
// ---------------------------------------------------------------------------
//
// Every threshold is justified in the text of the test that uses it. The
// hand-computed cases use the exact geometry of the chain model: the fold of a
// mirror tilted by θ about x turns the axis from +z to (0, sin 2θ, −cos 2θ).

// gClose fails when got differs from want by more than tol (a plain absolute
// comparison: the layout values are exact to the arithmetic of the trig
// functions, ~1e-16, so the tests use 1e-12).
func gClose(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.IsNaN(got) || math.Abs(got-want) > tol {
		t.Errorf("%s = %.15g, want %.15g (±%g)", name, got, want, tol)
	}
}

// gVecClose fails when the two vectors differ by more than tol.
func gVecClose(t *testing.T, name string, got, want geo.Vec3, tol float64) {
	t.Helper()
	if d := got.Sub(want).Norm(); d > tol {
		t.Errorf("%s = %v, want %v (|diff| %.3g > %g)", name, got, want, d, tol)
	}
}

// gWantGroup validates a group or fails the test.
func gWantGroup(t *testing.T, g *Group) *Group {
	t.Helper()
	if err := g.Validate(); err != nil {
		t.Fatalf("the test group is invalid: %v", err)
	}
	return g
}

// rotXVec applies the right-handed active rotation about x to v (the rotation
// of geo.CoordBreak, spelled out here to hand-compute the expectations).
func rotXVec(a float64, v geo.Vec3) geo.Vec3 {
	c, s := math.Cos(a), math.Sin(a)
	return geo.V(c*v.X-s*v.Z, v.Y, s*v.X+c*v.Z)
}

// TestGroupLayoutMatchesVertexZ checks the chain arithmetic against package
// geo: for a coaxial group (no poses, no folds) every layout vertex is the
// vertex System.VertexZ reports for the same surface and every axis is the
// global +z. VertexZ walks the frames of the geo trace, so the two must agree
// to the arithmetic of the chain: 1e-12 mm (they are in fact identical
// integers).
func TestGroupLayoutMatchesVertexZ(t *testing.T) {
	g := gWantGroup(t, &Group{
		Name:      "coaxial trio",
		ObjectZMM: math.Inf(-1),
		OriginZ:   7.5,
		Elements: []Element{
			{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: 5}, GapMM: 20},
			LensElement(103.361, -103.361, 4, "N-BK7", 12.7),
			MirrorElement(-400, -1, 10),
		},
	})
	g.Elements[1].GapMM = 20
	// The gap of the last element: the mirror has no next element, so this is
	// the nominal image distance (not part of the system).
	g.Elements[2].GapMM = 200
	sys, err := g.System()
	if err != nil {
		t.Fatal(err)
	}
	lay, err := g.Layout()
	if err != nil {
		t.Fatal(err)
	}
	if len(lay) != 4 {
		t.Fatalf("layout has %d surfaces, want 4", len(lay))
	}
	want := map[int]float64{ // surface index -> local axial position
		0: 0,         // the stop
		1: 20, 2: 24, // the lens: 20 mm behind the stop, 4 mm thick
		3: 44, // the mirror: 20 mm behind the lens rear
	}
	for _, l := range lay {
		gClose(t, fmt.Sprintf("global z of surface %d", l.Surface), l.Vertex.Z, g.OriginZ+want[l.Surface], 1e-12)
		gVecClose(t, fmt.Sprintf("axis of surface %d", l.Surface), l.Axis, geo.V(0, 0, 1), 1e-12)
		gClose(t, "folds of surface", float64(l.Folds), 0, 0)
		if got := sys.VertexZ(l.Surface); math.Abs(got-l.Vertex.Z) > 1e-12 {
			t.Errorf("surface %d: layout z %.15g, geo VertexZ %.15g", l.Surface, l.Vertex.Z, got)
		}
	}
	// The optional fields: nothing is decentred, so every vertex lies on the
	// axis and the normal is the axis.
	for _, l := range lay {
		gVecClose(t, "vertex on the axis", geo.V(l.Vertex.X, l.Vertex.Y, 0), geo.Vec3{}, 1e-12)
		gVecClose(t, "vertex normal", l.Normal, l.Axis, 1e-12)
	}
	// OriginZ and the axial positions of the elements.
	for i, want := range []float64{7.5, 27.5, 51.5} {
		z, err := g.PositionMM(i)
		if err != nil {
			t.Fatal(err)
		}
		gClose(t, "position of element", z, want, 0)
	}
}

// TestGroupLayoutSingleFold90Deg is the first hand-computed fold: a plane
// mirror tilted by 45 degrees about x turns the beam from +z to +y (the
// reflection of ẑ at the normal (0, −sin 45°, cos 45°)), so the following
// element sits at the mirror vertex plus the gap along +y and its axis is +y.
// Threshold 1e-12: the fold rotation is exact trig, the chain is two
// additions.
func TestGroupLayoutSingleFold90Deg(t *testing.T) {
	mirror := MirrorElement(0, 0, 20)
	mirror.TiltXDeg = 45
	mirror.GapMM = 50
	lens := LensElement(103.361, -103.361, 4, "N-BK7", 12.7)
	lens.GapMM = 0
	g := gWantGroup(t, &Group{Name: "one fold", ObjectZMM: math.Inf(-1), Elements: []Element{mirror, lens}})
	lay, err := g.Layout()
	if err != nil {
		t.Fatal(err)
	}
	if len(lay) != 3 {
		t.Fatalf("layout has %d surfaces, want 3", len(lay))
	}
	// The mirror itself keeps its own (tilted) axis, and starts the chain at
	// the origin.
	gVecClose(t, "mirror vertex", lay[0].Vertex, geo.V(0, 0, 0), 1e-12)
	gVecClose(t, "mirror axis", lay[0].Axis, geo.V(0, -math.Sin(math.Pi/4), math.Cos(math.Pi/4)), 1e-12)
	gClose(t, "mirror folds", float64(lay[0].Folds), 0, 0)
	// The lens sits 50 mm along the folded axis (+y).
	wantAxis := geo.V(0, 1, 0)
	gVecClose(t, "lens front vertex", lay[1].Vertex, geo.V(0, 50, 0), 1e-12)
	gVecClose(t, "lens front axis", lay[1].Axis, wantAxis, 1e-12)
	gVecClose(t, "lens rear vertex", lay[2].Vertex, geo.V(0, 54, 0), 1e-12)
	gVecClose(t, "lens rear axis", lay[2].Axis, wantAxis, 1e-12)
	gClose(t, "lens folds", float64(lay[1].Folds), 1, 0)
	// The exit frame is the vertex of the last surface and the folded axis.
	v, ax, err := g.ExitFrame()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "exit vertex", v, geo.V(0, 54, 0), 1e-12)
	gVecClose(t, "exit axis", ax, wantAxis, 1e-12)
	// geo walks the same geometry: the vertex z of the lens front is 0 (it
	// lies on the folded axis, at y = 50, z = 0).
	if got := sysVertexZ(t, g, 1); math.Abs(got) > 1e-12 {
		t.Errorf("geo VertexZ(1) = %g, want 0", got)
	}
}

// sysVertexZ returns System.VertexZ of a surface of the group.
func sysVertexZ(t *testing.T, g *Group, i int) float64 {
	t.Helper()
	sys, err := g.System()
	if err != nil {
		t.Fatal(err)
	}
	return sys.VertexZ(i)
}

// TestGroupLayoutTwoFolds90Deg is the second hand-computed case: two 45 degree
// mirrors in series. The first turns +z into +y; the second is tilted by 45
// degrees in the already folded frame, so its normal is (0, sin 45°, cos 45°)
// and it turns +y into −z. The third element therefore sits at
// (0, g1, 0) + g2·(0, 0, −1): the beam comes back along −z, parallel to the
// original axis. This is the periscope/reverted layout, and its image point is
// checked against the same hand geometry in the 3D test file. Threshold 1e-12.
func TestGroupLayoutTwoFolds90Deg(t *testing.T) {
	m1 := MirrorElement(0, 0, 20)
	m1.TiltXDeg = 45
	m1.GapMM = 40
	m2 := MirrorElement(0, 0, 20)
	m2.TiltXDeg = 45
	m2.GapMM = 30
	lens := LensElement(103.361, -103.361, 4, "N-BK7", 12.7)
	g := gWantGroup(t, &Group{Name: "two folds", ObjectZMM: math.Inf(-1), Elements: []Element{m1, m2, lens}})
	lay, err := g.Layout()
	if err != nil {
		t.Fatal(err)
	}
	byFace := func(e int, f Face) LayoutSurface {
		for _, l := range lay {
			if l.Element == e && l.Face == f {
				return l
			}
		}
		t.Fatalf("no layout entry for element %d %s", e, f)
		return LayoutSurface{}
	}
	gVecClose(t, "mirror 1 vertex", byFace(0, Front).Vertex, geo.V(0, 0, 0), 1e-12)
	gVecClose(t, "mirror 2 vertex", byFace(1, Front).Vertex, geo.V(0, 40, 0), 1e-12)
	gVecClose(t, "mirror 2 axis", byFace(1, Front).Axis, geo.V(0, math.Sin(math.Pi/4), math.Cos(math.Pi/4)), 1e-12)
	gVecClose(t, "lens vertex after two folds", byFace(2, Front).Vertex, geo.V(0, 40, -30), 1e-12)
	gVecClose(t, "lens axis after two folds", byFace(2, Front).Axis, geo.V(0, 0, -1), 1e-12)
	gVecClose(t, "lens rear after two folds", byFace(2, Rear).Vertex, geo.V(0, 40, -34), 1e-12)
	if n := byFace(2, Front).Folds; n != 2 {
		t.Errorf("folds in front of the lens = %d, want 2", n)
	}
	// The chain is rigid: the distance between the two mirrors is the gap, and
	// the lens sits exactly the second gap away from the second mirror.
	d := byFace(1, Front).Vertex.Sub(byFace(0, Front).Vertex).Norm()
	gClose(t, "mirror separation", d, 40, 1e-12)
	d2 := byFace(2, Front).Vertex.Sub(byFace(1, Front).Vertex).Norm()
	gClose(t, "mirror to lens distance", d2, 30, 1e-12)
}

// TestGroupLayout45DegTurn checks a fold that is not a right angle: a mirror
// tilted by 22.5 degrees turns the beam by 45 degrees, so the next element sits
// at the vertex plus the gap along (0, sin 45°, −cos 45°) — the general
// reflection ẑ → (0, sin 2θ, −cos 2θ) with 2θ = 45°. Threshold 1e-12.
func TestGroupLayout45DegTurn(t *testing.T) {
	m := MirrorElement(0, 0, 20)
	m.TiltXDeg = 22.5
	m.GapMM = 40
	s := StopElement(8)
	g := gWantGroup(t, &Group{Name: "45 degree turn", ObjectZMM: math.Inf(-1), Elements: []Element{m, s}})
	lay, err := g.Layout()
	if err != nil {
		t.Fatal(err)
	}
	// The reflected axis of a mirror tilted by θ about x is
	// (0, sin 2θ, −cos 2θ); here 2θ = 45°.
	want := geo.V(0, math.Sin(math.Pi/4), -math.Cos(math.Pi/4))
	gVecClose(t, "turn axis", lay[1].Axis, want, 1e-12)
	gVecClose(t, "next element vertex", lay[1].Vertex, geo.V(0, 40*math.Sin(math.Pi/4), -40*math.Cos(math.Pi/4)), 1e-12)
	// The deflection from the reversed incident direction is 2θ = 45°, the
	// classical mirror law (geo's pose_test checks the same on a surface).
	rev := geo.V(0, 0, -1)
	ang := math.Atan2(lay[1].Axis.Cross(rev).Norm(), lay[1].Axis.Dot(rev))
	gClose(t, "beam deflection from the reversed axis (rad)", ang, math.Pi/4, 1e-12)
}

// TestGroupSystemStructure checks the surface list of the converted geo.System:
// the element surfaces in beam order with the coordinate breaks between them,
// the media of a cemented interface, the gaps (including the zero gaps that
// close a posed element) and the stop index.
func TestGroupSystemStructure(t *testing.T) {
	doubletFront := Element{
		Name: "crown", Kind: ElementLens,
		Front: SurfaceSpec{RadiusMM: 44.8, SemiAperMM: 10},
		Rear:  &SurfaceSpec{RadiusMM: 0, SemiAperMM: 10},
		Glass: "N-BK7", GapGlass: "F2", ThicknessMM: 5, GapMM: 0,
	}
	doubletRear := SingleSurface(94.8, 10)
	doubletRear.GapGlass = ""
	doubletRear.GapMM = 16.7
	mirror := MirrorElement(0, 0, 20)
	mirror.TiltXDeg = 45
	mirror.GapMM = 40
	g := gWantGroup(t, &Group{
		Name:      "doublet and fold",
		ObjectZMM: math.Inf(-1),
		Elements:  []Element{doubletFront, doubletRear, mirror},
	})
	sys, err := g.System()
	if err != nil {
		t.Fatal(err)
	}
	// 3 surfaces of the doublet, then the frame break, the mirror and the fold
	// break.
	if n := sys.NumSurfaces(); n != 6 {
		t.Fatalf("system has %d surfaces, want 6 (3 doublet + frame + mirror + fold):\n%s",
			n, surfaceSummary(sys))
	}
	if !sys.Surfaces[3].Break || sys.Surfaces[3].Pose == nil {
		t.Errorf("surface 3 must be the element frame break: %+v", sys.Surfaces[3])
	}
	if sys.Surfaces[4].Profile == nil || !sys.Surfaces[4].Mirror {
		t.Errorf("surface 4 must be the mirror: %+v", sys.Surfaces[4])
	}
	if !sys.Surfaces[5].Break || sys.Surfaces[5].Pose == nil {
		t.Errorf("surface 5 must be the fold break: %+v", sys.Surfaces[5])
	}
	// Media: air before the front, BK7 inside the crown, F2 in the gap that
	// the GapGlass of the first element declares, air after the doublet.
	wantMed := []string{"air", "N-BK7", "F2", "air", "air", "air", "air"}
	for i, want := range wantMed {
		if got := sys.Mediums[i].Name(); got != want {
			t.Errorf("medium %d = %q, want %q", i, got, want)
		}
	}
	// Thicknesses: the doublet gaps, then a zero gap to the frame break (the
	// mirror's own gap is carried by the fold break).
	wantThick := []float64{5, 0, 16.7, 0, 0}
	for i, want := range wantThick {
		gClose(t, "thickness", sys.Thickness[i], want, 0)
	}
	if sys.StopIndex != 0 {
		t.Errorf("stop index %d, want 0", sys.StopIndex)
	}
	if sys.OriginZ != 0 || sys.Name != "doublet and fold" {
		t.Errorf("system origin/name %g/%q", sys.OriginZ, sys.Name)
	}
}

// surfaceSummary describes a system for error messages.
func surfaceSummary(sys *geo.System) string {
	var b strings.Builder
	for i, s := range sys.Surfaces {
		kind := "surface"
		if s.Break {
			kind = "break"
		}
		b.WriteString(kind)
		b.WriteString(" ")
		b.WriteString(s.Label())
		if i < len(sys.Surfaces)-1 {
			b.WriteString("; ")
		}
	}
	return b.String()
}

// TestGroupPrescriptionCassegrain checks the axial fold-back: a Cassegrain is
// two mirrors and an aperture stop, and its card form is the classical one —
// the return of the beam is a negative thickness (the library convention, see
// library/geo/mirror-cassegrain-f1000-d300-b300.json, whose numbers are used
// here). The group and the card prescription must produce the same system: the
// same surfaces, media, thicknesses and stop, and the same traced ray.
func TestGroupPrescriptionCassegrain(t *testing.T) {
	stop := StopElement(150)
	stop.Name = "aperture stop"
	stop.GapMM = 60
	primary := MirrorElement(-2000, 0, 153)
	primary.Name = "primary"
	primary.GapMM = -300
	secondary := MirrorElement(-646.1538461538461, 0, 117.5)
	secondary.Name = "secondary"
	secondary.GapMM = 600
	g := gWantGroup(t, &Group{Name: "cassegrain", ObjectZMM: math.Inf(-1), Elements: []Element{stop, primary, secondary}})

	p, err := g.Prescription()
	if err != nil {
		t.Fatal(err)
	}
	if len(p.Cards) != 3 {
		t.Fatalf("prescription has %d cards, want 3", len(p.Cards))
	}
	gClose(t, "card 0 radius", p.Cards[0].RadiusMM, 0, 0)
	gClose(t, "card 0 thickness", p.Cards[0].ThicknessMM, 60, 0)
	gClose(t, "card 0 semi-aperture", p.Cards[0].SemiAperMM, 150, 0)
	gClose(t, "card 1 radius", p.Cards[1].RadiusMM, -2000, 0)
	gClose(t, "card 1 thickness (the return)", p.Cards[1].ThicknessMM, -300, 0)
	if !p.Cards[1].Mirror || !p.Cards[2].Mirror {
		t.Error("the two mirrors must be mirror cards")
	}
	gClose(t, "card 2 thickness (the image distance)", p.Cards[2].ThicknessMM, 600, 0)
	if p.Stop != 0 || p.ObjectZMM > -1e300 {
		t.Errorf("stop %d, object %v", p.Stop, p.ObjectZMM)
	}
	// The two systems describe the same surfaces.
	sysA, err := g.System()
	if err != nil {
		t.Fatal(err)
	}
	sysB, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	if err := sameGeometry(sysA, sysB); err != nil {
		t.Errorf("the group and the prescription disagree: %v", err)
	}
	// The same ray: an axial ray, reflected by the primary and the secondary.
	r := geo.NewRay(geo.V(0, 20, -10), geo.V(0, 0, 1), geo.Wavelengthd)
	ra, rb := sysA.Trace(r), sysB.Trace(r)
	if !ra.OK() || !rb.OK() {
		t.Fatalf("trace status %v / %v", ra.Status, rb.Status)
	}
	gVecClose(t, "traced final point", ra.Ray.Origin, rb.Ray.Origin, 1e-12)
	gVecClose(t, "traced final direction", ra.Ray.Dir, rb.Ray.Dir, 1e-12)
	gClose(t, "optical path length", ra.Ray.OPL, rb.Ray.OPL, 1e-12)
}

// sameGeometry compares two systems surface by surface: count, profiles,
// mirrors, apertures, thicknesses and medium names.
func sameGeometry(a, b *geo.System) error {
	// The coordinate-break elements and the pose of a surface are a property of
	// the description, not of the geometry: compare the optical surfaces.
	realSurfaces := func(s *geo.System) []geo.Surface {
		out := make([]geo.Surface, 0, len(s.Surfaces))
		for _, sf := range s.Surfaces {
			if sf.Break {
				continue
			}
			out = append(out, sf)
		}
		return out
	}
	sa, sb := realSurfaces(a), realSurfaces(b)
	if len(sa) != len(sb) {
		return &geometryError{"surface count", len(sa), len(sb)}
	}
	for i := range sa {
		x, y := sa[i], sb[i]
		if (x.Profile == nil) != (y.Profile == nil) {
			return &geometryError{"nil profile", i, i}
		}
		if x.Profile != nil && x.Profile.String() != y.Profile.String() {
			return &geometryError{"profile " + x.Profile.String(), i, y.Profile.String()}
		}
		if x.Mirror != y.Mirror {
			return &geometryError{"mirror flag", i, i}
		}
		if x.MaxRadius() != y.MaxRadius() {
			return &geometryError{"aperture radius", x.MaxRadius(), y.MaxRadius()}
		}
	}
	// The gaps between the optical surfaces (the zero gaps that close a posed
	// element are bookkeeping of the chain, not geometry).
	gaps := func(s *geo.System) []float64 {
		out := []float64{}
		cur := 0.0
		first := true
		for i, sf := range s.Surfaces {
			if sf.Break {
				continue
			}
			if !first {
				out = append(out, cur)
			}
			first = false
			cur = 0
			if i < len(s.Thickness) {
				// accumulate until the next real surface
				for k := i; k < len(s.Thickness); k++ {
					cur += s.Thickness[k]
					if k+1 < len(s.Surfaces) && !s.Surfaces[k+1].Break {
						break
					}
				}
			}
		}
		return out
	}
	ga, gb := gaps(a), gaps(b)
	if len(ga) != len(gb) {
		return &geometryError{"gap count", len(ga), len(gb)}
	}
	for i := range ga {
		if math.Abs(ga[i]-gb[i]) > 1e-12 {
			return &geometryError{"gap", ga[i], gb[i]}
		}
	}
	// The stop is a surface: compare its position among the optical surfaces
	// (the break entries of a posed element shift the index).
	stopOrdinal := func(s *geo.System) int {
		n := 0
		for i := 0; i < s.StopIndex && i < len(s.Surfaces); i++ {
			if !s.Surfaces[i].Break {
				n++
			}
		}
		return n
	}
	if stopOrdinal(a) != stopOrdinal(b) {
		return &geometryError{"stop (optical surface ordinal)", stopOrdinal(a), stopOrdinal(b)}
	}
	// The media in front of and behind every optical surface.
	med := func(s *geo.System) []string {
		out := []string{}
		for i, sf := range s.Surfaces {
			if sf.Break {
				continue
			}
			if i < len(s.Mediums) {
				out = append(out, s.Mediums[i].Name())
			}
			if i+1 < len(s.Mediums) {
				out = append(out, s.Mediums[i+1].Name())
			}
		}
		return out
	}
	ma, mb := med(a), med(b)
	if len(ma) != len(mb) {
		return &geometryError{"medium count", len(ma), len(mb)}
	}
	for i := range ma {
		if ma[i] != mb[i] {
			return &geometryError{"medium", ma[i], mb[i]}
		}
	}
	return nil
}

// geometryError is the mismatch reported by sameGeometry.
type geometryError struct {
	what string
	a, b interface{}
}

func (e *geometryError) Error() string {
	return fmt.Sprintf("%s: %v vs %v", e.what, e.a, e.b)
}

// TestGroupPrescriptionRejects3DFold checks the documented limits of the card
// form: an out-of-plane fold with elements after it, and a tilted two-surface
// element, are refused with an error that names the element; the same fold as
// the last element is accepted (nothing follows it).
func TestGroupPrescriptionRejects3DFold(t *testing.T) {
	m := MirrorElement(0, 0, 20)
	m.TiltXDeg = 45
	m.GapMM = 40
	g := gWantGroup(t, &Group{Name: "fold", ObjectZMM: math.Inf(-1), Elements: []Element{m, StopElement(8)}})
	if _, err := g.Prescription(); err == nil {
		t.Error("a fold followed by an element must not convert to a card prescription")
	} else if !strings.Contains(err.Error(), "element 0") {
		t.Errorf("the error must name the folding element: %v", err)
	}
	// The same mirror as the last element converts: the beam simply leaves.
	g2 := gWantGroup(t, &Group{Name: "fold last", ObjectZMM: math.Inf(-1), Elements: []Element{StopElement(8), m}})
	if _, err := g2.Prescription(); err != nil {
		t.Errorf("a fold at the last element: %v", err)
	}
	// A tilted two-surface element is refused for the same reason (the card
	// model would measure its thickness along the chain axis).
	lens := LensElement(103.361, -103.361, 4, "N-BK7", 12.7)
	lens.TiltXDeg = 3
	g3 := gWantGroup(t, &Group{Name: "tilted lens", ObjectZMM: math.Inf(-1), Elements: []Element{lens}})
	if _, err := g3.Prescription(); err == nil {
		t.Error("a tilted two-surface element must not convert to cards")
	}
	// A decentred one converts: the two constructions agree exactly, because a
	// decentre does not rotate the element axis.
	lens2 := LensElement(103.361, -103.361, 4, "N-BK7", 12.7)
	lens2.DecenterXMM, lens2.DecenterYMM = 0.5, -0.25
	g4 := gWantGroup(t, &Group{Name: "decentred lens", ObjectZMM: math.Inf(-1), Elements: []Element{lens2}})
	p4, err := g4.Prescription()
	if err != nil {
		t.Fatalf("a decentred lens must convert: %v", err)
	}
	gClose(t, "decentred card x", p4.Cards[0].DecenterXMM, 0.5, 0)
	sysA, err := g4.System()
	if err != nil {
		t.Fatal(err)
	}
	sysB, err := p4.System()
	if err != nil {
		t.Fatal(err)
	}
	if err := sameGeometry(sysA, sysB); err != nil {
		t.Errorf("decentred lens: %v", err)
	}
}

// TestPrescriptionGroupRoundTrip converts a hand-written cemented doublet
// prescription (the layout of the library achromats) into elements and back,
// and checks that both systems describe the same geometry and trace the same
// rays. The doublet is [R1 R2 R3] with the glasses N-BK7 / F2, i.e. a lens
// [R1 R2] whose GapGlass is F2 at gap 0 followed by the single surface R3.
func TestPrescriptionGroupRoundTrip(t *testing.T) {
	p := &Prescription{
		Name:      "cemented doublet",
		ObjectZMM: math.Inf(-1),
		Stop:      0,
		Cards: []Card{
			{Name: "crown front", RadiusMM: 44.8, ThicknessMM: 5, Glass: "N-BK7", SemiAperMM: 10},
			{Name: "cement", RadiusMM: 0, ThicknessMM: 3, Glass: "F2", SemiAperMM: 10},
			{Name: "flint rear", RadiusMM: 94.8, ThicknessMM: 16.7, SemiAperMM: 10},
		},
	}
	g, err := p.Group()
	if err != nil {
		t.Fatal(err)
	}
	if len(g.Elements) != 2 {
		t.Fatalf("the doublet must give 2 elements, got %d", len(g.Elements))
	}
	e0, e1 := g.Elements[0], g.Elements[1]
	if e0.Kind != ElementLens || e0.Rear == nil || e0.Glass != "N-BK7" || e0.GapGlass != "F2" {
		t.Errorf("first element %+v", e0)
	}
	gClose(t, "crown thickness", e0.ThicknessMM, 5, 0)
	gClose(t, "cement gap", e0.GapMM, 3, 0)
	if e1.Kind != ElementLens || e1.Rear != nil {
		t.Errorf("second element must be the single flint surface: %+v", e1)
	}
	gClose(t, "single surface radius", e1.Front.RadiusMM, 94.8, 0)
	gClose(t, "flint gap", e1.GapMM, 16.7, 0)
	if g.Stop != 0 {
		t.Errorf("stop %d, want 0", g.Stop)
	}
	// The geometry is preserved.
	sysA, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	sysB, err := g.System()
	if err != nil {
		t.Fatal(err)
	}
	if err := sameGeometry(sysA, sysB); err != nil {
		t.Errorf("round trip: %v", err)
	}
	// … and so is the trace (three surfaces, same media).
	for _, y := range []float64{0, 2.5, -4} {
		r := geo.NewRay(geo.V(0.3, y, -20), geo.V(0, 0, 1), geo.Wavelengthd)
		ra, rb := sysA.Trace(r), sysB.Trace(r)
		if !ra.OK() || !rb.OK() {
			t.Fatalf("trace status %v / %v for y=%g", ra.Status, rb.Status, y)
		}
		gVecClose(t, "traced point", ra.Ray.Origin, rb.Ray.Origin, 1e-12)
		gClose(t, "traced OPL", ra.Ray.OPL, rb.Ray.OPL, 1e-12)
	}
	// Back to cards: the same three cards, field by field (an element carries
	// the name of the first card of its pair, so only that name survives).
	p2, err := g.Prescription()
	if err != nil {
		t.Fatal(err)
	}
	strip := func(cs []Card) []Card {
		out := make([]Card, len(cs))
		for i, c := range cs {
			c.Name = ""
			out[i] = c
		}
		return out
	}
	if !reflect.DeepEqual(strip(p.Cards), strip(p2.Cards)) {
		t.Errorf("cards changed in the round trip:\n%+v\n%+v", p.Cards, p2.Cards)
	}
	if p2.Cards[0].Name != "crown front" {
		t.Errorf("the front card name was lost: %q", p2.Cards[0].Name)
	}
}

// TestGroupJSONByteStable checks the JSON round trip of a group with every
// feature (a lens with an asphere and an aperture, a cemented interface, a
// folded mirror, a stop): marshal → unmarshal → marshal is byte-identical, the
// loaded group is deep-equal to the original, and the object position survives
// as "inf".
func TestGroupJSONByteStable(t *testing.T) {
	lens := Element{
		Name: "asphere", Kind: ElementLens,
		Front: SurfaceSpec{RadiusMM: 50, ConicK: -0.6, Coeffs: []float64{-1e-6, 2e-9}, SemiAperMM: 12},
		Rear:  &SurfaceSpec{RadiusMM: -80, SemiAperMM: 12},
		Glass: "N-BK7", GapGlass: "F2", ThicknessMM: 6, GapMM: 0,
		DecenterXMM: 0.25, TiltZDeg: 30,
	}
	rear := SingleSurface(0, 12)
	rear.GapGlass = ""
	rear.GapMM = 25
	m := MirrorElement(0, 0, 15)
	m.TiltXDeg = 45
	m.Name = "fold"
	m.GapMM = 30
	stop := StopElement(8)
	stop.GapMM = 10
	g := gWantGroup(t, &Group{Name: "json", ObjectZMM: math.Inf(-1), Elements: []Element{stop, lens, rear, m}})

	buf, err := json.Marshal(g)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(buf), `"object_z_mm":"inf"`) {
		t.Errorf("the object at infinity must be written as \"inf\": %s", buf)
	}
	var back Group
	if err := json.Unmarshal(buf, &back); err != nil {
		t.Fatal(err)
	}
	buf2, err := json.Marshal(&back)
	if err != nil {
		t.Fatal(err)
	}
	if string(buf) != string(buf2) {
		t.Errorf("marshal → unmarshal → marshal is not byte-stable:\n%s\n%s", buf, buf2)
	}
	if !reflect.DeepEqual(g.Elements, back.Elements) || back.Name != g.Name || back.Stop != g.Stop || back.OriginZ != g.OriginZ {
		t.Errorf("the loaded group differs:\n%+v\n%+v", g, back)
	}
	if !math.IsInf(back.ObjectZMM, -1) {
		t.Errorf("object position %v, want -Inf", back.ObjectZMM)
	}
	// The same group also survives the strong check of the deep clone.
	if !reflect.DeepEqual(g.Elements, g.Clone().Elements) {
		t.Error("Clone lost element data")
	}
	// A finite object position and a non-zero origin are written as numbers.
	g2 := g.Clone()
	g2.ObjectZMM, g2.OriginZ = -200, 12.5
	buf, err = json.Marshal(g2)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(buf), `"object_z_mm":-200`) || !strings.Contains(string(buf), `"origin_z_mm":12.5`) {
		t.Errorf("finite object/origin not written: %s", buf)
	}
	var back2 Group
	if err := json.Unmarshal(buf, &back2); err != nil {
		t.Fatal(err)
	}
	gClose(t, "object position", back2.ObjectZMM, -200, 0)
	// A hand-written group with a missing object position loads as an object at
	// infinity, and a bad value is refused.
	var hand Group
	if err := json.Unmarshal([]byte(`{"name":"hand","stop":0,"elements":[{"kind":"stop","front":{"radius_mm":0,"semi_aperture_mm":5},"gap_mm":0}]}`), &hand); err != nil {
		t.Fatal(err)
	}
	if !math.IsInf(hand.ObjectZMM, -1) || hand.Len() != 1 {
		t.Errorf("hand-written group: object %v, %d elements", hand.ObjectZMM, hand.Len())
	}
	if err := json.Unmarshal([]byte(`{"object_z_mm":"soon","elements":[]}`), &hand); err == nil {
		t.Error("a nonsense object position must be refused")
	}
}

// TestGroupPosedElementIsRigidAndLeavesTheChain checks the two claims the
// element model makes about a posed element: it is rigid (its surfaces move
// together, the thickness is measured along the element's own axis) and it does
// not move the following elements (the chain frame is restored, origin
// included, as the auto-return of a posed Card does). A lens tilted by 10
// degrees about y therefore has its rear vertex at t·(sin 10°, 0, cos 10°) —
// the hand-computed rigid displacement — while the stop behind it stays on the
// axis at the nominal gap. Threshold 1e-12 mm.
func TestGroupPosedElementIsRigidAndLeavesTheChain(t *testing.T) {
	const (
		thk = 5.0
		gap = 20.0
		ty  = 10.0
		dx  = 2.0
	)
	th := ty * math.Pi / 180
	lens := LensElement(103.361, -103.361, thk, "N-BK7", 12)
	lens.TiltYDeg = ty
	lens.DecenterXMM = dx
	lens.GapMM = gap
	g := gWantGroup(t, &Group{
		Name: "posed", ObjectZMM: math.Inf(-1),
		Elements: []Element{lens, StopElement(5)},
	})
	lay, err := g.Layout()
	if err != nil {
		t.Fatal(err)
	}
	if len(lay) != 3 {
		t.Fatalf("layout has %d surfaces, want 3", len(lay))
	}
	// Front: the decentre only (the tilt rotates about the decentred vertex).
	gVecClose(t, "posed front vertex", lay[0].Vertex, geo.V(dx, 0, 0), 1e-12)
	wantAxis := geo.V(math.Sin(th), 0, math.Cos(th))
	gVecClose(t, "posed front axis", lay[0].Axis, wantAxis, 1e-12)
	// Rear: the rigid translation along the element's own axis.
	gVecClose(t, "posed rear vertex", lay[1].Vertex, geo.V(dx, 0, 0).Add(wantAxis.Mul(thk)), 1e-12)
	gVecClose(t, "posed rear axis", lay[1].Axis, wantAxis, 1e-12)
	// The following element is back on the chain axis at the nominal gap.
	gVecClose(t, "the following element is not moved", lay[2].Vertex, geo.V(0, 0, gap), 1e-12)
	gVecClose(t, "the following element keeps the chain axis", lay[2].Axis, geo.V(0, 0, 1), 1e-12)
	if n := lay[2].Folds; n != 0 {
		t.Errorf("a tilted lens must not fold the chain: %d folds", n)
	}
	// PositionMM agrees with the layout for posed and unposed elements alike
	// (a posed element does not advance the chain by its thickness).
	z, err := g.PositionMM(1)
	if err != nil {
		t.Fatal(err)
	}
	gClose(t, "chain position of the following element", z, gap, 1e-12)
	unposed := g.Clone()
	if err := unposed.SetPose(0, 0, 0, 0, 0, 0); err != nil {
		t.Fatal(err)
	}
	z, err = unposed.PositionMM(1)
	if err != nil {
		t.Fatal(err)
	}
	gClose(t, "chain position with an unposed element", z, gap+thk, 1e-12)
	// The traced system follows the same layout: the ray hits the lens at the
	// decentred vertex and the stop on the axis.
	sys, err := g.System()
	if err != nil {
		t.Fatal(err)
	}
	res := sys.Trace(geo.NewRay(geo.V(dx, 0, -20), geo.V(0, 0, 1), geo.Wavelengthd))
	if !res.OK() || len(res.Steps) != 3 {
		t.Fatalf("trace %v with %d steps", res.Status, len(res.Steps))
	}
	gVecClose(t, "traced front hit", res.Steps[0].Point, geo.V(dx, 0, 0), 1e-12)
	if z := res.Steps[2].Point.Z; math.Abs(z-gap) > 1e-9 {
		t.Errorf("the ray hits the stop at z = %g, want %g", z, gap)
	}
}

// TestGroupFoldOffKeepsTheChain checks the FoldMode switch: with FoldOff a
// tilted mirror steers the beam out of the sequence (the chain keeps the
// incoming axis, so the following elements stay on it and the deflected beam
// misses them), while FoldAuto lays them along the reflected beam. The test
// pins the layout of both variants and the trace status of each.
func TestGroupFoldOffKeepsTheChain(t *testing.T) {
	mk := func(fold FoldMode) *Group {
		m := MirrorElement(0, 0, 20)
		m.TiltXDeg = 45
		m.Fold = fold
		m.GapMM = 30
		return gWantGroup(t, &Group{
			Name: "fold mode", ObjectZMM: math.Inf(-1),
			Elements: []Element{
				func() Element { e := StopElement(6); e.GapMM = 15; return e }(),
				m,
				planoElement(t, 100, 5, 8, "N-BK7"),
			},
		})
	}
	off := mk(FoldOff)
	layOff, err := off.Layout()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "FoldOff: the lens stays on the incoming axis", layOff[2].Vertex, geo.V(0, 0, 45), 1e-12)
	gVecClose(t, "FoldOff: the lens keeps the incoming axis", layOff[2].Axis, geo.V(0, 0, 1), 1e-12)
	if n := layOff[2].Folds; n != 0 {
		t.Errorf("FoldOff: %d folds, want 0", n)
	}
	// The exit frame is the frame of the last surface of the group. With
	// FoldOff the beam leaves the sequence at the mirror and never reaches that
	// surface, so the frame describes the unused path: the lens axis, +z.
	if _, ax, err := off.ExitFrame(); err != nil {
		t.Fatal(err)
	} else {
		gVecClose(t, "FoldOff: the exit frame is the axis of the last element", ax, geo.V(0, 0, 1), 1e-12)
	}
	sys, err := off.System()
	if err != nil {
		t.Fatal(err)
	}
	if r := sys.Trace(geo.NewRay(geo.V(0, 0, -20), geo.V(0, 0, 1), geo.Wavelengthd)); r.OK() {
		t.Errorf("FoldOff: the deflected beam must miss the lens behind it, got %v", r.Status)
	}
	on := mk(FoldAuto)
	layOn, err := on.Layout()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "FoldAuto: the lens follows the folded axis", layOn[2].Vertex, geo.V(0, 30, 15), 1e-12)
	gVecClose(t, "FoldAuto: the folded axis is +y", layOn[2].Axis, geo.V(0, 1, 0), 1e-12)
	if n := layOn[2].Folds; n != 1 {
		t.Errorf("FoldAuto: %d folds, want 1", n)
	}
	sysOn, err := on.System()
	if err != nil {
		t.Fatal(err)
	}
	if r := sysOn.Trace(geo.NewRay(geo.V(0, 0, -20), geo.V(0, 0, 1), geo.Wavelengthd)); !r.OK() {
		t.Errorf("FoldAuto: the beam must traverse the folded sequence, got %v", r.Status)
	}
}
