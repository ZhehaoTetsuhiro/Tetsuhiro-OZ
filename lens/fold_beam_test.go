package lens

// Beam bookkeeping of the sequential chain in return paths (the -z beam of a
// Cassegrain-style two-mirror system):
//
//   - a fold behind a mirror that returns the beam must fold the BEAM, which is
//     antiparallel to the chain axis there, not the chain axis: before the fix
//     the exit axis came out pointing backwards and the focus search ran on the
//     wrong plane through the image;
//   - a fold must reset the beam state for the following elements: before the
//     fix a mirror behind a fold reflected the stale direction of the fold
//     surface (the fold-then-return chain below exposed it as a 45° exit axis);
//   - a refracting surface must keep the propagation sense of the incoming
//     beam: before the fix a corrector in a return path silently flipped the
//     beam to +z and the mirror behind it sent the exit axis backwards.

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
)

// beamTestBundle is the on-axis hexapolar bundle the checks below share.
func beamTestBundle() Bundle {
	return Bundle{Object: geo.FieldTan(0, 0), Samples: geo.Hexapolar(3), LambdaNm: geo.Wavelengthd}
}

// TestGroupFoldFoldsTheBeamBehindAReturnMirror folds a converging return beam:
// a parabola of f = 2000 (R = -4000) returns the beam at -z and the 45 degree
// fold 900 mm into the return path sends it to -y, still converging to a point
// 1100 mm along -y from the fold vertex at (0, 0, -885). The exit axis and the
// image search axis are the folded beam (-y); the values are exact to the
// arithmetic (the parabola is perfect on axis).
func TestGroupFoldFoldsTheBeamBehindAReturnMirror(t *testing.T) {
	g := gWantGroup(t, &Group{
		Name: "return then fold", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(20); e.GapMM = 15; return e }(),
			func() Element { e := MirrorElement(-4000, -1, 0); e.GapMM = -900; return e }(),
			func() Element { e := MirrorElement(0, 0, 0); e.TiltXDeg = 45; e.GapMM = 20; return e }(),
		},
	})
	v, ax, err := g.ExitFrame()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "the fold vertex sits 900 mm into the return path", v, geo.V(0, 0, -885), 1e-9)
	gVecClose(t, "the exit axis is the folded beam (-y), not the chain axis", ax, geo.V(0, -1, 0), 1e-12)
	res, err := g.ImagePosition(beamTestBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "the image search runs along the folded beam", res.Axis, geo.V(0, -1, 0), 1e-12)
	gClose(t, "the focus lies 1100 mm along the folded beam", res.FocusMM, 1100, 1e-6)
	gVecClose(t, "the image point", res.Point, geo.V(0, -1100, -885), 1e-6)
	if res.Failed != 0 {
		t.Errorf("%d of the sampled rays were blocked", res.Failed)
	}
	if res.RMS > 1e-9 {
		t.Errorf("the on-axis spot RMS is %.3g mm, want < 1e-9 (a perfect parabola on axis)", res.RMS)
	}
}

// TestGroupFoldThenReturnMirror folds first and returns behind: the 45 degree
// fold lays the parabola (R = -400) on the folded axis at (0, 20, 15), the
// parabola returns the beam at -y and images at (0, -180, 15), 200 mm from the
// fold vertex. Before the fix the mirror behind the fold reflected the stale
// beam of the fold surface and the exit axis came out (0, 0.707, -0.707).
func TestGroupFoldThenReturnMirror(t *testing.T) {
	g := gWantGroup(t, &Group{
		Name: "fold then return", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(20); e.GapMM = 15; return e }(),
			func() Element { e := MirrorElement(0, 0, 0); e.TiltXDeg = 45; e.GapMM = 20; return e }(),
			MirrorElement(-400, -1, 0),
		},
	})
	lay, err := g.Layout()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "the return mirror lies on the folded axis", lay[2].Vertex, geo.V(0, 20, 15), 1e-9)
	if lay[2].Folds != 1 {
		t.Errorf("%d folds in front of the return mirror, want 1", lay[2].Folds)
	}
	_, ax, err := g.ExitFrame()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "the exit axis is the returned beam (-y)", ax, geo.V(0, -1, 0), 1e-12)
	res, err := g.ImagePosition(beamTestBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "the image search runs along the returned beam", res.Axis, geo.V(0, -1, 0), 1e-12)
	gClose(t, "the focus lies 200 mm along the returned beam", res.FocusMM, 200, 1e-6)
	gVecClose(t, "the image point", res.Point, geo.V(0, -180, 15), 1e-6)
	if res.RMS > 1e-9 {
		t.Errorf("the on-axis spot RMS is %.3g mm, want < 1e-9 (a perfect parabola on axis)", res.RMS)
	}
}

// TestGroupRefractorKeepsTheBeamSenseInAReturnPath puts a zero-power plate
// (focal infinity) in the return beam of a collecting mirror: the beam keeps
// its -z sense through the refracting surface, so the exit axis is -z and the
// image of the parabola (f = 200) forms 50 mm behind the stop at (0, 0, -50),
// 150 mm along -z from the plate at (0, 0, 100). Before the fix the plate
// flipped the beam sense to +z and the exit axis (and the image search) came
// out reversed.
func TestGroupRefractorKeepsTheBeamSenseInAReturnPath(t *testing.T) {
	g := gWantGroup(t, &Group{
		Name: "plate in return path", ObjectZMM: math.Inf(-1),
		Elements: []Element{
			func() Element { e := StopElement(20); e.GapMM = 150; return e }(),
			func() Element { e := MirrorElement(-400, -1, 0); e.GapMM = -50; return e }(),
			Element{Name: "plate", Kind: ElementLens, Front: SurfaceSpec{RadiusMM: 0}},
		},
	})
	v, ax, err := g.ExitFrame()
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "the plate sits 50 mm into the return path", v, geo.V(0, 0, 100), 1e-9)
	gVecClose(t, "the exit axis keeps the -z sense of the return beam", ax, geo.V(0, 0, -1), 1e-12)
	res, err := g.ImagePosition(beamTestBundle(), ImageOptions{})
	if err != nil {
		t.Fatal(err)
	}
	gVecClose(t, "the image search runs along the return beam", res.Axis, geo.V(0, 0, -1), 1e-12)
	gClose(t, "the focus lies 150 mm along -z from the plate", res.FocusMM, 150, 1e-6)
	gVecClose(t, "the image point", res.Point, geo.V(0, 0, -50), 1e-6)
	if res.Failed != 0 {
		t.Errorf("%d of the sampled rays were blocked", res.Failed)
	}
	if res.RMS > 1e-9 {
		t.Errorf("the on-axis spot RMS is %.3g mm, want < 1e-9 (a flat plate does not bend the on-axis bundle)", res.RMS)
	}
}
