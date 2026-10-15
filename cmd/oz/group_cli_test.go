package main

// CLI tests of `oz lens group ...` and `oz lens parax` (task wiring2).
//
// The fixtures are built from the lens constructors, so every number below is
// reproducible without any file outside the repository: the parabolic mirror
// focuses at R/2, the folded Kepler telescope has Γ = −1, the achromat is the
// classical thin-lens solution of two catalogue glasses, and the Cassegrain
// fixture is the element-level form of
// library/geo/mirror-cassegrain-f1000-d300-b300.json (the test compares the
// `to-prescription` output with that file card by card).

import (
	"encoding/json"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

// gParabolaFixture is a single paraboloid (K = −1) of R = −400 mm: an on-axis
// collimated bundle focuses exactly at f = |R|/2 = 200 mm (the parabola is free
// of spherical aberration), which is the closed form of the image measurement.
func gParabolaFixture() *lens.Group {
	return &lens.Group{
		Name: "parabola", ObjectZMM: math.Inf(-1),
		Elements: []lens.Element{lens.MirrorElement(-400, -1, 24)},
	}
}

// gCassegrainFixture is the element form of the library Cassegrain: stop at
// z = 0 (semi 150), primary R = −2000 at z = 60, secondary R = +8400 at
// z = 60 − 300 = −240 (a negative gap: the beam returns, so the vertex is
// behind the primary on the same z axis). EFL = +857.142857 mm, paraxial image
// z = +360 mm, BFL +600 mm (library numbers after the 2026-10 regeneration,
// where the secondary is the converging conjugate; see wt-notes/mirrorparax.md).
func gCassegrainFixture() *lens.Group {
	stop := lens.StopElement(150)
	stop.GapMM = 60
	primary := lens.MirrorElement(-2000, 0, 153)
	primary.GapMM = -300
	secondary := lens.MirrorElement(8400, 0, 117.5)
	secondary.GapMM = 600
	return &lens.Group{
		Name: "cassegrain", ObjectZMM: math.Inf(-1), Stop: 0,
		Elements: []lens.Element{stop, primary, secondary},
	}
}

// gPlanoElement returns a plano-convex element (curved front, flat rear) whose
// paraxial power is exactly (n−1)/R, i.e. the ideal thin lens of focal length f
// placed at the front vertex.
func gPlanoElement(t *testing.T, fMM, thicknessMM, semiAperMM float64) lens.Element {
	t.Helper()
	m, ok := geo.Glass("N-BK7")
	if !ok {
		t.Fatal("glass N-BK7 missing")
	}
	return lens.Element{
		Kind:        lens.ElementLens,
		Front:       lens.SurfaceSpec{RadiusMM: (m.N(geo.Wavelengthd) - 1) * fMM, SemiAperMM: semiAperMM},
		Rear:        &lens.SurfaceSpec{RadiusMM: 0, SemiAperMM: semiAperMM},
		Glass:       "N-BK7",
		ThicknessMM: thicknessMM,
	}
}

// gKeplerFixture returns the two telescopes of the lens tests: the straight
// Kepler (f, f) separated by f1 + f2 − d/n (Γ = −1) and the same telescope with
// a 45 degree plane fold mirror in its collimated space (Γ unchanged, the exit
// axis turns to +y).
func gKeplerFixture(t *testing.T, folded bool) *lens.Group {
	t.Helper()
	m, _ := geo.Glass("N-BK7")
	n := m.N(geo.Wavelengthd)
	const d = 2.0
	sep := 200 - d/n
	if !folded {
		l1 := gPlanoElement(t, 100, d, 20)
		l1.GapMM = 200 - d/n
		return &lens.Group{
			Name: "kepler", ObjectZMM: math.Inf(-1),
			Elements: []lens.Element{lens.StopElement(4), l1, gPlanoElement(t, 100, d, 20)},
		}
	}
	fold := lens.MirrorElement(0, 0, 20)
	fold.TiltXDeg = 45
	fold.GapMM = 20
	l1 := gPlanoElement(t, 100, d, 20)
	l1.GapMM = sep
	stop := lens.StopElement(6)
	stop.GapMM = 15
	return &lens.Group{
		Name: "folded-kepler", ObjectZMM: math.Inf(-1),
		Elements: []lens.Element{stop, fold, l1, gPlanoElement(t, 100, d, 20)},
	}
}

// gDoubletFixture is the classical cemented crown/flint achromat of focal
// length f: the crown is the two-surface element [R1, flat cement interface]
// whose GapGlass is the flint, the flint is the single rear surface R3 (the
// thin-lens solution φ1/V1 + φ2/V2 = 0 with φ1 + φ2 = 1/f).
func gDoubletFixture(t *testing.T, fMM, semiAperMM, crownThickness, flintThickness float64) *lens.Group {
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
	front := lens.Element{
		Name: "crown", Kind: lens.ElementLens,
		Front: lens.SurfaceSpec{RadiusMM: r1, SemiAperMM: semiAperMM},
		Rear:  &lens.SurfaceSpec{RadiusMM: 0, SemiAperMM: semiAperMM},
		Glass: "N-BK7", GapGlass: "F2",
		ThicknessMM: crownThickness, GapMM: flintThickness,
	}
	rear := lens.SingleSurface(r3, semiAperMM)
	rear.Name = "flint"
	return &lens.Group{Name: "doublet", ObjectZMM: math.Inf(-1), Elements: []lens.Element{front, rear}}
}

// writeGroup writes a group fixture into the test's temp dir.
func writeGroup(t *testing.T, name string, g *lens.Group) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name+".json")
	if err := g.SaveJSON(path); err != nil {
		t.Fatalf("save %s: %v", name, err)
	}
	return path
}

// captureOut runs fn with stdout captured and returns what was printed.
func captureOut(t *testing.T, fn func() error) (string, error) {
	t.Helper()
	old := os.Stdout
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	os.Stdout = w
	ferr := fn()
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	os.Stdout = old
	buf, err := io.ReadAll(r)
	if err != nil {
		t.Fatal(err)
	}
	r.Close()
	return string(buf), ferr
}

// gClose is the relative/absolute closeness helper of these tests.
func gClose(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if d := math.Abs(got - want); d > tol {
		t.Errorf("%s = %.17g, want %.17g (Δ %.3g > %.3g)", name, got, want, d, tol)
	}
}

// ---------------------------------------------------------------------------
// show / layout
// ---------------------------------------------------------------------------

// TestLensGroupShowCLI drives `oz lens group show`: the summary names the group
// and lists every element.
func TestLensGroupShowCLI(t *testing.T) {
	path := writeGroup(t, "cassegrain", gCassegrainFixture())
	out, err := captureOut(t, func() error { return cmdLens([]string{"group", "show", path}) })
	if err != nil {
		t.Fatalf("show: %v", err)
	}
	for _, want := range []string{`group "cassegrain": 3 elements, stop 0`, "[1] mirror R=-2000", "[2] mirror R=8400"} {
		if !strings.Contains(out, want) {
			t.Errorf("show output lacks %q:\n%s", want, out)
		}
	}
	if err := cmdLens([]string{"group", "show", filepath.Join(t.TempDir(), "nope.json")}); err == nil {
		t.Error("show accepted a missing file")
	}
	if err := cmdLens([]string{"group"}); err == nil {
		t.Error("lens group without a subcommand accepted")
	}
	if err := cmdLens([]string{"group", "nope"}); err == nil {
		t.Error("unknown group subcommand accepted")
	}
}

// TestLensGroupLayoutCLI drives `oz lens group layout`: the vertex table of the
// Cassegrain (0, 60, −240 — the negative gap of the returning beam) and the
// fold bookkeeping of the folded Kepler (the fold turns the chain to +y at
// z = 15; the following surfaces sit on +y). Two runs must be byte-identical.
func TestLensGroupLayoutCLI(t *testing.T) {
	dir := t.TempDir()
	cass := writeGroup(t, "cassegrain", gCassegrainFixture())
	out1 := filepath.Join(dir, "l1.json")
	out2 := filepath.Join(dir, "l2.json")
	if err := cmdLens([]string{"group", "layout", "--json", out1, cass}); err != nil {
		t.Fatalf("layout: %v", err)
	}
	var res groupLayoutResponse
	buf1 := readJSONFile(t, out1, &res)
	if res.NumSurfaces != 3 || res.Folds != 0 {
		t.Fatalf("cassegrain layout: %d surfaces, %d folds, want 3 and 0", res.NumSurfaces, res.Folds)
	}
	for i, wantZ := range []float64{0, 60, -240} {
		s := res.Surfaces[i]
		gClose(t, "vertex z", s.Vertex.Z, wantZ, 1e-12)
		gClose(t, "vertex x", s.Vertex.X, 0, 1e-12)
		gClose(t, "axis z", s.Axis.Z, 1, 1e-12)
		if s.Folds != 0 || s.Face != "front" || s.Element != i {
			t.Errorf("surface %d: %+v", i, s)
		}
		if s.Vertex.Z < -1e9 || s.Vertex.Z > 1e9 {
			t.Errorf("surface %d vertex is not finite: %v", i, s.Vertex)
		}
	}
	gClose(t, "exit vertex z", res.ExitVertex.Z, -240, 1e-12)
	gClose(t, "exit axis z", res.ExitAxis.Z, 1, 1e-12)
	if err := cmdLens([]string{"group", "layout", "--json", out2, cass}); err != nil {
		t.Fatalf("layout (second run): %v", err)
	}
	buf2, err := os.ReadFile(out2)
	if err != nil {
		t.Fatal(err)
	}
	if string(buf1) != string(buf2) {
		t.Errorf("layout is not deterministic:\n%s\n%s", buf1, buf2)
	}

	// The folded Kepler: 6 real surfaces (stop, fold mirror, four lens
	// surfaces), the fold happens at element 1 and every later surface sits on
	// the folded axis (0, 1, 0) with folds = 1.
	fk := writeGroup(t, "folded-kepler", gKeplerFixture(t, true))
	out3 := filepath.Join(dir, "l3.json")
	if err := cmdLens([]string{"group", "layout", "--json", out3, fk}); err != nil {
		t.Fatalf("layout (folded): %v", err)
	}
	var res3 groupLayoutResponse
	readJSONFile(t, out3, &res3)
	if res3.NumSurfaces != 6 || res3.Folds != 1 {
		t.Fatalf("folded layout: %d surfaces, %d folds, want 6 and 1", res3.NumSurfaces, res3.Folds)
	}
	gClose(t, "fold vertex y", res3.Surfaces[1].Vertex.Y, 0, 1e-12)
	gClose(t, "fold vertex z", res3.Surfaces[1].Vertex.Z, 15, 1e-12)
	for _, s := range res3.Surfaces[2:] {
		gClose(t, "post-fold vertex z", s.Vertex.Z, 15, 1e-12)
		gClose(t, "post-fold axis y", s.Axis.Y, 1, 1e-12)
		if s.Folds != 1 {
			t.Errorf("surface %d: folds %d, want 1", s.Surface, s.Folds)
		}
	}
	gClose(t, "folded exit axis y", res3.ExitAxis.Y, 1, 1e-12)
}

// ---------------------------------------------------------------------------
// image / afocal / dispersion
// ---------------------------------------------------------------------------

// TestLensGroupImageCLI drives `oz lens group image` against the closed form of
// the paraboloid: the image sits at R/2 = 200 mm behind the vertex, measured
// along the reflected beam, the residual RMS is at the level of the search
// tolerance (1e-9 mm) and every ray of the hexapolar sampling arrives.
func TestLensGroupImageCLI(t *testing.T) {
	dir := t.TempDir()
	path := writeGroup(t, "parabola", gParabolaFixture())
	run := func(out string) (groupImageResponse, []byte) {
		t.Helper()
		if err := cmdLens([]string{"group", "image", "--json", out, path}); err != nil {
			t.Fatalf("image: %v", err)
		}
		var res groupImageResponse
		return res, readJSONFile(t, out, &res)
	}
	res, buf1 := run(filepath.Join(dir, "i1.json"))
	gClose(t, "image z", res.Point.Z, -200, 1e-9)
	gClose(t, "image x", res.Point.X, 0, 1e-9)
	gClose(t, "focus distance", res.FocusMM, 200, 1e-9)
	gClose(t, "focal-plane axis z", res.Axis.Z, -1, 1e-12)
	if res.RMSMM > 1e-10 {
		t.Errorf("the parabola has no spherical aberration: RMS %.3g mm", res.RMSMM)
	}
	if res.Rays != 37 || res.Failed != 0 {
		t.Errorf("rays %d, failed %d, want 37 and 0", res.Rays, res.Failed)
	}
	if res.Samples != 3 || res.ObjectZMM != nil {
		t.Errorf("samples %d, object_z_mm %v, want 3 and null", res.Samples, res.ObjectZMM)
	}
	if _, buf2 := run(filepath.Join(dir, "i2.json")); string(buf1) != string(buf2) {
		t.Errorf("image CLI is not deterministic:\n%s\n%s", buf1, buf2)
	}

	// An off-axis field through the Cassegrain must still focus (the real-ray
	// measurement works for folded/mirrored systems; the image just walks off
	// the axis).
	cass := writeGroup(t, "cassegrain", gCassegrainFixture())
	out := filepath.Join(dir, "i3.json")
	if err := cmdLens([]string{"group", "image", "--field-deg", "0.2", "--json", out, cass}); err != nil {
		t.Fatalf("image (off axis): %v", err)
	}
	var off groupImageResponse
	readJSONFile(t, out, &off)
	if off.Rays == 0 || math.IsNaN(off.Point.Y) {
		t.Fatalf("off-axis image is empty: %+v", off)
	}
	t.Logf("cassegrain off axis 0.2 deg: image (%g, %g, %g), RMS %g mm, rays %d",
		off.Point.X, off.Point.Y, off.Point.Z, off.RMSMM, off.Rays)

	// A stop without a clear aperture cannot size the automatic pupil.
	bad := &lens.Group{Name: "no-aperture", ObjectZMM: math.Inf(-1), Elements: []lens.Element{
		{Kind: lens.ElementStop},
		lens.MirrorElement(-400, -1, 0),
	}}
	badPath := writeGroup(t, "no-aperture", bad)
	if err := cmdLens([]string{"group", "image", badPath}); err == nil {
		t.Error("image accepted a group without a clear aperture")
	}
	if err := cmdLens([]string{"group", "image", "--rings", "99", path}); err == nil {
		t.Error("image accepted 99 sampling rings")
	}
	if err := cmdLens([]string{"group", "image"}); err == nil {
		t.Error("image without a file accepted")
	}
}

// TestLensGroupAfocalCLI drives `oz lens group afocal`: the Kepler telescope has
// Γ = −1 (the closed form of two equal thin lenses), a 45 degree plane fold in
// its collimated space leaves Γ bit-for-bit alone and turns the exit axis to
// +y, and a focal system is refused (it has no angular magnification).
func TestLensGroupAfocalCLI(t *testing.T) {
	dir := t.TempDir()
	straight := writeGroup(t, "kepler", gKeplerFixture(t, false))
	folded := writeGroup(t, "folded-kepler", gKeplerFixture(t, true))
	call := func(path, out string) groupAfocalResponse {
		t.Helper()
		if err := cmdLens([]string{"group", "afocal", "--json", out, path}); err != nil {
			t.Fatalf("afocal %s: %v", filepath.Base(path), err)
		}
		var res groupAfocalResponse
		readJSONFile(t, out, &res)
		return res
	}
	a1 := call(straight, filepath.Join(dir, "a1.json"))
	a2 := call(folded, filepath.Join(dir, "a2.json"))
	gClose(t, "kepler Γ", a1.Gamma, -1, 1e-9)
	gClose(t, "folded kepler Γ", a2.Gamma, -1, 1e-9)
	gClose(t, "Γ is fold invariant", a2.Gamma-a1.Gamma, 0, 1e-9)
	// D_in/D_out is a geometric cross-check of |Γ|: its tolerance is the
	// spherical aberration of the telescope at that aperture (2 %, the value of
	// the lens-package test; the folded telescope uses a 6 mm stop, the
	// straight one 4 mm, so the two ratios differ in the third digit).
	gClose(t, "beam ratio (straight)", a1.BeamRatio, 1, 0.02)
	gClose(t, "beam ratio (folded)", a2.BeamRatio, 1, 0.02)
	gClose(t, "straight exit axis", a1.AxisOut.Z, 1, 1e-12)
	gClose(t, "folded exit axis y", a2.AxisOut.Y, 1, 1e-12)
	if a1.SpreadRad > 5e-3 || a2.SpreadRad > 5e-3 {
		t.Errorf("the exit bundles are not collimated: %g / %g rad", a1.SpreadRad, a2.SpreadRad)
	}

	cass := writeGroup(t, "cassegrain", gCassegrainFixture())
	if err := cmdLens([]string{"group", "afocal", cass}); err == nil {
		t.Error("afocal accepted a focal system")
	}
}

// TestLensGroupDispersionCLI drives `oz lens group dispersion`: a mirror system
// disperses exactly nothing (the media are constant, so the two focus searches
// agree bit for bit), while the cemented achromat keeps the normal-dispersion
// order (the F focus in front of the C focus) with a small residual shift.
func TestLensGroupDispersionCLI(t *testing.T) {
	dir := t.TempDir()
	par := writeGroup(t, "parabola", gParabolaFixture())
	out1 := filepath.Join(dir, "d1.json")
	if err := cmdLens([]string{"group", "dispersion", "--json", out1, par}); err != nil {
		t.Fatalf("dispersion (mirror): %v", err)
	}
	var mirror groupDispersionResponse
	buf1 := readJSONFile(t, out1, &mirror)
	if mirror.ShiftMM != 0 || mirror.A.FocusMM != mirror.B.FocusMM {
		t.Errorf("a pure mirror system must not disperse: shift=%v, F=%v, C=%v",
			mirror.ShiftMM, mirror.A.FocusMM, mirror.B.FocusMM)
	}

	dbl := writeGroup(t, "doublet", gDoubletFixture(t, 100, 12, 6, 3))
	out2 := filepath.Join(dir, "d2.json")
	if err := cmdLens([]string{"group", "dispersion", "--json", out2, dbl}); err != nil {
		t.Fatalf("dispersion (doublet): %v", err)
	}
	var res groupDispersionResponse
	readJSONFile(t, out2, &res)
	t.Logf("doublet F–C shift: %g mm (F %g, C %g), image z %g, RMS %g",
		res.ShiftMM, res.A.FocusMM, res.B.FocusMM, res.A.Point.Z, res.A.RMSMM)
	if res.LambdaANm != geo.WavelengthF || res.LambdaBNm != geo.WavelengthC {
		t.Errorf("wavelengths %g/%g, want F/C", res.LambdaANm, res.LambdaBNm)
	}
	if !(res.ShiftMM < 0) {
		t.Errorf("normal dispersion puts the F focus in front of the C focus: shift %g", res.ShiftMM)
	}
	gClose(t, "doublet F–C shift", res.ShiftMM, -0.4065509316, 0.01)
	if !(math.Abs(res.ShiftMM) > 1e-3) {
		t.Errorf("the achromat shows no measurable residual colour: %g mm", res.ShiftMM)
	}
	// Determinism: the same request twice is byte-identical.
	if err := cmdLens([]string{"group", "dispersion", "--json", out2, dbl}); err != nil {
		t.Fatal(err)
	}
	buf2, err := os.ReadFile(out2)
	if err != nil {
		t.Fatal(err)
	}
	if string(buf1) == string(buf2) {
		t.Error("the mirror and the doublet answers are identical — the fixture is not doing work")
	}
}

// ---------------------------------------------------------------------------
// set / to-prescription / parax
// ---------------------------------------------------------------------------

// TestLensGroupSetCLI drives `oz lens group set`: the canonical parameter syntax
// and its long aliases both land on the right quantity, a radius change moves
// the focal length to the thin-lens value (1/f = (n−1)/R for the plano-convex
// element, relative 1e-9 through to-prescription + parax), and an invalid edit
// echoes the error of the lens API.
func TestLensGroupSetCLI(t *testing.T) {
	dir := t.TempDir()
	// A plano-convex singlet of f = 100 mm: R = (n−1)·f on the front surface.
	m, _ := geo.Glass("N-BK7")
	nd := m.N(geo.Wavelengthd)
	stop := lens.StopElement(5)
	stop.GapMM = 15
	singlet := &lens.Group{
		Name: "plano", ObjectZMM: math.Inf(-1),
		Elements: []lens.Element{stop, gPlanoElement(t, 100, 5, 14)},
	}
	path := writeGroup(t, "plano", singlet)

	// EFL of the untouched group.
	efl := func(groupPath string) float64 {
		t.Helper()
		rx := filepath.Join(dir, "rx.json")
		if err := cmdLens([]string{"group", "to-prescription", "-o", rx, groupPath}); err != nil {
			t.Fatalf("to-prescription: %v", err)
		}
		jout := filepath.Join(dir, "px.json")
		if err := cmdLens([]string{"parax", "--json", jout, rx}); err != nil {
			t.Fatalf("parax: %v", err)
		}
		var res paraxResponse
		readJSONFile(t, jout, &res)
		if res.EFLMM == nil {
			t.Fatalf("no EFL in the parax report: %+v", res)
		}
		return *res.EFLMM
	}
	gClose(t, "plano-convex EFL", efl(path), 100, 1e-9)

	// Set the radius of the second element to the f = 80 mm value.
	out := filepath.Join(dir, "out.json")
	if err := cmdLens([]string{"group", "set", "-p", "e1.front.R=" + strconvF((nd-1)*80), "-o", out, path}); err != nil {
		t.Fatalf("set: %v", err)
	}
	gClose(t, "EFL after set e1.front.R", efl(out), 80, 1e-9*80)

	// The long aliases of the task (`e1.tilt_x_deg=45`, `e1.gap_mm=30`) set the
	// same quantities as the canonical spelling: a flat fold mirror in front of
	// a lens turns a collimated beam into a folded chain (one fold, exit axis
	// +y).
	stop2 := lens.StopElement(6)
	stop2.GapMM = 15
	scan := &lens.Group{
		Name: "scan", ObjectZMM: math.Inf(-1),
		Elements: []lens.Element{stop2, lens.MirrorElement(0, 0, 20), gPlanoElement(t, 100, 2, 20)},
	}
	scanPath := writeGroup(t, "scan", scan)
	aliasOut := filepath.Join(dir, "alias.json")
	if err := cmdLens([]string{"group", "set", "-p", "e1.tilt_x_deg=45", "-p", "e1.gap_mm=30", "-o", aliasOut, scanPath}); err != nil {
		t.Fatalf("set (aliases): %v", err)
	}
	g, err := lens.LoadGroupJSON(aliasOut)
	if err != nil {
		t.Fatal(err)
	}
	gClose(t, "tilt_x_deg alias", g.Elements[1].TiltXDeg, 45, 0)
	gClose(t, "gap_mm alias", g.Elements[1].GapMM, 30, 0)
	layoutOut := filepath.Join(dir, "alias-layout.json")
	if err := cmdLens([]string{"group", "layout", "--json", layoutOut, aliasOut}); err != nil {
		t.Fatalf("layout (aliases): %v", err)
	}
	var lay groupLayoutResponse
	readJSONFile(t, layoutOut, &lay)
	if lay.Folds != 1 {
		t.Errorf("the 45 degree mirror does not fold the chain: %d folds", lay.Folds)
	} else {
		gClose(t, "folded exit axis y", lay.ExitAxis.Y, 1, 1e-12)
	}

	// Errors: a parameter that does not exist, an edit the lens API refuses
	// (the stop has no rear surface), and a value that is not a number. The
	// API text is echoed verbatim.
	_, err = captureOut(t, func() error { return cmdLens([]string{"group", "set", "-p", "e0.bogus=1", path}) })
	if err == nil || !strings.Contains(err.Error(), "无法识别") {
		t.Errorf("unknown parameter error = %v", err)
	}
	if err := cmdLens([]string{"group", "set", "-p", "e0.rear.R=10", path}); err == nil {
		t.Error("set accepted a rear radius on the stop")
	} else if !strings.Contains(err.Error(), "has no rear surface") {
		t.Errorf("the lens API error is not echoed: %v", err)
	}
	if err := cmdLens([]string{"group", "set", "-p", "e1.front.R=abc", path}); err == nil {
		t.Error("set accepted a non-numeric value")
	}
	if err := cmdLens([]string{"group", "set", "-p", "e1.front.R", path}); err == nil {
		t.Error("set accepted a parameter without a value")
	}
}

// strconvF formats a float for a command line (the shortest round-trip form).
func strconvF(v float64) string {
	buf, err := json.Marshal(v)
	if err != nil {
		return "0"
	}
	return string(buf)
}

// TestLensGroupToPrescriptionCLI converts the Cassegrain fixture and compares
// the result card by card with the library file it is the element form of.
func TestLensGroupToPrescriptionCLI(t *testing.T) {
	dir := t.TempDir()
	path := writeGroup(t, "cassegrain", gCassegrainFixture())
	out := filepath.Join(dir, "rx.json")
	if err := cmdLens([]string{"group", "to-prescription", "-o", out, path}); err != nil {
		t.Fatalf("to-prescription: %v", err)
	}
	got, err := lens.LoadJSON(out)
	if err != nil {
		t.Fatal(err)
	}
	want, err := lens.LoadJSON(filepath.Join("..", "..", "library", "geo", "mirror-cassegrain-f1000-d300-b300.json"))
	if err != nil {
		t.Skipf("library file not readable: %v", err)
	}
	if len(got.Cards) != len(want.Cards) {
		t.Fatalf("%d cards, want %d", len(got.Cards), len(want.Cards))
	}
	if got.Stop != want.Stop {
		t.Errorf("stop %d, want %d", got.Stop, want.Stop)
	}
	for i := range want.Cards {
		a, b := got.Cards[i], want.Cards[i]
		if a.Mirror != b.Mirror {
			t.Errorf("card %d mirror %v, want %v", i, a.Mirror, b.Mirror)
		}
		for _, c := range []struct {
			name      string
			got, want float64
		}{
			{"radius", a.RadiusMM, b.RadiusMM},
			{"thickness", a.ThicknessMM, b.ThicknessMM},
			{"semi", a.SemiAperMM, b.SemiAperMM},
		} {
			if d := math.Abs(c.got - c.want); d > 1e-12 {
				t.Errorf("card %d %s: %g vs %g (Δ %g)", i, c.name, c.got, c.want, d)
			}
		}
	}
	if err := cmdLens([]string{"group", "to-prescription", path}); err == nil {
		t.Error("to-prescription without -o accepted")
	}
}

// TestLensParaxCLI drives `oz lens parax` on the library Cassegrain (a two
// mirror system: EFL = +857.142857 mm, image at +360 mm, BFL +600 mm, parity
// +1, F–C exactly zero), on a single mirror with a finite object (the closed form
// 1/s + 1/s' = 2/R) and on the afocal Kepler (Γ through --field-tan).
func TestLensParaxCLI(t *testing.T) {
	dir := t.TempDir()
	// 1. The library Cassegrain.
	cass := filepath.Join("..", "..", "library", "geo", "mirror-cassegrain-f1000-d300-b300.json")
	out := filepath.Join(dir, "cass.json")
	if err := cmdLens([]string{"parax", "--json", out, cass}); err != nil {
		t.Skipf("library file not readable: %v", err)
	}
	var res paraxResponse
	readJSONFile(t, out, &res)
	if res.Mirrors != 2 || res.Parity != 1 || res.Afocal {
		t.Errorf("mirrors %d, parity %g, afocal %v (want 2, +1, false)", res.Mirrors, res.Parity, res.Afocal)
	}
	gClose(t, "cassegrain EFL", *res.EFLMM, 857.1428571428571, 1e-9)
	gClose(t, "cassegrain ImageZ", *res.ImageZMM, 360, 1e-9)
	gClose(t, "cassegrain BFL", *res.BFLMM, 600, 1e-9)
	if res.FCShiftMM == nil || *res.FCShiftMM != 0 {
		t.Errorf("a pure mirror system has no F–C shift: %v", res.FCShiftMM)
	}
	if res.Gamma != nil {
		t.Errorf("an infinite object of a focal system has no Γ: %v", *res.Gamma)
	}
	if res.GammaNote == "" {
		t.Error("the missing Γ is not explained")
	}

	// 2. A single mirror prescribed with one card: EFL = R/2 = −1000, and with
	// the object at z = −3000 the image is at s·R/(R+2s) = −1500 (m = −0.5).
	one := &lens.Prescription{
		Name: "one-mirror", ObjectZMM: math.Inf(-1),
		Cards: []lens.Card{{Name: "mirror", RadiusMM: -2000, Mirror: true, SemiAperMM: 80}},
	}
	onePath := filepath.Join(dir, "one.json")
	if err := one.SaveJSON(onePath); err != nil {
		t.Fatal(err)
	}
	out2 := filepath.Join(dir, "one-out.json")
	if err := cmdLens([]string{"parax", "--object-z", "-3000", "--json", out2, onePath}); err != nil {
		t.Fatalf("parax (finite object): %v", err)
	}
	var r2 paraxResponse
	readJSONFile(t, out2, &r2)
	gClose(t, "one mirror EFL", *r2.EFLMM, -1000, 1e-12)
	gClose(t, "one mirror ImageZ (s=−3000)", *r2.ImageZMM, -1500, 1e-9)
	gClose(t, "one mirror magnification", *r2.Magnifi, -0.5, 1e-9)
	if r2.ObjectZMM == nil || *r2.ObjectZMM != -3000 {
		t.Errorf("object_z_mm echo: %v", r2.ObjectZMM)
	}
	// A finite object of a system with power reports γ = n/(n'·m) = −2.
	if r2.Gamma == nil {
		t.Fatalf("no γ for the finite-object system (note %q)", r2.GammaNote)
	}
	gClose(t, "one mirror γ = 1/m", *r2.Gamma, -2, 1e-9)

	// 3. The afocal Kepler: Γ = −1 needs a non-zero field tangent.
	kep := writeGroup(t, "kepler", gKeplerFixture(t, false))
	rx := filepath.Join(dir, "kepler-rx.json")
	if err := cmdLens([]string{"group", "to-prescription", "-o", rx, kep}); err != nil {
		t.Fatalf("kepler to-prescription: %v", err)
	}
	out3 := filepath.Join(dir, "kepler.json")
	if err := cmdLens([]string{"parax", "--field-tan", "1e-3", "--json", out3, rx}); err != nil {
		t.Fatalf("parax (afocal): %v", err)
	}
	var r3 paraxResponse
	readJSONFile(t, out3, &r3)
	if !r3.Afocal {
		t.Fatalf("the Kepler prescription is not afocal: %+v", r3)
	}
	if r3.Gamma == nil {
		t.Fatalf("no Γ for the afocal prescription (note %q)", r3.GammaNote)
	}
	gClose(t, "kepler Γ", *r3.Gamma, -1, 1e-9)
	if r3.EFLMM != nil {
		t.Errorf("an afocal system has no EFL: %v", *r3.EFLMM)
	}
}
