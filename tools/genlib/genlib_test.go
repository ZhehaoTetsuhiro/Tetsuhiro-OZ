package main

import (
	"fmt"
	"math"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// Regression tests of the three generator failure classes. They are small
// on purpose: each builds one or two candidates, not the library.

func resetRecords() {
	buildFails, candFails, merges, renames = nil, nil, nil, nil
}

// TestTrimNumKeepsIntegerDigits: trimNum used to strip the zeros of the
// integer part (100 -> "1"), which made mirror-parabola-f100-f4 and
// mirror-parabola-f1000-f4 collide on one file name.
func TestTrimNumKeepsIntegerDigits(t *testing.T) {
	cases := []struct {
		x    float64
		dec  int
		want string
	}{
		{100, 0, "100"}, {1000, 0, "1000"}, {200, 0, "200"}, {2.5, 2, "2.5"},
		{2, 1, "2"}, {0.10, 2, "0.1"}, {-25, 0, "-25"}, {1.40, 2, "1.4"},
	}
	for _, c := range cases {
		if got := trimNum(c.x, c.dec); got != c.want {
			t.Errorf("trimNum(%v, %d) = %q, want %q", c.x, c.dec, got, c.want)
		}
	}
}

// TestDuplicateFileNameDoesNotPanic reproduces the original crash
// "panic: duplicate file name mirror-parabola-f1-f4.json" (layout.go:40):
// two different parabolas must get two distinct file names; the same design
// twice is merged; a clash of file names with different optics gets a suffix.
func TestDuplicateFileNameDoesNotPanic(t *testing.T) {
	resetRecords()
	l := newLib()
	parabola(l, 100, 4)
	parabola(l, 1000, 4)
	if len(l.items) != 2 {
		t.Fatalf("two different parabolas: %d items, want 2", len(l.items))
	}
	if l.items[0].file == l.items[1].file {
		t.Fatalf("file names collide: %s", l.items[0].file)
	}
	if l.items[0].file != "mirror-parabola-f100-f4.json" || l.items[1].file != "mirror-parabola-f1000-f4.json" {
		t.Errorf("unexpected names %s, %s", l.items[0].file, l.items[1].file)
	}

	// the same optical design twice is merged, not added
	parabola(l, 100, 4)
	if len(l.items) != 2 || len(merges) != 1 {
		t.Errorf("identical design: items=%d merges=%d, want 2 and 1", len(l.items), len(merges))
	}

	// same file name, different optics: readable suffix, no panic
	p1, err := buildLayout("a", []grp{{name: "L1", curv: seedBic("N-BK7", 50), glasses: []string{"N-BK7", ""}, thick: 3, after: 30}}, 2.5, 50, 5, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	p2, err := buildLayout("b", []grp{{name: "L1", curv: seedBic("N-BK7", 80), glasses: []string{"N-BK7", ""}, thick: 3, after: 30}}, 2.5, 80, 6, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	l.add("same-name.json", p1, "singlet", 1, []string{"a", "b", "c"}, "")
	l.add("same-name.json", p2, "singlet", 1, []string{"a", "b", "c"}, "")
	if !l.seen["same-name-alt2.json"] {
		t.Errorf("clashing name not suffixed: seen=%v", l.seen)
	}
	if len(renames) != 1 {
		t.Errorf("renames=%d, want 1", len(renames))
	}
}

// TestSolveEFLNoSignChangeReturnsError: a plane-parallel plate has no finite
// focal length, so no radius scale can reach the target. The solver must
// return an error and leave the radii untouched (no residual-1e16 system).
func TestSolveEFLNoSignChangeReturnsError(t *testing.T) {
	p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
		{Name: "stop", RadiusMM: math.Inf(1), ThicknessMM: 5, SemiAperMM: 5},
		{Name: "plate", RadiusMM: 0, ThicknessMM: 3, Glass: "N-BK7", SemiAperMM: 5},
		{Name: "plate", RadiusMM: 0, ThicknessMM: 10, SemiAperMM: 5},
	}}
	before := append([]lens.Card(nil), p.Cards...)
	err := solveEFL(p, radiiIdx(p.Cards), 100)
	if err == nil {
		t.Fatal("solveEFL on a plane-parallel plate returned nil")
	}
	for i := range before {
		if p.Cards[i].RadiusMM != before[i].RadiusMM {
			t.Errorf("card %d radius changed by a failed solve", i)
		}
	}
}

// checkSolvedOrReported asserts the contract of finishLens: either a
// prescription whose EFL meets the target to 1e-6 relative, or an error.
func checkSolvedOrReported(t *testing.T, name string, p *lens.Prescription, err error, target float64) {
	t.Helper()
	if err != nil {
		if !strings.Contains(err.Error(), "solveEFL") && !strings.Contains(err.Error(), "EFL residual") &&
			!strings.Contains(err.Error(), "final solve") {
			t.Errorf("%s: unexpected error text: %v", name, err)
		}
		return
	}
	par, perr := p.Paraxial(geo.Wavelengthd)
	if perr != nil {
		t.Fatalf("%s: accepted prescription has no paraxial solution: %v", name, perr)
	}
	if rel := math.Abs(par.EFL-target) / math.Abs(target); !(rel <= 1e-6) {
		t.Errorf("%s: accepted with EFL residual %.3g (target %g, got %g)", name, rel, target, par.EFL)
	}
}

// TestBendFamilySolveSolvesOrReports: the bend family k=1.3 and k=2 failed
// with "no sign change on [0.001,1000]". Each must now either solve to the
// target or report an error — never a silent garbage system.
func TestBendFamilySolveSolvesOrReports(t *testing.T) {
	for _, k := range []float64{1.3, 2} {
		target, sign := 100.0, 1.0
		if k > 1 {
			target, sign = -100, -1
		}
		g := grp{name: "L1", curv: []float64{sign, sign * k}, glasses: []string{"N-BK7", ""}, thick: thickFor(100), after: 30}
		p, err := buildLayout("bend", []grp{g}, 2.5, target, 10, math.Inf(-1))
		checkSolvedOrReported(t, fmt.Sprintf("bend k=%g", k), p, err, target)
	}
}

// TestMacroSolveSolvesOrReports: the macro layout produced residuals of
// 2.3e16 (f=50) and 2.3e16 (f=100). It must now solve or report.
func TestMacroSolveSolvesOrReports(t *testing.T) {
	for _, f := range []float64{50, 100} {
		fnum := fnumFor(f, 4)
		p, err := buildLayout("macro", layoutMacro(f), 2.5, f, f/(2*fnum), math.Inf(-1))
		checkSolvedOrReported(t, fmt.Sprintf("macro f=%g", f), p, err, f)
	}
}

// microLayout is the doublet used by the 60x microscope objective.
func microLayout(f float64) []grp {
	return []grp{{name: "doublet", curv: seedAchromat("N-BK7", "F2", f), glasses: []string{"N-BK7", "F2", ""}, thick: 0.08 * f, after: 1.2 * f}}
}

// TestObjectDistanceClippedToSurface0: the 60x objective solve used to bracket
// the object distance across surface 0 ("object at z=0.359 is not in front of
// surface 0"). The solved object must lie in front of surface 0 and reproduce
// the magnification.
func TestObjectDistanceClippedToSurface0(t *testing.T) {
	const f, mag = 2.8, 60.0
	p, err := buildLayout("micro", microLayout(f), 0.2, f, 0.2, -f*(1+1/mag)*1.3)
	if err != nil {
		t.Fatal(err)
	}
	z, err := solveObjectDistance(p, -mag, f)
	if err != nil {
		t.Fatalf("solve: %v", err)
	}
	sys, _ := p.System()
	if !(z < sys.VertexZ(0)) {
		t.Fatalf("object z=%g is not in front of surface 0 (z=%g)", z, sys.VertexZ(0))
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if rel := math.Abs(par.Magnification+mag) / mag; rel > 1e-6 {
		t.Errorf("magnification %g, want %g (residual %.3g)", par.Magnification, -mag, rel)
	}
}

// TestObjectDistanceOutOfRangeReturnsError: a demagnification of 1e-6 needs an
// object about a million focal lengths away, outside the scanned region. The
// solve must report that instead of returning a distant or clipped object.
func TestObjectDistanceOutOfRangeReturnsError(t *testing.T) {
	const f = 2.8
	p, err := buildLayout("micro", microLayout(f), 2.5, f, 0.2, -f*2*1.3)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := solveObjectDistance(p, -1e-6, f); err == nil {
		t.Fatal("out-of-range magnification was solved")
	}
}

// TestCandidateFailureIsRecordedNotWritten: a failing construction reaches
// the failure list exactly once and adds no item.
func TestCandidateFailureIsRecordedNotWritten(t *testing.T) {
	resetRecords()
	l := newLib()
	before := len(buildFails)
	failf("singlet bend k=%g: %v", 1.3, fmt.Errorf("solveEFL: no sign change"))
	if len(buildFails) != before+1 || len(l.items) != 0 {
		t.Fatalf("construction failure not recorded once: fails=%d items=%d", len(buildFails), len(l.items))
	}
	if !strings.Contains(buildFails[len(buildFails)-1], "k=1.3") {
		t.Errorf("failure text lacks its parameters: %q", buildFails[len(buildFails)-1])
	}
}
