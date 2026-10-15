// The exit pupil of a system is the image of the aperture stop formed by the
// elements after it, evaluated in the image space. The transfer from the stop
// plane to a plane a gap L beyond the last surface left-multiplies by the
// appended transfer [[1, L/nk],[0,1]], so the conjugate condition B = 0 reads
// B + L·D/nk = 0 — the appended gap is governed by the D entry, exactly as the
// prepended gap of the entrance pupil is governed by the A entry.
//
// These tests pin the exit pupil of two systems against the real tracer:
//   - a singlet with the stop 50 mm in front, whose exit pupil is virtual at
//     z ≈ −51 (the old formula, which used the A entry, reported z ≈ +0.7,
//     i.e. the lens plane itself);
//   - a concave mirror of R = −2000 with the stop 60 mm in front, whose exit
//     pupil is virtual at z = 60 + 3000/47 with 50/47 of the stop radius.
package geo

import (
	"math"
	"testing"
)

// exitPupilCrossing traces two real rays from one point of the stop plane and
// returns the intersection of their emergent lines and the height there. It is
// the real-trace image of that stop point: for the stop centre it must be
// (ExitPupilZ, 0); for the stop edge (semi, ExitPupilZ, |ExitPupilMag|·semi).
func exitPupilCrossing(t *testing.T, sys *System, y0, u1, u2 float64) (z, y float64) {
	t.Helper()
	r1 := sys.Trace(NewRay(V(0, y0, 0), V(0, u1, 1), Wavelengthd))
	r2 := sys.Trace(NewRay(V(0, y0, 0), V(0, u2, 1), Wavelengthd))
	if !r1.OK() || !r2.OK() {
		t.Fatalf("ray trace failed: %v / %v", r1.Error(), r2.Error())
	}
	s1 := r1.Steps[len(r1.Steps)-1]
	s2 := r2.Steps[len(r2.Steps)-1]
	z = pairCrossZTrue(s1.Point, s1.DirOut, s2.Point, s2.DirOut)
	if math.IsNaN(z) {
		t.Fatalf("the emergent lines do not cross")
	}
	y = s1.Point.Y + (z-s1.Point.Z)*s1.DirOut.Y/s1.DirOut.Z
	return z, y
}

// TestExitPupilSingletConjugation checks [stop, gap 50, biconvex N-BK7]: the
// exit pupil is the image of the stop by the lens, virtual in front of the
// system (z ≈ −51 for f ≈ 100), with about twice the stop radius.
func TestExitPupilSingletConjugation(t *testing.T) {
	n := N_BK7.N(Wavelengthd)
	R := 2 * (n - 1) * 100 // equiconvex N-BK7, f ≈ 100 mm
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(5)}, // stop at z = 0
			{Profile: Sphere{R: R}, Aperture: CircleAperture(12)},
			{Profile: Sphere{R: -R}, Aperture: CircleAperture(12)},
		},
		[]float64{50, 6}, []Medium{Air{}, Air{}, N_BK7, Air{}})
	sys.StopIndex = 0
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	// Thin-lens baseline: a stop s = 50 in front of f ≈ 100 images to
	// z = −s·f/(s−f) = −50 (virtual, in front of the lens); the thick lens
	// shifts it by a few mm. The old A-entry formula returned z ≈ 0 (the
	// lens plane) and about half the true radius, so both bounds below also
	// pin the fix. The exact values are the frozen engine output.
	if par.ExitPupilZ > -45 || par.ExitPupilZ < -60 {
		t.Errorf("exit pupil z = %g mm, want in [−60, −45] (virtual, in front of the lens)", par.ExitPupilZ)
	}
	checkRel(t, "singlet exit pupil z (frozen)", par.ExitPupilZ, -53.172441692855926, 1e-12)
	checkRel(t, "singlet exit pupil radius (frozen)", par.ExitPupilRadius, 10.305748056428532, 1e-12)
	checkRel(t, "singlet exit pupil mag (frozen)", par.ExitPupilMag, 2.0611496112857064, 1e-12)
	if math.Abs(par.ExitPupilRadius-2*5) > 0.5 {
		t.Errorf("exit pupil radius = %g mm, want ≈ 2×5 (|m| = 2)", par.ExitPupilRadius)
	}
	if math.Abs(par.ExitPupilMag-2) > 0.1 {
		t.Errorf("exit pupil magnification = %g, want ≈ 2", par.ExitPupilMag)
	}
	// Real rays from the stop centre and from the stop edge: the emergent
	// lines must cross at the exit pupil plane, the edge two times higher.
	// The tolerances allow the third-order aberration of the traced rays.
	const (
		u1, u2 = 0.005, 0.015
	)
	zc, yc := exitPupilCrossing(t, sys, 0, u1, u2)
	if math.Abs(zc-par.ExitPupilZ) > 0.1 {
		t.Errorf("real-ray centre crossing z = %.9f, paraxial exit pupil z = %.9f", zc, par.ExitPupilZ)
	}
	if math.Abs(yc) > 0.01 {
		t.Errorf("the image of the stop centre is at height %g, want 0", yc)
	}
	ze, ye := exitPupilCrossing(t, sys, 5, u1, u2)
	if math.Abs(ze-par.ExitPupilZ) > 1.5 {
		t.Errorf("real-ray edge crossing z = %.9f, paraxial exit pupil z = %.9f (a plane object images to one plane)", ze, par.ExitPupilZ)
	}
	if rel := math.Abs(ye-par.ExitPupilRadius) / par.ExitPupilRadius; rel > 0.03 {
		t.Errorf("real-ray edge height %g vs exit pupil radius %g (rel %.3g)", ye, par.ExitPupilRadius, rel)
	}
	t.Logf("singlet: exit pupil z %.17g (real ray %.9f), radius %.17g (real ray %.17g), mag %.17g",
		par.ExitPupilZ, zc, par.ExitPupilRadius, ye, par.ExitPupilMag)
}

// TestExitPupilMirrorConjugation checks [plane stop, gap 60, concave mirror
// R = −2000]. The mirror images the stop to z = 60 + s·R/(R+2s) = 60+3000/47
// with magnification |R/(R+2s)| = 50/47 — both closed forms are exact in the
// paraxial engine and are reproduced by the real tracer to the round-off.
func TestExitPupilMirrorConjugation(t *testing.T) {
	sys := stopMirrorPair(t, 60, 10, -2000, 0, 50)
	par, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	const (
		wantZ = 60 + 3000.0/47.0
		wantR = 10 * (50.0 / 47.0)
		wantM = 50.0 / 47.0
	)
	checkRel(t, "mirror exit pupil z", par.ExitPupilZ, wantZ, 1e-12)
	checkRel(t, "mirror exit pupil radius", par.ExitPupilRadius, wantR, 1e-12)
	checkRel(t, "mirror exit pupil mag", par.ExitPupilMag, wantM, 1e-12)
	// The entrance pupil of a stop at surface 0 is the stop itself.
	checkRel(t, "mirror entrance pupil z", par.EntrancePupilZ, 0, 1e-15)
	checkRel(t, "mirror entrance pupil radius", par.EntrancePupilRadius, 10, 1e-15)
	// Real rays: the emergent lines of the reflected bundle cross at the
	// virtual pupil behind the mirror (the lines extend through it). The
	// tolerance allows the third-order coma of the 9.5° field of the edge
	// point; the centre crossing is exact to the trace precision.
	zc, yc := exitPupilCrossing(t, sys, 0, 0.001, 0.002)
	if math.Abs(zc-wantZ) > 1e-4 {
		t.Errorf("real-ray crossing z = %.9f, want %.9f", zc, wantZ)
	}
	if math.Abs(yc) > 1e-3 {
		t.Errorf("the image of the stop centre is at height %g, want 0", yc)
	}
	ze, ye := exitPupilCrossing(t, sys, 10, 0.001, 0.002)
	if math.Abs(ze-wantZ) > 0.15 {
		t.Errorf("real-ray edge crossing z = %.9f, want %.9f", ze, wantZ)
	}
	if math.Abs(ye-wantR) > 1e-3 {
		t.Errorf("real-ray edge height %.9f, want %.9f", ye, wantR)
	}
	t.Logf("mirror: exit pupil z %.9f, radius %.9f, mag %.9f (real ray z %.9f, y %.9f)",
		par.ExitPupilZ, par.ExitPupilRadius, par.ExitPupilMag, ze, ye)
}
