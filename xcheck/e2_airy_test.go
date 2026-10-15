package xcheck

import (
	"math"
	"testing"
)

// besselJ1Zero1 is the first positive zero j_{1,1} of the Bessel function J1
// (the first dark ring of an Airy pattern is at j_{1,1}/π · λF# = 1.21967·λF#).
const besselJ1Zero1 = 3.8317059702075123

// E2 (plan §6.4): the first dark ring of an ideal circular aperture is at
// r1 = 1.2197·λ·F#.
//
// Setup: aperture diameter D = 20 mm (R = 128 pixels across the radius),
// focal length f = 100 mm, so F# = f/D = 5; the optics Fraunhofer focal field on
// an N = 2048 grid has pixel p = λ·f/(N·dx) = λF#/8 with dx = D/(2R) (the
// test reports p in µm). The ring is taken from the mean of the four
// axis-aligned cuts through the centre, as the first local minimum with a
// three-point parabolic sub-pixel refinement.
//
// Criterion: |r_measured − j_{1,1}/π·λF#| ≤ 1 pixel. The analytic constant is
// the Airy first zero, r1 = j_{1,1}/π · λF# with j_{1,1} = 3.8317059702.
func TestE2AiryFirstDarkRing(t *testing.T) {
	const (
		focalMM    = 100.0
		apertureMM = 20.0
		R          = 128
		N          = 2048
	)
	lamM := dLineM
	dxM := apertureMM * 1e-3 / float64(2*R)
	lat := newLattice(N, R)
	f := lat.fraunhofer(lamM, dxM, focalMM*1e-3, nil)

	Fnum := focalMM / apertureMM
	lamF := lamM * Fnum // λF# (m)
	p := f.DX           // image pixel (m)
	c := N / 2
	ic := f.Intensity(c*N + c)

	// Radial profile: mean of the four axis cuts, index r−1 = radius r px.
	prof := make([]float64, N/2-1)
	for r := 1; r < N/2; r++ {
		s := f.Intensity(c*N+c+r) + f.Intensity(c*N+c-r) +
			f.Intensity((c+r)*N+c) + f.Intensity((c-r)*N+c)
		prof[r-1] = s / 4 / ic
	}
	// First local minimum (after the central lobe), parabolic refinement.
	minIdx := -1
	for i := 2; i < len(prof)-1; i++ {
		if prof[i] < prof[i-1] && prof[i] <= prof[i+1] {
			minIdx = i
			break
		}
	}
	if minIdx < 0 {
		t.Fatal("E2: no dark ring found in the radial profile")
	}
	y0, y1, y2 := prof[minIdx-1], prof[minIdx], prof[minIdx+1]
	sub := float64(minIdx+1) + 0.5*(y0-y2)/(y0-2*y1+y2) // radius in pixels (index i ↔ radius i+1)

	theoPx := besselJ1Zero1 / math.Pi * lamF / p
	dev := sub - theoPx
	t.Logf("E2 aperture D = %.1f mm, f = %.1f mm, F# = %.3f, λF# = %.4f µm", apertureMM, focalMM, Fnum, lamF*1e6)
	t.Logf("E2 pixel p = %.4f µm = λF#/%.3f, N = %d, pupil radius %d px", p*1e6, lamF/p, N, R)
	t.Logf("E2 first dark ring: theory %.4f px (= %.4f µm = 1.21967·λF#); sampled minimum at r = %d px; parabolic %.4f px; deviation %+.3f px",
		theoPx, theoPx*p*1e6, minIdx+1, sub, dev)
	if math.Abs(dev) > 1 {
		t.Errorf("E2: first dark ring at %.3f px, theory %.3f px: |Δ| = %.3f px > 1 px", sub, theoPx, math.Abs(dev))
	}
}
