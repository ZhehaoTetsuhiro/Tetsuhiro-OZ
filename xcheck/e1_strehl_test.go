package xcheck

import (
	"testing"
)

// E1 (plan §6.4): Strehl ratio of the real achromat at its geometric best
// focus, computed from the traced wavefront on two sides.
//
//	geometric side: the diffraction integral |∫ exp(i·2πW) dA|² / (∫ dA)² of the
//	  traced OPD (exact product rule on the pupil: 64 Gauss–Legendre radial nodes
//	  × 128 equispaced azimuths), and the Maréchal approximation exp(−(2πσ)²) of
//	  the same OPD's RMS σ;
//	wave side: the peak of the optics Fraunhofer focal field of the same OPD on
//	  a 256-pixel pupil (N = 1024), divided by the peak of the unaberrated field
//	  on the same grid (the peak ratio is the Strehl ratio).
//
// Criterion: |S_wave − S_exact| / S_exact ≤ 1 %. The Maréchal value is reported
// with its deviation from the exact integral and asserted at the same 1 %.
//
// Why Fraunhofer and not the lens element + ASM: the lens phase
// exp(−ik r²/2f) has local spatial frequency r/(λf); at the pupil edge
// (r = 10 mm, f = 100 mm) that is 170 cycles/mm, against a Nyquist limit of
// 6.4 cycles/mm for this 256-pixel pupil (pixel 0.078 mm). Propagating over
// the full pupil with the lens element would therefore alias; the Fraunhofer
// kernel is the exact lens-focused field, the quadratic phase having been
// absorbed by the lens.
func TestE1StrehlGeoVsWave(t *testing.T) {
	sys, par, pup := achromatDesign(t)
	zB, rmsB := bestFocusZ(t, sys, pup, par.ImageZ)

	// Geometric side (product rule on the pupil).
	qp, qw := diskQuadrature(64, 128)
	gq := rayOPD(t, sys, pup, qp, zB)
	if nan := countNaN(gq.OPD); nan > 0 {
		t.Fatalf("E1: %d of %d quadrature rays are blocked", nan, len(gq.OPD))
	}
	sExact, rmsQ := integralStrehl(wavesOf(gq.OPD), qw)
	sMar := marechal(rmsQ)

	// Wave side (lattice pupil through the optics kernel).
	lat := newLattice(1024, 128)
	gl := rayOPD(t, sys, pup, lat.points(), zB)
	if nan := countNaN(gl.OPD); nan > 0 {
		t.Fatalf("E1: %d of %d lattice rays are blocked", nan, len(gl.OPD))
	}
	zM := gl.SphereRadius * 1e-3              // focal distance: pupil plane to image point C (m)
	dxM := pup.Radius / float64(lat.R) * 1e-3 // pupil pixel (m)
	ideal := lat.fraunhofer(dLineM, dxM, zM, nil)
	aber := lat.fraunhofer(dLineM, dxM, zM, opdToM(gl.OPD))
	sWave := peak(aber) / peak(ideal)

	dWave := relErr(sWave, sExact)
	dMar := relErr(sMar, sExact)
	t.Logf("E1 best focus z_g = %.7f mm (paraxial ImageZ %.6f), RMS = %.6f λ", zB, par.ImageZ, rmsB)
	t.Logf("E1 geometric exact integral S = %.7f (product rule 64×128, RMS %.6f λ)", sExact, rmsQ)
	t.Logf("E1 geometric Maréchal    S = %.7f  |S_M − S_exact|/S_exact = %.3e", sMar, dMar)
	t.Logf("E1 wave peak ratio (optics Fraunhofer, N=%d, %d px pupil, pixel %.4f µm = λF#/%.2f, F#=%.4f) S = %.7f  |S_W − S_exact|/S_exact = %.3e",
		lat.N, lat.R, aber.DX*1e6, dLineM*(zM/(2*dxM*float64(lat.R)))/aber.DX, zM/(2*pup.Radius*1e-3), sWave, dWave)
	if dWave > 0.01 {
		t.Errorf("E1: wave Strehl %.6f vs exact geometric %.6f: relative difference %.3e > 1%%", sWave, sExact, dWave)
	}
	if dMar > 0.01 {
		t.Errorf("E1: Maréchal %.6f vs exact geometric %.6f: relative difference %.3e > 1%%", sMar, sExact, dMar)
	}
}

// countNaN returns the number of blocked (NaN) OPD samples.
func countNaN(v []float64) int {
	n := 0
	for _, x := range v {
		if x != x {
			n++
		}
	}
	return n
}
