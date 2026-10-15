package xcheck

import (
	"math"
	"testing"
)

// E3 (plan §6.4): Strehl ratio of a pure defocus W020 = 0.25 λ, against an
// exact Fourier integral.
//
// Wavefront: W(ρ) = w·ρ² waves over the unit pupil, w = 0.25 (W020 in the
// usual Zernike/Seidel convention, the edge OPD being w·λ).
//
// Exact reference. For the unit disc,
//
//	S = |∫_0^1 2ρ·exp(i·2π·w·ρ²) dρ|² = (sin(πw)/(πw))²,
//
// evaluated in closed form, and also as a numerical Fourier integral with a
// 200 × 256 product rule (checked against the closed form to 1e-10). The
// Maréchal value exp(−(2πσ)²), σ = w/√12, is 0.8141450; it is NOT an exact
// reference (its deviation from the integral is reported).
//
// Wave side: the optics Fraunhofer focal field on lattice pupils of radius
// R = 32, 64, 128 pixels (N = 8R, pixel λF#/4 for F# = 5): S_R = peak ratio to
// the unaberrated peak. The lattice (pixel centres inside the circle) is a
// discretisation of the disc, so S_R converges to the exact value as R grows.
//
// Criteria: |S_128 − S_exact| / S_exact ≤ 1e-3, and the error decreases
// monotonically over R = 32, 64, 128. The 1e-3 bound is 3× the error measured at
// R = 128; it is a threshold set on this sampling, not the plan's 0.3 %, which
// was against the Maréchal value (see the report).
func TestE3DefocusStrehlAgainstIntegral(t *testing.T) {
	const w = 0.25
	closed := math.Pow(math.Sin(math.Pi*w)/(math.Pi*w), 2)

	qp, qw := diskQuadrature(200, 256)
	wv := make([]float64, len(qp))
	for i, p := range qp {
		wv[i] = w * (p.PX*p.PX + p.PY*p.PY)
	}
	sQuad, _ := integralStrehl(wv, qw)
	if e := relErr(sQuad, closed); e > 1e-10 {
		t.Errorf("E3: quadrature %.12f vs closed form %.12f: relative difference %.3e", sQuad, closed, e)
	}

	var errs []float64
	for _, R := range []int{32, 64, 128} {
		lat := newLattice(8*R, R)
		opdM := make([]float64, len(lat.px))
		for i := range lat.px {
			r2 := lat.px[i]*lat.px[i] + lat.py[i]*lat.py[i]
			opdM[i] = w * r2 * dLineM
		}
		dxM := 10e-3 / float64(R) // pupil radius 10 mm
		zM := 100e-3              // focal length 100 mm (F/5)
		ideal := lat.fraunhofer(dLineM, dxM, zM, nil)
		aber := lat.fraunhofer(dLineM, dxM, zM, opdM)
		s := peak(aber) / peak(ideal)
		errs = append(errs, relErr(s, closed))
		t.Logf("E3 R = %3d px (N = %4d, %6d samples): wave S = %.7f, |S − S_exact|/S_exact = %.3e",
			R, lat.N, len(lat.px), s, errs[len(errs)-1])
	}
	t.Logf("E3 exact S = %.10f (closed form (sin πw/πw)², w = %.2f); quadrature agrees to %.1e",
		closed, w, relErr(sQuad, closed))
	t.Logf("E3 Maréchal S = %.7f: |S_M − S_exact|/S_exact = %.3e (not the reference)",
		marechal(w/math.Sqrt(12)), relErr(marechal(w/math.Sqrt(12)), closed))
	if errs[2] > 1e-3 {
		t.Errorf("E3: wave Strehl error %.3e at R = 128 exceeds 1e-3", errs[2])
	}
	if !(errs[2] < errs[1] && errs[1] < errs[0]) {
		t.Errorf("E3: discretisation error does not decrease with R: %.3e, %.3e, %.3e", errs[0], errs[1], errs[2])
	}
}
