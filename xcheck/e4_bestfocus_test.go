package xcheck

import (
	"math"
	"testing"

	"tetsuhiro-oz/optics"
)

// zScanPeak finds the defocus dz (µm) that maximises the peak intensity of the
// focal-plane field base after an angular-spectrum propagation (optics
// MethodASM) by dz. The scan is coarse (2 µm steps over ±40 µm) and then refined
// by golden-section search on the peak to 1 nm. It returns dz and the peak
// intensity there.
func zScanPeak(base *optics.Field, lamM float64) (float64, float64) {
	ctx := func() *optics.Context { return &optics.Context{Wavelength: lamM, Evanescent: "decay"} }
	peakAt := func(dzUM float64) float64 {
		c := base.Clone()
		if err := optics.Propagate(c, dzUM*1e-6, optics.MethodASM, ctx()); err != nil {
			panic(err)
		}
		return peak(c)
	}
	bestI, bestP := 0, math.Inf(-1)
	for i := -20; i <= 20; i++ {
		if v := peakAt(float64(i) * 2); v > bestP {
			bestI, bestP = i, v
		}
	}
	dz := goldenMin(func(x float64) float64 { return -peakAt(x) }, float64(bestI)*2-2, float64(bestI)*2+2, 1e-3)
	return dz, peakAt(dz)
}

// E4 (plan §6.4): the geometric best focus and the wave best focus of the
// real achromat agree within one pixel, and the wave Strehl is ≥ 0.98.
//
//	geometric best focus z_g: the image plane minimising the RMS OPD (golden
//	  section, 32 × 64 pupil product rule);
//	wave best focus z_w: the defocus maximising the peak intensity of the wave
//	  field. The field is the optics Fraunhofer focal field (N = 1024, 256-px
//	  pupil) of the OPD traced to z_g; the defocus scan is an angular-spectrum
//	  propagation (optics MethodASM) of that focal field.
//
// Criteria: |z_w − z_g| ≤ 1 pixel, the pixel being the transverse pixel p of
// the focal grid (λF#/4, about 0.72 µm), with z_w and z_g both absolute image
// positions; the wave Strehl (peak ratio at z_w) must be ≥ 0.98. Controls: the
// same scan on the unaberrated pupil must reproduce the analytic diffraction
// focal shift −48·f³/(k²a⁴) of an F/5 pupil (5 %), and the reference-sphere
// radius must equal z_g − zQ (zQ: vertex of the last surface).
func TestE4BestFocusGeoVsWave(t *testing.T) {
	sys, par, pup := achromatDesign(t)
	zG, rmsG := bestFocusZ(t, sys, pup, par.ImageZ)

	lat := newLattice(1024, 128)
	gl := rayOPD(t, sys, pup, lat.points(), zG)
	if nan := countNaN(gl.OPD); nan > 0 {
		t.Fatalf("E4: %d of %d lattice rays are blocked", nan, len(gl.OPD))
	}
	// Absolute geometry. The chief ray leaves the last surface on axis at its
	// vertex zQ, and the reference sphere of the OPD is centred on the chief
	// image point C = (0, 0, z_g) and passes through that point, so its radius
	// is z_g − zQ. The wave model puts its pupil plane at zQ and its focal
	// plane at the distance d = R_sphere: the wave focus is then at zQ + d
	// before the defocus scan.
	zQ := sys.VertexZ(len(sys.Surfaces) - 1)
	if d := math.Abs(gl.SphereRadius - (zG - zQ)); d > 1e-9 {
		t.Fatalf("E4: reference sphere radius %.9f mm ≠ z_g − zQ = %.9f mm", gl.SphereRadius, zG-zQ)
	}
	zM := gl.SphereRadius * 1e-3
	dxM := pup.Radius / float64(lat.R) * 1e-3
	focus := lat.fraunhofer(dLineM, dxM, zM, opdToM(gl.OPD))
	ideal := lat.fraunhofer(dLineM, dxM, zM, nil)
	pIdeal := peak(ideal)

	dzW, pW := zScanPeak(focus, dLineM)
	strehlW := pW / pIdeal
	pixUM := focus.DX * 1e6
	zWaveAbs := zQ + zM*1e3 + dzW*1e-3 // mm, absolute wave best focus
	dzI, _ := zScanPeak(ideal, dLineM)

	// Geometric Strehl at z_g for context.
	qp, qw := diskQuadrature(64, 128)
	gq := rayOPD(t, sys, pup, qp, zG)
	sExact, _ := integralStrehl(wavesOf(gq.OPD), qw)

	t.Logf("E4 geometric best focus z_g = %.7f mm (RMS %.6f λ); geometric Strehl at z_g: exact %.6f, Maréchal %.6f",
		zG, rmsG, sExact, marechal(rmsG))
	t.Logf("E4 wave best focus: z_w = %.7f mm (zQ + d + Δz_w, Δz_w = %+.4f µm, peak ratio %.6f); z_w − z_g = %+.4f µm; pixel p = %.4f µm; |z_w − z_g|/p = %.4f; 2F#·p = %.3f µm",
		zWaveAbs, dzW, strehlW, (zWaveAbs-zG)*1e3, pixUM, math.Abs(zWaveAbs-zG)*1e3/pixUM, 2*(zM/(2*pup.Radius*1e-3))*pixUM)
	// Analytic control of the scan. On axis, a uniform aperture of radius a
	// focused to f has I(z) ∝ z⁻²·sinc²(u), u = k·a²·(1/z − 1/f)/4. Its maximum
	// satisfies cot u − 1/u = −4z/(k·a²), i.e. u ≈ 12z/(k·a²), which gives
	// z_peak − f ≈ −48·f³/(k²·a⁴) (the diffraction focal shift; no aberration).
	k := 2 * math.Pi / dLineM
	a := pup.Radius * 1e-3
	shiftAna := -48 * zM * zM * zM / (k * k * a * a * a * a) * 1e6 // µm
	t.Logf("E4 control (unaberrated pupil): scan peak at %+.4f µm; analytic diffraction focal shift −48f³/(k²a⁴) = %+.4f µm (|Δz|/p = %.4f)",
		dzI, shiftAna, math.Abs(dzI)/pixUM)
	if math.Abs(dzI-shiftAna) > 0.05*math.Abs(shiftAna) {
		t.Errorf("E4: unaberrated scan peak %+.4f µm vs analytic focal shift %+.4f µm (tolerance 5%%)", dzI, shiftAna)
	}
	if math.Abs(zWaveAbs-zG)*1e3 > pixUM {
		t.Errorf("E4: wave best focus z_w = %.7f mm differs from geometric z_g = %.7f mm by %.4f µm (> one pixel, %.4f µm)",
			zWaveAbs, zG, (zWaveAbs-zG)*1e3, pixUM)
	}
	if strehlW < 0.98 {
		t.Errorf("E4: wave Strehl %.6f < 0.98", strehlW)
	}
}
