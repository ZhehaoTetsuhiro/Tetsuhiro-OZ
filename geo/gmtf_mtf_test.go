package geo

import (
	"math"
	"testing"
)

// gmtfUniformBundle returns the LSF of the uniform circular bundle of the
// geometric MTF test: nRays rays on the equal-area disc, aimed at an
// aberration-free focus at z=0 from a pupil of radius a at distance L, observed
// on a plane a distance delta beyond the focus, binned into nBins bins over
// [−span·R, span·R). The footprint radius is R = a·delta/L.
func gmtfUniformBundle(t *testing.T, nRays, nBins int, delta, span float64) (LSF, float64) {
	t.Helper()
	const L, a = 100.0, 5.0
	R := a * delta / L
	det, err := NewSystem(Air{}, []Surface{{Profile: Plane{}, Name: "detector"}}, nil, nil)
	if err != nil {
		t.Fatalf("detector system: %v", err)
	}
	det.OriginZ = delta
	hits, err := det.ConvergingBundleHits(0, L, a, EqualAreaDisc(nRays), Wavelengthd)
	if err != nil {
		t.Fatalf("ConvergingBundleHits: %v", err)
	}
	x := make([]float64, len(hits))
	for i, h := range hits {
		x[i] = h.X
	}
	lsf, err := BinLSF(x, -span*R, span*R, nBins)
	if err != nil {
		t.Fatalf("BinLSF: %v", err)
	}
	return lsf, R
}

// TestGeometricMTFUniformBundleAnalytic is the analytic validation of the
// geometric MTF pipeline (rays → LSF → MTF) against the uniform circular
// bundle, whose MTF is |2·J1(x)/x| with x = 2π·f·R (math.J1).
//
// Configuration: pupil radius a = 5 mm at L = 100 mm before an aberration-free
// focus, observed at delta = 2 mm beyond it, so R = 0.1 mm. Rays: 10^6 on an
// equal-area disc (EqualAreaDisc). Bins: 210 over [−1.05R, 1.05R), so the bin
// width is Δ = R/100 = 0.001 mm and the Nyquist cutoff 1/(2Δ) = 500 cycles/mm.
// Frequencies: 2000 points in 0 < f ≤ 1/(2Δ).
//
// Errors:
//   - absolute error max |M_num − M_ana| over the whole range; M(0) = 1, so
//     this is the error relative to the full scale;
//   - pointwise relative error |M_num − M_ana|/M_ana where M_ana ≥ 0.1, the
//     range where a pointwise relative error is meaningful (x ≲ 3.0).
//
// The first zero f0 = j₁,₁/(2πR) = 6.098 cycles/mm is also checked: the
// numerical MTF must vanish there to within the absolute tolerance.
//
// Thresholds: both errors ≤ 1 %, and |M_num(f0)| ≤ 1 %. The measured values
// are in the test log.
func TestGeometricMTFUniformBundleAnalytic(t *testing.T) {
	const (
		nRays  = 1000000
		nBins  = 210
		delta  = 2.0
		span   = 1.05
		nFreqs = 2000
		tolAbs = 1e-2
		tolRel = 1e-2
		minRef = 0.1
	)
	lsf, R := gmtfUniformBundle(t, nRays, nBins, delta, span)
	fN := 1 / (2 * lsf.Width)
	freqs := make([]float64, nFreqs)
	for k := range freqs {
		freqs[k] = fN * float64(k+1) / nFreqs
	}
	mNum := lsf.MTF(freqs)
	maxAbs, fAbs := 0.0, 0.0
	maxRel, fRel := 0.0, 0.0
	for k, f := range freqs {
		ref := DiscMTF(f, R)
		if d := math.Abs(mNum[k] - ref); d > maxAbs {
			maxAbs, fAbs = d, f
		}
		if ref >= minRef {
			if d := math.Abs(mNum[k]-ref) / ref; d > maxRel {
				maxRel, fRel = d, f
			}
		}
	}
	f0 := 3.8317059702075125 / (2 * math.Pi * R) // first zero of J1
	mZero := lsf.MTF([]float64{f0})[0]
	t.Logf("uniform bundle: N=%d rays, %d bins, R=%.6g mm, Δ=%.6g mm (R/%.0f), Nyquist 1/(2Δ)=%.6g cycles/mm",
		lsf.N, nBins, R, lsf.Width, R/lsf.Width, fN)
	t.Logf("M(0)=%.15g; max|M_num-M_ana|=%.3e at f=%.4g (full scale); max relative error (M_ana>=%.1f)=%.3e at f=%.4g",
		lsf.MTF([]float64{0})[0], maxAbs, fAbs, minRef, maxRel, fRel)
	t.Logf("first zero f0=%.6g cycles/mm: |M_num(f0)|=%.3e", f0, mZero)
	if d := math.Abs(lsf.MTF([]float64{0})[0] - 1); d > 1e-12 {
		t.Errorf("M(0) differs from 1 by %.3e", d)
	}
	if maxAbs > tolAbs {
		t.Errorf("max |M_num-M_ana| = %.3e > %.0e over 0 < f <= %.4g cycles/mm", maxAbs, tolAbs, fN)
	}
	if maxRel > tolRel {
		t.Errorf("max relative error %.3e > %.0e where M_ana >= %.1f", maxRel, tolRel, minRef)
	}
	if mZero > tolAbs {
		t.Errorf("|M_num(f0)| = %.3e > %.0e at the first zero", mZero, tolAbs)
	}
}

// TestGeometricMTFPerfectFocusIsUnity checks the aberration-free limit. With
// delta = 0 every ray reaches the focal point to round-off, so the spot is a
// point and its MTF is 1 at every frequency (the analytic limit of 2·J1(x)/x
// as R → 0).
func TestGeometricMTFPerfectFocusIsUnity(t *testing.T) {
	det, err := NewSystem(Air{}, []Surface{{Profile: Plane{}, Name: "detector"}}, nil, nil)
	if err != nil {
		t.Fatal(err)
	}
	hits, err := det.ConvergingBundleHits(0, 100, 5, EqualAreaDisc(20000), Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	worst := 0.0
	x := make([]float64, len(hits))
	for i, h := range hits {
		worst = math.Max(worst, math.Hypot(h.X, h.Y))
		x[i] = h.X
	}
	if worst > 1e-12 {
		t.Fatalf("perfect focus: worst ray distance from the focus %.3e mm > 1e-12", worst)
	}
	// The spot lies within the single bin that straddles the focus at its
	// centre, so the MTF is the modulus of one unit phasor: 1 for every f.
	lsf, err := BinLSF(x, -0.5, 0.5, 3)
	if err != nil {
		t.Fatal(err)
	}
	for _, f := range []float64{1, 10, 50, 500} {
		if m := lsf.MTF([]float64{f})[0]; math.Abs(m-1) > 1e-12 {
			t.Errorf("point spot MTF(%g) = %.15g, want 1", f, m)
		}
	}
	t.Logf("perfect focus: worst ray distance %.3e mm over %d rays", worst, len(hits))
}

// TestBinLSFRejectsOutOfRange documents the contract of BinLSF: a coordinate
// outside the binned range is an error, not a silently dropped sample.
func TestBinLSFRejectsOutOfRange(t *testing.T) {
	if _, err := BinLSF([]float64{0, 0.5, 1.5}, -1, 1, 10); err == nil {
		t.Fatal("ray at 1.5 outside [-1,1): no error")
	}
	if _, err := BinLSF(nil, -1, 1, 10); err == nil {
		t.Fatal("empty ray set accepted")
	}
	if _, err := BinLSF([]float64{0}, 1, 1, 10); err == nil {
		t.Fatal("empty range accepted")
	}
}
