package wave

import (
	"math"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// syntheticGrid builds an OPD grid directly (no ray tracing): the wavefront
// function is evaluated in units of the wavelength (waves), the pupil disc
// is the unit disc.
func syntheticGrid(np int, lambdaNm float64, waves func(px, py float64) float64) *geo.OPDGrid {
	pts := SquareGrid(np)
	lamMM := lambdaNm * 1e-6
	opd := make([]float64, len(pts))
	for i, p := range pts {
		opd[i] = waves(p.PX, p.PY) * lamMM
	}
	return &geo.OPDGrid{
		Lambda: lambdaNm,
		Pupil:  geo.Pupil{Z: 0, Radius: 10}, // mm; F# = f/(2R)
		Pupils: pts,
		OPD:    opd,
	}
}

// TestPSFAiry checks the diffraction core against the analytic Airy result:
// the first minimum of |MF|² of a circular pupil sits at 1.2197·λ·F#, and
// the first bright ring carries ~1.75% of the peak.
func TestPSFAiry(t *testing.T) {
	const lambdaNm = 550.0
	const focal = 100.0 // mm; R = 10 mm -> F# = 5
	np, pad := 128, 4
	g := syntheticGrid(np, lambdaNm, func(px, py float64) float64 { return 0 })
	p, err := PSFFromGrid(g, focal, np, pad)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(p.Strehl-1) > 1e-12 {
		t.Fatalf("zero-OPD Strehl = %.15f, want 1", p.Strehl)
	}

	// First local minimum along the centre row, starting a few pixels out
	// so that the noise floor near the array edges cannot trigger first.
	half := p.N / 2
	row := p.Intensity[half*p.N : (half+1)*p.N]
	minIdx := -1
	for c := half + 3; c < p.N-1 && c < half+3*p.N/8; c++ {
		if row[c] < row[c-1] && row[c] <= row[c+1] {
			minIdx = c
			break
		}
	}
	if minIdx < 0 {
		t.Fatal("no minimum found along the centre row")
	}
	x := float64(minIdx-half) * p.PixelMM
	want := p.AiryRadiusMM()
	if math.Abs(x-want) > 0.9*p.PixelMM {
		t.Fatalf("Airy first minimum at %.4g mm, want %.4g mm (±%.2g)",
			x, want, 0.9*p.PixelMM)
	}

	// First bright ring energy ratio.
	ringMax := 0.0
	for r := 0; r < p.N; r++ {
		for c := 0; c < p.N; c++ {
			dx := float64(c-half) * p.PixelMM
			dy := float64(r-half) * p.PixelMM
			rad := math.Hypot(dx, dy)
			if rad > 1.5*want && rad < 2.6*want {
				if v := p.Intensity[r*p.N+c]; v > ringMax {
					ringMax = v
				}
			}
		}
	}
	if ratio := ringMax / p.Intensity[half*p.N+half]; ratio < 0.008 || ratio > 0.032 {
		t.Fatalf("first bright ring ratio = %.4f, want ~0.0175", ratio)
	}
}

// TestPSFMTFvsAiry checks the MTF of the same perfect aperture against the
// analytic aberration-free MTF (2/π)(acos x − x·sqrt(1−x²)) of the relative
// frequency x.
func TestPSFMTFvsAiry(t *testing.T) {
	const lambdaNm = 550.0
	const focal = 100.0
	np, pad := 128, 4
	g := syntheticGrid(np, lambdaNm, func(px, py float64) float64 { return 0 })
	p, err := PSFFromGrid(g, focal, np, pad)
	if err != nil {
		t.Fatal(err)
	}
	vc := p.CutoffCyclesPerMM()
	for _, x := range []float64{0.2, 0.4, 0.6} {
		got := p.MTFAt(x * vc)
		want := AiryMTF(x)
		if math.Abs(got-want) > 0.025 {
			t.Fatalf("MTF(%.1f·νc) = %.4f, want %.4f (analytic)", x, got, want)
		}
	}
	if v := p.MTFAt(0); math.Abs(v-1) > 1e-9 {
		t.Fatalf("MTF(0) = %.12f, want 1", v)
	}
	if v := p.MTFAt(1.05 * vc); v > 0.10 {
		t.Fatalf("MTF above cutoff = %.4f, want ~0 (clamped at Nyquist)", v)
	}
}

// TestPSFDefocusMarechal checks the Strehl ratio of a pure defocus against
// the Maréchal approximation exp(−(2πσ)²) with σ the RMS of W020·ρ² over
// the pupil disc, σ = W020/√12.
func TestPSFDefocusMarechal(t *testing.T) {
	const lambdaNm = 550.0
	const w020 = 0.25 // waves of defocus
	np, pad := 128, 4
	g := syntheticGrid(np, lambdaNm, func(px, py float64) float64 {
		return w020 * (px*px + py*py)
	})
	p, err := PSFFromGrid(g, 100, np, pad)
	if err != nil {
		t.Fatal(err)
	}
	sigma := w020 / math.Sqrt(12)
	marechal := math.Exp(-(2 * math.Pi * sigma) * (2 * math.Pi * sigma))
	if math.Abs(p.Strehl-marechal) > 0.015 {
		t.Fatalf("defocus Strehl = %.5f, Maréchal %.5f (σ = %.4f λ)", p.Strehl, marechal, sigma)
	}
	t.Logf("defocus: Strehl = %.5f, Maréchal = %.5f", p.Strehl, marechal)
}

// TestPSFRealLens is the geometric→wave bridge end to end: design the
// achromat with the lens package, trace its wavefront on a uniform pupil
// grid and compute the PSF; the d-line Strehl of a well-corrected f/5
// doublet must be high, with the peak at the chief ray.
func TestPSFRealLens(t *testing.T) {
	p, _, err := lens.DesignAchromat(100, "N-BK7", "F2", 10, nil)
	if err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	fld := geo.Field{}
	pup, err := sys.EntrancePupil(fld, geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	np, pad := 64, 2
	pts := SquareGrid(np)
	// Scan a small focus range for the wavefront's best focus (the paraxial
	// image plane of a real design carries a little residual defocus), then
	// compute the PSF there.
	bestZ, bestRMS := par.ImageZ, math.Inf(1)
	var g *geo.OPDGrid
	for k := -4; k <= 4; k++ {
		z := par.ImageZ + float64(k)*par.ImageZ*5e-4
		cand, err := sys.OPDGrid(fld, pup, pts, z, geo.Wavelengthd, geo.OPDOptions{})
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("focus scan k=%+d z=%.5f mm: RMS = %.6g λ", k, z, cand.RMSWaves())
		if cand.RMSWaves() < bestRMS {
			bestRMS, bestZ, g = cand.RMSWaves(), z, cand
		}
	}
	psf, err := PSFFromGrid(g, par.EFL, np, pad)
	if err != nil {
		t.Fatal(err)
	}
	t.Logf("real lens: best focus z = %.5f mm, OPD RMS = %.4g λ, PV = %.4g λ, Strehl = %.5f, F# = %.2f, dx = %.4g mm",
		bestZ, bestRMS, g.PVWaves(), psf.Strehl, psf.FNumber, psf.PixelMM)
	if psf.Strehl < 0.6 {
		t.Fatalf("Strehl = %.4f too low for the designed achromat at best focus (OPD RMS %.4g λ)", psf.Strehl, bestRMS)
	}
	peak, pr, pc := psf.Peak()
	if peak <= 0 {
		t.Fatal("empty PSF")
	}
	half := psf.N / 2
	if pr < half-2 || pr > half+2 || pc < half-2 || pc > half+2 {
		t.Fatalf("PSF peak at (%d,%d), want the grid centre (%d,%d)", pr, pc, half, half)
	}
	if m := psf.MTFAt(50); m < 0.1 || m > 0.999 {
		t.Fatalf("MTF(50 cycles/mm) = %.4f out of range", m)
	}
}

// TestPSFDeterministic pins bit-identical repeatability.
func TestPSFDeterministic(t *testing.T) {
	np, pad := 64, 2
	g := syntheticGrid(np, 550, func(px, py float64) float64 {
		return 0.1 * (px*px + py*py) * (px*px + py*py) // some spherical aberration
	})
	a, err := PSFFromGrid(g, 100, np, pad)
	if err != nil {
		t.Fatal(err)
	}
	b, err := PSFFromGrid(g, 100, np, pad)
	if err != nil {
		t.Fatal(err)
	}
	if a.Strehl != b.Strehl {
		t.Fatalf("Strehl differs: %g vs %g", a.Strehl, b.Strehl)
	}
	for i := range a.Intensity {
		if a.Intensity[i] != b.Intensity[i] {
			t.Fatalf("intensity differs at %d", i)
		}
	}
}
