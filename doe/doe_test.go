package doe

import (
	"image/png"
	"math"
	"math/cmplx"
	"os"
	"path/filepath"
	"strconv"
	"testing"
)

// Physical parameters of the focusing tests (µm). N_F = w0²/(λf) = 32 gives a
// 1/e amplitude radius w0 = 450 µm on the DOE and a focal spot radius
// w_f = λf/(πw0) ≈ 4.48 µm. The grid is M = 2048 with pitch dx = 8w0/M ≈ 1.76 µm,
// so w_f/dx = M/(8πN_F) ≈ 2.55 (the tests print the exact values).
const (
	testLambda = 0.6328  // µm (He–Ne)
	testF      = 10000.0 // µm (10 mm)
	testNF     = 32.0    // Fresnel number w0²/(λf)
	testM      = 2048    // grid size (power of two)
	testWinW   = 2.5     // efficiency window radius in focal spot radii
	gratingM   = 1024    // samples per grating period in the FFT test
	effTol     = 1e-3    // acceptance threshold on |η − η_N| (brief)
)

// TestEtaNAnalytic checks the closed form η_N = [sin(π/N)/(π/N)]² at the
// values quoted in the brief and its limits.
func TestEtaNAnalytic(t *testing.T) {
	eta2, err := EtaN(2)
	if err != nil {
		t.Fatal(err)
	}
	if d := math.Abs(eta2 - 4/(math.Pi*math.Pi)); d > 1e-15 {
		t.Errorf("η_2 = %.15g, want 4/π² = %.15g", eta2, 4/(math.Pi*math.Pi))
	}
	if math.Abs(eta2-0.405) > 5e-4 {
		t.Errorf("η_2 = %.6f, brief quotes ≈ 0.405", eta2)
	}
	prev := 0.0
	for _, n := range []int{2, 4, 8, 16, 32, 64} {
		e, _ := EtaN(n)
		if e <= prev || e > 1 {
			t.Errorf("η_%d = %.6g is not increasing towards 1", n, e)
		}
		prev = e
	}
	if e, _ := EtaN(4096); math.Abs(e-1) > 1e-6 {
		t.Errorf("η_4096 = %.9f, want → 1", e)
	}
	if _, err := EtaN(1); err == nil {
		t.Error("EtaN(1) must be rejected (no modulation)")
	}
}

// TestQuantizeLevels: the quantiser returns exactly N distinct levels 2πk/N on
// a fine sweep of the phase, and maps the wrapped interval correctly.
func TestQuantizeLevels(t *testing.T) {
	for _, n := range []int{2, 4, 8, 16} {
		seen := map[int]bool{}
		for i := 0; i < 100000; i++ {
			phi := -20*math.Pi + float64(i)*40*math.Pi/100000 // covers several 2π periods, both signs
			q, err := Quantize(phi, n)
			if err != nil {
				t.Fatal(err)
			}
			k := q * float64(n) / (2 * math.Pi)
			ki := math.Round(k)
			if math.Abs(k-ki) > 1e-9 || ki < 0 || ki >= float64(n) {
				t.Fatalf("N=%d: level %.12g is not 2πk/N", n, q)
			}
			seen[int(ki)] = true
		}
		if len(seen) != n {
			t.Errorf("N=%d: %d distinct levels produced, want %d", n, len(seen), n)
		}
	}
	if q, _ := Quantize(math.Pi*1.2, 2); math.Abs(q-math.Pi) > 1e-12 {
		t.Errorf("Quantize(1.2π, 2) = %v, want π", q)
	}
	if _, err := Quantize(0.1, 1); err == nil {
		t.Error("Quantize with N < 2 must fail")
	}
}

// TestLensPhaseFormula checks the target phase against the defining formula
// (mod 2π) and its paraxial limit −k r²/(2f).
func TestLensPhaseFormula(t *testing.T) {
	k := 2 * math.Pi / testLambda
	for _, xy := range [][2]float64{{0, 0}, {300, -120}, {-2500, 1700}} {
		x, y := xy[0], xy[1]
		got := LensPhase(x, y, testF, testLambda)
		want := -k * (math.Sqrt(x*x+y*y+testF*testF) - testF)
		d := math.Mod(got-want, 2*math.Pi)
		if d > math.Pi {
			d -= 2 * math.Pi
		} else if d < -math.Pi {
			d += 2 * math.Pi
		}
		if math.Abs(d) > 1e-9 {
			t.Errorf("LensPhase(%g,%g) = %.12g, defining formula %.12g (mod 2π diff %.2g)", x, y, got, want, d)
		}
		if got < 0 || got >= 2*math.Pi {
			t.Errorf("LensPhase(%g,%g) = %g outside [0, 2π)", x, y, got)
		}
	}
	// Paraxial limit at r = 20 µm: exact − (−k r²/2f) ≈ k r⁴/(8f³) ≈ 1e-6 rad.
	x := 20.0
	par := wrap2Pi(-k * x * x / (2 * testF))
	ex := LensPhase(x, 0, testF, testLambda)
	d := math.Mod(ex-par, 2*math.Pi)
	if d > math.Pi {
		d -= 2 * math.Pi
	} else if d < -math.Pi {
		d += 2 * math.Pi
	}
	if math.Abs(d) > 1e-6 {
		t.Errorf("paraxial phase error %.3g rad at r=20 µm, want ≤ 1e-6", d)
	}
}

// TestFFT2RoundTripAndParseval checks the 2-D FFT against the inverse and
// Parseval's theorem on a deterministic test field.
func TestFFT2RoundTripAndParseval(t *testing.T) {
	const n = 64
	a := make([]complex128, n*n)
	for i := range a {
		x := float64(i%n) / n
		y := float64(i/n) / n
		a[i] = complex(math.Sin(7*x+3*y)+0.25*math.Cos(11*x*y), math.Cos(5*x-2*y))
	}
	orig := append([]complex128(nil), a...)
	var eIn float64
	for _, v := range orig {
		eIn += real(v)*real(v) + imag(v)*imag(v)
	}
	if err := FFT2(a, n, false); err != nil {
		t.Fatal(err)
	}
	var eOut float64
	for _, v := range a {
		eOut += real(v)*real(v) + imag(v)*imag(v)
	}
	if rel := math.Abs(eOut-float64(n*n)*eIn) / (float64(n*n) * eIn); rel > 1e-12 {
		t.Errorf("Parseval: rel error %.3g", rel)
	}
	if err := FFT2(a, n, true); err != nil {
		t.Fatal(err)
	}
	var maxErr float64
	for i := range a {
		maxErr = math.Max(maxErr, cmplx.Abs(a[i]-orig[i]))
	}
	if maxErr > 1e-12 {
		t.Errorf("FFT2 round trip max error %.3g, want ≤ 1e-12", maxErr)
	}
	if err := FFT2(make([]complex128, 12), 3, false); err == nil {
		t.Error("FFT2 must reject a non power-of-two size")
	}
}

// TestASMRoundTripAndEnergy: a Gaussian beam propagated by +z and back by −z
// returns to itself, and the power is conserved (the band-limit and the
// evanescent cut remove nothing for this beam width).
func TestASMRoundTripAndEnergy(t *testing.T) {
	const m = 256
	const dx = 2.0
	const w = 60.0
	u := make([]complex128, m*m)
	for r := 0; r < m; r++ {
		y := Coord(r, m, dx)
		for c := 0; c < m; c++ {
			x := Coord(c, m, dx)
			u[r*m+c] = complex(math.Exp(-(x*x+y*y)/(w*w)), 0)
		}
	}
	const z = 2000.0
	fwd, err := Propagate(u, m, dx, testLambda, z)
	if err != nil {
		t.Fatal(err)
	}
	back, err := Propagate(fwd, m, dx, testLambda, -z)
	if err != nil {
		t.Fatal(err)
	}
	var eu, ef, maxErr float64
	for i := range u {
		eu += real(u[i])*real(u[i]) + imag(u[i])*imag(u[i])
		ef += real(fwd[i])*real(fwd[i]) + imag(fwd[i])*imag(fwd[i])
		maxErr = math.Max(maxErr, cmplx.Abs(back[i]-u[i]))
	}
	if rel := math.Abs(ef-eu) / eu; rel > 1e-9 {
		t.Errorf("power after 2 mm: relative change %.3g, want ≤ 1e-9", rel)
	}
	if maxErr > 1e-9 {
		t.Errorf("ASM round trip (+z then −z) max error %.3g, want ≤ 1e-9", maxErr)
	}
}

// TestAxialFieldMatchesPropagate: the on-axis field from the spectral sum
// equals the centre sample of the propagated field.
func TestAxialFieldMatchesPropagate(t *testing.T) {
	const m = 128
	const dx = 4.0
	lp, err := GaussianLens(m, dx, testF, testLambda, 120, 8)
	if err != nil {
		t.Fatal(err)
	}
	const z = 9000.0
	field, err := Propagate(lp.T, m, dx, testLambda, z)
	if err != nil {
		t.Fatal(err)
	}
	spec, err := Spectrum(lp.T, m)
	if err != nil {
		t.Fatal(err)
	}
	ax, err := AxialField(spec, m, dx, testLambda, []float64{z})
	if err != nil {
		t.Fatal(err)
	}
	peak := cmplx.Abs(field[0])
	if d := cmplx.Abs(ax[0]-field[0]) / peak; d > 1e-9 {
		t.Errorf("axial field %v vs propagated centre %v (rel %.3g)", ax[0], field[0], d)
	}
}

// TestGratingFirstOrderFFT is the primary numerical criterion of the brief:
// the first-order efficiency of the full-depth N-level quantised phase grating,
// computed by the scalar Fraunhofer model (FFT of the transmittance over one
// period of M = 1024 samples), against η_N, tolerance 1e-3. The sampling is
// exact (N divides M), so the discretisation error is only (π/M)²/3 ≈ 3e-6.
func TestGratingFirstOrderFFT(t *testing.T) {
	for _, n := range []int{2, 4, 8, 16} {
		eta, total, err := GratingEta(n, gratingM)
		if err != nil {
			t.Fatal(err)
		}
		want, _ := EtaN(n)
		if d := math.Abs(eta - want); d > effTol {
			t.Errorf("N=%d: FFT η_1 = %.9f, analytic %.9f (|diff| %.3g > %g)", n, eta, want, d, effTol)
		}
		exact, _ := GratingEtaDiscreteExact(n, gratingM)
		if d := math.Abs(eta - exact); d > 1e-12 {
			t.Errorf("N=%d: FFT η_1 %.15g vs discrete closed form %.15g", n, eta, exact)
		}
		if d := math.Abs(total - 1); d > 1e-12 {
			t.Errorf("N=%d: Σ_m |c_m|² = %.15g, want 1 (lossless)", n, total)
		}
		t.Logf("N=%d M=%d: η_1(FFT)=%.9f  η_N=%.9f  |diff|=%.2e", n, gratingM, eta, want, math.Abs(eta-want))
	}
}

// TestFocusDOE simulates the Gaussian-apodised focusing DOE with the band-
// limited ASM on a 2048² grid (dx ≈ 1.76 µm, N_F = 32) and checks, for each N:
//   - the first-order efficiency in the focal window against η_N (≤ 1e-3);
//   - the axial intensity peak against the design focal length f (≤ 1 pixel);
//   - the transverse intensity peak against the axis (≤ 1 pixel).
func TestFocusDOE(t *testing.T) {
	w0 := math.Sqrt(testNF * testLambda * testF)
	dx := 8 * w0 / testM
	for _, n := range []int{2, 4, 8, 16} {
		res, err := SimulateFocus(FocusSpec{
			Lambda: testLambda, F: testF, Levels: n, W0: w0, M: testM, Dx: dx, WinW: testWinW,
		})
		if err != nil {
			t.Fatal(err)
		}
		if d := math.Abs(res.Eta - res.EtaAnalytic); d > effTol {
			t.Errorf("N=%d: focal η = %.6f, η_N = %.6f (|diff| %.3g > %g)", n, res.Eta, res.EtaAnalytic, d, effTol)
		}
		if d := math.Abs(res.AxialZ - testF); d > dx {
			t.Errorf("N=%d: axial peak %.3f µm, f = %.1f µm (|diff| %.3g µm > 1 px = %.3g µm)", n, res.AxialZ, testF, d, dx)
		}
		if res.PeakOffsetPx > 1 {
			t.Errorf("N=%d: transverse peak %.2f px off axis", n, res.PeakOffsetPx)
		}
		t.Logf("N=%d: η=%.6f η_N=%.6f |diff|=%.2e; axial peak %.3f µm (f %+.3f µm = %+.2f px); transverse offset %.2f px; w_f=%.2f µm = %.2f px; P_in=%.4g",
			n, res.Eta, res.EtaAnalytic, math.Abs(res.Eta-res.EtaAnalytic),
			res.AxialZ, res.AxialZ-testF, (res.AxialZ-testF)/dx, res.PeakOffsetPx, res.Wf, res.Wf/dx, res.PowerIn)
	}
}

// TestPhaseMapOutput writes the quantised phase maps of the lens (N = 2, 4,
// 8, 16) as PNG files when DOE_PHASE_DIR is set; it always checks that the
// image decodes with exactly N grey levels.
func TestPhaseMapOutput(t *testing.T) {
	w0 := math.Sqrt(testNF * testLambda * testF)
	dx := 8 * w0 / 512
	dir := os.Getenv("DOE_PHASE_DIR")
	for _, n := range []int{2, 4, 8, 16} {
		lp, err := GaussianLens(512, dx, testF, testLambda, w0, n)
		if err != nil {
			t.Fatal(err)
		}
		path := filepath.Join(t.TempDir(), "phase.png")
		if dir != "" {
			path = filepath.Join(dir, "lens_phase_N"+strconv.Itoa(n)+".png")
		}
		if err := WritePhasePNG(path, lp); err != nil {
			t.Fatal(err)
		}
		f, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		img, err := png.Decode(f)
		f.Close()
		if err != nil {
			t.Fatal(err)
		}
		levels := map[uint32]bool{}
		b := img.Bounds()
		for y := b.Min.Y; y < b.Max.Y; y++ {
			for x := b.Min.X; x < b.Max.X; x++ {
				r, _, _, _ := img.At(x, y).RGBA()
				levels[r] = true
			}
		}
		if len(levels) != n {
			t.Errorf("N=%d: phase image has %d grey levels, want %d", n, len(levels), n)
		}
	}
}
