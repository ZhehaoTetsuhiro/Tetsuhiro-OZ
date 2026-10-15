package wave

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
)

// naiveDFT is the O(n²) reference transform. The twiddle exp(±2πi·m/n) is
// tabulated for m = (k·j) mod n, so every angle lies in [0, 2π) and carries
// no argument-reduction error (k·j/n can reach ~10³ turns for n = 1024).
func naiveDFT(a []complex128, inverse bool) []complex128 {
	n := len(a)
	out := make([]complex128, n)
	sign := -1.0
	if inverse {
		sign = 1.0
	}
	tw := make([]complex128, n)
	for m := range tw {
		tw[m] = cmplx.Exp(complex(0, sign*2*math.Pi*float64(m)/float64(n)))
	}
	for k := 0; k < n; k++ {
		var s complex128
		for j := 0; j < n; j++ {
			s += a[j] * tw[(k*j)%n]
		}
		if inverse {
			s /= complex(float64(n), 0)
		}
		out[k] = s
	}
	return out
}

// maxRelErr is max_k |got_k − want_k| / max_k |want_k|: the error relative to
// the largest reference coefficient (per-bin relative errors are ill-defined
// where the spectrum is near zero).
func maxRelErr(got, want []complex128) float64 {
	var d, scale float64
	for i := range want {
		if v := cmplx.Abs(got[i] - want[i]); v > d {
			d = v
		}
		if v := cmplx.Abs(want[i]); v > scale {
			scale = v
		}
	}
	return d / scale
}

// TestFFT1DNaiveDFTPowersOfTwo checks the transform against the naive DFT for
// N = 8, 16, …, 1024 on a fixed-seed random complex input, in both
// directions. Acceptance: max relative error ≤ 1e-12.
func TestFFT1DNaiveDFTPowersOfTwo(t *testing.T) {
	rng := rand.New(rand.NewSource(20261015)) // fixed seed: reproducible input
	worst := 0.0
	for n := 8; n <= 1024; n <<= 1 {
		a := make([]complex128, n)
		for i := range a {
			a[i] = complex(2*rng.Float64()-1, 2*rng.Float64()-1)
		}
		ef := maxRelErr(FFT1D(a, false), naiveDFT(a, false))
		ei := maxRelErr(FFT1D(a, true), naiveDFT(a, true))
		t.Logf("N=%4d: max relative error forward %.3e, inverse %.3e", n, ef, ei)
		if ef > 1e-12 {
			t.Errorf("N=%d forward: max relative error %.3g > 1e-12", n, ef)
		}
		if ei > 1e-12 {
			t.Errorf("N=%d inverse: max relative error %.3g > 1e-12", n, ei)
		}
		worst = math.Max(worst, math.Max(ef, ei))
	}
	t.Logf("worst over N = 8…1024: %.3e (threshold 1e-12)", worst)
}

func TestFFT1DMatchesDFT(t *testing.T) {
	for _, n := range []int{2, 4, 8, 16, 64} {
		a := make([]complex128, n)
		for i := range a {
			a[i] = complex(math.Sin(float64(i)*0.7)+0.3, math.Cos(float64(i)*0.31)-0.2)
		}
		got := FFT1D(a, false)
		want := naiveDFT(a, false)
		for i := range got {
			if d := cmplx.Abs(got[i] - want[i]); d > 1e-12*float64(n) {
				t.Fatalf("n=%d forward[%d]: |diff| = %g", n, i, d)
			}
		}
		back := FFT1D(got, true)
		for i := range back {
			if d := cmplx.Abs(back[i] - a[i]); d > 1e-12 {
				t.Fatalf("n=%d round trip[%d]: |diff| = %g", n, i, d)
			}
		}
	}
}

func TestFFT2DMatchesSeparableAndRoundTrip(t *testing.T) {
	n := 16
	a := make([]complex128, n*n)
	for r := 0; r < n; r++ {
		for c := 0; c < n; c++ {
			a[r*n+c] = complex(math.Sin(float64(r)*0.3+float64(c)*0.11), math.Cos(float64(r*c)*0.07))
		}
	}
	// Reference: naive rows, then naive columns.
	want := make([]complex128, n*n)
	for r := 0; r < n; r++ {
		row := naiveDFT(a[r*n:(r+1)*n], false)
		copy(want[r*n:(r+1)*n], row)
	}
	for c := 0; c < n; c++ {
		col := make([]complex128, n)
		for r := 0; r < n; r++ {
			col[r] = want[r*n+c]
		}
		col = naiveDFT(col, false)
		for r := 0; r < n; r++ {
			want[r*n+c] = col[r]
		}
	}
	got := append([]complex128(nil), a...)
	FFT2D(got, n, false)
	var maxd float64
	for i := range got {
		if d := cmplx.Abs(got[i] - want[i]); d > maxd {
			maxd = d
		}
	}
	if maxd > 1e-11 {
		t.Fatalf("2-D forward: max |diff| = %g", maxd)
	}
	FFT2D(got, n, true)
	for i := range got {
		if d := cmplx.Abs(got[i] - a[i]); d > 1e-12 {
			t.Fatalf("2-D round trip[%d]: |diff| = %g", i, d)
		}
	}
}
