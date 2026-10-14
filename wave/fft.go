// Package wave adds the diffraction layer of Tetsuhiro-OZ: it turns the
// wavefront (optical path difference) of a lens system — as computed by the
// geo package's ray tracing — into the diffraction point-spread function,
// the Strehl ratio and the modulation transfer function.
//
// This is the "geometric → wave" bridge that closes the loop with the
// wave-optics kernel of Tetsuhiro-WOS: the geometric layer answers where the
// rays go, this package answers what the diffraction pattern at the focus
// looks like, and the two are validated against each other and against the
// analytic Airy/MTF results in the test suite.
//
// Conventions: lengths in mm (matching geo/lens), wavelengths in nm where
// carried from geo, but every function here takes lengths and wavelengths in
// the same unit (millimetre) unless a name says otherwise; frequencies in
// cycles/mm. The FFT is a compact, dependency-free, deterministic
// implementation (radix-2, sizes are powers of two).
package wave

import "math"

// fft1D transforms a in place (radix-2 Cooley–Tukey). len(a) must be a power
// of two. inverse applies the 1/N scaling, so fft(fft(x)) == x.
func fft1D(a []complex128, inverse bool) {
	n := len(a)
	if n <= 1 {
		return
	}
	// Bit-reversal permutation.
	for i, j := 1, 0; i < n; i++ {
		bit := n >> 1
		for ; j&bit != 0; bit >>= 1 {
			j ^= bit
		}
		j |= bit
		if i < j {
			a[i], a[j] = a[j], a[i]
		}
	}
	for size := 2; size <= n; size <<= 1 {
		half := size >> 1
		wm := -math.Pi / float64(half) // forward: exp(-2πi·j/size)
		if inverse {
			wm = -wm
		}
		for i := 0; i < n; i += size {
			wr, wi := 1.0, 0.0
			cr, ci := math.Cos(wm), math.Sin(wm)
			for j := 0; j < half; j++ {
				ur, ui := real(a[i+j]), imag(a[i+j])
				vr, vi := real(a[i+j+half]), imag(a[i+j+half])
				tr, ti := vr*wr-vi*wi, vr*wi+vi*wr
				a[i+j] = complex(ur+tr, ui+ti)
				a[i+j+half] = complex(ur-tr, ui-ti)
				wr, wi = wr*cr-wi*ci, wr*ci+wi*cr
			}
		}
	}
	if inverse {
		s := 1 / float64(n)
		for i := range a {
			a[i] = complex(real(a[i])*s, imag(a[i])*s)
		}
	}
}

// isPowerOfTwo reports whether n is a positive power of two.
func isPowerOfTwo(n int) bool { return n > 0 && n&(n-1) == 0 }

// FFT1D returns the forward or inverse transform of a (a copy; the input is
// not modified). Panics for non-power-of-two lengths (a programming error the
// callers in this package exclude by construction).
func FFT1D(a []complex128, inverse bool) []complex128 {
	if !isPowerOfTwo(len(a)) {
		panic("wave: FFT length must be a power of two")
	}
	out := append([]complex128(nil), a...)
	fft1D(out, inverse)
	return out
}

// FFT2D transforms an n×n row-major array in place. n must be a power of two.
// The forward transform maps f(x,y) -> sum f·exp(-i·2π·(fx·x + fy·y)/n).
func FFT2D(a []complex128, n int, inverse bool) {
	if !isPowerOfTwo(n) || len(a) != n*n {
		panic("wave: FFT2D needs an n×n array with n a power of two")
	}
	// Rows.
	for r := 0; r < n; r++ {
		fft1D(a[r*n:(r+1)*n], inverse)
	}
	// Columns.
	col := make([]complex128, n)
	for c := 0; c < n; c++ {
		for r := 0; r < n; r++ {
			col[r] = a[r*n+c]
		}
		fft1D(col, inverse)
		for r := 0; r < n; r++ {
			a[r*n+c] = col[r]
		}
	}
}
