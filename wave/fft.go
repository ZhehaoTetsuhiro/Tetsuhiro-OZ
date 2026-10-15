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
// cycles/mm. The FFT is the wave-optics kernel's own transform (package
// optics), so the bridge and the kernel share one implementation; sizes here
// are powers of two.
package wave

import "tetsuhiro-oz/optics"

// isPowerOfTwo reports whether n is a positive power of two.
func isPowerOfTwo(n int) bool { return n > 0 && n&(n-1) == 0 }

// FFT1D returns the forward or inverse transform of a (a copy; the input is
// not modified). The forward transform is unnormalised, X_k = Σ_j a_j·
// exp(−2πi·jk/n); the inverse carries the 1/n factor, so
// FFT1D(FFT1D(a, false), true) == a. Panics for non-power-of-two lengths (a
// programming error the callers in this package exclude by construction).
func FFT1D(a []complex128, inverse bool) []complex128 {
	if !isPowerOfTwo(len(a)) {
		panic("wave: FFT length must be a power of two")
	}
	out := append([]complex128(nil), a...)
	optics.FFT1D(out, inverse)
	return out
}

// FFT2D transforms an n×n row-major array in place. n must be a power of two.
// The forward transform maps f(x,y) -> sum f·exp(-i·2π·(fx·x + fy·y)/n); the
// inverse carries 1/n² (1/n per dimension).
func FFT2D(a []complex128, n int, inverse bool) {
	if !isPowerOfTwo(n) || len(a) != n*n {
		panic("wave: FFT2D needs an n×n array with n a power of two")
	}
	optics.FFT2D(a, n, inverse)
}
