package optics

// Exported entry points to the kernel's FFT for other packages (the geometric
// → wave bridge in package wave). These are thin wrappers: they dispatch
// exactly like the package-private entry points (cuFFT/MUSA when the GPU
// backend is enabled, the pure-Go transform otherwise) and keep their
// normalisation, so there is a single FFT implementation in the repository.

// FFT1D transforms a in place. The forward transform is unnormalised,
//
//	X_k = Σ_j x_j · exp(−2πi·j·k/n),
//
// and the inverse applies the 1/n factor, so FFT1D(FFT1D(x, false), true) == x.
// Power-of-two lengths use radix-2; other lengths use Bluestein's algorithm.
func FFT1D(a []complex128, inverse bool) {
	fft1DAny(a, inverse)
}

// FFT2D transforms an n×n row-major array in place with the same convention as
// FFT1D applied to every row and then every column (the inverse carries 1/n²).
// len(a) must equal n*n.
func FFT2D(a []complex128, n int, inverse bool) {
	fft2D(a, n, inverse)
}
