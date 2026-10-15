// Package doe implements the diffractive optical element (DOE) layer of
// Tetsuhiro-OZ (G8): focusing phase profiles quantised to N levels, the
// scalar diffraction model used to evaluate them, and the first-order
// diffraction efficiency compared with the analytic full-depth N-level value.
//
// Scalar model. A thin, lossless, phase-only element multiplies the incident
// field by t(x, y) = a(x, y)·exp(i φ_q(x, y)), where a is the amplitude
// (unity for a uniform beam) and φ_q the quantised phase. The field at a
// distance z is obtained with the band-limited angular spectrum method (ASM;
// K. Matsushima and T. Shimobaba, Opt. Express 17(22), 19662–19673 (2009),
// doi:10.1364/OE.17.019662) computed with a 2-D FFT:
//
//	U(fx, fy; z) = U0(fx, fy)·H(fx, fy; z),
//	H = exp(i 2π z sqrt(1/λ² − fx² − fy²))   for |fx| < f_lim, |fy| < f_lim,
//	    and 0 otherwise (evanescent part and band-limit removed),
//	f_lim = 1/(λ sqrt((2 Δf z)² + 1)),   Δf = 1/(M·dx).
//
// The band limit suppresses aliasing of the transfer function on the
// periodic M×M grid; the design-order focus lies well inside it.
//
// Grid convention: a field of M×M samples (M a power of two) with pitch dx
// has its origin at index 0 and sample x_j = dx·j for j in [−M/2, M/2), so no
// fftshift is needed. Lengths and λ use the same unit (µm in the tests).
//
// Analytic reference (full depth, N levels). For a phase that wraps through a
// full 2π and is quantised to N equal levels (phase bins of width 2π/N, each
// level transmitting with unit amplitude), the Fourier coefficients of the
// periodic (grating) profile are
//
//	c_m = (N/π m) sin(π m/N) exp(−iπ m/N)   for m ≡ 1 (mod N),   c_m = 0 otherwise,
//
// so the first-order efficiency is
//
//	η_N = |c_1|² = [sin(π/N)/(π/N)]² ,   η_2 = 4/π² ≈ 0.405285 .
//
// Validity: (i) phase depth exactly 2π (no residual modulation); (ii) equal
// level widths; (iii) unit amplitude transmittance; (iv) scalar theory, i.e.
// the local grating period is much larger than λ (the diffraction angles are
// small); (v) the sampling resolves the level steps (several pixels per step)
// for the numerical check; (vi) the value is the diffraction efficiency of the
// first order alone, which for the grating means the energy in order m = 1
// with respect to the incident power. For a focusing (kinoform) element the
// same η_N applies locally to every zone, because the local period changes
// slowly over the aperture; the other orders m = 1 + jN carry the energy
// η_N/(1 + jN)² and focus elsewhere (they are defocused at the design focus).
// The sum over all orders equals 1 (Parseval), which the tests check.
//
// Literature. The derivation above is the Fourier series of the N-level
// staircase and is self-contained. Continuous (sawtooth) kinoform efficiencies
// of the form sinc²(·) are discussed in the applied-optics literature (e.g.
// "Optical performance of holographic kinoforms", Applied Optics 28, 976 —
// URL https://www.apollooptical.com/hubfs/articles/ao28_976.pdf, read for this
// revision); the N-level formula itself was not found in a retrieved source and
// is 待核对 against a primary reference.
package doe
