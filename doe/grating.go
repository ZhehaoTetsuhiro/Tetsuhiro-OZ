package doe

import (
	"fmt"
	"math"
)

// EtaN returns the analytic first-order efficiency of a full-depth N-level
// quantised lossless phase grating, η_N = [sin(π/N)/(π/N)]² (see doc.go for
// the validity conditions). N must be at least 2; η_2 = 4/π² ≈ 0.405285.
func EtaN(n int) (float64, error) {
	if n < 2 {
		return 0, fmt.Errorf("doe: EtaN needs N >= 2 (N=%d)", n)
	}
	x := math.Pi / float64(n)
	s := math.Sin(x) / x
	return s * s, nil
}

// GratingOrders returns the diffraction orders c_m, m = 0…M−1 in DFT order,
// of the one-dimensional full-depth N-level quantised sawtooth sampled on M
// points per period. The transmittance of sample j is exp(i·2π·k_j/N) with
// the integer level k_j = ⌊jN/M⌋ (exact when N divides M, so every bin edge
// lies on a sample). The orders are c_m = (1/M) Σ_j t_j exp(−2πi jm/M), i.e.
// the scalar Fraunhofer diffraction of the period, computed with the FFT.
func GratingOrders(n, m int) ([]complex128, error) {
	if n < 2 {
		return nil, fmt.Errorf("doe: grating needs N >= 2 (N=%d)", n)
	}
	if m < n || m%n != 0 {
		return nil, fmt.Errorf("doe: samples per period M=%d must be a multiple of N=%d", m, n)
	}
	t := make([]complex128, m)
	for j := 0; j < m; j++ {
		k := (j * n) / m // exact integer level
		ph := 2 * math.Pi * float64(k) / float64(n)
		s, c := math.Sincos(ph)
		t[j] = complex(c, s)
	}
	if err := fft1(t, false); err != nil {
		return nil, fmt.Errorf("doe: samples per period M=%d: %w", m, err)
	}
	inv := 1 / float64(m)
	for i := range t {
		t[i] *= complex(inv, 0)
	}
	return t, nil
}

// GratingEta returns the first-order efficiency |c_1|² and the total power in
// all orders Σ_m |c_m|² (1 for a lossless grating) of the sampled grating.
func GratingEta(n, m int) (eta1, total float64, err error) {
	c, err := GratingOrders(n, m)
	if err != nil {
		return 0, 0, err
	}
	for _, v := range c {
		total += real(v)*real(v) + imag(v)*imag(v)
	}
	eta1 = real(c[1])*real(c[1]) + imag(c[1])*imag(c[1])
	return eta1, total, nil
}

// GratingEtaDiscreteExact is the closed form of |c_1|² for the sampled grating
// of GratingEta (M a multiple of N): |c_1| = (N/M)·sin(π/N)/sin(π/M). It is
// used as an internal-consistency reference for the FFT value; as M → ∞ it
// tends to EtaN(N).
func GratingEtaDiscreteExact(n, m int) (float64, error) {
	if n < 2 || m < n || m%n != 0 {
		return 0, fmt.Errorf("doe: need N >= 2 and M a multiple of N")
	}
	a := float64(n) / float64(m) * math.Sin(math.Pi/float64(n)) / math.Sin(math.Pi/float64(m))
	return a * a, nil
}
