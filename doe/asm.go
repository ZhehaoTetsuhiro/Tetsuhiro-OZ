package doe

import (
	"fmt"
	"math"
)

// signedIndex maps FFT output index k (0…m−1) to the signed frequency index.
func signedIndex(k, m int) int {
	if k >= m/2 {
		return k - m
	}
	return k
}

// asmTransfer returns the band-limited angular-spectrum transfer function H
// (see doc.go) on the m×m FFT grid for distance z. Entries outside the
// propagating, band-limited region are zero.
func asmTransfer(m int, dx, lambda, z float64) []complex128 {
	h := make([]complex128, m*m)
	df := 1 / (float64(m) * dx)
	// Band limit of the transfer function on this grid (Matsushima & Shimobaba).
	fLim := 1 / (lambda * math.Sqrt(math.Pow(2*df*z, 2)+1))
	inv2 := 1 / (lambda * lambda)
	for r := 0; r < m; r++ {
		fy := float64(signedIndex(r, m)) * df
		if math.Abs(fy) >= fLim {
			continue
		}
		for c := 0; c < m; c++ {
			fx := float64(signedIndex(c, m)) * df
			if math.Abs(fx) >= fLim {
				continue
			}
			arg := inv2 - fx*fx - fy*fy
			if arg <= 0 {
				continue // evanescent
			}
			s, co := math.Sincos(2 * math.Pi * z * math.Sqrt(arg))
			h[r*m+c] = complex(co, s)
		}
	}
	return h
}

// Propagate returns the complex field u (m×m, origin at index 0, pitch dx)
// after propagation over distance z in a homogeneous medium of wavelength
// lambda, by the band-limited angular spectrum method. u is not modified.
func Propagate(u []complex128, m int, dx, lambda, z float64) ([]complex128, error) {
	if len(u) != m*m {
		return nil, fmt.Errorf("doe: field has %d samples, want %d", len(u), m*m)
	}
	if !(dx > 0 && lambda > 0) {
		return nil, fmt.Errorf("doe: dx and lambda must be positive")
	}
	spec := append([]complex128(nil), u...)
	if err := FFT2(spec, m, false); err != nil {
		return nil, err
	}
	h := asmTransfer(m, dx, lambda, z)
	for i := range spec {
		spec[i] *= h[i]
	}
	if err := FFT2(spec, m, true); err != nil {
		return nil, err
	}
	return spec, nil
}

// Spectrum returns the unnormalised forward FFT of the field u (see FFT2),
// which AxialField needs.
func Spectrum(u []complex128, m int) ([]complex128, error) {
	if len(u) != m*m {
		return nil, fmt.Errorf("doe: field has %d samples, want %d", len(u), m*m)
	}
	spec := append([]complex128(nil), u...)
	if err := FFT2(spec, m, false); err != nil {
		return nil, err
	}
	return spec, nil
}

// AxialField returns the complex field on the optical axis, u(0, 0, z), for
// each distance in zs, from the spectrum of the input field (see Spectrum).
// Using u(0,0,z) = (1/M²) Σ_k U_k H_k(z), it needs no FFT per plane and is
// the same band-limited ASM as Propagate.
func AxialField(spec []complex128, m int, dx, lambda float64, zs []float64) ([]complex128, error) {
	if len(spec) != m*m {
		return nil, fmt.Errorf("doe: spectrum has %d samples, want %d", len(spec), m*m)
	}
	out := make([]complex128, len(zs))
	inv := 1 / float64(m*m)
	for i, z := range zs {
		h := asmTransfer(m, dx, lambda, z)
		var s complex128
		for k := range spec {
			s += spec[k] * h[k]
		}
		out[i] = s * complex(inv, 0)
	}
	return out, nil
}
