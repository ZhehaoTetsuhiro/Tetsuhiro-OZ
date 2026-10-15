package doe

import (
	"fmt"
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
)

// LensPhase is the target phase (rad, in [0, 2π)) of an ideal focusing
// element of focal length f at transverse position (x, y) for wavelength
// lambda (same length unit):
//
//	φ(x, y) = −(2π/λ) (sqrt(x² + y² + f²) − f)   mod 2π .
//
// The sign makes the element converging for a plane wave incident from the
// −z side: the phase retarded by the element decreases away from the axis.
func LensPhase(x, y, f, lambda float64) float64 {
	phi := -2 * math.Pi / lambda * (math.Sqrt(x*x+y*y+f*f) - f)
	return wrap2Pi(phi)
}

// wrap2Pi maps a phase to [0, 2π).
func wrap2Pi(phi float64) float64 {
	w := math.Mod(phi, 2*math.Pi)
	if w < 0 {
		w += 2 * math.Pi
	}
	if w >= 2*math.Pi { // guard against rounding up to 2π
		w = 0
	}
	return w
}

// Quantize maps a phase to one of n equally spaced levels 2πk/n, k = 0..n−1,
// using the bin [2πk/n, 2π(k+1)/n) (floor). A constant offset of the level set
// multiplies the first-order amplitude by a unit phase and does not change the
// efficiency, so the floor convention is as good as rounding; the floor keeps
// the bin boundaries at 2πk/n, which lie on sample points whenever n divides
// the grid, so the sampled grating is exact (see GratingEtaFFT).
func Quantize(phi float64, n int) (float64, error) {
	if n < 2 {
		return 0, fmt.Errorf("doe: quantisation needs at least 2 levels (n=%d)", n)
	}
	w := wrap2Pi(phi)
	k := int(math.Floor(w * float64(n) / (2 * math.Pi)))
	if k >= n { // w is just below 2π after rounding
		k = n - 1
	}
	return 2 * math.Pi * float64(k) / float64(n), nil
}

// LensProfile holds the complex transmittance of an N-level quantised lens on
// an m×m grid (row-major, origin at index 0, pitch dx).
type LensProfile struct {
	M      int
	Dx     float64
	Levels int
	T      []complex128 // t = a·exp(i φ_q), a = amplitude
}

// Coord returns the transverse coordinate of grid index i (0 ≤ i < M) in the
// convention x_j = dx·j for j in [−M/2, M/2).
func Coord(i, m int, dx float64) float64 {
	j := i
	if j >= m/2 {
		j -= m
	}
	return float64(j) * dx
}

// GaussianLens builds the transmittance of a Gaussian-apodised lens: amplitude
// a(x, y) = exp(−(x²+y²)/w0²) (intensity radius w0, power-normalised later)
// times exp(i φ_q), φ_q the N-level quantised LensPhase of focal length f.
// The grid is m×m with pitch dx; the lens centre is grid index (0, 0).
func GaussianLens(m int, dx, f, lambda, w0 float64, n int) (LensProfile, error) {
	if m < 2 || m&(m-1) != 0 {
		return LensProfile{}, fmt.Errorf("doe: grid size %d must be a power of two", m)
	}
	if !(dx > 0 && f > 0 && lambda > 0 && w0 > 0) {
		return LensProfile{}, fmt.Errorf("doe: dx, f, lambda, w0 must be positive")
	}
	if n < 2 {
		return LensProfile{}, fmt.Errorf("doe: need at least 2 levels")
	}
	t := make([]complex128, m*m)
	for r := 0; r < m; r++ {
		y := Coord(r, m, dx)
		for c := 0; c < m; c++ {
			x := Coord(c, m, dx)
			phi, err := Quantize(LensPhase(x, y, f, lambda), n)
			if err != nil {
				return LensProfile{}, err
			}
			amp := math.Exp(-(x*x + y*y) / (w0 * w0))
			s, co := math.Sincos(phi)
			t[r*m+c] = complex(amp*co, amp*s)
		}
	}
	return LensProfile{M: m, Dx: dx, Levels: n, T: t}, nil
}

// Power returns Σ|t|²·dx², the power transmitted by the element (the element
// is lossless, so this is the incident power on the aperture).
func (p LensProfile) Power() float64 {
	var s float64
	for _, v := range p.T {
		s += real(v)*real(v) + imag(v)*imag(v)
	}
	return s * p.Dx * p.Dx
}

// PhaseImage returns the quantised phase of p as an 8-bit grey image of
// the m×m grid, with the grid origin at the image centre (x to the right,
// y downwards; level k maps to the grey value round(255·k/N)). The image is
// only for inspection; the numerical model uses the complex transmittance.
func PhaseImage(p LensProfile) *image.Gray {
	m := p.M
	img := image.NewGray(image.Rect(0, 0, m, m))
	for r := 0; r < m; r++ {
		for c := 0; c < m; c++ {
			// Reorder so that the origin (grid index 0) is at the image centre.
			rr := (r + m/2) % m
			cc := (c + m/2) % m
			v := p.T[rr*m+cc]
			ph := wrap2Pi(math.Atan2(imag(v), real(v)))
			k := math.Round(ph * float64(p.Levels) / (2 * math.Pi))
			if int(k) >= p.Levels {
				k = 0
			}
			g := uint8(math.Round(255 * k / float64(p.Levels)))
			img.SetGray(c, r, color.Gray{Y: g})
		}
	}
	return img
}

// WritePhasePNG writes PhaseImage(p) to a PNG file.
func WritePhasePNG(path string, p LensProfile) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	if err := png.Encode(f, PhaseImage(p)); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}
