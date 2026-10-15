package doe

import (
	"fmt"
	"math"
)

// FocusSpec describes a focusing N-level DOE simulated with the scalar model.
// All lengths share one unit (the tests use µm).
type FocusSpec struct {
	Lambda float64 // wavelength
	F      float64 // design focal length
	Levels int     // number of phase levels N (full 2π depth)
	W0     float64 // 1/e amplitude radius of the Gaussian illumination on the DOE
	M      int     // grid size (power of two)
	Dx     float64 // pixel pitch of the DOE and of the focal plane
	WinW   float64 // efficiency window radius, in units of the focal spot radius
}

// FocusResult is the outcome of SimulateFocus.
type FocusResult struct {
	Eta          float64 // first-order efficiency: power in the focal window / incident power
	EtaAnalytic  float64 // EtaN(Levels)
	Wf           float64 // focal spot radius λF/(πW0) (1/e intensity)
	PeakOffsetPx float64 // distance of the focal-plane intensity maximum from the axis (pixels)
	AxialZ       float64 // axial intensity maximum, parabolic interpolation (length unit)
	AxialStep    float64 // axial scan pitch used for AxialZ (the axial "pixel")
	PowerIn      float64 // incident power Σ|t|²dx²
	WindowFrac   float64 // fraction of the design-order power inside the window, 1 − exp(−2 (WinW)²) for an ideal Gaussian focus
}

// SimulateFocus runs the scalar model of a Gaussian-apodised N-level DOE:
//
//  1. the quantised lens transmittance on an M×M grid (GaussianLens);
//  2. band-limited ASM propagation to the design focal plane z = F; the
//     first-order efficiency is the power inside a window of radius
//     WinW·w_f around the axis divided by the incident power;
//  3. the on-axis intensity, from the spectrum of the lens field: a coarse
//     scan over ±½ Rayleigh range, then a one-pixel (dx) scan around its
//     maximum; the axial peak AxialZ is refined by parabolic interpolation.
//
// The analytic reference is EtaN(Levels) (see doc.go).
func SimulateFocus(s FocusSpec) (FocusResult, error) {
	if !(s.WinW > 0) {
		return FocusResult{}, fmt.Errorf("doe: window radius must be positive")
	}
	lp, err := GaussianLens(s.M, s.Dx, s.F, s.Lambda, s.W0, s.Levels)
	if err != nil {
		return FocusResult{}, err
	}
	eta, err := EtaN(s.Levels)
	if err != nil {
		return FocusResult{}, err
	}
	res := FocusResult{EtaAnalytic: eta, PowerIn: lp.Power()}
	res.Wf = s.Lambda * s.F / (math.Pi * s.W0)

	// (2) Focal plane.
	u, err := Propagate(lp.T, s.M, s.Dx, s.Lambda, s.F)
	if err != nil {
		return FocusResult{}, err
	}
	rWin := s.WinW * res.Wf
	var inWin, peak float64
	peakR := 0
	peakC := 0
	for r := 0; r < s.M; r++ {
		y := Coord(r, s.M, s.Dx)
		for c := 0; c < s.M; c++ {
			x := Coord(c, s.M, s.Dx)
			v := u[r*s.M+c]
			I := real(v)*real(v) + imag(v)*imag(v)
			if math.Hypot(x, y) <= rWin {
				inWin += I
			}
			if I > peak {
				peak = I
				peakR, peakC = r, c
			}
		}
	}
	res.Eta = inWin * s.Dx * s.Dx / res.PowerIn
	dr := float64(signedIndex(peakR, s.M))
	dc := float64(signedIndex(peakC, s.M))
	res.PeakOffsetPx = math.Hypot(dr, dc)
	res.WindowFrac = 1 - math.Exp(-2*s.WinW*s.WinW)

	// (3) Axial peak: the maximum of the on-axis intensity on ±½ Rayleigh range
	// (unimodal there), by golden-section search to half a pixel, then a
	// parabolic refinement over three points spaced by half a pixel.
	spec, err := Spectrum(lp.T, s.M)
	if err != nil {
		return FocusResult{}, err
	}
	zR := math.Pi * res.Wf * res.Wf / s.Lambda
	zAx, err := axialPeak(spec, s.M, s.Dx, s.Lambda, s.F-0.5*zR, s.F+0.5*zR, 0.5*s.Dx)
	if err != nil {
		return FocusResult{}, err
	}
	res.AxialStep = s.Dx
	res.AxialZ = zAx
	return res, nil
}

func sqr(v complex128) float64 { return real(v)*real(v) + imag(v)*imag(v) }

// axialPeak returns the axial position in [a, b] of the maximum of the on-axis
// intensity |u(0,0,z)|², by golden-section search to bracket width tol and a
// parabolic refinement with points spaced tol apart. The intensity is assumed
// unimodal on [a, b]; the caller chooses the interval.
func axialPeak(spec []complex128, m int, dx, lambda, a, b, tol float64) (float64, error) {
	intensity := func(z float64) (float64, error) {
		v, err := AxialField(spec, m, dx, lambda, []float64{z})
		if err != nil {
			return 0, err
		}
		return sqr(v[0]), nil
	}
	const gr = 0.6180339887498949 // (√5 − 1)/2
	c := b - gr*(b-a)
	d := a + gr*(b-a)
	fc, err := intensity(c)
	if err != nil {
		return 0, err
	}
	fd, err := intensity(d)
	if err != nil {
		return 0, err
	}
	for b-a > tol {
		if fc > fd {
			b, d, fd = d, c, fc
			c = b - gr*(b-a)
			if fc, err = intensity(c); err != nil {
				return 0, err
			}
		} else {
			a, c, fc = c, d, fd
			d = a + gr*(b-a)
			if fd, err = intensity(d); err != nil {
				return 0, err
			}
		}
	}
	z := 0.5 * (a + b)
	y0, err := intensity(z - tol)
	if err != nil {
		return 0, err
	}
	y1, err := intensity(z)
	if err != nil {
		return 0, err
	}
	y2, err := intensity(z + tol)
	if err != nil {
		return 0, err
	}
	if den := y0 - 2*y1 + y2; den < 0 {
		z += 0.5 * (y0 - y2) / den * tol
	}
	return z, nil
}

// argmaxSq returns the index of the largest |v|².
func argmaxSq(v []complex128) int {
	best := 0
	for i := range v {
		if sqr(v[i]) > sqr(v[best]) {
			best = i
		}
	}
	return best
}
