package wave

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// SquareGrid returns the n×n uniform grid of normalised pupil points over
// [−1, 1]² with PX (fast) and PY (slow) running from −1 to +1 in n steps.
// It is the sampling PSFFromGrid expects: build a geo.OPDGrid with these
// points (geo.System.OPDGrid traces them), then hand the grid over.
func SquareGrid(n int) []geo.PupilPoint {
	if n < 2 {
		return []geo.PupilPoint{{}}
	}
	out := make([]geo.PupilPoint, 0, n*n)
	for j := 0; j < n; j++ {
		py := 2*float64(j)/float64(n-1) - 1
		for i := 0; i < n; i++ {
			px := 2*float64(i)/float64(n-1) - 1
			out = append(out, geo.PupilPoint{PX: px, PY: py})
		}
	}
	return out
}

// PSF is a diffraction point-spread function computed from a sampled
// wavefront, normalised so the aberrated peak equals the Strehl ratio (the
// peak of the same aperture without aberration is 1).
type PSF struct {
	N             int       // grid side (power of two)
	PixelMM       float64   // image-plane pixel spacing (mm)
	Intensity     []float64 // N×N, centred (row-major; the peak sits near [N/2, N/2])
	Strehl        float64   // peak / aberration-free peak of the same aperture
	WavelengthMM  float64   // wavelength used (mm)
	PupilRadiusMM float64   // pupil semi-diameter (mm)
	FocalMM       float64   // focal length mapping (mm)
	FNumber       float64   // Focal / (2·PupilRadius)
}

// Peak returns the maximum intensity and its grid position (row, column).
func (p *PSF) Peak() (float64, int, int) {
	best, br, bc := 0.0, 0, 0
	for r := 0; r < p.N; r++ {
		for c := 0; c < p.N; c++ {
			if v := p.Intensity[r*p.N+c]; v > best {
				best, br, bc = v, r, c
			}
		}
	}
	return best, br, bc
}

// PSFFromGrid computes the Fraunhofer PSF of the wavefront described by an
// OPD grid. g must be sampled on the uniform np×np square grid of
// SquareGrid(np) (as produced by geo.System.OPDGrid with those pupil
// points). The pupil samples are zero-padded into an FFT array of
// np·pad × np·pad, which oversamples the image plane by the same factor
// pad (pad = 1 gives one image sample per λ·F#; the Airy core of a typical
// system then spans only ~1.2 samples, so pad ≥ 4 is recommended for
// profile work). focalMM is the image-side effective focal length mapping
// pupil frequencies to image coordinates. The pupil amplitude is 1 inside
// the unit disc and 0 outside; samples whose ray did not reach the image
// plane (NaN OPD) are blocked as well.
func PSFFromGrid(g *geo.OPDGrid, focalMM float64, np, pad int) (*PSF, error) {
	if g == nil {
		return nil, fmt.Errorf("wave: nil OPD grid")
	}
	if !isPowerOfTwo(np) || np < 8 {
		return nil, fmt.Errorf("wave: pupil sampling %d must be a power of two >= 8", np)
	}
	if pad < 1 || pad > 64 {
		return nil, fmt.Errorf("wave: padding factor %d must be in [1, 64]", pad)
	}
	if len(g.Pupils) != np*np {
		return nil, fmt.Errorf("wave: OPD grid has %d samples, want %d×%d (sample it with wave.SquareGrid(%d))",
			len(g.Pupils), np, np, np)
	}
	if !(g.Pupil.Radius > 0) {
		return nil, fmt.Errorf("wave: pupil radius %g is not positive", g.Pupil.Radius)
	}
	lam := g.Lambda * 1e-6 // nm -> mm
	if !(lam > 0) {
		return nil, fmt.Errorf("wave: wavelength %g nm is not positive", g.Lambda)
	}
	if !(focalMM > 0) {
		return nil, fmt.Errorf("wave: focal length %g is not positive", focalMM)
	}

	nf := np * pad
	off := (nf - np) / 2
	pupil := make([]complex128, nf*nf)
	countA := 0
	var sumP complex128
	for j := 0; j < np; j++ {
		for i := 0; i < np; i++ {
			p := g.Pupils[j*np+i]
			if p.PX*p.PX+p.PY*p.PY > 1+1e-12 {
				continue // outside the aperture: blocked
			}
			opd := g.OPD[j*np+i]
			if math.IsNaN(opd) {
				continue // ray did not reach the image plane: blocked
			}
			ph := 2 * math.Pi * opd / lam
			c := complex(math.Cos(ph), math.Sin(ph))
			pupil[(j+off)*nf+(i+off)] = c
			countA++
			sumP += c
		}
	}
	if countA == 0 {
		return nil, fmt.Errorf("wave: the pupil is empty (every sample blocked or outside the aperture)")
	}

	FFT2D(pupil, nf, false)
	psf := make([]float64, nf*nf)
	peakAbs := float64(countA) * float64(countA) // |Σ A|², the aberration-free peak
	half := nf / 2
	for r := 0; r < nf; r++ {
		for c := 0; c < nf; c++ {
			v := pupil[r*nf+c]
			i2 := real(v)*real(v) + imag(v)*imag(v)
			psf[((r+half)%nf)*nf+((c+half)%nf)] = i2 / peakAbs
		}
	}
	strehl := (real(sumP)*real(sumP) + imag(sumP)*imag(sumP)) / peakAbs

	// Image-plane pixel spacing: pupil samples sit du = 2R/np apart, the
	// zero-padded FFT has bins Δfx = 1/(nf·du) and the Fraunhofer image
	// coordinate is x = λ·f·fx, so dx = λ·f·np/(2R·nf) = (np/nf)·λ·F#.
	dx := lam * focalMM * float64(np) / (2 * g.Pupil.Radius * float64(nf))
	return &PSF{
		N: nf, PixelMM: dx, Intensity: psf, Strehl: strehl,
		WavelengthMM: lam, PupilRadiusMM: g.Pupil.Radius, FocalMM: focalMM,
		FNumber: focalMM / (2 * g.Pupil.Radius),
	}, nil
}

// CutoffCyclesPerMM returns the diffraction cutoff frequency 1/(λ·F#) of the
// aperture (cycles/mm).
func (p *PSF) CutoffCyclesPerMM() float64 {
	return 1 / (p.WavelengthMM * p.FNumber)
}

// MTFAt returns the modulation transfer function at the spatial frequency
// freq (cycles/mm), evaluated along the x axis of the PSF grid as the
// normalised modulus of the 2-D optical transfer function at (freq, 0):
//
//	MTF(f) = |Σ psf(x,y)·exp(−2πi f x)| / Σ psf .
//
// The evaluation is direct (no interpolation), deterministic and exact for
// the sampled PSF; frequencies above the sampling Nyquist limit 1/(2·dx)
// are clamped to the Nyquist value.
func (p *PSF) MTFAt(freq float64) float64 {
	if freq < 0 {
		freq = -freq
	}
	nyq := 1 / (2 * p.PixelMM)
	if freq > nyq {
		freq = nyq
	}
	half := p.N / 2
	var re, im, total float64
	w := 2 * math.Pi * freq * p.PixelMM
	for r := 0; r < p.N; r++ {
		for c := 0; c < p.N; c++ {
			v := p.Intensity[r*p.N+c]
			if v == 0 {
				continue
			}
			x := float64(c - half) // pixel offset from the peak centre
			ph := w * x
			re += v * math.Cos(ph)
			im -= v * math.Sin(ph)
			total += v
		}
	}
	if total == 0 {
		return 0
	}
	return math.Hypot(re, im) / total
}

// MTFCurve returns the MTF at every frequency of freqs.
func (p *PSF) MTFCurve(freqs []float64) []float64 {
	out := make([]float64, len(freqs))
	for i, f := range freqs {
		out[i] = p.MTFAt(f)
	}
	return out
}

// RadialProfile returns the azimuthally averaged intensity as a function of
// the image-plane radius: the mean intensity in the annulus
// [(k−.5)·dr, (k+.5)·dr) for k = 0..bins−1, with dr = maxRadiusMM/bins.
// It is the standard way to compare a computed PSF with the analytic Airy
// profile.
func (p *PSF) RadialProfile(maxRadiusMM float64, bins int) []float64 {
	if bins <= 0 || !(maxRadiusMM > 0) {
		return nil
	}
	sum := make([]float64, bins)
	cnt := make([]int, bins)
	dr := maxRadiusMM / float64(bins)
	half := p.N / 2
	for r := 0; r < p.N; r++ {
		for c := 0; c < p.N; c++ {
			x := float64(c-half) * p.PixelMM
			y := float64(r-half) * p.PixelMM
			rad := math.Hypot(x, y)
			b := int(rad / dr)
			if b >= bins {
				continue
			}
			sum[b] += p.Intensity[r*p.N+c]
			cnt[b]++
		}
	}
	out := make([]float64, bins)
	for b := range out {
		if cnt[b] > 0 {
			out[b] = sum[b] / float64(cnt[b])
		}
	}
	return out
}

// AiryRadiusMM returns the first dark-ring radius of the aberration-free
// circular pupil of this aperture: 1.2197·λ·F# (mm). It is the analytic
// reference the PSF radial profile is checked against.
func (p *PSF) AiryRadiusMM() float64 {
	return 1.21967 * p.WavelengthMM * p.FNumber
}

// SPHERICAL/analytic reference -------------------------------------------------

// AiryMTF returns the analytic MTF of a perfect circular pupil at the
// relative frequency x = freq/cutoff (0 <= x <= 1):
//
//	MTF(x) = (2/π)·(acos(x) − x·sqrt(1−x²)) .
func AiryMTF(x float64) float64 {
	if x <= 0 {
		return 1
	}
	if x >= 1 {
		return 0
	}
	return (2 / math.Pi) * (math.Acos(x) - x*math.Sqrt(1-x*x))
}
