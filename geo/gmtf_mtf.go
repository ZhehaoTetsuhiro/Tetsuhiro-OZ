package geo

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// Geometric MTF from a ray distribution (G11)
// ---------------------------------------------------------------------------
//
// The geometric MTF is obtained from the ray footprint on a plane: the rays
// are binned into a line-spread function (LSF) along one transverse axis, and
// the modulus of its Fourier transform gives the MTF,
//
//	MTF(f) = | Σ_k p_k · exp(−2πi f x_k) | / Σ_k p_k ,
//
// where x_k are the bin centres and p_k the fraction of rays in bin k. The
// sum is evaluated directly at the requested frequencies (no FFT is needed
// for a few hundred frequencies, and the frequencies are not restricted to a
// grid). The two-dimensional PSF gives the same radial MTF through the
// projection-slice theorem, so the LSF is sufficient for the one-dimensional
// MTF of an axis.

// LSF is the line-spread function of a ray distribution: the fraction of rays
// per equal-width bin along one transverse axis.
type LSF struct {
	Lo    float64   // left edge of the first bin (mm)
	Width float64   // bin width Δ (mm)
	Prob  []float64 // fraction of the rays in each bin (sums to 1)
	N     int       // number of rays binned
}

// BinLSF bins the transverse coordinates x into nBins equal bins covering
// [lo, hi). The range must contain every coordinate: a ray outside it is an
// error rather than a silently dropped sample, because a dropped tail changes
// the MTF at every frequency.
func BinLSF(x []float64, lo, hi float64, nBins int) (LSF, error) {
	if len(x) == 0 {
		return LSF{}, fmt.Errorf("geo: LSF: no rays")
	}
	if nBins <= 0 || !(hi > lo) || math.IsInf(lo, 0) || math.IsInf(hi, 0) {
		return LSF{}, fmt.Errorf("geo: LSF: invalid bin range [%g, %g) with %d bins", lo, hi, nBins)
	}
	width := (hi - lo) / float64(nBins)
	count := make([]float64, nBins)
	for i, v := range x {
		if !(v >= lo && v < hi) {
			return LSF{}, fmt.Errorf("geo: LSF: ray %d at %g lies outside [%g, %g)", i, v, lo, hi)
		}
		k := int((v - lo) / width)
		if k >= nBins { // v just below hi, rounding up
			k = nBins - 1
		}
		count[k]++
	}
	prob := make([]float64, nBins)
	n := float64(len(x))
	for k, c := range count {
		prob[k] = c / n
	}
	return LSF{Lo: lo, Width: width, Prob: prob, N: len(x)}, nil
}

// Center returns the centre of bin k (mm).
func (l LSF) Center(k int) float64 { return l.Lo + (float64(k)+0.5)*l.Width }

// MTF returns the geometric MTF at the spatial frequencies freqs (cycles/mm),
// with the modulus of the bin-centre sum of the LSF (see the file comment). The
// value is exactly 1 at f = 0 by construction.
//
// Binning is a box average, whose Fourier transform is the box window
// sinc(fΔ) multiplying the bin-centre sum. The bin-centre sum is returned
// without that factor, so the estimate is the MTF of the ray distribution; the
// neglected difference between the bin-centre sum and the true transform of
// the ray density is of relative order (πfΔ)²/6. BoxWindow gives the factor
// for users who want the detector-integrated MTF instead.
func (l LSF) MTF(freqs []float64) []float64 {
	out := make([]float64, len(freqs))
	for i, f := range freqs {
		var re, im float64
		for k, p := range l.Prob {
			if p == 0 {
				continue
			}
			s, c := math.Sincos(-2 * math.Pi * f * l.Center(k))
			re += p * c
			im += p * s
		}
		out[i] = math.Hypot(re, im)
	}
	return out
}

// BoxWindow returns sinc(fΔ) = sin(πfΔ)/(πfΔ), the Fourier transform of a box
// of width Δ normalised to unit area; sinc(0) = 1.
func (l LSF) BoxWindow(f float64) float64 {
	u := math.Pi * f * l.Width
	if u == 0 {
		return 1
	}
	return math.Sin(u) / u
}

// DiscOTF returns the optical transfer function of a uniformly filled disc of
// radius R (mm): 2·J1(x)/x with x = 2π·f·R (math.J1). It is real and changes
// sign at its zeros x = j₁,ₖ (first zero x ≈ 3.8317).
func DiscOTF(f, R float64) float64 {
	x := 2 * math.Pi * f * R
	if x == 0 {
		return 1
	}
	return 2 * math.J1(x) / x
}

// DiscMTF returns the geometric MTF of a uniformly filled disc of radius R
// (mm), the modulus |2·J1(x)/x| with x = 2π·f·R. This is the MTF of a uniform
// circular ray bundle whose footprint is the disc, e.g. an aberration-free
// bundle observed a distance δ from its focus, R = a·|δ|/L (a: pupil radius,
// L: pupil-to-focus distance).
func DiscMTF(f, R float64) float64 { return math.Abs(DiscOTF(f, R)) }

// EqualAreaDisc returns n points of the unit disc with equal area per point
// (Vogel's sunflower spiral): point k sits at radius sqrt((k+½)/n) and azimuth
// k·(π(3−√5)). Each point represents the same area π/n, so the set is an
// equal-area, deterministic sampling of the disc with no random noise. Its
// quality is measured, not assumed: the MTF test reports the error it produces.
func EqualAreaDisc(n int) []PupilPoint {
	if n <= 0 {
		return nil
	}
	golden := math.Pi * (3 - math.Sqrt(5))
	out := make([]PupilPoint, n)
	for k := range out {
		r := math.Sqrt((float64(k) + 0.5) / float64(n))
		a := golden * float64(k)
		s, c := math.Sincos(a)
		out[k] = PupilPoint{PX: r * c, PY: r * s}
	}
	return out
}

// ConvergingBundleHits returns where the rays of an aberration-free focusing
// bundle land on the plane z = imageZ. The bundle has its pupil on the plane
// z = focusZ − L, radius a, filled according to pts (unit-disc coordinates),
// and every ray is aimed at the focal point (0, 0, focusZ): this is what a
// perfect lens produces, so the rays have no aberration by construction. The
// rays are traced by System.Trace through the system given (its surfaces must
// lie between the pupil plane and the image plane; for a plane detector this
// is a single geo.Plane surface at imageZ), and the returned points are the
// final ray origins on the image plane.
//
// For a focal point at the origin the footprint of a defocused bundle is the
// disc of radius a·|imageZ − focusZ|/L, exactly uniform, which is the geometry
// behind DiscMTF.
func (s *System) ConvergingBundleHits(focusZ, L, a float64, pts []PupilPoint, lambdaNm float64) ([]Vec3, error) {
	if !(L > 0) || !(a > 0) {
		return nil, fmt.Errorf("geo: converging bundle: L=%g and a=%g must be positive", L, a)
	}
	pupilZ := focusZ - L
	focus := V(0, 0, focusZ)
	hits := make([]Vec3, len(pts))
	for i, p := range pts {
		origin := V(a*p.PX, a*p.PY, pupilZ)
		r := RayThrough(origin, focus, lambdaNm)
		if !r.Valid() {
			return nil, fmt.Errorf("geo: converging bundle: degenerate ray %d", i)
		}
		res := s.Trace(r)
		if !res.OK() {
			return nil, fmt.Errorf("geo: converging bundle: ray %d: %v", i, res.Error())
		}
		hits[i] = res.Ray.Origin
	}
	return hits, nil
}
