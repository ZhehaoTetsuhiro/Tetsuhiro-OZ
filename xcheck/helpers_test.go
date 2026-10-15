package xcheck

import (
	"math"
	"math/cmplx"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/optics"
)

// Units. The geometric layer (geo, lens) works in mm; the wave kernel (optics)
// works in SI (m). Every conversion is written out at its use.

// dLineMM and dLineM are the d-line wavelength in mm and m.
const (
	dLineMM = geo.Wavelengthd * 1e-6
	dLineM  = geo.Wavelengthd * 1e-9
)

// achromatDesign returns the f = 100 mm, F/5 cemented doublet of the lens
// package (N-BK7 + F2, auto-designed and DLS-polished), its paraxial data at
// the d-line and its entrance pupil.
func achromatDesign(t *testing.T) (*geo.System, *geo.Paraxial, geo.Pupil) {
	t.Helper()
	p, _, err := lens.DesignAchromat(100, "N-BK7", "F2", 10, nil)
	if err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	pup, err := sys.EntrancePupil(geo.Field{}, geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	return sys, par, pup
}

// rayOPD traces the pupil points to the image plane z (mm) and returns the
// OPD grid (on-axis field, d-line). Blocked rays carry NaN OPD.
func rayOPD(t *testing.T, sys *geo.System, pup geo.Pupil, pts []geo.PupilPoint, z float64) *geo.OPDGrid {
	t.Helper()
	g, err := sys.OPDGrid(geo.Field{}, pup, pts, z, geo.Wavelengthd, geo.OPDOptions{})
	if err != nil {
		t.Fatal(err)
	}
	return g
}

// wavesOf converts an OPD in mm to waves at the d-line (blocked stays NaN).
func wavesOf(opdMM []float64) []float64 {
	out := make([]float64, len(opdMM))
	for i, v := range opdMM {
		out[i] = v / dLineMM
	}
	return out
}

// opdToM converts an OPD in mm to m (blocked stays NaN).
func opdToM(opdMM []float64) []float64 {
	out := make([]float64, len(opdMM))
	for i, v := range opdMM {
		out[i] = v * 1e-3
	}
	return out
}

// goldenMin returns the minimiser of a unimodal f on [a, b] to tolerance tol.
func goldenMin(f func(float64) float64, a, b, tol float64) float64 {
	gr := (math.Sqrt(5) - 1) / 2
	c, d := b-gr*(b-a), a+gr*(b-a)
	fc, fd := f(c), f(d)
	for b-a > tol {
		if fc < fd {
			b, d, fd = d, c, fc
			c = b - gr*(b-a)
			fc = f(c)
		} else {
			a, c, fc = c, d, fd
			d = a + gr*(b-a)
			fd = f(d)
		}
	}
	return (a + b) / 2
}

// bestFocusZ returns the image-plane position z (mm) that minimises the RMS
// wavefront error of the on-axis field, and that RMS in waves. The search is a
// golden-section minimisation of the RMS over [z0−0.2, z0+0.2] mm, evaluated
// with a 32×64 product rule on the pupil.
func bestFocusZ(t *testing.T, sys *geo.System, pup geo.Pupil, z0 float64) (float64, float64) {
	t.Helper()
	qp, qw := diskQuadrature(32, 64)
	rms := func(z float64) float64 {
		_, r := integralStrehl(wavesOf(rayOPD(t, sys, pup, qp, z).OPD), qw)
		return r
	}
	z := goldenMin(rms, z0-0.2, z0+0.2, 1e-7)
	return z, rms(z)
}

// legendre returns P_n(z) and P_n'(z) by the three-term recurrence.
func legendre(n int, z float64) (p, dp float64) {
	if n == 0 {
		return 1, 0
	}
	p0, p1 := 1.0, z
	for k := 2; k <= n; k++ {
		p0, p1 = p1, ((2*float64(k)-1)*z*p1-(float64(k)-1)*p0)/float64(k)
	}
	return p1, float64(n) * (z*p1 - p0) / (z*z - 1)
}

// gaussLegendre returns the n-point Gauss–Legendre nodes and weights on [−1, 1]
// (Newton iteration on P_n from the Chebyshev initial guess).
func gaussLegendre(n int) (x, w []float64) {
	x = make([]float64, n)
	w = make([]float64, n)
	for i := 0; i < n; i++ {
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		for it := 0; it < 100; it++ {
			p, dp := legendre(n, z)
			dz := p / dp
			z -= dz
			if math.Abs(dz) < 1e-16 {
				break
			}
		}
		_, dp := legendre(n, z)
		x[i] = z
		w[i] = 2 / ((1 - z*z) * dp * dp)
	}
	return x, w
}

// diskQuadrature returns a product rule on the unit disc: M Gauss–Legendre
// nodes in ρ on (0, 1) times K equispaced azimuths. The weights sum to π, and
// ∫_disc f dA ≈ Σ_i w_i f(p_i) is exact for polynomials in ρ of degree ≤ 2M−1
// and for trigonometric degree < K in θ.
func diskQuadrature(M, K int) ([]geo.PupilPoint, []float64) {
	x, w := gaussLegendre(M)
	pts := make([]geo.PupilPoint, 0, M*K)
	wts := make([]float64, 0, M*K)
	for m := 0; m < M; m++ {
		rho := (x[m] + 1) / 2
		wr := w[m] / 2 * rho
		for k := 0; k < K; k++ {
			th := 2 * math.Pi * float64(k) / float64(K)
			pts = append(pts, geo.PupilPoint{PX: rho * math.Cos(th), PY: rho * math.Sin(th)})
			wts = append(wts, wr*2*math.Pi/float64(K))
		}
	}
	return pts, wts
}

// integralStrehl returns the diffraction-integral Strehl ratio
// |∫ exp(i·2π·W) dA|² / (∫ dA)² of a wavefront W (waves) sampled with weights
// w, and the weighted RMS of W about its mean (waves). Blocked samples (NaN)
// must be removed by the caller.
func integralStrehl(wavesAt []float64, w []float64) (strehl, rms float64) {
	var sw, mean float64
	for i, v := range wavesAt {
		sw += w[i]
		mean += w[i] * v
	}
	mean /= sw
	var sumE complex128
	var var2 float64
	for i, v := range wavesAt {
		sumE += complex(w[i], 0) * cmplx.Exp(complex(0, 2*math.Pi*v))
		d := v - mean
		var2 += w[i] * d * d
	}
	strehl = (real(sumE)*real(sumE) + imag(sumE)*imag(sumE)) / (sw * sw)
	return strehl, math.Sqrt(var2 / sw)
}

// marechal returns the Maréchal approximation exp(−(2πσ)²) of the Strehl
// ratio for an RMS wavefront error σ in waves.
func marechal(rmsWaves float64) float64 {
	return math.Exp(-(2 * math.Pi * rmsWaves) * (2 * math.Pi * rmsWaves))
}

// relErr is |got − want| / |want|.
func relErr(got, want float64) float64 { return math.Abs(got-want) / math.Abs(want) }

// lattice is a circular pupil on the centred N×N grid of the wave kernel: the
// pixel centres (i−N/2, j−N/2) with radius < R pixels. Pixel k of the lattice
// has normalised pupil coordinates (px, py) = offset/R, so the pupil radius
// is R pixels whatever the physical pixel size.
type lattice struct {
	N, R   int
	idx    []int     // row-major field index of each pupil pixel
	px, py []float64 // normalised pupil coordinates of each pixel
}

func newLattice(N, R int) lattice {
	l := lattice{N: N, R: R}
	c := N / 2
	for j := 0; j < N; j++ {
		for i := 0; i < N; i++ {
			dx, dy := float64(i-c), float64(j-c)
			if dx*dx+dy*dy < float64(R*R)-1e-9 {
				l.idx = append(l.idx, j*N+i)
				l.px = append(l.px, dx/float64(R))
				l.py = append(l.py, dy/float64(R))
			}
		}
	}
	return l
}

// points returns the normalised pupil points of the lattice for geo tracing.
func (l lattice) points() []geo.PupilPoint {
	pts := make([]geo.PupilPoint, len(l.px))
	for i := range l.px {
		pts[i] = geo.PupilPoint{PX: l.px[i], PY: l.py[i]}
	}
	return pts
}

// fraunhofer returns the focal-plane field of the lattice pupil. opdM holds the
// wavefront in metres at each lattice pixel (nil = unaberrated; NaN = blocked
// ray, which carries no amplitude). The field is propagated by the kernel's
// Fraunhofer method over zM metres: this is the focal-plane field of an ideal
// lens of focal length zM (the lens's quadratic phase cancels the Fresnel
// kernel's quadratic phase over that distance). See TestE1StrehlGeoVsWave for
// why the lens element with an ASM propagation is not used on this grid.
func (l lattice) fraunhofer(lamM, dxM, zM float64, opdM []float64) *optics.Field {
	f := optics.NewField(l.N, dxM, false)
	k := 2 * math.Pi / lamM
	for n, id := range l.idx {
		if opdM == nil {
			f.Ex[id] = 1
			continue
		}
		if math.IsNaN(opdM[n]) {
			continue
		}
		ph := k * opdM[n]
		f.Ex[id] = complex(math.Cos(ph), math.Sin(ph))
	}
	ctx := &optics.Context{Wavelength: lamM, Evanescent: "decay"}
	if err := optics.Propagate(f, zM, optics.MethodFraunhofer, ctx); err != nil {
		panic(err)
	}
	return f
}

// peak returns the maximum intensity of a field.
func peak(f *optics.Field) float64 {
	m := 0.0
	for i := range f.Ex {
		if v := f.Intensity(i); v > m {
			m = v
		}
	}
	return m
}
