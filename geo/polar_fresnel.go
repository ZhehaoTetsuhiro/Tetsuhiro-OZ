package geo

import (
	"fmt"
	"math"
	"math/cmplx"
)

// ---------------------------------------------------------------------------
// Polarization of rays: Fresnel amplitude coefficients of one surface.
//
// Conventions shared by the polar_*.go files (they match optics/polar.go):
//
//   - Time dependence exp(-i·ω·t), so a plane wave is exp(i(k·r - ωt)); a
//     complex index is n = n' + i·κ with κ ≥ 0 for absorption (as in package
//     coating). Fresnel coefficients are real for lossless media below the
//     critical angle.
//   - A Jones vector (Ex, Ey) lives in the transverse frame of the ray, see
//     TransverseFrame: ex × ey = d, the propagation direction. On axis this is
//     the global (x, y) frame of optics/polar.go.
//   - s is the unit vector normal to the plane of incidence, s = d×n/|d×n|,
//     recomputed at every surface. The p vectors are p_in = s×d_in for the
//     incident wave, p_out = s×d_out for a transmitted wave and p_out = d_out×s
//     for a reflected wave (the unified convention below).
//   - The amplitude coefficients are electric-field amplitude ratios in the
//     unified convention, in which s and p coincide at normal incidence:
//
//     r_s = (n1·cosθ1 − n2·cosθ2)/(n1·cosθ1 + n2·cosθ2)
//     t_s = 2·n1·cosθ1/(n1·cosθ1 + n2·cosθ2)
//     r_p = (n1·cosθ2 − n2·cosθ1)/(n1·cosθ2 + n2·cosθ1)   = −r_p^BW
//     t_p = 2·n1·cosθ1/(n1·cosθ2 + n2·cosθ1)              =  t_p^BW
//
//     where r_p^BW, t_p^BW are the Born & Wolf coefficients (p vectors s×d for
//     every wave). Only the sign of r_p differs. Check: a perfectly conducting
//     mirror has r_s = r_p = −1 in this convention, see MetalMirrorAmplitudes.
// ---------------------------------------------------------------------------

// Amplitudes are the complex amplitude coefficients of one surface for the s
// and p eigenstates. R* are reflected and T* transmitted electric-field
// amplitude ratios in the convention above. For a perfect mirror only the R*
// are used.
type Amplitudes struct {
	RS, TS complex128
	RP, TP complex128
}

// AmplitudeFunc returns the amplitude coefficients of the interaction described
// by st (one traversed surface of a trace: incident and outgoing directions,
// the +z oriented normal, the indices before and after the surface) at the
// vacuum wavelength wavelengthNm. It is the hook through which coated surfaces
// are supplied (see CoatedAmplitudes); FresnelStep is the default.
type AmplitudeFunc func(st Step, wavelengthNm float64) (Amplitudes, error)

// FresnelStep is the default AmplitudeFunc. Refracting surfaces use the plane
// interface Fresnel coefficients at the actual incidence angle of the step, and
// Mirror surfaces use MetalMirrorAmplitudes.
func FresnelStep(st Step, wavelengthNm float64) (Amplitudes, error) {
	if st.Mirror {
		return MetalMirrorAmplitudes(), nil
	}
	cosI := math.Abs(st.DirIn.Dot(st.Normal))
	a, _, err := FresnelAmplitudes(complex(st.NBefore, 0), complex(st.NAfter, 0), cosI)
	return a, err
}

// MetalMirrorAmplitudes returns the coefficients of an ideal conductor. The
// tangential electric field vanishes at the surface, which gives r_s = r_p = −1
// for every angle (the limit n2 → ∞ of FresnelAmplitudes in the unified
// convention) and no transmitted wave.
func MetalMirrorAmplitudes() Amplitudes {
	return Amplitudes{RS: -1, RP: -1}
}

// FresnelAmplitudes returns the amplitude coefficients of a plane interface from
// index n1 (incidence side) into index n2 for an incidence cosine
// cosI = |d·n| in (0, 1], together with the transmitted cosine cosT.
//
// The transmitted cosine follows Snell's law on the physical branch:
// cosT = q/n2 with q = sqrt(n2² − (n1·sinθ1)²) and Im q ≥ 0. For real indices
// below the critical angle everything is real. Beyond it cosT is imaginary,
// |r_s| = |r_p| = 1 and the transmitted wave is evanescent (a ray cannot follow
// it; TracePolarized stops at the TIR status of the trace).
//
// With the tilted admittances η_s = n·cosθ and η_p = n/cosθ the coefficients
// are r = (η1−η2)/(η1+η2) and t_s = 2η1s/(η1s+η2s), while t_p is the
// tangential ratio 2η1p/(η1p+η2p) multiplied by cosθ1/cosθ2 to give the
// electric-field ratio.
func FresnelAmplitudes(n1, n2 complex128, cosI float64) (Amplitudes, complex128, error) {
	var a Amplitudes
	if !(cosI > 0 && cosI <= 1) {
		return a, 0, fmt.Errorf("geo: polarization: incidence cosine %g outside (0,1]", cosI)
	}
	if n1 == 0 || n2 == 0 || cmplx.IsNaN(n1) || cmplx.IsNaN(n2) || cmplx.IsInf(n1) || cmplx.IsInf(n2) {
		return a, 0, fmt.Errorf("geo: polarization: refractive index must be finite and non-zero")
	}
	sinI := math.Sqrt(1 - cosI*cosI)
	alpha := n1 * complex(sinI, 0)
	q := polSqrt(n2*n2 - alpha*alpha)
	cosT := q / n2

	// s-polarization: η_s = n·cosθ.
	es1 := n1 * complex(cosI, 0)
	es2 := q
	a.RS = (es1 - es2) / (es1 + es2)
	a.TS = 2 * es1 / (es1 + es2)

	// p-polarization: η_p = n/cosθ = n²/η_s.
	ep1 := n1 / complex(cosI, 0)
	ep2 := n2 / cosT
	a.RP = (ep1 - ep2) / (ep1 + ep2)
	a.TP = 2 * ep1 / (ep1 + ep2) * complex(cosI, 0) / cosT
	return a, cosT, nil
}

// polSqrt is the square root on the physical branch of the e^{-iωt} convention:
// Im ≥ 0, so that a forward wave exp(i·k_z·z) decays into an absorbing or
// evanescent medium. For a positive real argument it is the ordinary root.
func polSqrt(q complex128) complex128 {
	r := cmplx.Sqrt(q)
	if imag(r) < 0 {
		r = -r
	}
	return r
}
