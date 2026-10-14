package coating

import (
	"errors"
	"fmt"
	"math"
	"math/cmplx"
)

// Polarization selects the polarization state used to evaluate a stack.
type Polarization uint8

const (
	// PolS is s-polarization: E perpendicular to the plane of incidence, so the
	// tangential electric field is Ey and the tilted admittance is n·cos(theta).
	PolS Polarization = iota
	// PolP is p-polarization: E inside the plane of incidence, so the tangential
	// electric field is Ex and the tilted admittance is n/cos(theta).
	PolP
	// PolUnpolarized averages the s and p results equally; use it for natural
	// light or when the polarization of the beam is not defined.
	PolUnpolarized
)

// String implements fmt.Stringer.
func (p Polarization) String() string {
	switch p {
	case PolS:
		return "s"
	case PolP:
		return "p"
	case PolUnpolarized:
		return "unpolarized"
	default:
		return fmt.Sprintf("Polarization(%d)", uint8(p))
	}
}

// ParsePolarization converts a short label into a Polarization: "s", "p",
// "unpolarized"/"u"/"avg"/"average"/"" (empty means unpolarized). The
// comparison ignores case.
func ParsePolarization(label string) (Polarization, error) {
	switch lower(label) {
	case "s", "te", "perpendicular":
		return PolS, nil
	case "p", "tm", "parallel":
		return PolP, nil
	case "", "u", "avg", "average", "unpolarized", "unpolarised", "natural":
		return PolUnpolarized, nil
	}
	return PolUnpolarized, fmt.Errorf("coating: unknown polarization %q", label)
}

// lower is a tiny ASCII lower-caser, kept local so the package keeps its
// dependency list at the standard library minimum.
func lower(s string) string {
	b := []byte(s)
	for i := range b {
		if b[i] >= 'A' && b[i] <= 'Z' {
			b[i] += 'a' - 'A'
		}
	}
	return string(b)
}

// Layer is a single coating layer: a material and a physical thickness in nm.
type Layer struct {
	Material    Material
	ThicknessNm float64
}

// Stack is a complete thin-film system: the incident medium, the layer
// sequence (front surface first) and a semi-infinite substrate.
//
// The substrate is assumed semi-infinite and the whole stack is treated
// coherently. TODO: incoherent back-surface reflection of a finite plane
// parallel substrate is not implemented - see the package documentation.
type Stack struct {
	Incident  Material
	Layers    []Layer
	Substrate Material
}

// NewStack builds a stack from an incidence medium, a substrate and the layers.
func NewStack(incident, substrate Material, layers ...Layer) Stack {
	return Stack{Incident: incident, Substrate: substrate, Layers: layers}
}

// TotalThicknessNm returns the summed physical thickness of all layers.
func (s Stack) TotalThicknessNm() float64 {
	t := 0.0
	for _, l := range s.Layers {
		t += l.ThicknessNm
	}
	return t
}

// Clone returns a deep copy of the stack (the layer slice is copied).
func (s Stack) Clone() Stack {
	out := s
	out.Layers = append([]Layer(nil), s.Layers...)
	return out
}

// ThicknessesNm returns the layer thicknesses as a slice.
func (s Stack) ThicknessesNm() []float64 {
	out := make([]float64, len(s.Layers))
	for i, l := range s.Layers {
		out[i] = l.ThicknessNm
	}
	return out
}

// WithThicknesses returns a copy of the stack in which the layer thicknesses
// are replaced by t. It is the natural inverse of ThicknessesNm and is used by
// the optimiser and by callers that keep their own parameter vector.
func (s Stack) WithThicknesses(t []float64) Stack {
	out := s.Clone()
	for i := range out.Layers {
		if i < len(t) {
			out.Layers[i].ThicknessNm = t[i]
		}
	}
	return out
}

// Reversed returns a stack with the layer order flipped and the substrate and
// incidence media exchanged: the reflectance of the reversed coating seen from
// the "air" side equals the transmittance behaviour of the original one. It is
// useful when a design has to be deposited in the opposite order.
func (s Stack) Reversed() Stack {
	out := Stack{Incident: s.Substrate, Substrate: s.Incident}
	out.Layers = make([]Layer, len(s.Layers))
	for i, l := range s.Layers {
		out.Layers[len(s.Layers)-1-i] = l
	}
	return out
}

// Validate reports whether the stack can be evaluated: both surrounding media
// and every layer need a dispersion model, and every thickness must be finite
// and non-negative.
func (s Stack) Validate() error {
	if !s.Incident.Valid() {
		return errors.New("coating: stack has no valid incidence medium")
	}
	if !s.Substrate.Valid() {
		return errors.New("coating: stack has no valid substrate")
	}
	for i, l := range s.Layers {
		if !l.Material.Valid() {
			return fmt.Errorf("coating: layer %d has no valid material", i)
		}
		if math.IsNaN(l.ThicknessNm) || math.IsInf(l.ThicknessNm, 0) || l.ThicknessNm < 0 {
			return fmt.Errorf("coating: layer %d has invalid thickness %g nm", i, l.ThicknessNm)
		}
	}
	return nil
}

// Sample is one evaluated operating point of a stack: a wavelength, an
// incidence angle and a polarization, together with the resulting power ratios
// and complex amplitudes.
type Sample struct {
	LambdaNm float64
	AngleDeg float64
	Pol      Polarization

	// R, T, A are the power reflectance, transmittance and absorptance of the
	// semi-infinite substrate system; R + T + A = 1.
	R float64 // reflected fraction
	T float64 // transmitted fraction
	A float64 // absorbed fraction, 1 - R - T

	// Rc and Tc are the complex amplitude reflection and transmission
	// coefficients with |Rc|^2 = R (per polarization). The transmission
	// coefficient is the tangential (in-plane) field amplitude at the substrate
	// interface for an incident amplitude of 1; for PolUnpolarized both are
	// zeroed because the s and p phases are not averaged.
	Rc, Tc complex128
}

// PhaseRdeg returns the phase of the reflected amplitude in degrees, in
// (-180, 180].
func (s Sample) PhaseRdeg() float64 { return phaseDeg(s.Rc) }

// PhaseTdeg returns the phase of the transmitted amplitude in degrees, in
// (-180, 180].
func (s Sample) PhaseTdeg() float64 { return phaseDeg(s.Tc) }

// phaseDeg converts a complex amplitude into its phase in degrees.
func phaseDeg(c complex128) float64 {
	if c == 0 {
		return 0
	}
	return math.Atan2(imag(c), real(c)) * 180 / math.Pi
}

// ---------------------------------------------------------------------------
// Core evaluation
// ---------------------------------------------------------------------------

// principalSqrt returns the square root of q = n^2 - alpha^2 on the physical
// branch for the exp(-i·omega·t) convention used throughout this package: the
// imaginary part is non-negative, so the forward wave exp(i·k_z·z) decays when
// it enters an absorbing medium (n = n' + i·k, k >= 0) and is evanescent when
// the transverse wavevector exceeds the medium's index. For a real positive q
// the result is the ordinary positive square root.
func principalSqrt(q complex128) complex128 {
	r := cmplx.Sqrt(q)
	if imag(r) < 0 || (imag(r) == 0 && real(r) < 0) {
		r = -r
	}
	return r
}

// admittance returns the tilted optical admittances of a medium of complex
// index n for a plane wave with Snell invariant alpha = n0·sin(theta0):
//
//	etaS = n·cos(theta) = sqrt(n^2 - alpha^2),   etaP = n/cos(theta) = n^2/etaS.
//
// cosTheta = sqrt(n^2-alpha^2)/n is returned as well (used for phase
// thicknesses and field profiles). A degenerate medium (n·cos(theta) = 0)
// yields an error rather than a division by zero.
func admittance(n, alpha complex128) (etaS, etaP, cosTheta complex128, err error) {
	if math.IsNaN(real(n)) || math.IsNaN(imag(n)) || math.IsInf(real(n), 0) || math.IsInf(imag(n), 0) {
		return 0, 0, 0, errors.New("coating: non-finite refractive index")
	}
	q := n*n - alpha*alpha
	if q == 0 {
		return 0, 0, 0, errors.New("coating: degenerate propagation (n·cos(theta) = 0)")
	}
	c := principalSqrt(q)
	return c, n * n / c, c / n, nil
}

// polCoefficients computes the complex amplitudes and power ratios of a stack
// for one single polarization. It is the hot path of the package: no
// allocations, no maps, one pass over the layers.
func (s Stack) polCoefficients(lambdaNm, angleDeg float64, pol Polarization) (rc, tc complex128, R, T float64, err error) {
	if !(lambdaNm > 0) || math.IsInf(lambdaNm, 0) {
		err = fmt.Errorf("coating: wavelength %g nm must be finite and positive", lambdaNm)
		return
	}
	if math.IsNaN(angleDeg) || math.Abs(angleDeg) >= 90 {
		err = fmt.Errorf("coating: incidence angle %g deg must satisfy |theta| < 90", angleDeg)
		return
	}
	useP := pol == PolP
	k0 := 2 * math.Pi / lambdaNm
	theta0 := angleDeg * math.Pi / 180

	n0 := s.Incident.IndexNm(lambdaNm)
	ns := s.Substrate.IndexNm(lambdaNm)
	alpha := n0 * complex(math.Sin(theta0), 0)

	eta0S, eta0P, _, err := admittance(n0, alpha)
	if err != nil {
		return 0, 0, 0, 0, fmt.Errorf("coating: incidence medium %s at %g nm: %w", s.Incident, lambdaNm, err)
	}
	eta := eta0S
	if useP {
		eta = eta0P
	}
	if real(eta) == 0 {
		err = fmt.Errorf("coating: incidence medium %s has zero admittance at %g nm", s.Incident, lambdaNm)
		return
	}

	// Ordered product of the characteristic matrices, accumulated scalar-wise
	// as (m11, m12; m21, m22) starting from the identity.
	m11, m12 := complex(1, 0), complex(0, 0)
	m21, m22 := complex(0, 0), complex(1, 0)
	for i := range s.Layers {
		layer := &s.Layers[i]
		if layer.ThicknessNm < 0 || math.IsNaN(layer.ThicknessNm) {
			err = fmt.Errorf("coating: layer %d has invalid thickness %g nm", i, layer.ThicknessNm)
			return
		}
		if layer.ThicknessNm == 0 {
			continue
		}
		n := layer.Material.IndexNm(lambdaNm)
		etaS, etaP, _, aerr := admittance(n, alpha)
		if aerr != nil {
			return 0, 0, 0, 0, fmt.Errorf("coating: layer %d (%s) at %g nm: %w", i, layer.Material, lambdaNm, aerr)
		}
		e := etaS
		if useP {
			e = etaP
		}
		delta := complex(k0*layer.ThicknessNm, 0) * etaS
		cd, sd := cmplx.Cos(delta), cmplx.Sin(delta)
		mj12 := complex(0, -1) * sd / e
		mj21 := complex(0, -1) * e * sd
		// M *= M_j
		m11, m12, m21, m22 = m11*cd+m12*mj21, m11*mj12+m12*cd, m21*cd+m22*mj21, m21*mj12+m22*cd
	}

	etaSubS, etaSubP, _, err := admittance(ns, alpha)
	if err != nil {
		return 0, 0, 0, 0, fmt.Errorf("coating: substrate %s at %g nm: %w", s.Substrate, lambdaNm, err)
	}
	etaSub := etaSubS
	if useP {
		etaSub = etaSubP
	}

	den := m11 + m12*etaSub
	if den == 0 {
		err = fmt.Errorf("coating: singular characteristic matrix at %g nm", lambdaNm)
		return
	}
	Y := (m21 + m22*etaSub) / den
	d := eta + Y
	if d == 0 {
		err = fmt.Errorf("coating: singular admittance at %g nm", lambdaNm)
		return
	}
	rc = (eta - Y) / d
	tc = 2 * eta / (d * den)

	R = real(rc * cmplx.Conj(rc))
	switch {
	case real(etaSub) > 0:
		T = (real(etaSub) / real(eta)) * real(tc*cmplx.Conj(tc))
	case real(etaSub) == 0:
		T = 0 // evanescent substrate: no propagating flux
	default:
		// An evanescent (or gain) substrate would flip the sign of the flux;
		// passive, propagating media are the supported case.
		T = 0
	}
	// Trim rounding noise: |t|^2 may pick up a relative error of a few 1e-16.
	if T < 0 && T > -1e-12 {
		T = 0
	}
	return rc, tc, R, T, nil
}

// Evaluate computes one operating point of the stack. For PolUnpolarized the
// power ratios are the arithmetic mean of the s and p results and the complex
// amplitudes are zeroed (there is no single meaningful phase).
func (s Stack) Evaluate(lambdaNm, angleDeg float64, pol Polarization) (Sample, error) {
	out := Sample{LambdaNm: lambdaNm, AngleDeg: angleDeg, Pol: pol}
	switch pol {
	case PolS, PolP:
		rc, tc, R, T, err := s.polCoefficients(lambdaNm, angleDeg, pol)
		if err != nil {
			return Sample{}, err
		}
		out.Rc, out.Tc = rc, tc
		out.R, out.T = R, T
		out.A = absorptance(R, T)
		return out, nil
	case PolUnpolarized:
		_, _, Rs, Ts, err := s.polCoefficients(lambdaNm, angleDeg, PolS)
		if err != nil {
			return Sample{}, err
		}
		_, _, Rp, Tp, err := s.polCoefficients(lambdaNm, angleDeg, PolP)
		if err != nil {
			return Sample{}, err
		}
		out.R = 0.5 * (Rs + Rp)
		out.T = 0.5 * (Ts + Tp)
		out.A = absorptance(out.R, out.T)
		return out, nil
	default:
		return Sample{}, fmt.Errorf("coating: invalid polarization value %d", uint8(pol))
	}
}

// absorptance returns 1 - R - T, clamping pure rounding noise to zero so that
// a lossless stack reports A = 0 exactly.
func absorptance(R, T float64) float64 {
	a := 1 - R - T
	if a < 0 && a > -1e-12 {
		return 0
	}
	return a
}

// Reflectance is a convenience wrapper returning only the reflected fraction
// of Evaluate.
func (s Stack) Reflectance(lambdaNm, angleDeg float64, pol Polarization) (float64, error) {
	smp, err := s.Evaluate(lambdaNm, angleDeg, pol)
	return smp.R, err
}

// Transmittance is a convenience wrapper returning only the transmitted
// fraction of Evaluate.
func (s Stack) Transmittance(lambdaNm, angleDeg float64, pol Polarization) (float64, error) {
	smp, err := s.Evaluate(lambdaNm, angleDeg, pol)
	return smp.T, err
}

// ---------------------------------------------------------------------------
// Scans
// ---------------------------------------------------------------------------

// Spectrum evaluates the stack at every wavelength in lambdas, at the given
// incidence angle and polarization.
func (s Stack) Spectrum(lambdas []float64, angleDeg float64, pol Polarization) ([]Sample, error) {
	out := make([]Sample, len(lambdas))
	for i, l := range lambdas {
		smp, err := s.Evaluate(l, angleDeg, pol)
		if err != nil {
			return nil, err
		}
		out[i] = smp
	}
	return out, nil
}

// SpectrumRange evaluates the stack from fromNm to toNm in steps of stepNm
// (inclusive of the endpoints; the last step is shortened when the range is not
// an integer number of steps). It is the wavelength scan used to plot a design.
func (s Stack) SpectrumRange(fromNm, toNm, stepNm, angleDeg float64, pol Polarization) ([]Sample, error) {
	lambdas, err := grid(fromNm, toNm, stepNm)
	if err != nil {
		return nil, err
	}
	return s.Spectrum(lambdas, angleDeg, pol)
}

// AngleScan evaluates the stack at every incidence angle in anglesDeg, for the
// given wavelength and polarization.
func (s Stack) AngleScan(lambdaNm float64, anglesDeg []float64, pol Polarization) ([]Sample, error) {
	out := make([]Sample, len(anglesDeg))
	for i, a := range anglesDeg {
		smp, err := s.Evaluate(lambdaNm, a, pol)
		if err != nil {
			return nil, err
		}
		out[i] = smp
	}
	return out, nil
}

// AngleScanRange evaluates the stack from fromDeg to toDeg in steps of stepDeg,
// inclusive of the endpoints. Angles may be negative (they are symmetric about
// the normal for isotropic stacks).
func (s Stack) AngleScanRange(lambdaNm, fromDeg, toDeg, stepDeg float64, pol Polarization) ([]Sample, error) {
	angles, err := grid(fromDeg, toDeg, stepDeg)
	if err != nil {
		return nil, err
	}
	return s.AngleScan(lambdaNm, angles, pol)
}

// grid builds an inclusive sample grid from a start value, an end value and a
// positive step.
func grid(from, to, step float64) ([]float64, error) {
	if !(step > 0) || math.IsInf(step, 0) {
		return nil, fmt.Errorf("coating: step %g must be positive and finite", step)
	}
	if to < from {
		return nil, fmt.Errorf("coating: range end %g is below start %g", to, from)
	}
	if math.IsNaN(from) || math.IsNaN(to) {
		return nil, errors.New("coating: range endpoints must be finite")
	}
	n := int(math.Floor((to-from)/step)) + 1
	out := make([]float64, 0, n+1)
	for i := 0; i < n; i++ {
		v := from + float64(i)*step
		if v > to {
			break
		}
		out = append(out, v)
	}
	if len(out) == 0 || out[len(out)-1] < to {
		out = append(out, to)
	}
	return out, nil
}

// BandAverage returns the mean of the requested quantity over the spectral band
// [fromNm, toNm], sampled every stepNm at the given angle and polarization.
func (s Stack) BandAverage(fromNm, toNm, stepNm, angleDeg float64, pol Polarization, qty Quantity) (float64, error) {
	samples, err := s.SpectrumRange(fromNm, toNm, stepNm, angleDeg, pol)
	if err != nil {
		return 0, err
	}
	sum := 0.0
	for _, smp := range samples {
		switch qty {
		case Reflectance:
			sum += smp.R
		case Transmittance:
			sum += smp.T
		case Absorptance:
			sum += smp.A
		default:
			return 0, fmt.Errorf("coating: invalid quantity %d", uint8(qty))
		}
	}
	return sum / float64(len(samples)), nil
}

// BandAverageR returns the mean reflectance over [fromNm, toNm].
func (s Stack) BandAverageR(fromNm, toNm, stepNm, angleDeg float64, pol Polarization) (float64, error) {
	return s.BandAverage(fromNm, toNm, stepNm, angleDeg, pol, Reflectance)
}

// BandAverageT returns the mean transmittance over [fromNm, toNm].
func (s Stack) BandAverageT(fromNm, toNm, stepNm, angleDeg float64, pol Polarization) (float64, error) {
	return s.BandAverage(fromNm, toNm, stepNm, angleDeg, pol, Transmittance)
}

// ---------------------------------------------------------------------------
// Internal field profile
// ---------------------------------------------------------------------------

// mat2 is a 2x2 complex matrix used for the interface-level bookkeeping of the
// field profile.
type mat2 [2][2]complex128

func (a mat2) mul(b mat2) mat2 {
	return mat2{
		{a[0][0]*b[0][0] + a[0][1]*b[1][0], a[0][0]*b[0][1] + a[0][1]*b[1][1]},
		{a[1][0]*b[0][0] + a[1][1]*b[1][0], a[1][0]*b[0][1] + a[1][1]*b[1][1]},
	}
}

// layerMatrix returns the characteristic matrix of a layer with phase thickness
// delta and tilted admittance eta, in the exp(-i·omega·t) convention used by
// this package:
//
//	M = [ cos(delta)        -i·sin(delta)/eta ]
//	    [ -i·eta·sin(delta)   cos(delta)       ].
//
// The sign of the imaginary parts is what makes an absorbing layer
// (n = n' + i·k, k > 0) attenuate instead of amplify; with the opposite
// convention the same physical layer would show a negative absorptance.
func layerMatrix(delta, eta complex128) mat2 {
	cd, sd := cmplx.Cos(delta), cmplx.Sin(delta)
	return mat2{{cd, complex(0, -1) * sd / eta}, {complex(0, -1) * eta * sd, cd}}
}

// FieldProfile returns the complex tangential electric-field amplitude inside
// the coating, normalised to an incident amplitude of 1: Ey for PolS and Ex for
// PolP. Depths are in nm and measured from the front surface of the first
// layer, so negative depths lie in the incidence medium and depths beyond the
// total thickness lie in the substrate.
//
// The profile is the coherent standing wave of the stack; for a substrate with
// absorption it decays with the substrate extinction. It exposes exactly the
// quantity that coating engineers plot and is used by the test suite to verify
// energy conservation from the field rather than from the matrix bookkeeping.
func (s Stack) FieldProfile(lambdaNm, angleDeg float64, pol Polarization, depthsNm []float64) ([]complex128, error) {
	if pol != PolS && pol != PolP {
		return nil, fmt.Errorf("coating: FieldProfile needs PolS or PolP, got %s", pol)
	}
	if !(lambdaNm > 0) || math.IsInf(lambdaNm, 0) {
		return nil, fmt.Errorf("coating: wavelength %g nm must be finite and positive", lambdaNm)
	}
	if math.IsNaN(angleDeg) || math.Abs(angleDeg) >= 90 {
		return nil, fmt.Errorf("coating: incidence angle %g deg must satisfy |theta| < 90", angleDeg)
	}
	useP := pol == PolP
	k0 := 2 * math.Pi / lambdaNm
	theta0 := angleDeg * math.Pi / 180

	n0 := s.Incident.IndexNm(lambdaNm)
	ns := s.Substrate.IndexNm(lambdaNm)
	alpha := n0 * complex(math.Sin(theta0), 0)

	_, _, cos0, err := admittance(n0, alpha)
	if err != nil {
		return nil, fmt.Errorf("coating: incidence medium %s at %g nm: %w", s.Incident, lambdaNm, err)
	}
	etaSubS, etaSubP, cosSub, err := admittance(ns, alpha)
	if err != nil {
		return nil, fmt.Errorf("coating: substrate %s at %g nm: %w", s.Substrate, lambdaNm, err)
	}
	etaSub := etaSubS
	if useP {
		etaSub = etaSubP
	}

	// Layer matrices, their phase factors and the admittances (n·cos(theta)).
	n := len(s.Layers)
	mats := make([]mat2, n)
	etas := make([]complex128, n)   // n·cos(theta): phase thickness factor
	admits := make([]complex128, n) // tilted admittance of the chosen polarization
	starts := make([]float64, n+1)
	for i, l := range s.Layers {
		starts[i+1] = starts[i] + l.ThicknessNm
		idx := l.Material.IndexNm(lambdaNm)
		etaS, etaP, _, aerr := admittance(idx, alpha)
		if aerr != nil {
			return nil, fmt.Errorf("coating: layer %d (%s) at %g nm: %w", i, l.Material, lambdaNm, aerr)
		}
		e := etaS
		if useP {
			e = etaP
		}
		delta := complex(k0*l.ThicknessNm, 0) * etaS
		mats[i] = layerMatrix(delta, e)
		etas[i] = etaS
		admits[i] = e
	}
	// Suffix products: suffix[i] = M_i·M_{i+1}···M_{n-1}, suffix[n] = I.
	suffix := make([]mat2, n+1)
	suffix[n] = mat2{{1, 0}, {0, 1}}
	for i := n - 1; i >= 0; i-- {
		suffix[i] = mats[i].mul(suffix[i+1])
	}

	rc, tc, _, _, err := s.polCoefficients(lambdaNm, angleDeg, pol)
	if err != nil {
		return nil, err
	}

	// Field vector (E, H) at the exit of layer i: suffix[i+1] applied to the
	// substrate state (tc, tc·etaSub).
	stateAt := func(i int) (complex128, complex128) {
		e := suffix[i][0][0]*tc + suffix[i][0][1]*tc*etaSub
		h := suffix[i][1][0]*tc + suffix[i][1][1]*tc*etaSub
		return e, h
	}
	// fieldInLayer returns the tangential E amplitude at depth z measured from
	// the front surface of layer i. The characteristic matrix of the remaining
	// distance (d_i - z) relates that plane to the back of the layer.
	fieldInLayer := func(i int, z float64) complex128 {
		eOut, hOut := stateAt(i + 1)
		delta := complex(k0*(s.Layers[i].ThicknessNm-z), 0) * etas[i]
		m := layerMatrix(delta, admits[i])
		return m[0][0]*eOut + m[0][1]*hOut
	}

	total := starts[n]
	out := make([]complex128, len(depthsNm))
	for k, z := range depthsNm {
		if math.IsNaN(z) || math.IsInf(z, 0) {
			return nil, fmt.Errorf("coating: depth %g nm must be finite", z)
		}
		switch {
		case z < 0:
			// Incidence medium: the unit incident wave exp(+i·k_z·z) plus the
			// reflected one exp(-i·k_z·z).
			ph := complex(k0*z, 0) * cos0
			out[k] = cmplx.Exp(complex(0, 1)*ph) + rc*cmplx.Exp(-complex(0, 1)*ph)
		case z >= total:
			ph := complex(k0*(z-total), 0) * cosSub
			out[k] = tc * cmplx.Exp(complex(0, 1)*ph)
		default:
			i := len(s.Layers) - 1
			for j := range s.Layers {
				if z < starts[j+1] {
					i = j
					break
				}
			}
			out[k] = fieldInLayer(i, z-starts[i])
		}
	}
	return out, nil
}
