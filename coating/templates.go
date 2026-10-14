package coating

import (
	"fmt"
	"math"
)

// This file contains the built-in template designers: each function states a
// physical target set, builds an initial layer stack from the material library
// and returns the optimised design together with the optimiser report. They
// are the entry points used by the CLI, by the GUI and by the test suite, so
// their promises are written down in the doc comments and verified in the
// tests.

// ARCandidates are the coating materials a single-layer antireflection design
// may choose from, ordered by index. Films that evaporate well and are
// transparent over the visible are preferred.
func ARCandidates() []Material {
	return []Material{MgF2, SiO2, Al2O3, Y2O3, HfO2, ZrO2, Ta2O5, TiO2}
}

// HRCandidates are the high-index materials of a quarter-wave reflector.
func HRCandidates() []Material {
	return []Material{TiO2, Ta2O5, ZrO2, HfO2, Y2O3, Al2O3}
}

// QWOTNm returns the quarter-wave optical thickness lambda/(4·n') of a material
// at lambdaNm, in nm: the layer thickness that turns the material into a
// quarter-wave layer, the building block of reflectors and V-coats.
func QWOTNm(m Material, lambdaNm float64) (float64, error) {
	n := m.N(lambdaNm)
	if !(n > 0) || math.IsInf(n, 0) {
		return 0, fmt.Errorf("coating: material %s has non-physical index %g at %g nm", m, n, lambdaNm)
	}
	return lambdaNm / (4 * n), nil
}

// closestMaterial returns the candidate whose real index at lambdaNm is
// nearest to targetN.
func closestMaterial(candidates []Material, targetN, lambdaNm float64) Material {
	best := candidates[0]
	bestErr := math.Inf(1)
	for _, m := range candidates {
		n := m.N(lambdaNm)
		if !(n > 0) {
			continue
		}
		if e := math.Abs(n - targetN); e < bestErr {
			best, bestErr = m, e
		}
	}
	return best
}

// SingleLayerAR designs the classic single quarter-wave antireflection coating
// for the incident medium / substrate pair at lambdaNm. The analytic optimum
// for a lossless layer at normal incidence is n1 = sqrt(n0·ns) with thickness
// lambda/(4·n1), which would make the reflectance vanish; the material closest
// to that index is taken from ARCandidates and its thickness is then refined by
// a one-variable damped least squares fit, so the reported design is the best
// single layer available from the library (for MgF2 on BK7 the residual
// reflectance is about 1.2%, the well-known limit of a one-layer AR).
func SingleLayerAR(inc, substrate Material, lambdaNm float64, opts Options) (Stack, *Result, error) {
	opts = opts.withDefaults()
	n0, ns := inc.N(lambdaNm), substrate.N(lambdaNm)
	if !(n0 > 0) || !(ns > 0) {
		return Stack{}, nil, fmt.Errorf("coating: SingleLayerAR needs lossless media with positive index (n0=%g, ns=%g)", n0, ns)
	}
	m := closestMaterial(ARCandidates(), math.Sqrt(n0*ns), lambdaNm)
	d, err := QWOTNm(m, lambdaNm)
	if err != nil {
		return Stack{}, nil, err
	}
	start := NewStack(inc, substrate, Layer{Material: m, ThicknessNm: d})
	targets := []Target{RTarget(lambdaNm, 1e-6, 1)}
	sub := opts
	sub.Restarts = 0
	sub.MaxIter = 100
	res, err := Optimize(start, targets, sub)
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// SingleLayerARExact returns the ideal single-layer antireflection coating: a
// lossless synthetic layer with the exact optimum index n1 = sqrt(n0·ns) and
// quarter-wave thickness lambda/(4·n1). At normal incidence the amplitude
// reflection vanishes identically, so R at lambdaNm is zero to machine
// precision; use it as a reference design (and in tests), not as a coating
// that can be manufactured.
func SingleLayerARExact(inc, substrate Material, lambdaNm float64) (Stack, error) {
	n0, ns := inc.N(lambdaNm), substrate.N(lambdaNm)
	if !(n0 > 0) || !(ns > 0) {
		return Stack{}, fmt.Errorf("coating: SingleLayerARExact needs positive indices (n0=%g, ns=%g)", n0, ns)
	}
	n1 := math.Sqrt(n0 * ns)
	return NewStack(inc, substrate, Layer{
		Material:    NewMaterial("ideal n1=sqrt(n0*ns)", ConstIndex{N: n1}),
		ThicknessNm: lambdaNm / (4 * n1),
	}), nil
}

// quarterWaveBandHalfWidth returns the half-width dg of the high-reflectance
// zone of a quarter-wave stack of indices nH > nL, in the normalised frequency
// g = lambda0/lambda: dg = (2/pi)·asin((nH-nL)/(nH+nL)). The zone spans
// g = 1-+dg, i.e. the wavelengths lambda0/(1+dg) .. lambda0/(1-dg).
func quarterWaveBandHalfWidth(nH, nL float64) float64 {
	if !(nH > 0) || !(nL > 0) || nH <= nL {
		return 0
	}
	return (2 / math.Pi) * math.Asin((nH-nL)/(nH+nL))
}

// VCoat designs a two-layer "V-coat" whose reflectance dips to (almost) zero
// at lambdaNm for normal incidence. Two quarter-wave layers with indices n1
// (outside) and n2 (against the substrate) reflect nothing when n1^2·ns =
// n0·n2^2; the pair of library materials closest to that condition is used as
// the starting design and both thicknesses are optimised against a
// logarithmic zero-reflectance target plus two low-reflectance side targets
// that keep the dip smooth. The returned design meets R < 0.3% at lambdaNm,
// the specification verified by the test suite.
func VCoat(inc, substrate Material, lambdaNm float64, opts Options) (Stack, *Result, error) {
	opts = opts.withDefaults()
	n0, ns := inc.N(lambdaNm), substrate.N(lambdaNm)
	if !(n0 > 0) || !(ns > 0) {
		return Stack{}, nil, fmt.Errorf("coating: VCoat needs positive indices (n0=%g, ns=%g)", n0, ns)
	}
	low, high := vCoatPair(n0, ns, lambdaNm)
	dl, err := QWOTNm(low, lambdaNm)
	if err != nil {
		return Stack{}, nil, err
	}
	dh, err := QWOTNm(high, lambdaNm)
	if err != nil {
		return Stack{}, nil, err
	}
	start := NewStack(inc, substrate,
		Layer{Material: low, ThicknessNm: dl},
		Layer{Material: high, ThicknessNm: dh},
	)
	targets := []Target{
		RTarget(lambdaNm, 1e-5, 1.0),
		RTarget(lambdaNm*0.95, 0.005, 0.15),
		RTarget(lambdaNm*1.05, 0.005, 0.15),
	}
	sub := opts
	if sub.Restarts == 0 {
		sub.Restarts = 3
	}
	res, err := Optimize(start, targets, sub)
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// vCoatPair picks the low/high index material pair whose quarter-wave
// reflectance is closest to zero, i.e. which minimises the relative mismatch of
// n1^2·ns = n0·n2^2.
func vCoatPair(n0, ns, lambdaNm float64) (low, high Material) {
	cands := ARCandidates()
	low, high = cands[0], cands[len(cands)-1]
	best := math.Inf(1)
	for _, a := range cands {
		na := a.N(lambdaNm)
		for _, b := range cands {
			nb := b.N(lambdaNm)
			if !(na > 0 && nb > 0) || na >= nb {
				continue
			}
			den := n0 * nb * nb
			if den <= 0 {
				continue
			}
			score := math.Abs(na*na*ns-den) / den
			if score < best {
				best, low, high = score, a, b
			}
		}
	}
	return low, high
}

// BroadbandAR designs a multilayer antireflection coating for the spectral band
// [fromNm, toNm] at normal incidence. The starting design is a six-layer
// Ta2O5/SiO2 stack graded from the substrate (high index) to air (low index);
// the targets ask for a small logarithmic reflectance on a grid across the band
// and the optimiser is run with random restarts. The promised specification,
// verified by the tests on BK7 over 420-680 nm, is a mean reflectance below
// 1%.
func BroadbandAR(inc, substrate Material, fromNm, toNm float64, opts Options) (Stack, *Result, error) {
	opts = opts.withDefaults()
	if !(fromNm > 0) || !(toNm > fromNm) {
		return Stack{}, nil, fmt.Errorf("coating: BroadbandAR needs 0 < from < to, got %g..%g nm", fromNm, toNm)
	}
	lh, err := QWOTNm(Ta2O5, fromNm)
	if err != nil {
		return Stack{}, nil, err
	}
	ll, err := QWOTNm(SiO2, fromNm)
	if err != nil {
		return Stack{}, nil, err
	}
	lo, err := QWOTNm(MgF2, fromNm)
	if err != nil {
		return Stack{}, nil, err
	}
	start := NewStack(inc, substrate,
		Layer{Material: MgF2, ThicknessNm: lo * 1.5},
		Layer{Material: Ta2O5, ThicknessNm: lh * 1.2},
		Layer{Material: SiO2, ThicknessNm: ll * 0.6},
		Layer{Material: Ta2O5, ThicknessNm: lh * 0.5},
		Layer{Material: SiO2, ThicknessNm: ll * 1.5},
		Layer{Material: Ta2O5, ThicknessNm: lh * 2.5},
	)
	step := (toNm - fromNm) / 12
	targets, err := BandTargets(fromNm, toNm, step, 0, PolUnpolarized, Reflectance, 1e-4, 1, true)
	if err != nil {
		return Stack{}, nil, err
	}
	sub := opts
	if sub.Restarts == 0 {
		sub.Restarts = 4
	}
	res, err := Optimize(start, targets, sub)
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// HighReflector designs a quarter-wave high-reflectance stack (H L)^pairs at
// lambdaNm using the default high/low material pair for that wavelength
// (Ta2O5/SiO2 in the near infrared, TiO2/SiO2 in the visible). At the design
// wavelength the admittance transformation of a quarter-wave stack gives
// Y = (nH/nL)^(2·pairs)·ns and therefore R = ((n0-Y)/(n0+Y))^2, so 14 pairs on
// BK7 at 1064 nm reach R > 99.99% (the tests require the promised R > 99.9%).
func HighReflector(inc, substrate Material, lambdaNm float64, pairs int, opts Options) (Stack, *Result, error) {
	high, low := defaultHRPair(lambdaNm)
	return HighReflectorPair(inc, substrate, lambdaNm, pairs, high, low, opts)
}

// HighReflectorPair is HighReflector with an explicit material pair.
func HighReflectorPair(inc, substrate Material, lambdaNm float64, pairs int, high, low Material, opts Options) (Stack, *Result, error) {
	opts = opts.withDefaults()
	if pairs < 1 {
		return Stack{}, nil, fmt.Errorf("coating: HighReflector needs at least one pair, got %d", pairs)
	}
	dh, err := QWOTNm(high, lambdaNm)
	if err != nil {
		return Stack{}, nil, err
	}
	dl, err := QWOTNm(low, lambdaNm)
	if err != nil {
		return Stack{}, nil, err
	}
	layers := make([]Layer, 0, 2*pairs)
	for i := 0; i < pairs; i++ {
		layers = append(layers, Layer{Material: high, ThicknessNm: dh}, Layer{Material: low, ThicknessNm: dl})
	}
	start := NewStack(inc, substrate, layers...)
	// High-reflectance goals are linear (the log scale is flat as R -> 1). The
	// design wavelength carries the weight; the two side targets only ask that
	// the stopband top does not become a needle.
	targets := []Target{
		RTarget(lambdaNm, 0.9999, 1.0).LinearScale(),
		RTarget(lambdaNm*0.95, 0.999, 0.2).LinearScale(),
		RTarget(lambdaNm*1.05, 0.999, 0.2).LinearScale(),
	}
	sub := opts
	sub.MaxIter = 80
	sub.Restarts = 0
	res, err := Optimize(start, targets, sub)
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// defaultHRPair returns the standard high/low index pair for a wavelength:
// titanium dioxide and silica in the visible, tantalum pentoxide and silica in
// the near infrared (Ta2O5 has the lower loss and the wider transparency range
// there).
func defaultHRPair(lambdaNm float64) (high, low Material) {
	if lambdaNm >= 900 {
		return Ta2O5, SiO2
	}
	return TiO2, SiO2
}

// edgeTargets builds the target set of an edge filter: blocked band targets and
// passing band targets. blocked is the spectral region that must not be
// transmitted, passing the region that must be.
func edgeTargets(blocked0, blocked1, passing0, passing1 float64) ([]Target, error) {
	blocked, err := BandTargets(blocked0, blocked1, (blocked1-blocked0)/7, 0, PolUnpolarized, Transmittance, 1e-3, 1, true)
	if err != nil {
		return nil, err
	}
	passing, err := BandTargets(passing0, passing1, (passing1-passing0)/6, 0, PolUnpolarized, Transmittance, 0.99, 2, false)
	if err != nil {
		return nil, err
	}
	return append(blocked, passing...), nil
}

// edgeStack builds the starting design of an edge filter: a quarter-wave stack
// (H L)^pairs whose high-reflectance zone ends (long pass) or starts (short
// pass) at the cutoff, followed by antireflection matching on the air side.
func edgeStack(inc, substrate Material, lambdaNm float64, pairs int, high, low Material, matchLow Material) (Stack, error) {
	dh, err := QWOTNm(high, lambdaNm)
	if err != nil {
		return Stack{}, err
	}
	dl, err := QWOTNm(low, lambdaNm)
	if err != nil {
		return Stack{}, err
	}
	layers := make([]Layer, 0, 2*pairs+2)
	if matchLow != (Material{}) {
		dm, err := QWOTNm(matchLow, lambdaNm)
		if err != nil {
			return Stack{}, err
		}
		layers = append(layers, Layer{Material: matchLow, ThicknessNm: 0.5 * dm})
	}
	for i := 0; i < pairs; i++ {
		layers = append(layers, Layer{Material: high, ThicknessNm: dh}, Layer{Material: low, ThicknessNm: dl})
	}
	return NewStack(inc, substrate, layers...), nil
}

// LongWavePass designs a long-wave-pass (short-wavelength blocking) edge
// filter: wavelengths below cutoffNm are reflected (T < 1% over
// [0.8·cutoff, 0.98·cutoff]) and wavelengths above it pass (T > 90% over
// [1.03·cutoff, 1.25·cutoff]). The starting design is a quarter-wave stack
// whose high-reflectance zone ends exactly at the cutoff plus a matching layer
// on the air side, and all thicknesses are then optimised.
func LongWavePass(inc, substrate Material, cutoffNm float64, opts Options) (Stack, *Result, error) {
	opts = opts.withDefaults()
	high, low := defaultHRPair(cutoffNm)
	nH, nL := high.N(cutoffNm), low.N(cutoffNm)
	lam0 := cutoffNm * (1 - quarterWaveBandHalfWidth(nH, nL))
	if !(lam0 > 0) {
		return Stack{}, nil, fmt.Errorf("coating: LongWavePass cannot place the stopband for cutoff %g nm", cutoffNm)
	}
	// Scale the quarter-wave stack from the cutoff to the stopband centre.
	scale := lam0 / cutoffNm
	start, err := edgeStack(inc, substrate, cutoffNm, 10, high, low, low)
	if err != nil {
		return Stack{}, nil, err
	}
	for i := range start.Layers {
		if i > 0 { // keep the matching layer at the cutoff scale
			start.Layers[i].ThicknessNm *= scale
		}
	}
	targets, err := edgeTargets(0.80*cutoffNm, 0.98*cutoffNm, 1.03*cutoffNm, 1.25*cutoffNm)
	if err != nil {
		return Stack{}, nil, err
	}
	sub := opts
	if sub.Restarts == 0 {
		sub.Restarts = 2
	}
	res, err := Optimize(start, targets, sub)
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// ShortWavePass designs a short-wave-pass (long-wavelength blocking) edge
// filter: wavelengths above cutoffNm are reflected (T < 1% over
// [1.02·cutoff, 1.25·cutoff]) and wavelengths below it pass (T > 90% over
// [0.80·cutoff, 0.97·cutoff]). The starting design is a quarter-wave stack
// whose high-reflectance zone starts at the cutoff plus a matching layer on the
// air side, and all thicknesses are then optimised.
func ShortWavePass(inc, substrate Material, cutoffNm float64, opts Options) (Stack, *Result, error) {
	opts = opts.withDefaults()
	high, low := defaultHRPair(cutoffNm)
	nH, nL := high.N(cutoffNm), low.N(cutoffNm)
	lam0 := cutoffNm * (1 + quarterWaveBandHalfWidth(nH, nL))
	if !(lam0 > 0) {
		return Stack{}, nil, fmt.Errorf("coating: ShortWavePass cannot place the stopband for cutoff %g nm", cutoffNm)
	}
	scale := lam0 / cutoffNm
	start, err := edgeStack(inc, substrate, cutoffNm, 10, high, low, low)
	if err != nil {
		return Stack{}, nil, err
	}
	for i := range start.Layers {
		if i > 0 {
			start.Layers[i].ThicknessNm *= scale
		}
	}
	targets, err := edgeTargets(1.02*cutoffNm, 1.25*cutoffNm, 0.80*cutoffNm, 0.97*cutoffNm)
	if err != nil {
		return Stack{}, nil, err
	}
	sub := opts
	if sub.Restarts == 0 {
		sub.Restarts = 2
	}
	res, err := Optimize(start, targets, sub)
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// Splitter50 designs a 50/50 dielectric beam-splitter coating for use at normal
// incidence: the reflectance is driven to 0.5 over a band around lambdaNm
// (+-10% of the wavelength are targeted). The starting design is a
// Ta2O5/SiO2 stack of six layers at quarter-wave and half-wave optical
// thickness; the tests verify |R-0.5| < 0.02 at the design wavelength.
func Splitter50(inc, substrate Material, lambdaNm float64, opts Options) (Stack, *Result, error) {
	opts = opts.withDefaults()
	if !(lambdaNm > 0) {
		return Stack{}, nil, fmt.Errorf("coating: Splitter50 needs a positive wavelength, got %g", lambdaNm)
	}
	dh, err := QWOTNm(Ta2O5, lambdaNm)
	if err != nil {
		return Stack{}, nil, err
	}
	dl, err := QWOTNm(SiO2, lambdaNm)
	if err != nil {
		return Stack{}, nil, err
	}
	start := NewStack(inc, substrate,
		Layer{Material: Ta2O5, ThicknessNm: 0.5 * dh},
		Layer{Material: SiO2, ThicknessNm: 0.5 * dl},
		Layer{Material: Ta2O5, ThicknessNm: 0.5 * dh},
		Layer{Material: SiO2, ThicknessNm: 0.75 * dl},
		Layer{Material: Ta2O5, ThicknessNm: 0.75 * dh},
		Layer{Material: SiO2, ThicknessNm: 0.5 * dl},
	)
	targets := []Target{RTarget(lambdaNm, 0.5, 1).LinearScale()}
	for _, f := range []float64{0.9, 0.95, 1.05, 1.1} {
		targets = append(targets, RTarget(lambdaNm*f, 0.5, 0.5).LinearScale())
	}
	sub := opts
	if sub.Restarts == 0 {
		sub.Restarts = 4
	}
	res, err := Optimize(start, targets, sub)
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}
