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
//
// Every template returns a manufacturable design: no layer is thinner than
// templateMinThicknessNm. The optimiser itself is allowed to annihilate a
// layer (Options.MinThicknessNm = 0, what the zero value selects), which the
// first version of the beam-splitter template exploited and returned a 0 nm
// layer; a template that ends up there is re-optimised with the floor instead
// (see optimizeDesign). Search budgets and target grids are fixed by the
// templates and are documented with the specification each one promises.

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
// single layer available from the library. Its reflectance is the analytic
// residual of that material, R = ((n1^2-n0·ns)/(n1^2+n0·ns))^2, which the test
// suite checks exactly: 1.26% for MgF2 on BK7 (the well-known limit of a
// one-layer AR, since the ideal index 1.2323 is not available), 0.098% on
// SF11 and 0.25% for a cemented 1.52/1.62 glass pair. SingleLayerARExact is
// the machine-precision-zero reference design with the ideal index.
func SingleLayerAR(inc, substrate Material, lambdaNm float64, opts Options) (Stack, *Result, error) {
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
	res, err := optimizeDesign(start, targets, withTemplateDefaults(opts, 100, 0))
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
// that keep the dip smooth.
//
// The returned design has R < 0.1% at the design wavelength across the visible
// (measured 2.5e-5 on BK7 at 550 nm, 9.1e-5 at 420 nm and 2.0e-4 at 680 nm);
// the test suite runs the template at 420, 550 and 680 nm and verifies
// R < 0.1% at each design wavelength. The dip is narrow by construction (the V
// shape): 5% away from the design wavelength the reflectance is already of the
// order of 1%, and far from the dip a two-layer design reflects more than the
// bare substrate - use BroadbandAR when a low average reflectance over a band
// is wanted.
func VCoat(inc, substrate Material, lambdaNm float64, opts Options) (Stack, *Result, error) {
	if !(lambdaNm > 0) || math.IsInf(lambdaNm, 0) {
		return Stack{}, nil, fmt.Errorf("coating: VCoat needs a positive design wavelength, got %g nm", lambdaNm)
	}
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
	res, err := optimizeDesign(start, targets, withTemplateDefaults(opts, 0, 3))
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
// the targets ask for a small logarithmic reflectance on a grid of roughly
// 10 nm across the band and the optimiser is run with random restarts. The
// promised specification, verified by the tests on BK7 over 420-680 nm, is a
// mean reflectance below 1% (measured 0.70% on a 1 nm grid, with no point
// above 5%).
//
// The target grid is deliberately denser than the ~20 nm grid of a first
// version: a coarse grid lets the optimiser leave narrow reflectance spikes
// between the sample points, which then show up in the measured spectrum.
func BroadbandAR(inc, substrate Material, fromNm, toNm float64, opts Options) (Stack, *Result, error) {
	if !(fromNm > 0) || !(toNm > fromNm) || math.IsInf(toNm, 0) {
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
	targets, err := BandTargets(fromNm, toNm, (toNm-fromNm)/25, 0, PolUnpolarized, Reflectance, 1e-4, 1, true)
	if err != nil {
		return Stack{}, nil, err
	}
	res, err := optimizeDesign(start, targets, withTemplateDefaults(opts, 0, 4))
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// HighReflector designs a quarter-wave high-reflectance stack (H L)^pairs at
// lambdaNm using the default high/low material pair for that wavelength
// (Ta2O5/SiO2 in the near infrared, TiO2/SiO2 in the visible). At the design
// wavelength the admittance transformation of a quarter-wave stack gives
// Y = (nH/nL)^(2·pairs)·ns and therefore R = ((n0-Y)/(n0+Y))^2: 14 pairs of
// Ta2O5/SiO2 on BK7 at 1064 nm reach an analytic peak of 99.975%, of which
// the optimised design keeps more than 99.9% (the test suite verifies the
// promised R > 99.9% at the design wavelength and at +-5% of it).
func HighReflector(inc, substrate Material, lambdaNm float64, pairs int, opts Options) (Stack, *Result, error) {
	high, low := defaultHRPair(lambdaNm)
	return HighReflectorPair(inc, substrate, lambdaNm, pairs, high, low, opts)
}

// HighReflectorPair is HighReflector with an explicit material pair.
func HighReflectorPair(inc, substrate Material, lambdaNm float64, pairs int, high, low Material, opts Options) (Stack, *Result, error) {
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
	// peak target asks for R = 1, which a linear residual treats as "as close
	// to unity as possible": a target of 0.9999 would instead penalise a design
	// for being better than 0.9999 and pull the peak back down. The side
	// targets only ask that the stopband top does not collapse away from the
	// design wavelength.
	targets := []Target{
		RTarget(lambdaNm, 1.0, 1.0).LinearScale(),
		RTarget(lambdaNm*0.95, 0.9999, 0.1).LinearScale(),
		RTarget(lambdaNm*1.05, 0.9999, 0.1).LinearScale(),
	}
	res, err := optimizeDesign(start, targets, withTemplateDefaults(opts, 80, 0))
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

// The edge-filter templates (LongWavePass, ShortWavePass) fix their own search
// budget and target spacing; they are validated on the visible band, where the
// quarter-wave zone of the TiO2/SiO2 pair is wide enough for the promised
// blocked band (see quarterWaveZoneCentre).
const (
	// edgePairs is the number of quarter-wave (H L) pairs of the starting
	// design.
	edgePairs = 10
	// lwpStepNm and swpStepNm are the wavelength spacings of the target grids
	// of LongWavePass and ShortWavePass. Both have to resolve the passband
	// ripple of the quarter-wave stack (a period of a few tens of nanometres);
	// the spacings were calibrated by scanning 6-20 nm over the target recipe:
	// 10 nm (LWP) and 8 nm (SWP) gave the best measured passband (mean
	// transmittance 95.3% / 95.1% at 550 nm), while the 20 nm grid of the first
	// version let the optimiser hide ripple dips between the sample points and
	// reached only 87% / 88%.
	lwpStepNm = 10.0
	swpStepNm = 8.0
	// edgePassWeight is the weight of the passband transmittance goals (the
	// blocked-band goals keep weight 1).
	edgePassWeight = 2.0
	// edgeMaxIter and edgeRestarts are the default search budget of the edge
	// filters; pass MaxIter/Restarts in opts to override.
	edgeMaxIter  = 100
	edgeRestarts = 1
)

// templateMinThicknessNm is the thinnest layer a template may return. A layer
// of a few tenths of a nanometre is not a manufacturable film and a layer of
// exactly zero is not a layer at all, so the templates never hand one back.
const templateMinThicknessNm = 1.0

// withTemplateDefaults applies the template's search budget to every field of
// opts the caller left unset (zero or negative), then fills the remaining
// package defaults. A non-positive maxIter/restarts means "keep the package
// default", so a template without an opinion about the budget passes 0.
func withTemplateDefaults(opts Options, maxIter, restarts int) Options {
	if maxIter > 0 && opts.MaxIter <= 0 {
		opts.MaxIter = maxIter
	}
	if restarts >= 0 && opts.Restarts <= 0 {
		opts.Restarts = restarts
	}
	return opts.withDefaults()
}

// thinnestLayerNm returns the thickness of the thinnest layer of a stack, or
// +inf for a stack without layers.
func thinnestLayerNm(s Stack) float64 {
	min := math.Inf(1)
	for _, l := range s.Layers {
		if l.ThicknessNm < min {
			min = l.ThicknessNm
		}
	}
	return min
}

// optimizeDesign runs the optimiser for a template and guarantees that the
// returned design is manufacturable. The zero value of Options.MinThicknessNm
// lets the optimiser annihilate a layer (its documented meaning), and the
// unconstrained optimum of the beam-splitter template really did contain a
// layer of 0 nm; rather than forbidding the freedom up front - which costs
// design quality, the long-wave-pass filter loses most of its passband
// flatness when the floor is applied to the search - the unconstrained design
// is accepted when every layer is at least templateMinThicknessNm thick, and
// otherwise the design is re-optimised from the same start with the floor in
// place. Both paths are deterministic, so the template stays reproducible.
func optimizeDesign(start Stack, targets []Target, opts Options) (*Result, error) {
	res, err := Optimize(start, targets, opts)
	if err != nil {
		return nil, err
	}
	if thinnestLayerNm(res.Stack) >= templateMinThicknessNm {
		return res, nil
	}
	if opts.MinThicknessNm < templateMinThicknessNm {
		opts.MinThicknessNm = templateMinThicknessNm
	}
	return Optimize(start, targets, opts)
}

// quarterWaveZoneCentre converts the requested stopband edge cutoffNm into the
// design wavelength of the quarter-wave stack whose high-reflectance zone ends
// (longPass) or starts (shortPass) at the cutoff, and refuses material pairs
// whose zone cannot cover the promised blocked band. The zone of a stack of
// indices nH > nL has the half-width dg = (2/pi)·asin((nH-nL)/(nH+nL)) around
// its centre, so covering the blocked band [0.80, 0.98]·cutoff (long pass) or
// [1.02, 1.25]·cutoff (short pass) requires dg >= 1/9; TiO2/SiO2 satisfies
// that throughout the visible (dg = 0.15 at 550 nm) whereas the near-infrared
// pair Ta2O5/SiO2 reaches only dg = 0.105 at 1064 nm.
func quarterWaveZoneCentre(cutoffNm float64, high, low Material, longPass bool) (float64, error) {
	nH, nL := high.N(cutoffNm), low.N(cutoffNm)
	dg := quarterWaveBandHalfWidth(nH, nL)
	if dg < 1.0/9.0 {
		return 0, fmt.Errorf("coating: the %s/%s quarter-wave zone half-width dg = %.4f is too narrow to cover the blocked band at %g nm (dg >= 1/9 is required); choose a material pair with a larger index ratio", high, low, dg, cutoffNm)
	}
	if longPass {
		return cutoffNm * (1 - dg), nil
	}
	return cutoffNm * (1 + dg), nil
}

// edgeTargets builds the target set of an edge filter: blocked is the spectral
// region that must not be transmitted, passing the region that must be, both
// sampled every stepNm. The blocked goals use the logarithmic scale so that
// the relative rejection matters (1e-3 and 1e-4 are one decade apart), the
// passing goals are linear around T = 0.99 with the given weight.
func edgeTargets(blocked0, blocked1, passing0, passing1, stepNm, passWeight float64) ([]Target, error) {
	blocked, err := BandTargets(blocked0, blocked1, stepNm, 0, PolUnpolarized, Transmittance, 1e-3, 1, true)
	if err != nil {
		return nil, err
	}
	passing, err := BandTargets(passing0, passing1, stepNm, 0, PolUnpolarized, Transmittance, 0.99, passWeight, false)
	if err != nil {
		return nil, err
	}
	return append(blocked, passing...), nil
}

// edgeStack builds the starting design of an edge filter: a half-wave matching
// layer of the low-index material on the air side, followed by a quarter-wave
// stack (H L)^pairs designed at lambdaNm. The caller scales the stack layers to
// the stopband centre while keeping the matching layer at the cutoff scale.
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
// filter for normal incidence: wavelengths below cutoffNm are reflected,
// wavelengths above it pass. The starting design is a quarter-wave stack of
// (H L) pairs whose high-reflectance zone ends at the cutoff plus a matching
// layer on the air side; all thicknesses are then optimised against dense
// target grids (lwpStepNm) in both bands.
//
// The returned design satisfies, and the test suite verifies on a 1 nm grid:
//
//	blocked band [0.80, 0.98]·cutoff: mean T < 1%, T < 3% everywhere,
//	pass band [1.03, 1.25]·cutoff: mean T > 90%, T > 85% everywhere at
//	550 nm (measured 95.3% / 91.2%); at 632.8 nm the pass band still keeps
//	mean T > 85% and T > 60% everywhere (measured 88.4% / 67.5%), which is
//	the promise verified for that cutoff.
//
// The blocked band has to lie inside the quarter-wave high-reflectance zone
// of the material pair (dg >= 1/9), which TiO2/SiO2 satisfies throughout the
// visible; a cutoff whose pair is too narrow (e.g. the Ta2O5/SiO2 far-infrared
// default) is refused with an error instead of returning a leaky design.
func LongWavePass(inc, substrate Material, cutoffNm float64, opts Options) (Stack, *Result, error) {
	if !(cutoffNm > 0) || math.IsInf(cutoffNm, 0) {
		return Stack{}, nil, fmt.Errorf("coating: LongWavePass needs a positive cutoff wavelength, got %g nm", cutoffNm)
	}
	high, low := defaultHRPair(cutoffNm)
	lam0, err := quarterWaveZoneCentre(cutoffNm, high, low, true)
	if err != nil {
		return Stack{}, nil, err
	}
	// Scale the quarter-wave stack from the cutoff to the stopband centre,
	// keeping the matching layer at the cutoff scale.
	scale := lam0 / cutoffNm
	start, err := edgeStack(inc, substrate, cutoffNm, edgePairs, high, low, low)
	if err != nil {
		return Stack{}, nil, err
	}
	for i := range start.Layers {
		if i > 0 {
			start.Layers[i].ThicknessNm *= scale
		}
	}
	targets, err := edgeTargets(0.80*cutoffNm, 0.98*cutoffNm, 1.03*cutoffNm, 1.25*cutoffNm, lwpStepNm, edgePassWeight)
	if err != nil {
		return Stack{}, nil, err
	}
	res, err := optimizeDesign(start, targets, withTemplateDefaults(opts, edgeMaxIter, edgeRestarts))
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// ShortWavePass designs a short-wave-pass (long-wavelength blocking) edge
// filter for normal incidence: wavelengths above cutoffNm are reflected,
// wavelengths below it pass. The starting design is a quarter-wave stack of
// (H L) pairs whose high-reflectance zone starts at the cutoff plus a matching
// layer on the air side; all thicknesses are then optimised against dense
// target grids in both bands.
//
// The returned design satisfies, and the test suite verifies on a 1 nm grid:
//
//	pass band [0.80, 0.97]·cutoff: mean T > 90%, T > 70% everywhere,
//	blocked band [1.02, 1.25]·cutoff: mean T < 1%, T < 2% everywhere
//	(measured 95.1% / 76.0% and 0.15% / 0.77% at 550 nm; 95.0% / 84.0%
//	and 0.17% / 1.14% at 632.8 nm).
//
// The blocked band has to lie inside the quarter-wave high-reflectance zone of
// the material pair (dg >= 1/9, see quarterWaveZoneCentre); a cutoff whose
// pair is too narrow is refused with an error.
func ShortWavePass(inc, substrate Material, cutoffNm float64, opts Options) (Stack, *Result, error) {
	if !(cutoffNm > 0) || math.IsInf(cutoffNm, 0) {
		return Stack{}, nil, fmt.Errorf("coating: ShortWavePass needs a positive cutoff wavelength, got %g nm", cutoffNm)
	}
	high, low := defaultHRPair(cutoffNm)
	lam0, err := quarterWaveZoneCentre(cutoffNm, high, low, false)
	if err != nil {
		return Stack{}, nil, err
	}
	scale := lam0 / cutoffNm
	start, err := edgeStack(inc, substrate, cutoffNm, edgePairs, high, low, low)
	if err != nil {
		return Stack{}, nil, err
	}
	for i := range start.Layers {
		if i > 0 {
			start.Layers[i].ThicknessNm *= scale
		}
	}
	targets, err := edgeTargets(1.02*cutoffNm, 1.25*cutoffNm, 0.80*cutoffNm, 0.97*cutoffNm, swpStepNm, edgePassWeight)
	if err != nil {
		return Stack{}, nil, err
	}
	res, err := optimizeDesign(start, targets, withTemplateDefaults(opts, edgeMaxIter, edgeRestarts))
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}

// Splitter50 designs a 50/50 dielectric beam-splitter coating for use at normal
// incidence: the reflectance is driven to 0.5 over a band around lambdaNm
// (+-10% of the wavelength are targeted). The starting design is a
// Ta2O5/SiO2 stack of six layers at quarter-wave and half-wave optical
// thickness. The test suite verifies, on a 1 nm grid, |R-0.5| < 0.02 at the
// design wavelength (measured 0.0072 at 550 nm and 0.0145 at 632.8 nm), a band
// mean of [0.9, 1.1]·lambda within 1% of 0.5 (measured 0.4995 / 0.5017) and R
// within 0.5+-0.05 everywhere in that band (measured 0.4910..0.5076 at 550 nm
// and 0.4731..0.5154 at 632.8 nm).
func Splitter50(inc, substrate Material, lambdaNm float64, opts Options) (Stack, *Result, error) {
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
	res, err := optimizeDesign(start, targets, withTemplateDefaults(opts, 0, 4))
	if err != nil {
		return Stack{}, nil, err
	}
	return res.Stack, res, nil
}
