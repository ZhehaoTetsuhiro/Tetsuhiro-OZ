package coating

import (
	"math"
	"testing"
)

// This file holds the end-to-end tests of the built-in template designers.
// Every threshold is the measured performance of the template with its default
// options, moved away from the measurement by a comfortable margin; the number
// measured when the threshold was set is quoted in the failure message, and
// the same specification is written down in the doc comment of the template.
// The optimiser is deterministic (fixed seed, no wall-clock input), so these
// are exact reproductions, not statistical statements.

// spectrumStats summarises a normal-incidence spectral scan of the requested
// quantity over [fromNm, toNm] on a stepNm grid: the mean, the minimum and the
// maximum sample.
func spectrumStats(t *testing.T, st Stack, fromNm, toNm, stepNm float64, qty Quantity) (mean, min, max float64) {
	t.Helper()
	samples, err := st.SpectrumRange(fromNm, toNm, stepNm, 0, PolUnpolarized)
	if err != nil {
		t.Fatalf("spectrum %.4g-%.4g nm: %v", fromNm, toNm, err)
	}
	min, max = math.Inf(1), math.Inf(-1)
	sum := 0.0
	for _, s := range samples {
		var v float64
		switch qty {
		case Reflectance:
			v = s.R
		case Transmittance:
			v = s.T
		case Absorptance:
			v = s.A
		default:
			t.Fatalf("invalid quantity %d", uint8(qty))
		}
		sum += v
		if v < min {
			min = v
		}
		if v > max {
			max = v
		}
	}
	return sum / float64(len(samples)), min, max
}

// checkEnergyConserved asserts R + T + A = 1 (to rounding noise) over a scan,
// for the two polarizations and their mean, and that the absorptance of a
// lossless design never goes negative. It is the energy-conservation contract
// of the package applied to the designs the templates return.
func checkEnergyConserved(t *testing.T, st Stack, fromNm, toNm, stepNm float64) {
	t.Helper()
	for _, pol := range []Polarization{PolS, PolP, PolUnpolarized} {
		samples, err := st.SpectrumRange(fromNm, toNm, stepNm, 0, pol)
		if err != nil {
			t.Fatalf("energy scan %.4g-%.4g nm (%s): %v", fromNm, toNm, pol, err)
		}
		for _, s := range samples {
			if d := math.Abs(s.R + s.T + s.A - 1); d > 1e-12 {
				t.Errorf("%s at %.4g nm: R+T+A-1 = %g, want within 1e-12", pol, s.LambdaNm, s.R+s.T+s.A-1)
				break
			}
			if s.A < -1e-12 {
				t.Errorf("%s at %.4g nm: A = %g is negative for a lossless design", pol, s.LambdaNm, s.A)
				break
			}
		}
	}
}

// TestSingleLayerARTemplate checks the library-material single-layer
// antireflection coating. The best available material is selected
// automatically (MgF2 for air/BK7, n = 1.3785 against the ideal sqrt(n0·ns) =
// 1.2323) and the optimiser must land on that material's analytic optimum
// R = ((n1^2 - n0·ns)/(n1^2 + n0·ns))^2 - the single-layer reflectance can not
// go below it, which is why the promise is "R < 2% at the design wavelength"
// and not the machine-precision zero of SingleLayerARExact.
func TestSingleLayerARTemplate(t *testing.T) {
	t.Parallel()
	const lambda = 550.0
	cases := []struct {
		name    string
		n0, ns  float64
		wantMat string
		maxR    float64
	}{
		{"air/BK7", 1.0, 1.5168, "MgF2", 0.02},         // measured R = 0.0125942
		{"air/SF11", 1.0, 1.7847, "MgF2", 0.002},       // measured R = 0.000983742
		{"cemented glass", 1.52, 1.62, "Al2O3", 0.005}, // measured R = 0.00251652
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			inc := NewMaterial("incident", ConstIndex{N: tc.n0})
			sub := NewMaterial("substrate", ConstIndex{N: tc.ns})
			st, res, err := SingleLayerAR(inc, sub, lambda, Options{})
			if err != nil {
				t.Fatalf("SingleLayerAR: %v", err)
			}
			if len(st.Layers) != 1 {
				t.Fatalf("design has %d layers, want 1", len(st.Layers))
			}
			if got := st.Layers[0].Material.Name; got != tc.wantMat {
				t.Errorf("selected material %q, want %q (the candidate nearest to sqrt(n0·ns))", got, tc.wantMat)
			}
			// The thickness must be the quarter-wave value of the material and
			// the reflectance must equal the analytic optimum.
			qwot, err := QWOTNm(st.Layers[0].Material, lambda)
			if err != nil {
				t.Fatal(err)
			}
			closeFloat(t, "quarter-wave thickness", st.Layers[0].ThicknessNm, qwot, 1e-6)
			n1 := st.Layers[0].Material.N(lambda)
			analytic := math.Pow((n1*n1-tc.n0*tc.ns)/(n1*n1+tc.n0*tc.ns), 2)
			got, err := st.Reflectance(lambda, 0, PolUnpolarized)
			if err != nil {
				t.Fatal(err)
			}
			closeFloat(t, "R at the design wavelength", got, analytic, 1e-9)
			if got > tc.maxR {
				t.Errorf("R = %g at the design wavelength, want <= %g", got, tc.maxR)
			}
			// The quarter-wave AR is narrow: half a design wavelength away the
			// layer is no longer matched and the reflectance is back to a few
			// per cent.
			far, err := st.Reflectance(0.5*lambda, 0, PolUnpolarized)
			if err != nil {
				t.Fatal(err)
			}
			if far < 1e-3 {
				t.Errorf("R at 0.5·lambda0 = %g, want a clearly non-zero value (the dip is narrow)", far)
			}
			// The reported merit must be the merit of the returned design.
			m, err := Merit(st, res.Targets)
			if err != nil {
				t.Fatal(err)
			}
			closeFloat(t, "reported merit", m, res.Merit, 1e-12)
			checkEnergyConserved(t, st, lambda, lambda, 1)
		})
	}
}

// TestVCoatTemplate runs the V-coat at three design wavelengths across the
// visible. Each design must null the reflectance below 0.1% at its own design
// wavelength (measured 9.1e-5 at 420 nm, 2.5e-5 at 550 nm, 2.0e-4 at 680 nm),
// keep the V shape (the shoulders at +-5% an order of magnitude above the dip)
// and stay below 15% anywhere in 400-700 nm (measured maximum 12.7%).
func TestVCoatTemplate(t *testing.T) {
	t.Parallel()
	for _, lam := range []float64{420, 550, 680} {
		st, res, err := VCoat(Air, BK7, lam, Options{})
		if err != nil {
			t.Fatalf("VCoat at %g nm: %v", lam, err)
		}
		if len(st.Layers) != 2 {
			t.Fatalf("VCoat at %g nm: %d layers, want 2", lam, len(st.Layers))
		}
		r0, err := st.Reflectance(lam, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if r0 > 1e-3 {
			t.Errorf("R(%.0f nm) = %g, want < 0.1%% (a two-layer V-coat must reach a deep null)", lam, r0)
		}
		for _, f := range []float64{0.95, 1.05} {
			r, err := st.Reflectance(lam*f, 0, PolUnpolarized)
			if err != nil {
				t.Fatal(err)
			}
			if r <= 10*r0 {
				t.Errorf("V shape at %.0f nm: R(%.2f·lambda0) = %g is not clearly above the dip %g", lam, f, r, r0)
			}
			if r > 0.02 {
				t.Errorf("shoulder at %.0f nm: R(%.2f·lambda0) = %g, want < 2%% (measured 0.45%%-1.2%%)", lam, f, r)
			}
		}
		_, _, mx := spectrumStats(t, st, 400, 700, 5, Reflectance)
		if mx > 0.15 {
			t.Errorf("V-coat at %.0f nm: maximum R over 400-700 nm = %g, want < 15%% (the two-layer V-coat is worse than the bare substrate away from the dip)", lam, mx)
		}
		m, err := Merit(st, res.Targets)
		if err != nil {
			t.Fatal(err)
		}
		closeFloat(t, "reported merit", m, res.Merit, 1e-12)
		checkEnergyConserved(t, st, 400, 700, 20)
	}
}

// TestBroadbandARTemplate checks the six-layer broadband antireflection design
// on BK7 over 420-680 nm. The promised mean reflectance is below 1%: the
// design measures 0.70% on a 1 nm grid (0.74% on a 5 nm grid) with a worst
// point of 4.6%, which is the reason the test also bounds the peak.
func TestBroadbandARTemplate(t *testing.T) {
	t.Parallel()
	st, res, err := BroadbandAR(Air, BK7, 420, 680, Options{})
	if err != nil {
		t.Fatalf("BroadbandAR: %v", err)
	}
	if len(st.Layers) != 6 {
		t.Errorf("design has %d layers, want the documented 6", len(st.Layers))
	}
	mean1, _, max1 := spectrumStats(t, st, 420, 680, 1, Reflectance)
	if mean1 >= 0.01 {
		t.Errorf("mean R over 420-680 nm = %.5f, want < 1%%", mean1)
	}
	if max1 > 0.05 {
		t.Errorf("peak R over 420-680 nm = %.5f, want < 5%% (measured 4.59%%)", max1)
	}
	mean5, _, _ := spectrumStats(t, st, 420, 680, 5, Reflectance)
	if mean5 >= 0.01 {
		t.Errorf("mean R over 420-680 nm (5 nm grid) = %.5f, want < 1%%", mean5)
	}
	if tot := st.TotalThicknessNm(); tot > 1000 {
		t.Errorf("total thickness %.1f nm, want a practical design below 1000 nm", tot)
	}
	if res.Starts != 5 {
		t.Errorf("Starts = %d, want 5 (one initial design plus the four documented restarts)", res.Starts)
	}
	if res.Merit >= res.InitialMerit {
		t.Errorf("merit did not improve: initial %g, final %g", res.InitialMerit, res.Merit)
	}
	m, err := Merit(st, res.Targets)
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "reported merit", m, res.Merit, 1e-12)
	checkEnergyConserved(t, st, 420, 680, 10)
}

// TestHighReflectorTemplate checks the 1064 nm quarter-wave high reflector on
// BK7 (14 Ta2O5/SiO2 pairs). The pure quarter-wave stack has the analytic peak
// R = ((n0-Y)/(n0+Y))^2 with Y = (nH/nL)^28·ns = 99.975%; the optimised design
// must keep more than 99.9% at the design wavelength and at +-5% of it
// (measured 99.975%, 99.935% and 99.940%).
func TestHighReflectorTemplate(t *testing.T) {
	t.Parallel()
	const lam = 1064.0
	st, res, err := HighReflector(Air, BK7, lam, 14, Options{})
	if err != nil {
		t.Fatalf("HighReflector: %v", err)
	}
	if len(st.Layers) != 28 {
		t.Fatalf("design has %d layers, want 28 (14 pairs)", len(st.Layers))
	}
	nH, nL, ns := Ta2O5.N(lam), SiO2.N(lam), BK7.N(lam)
	Y := math.Pow(nH/nL, 28) * ns
	analytic := math.Pow((1-Y)/(1+Y), 2)
	got, err := st.Reflectance(lam, 0, PolUnpolarized)
	if err != nil {
		t.Fatal(err)
	}
	if got <= 0.999 {
		t.Errorf("R(%.0f nm) = %.6f, want > 99.9%%", lam, got)
	}
	if got < analytic-1e-4 {
		t.Errorf("optimised R = %.8f is more than 1e-4 below the analytic quarter-wave peak %.8f", got, analytic)
	}
	for _, f := range []float64{0.95, 1.05} {
		r, err := st.Reflectance(lam*f, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if r <= 0.999 {
			t.Errorf("stopband: R(%.2f·lambda0) = %.6f, want > 99.9%% (measured 99.935%%/99.940%%)", f, r)
		}
	}
	// The default pair switches to the visible pair below 900 nm.
	vis, _, err := HighReflector(Air, BK7, 550, 4, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if vis.Layers[0].Material.Name != "TiO2" || vis.Layers[1].Material.Name != "SiO2" {
		t.Errorf("visible design uses %s/%s, want TiO2/SiO2", vis.Layers[0].Material.Name, vis.Layers[1].Material.Name)
	}
	rv, err := vis.Reflectance(550, 0, PolUnpolarized)
	if err != nil {
		t.Fatal(err)
	}
	if rv < 0.9 {
		t.Errorf("4-pair visible reflector: R = %.6f, want > 90%%", rv)
	}
	m, err := Merit(st, res.Targets)
	if err != nil {
		t.Fatal(err)
	}
	closeFloat(t, "reported merit", m, res.Merit, 1e-12)
	checkEnergyConserved(t, st, 1000, 1130, 10)
}

// TestLongWavePassTemplate checks the long-wave-pass (short-wavelength
// blocking) edge filter at two cutoffs: the nominal 550 nm design and 632.8 nm
// (the HeNe line), which demonstrates that the template works across the
// visible. On a 1 nm grid the blocked band [0.80, 0.98]·cutoff must keep a
// mean transmittance below 1% (measured 0.27%) and the pass band
// [1.03, 1.25]·cutoff a mean above the threshold given below (measured 95.3%
// at 550 nm and 88.4% at 632.8 nm, hence the different promise for the second
// cutoff). The transition itself is checked for position: T < 2% just below
// the cutoff and T > 85% well above it.
func TestLongWavePassTemplate(t *testing.T) {
	t.Parallel()
	cases := []struct {
		cutoff                   float64
		stopMeanMax, stopPeakMax float64
		passMeanMin, passMinMin  float64
	}{
		{550, 0.01, 0.03, 0.90, 0.85},    // measured 0.0027/0.0205 and 0.9529/0.9115
		{632.8, 0.01, 0.025, 0.85, 0.60}, // measured 0.0028/0.0171 and 0.8838/0.6750
	}
	for _, tc := range cases {
		st, res, err := LongWavePass(Air, BK7, tc.cutoff, Options{})
		if err != nil {
			t.Fatalf("LongWavePass at %g nm: %v", tc.cutoff, err)
		}
		stopMean, _, stopMax := spectrumStats(t, st, 0.80*tc.cutoff, 0.98*tc.cutoff, 1, Transmittance)
		passMean, passMin, _ := spectrumStats(t, st, 1.03*tc.cutoff, 1.25*tc.cutoff, 1, Transmittance)
		if stopMean >= tc.stopMeanMax {
			t.Errorf("cutoff %.1f nm: mean T over the blocked band = %.5f, want < %g", tc.cutoff, stopMean, tc.stopMeanMax)
		}
		if stopMax >= tc.stopPeakMax {
			t.Errorf("cutoff %.1f nm: maximum T over the blocked band = %.5f, want < %g", tc.cutoff, stopMax, tc.stopPeakMax)
		}
		if passMean <= tc.passMeanMin {
			t.Errorf("cutoff %.1f nm: mean T over the pass band = %.5f, want > %g", tc.cutoff, passMean, tc.passMeanMin)
		}
		if passMin <= tc.passMinMin {
			t.Errorf("cutoff %.1f nm: minimum T over the pass band = %.5f, want > %g", tc.cutoff, passMin, tc.passMinMin)
		}
		// Transition position: blocked just below, open well above the cutoff.
		below, err := st.Transmittance(0.95*tc.cutoff, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if below >= 0.02 {
			t.Errorf("cutoff %.1f nm: T at 0.95·cutoff = %.5f, want < 2%% (measured 0.4%%)", tc.cutoff, below)
		}
		above, err := st.Transmittance(1.20*tc.cutoff, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if above <= 0.85 {
			t.Errorf("cutoff %.1f nm: T at 1.20·cutoff = %.5f, want > 85%% (measured 96.2%%/91.8%%)", tc.cutoff, above)
		}
		m, err := Merit(st, res.Targets)
		if err != nil {
			t.Fatal(err)
		}
		closeFloat(t, "reported merit", m, res.Merit, 1e-12)
		checkEnergyConserved(t, st, 0.75*tc.cutoff, 1.3*tc.cutoff, 5)
	}
	// A cutoff whose material pair has too narrow a quarter-wave zone must be
	// refused rather than silently designed into a leaky filter.
	if _, _, err := LongWavePass(Air, BK7, 1064, Options{}); err == nil {
		t.Error("LongWavePass at 1064 nm (Ta2O5/SiO2 zone dg = 0.105 < 1/9): expected an error")
	}
	if _, _, err := LongWavePass(Air, BK7, 0, Options{}); err == nil {
		t.Error("LongWavePass with zero cutoff: expected an error")
	}
}

// TestShortWavePassTemplate is the mirror of TestLongWavePassTemplate for the
// short-wave-pass (long-wavelength blocking) edge filter. The blocked band
// [1.02, 1.25]·cutoff must stay below 2% everywhere (measured 0.77% at 550 nm
// and 1.14% at 632.8 nm) with a mean below 1%, and the pass band
// [0.80, 0.97]·cutoff above a 90% mean (measured 95.1% / 95.0%) with a
// minimum above 70% (measured 76.0% / 84.0%).
func TestShortWavePassTemplate(t *testing.T) {
	t.Parallel()
	for _, cutoff := range []float64{550, 632.8} {
		st, res, err := ShortWavePass(Air, BK7, cutoff, Options{})
		if err != nil {
			t.Fatalf("ShortWavePass at %g nm: %v", cutoff, err)
		}
		passMean, passMin, _ := spectrumStats(t, st, 0.80*cutoff, 0.97*cutoff, 1, Transmittance)
		stopMean, _, stopMax := spectrumStats(t, st, 1.02*cutoff, 1.25*cutoff, 1, Transmittance)
		if passMean <= 0.90 {
			t.Errorf("cutoff %.1f nm: mean T over the pass band = %.5f, want > 90%%", cutoff, passMean)
		}
		if passMin <= 0.70 {
			t.Errorf("cutoff %.1f nm: minimum T over the pass band = %.5f, want > 70%%", cutoff, passMin)
		}
		if stopMean >= 0.01 {
			t.Errorf("cutoff %.1f nm: mean T over the blocked band = %.5f, want < 1%%", cutoff, stopMean)
		}
		if stopMax >= 0.02 {
			t.Errorf("cutoff %.1f nm: maximum T over the blocked band = %.5f, want < 2%%", cutoff, stopMax)
		}
		above, err := st.Transmittance(1.05*cutoff, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if above >= 0.02 {
			t.Errorf("cutoff %.1f nm: T at 1.05·cutoff = %.5f, want < 2%% (measured 0.2%%/0.25%%)", cutoff, above)
		}
		below, err := st.Transmittance(0.85*cutoff, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if below <= 0.85 {
			t.Errorf("cutoff %.1f nm: T at 0.85·cutoff = %.5f, want > 85%% (measured 96.4%%/97.4%%)", cutoff, below)
		}
		m, err := Merit(st, res.Targets)
		if err != nil {
			t.Fatal(err)
		}
		closeFloat(t, "reported merit", m, res.Merit, 1e-12)
		checkEnergyConserved(t, st, 0.75*cutoff, 1.3*cutoff, 5)
	}
	if _, _, err := ShortWavePass(Air, BK7, 1064, Options{}); err == nil {
		t.Error("ShortWavePass at 1064 nm (Ta2O5/SiO2 zone dg = 0.105 < 1/9): expected an error")
	}
}

// TestSplitter50Template checks the 50/50 beam-splitter coating at 550 nm and
// 632.8 nm: |R - 0.5| < 0.02 at the design wavelength (measured 0.0106 and
// 0.0145), the mean over +-10% within 1% of 0.5 (measured 0.0015 and 0.0017)
// and no point of that band further than 0.05 from 0.5 (measured 0.4706..0.5116
// at 550 nm).
func TestSplitter50Template(t *testing.T) {
	t.Parallel()
	for _, lam := range []float64{550, 632.8} {
		st, res, err := Splitter50(Air, BK7, lam, Options{})
		if err != nil {
			t.Fatalf("Splitter50 at %g nm: %v", lam, err)
		}
		r0, err := st.Reflectance(lam, 0, PolUnpolarized)
		if err != nil {
			t.Fatal(err)
		}
		if d := math.Abs(r0 - 0.5); d > 0.02 {
			t.Errorf("at %.1f nm: |R-0.5| = %.4f, want <= 0.02", lam, d)
		}
		mean, min, max := spectrumStats(t, st, 0.9*lam, 1.1*lam, 1, Reflectance)
		if d := math.Abs(mean - 0.5); d > 0.01 {
			t.Errorf("at %.1f nm: band mean R = %.5f, want within 1%% of 0.5", lam, mean)
		}
		if min < 0.45 || max > 0.55 {
			t.Errorf("at %.1f nm: band R spans %.4f..%.4f, want inside 0.5+-0.05", lam, min, max)
		}
		// The unconstrained optimum of this template annihilates a layer; the
		// template must return the floored re-run instead.
		for i, l := range st.Layers {
			if l.ThicknessNm < templateMinThicknessNm {
				t.Errorf("at %.1f nm: layer %d is %.4f nm thick, want at least %g nm", lam, i, l.ThicknessNm, templateMinThicknessNm)
			}
		}
		m, err := Merit(st, res.Targets)
		if err != nil {
			t.Fatal(err)
		}
		closeFloat(t, "reported merit", m, res.Merit, 1e-12)
		checkEnergyConserved(t, st, 0.9*lam, 1.1*lam, 5)
	}
}
