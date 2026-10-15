package lens

import (
	"math"
	"math/rand"
	"testing"

	"tetsuhiro-oz/geo"
)

// thinBiconvex returns a zero-thickness biconvex singlet in air (object at
// infinity) of the given glass: R1 = +80 mm, R2 = −80 mm.
func thinBiconvex(glass string) *Prescription {
	return &Prescription{
		Name:      "thin " + glass,
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{RadiusMM: 80, ThicknessMM: 0, Glass: glass, SemiAperMM: 10},
			{RadiusMM: -80, ThicknessMM: 0, Glass: "", SemiAperMM: 10},
		},
	}
}

// thickBiconvex is a 8 mm thick biconvex singlet, R1 = +60, R2 = −120 mm.
func thickBiconvex(glass string) *Prescription {
	return &Prescription{
		Name:      "thick " + glass,
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{RadiusMM: 60, ThicknessMM: 8, Glass: glass, SemiAperMM: 12},
			{RadiusMM: -120, ThicknessMM: 0, Glass: "", SemiAperMM: 12},
		},
	}
}

func relDiff(a, b float64) float64 {
	return math.Abs(a-b) / math.Max(1, math.Abs(b))
}

// TestThermalZeroIdentity: at ΔT = 0 the thermal prescription has a system
// identical to the original, so every paraxial quantity and a real trace are
// unchanged to 1e-15 (relative, or absolute below 1).
func TestThermalZeroIdentity(t *testing.T) {
	worst := 0.0
	defer func() { t.Logf("max relative difference at ΔT=0 = %.3g (bound 1e-15)", worst) }()
	achro, _, err := DesignAchromat(100, "N-BK7", "F2", 10, &DesignOptions{SkipPolish: true})
	if err != nil {
		t.Fatal(err)
	}
	asph := &Prescription{
		Name:      "asphere",
		ObjectZMM: -300,
		Cards: []Card{
			{RadiusMM: 45, ConicK: -0.5, Coeffs: []float64{2e-6, -3e-9}, ThicknessMM: 6, Glass: "N-SK16", SemiAperMM: 9},
			{RadiusMM: -90, ThicknessMM: 12, SemiAperMM: 9},
			{RadiusMM: 0, ThicknessMM: 0, Glass: "fused silica", SemiAperMM: 9},
			{RadiusMM: -60, ThicknessMM: 0, SemiAperMM: 9},
		},
	}
	for _, p := range []*Prescription{thinBiconvex("N-BK7"), thickBiconvex("F2"), achro, asph} {
		hot, err := p.Thermal(0, nil)
		if err != nil {
			t.Fatalf("%s: %v", p.Name, err)
		}
		a, err := p.Paraxial(geo.Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		b, err := hot.Paraxial(geo.Wavelengthd)
		if err != nil {
			t.Fatal(err)
		}
		pairs := map[string][2]float64{
			"EFL": {a.EFL, b.EFL}, "BFL": {a.BFL, b.BFL}, "FFL": {a.FFL, b.FFL},
			"EPZ": {a.EntrancePupilZ, b.EntrancePupilZ}, "EPR": {a.EntrancePupilRadius, b.EntrancePupilRadius},
			"F#": {a.FNumber, b.FNumber}, "M": {a.Magnification, b.Magnification},
		}
		for k, v := range pairs {
			worst = math.Max(worst, relDiff(v[1], v[0]))
			if relDiff(v[1], v[0]) > 1e-15 {
				t.Errorf("%s: %s changed at ΔT=0: %.17g → %.17g", p.Name, k, v[0], v[1])
			}
		}
		sa, err := p.System()
		if err != nil {
			t.Fatal(err)
		}
		sb, err := hot.System()
		if err != nil {
			t.Fatal(err)
		}
		ra := sa.Trace(geo.NewRay(geo.V(0, 3, -50), geo.V(0, 0, 1), geo.Wavelengthd))
		rb := sb.Trace(geo.NewRay(geo.V(0, 3, -50), geo.V(0, 0, 1), geo.Wavelengthd))
		if ra.Status != rb.Status {
			t.Fatalf("%s: trace status %v → %v", p.Name, ra.Status, rb.Status)
		}
		if ra.OK() {
			for _, pair := range [][2]float64{
				{ra.Final().X, rb.Final().X}, {ra.Final().Y, rb.Final().Y}, {ra.Final().Z, rb.Final().Z},
				{ra.Direction().Y, rb.Direction().Y}, {ra.OPL(), rb.OPL()},
			} {
				worst = math.Max(worst, relDiff(pair[1], pair[0]))
				if relDiff(pair[1], pair[0]) > 1e-15 {
					t.Errorf("%s: trace changed at ΔT=0: %.17g → %.17g", p.Name, pair[0], pair[1])
				}
			}
		}
	}
}

// TestThinLensAnalyticDrift: for a thin lens in air, dEFL/dT from the closed
// form EFL·(α − n_rel'/(n_rel − 1)) — with the index derivative and the CTE
// of the thermal model — agrees with a central finite difference of the
// paraxial EFL to 1e-6 relative, at several temperatures and for every
// catalogue glass.
func TestThinLensAnalyticDrift(t *testing.T) {
	const h = 0.02 // K
	const lam = geo.Wavelengthd
	worst := 0.0
	defer func() { t.Logf("max relative |FD − analytic| of dEFL/dT = %.3g (bound 1e-6)", worst) }()
	for _, glass := range geo.GlassNames() {
		g, ok := geo.ThermalGlassByName(glass)
		if !ok {
			t.Fatalf("%s: no thermal data", glass)
		}
		base, _ := geo.Glass(glass)
		p := thinBiconvex(glass)
		for _, dT := range []float64{-40, 0, 25, 60} {
			efl := func(x float64) float64 {
				q, err := p.Thermal(x, nil)
				if err != nil {
					t.Fatal(err)
				}
				par, err := q.Paraxial(lam)
				if err != nil {
					t.Fatal(err)
				}
				return par.EFL
			}
			fd := (efl(dT+h) - efl(dT-h)) / (2 * h)
			n0 := base.N(lam)
			nRel := n0 + g.DeltaNRel(lam, dT, n0)
			dn := g.DnRelDT(lam, dT, n0)
			analytic := geo.ThinLensEFLDrift(efl(dT), nRel, dn, g.CTE, dT)
			rel := math.Abs(fd-analytic) / math.Abs(analytic)
			worst = math.Max(worst, rel)
			if rel > 1e-6 {
				t.Errorf("%s ΔT=%g: dEFL/dT FD %.10g vs analytic %.10g (rel %.2g)", glass, dT, fd, analytic, rel)
			}
		}
	}
}

// TestThickLensExactTrace: the paraxial EFL and BFL of the thermal system
// agree with an exact trace of a near-axial ray (height 0.1 mm) through the
// same thermal system to 1e-4 relative, at several temperatures.
func TestThickLensExactTrace(t *testing.T) {
	const lam = geo.Wavelengthd
	const h = 0.1 // mm
	worst := 0.0
	defer func() { t.Logf("max relative EFL/BFL exact-vs-paraxial difference = %.3g (bound 1e-4)", worst) }()
	for _, dT := range []float64{-30, 0, 40} {
		p := thickBiconvex("N-BK7")
		hot, err := p.Thermal(dT, nil)
		if err != nil {
			t.Fatal(err)
		}
		par, err := hot.Paraxial(lam)
		if err != nil {
			t.Fatal(err)
		}
		sys, err := hot.System()
		if err != nil {
			t.Fatal(err)
		}
		res := sys.Trace(geo.NewRay(geo.V(0, h, -50), geo.V(0, 0, 1), lam))
		if !res.OK() {
			t.Fatalf("ΔT=%g: trace failed: %v", dT, res.Error())
		}
		d := res.Direction()
		u := d.Y / d.Z                            // exit slope dy/dz
		eflExact := -h / u                        // paraxial relation EFL = −y/u'
		zCross := res.Final().Z - res.Final().Y/u // axis crossing after the last surface
		bflExact := zCross - res.Final().Z        // distance from the last surface
		worst = math.Max(worst, math.Abs(eflExact-par.EFL)/math.Abs(par.EFL))
		worst = math.Max(worst, math.Abs(bflExact-par.BFL)/math.Abs(par.BFL))
		if rel := math.Abs(eflExact-par.EFL) / math.Abs(par.EFL); rel > 1e-4 {
			t.Errorf("ΔT=%g: EFL exact %.8f vs paraxial %.8f (rel %.2g)", dT, eflExact, par.EFL, rel)
		}
		if rel := math.Abs(bflExact-par.BFL) / math.Abs(par.BFL); rel > 1e-4 {
			t.Errorf("ΔT=%g: BFL exact %.8f vs paraxial %.8f (rel %.2g)", dT, bflExact, par.BFL, rel)
		}
	}
}

// TestThermalSystemScaling: the geometry follows the CTE rules of the model.
// A glass element scales its radii, its gap and its aperture by 1 + α ΔT; a
// surface between two glasses takes the mean of the two CTEs; an air gap
// scales by the mount CTE only (zero for a rigid mount); the index of each
// glass rises with ΔT > 0 (relative to air).
func TestThermalSystemScaling(t *testing.T) {
	const dT = 25.0
	bk7, _ := geo.ThermalGlassByName("N-BK7")
	f2, _ := geo.ThermalGlassByName("F2")

	p := thickBiconvex("N-BK7")
	hot, err := p.Thermal(dT, nil)
	if err != nil {
		t.Fatal(err)
	}
	s := 1 + bk7.CTE*dT
	checkEq := func(name string, got, want float64) {
		t.Helper()
		if relDiff(got, want) > 1e-14 {
			t.Errorf("%s = %.15g, want %.15g", name, got, want)
		}
	}
	checkEq("glass thickness", hot.Cards[0].ThicknessMM, p.Cards[0].ThicknessMM*s)
	checkEq("front radius", hot.Cards[0].RadiusMM, p.Cards[0].RadiusMM*s)
	checkEq("rear radius (glass→air)", hot.Cards[1].RadiusMM, p.Cards[1].RadiusMM*s)
	checkEq("front aperture", hot.Cards[0].SemiAperMM, p.Cards[0].SemiAperMM*s)

	// Cemented surface: mean of the two CTEs; air gap: mount CTE (rigid = 0).
	cem := &Prescription{
		ObjectZMM: math.Inf(-1),
		Cards: []Card{
			{RadiusMM: 80, ThicknessMM: 5, Glass: "N-BK7", SemiAperMM: 10},
			{RadiusMM: -40, ThicknessMM: 3, Glass: "F2", SemiAperMM: 10},
			{RadiusMM: -150, ThicknessMM: 40, SemiAperMM: 10},
			{RadiusMM: 0, ThicknessMM: 0, SemiAperMM: 10},
		},
	}
	ch, err := cem.Thermal(dT, nil)
	if err != nil {
		t.Fatal(err)
	}
	checkEq("cemented radius (mean CTE)", ch.Cards[1].RadiusMM, -40*(1+0.5*(bk7.CTE+f2.CTE)*dT))
	checkEq("F2 gap", ch.Cards[1].ThicknessMM, 3*(1+f2.CTE*dT))
	checkEq("air gap, rigid mount", ch.Cards[2].ThicknessMM, 40)
	mounted, err := cem.Thermal(dT, &ThermalOptions{MountCTE: 1e-5})
	if err != nil {
		t.Fatal(err)
	}
	checkEq("air gap, mount CTE 1e-5", mounted.Cards[2].ThicknessMM, 40*(1+1e-5*dT))
	// Card 2's front surface touches F2 (glass before it): F2's CTE. Card 3's
	// surface touches only air on both sides: the mount CTE.
	checkEq("surface after F2 takes F2 CTE", ch.Cards[2].RadiusMM, -150*(1+f2.CTE*dT))
	air := &Prescription{ObjectZMM: math.Inf(-1), Cards: []Card{
		{RadiusMM: 80, ThicknessMM: 5, Glass: "N-BK7", SemiAperMM: 10},
		{RadiusMM: -60, ThicknessMM: 20, SemiAperMM: 10},
		{RadiusMM: -60, ThicknessMM: 0, SemiAperMM: 10},
	}}
	am, err := air.Thermal(dT, &ThermalOptions{MountCTE: 1e-5})
	if err != nil {
		t.Fatal(err)
	}
	checkEq("air-only surface, mount CTE", am.Cards[2].RadiusMM, -60*(1+1e-5*dT))
	checkEq("air gap, mount CTE (short)", am.Cards[1].ThicknessMM, 20*(1+1e-5*dT))

	sys, err := hot.System()
	if err != nil {
		t.Fatal(err)
	}
	if !(sys.Mediums[1].N(geo.Wavelengthd) > geo.N_BK7.N(geo.Wavelengthd)) {
		t.Errorf("BK7 index should rise with ΔT=+25 K (relative to air)")
	}
	if _, ok := sys.Mediums[1].(geo.ThermalMedium); !ok {
		t.Errorf("glass medium is %T, want geo.ThermalMedium", sys.Mediums[1])
	}
}

// TestThermalRejectsBadInput: non-finite inputs are errors, not silent results.
func TestThermalRejectsBadInput(t *testing.T) {
	p := thinBiconvex("N-BK7")
	if _, err := p.Thermal(math.NaN(), nil); err == nil {
		t.Error("NaN ΔT accepted")
	}
	if _, err := p.Thermal(math.Inf(1), nil); err == nil {
		t.Error("infinite ΔT accepted")
	}
	if _, err := p.Thermal(10, &ThermalOptions{MountCTE: math.NaN()}); err == nil {
		t.Error("NaN mount CTE accepted")
	}
	if _, err := p.ThermalDriftOf(geo.Wavelengthd, 0, nil); err == nil {
		t.Error("zero ΔT drift accepted")
	}
}

// TestThermalDriftIntegratesAnalytic: the EFL drift returned by ThermalDriftOf
// over [0, ΔT] equals the integral of the analytic thin-lens derivative over
// the same interval (Simpson's rule on 17 points; the integrand is smooth, so
// the quadrature error is far below the 1e-6 bound).
func TestThermalDriftIntegratesAnalytic(t *testing.T) {
	const lam = geo.Wavelengthd
	for _, glass := range []string{"N-BK7", "fused silica"} {
		p := thinBiconvex(glass)
		base, _ := geo.Glass(glass)
		g, _ := geo.ThermalGlassByName(glass)
		const dT = 30.0
		drift, err := p.ThermalDriftOf(lam, dT, nil)
		if err != nil {
			t.Fatal(err)
		}
		deriv := func(x float64) float64 {
			q, err := p.Thermal(x, nil)
			if err != nil {
				t.Fatal(err)
			}
			par, err := q.Paraxial(lam)
			if err != nil {
				t.Fatal(err)
			}
			n0 := base.N(lam)
			return geo.ThinLensEFLDrift(par.EFL, n0+g.DeltaNRel(lam, x, n0), g.DnRelDT(lam, x, n0), g.CTE, x)
		}
		const nSeg = 16
		h := dT / nSeg
		sum := deriv(0) + deriv(dT)
		for i := 1; i < nSeg; i++ {
			w := 2.0
			if i%2 == 1 {
				w = 4.0
			}
			sum += w * deriv(float64(i)*h)
		}
		integral := sum * h / 3
		if rel := math.Abs(drift.EFL1-drift.EFL0-integral) / math.Abs(integral); rel > 1e-6 {
			t.Errorf("%s: ΔEFL over %g K = %.10g mm, ∫ analytic derivative = %.10g mm (rel %.2g)",
				glass, dT, drift.EFL1-drift.EFL0, integral, rel)
		}
		t.Logf("%s thin singlet f=%.4f mm: EFL drift %+.6f µm/K (secant over %g K), BFL %+.6f µm/K",
			glass, drift.EFL0, drift.EFLPerK*1e3, dT, drift.BFLPerK*1e3)
	}
}

// TestThermalComposesWithMeltHook: the thermal shift is applied on top of the
// glass-melt index offset installed by the tolerance analysis (both use the
// mediumWrap hook): N_total(ΔT) = N_melt + Δn_rel(λ, ΔT, N_melt).
func TestThermalComposesWithMeltHook(t *testing.T) {
	const lam = geo.Wavelengthd
	const dT = 30.0
	p := thinBiconvex("N-BK7")
	// A melt-only prescription: the tolerance hook with a fixed Δn = +0.0005.
	melt, err := perturbPrescription(p, ToleranceSpec{IndexDelta: 0.0005}, rand.New(rand.NewSource(7)))
	if err != nil {
		t.Fatal(err)
	}
	sysMelt, err := melt.System()
	if err != nil {
		t.Fatal(err)
	}
	both, err := melt.Thermal(dT, nil)
	if err != nil {
		t.Fatal(err)
	}
	sysBoth, err := both.System()
	if err != nil {
		t.Fatal(err)
	}
	g, _ := geo.ThermalGlassByName("N-BK7")
	nMelt := sysMelt.Mediums[1].N(lam)
	want := nMelt + g.DeltaNRel(lam, dT, nMelt)
	if got := sysBoth.Mediums[1].N(lam); math.Abs(got-want) > 1e-14 {
		t.Errorf("composed index %.15g, want %.15g", got, want)
	}
}
