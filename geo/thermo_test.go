package geo

import (
	"math"
	"testing"
)

// Datasheet values used as external cross-checks. Each is the coefficient
// averaged over the temperature interval +20…+40 °C (ΔT = 0…20 K), in 1e-6/K,
// as printed in the cited datasheet (one decimal; F2 from [S1] Tab. 2, two
// decimals). Wavelengths: 1060.0 nm, e-line 546.1 nm, g-line 435.8 nm. A value
// of -1 means "not used" (not printed in the excerpt that was read).
type dsRow struct {
	name       string
	relNIR     float64 // Δn_rel/ΔT at 1060.0 nm
	relE, relG float64
	absNIR     float64 // Δn_abs/ΔT at 1060.0 nm (-1 = not given)
	absE, absG float64
	source     string
}

var datasheetRows = []dsRow{
	{"N-BK7", 2.4, 3.0, 3.5, 1.1, 1.6, 2.1, "[S2] N-BK7 datasheet, Color Code table, +20/+40"},
	{"N-SF5", 1.8, 3.4, 5.5, 0.4, 2.0, 4.0, "[S3] N-SF5 datasheet, Color Code table, +20/+40"},
	{"N-SK16", 1.7, 2.3, 2.9, 0.3, 0.9, 1.4, "[S4] N-SK16 datasheet, Color Code table, +20/+40"},
	{"F2", -1, 4.40, -1, -1, -1, -1, "[S1] TIE-19 Tab. 2, F2 e-line +20/+40 (relative); other entries not printed in the excerpt (-1 = not used)"},
}

// TestThermalGlassTable checks the table itself: every entry resolves through
// the Sellmeier catalogue, carries a source, and the unverified entries are
// labelled 待核对.
func TestThermalGlassTable(t *testing.T) {
	names := ThermalGlassNames()
	if len(names) != len(GlassNames()) {
		t.Fatalf("thermal table has %d glasses, catalogue has %d", len(names), len(GlassNames()))
	}
	for _, n := range GlassNames() {
		g, ok := ThermalGlassByName(n)
		if !ok {
			t.Errorf("no thermo-optic data for catalogue glass %q", n)
			continue
		}
		if g.Source == "" {
			t.Errorf("%s: empty source", n)
		}
		if !(g.CTE > 0) || math.IsNaN(g.D0+g.D1+g.D2+g.E0+g.E1+g.LambdaTK) {
			t.Errorf("%s: non-physical constants %+v", n, g)
		}
		if !g.Verified && len(g.Note) < len("待核对") {
			t.Errorf("%s: unverified entry without a 待核对 note", n)
		}
	}
	if g, ok := ThermalGlassByName("BK7"); !ok || g.Name != "N-BK7" {
		t.Errorf("alias BK7 resolves to %+v, want N-BK7", g)
	}
	if _, ok := ThermalGlassByName("no-such-glass"); ok {
		t.Error("unknown glass must not resolve")
	}
	if g, _ := ThermalGlassByName("fused silica"); g.Verified {
		t.Error("fused silica is expected to stay 待核对 (grade and sample mismatch)")
	}
}

// TestDatasheetCoefficients compares the TIE-19 formula with the datasheet
// values at 1060 nm, e and g. The tolerance is the datasheet's own rounding
// (±0.05e-6/K for one decimal); F2 is given to two decimals in TIE-19.
func TestDatasheetCoefficients(t *testing.T) {
	const tol = 0.06e-6 // datasheet rounding ±0.05e-6/K plus margin
	lines := []float64{1060.0, Wavelengthe, Wavelengthg}
	worstRel, worstAbs := 0.0, 0.0
	defer func() {
		t.Logf("datasheet cross-check: max |model − datasheet| rel %.3g, abs %.3g (×1e-6 /K; datasheet rounding ±0.05)", worstRel*1e6, worstAbs*1e6)
	}()
	for _, row := range datasheetRows {
		g, ok := ThermalGlassByName(row.name)
		if !ok {
			t.Fatalf("%s: no thermal data", row.name)
		}
		base, _ := Glass(row.name)
		for k, lam := range lines {
			n0 := base.N(lam)
			// Average over ΔT = 0…20 K: (Δn(20 K) − Δn(0)) / 20 K.
			relWant := []float64{row.relNIR, row.relE, row.relG}[k] * 1e-6
			relGot := g.DeltaNRel(lam, 20, n0) / 20
			if relWant > -1e-9 {
				worstRel = math.Max(worstRel, math.Abs(relGot-relWant))
			}
			if relWant > -1e-9 && math.Abs(relGot-relWant) > tol {
				t.Errorf("%s λ=%.1f: Δn_rel/ΔT(+20/+40) = %.4g, datasheet %.4g (%s)",
					row.name, lam, relGot*1e6, relWant*1e6, row.source)
			}
			absWant := []float64{row.absNIR, row.absE, row.absG}[k] * 1e-6
			absGot := g.deltaNAbs(lam, 20, n0) / 20
			if absWant > -1e-9 {
				worstAbs = math.Max(worstAbs, math.Abs(absGot-absWant))
			}
			if absWant > -1e-9 && math.Abs(absGot-absWant) > tol {
				t.Errorf("%s λ=%.1f: Δn_abs/ΔT(+20/+40) = %.4g, datasheet %.4g (%s)",
					row.name, lam, absGot*1e6, absWant*1e6, row.source)
			}
		}
	}
}

// TestHeraeusFusedSilicaRelative compares the relative coefficient of the
// synthetic fused silica constants with the Heraeus table [S6]: 10.3 ppm/K at
// 546.1 nm for +20…+40 °C. The table has one decimal (±0.05 ppm/K); the
// tolerance adds the same margin as the datasheet test.
func TestHeraeusFusedSilicaRelative(t *testing.T) {
	g, ok := ThermalGlassByName("fused silica")
	if !ok {
		t.Fatal("fused silica: no thermal data")
	}
	base, _ := Glass("fused silica")
	got := g.DeltaNRel(Wavelengthe, 20, base.N(Wavelengthe)) / 20
	const want = 10.3e-6
	if math.Abs(got-want) > 0.1e-6 {
		t.Errorf("fused silica Δn_rel/ΔT(546.1 nm, +20/+40) = %.3f ppm/K, Heraeus table 10.3", got*1e6)
	}
}

// TestDnRelDTFiniteDifference checks the analytic derivative DnRelDT against a
// central difference of DeltaNRel in ΔT, for several glasses, wavelengths and
// temperature changes. Tolerance 1e-6 relative (the bound stated for the
// thermal model); the step is 0.01 K.
func TestDnRelDTFiniteDifference(t *testing.T) {
	const h = 0.01
	worst := 0.0
	defer func() { t.Logf("max relative |FD − analytic| of dn_rel/dT = %.3g (bound 1e-6)", worst) }()
	for _, name := range ThermalGlassNames() {
		g, _ := ThermalGlassByName(name)
		base, _ := Glass(name)
		for _, lam := range []float64{365.0146, Wavelengthe, Wavelengthd, 1060} {
			n0 := base.N(lam)
			for _, dT := range []float64{-60, -10, 0, 25, 100} {
				fd := (g.DeltaNRel(lam, dT+h, n0) - g.DeltaNRel(lam, dT-h, n0)) / (2 * h)
				an := g.DnRelDT(lam, dT, n0)
				rel := math.Abs(fd-an) / math.Abs(an)
				worst = math.Max(worst, rel)
				if rel > 1e-6 {
					t.Errorf("%s λ=%.1f ΔT=%g: FD %.10g vs analytic %.10g (rel %.2g)",
						name, lam, dT, fd, an, rel)
				}
			}
		}
	}
}

// TestThermalMediumIdentityAtZero: a ThermalMedium at ΔT = 0 returns the
// catalogue index bit for bit, and a non-zero ΔT changes it.
func TestThermalMediumIdentityAtZero(t *testing.T) {
	g, _ := ThermalGlassByName("N-BK7")
	for _, lam := range []float64{365, 486.1, 546.074, 587.56, 656.27, 1060} {
		m := ThermalMedium{Base: N_BK7, Glass: g, DT: 0}
		if got, want := m.N(lam), N_BK7.N(lam); got != want {
			t.Errorf("λ=%g: ΔT=0 index %.17g != base %.17g", lam, got, want)
		}
		hot := ThermalMedium{Base: N_BK7, Glass: g, DT: 30}
		if hot.N(lam) <= N_BK7.N(lam) {
			t.Errorf("λ=%g: BK7 index should rise with ΔT>0 (relative to air)", lam)
		}
	}
}

// TestScaleProfileIdentityAndScaling: s = 1 is the identity, and the scaled
// sag obeys z'(r) = s·z(r/s) for sphere, conic and even asphere.
func TestScaleProfileIdentityAndScaling(t *testing.T) {
	profiles := []Profile{
		Sphere{R: 80},
		Conic{R: -55, K: -0.7},
		EvenAsphere{R: 40, K: -1.2, Coeffs: []float64{3e-6, -2e-9, 1e-12}},
		Plane{},
	}
	for _, p := range profiles {
		id, err := ScaleProfile(p, 1)
		if err != nil {
			t.Fatalf("%v: %v", p, err)
		}
		for _, r := range []float64{0, 1, 5, 12.5} {
			if a, b := id.Sag(r), p.Sag(r); a != b && !(math.IsNaN(a) && math.IsNaN(b)) {
				t.Errorf("%v: s=1 changes the sag at r=%g", p, r)
			}
		}
		const s = 1.000713
		sc, err := ScaleProfile(p, s)
		if err != nil {
			t.Fatal(err)
		}
		for _, r := range []float64{0.5, 3, 9} {
			want := s * p.Sag(r/s)
			checkRel(t, "scaled sag "+p.String(), sc.Sag(r), want, 1e-12)
		}
	}
	if _, err := ScaleProfile(Sphere{R: 10}, 0); err == nil {
		t.Error("ScaleProfile must reject a zero factor")
	}
}
