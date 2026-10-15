package geo

import "fmt"

// Thermo-optic and expansion data of the catalogue glasses (G9).
//
// This file is new data; the Sellmeier catalogue in medium.go is untouched.
// Every constant below was read from the cited document (URL and version in
// Source). An entry whose constants could not all be checked against a cited
// source carries Verified = false and a Note starting with 待核对; such
// entries are usable but must be reported as unverified.
//
// Sources (all retrieved 2026-10-15 for this revision):
//
//	[S1] SCHOTT TIE-19 "Temperature Coefficient of the Refractive Index",
//	     version August 2023 — formula (2)/(3) and Tab. 2 (N-BK7, F2 e-line
//	     coefficients used as cross-checks).
//	     https://media.schott.com/api/public/content/63b62a0d73a640efa8c6fde92f2ad5c2?v=493c707e&download=true
//	[S2] SCHOTT datasheet N-BK7 (517642.251, "As of 01-Dec-2023"): dn/dT
//	     constants, Δn_rel/ΔT and Δn_abs/ΔT table, α(−30/+70) and α(+20/+300).
//	     https://media.schott.com/api/public/content/41e799d0bf874807a0bb8e702fbb75b5?v=54856406
//	[S3] SCHOTT datasheet N-SF5 (673323.286): as [S2].
//	     https://media.schott.com/api/public/content/f8d4ac6649594677ab1b17b05ebd6a8f?v=c75caf58
//	[S4] SCHOTT datasheet N-SK16 (620603.358): as [S2].
//	     https://media.schott.com/api/public/content/5e904e6952b54f05be70518c4b669395?v=47550206
//	[S5] SCHOTT datasheet F2 ("Stand 01.02.2014", German original): constants
//	     and α values (the F2 nd = 1.62004, vd = 36.37 match Glass("F2")).
//	     https://media.schott.com/api/public/content/d7e653d4b846492b83434c7d4b3ddcfd?v=40f5a5d0
//	[S6] Heraeus Quarzglas, "Fused Silica for Optical Applications — Data and
//	     Properties" (2026-09 brochure): dn/dT constants of Infrasil 301,
//	     Suprasil 311, Suprasil 3001 and the relative temperature coefficient
//	     table (ppm/K, to ambient air).
//	     https://www.heraeus-covantics.com/media/Media/Documents/Products_and_Solutions/OPT/EN/Heraeus_brochure_Silica_for_Optical_Applications_2026-09_EN_WEB.pdf
//	[S7] Heraeus Quarzglas, "Quartz Glass for Optics — Data and Properties"
//	     (mirror copy): mean thermal expansion coefficient table.
//	     https://lightmachinery.com/media/1560/heraeus-fused-silica.pdf
//
// Model conventions (see thermo_model.go): ΔT = T − 20 °C; the index values
// of Glass() are the catalogue indices relative to air, so the thermal
// correction is applied to the relative index (lens in air).

// ThermalGlass is the thermo-optic and expansion description of one glass.
//
// Index: the TIE-19 form of the absolute index temperature coefficient
//
//	dn_abs/dT(λ, T) = (n²(λ,T0) − 1) / (2 n(λ,T0)) ·
//	                  [ D0 + 2 D1 ΔT + 3 D2 ΔT² + (E0 + 2 E1 ΔT) / (λ² − λTK²) ]
//
// with λ in µm, ΔT in K, T0 = 20 °C. The integral Δn_abs(λ, ΔT) is the
// TIE-19 eq. (3). Validity stated by the source: −100 °C … +140 °C and
// 0.365 µm … 1.014 µm.
//
// Expansion: CTE is the linear expansion coefficient (1/K) used for every
// dimension (radii, thicknesses, apertures) of an element made of this glass.
// The value is the mean α(−30 °C … +70 °C) of the datasheet, the mean that is
// centred on 20 °C; the α(+20 … +300 °C) value is a high-temperature mean and
// is not used (it would bias the 20 °C linearisation by about 1e-6/K).
type ThermalGlass struct {
	Name     string  // canonical catalogue name, equal to the Sellmeier Label
	D0       float64 // 1/K
	D1       float64 // 1/K²
	D2       float64 // 1/K³
	E0       float64 // µm²/K
	E1       float64 // µm²/K²
	LambdaTK float64 // µm
	CTE      float64 // 1/K, mean α(−30/+70 °C)
	Source   string  // document(s) and table(s) the constants were read from
	Verified bool    // true: all constants read from a cited source
	Note     string  // caveats; starts with 待核对 when Verified is false
}

// thermalGlasses is the table, in a fixed order (no map iteration).
var thermalGlasses = []ThermalGlass{
	{
		Name: "N-BK7",
		D0:   1.86e-6, D1: 1.31e-8, D2: -1.37e-11,
		E0: 4.34e-7, E1: 6.27e-10, LambdaTK: 0.170,
		CTE:      7.1e-6,
		Source:   "[S2] N-BK7 datasheet: constants of formula for dn/dT; α−30/+70 = 7.1e-6 /K (the α+20/+300 = 8.3e-6 /K is not used); cross-checked against [S1] Tab. 2 (e-line, +20/+40 °C: 1.6e-6 abs / 3.0e-6 rel per K ≈ table value).",
		Verified: true,
	},
	{
		Name: "N-SF5",
		D0:   -2.51e-7, D1: 1.07e-8, D2: -2.40e-11,
		E0: 7.85e-7, E1: 1.15e-9, LambdaTK: 0.278,
		CTE:      7.9e-6,
		Source:   "[S3] N-SF5 datasheet (ID 673323.286): constants of formula for dn/dT; α−30/+70 = 7.9e-6 /K. Note: the TIE-19 Tab. 1 row 'SF5' is a different glass and is not used.",
		Verified: true,
	},
	{
		Name: "N-SK16",
		D0:   -2.37e-8, D1: 1.32e-8, D2: -1.29e-11,
		E0: 4.09e-7, E1: 5.17e-10, LambdaTK: 0.170,
		CTE:      6.3e-6,
		Source:   "[S4] N-SK16 datasheet (ID 620603.358): constants of formula for dn/dT; α−30/+70 = 6.3e-6 /K.",
		Verified: true,
	},
	{
		Name: "F2",
		D0:   1.51e-6, D1: 1.56e-8, D2: -2.78e-11,
		E0: 9.34e-7, E1: 1.04e-9, LambdaTK: 0.250,
		CTE:      8.2e-6,
		Source:   "[S5] F2 datasheet (Stand 01.02.2014): constants of formula for dn/dT; α−30/+70 = 8.2e-6 /K. Cross-checked against [S1] Tab. 2 (F2 e-line +20/+40 °C relative 4.40e-6 /K).",
		Verified: true,
		Note:     "Constants read from the 2014 issue. The 2023 TIE-19 e-line value (4.40e-6 /K relative) agrees within rounding; the current F2 datasheet issue was not retrieved (待核对).",
	},
	{
		Name: "fused silica",
		D0:   2.18e-5, D1: 2.45e-8, D2: -2.72e-11,
		E0: 2.31e-7, E1: 2.21e-10, LambdaTK: 0.235,
		CTE:      5.1e-7,
		Source:   "index: [S6] Suprasil 311 constants (synthetic fused silica); CTE: [S7] mean 0…100 °C = 5.1e-7 /K.",
		Verified: false,
		Note:     "待核对: the index constants are the Heraeus Suprasil 311 sample, not the Malitson sample used by the Sellmeier entry 'fused silica'; the CTE is a 0…100 °C mean of an unspecified grade, not the 20 °C value. The mapping of the Heraeus constant columns to the relative-coefficient columns is not stated in the excerpt.",
	},
}

// ThermalGlassNames lists the glasses that have thermo-optic data, in table
// order.
func ThermalGlassNames() []string {
	out := make([]string, len(thermalGlasses))
	for i, g := range thermalGlasses {
		out[i] = g.Name
	}
	return out
}

// ThermalGlassByName returns the thermo-optic data of a catalogue glass. The
// name may be any alias accepted by Glass (e.g. "BK7", "silica"); the lookup
// goes through the Sellmeier catalogue so the two tables cannot disagree.
func ThermalGlassByName(name string) (ThermalGlass, bool) {
	s, ok := Glass(name)
	if !ok {
		return ThermalGlass{}, false
	}
	for _, g := range thermalGlasses {
		if g.Name == s.Label {
			return g, true
		}
	}
	return ThermalGlass{}, false
}

// String implements fmt.Stringer.
func (g ThermalGlass) String() string {
	state := "verified"
	if !g.Verified {
		state = "待核对"
	}
	return fmt.Sprintf("thermal[%s α=%.3g/K %s]", g.Name, g.CTE, state)
}
