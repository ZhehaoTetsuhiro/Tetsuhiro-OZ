package coating

import (
	"fmt"
	"math"
	"sort"
	"strings"
)

// Dispersion maps a wavelength (nm) to the complex refractive index of a
// medium, using the convention n = n' + i·k with k >= 0 absorbing. All models
// are pure functions of the wavelength and are safe for concurrent use; a
// wavelength of zero or less, or a wavelength inside a model's pole, yields a
// non-finite index, which the optics entry points reject with an error.
type Dispersion interface {
	// IndexNm returns the complex refractive index at lambdaNm nanometres.
	IndexNm(lambdaNm float64) complex128
	// Lossless reports whether the extinction coefficient is identically zero,
	// so callers can expect R + T = 1 exactly.
	Lossless() bool
	// Model returns a short human-readable description of the model and its
	// validity range.
	Model() string
}

// ConstIndex is a non-dispersive medium, n = N + i·K. It is used for air and
// for synthetic "ideal glass" layers in analytic designs.
type ConstIndex struct {
	N float64 // real refractive index
	K float64 // extinction coefficient, >= 0 for an absorbing medium
}

// IndexNm implements Dispersion.
func (c ConstIndex) IndexNm(float64) complex128 { return complex(c.N, c.K) }

// Lossless implements Dispersion.
func (c ConstIndex) Lossless() bool { return c.K == 0 }

// Model implements Dispersion.
func (c ConstIndex) Model() string {
	if c.K == 0 {
		return fmt.Sprintf("constant n = %g (lossless)", c.N)
	}
	return fmt.Sprintf("constant n = %g + %gi", c.N, c.K)
}

// Cauchy is the Cauchy series in which the wavelength is expressed in
// micrometres,
//
//	n(lambda) = A + B/lambda^2 + C/lambda^4,
//	k(lambda) = Ka + Kb/lambda^2,
//
// the usual analytic fit for evaporated thin films over the visible and near
// infrared. Set Ka = Kb = 0 for a lossless model.
type Cauchy struct {
	A, B, C  float64 // refractive index coefficients (lambda in micrometres)
	Ka, Kb   float64 // extinction coefficients (lambda in micrometres)
}

// IndexNm implements Dispersion.
func (c Cauchy) IndexNm(lambdaNm float64) complex128 {
	if !(lambdaNm > 0) {
		return complex(math.NaN(), math.NaN())
	}
	l2 := lambdaNm * lambdaNm * 1e-6 // (nm -> um)^2
	if l2 == 0 {
		return complex(math.NaN(), math.NaN())
	}
	n := c.A + c.B/l2 + c.C/(l2*l2)
	k := c.Ka + c.Kb/l2
	return complex(n, k)
}

// Lossless implements Dispersion.
func (c Cauchy) Lossless() bool { return c.Ka == 0 && c.Kb == 0 }

// Model implements Dispersion.
func (c Cauchy) Model() string {
	s := fmt.Sprintf("Cauchy n = %g + %g/l^2 + %g/l^4 (l in um)", c.A, c.B, c.C)
	if !c.Lossless() {
		s += fmt.Sprintf(", k = %g + %g/l^2", c.Ka, c.Kb)
	}
	return s
}

// Sellmeier is the three-term Sellmeier model, with the wavelength in
// micrometres and the poles C_i in micrometres squared:
//
//	n(lambda)^2 = 1 + sum_i B_i·lambda^2 / (lambda^2 - C_i),
//	k(lambda)   = Ka + Kb/lambda^2.
//
// A zero denominator makes the index non-finite, which the optics entry points
// report as an error naming the offending material.
type Sellmeier struct {
	B  [3]float64 // oscillator strengths
	C  [3]float64 // resonance wavelengths squared, um^2
	Ka float64    // extinction offset
	Kb float64    // extinction 1/lambda^2 coefficient
}

// IndexNm implements Dispersion.
func (s Sellmeier) IndexNm(lambdaNm float64) complex128 {
	if !(lambdaNm > 0) {
		return complex(math.NaN(), math.NaN())
	}
	l2 := lambdaNm * lambdaNm * 1e-6
	n2 := 1.0
	for i := range s.B {
		d := l2 - s.C[i]
		if d == 0 {
			return complex(math.NaN(), math.NaN())
		}
		n2 += s.B[i] * l2 / d
	}
	if n2 < 0 {
		return complex(math.NaN(), math.NaN())
	}
	k := 0.0
	if s.Ka != 0 || s.Kb != 0 {
		k = s.Ka + s.Kb/l2
	}
	return complex(math.Sqrt(n2), k)
}

// Lossless implements Dispersion.
func (s Sellmeier) Lossless() bool { return s.Ka == 0 && s.Kb == 0 }

// Model implements Dispersion.
func (s Sellmeier) Model() string {
	str := fmt.Sprintf("Sellmeier n^2 = 1 + %g l^2/(l^2-%g) + %g l^2/(l^2-%g) + %g l^2/(l^2-%g) (l in um)",
		s.B[0], s.C[0], s.B[1], s.C[1], s.B[2], s.C[2])
	if !s.Lossless() {
		str += fmt.Sprintf(", k = %g + %g/l^2", s.Ka, s.Kb)
	}
	return str
}

// Tabulated is a measured n(lambda), k(lambda) table, linearly interpolated in
// wavelength (nm). The wavelengths must be strictly increasing. The table is
// clamped outside its range (the first/last value is repeated), so it should
// only be queried inside the range documented by the material using it; a
// table with fewer than two samples, or with non-increasing wavelengths, fails
// validation.
type Tabulated struct {
	Nm []float64 // wavelengths in nm, strictly increasing
	N  []float64 // real refractive index
	K  []float64 // extinction coefficient; nil means k = 0 everywhere
}

// IndexNm implements Dispersion.
func (t Tabulated) IndexNm(lambdaNm float64) complex128 {
	n := len(t.Nm)
	if n == 0 {
		return complex(math.NaN(), math.NaN())
	}
	if lambdaNm <= t.Nm[0] {
		return complex(t.N[0], t.k(0))
	}
	if lambdaNm >= t.Nm[n-1] {
		return complex(t.N[n-1], t.k(n-1))
	}
	i := sort.SearchFloat64s(t.Nm, lambdaNm)
	if i == 0 {
		return complex(t.N[0], t.k(0))
	}
	// t.Nm[i-1] < lambdaNm <= t.Nm[i]
	span := t.Nm[i] - t.Nm[i-1]
	w := 0.0
	if span > 0 {
		w = (lambdaNm - t.Nm[i-1]) / span
	}
	nv := t.N[i-1] + w*(t.N[i]-t.N[i-1])
	kv := t.k(i-1) + w*(t.k(i)-t.k(i-1))
	return complex(nv, kv)
}

// k returns the i-th extinction sample, treating a nil K slice as zero.
func (t Tabulated) k(i int) float64 {
	if len(t.K) == 0 {
		return 0
	}
	return t.K[i]
}

// Lossless implements Dispersion.
func (t Tabulated) Lossless() bool {
	for _, k := range t.K {
		if k != 0 {
			return false
		}
	}
	return true
}

// Model implements Dispersion.
func (t Tabulated) Model() string {
	if len(t.Nm) == 0 {
		return "empty table"
	}
	return fmt.Sprintf("tabulated n,k over %.4g-%.4g nm (%d samples, linear interpolation)",
		t.Nm[0], t.Nm[len(t.Nm)-1], len(t.Nm))
}

// Validate reports whether the table can be interpolated: at least two samples,
// strictly increasing wavelengths, and matched slice lengths.
func (t Tabulated) Validate() error {
	if len(t.Nm) < 2 {
		return fmt.Errorf("coating: tabulated dispersion needs at least 2 samples, got %d", len(t.Nm))
	}
	if len(t.N) != len(t.Nm) {
		return fmt.Errorf("coating: tabulated dispersion has %d wavelengths but %d index samples", len(t.Nm), len(t.N))
	}
	if len(t.K) != 0 && len(t.K) != len(t.Nm) {
		return fmt.Errorf("coating: tabulated dispersion has %d wavelengths but %d extinction samples", len(t.Nm), len(t.K))
	}
	for i := 1; i < len(t.Nm); i++ {
		if !(t.Nm[i] > t.Nm[i-1]) {
			return fmt.Errorf("coating: tabulated dispersion wavelengths must increase: %.6g at %d is not > %.6g", t.Nm[i], i, t.Nm[i-1])
		}
	}
	return nil
}

// Material is a named optical medium: a dispersion model plus the name used in
// reports. The zero Material has no dispersion model; the optics entry points
// reject it.
type Material struct {
	Name       string
	Dispersion Dispersion
}

// NewMaterial pairs a name with a dispersion model.
func NewMaterial(name string, d Dispersion) Material { return Material{Name: name, Dispersion: d} }

// IndexNm returns the complex refractive index at lambdaNm nanometres. An
// invalid (zero-value) Material returns NaN, which the optics entry points
// report as an error.
func (m Material) IndexNm(lambdaNm float64) complex128 {
	if m.Dispersion == nil {
		return complex(math.NaN(), math.NaN())
	}
	return m.Dispersion.IndexNm(lambdaNm)
}

// N returns the real refractive index n' at lambdaNm.
func (m Material) N(lambdaNm float64) float64 { return real(m.IndexNm(lambdaNm)) }

// K returns the extinction coefficient k at lambdaNm.
func (m Material) K(lambdaNm float64) float64 { return imag(m.IndexNm(lambdaNm)) }

// Lossless reports whether the medium has no absorption anywhere.
func (m Material) Lossless() bool { return m.Dispersion != nil && m.Dispersion.Lossless() }

// Valid reports whether the material carries a dispersion model and a name.
func (m Material) Valid() bool { return m.Dispersion != nil && m.Name != "" }

// String implements fmt.Stringer.
func (m Material) String() string {
	if m.Dispersion == nil {
		return m.Name + "(invalid)"
	}
	return m.Name
}

// siliconDispersion implements the Salzberg & Villa silicon model, whose
// constant and 1/lambda^2 terms sit outside the three-term Sellmeier form:
//
//	n^2 = 11.6858 + 0.939816/l^2 + 0.000993358/(l^2 - 1.22567^2),  l in um.
type siliconDispersion struct{}

func (siliconDispersion) IndexNm(lambdaNm float64) complex128 {
	if !(lambdaNm > 0) {
		return complex(math.NaN(), math.NaN())
	}
	l2 := lambdaNm * lambdaNm * 1e-6 // (nm -> um)^2
	d := l2 - 1.22567*1.22567
	if d == 0 || l2 == 0 {
		return complex(math.NaN(), math.NaN())
	}
	n2 := 11.6858 + 0.939816/l2 + 0.000993358/d
	if n2 < 0 {
		return complex(math.NaN(), math.NaN())
	}
	return complex(math.Sqrt(n2), 0)
}

func (siliconDispersion) Lossless() bool { return true }

func (siliconDispersion) Model() string {
	return "Sellmeier (Salzberg & Villa 1957) n^2 = 11.6858 + 0.939816/l^2 + 0.000993358/(l^2-1.22567^2), l in um"
}

// ---------------------------------------------------------------------------
// Built-in material library
//
// Every entry documents its source and validity range. Bulk single-crystal
// data is used where a reliable model exists (fluorides, silica, sapphire,
// semiconductors); the deposition-dependent oxides are typical evaporated /
// ion-assisted values fitted with the Cauchy series, because the index of a
// real film depends on the plant, the pressure and the substrate temperature.
// ---------------------------------------------------------------------------

var (
	// Air (vacuum), n = 1 exactly. Used as the default incidence medium.
	Air = Material{Name: "Air", Dispersion: ConstIndex{N: 1}}

	// MgF2: magnesium fluoride single crystal, ordinary ray. Sellmeier after
	// M. J. Dodge, "Refractive properties of magnesium fluoride", Appl. Opt. 23,
	// 1980-1985 (1984); valid 0.2-7.0 um, n(550 nm) = 1.3785. Evaporated MgF2
	// films are close to the crystal value (1.37-1.39) and slightly porous.
	MgF2 = Material{Name: "MgF2", Dispersion: Sellmeier{
		B: [3]float64{0.48755108, 0.39875031, 2.3120353},
		C: [3]float64{0.0018822238, 0.0089518950, 566.13567}, // (0.04338408, 0.09461442, 23.793604)^2
	}}

	// SiO2: fused silica. Sellmeier after I. H. Malitson, "Interspecimen
	// comparison of the refractive index of fused silica", J. Opt. Soc. Am. 55,
	// 1205-1209 (1965); valid 0.21-3.71 um, n(587.6 nm) = 1.4585, n(550 nm) =
	// 1.4601. Electron-beam evaporated SiO2 films run 0.5-1.5% lower
	// (1.44-1.47); use a Cauchy model for a specific plant.
	SiO2 = Material{Name: "SiO2", Dispersion: Sellmeier{
		B: [3]float64{0.6961663, 0.4079426, 0.8974794},
		C: [3]float64{0.00467914826, 0.0135120631, 97.9340025}, // (0.0684043, 0.1162414, 9.896161)^2
	}}

	// FusedSilica is the glass-substrate spelling of SiO2 (same model): the
	// standard UV/visible substrate, 0.21-3.7 um.
	FusedSilica = Material{Name: "FusedSilica", Dispersion: SiO2.Dispersion}

	// Sapphire: single-crystal Al2O3, ordinary ray, Sellmeier after I. H.
	// Malitson (as tabulated in the literature, e.g. the Sellmeier coefficient
	// table of refractiveindex.info, "sapphire, ordinary wave"); valid
	// 0.2-5 um, n(587.6 nm) = 1.7682, n(550 nm) = 1.7717. Sintered or
	// evaporated Al2O3 films are much lower - see Al2O3.
	Sapphire = Material{Name: "Sapphire", Dispersion: Sellmeier{
		B: [3]float64{1.43134930, 0.65054713, 5.3414021},
		C: [3]float64{0.00527992610, 0.0142382647, 325.017834},
	}}

	// Al2O3: alumina coating material, typical evaporated film. Cauchy fit to
	// n(550 nm) = 1.65; the Umicore thin-film catalogue gives 1.65-1.74 at
	// 550 nm depending on ion assistance and substrate temperature, with a
	// transparency range of roughly 200 nm - 8 um.
	Al2O3 = Material{Name: "Al2O3", Dispersion: Cauchy{A: 1.6000, B: 0.015125}}

	// Y2O3: yttrium oxide coating material, typical film. Cauchy fit to
	// n(550 nm) = 1.79 (Umicore range 1.76-1.95 at 550 nm, transparent from
	// about 250 nm to 12 um). Sometimes used as an intermediate-index layer.
	Y2O3 = Material{Name: "Y2O3", Dispersion: Cauchy{A: 1.7200, B: 0.021189}}

	// HfO2: hafnia coating material, typical ion-assisted film. Cauchy fit to
	// n(550 nm) = 2.00 (Umicore range 1.90-2.07 at 550 nm); transparent from
	// the near UV to about 8 um and the standard high-index material for
	// laser-damage-resistant coatings and 1064 nm mirrors.
	HfO2 = Material{Name: "HfO2", Dispersion: Cauchy{A: 1.9100, B: 0.027219}}

	// ZrO2: zirconia coating material, typical film. Cauchy fit to
	// n(550 nm) = 2.10; evaporated ZrO2 films range from about 2.0 (porous,
	// unheated) to 2.15 (dense), transparent from ~250 nm to ~8 um.
	ZrO2 = Material{Name: "ZrO2", Dispersion: Cauchy{A: 2.0000, B: 0.030250}}

	// Ta2O5: tantalum pentoxide coating material, typical ion-assisted film.
	// Cauchy fit to n(550 nm) = 2.10 (Umicore range 2.04-2.30 at 550 nm);
	// low absorption from ~350 nm to ~9 um, the standard high-index material
	// for near-infrared laser optics.
	Ta2O5 = Material{Name: "Ta2O5", Dispersion: Cauchy{A: 1.9900, B: 0.033247}}

	// TiO2: titanium dioxide coating material, typical electron-beam film.
	// Cauchy fit to n(550 nm) = 2.35 (Umicore range 2.07-2.50 at 550 nm);
	// the highest index available in the visible, transparent from ~400 nm
	// (fully oxidised) to ~7 um. Dense IAD films reach 2.4-2.5.
	TiO2 = Material{Name: "TiO2", Dispersion: Cauchy{A: 2.2000, B: 0.045374}}

	// Si: crystalline silicon. Sellmeier after C. D. Salzberg and J. J. Villa,
	// "Infrared refractive indexes of silicon, germanium and modified selenium
	// glass", J. Opt. Soc. Am. 47, 244 (1957); the model is lossless and valid
	// roughly 1.36-11 um, n(1.55 um) = 3.475 (measured 3.4777, Aspnes 1983).
	// Loss is deliberately not modelled: real silicon absorbs strongly below
	// ~1.1 um and weakly up to its IR phonon edge.
	Si = Material{Name: "Si", Dispersion: siliconDispersion{}}

	// Ge: germanium, mid-infrared. Tabulated n after the classic dispersion
	// tables (Tydex/HBOC germanium data, agreeing with Li 1980 and Palik);
	// transparent from ~2 to ~14 um, k is below 1e-4 across that band and is
	// modelled as zero. Below 2 um germanium is opaque - the absorption is NOT
	// modelled, so stay inside 2-14 um.
	Ge = Material{Name: "Ge", Dispersion: Tabulated{
		Nm: []float64{2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 10600, 11000, 12000, 13000, 14000},
		N:  []float64{4.1079, 4.0653, 4.0446, 4.0255, 4.0170, 4.0122, 4.0092, 4.0074, 4.0061, 4.0052, 4.0048, 4.0045, 4.0039, 4.0035, 4.0032},
	}}

	// Aluminum: metal mirror material. Tabulated n, k after E. D. Palik (ed.),
	// "Handbook of Optical Constants of Solids" (1985); sampled 400-800 nm,
	// linear interpolation in between. The table is clamped outside that range,
	// so use it in the visible only - the far-UV and infrared behaviour of
	// aluminium is very different.
	Aluminum = Material{Name: "Aluminum", Dispersion: Tabulated{
		Nm: []float64{400, 450, 500, 550, 600, 650, 700, 800},
		N:  []float64{0.490, 0.618, 0.855, 1.000, 1.130, 1.300, 1.497, 1.827},
		K:  []float64{4.860, 5.390, 6.040, 6.590, 7.150, 7.580, 7.910, 8.438},
	}}

	// BK7: SCHOTT N-BK7 borosilicate crown glass. Sellmeier from the SCHOTT
	// data sheet (nd = 1.5168, Abbe 64.17); valid 0.30-2.5 um, n(587.6 nm) =
	// 1.51673. The standard visible substrate.
	BK7 = Material{Name: "BK7", Dispersion: Sellmeier{
		B: [3]float64{1.03961212, 0.231792344, 1.01046945},
		C: [3]float64{0.00600069867, 0.0200179144, 103.560653},
	}}

	// SF5: SCHOTT N-SF5 dense flint glass. Sellmeier from the SCHOTT data
	// sheet / glassbank (nd = 1.67270); valid 0.365-2.33 um. High-index
	// substrate for wide-angle and immersed designs.
	SF5 = Material{Name: "SF5", Dispersion: Sellmeier{
		B: [3]float64{1.52481889, 0.187085527, 1.42729015},
		C: [3]float64{0.011254756, 0.0588995392, 129.141675},
	}}

	// SF11: SCHOTT N-SF11 dense flint glass. Sellmeier from the SCHOTT data
	// sheet (nd = 1.7847); valid 0.37-2.5 um. The densest common flint, used
	// where a very high substrate index is wanted.
	SF11 = Material{Name: "SF11", Dispersion: Sellmeier{
		B: [3]float64{1.73759695, 0.313747346, 1.89878101},
		C: [3]float64{0.013188707, 0.0623068142, 155.23629},
	}}
)

// Library returns the built-in material catalogue: coating materials first
// (MgF2, SiO2, Al2O3, Y2O3, HfO2, ZrO2, Ta2O5, TiO2, Si, Ge, Aluminum), then
// substrate glasses (BK7, FusedSilica, SF5, SF11, Sapphire) and Air.
func Library() []Material {
	return []Material{
		Air,
		MgF2, SiO2, Al2O3, Y2O3, HfO2, ZrO2, Ta2O5, TiO2, Si, Ge, Aluminum,
		BK7, FusedSilica, SF5, SF11, Sapphire,
	}
}

// materialAliases lets LookupMaterial accept the common alternative spellings.
var materialAliases = map[string]string{
	"fused silica": "FusedSilica",
	"fused-silica": "FusedSilica",
	"silica":       "SiO2",
	"alumina":      "Al2O3",
	"hafnia":       "HfO2",
	"zirconia":     "ZrO2",
	"tantala":      "Ta2O5",
	"titania":      "TiO2",
	"silicon":      "Si",
	"germanium":    "Ge",
	"aluminium":    "Aluminum",
	"al":           "Aluminum",
	"n-bk7":        "BK7",
	"n-sf5":        "SF5",
	"n-sf11":       "SF11",
	"sapphire":     "Sapphire",
	"vacuum":       "Air",
}

// LookupMaterial finds a built-in material by name. The comparison ignores
// case and a few common synonyms ("silica", "alumina", "n-bk7", ...).
func LookupMaterial(name string) (Material, bool) {
	key := strings.TrimSpace(name)
	if canonical, ok := materialAliases[strings.ToLower(key)]; ok {
		key = canonical
	}
	for _, m := range Library() {
		if strings.EqualFold(m.Name, key) {
			return m, true
		}
	}
	return Material{}, false
}
