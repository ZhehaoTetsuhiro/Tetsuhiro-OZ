package geo

import (
	"fmt"
	"math"
)

// Spectral line wavelengths in air (nm), from the standard Fraunhofer set. They
// are the reference wavelengths for refractive index data and colour
// aberration analysis (F and C define the usual Achromat doublet condition,
// d/e are the design wavelengths).
const (
	WavelengthUV   = 365.0146 // i-line (mercury)
	Wavelengthg    = 435.8343 // g-line (mercury)
	WavelengthF    = 486.1327 // F-line (hydrogen, blue)
	Wavelengthe    = 546.0740 // e-line (mercury)
	Wavelengthd    = 587.5618 // d-line (helium, yellow)
	WavelengthC    = 656.2725 // C-line (hydrogen, red)
	WavelengthIR   = 1013.98  // t-line (mercury)
	WavelengthHeNe = 632.8    // helium-neon laser
	WavelengthNdY  = 1064.0   // Nd:YAG laser
)

// Medium is a homogeneous, isotropic, non-absorbing optical medium. N returns
// the real refractive index at the given vacuum wavelength in nm; Name returns
// a short human readable label used in reports and error messages.
//
// Implementations must be safe for concurrent use and must not mutate their
// receivers.
type Medium interface {
	N(wavelengthNm float64) float64
	Name() string
}

// Air is vacuum/air: n = 1 at every wavelength. This is the default medium of
// a System and the zero-cost option for the spaces between lenses.
type Air struct{}

// N returns 1.
func (Air) N(wavelengthNm float64) float64 { return 1 }

// Name implements Medium.
func (Air) Name() string { return "air" }

// Constant is a non-dispersive medium of a fixed refractive index. The zero
// value is air (n = 1); if Label is empty a default name is generated.
type Constant struct {
	Index float64 // refractive index (dimensionless)
	Label string  // optional name
}

// ConstantIndex returns a non-dispersive medium with index n.
func ConstantIndex(n float64) Constant { return Constant{Index: n} }

// N returns the constant index.
func (m Constant) N(wavelengthNm float64) float64 { return m.Index }

// Name implements Medium.
func (m Constant) Name() string {
	if m.Label != "" {
		return m.Label
	}
	return fmt.Sprintf("n=%.5g", m.Index)
}

// Sellmeier is a medium described by the three-term Sellmeier dispersion
// equation (λ in µm, the returned index in the visible is the real refractive
// index of the glass):
//
//	n²(λ) = 1 + Σ_i B_i λ² / (λ² - C_i) ,  i = 1..3
//
// B is the dimensionless numerator and C the squared resonance wavelengths in
// µm² (C_i = λ_i²). The formula reproduces the catalogue index of Schott and
// fused-silica glasses to better than 1e-5 in the visible and near IR. The zero
// value is vacuum.
type Sellmeier struct {
	Label string
	B     [3]float64
	C     [3]float64
}

// N evaluates the Sellmeier equation at wavelengthNm (nm).
func (m Sellmeier) N(wavelengthNm float64) float64 {
	l := wavelengthNm / 1000 // µm
	l2 := l * l
	sum := 1.0
	for i := 0; i < 3; i++ {
		sum += m.B[i] * l2 / (l2 - m.C[i])
	}
	if sum < 0 {
		return math.NaN()
	}
	return math.Sqrt(sum)
}

// Name implements Medium.
func (m Sellmeier) Name() string {
	if m.Label == "" {
		return "sellmeier"
	}
	return m.Label
}

// Nd returns the refractive index at the helium d-line (587.5618 nm).
func (m Sellmeier) Nd() float64 { return m.N(Wavelengthd) }

// AbbeNumber returns the Abbe number Vd = (nd - 1) / (nF - nC) computed from
// this dispersion model. The result is +Inf when the model is dispersionless.
func (m Sellmeier) AbbeNumber() float64 {
	nd := m.N(Wavelengthd)
	den := m.N(WavelengthF) - m.N(WavelengthC)
	if den == 0 {
		return math.Inf(1)
	}
	return (nd - 1) / den
}

// Sellmeier coefficient sets. Sources: Schott Optical Glass Data Sheets
// (N-BK7, N-SF5, N-SK16, F2) and I. H. Malitson, J. Opt. Soc. Am. 55, 1205
// (1965) for fused silica. C values are the squares of the resonance
// wavelengths in µm (e.g. 0.0684043² = 0.004679148 for fused silica).
var (
	// N_BK7 is Schott N-BK7 borosilicate crown glass (nd = 1.51680, Vd = 64.17).
	N_BK7 = Sellmeier{
		Label: "N-BK7",
		B:     [3]float64{1.03961212, 0.231792344, 1.01046945},
		C:     [3]float64{0.00600069867, 0.0200179144, 103.560653},
	}

	// FusedSilica is synthetic fused silica / SiO₂ (Malitson 1965,
	// nd = 1.45846, Vd = 67.8), valid from the deep UV to ~3.5 µm.
	FusedSilica = Sellmeier{
		Label: "fused silica",
		B:     [3]float64{0.6961663, 0.4079426, 0.8974794},
		C:     [3]float64{0.004679148, 0.013512064, 97.934003},
	}

	// N_SF5 is Schott N-SF5 dense flint glass (nd = 1.67271, Vd = 32.25).
	N_SF5 = Sellmeier{
		Label: "N-SF5",
		B:     [3]float64{1.52481889, 0.187085527, 1.42729015},
		C:     [3]float64{0.011254756, 0.0588995392, 129.141675},
	}

	// N_SK16 is Schott N-SK16 lanthanum crown glass (nd = 1.62041, Vd = 60.32).
	N_SK16 = Sellmeier{
		Label: "N-SK16",
		B:     [3]float64{1.34317774, 0.241144399, 0.994317969},
		C:     [3]float64{0.00704687339, 0.0229005, 92.7508526},
	}

	// F2 is Schott F2 flint glass (nd = 1.62004, Vd = 36.37).
	F2 = Sellmeier{
		Label: "F2",
		B:     [3]float64{1.34533359, 0.209073176, 0.937357162},
		C:     [3]float64{0.00997743871, 0.0470450767, 111.886764},
	}
)

// Glass returns the Sellmeier model of a catalogue glass by name. Recognised
// names (case-sensitive): "N-BK7", "BK7", "fused silica", "silica", "N-SF5",
// "SF5", "N-SK16", "SK16", "F2".
func Glass(name string) (Sellmeier, bool) {
	switch name {
	case "N-BK7", "BK7", "nbk7", "bk7":
		return N_BK7, true
	case "fused silica", "FusedSilica", "silica", "SiO2":
		return FusedSilica, true
	case "N-SF5", "SF5", "nsf5", "sf5":
		return N_SF5, true
	case "N-SK16", "SK16", "nsk16", "sk16":
		return N_SK16, true
	case "F2", "f2":
		return F2, true
	}
	return Sellmeier{}, false
}

// GlassNames lists the catalogue glasses known to Glass, in a stable order.
func GlassNames() []string {
	return []string{"N-BK7", "fused silica", "N-SF5", "N-SK16", "F2"}
}
