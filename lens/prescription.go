package lens

import (
	"encoding/json"
	"fmt"
	"math"
	"os"

	"tetsuhiro-oz/geo"
)

// VarKind selects which quantity of a card a design variable refers to.
type VarKind int

// Variable kinds.
const (
	VarRadius VarKind = iota
	VarThickness
	VarConic
	VarCoeff // Coeffs[Order], see VarSpec.Order
)

// VarSpec points at one optimizable quantity inside a prescription: the
// radius, thickness, conic constant or one aspheric coefficient of card i.
// It is the bridge between the prescription text and the optimizer's flat
// variable vectors.
type VarSpec struct {
	Card  int
	Kind  VarKind
	Order int // Coeffs index for VarCoeff (0 = A4)
}

// String labels the variable for reports, e.g. "s2.R" or "s1.thk".
func (v VarSpec) String() string {
	switch v.Kind {
	case VarRadius:
		return fmt.Sprintf("s%d.R", v.Card)
	case VarThickness:
		return fmt.Sprintf("s%d.thk", v.Card)
	case VarConic:
		return fmt.Sprintf("s%d.K", v.Card)
	case VarCoeff:
		return fmt.Sprintf("s%d.A%d", v.Card, 2*v.Order+4)
	}
	return fmt.Sprintf("s%d.?", v.Card)
}

// Card is one surface of a sequential prescription.
//
// RadiusMM is the vertex radius of curvature in mm; RadiusMM == 0 or ±Inf
// means a flat surface. ConicK is the conic constant (0 = sphere, -1 =
// paraboloid, ...). Coeffs are the even aspheric coefficients A4, A6, A8...
// (empty for a pure conic/sphere). ThicknessMM is the axial distance to the
// next card (ignored on the last card). Glass names the medium that follows
// this surface; "" or "air" is air, otherwise a catalogue name accepted by
// geo.Glass (N-BK7, fused silica, N-SF5, N-SK16, F2). SemiAperMM > 0 sets a
// circular clear aperture of that radius; 0 means unlimited. Mirror makes
// the surface reflective (the medium does not change across it).
type Card struct {
	Name        string    `json:"name,omitempty"`
	RadiusMM    float64   `json:"radius_mm"`
	ConicK      float64   `json:"conic_k,omitempty"`
	Coeffs      []float64 `json:"coeffs,omitempty"`
	ThicknessMM float64   `json:"thickness_mm"`
	Glass       string    `json:"glass,omitempty"`
	SemiAperMM  float64   `json:"semi_aperture_mm,omitempty"`
	Mirror      bool      `json:"mirror,omitempty"`
}

// profile converts the card into a geo surface profile.
func (c Card) profile() (geo.Profile, error) {
	if math.IsInf(c.RadiusMM, 0) || c.RadiusMM == 0 {
		if len(c.Coeffs) > 0 {
			// Pure polynomial asphere on a flat base.
			return geo.EvenAsphere{R: math.Inf(1), K: 0, Coeffs: append([]float64(nil), c.Coeffs...)}, nil
		}
		return geo.Plane{}, nil
	}
	if len(c.Coeffs) > 0 {
		return geo.EvenAsphere{R: c.RadiusMM, K: c.ConicK, Coeffs: append([]float64(nil), c.Coeffs...)}, nil
	}
	if c.ConicK != 0 {
		return geo.Conic{R: c.RadiusMM, K: c.ConicK}, nil
	}
	return geo.Sphere{R: c.RadiusMM}, nil
}

// Prescription is a complete sequential lens description: an object
// position, an ordered list of cards and the stop. It converts to a
// geo.System for tracing and paraxial analysis, carries JSON I/O, and
// exposes its quantities as optimizer variables.
type Prescription struct {
	Name      string  `json:"name,omitempty"`
	ObjectZMM float64 `json:"object_z_mm"` // math.Inf(-1) for an object at infinity
	Cards     []Card  `json:"cards"`
	// Stop is the index of the aperture-stop card; the zero value (0) is
	// the common front-stop default. Negative or out-of-range values are
	// clamped by System().
	Stop int `json:"stop"`
}

// Air / glass helpers ------------------------------------------------------

// resolveGlass maps a card's glass name to a geo medium.
func resolveGlass(name string) (geo.Medium, error) {
	switch name {
	case "", "air", "Air", "AIR", "vacuum":
		return geo.Air{}, nil
	}
	if g, ok := geo.Glass(name); ok {
		return g, nil
	}
	return nil, fmt.Errorf("lens: unknown glass %q (known: %v)", name, geo.GlassNames())
}

// Validate checks structural consistency and returns a descriptive error.
func (p *Prescription) Validate() error {
	if len(p.Cards) == 0 {
		return fmt.Errorf("lens: prescription has no cards")
	}
	for i, c := range p.Cards {
		if math.IsNaN(c.RadiusMM) {
			return fmt.Errorf("lens: card %d radius is NaN", i)
		}
		if math.IsNaN(c.ThicknessMM) || math.IsInf(c.ThicknessMM, 0) {
			return fmt.Errorf("lens: card %d thickness is not finite (%v)", i, c.ThicknessMM)
		}
		if _, err := resolveGlass(c.Glass); err != nil {
			return fmt.Errorf("lens: card %d: %w", i, err)
		}
		if c.SemiAperMM < 0 {
			return fmt.Errorf("lens: card %d has a negative semi-aperture", i)
		}
		for j, a := range c.Coeffs {
			if math.IsNaN(a) || math.IsInf(a, 0) {
				return fmt.Errorf("lens: card %d coefficient A%d is not finite", i, 2*j+4)
			}
		}
	}
	return nil
}

// System converts the prescription into a geometry-ready geo.System. The
// object-side medium is air; every card's Glass fills the gap after it.
func (p *Prescription) System() (*geo.System, error) {
	if err := p.Validate(); err != nil {
		return nil, err
	}
	n := len(p.Cards)
	surfs := make([]geo.Surface, n)
	meds := make([]geo.Medium, n+1)
	meds[0] = geo.Air{}
	thick := make([]float64, 0, n-1)
	for i, c := range p.Cards {
		prof, err := c.profile()
		if err != nil {
			return nil, fmt.Errorf("lens: card %d: %w", i, err)
		}
		sf := geo.Surface{Profile: prof, Mirror: c.Mirror, Name: c.Name}
		if c.SemiAperMM > 0 {
			sf.Aperture = geo.CircleAperture(c.SemiAperMM)
		}
		surfs[i] = sf
		m, err := resolveGlass(c.Glass)
		if err != nil {
			return nil, fmt.Errorf("lens: card %d: %w", i, err)
		}
		meds[i+1] = m
		if i < n-1 {
			thick = append(thick, c.ThicknessMM)
		}
	}
	stop := p.Stop
	if stop < 0 {
		stop = 0
	}
	if stop > n-1 {
		stop = n - 1
	}
	sys := &geo.System{
		Surfaces:  surfs,
		Thickness: thick,
		Mediums:   meds,
		StopIndex: stop,
		Name:      p.Name,
	}
	if err := sys.Validate(); err != nil {
		return nil, err
	}
	return sys, nil
}

// Paraxial runs the geo first-order analysis of the prescription at
// wavelength lambdaNm (nm).
func (p *Prescription) Paraxial(lambdaNm float64) (*geo.Paraxial, error) {
	sys, err := p.System()
	if err != nil {
		return nil, err
	}
	return sys.Paraxial(lambdaNm, p.ObjectZMM)
}

// Variables ----------------------------------------------------------------

// GetVar reads one variable quantity.
func (p *Prescription) GetVar(s VarSpec) (float64, error) {
	if s.Card < 0 || s.Card >= len(p.Cards) {
		return 0, fmt.Errorf("lens: variable %s: card out of range", s)
	}
	c := &p.Cards[s.Card]
	switch s.Kind {
	case VarRadius:
		return c.RadiusMM, nil
	case VarThickness:
		return c.ThicknessMM, nil
	case VarConic:
		return c.ConicK, nil
	case VarCoeff:
		if s.Order < 0 || s.Order >= len(c.Coeffs) {
			return 0, fmt.Errorf("lens: variable %s: coefficient out of range", s)
		}
		return c.Coeffs[s.Order], nil
	}
	return 0, fmt.Errorf("lens: variable %s: unknown kind", s)
}

// SetVar writes one variable quantity.
func (p *Prescription) SetVar(s VarSpec, v float64) error {
	if s.Card < 0 || s.Card >= len(p.Cards) {
		return fmt.Errorf("lens: variable %s: card out of range", s)
	}
	c := &p.Cards[s.Card]
	switch s.Kind {
	case VarRadius:
		c.RadiusMM = v
	case VarThickness:
		c.ThicknessMM = v
	case VarConic:
		c.ConicK = v
	case VarCoeff:
		if s.Order < 0 || s.Order >= len(c.Coeffs) {
			return fmt.Errorf("lens: variable %s: coefficient out of range", s)
		}
		c.Coeffs[s.Order] = v
	default:
		return fmt.Errorf("lens: variable %s: unknown kind", s)
	}
	return nil
}

// Vars reads a list of quantities into a flat vector (optimizer convention:
// one slice, same order as the specs).
func (p *Prescription) Vars(specs []VarSpec) ([]float64, error) {
	out := make([]float64, len(specs))
	for i, s := range specs {
		v, err := p.GetVar(s)
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// SetVars writes a flat vector back into the prescription.
func (p *Prescription) SetVars(specs []VarSpec, x []float64) error {
	if len(specs) != len(x) {
		return fmt.Errorf("lens: %d specs for %d values", len(specs), len(x))
	}
	for i, s := range specs {
		if err := p.SetVar(s, x[i]); err != nil {
			return err
		}
	}
	return nil
}

// Clone returns a deep copy (Coeffs slices are duplicated).
func (p *Prescription) Clone() *Prescription {
	out := *p
	out.Cards = make([]Card, len(p.Cards))
	copy(out.Cards, p.Cards)
	for i := range out.Cards {
		if p.Cards[i].Coeffs != nil {
			out.Cards[i].Coeffs = append([]float64(nil), p.Cards[i].Coeffs...)
		}
	}
	return &out
}

// JSON I/O ------------------------------------------------------------------

// prescJSON is the wire form of Prescription: like the struct itself, except
// that the object position may be the string "inf" for an object at
// infinity (JSON has no infinity). A missing field also means infinity.
type prescJSON struct {
	Name      string          `json:"name,omitempty"`
	ObjectZMM json.RawMessage `json:"object_z_mm,omitempty"`
	Cards     []Card          `json:"cards"`
	Stop      int             `json:"stop"`
}

// MarshalJSON implements json.Marshaler. An object at infinity (ObjectZMM =
// math.Inf(-1)) is written as the string "inf".
func (p *Prescription) MarshalJSON() ([]byte, error) {
	w := prescJSON{Name: p.Name, Cards: p.Cards, Stop: p.Stop}
	switch {
	case math.IsInf(p.ObjectZMM, -1):
		w.ObjectZMM = json.RawMessage(`"inf"`)
	case math.IsNaN(p.ObjectZMM) || math.IsInf(p.ObjectZMM, 1):
		return nil, fmt.Errorf("lens: non-finite object position %v", p.ObjectZMM)
	default:
		buf, err := json.Marshal(p.ObjectZMM)
		if err != nil {
			return nil, err
		}
		w.ObjectZMM = buf
	}
	return json.Marshal(w)
}

// UnmarshalJSON implements json.Unmarshaler; a missing object position means
// an object at infinity.
func (p *Prescription) UnmarshalJSON(data []byte) error {
	var w prescJSON
	if err := json.Unmarshal(data, &w); err != nil {
		return err
	}
	p.Name = w.Name
	p.Cards = w.Cards
	p.Stop = w.Stop
	p.ObjectZMM = math.Inf(-1)
	if len(w.ObjectZMM) > 0 && string(w.ObjectZMM) != "null" {
		var f float64
		if err := json.Unmarshal(w.ObjectZMM, &f); err == nil {
			p.ObjectZMM = f
		} else {
			var s string
			if err := json.Unmarshal(w.ObjectZMM, &s); err != nil {
				return fmt.Errorf("lens: object_z_mm must be a number or \"inf\"")
			}
			switch s {
			case "inf", "-inf":
				p.ObjectZMM = math.Inf(-1)
			default:
				return fmt.Errorf("lens: unrecognised object_z_mm %q (use a number or \"inf\")", s)
			}
		}
	}
	return nil
}

// SaveJSON writes the prescription as pretty JSON.
func (p *Prescription) SaveJSON(path string) error {
	buf, err := json.MarshalIndent(p, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, append(buf, '\n'), 0o644)
}

// LoadJSON reads a prescription from a JSON file.
func LoadJSON(path string) (*Prescription, error) {
	buf, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var p Prescription
	if err := json.Unmarshal(buf, &p); err != nil {
		return nil, fmt.Errorf("lens: %s: %w", path, err)
	}
	if err := p.Validate(); err != nil {
		return nil, err
	}
	return &p, nil
}

// Prescription returns a one-line human-readable summary, useful in CLIs.
func (p *Prescription) String() string {
	s := fmt.Sprintf("prescription %q: %d cards", p.Name, len(p.Cards))
	for i, c := range p.Cards {
		r := "plane"
		if !math.IsInf(c.RadiusMM, 0) && c.RadiusMM != 0 {
			r = fmt.Sprintf("R=%.4g", c.RadiusMM)
		}
		k := ""
		if c.ConicK != 0 {
			k = fmt.Sprintf(" K=%.3g", c.ConicK)
		}
		asp := ""
		if len(c.Coeffs) > 0 {
			asp = fmt.Sprintf(" %d asphere coeffs", len(c.Coeffs))
		}
		mir := ""
		if c.Mirror {
			mir = " [mirror]"
		}
		g := c.Glass
		if g == "" {
			g = "air"
		}
		s += fmt.Sprintf("\n  [%d] %s%s%s%s → %.4gmm into %s", i, r, k, asp, mir, c.ThicknessMM, g)
	}
	return s
}
