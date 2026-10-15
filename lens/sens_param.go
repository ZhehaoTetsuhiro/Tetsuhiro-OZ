package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// ParamKind selects the prescription quantity a sensitivity Param perturbs.
type ParamKind int

// Parameter kinds. All kinds except ParamIndex refer to a quantity of one
// card; ParamIndex refers to the glass that follows the card.
const (
	ParamRadius    ParamKind = iota // vertex radius R (mm); the offset is added to R
	ParamCurvature                  // vertex curvature c = 1/R (1/mm); R = 1/(c0 + offset)
	ParamThickness                  // axial distance to the next card (mm)
	ParamConic                      // conic constant K
	ParamCoeff                      // aspheric coefficient Coeffs[Order]
	ParamIndex                      // refractive index of the medium after the card (offset Δn)
)

// Param is one perturbable quantity of a prescription. An offset x moves the
// quantity from its nominal value v0 to v0 + x; for ParamCurvature the offset
// acts on the curvature (R = 1/(c0 + x)), and for ParamIndex the nominal
// offset is zero and x is the index change Δn of the glass.
//
// Use ParamCurvature rather than ParamRadius for a plane surface, whose
// radius is not a usable coordinate.
type Param struct {
	Kind  ParamKind
	Card  int     // card (surface) index
	Order int     // ParamCoeff only: index into Coeffs (0 = A4)
	Scale float64 // characteristic size for the finite-difference step; 0 = default of the kind
}

// String labels the parameter for reports, e.g. "s1.R", "s0.c", "s1.thk", "s0.dn".
func (q Param) String() string {
	switch q.Kind {
	case ParamRadius:
		return fmt.Sprintf("s%d.R", q.Card)
	case ParamCurvature:
		return fmt.Sprintf("s%d.c", q.Card)
	case ParamThickness:
		return fmt.Sprintf("s%d.thk", q.Card)
	case ParamConic:
		return fmt.Sprintf("s%d.K", q.Card)
	case ParamCoeff:
		return fmt.Sprintf("s%d.A%d", q.Card, 2*q.Order+4)
	case ParamIndex:
		return fmt.Sprintf("s%d.dn", q.Card)
	}
	return fmt.Sprintf("s%d.?", q.Card)
}

// nominal returns the nominal value of the quantity in p (zero for ParamIndex).
func (q Param) nominal(p *Prescription) (float64, error) {
	if q.Card < 0 || q.Card >= len(p.Cards) {
		return 0, fmt.Errorf("lens: parameter %s: card out of range", q)
	}
	c := &p.Cards[q.Card]
	switch q.Kind {
	case ParamRadius:
		if c.RadiusMM == 0 || math.IsInf(c.RadiusMM, 0) {
			return 0, fmt.Errorf("lens: parameter %s: surface is plane, use ParamCurvature", q)
		}
		return c.RadiusMM, nil
	case ParamCurvature:
		if c.RadiusMM == 0 || math.IsInf(c.RadiusMM, 0) {
			return 0, nil
		}
		return 1 / c.RadiusMM, nil
	case ParamThickness:
		if q.Card >= len(p.Cards)-1 {
			return 0, fmt.Errorf("lens: parameter %s: the last card's thickness does not enter the system", q)
		}
		return c.ThicknessMM, nil
	case ParamConic:
		return c.ConicK, nil
	case ParamCoeff:
		if q.Order < 0 || q.Order >= len(c.Coeffs) {
			return 0, fmt.Errorf("lens: parameter %s: coefficient out of range", q)
		}
		return c.Coeffs[q.Order], nil
	case ParamIndex:
		return 0, nil
	}
	return 0, fmt.Errorf("lens: parameter %s: unknown kind %d", q, q.Kind)
}

// scale returns the characteristic size of the parameter at nominal value v0;
// the finite-difference step of Sensitivity is RelStep times this value.
func (q Param) scale(v0 float64) float64 {
	if q.Scale > 0 {
		return q.Scale
	}
	switch q.Kind {
	case ParamRadius, ParamThickness, ParamConic:
		return math.Max(math.Abs(v0), 1) // mm (conic: unit)
	case ParamCurvature:
		return math.Max(math.Abs(v0), 1e-2) // 1/mm
	case ParamCoeff:
		return math.Max(math.Abs(v0), 1e-3)
	}
	return 1 // index change
}

// applyParams returns a clone of p in which params[i] is moved by offsets[i]
// relative to its nominal value in p. The parameters must be distinct (see
// checkParams). The clone keeps every property of p, including the glass
// hook of the tolerance analysis; index offsets are composed on top of it.
func applyParams(p *Prescription, params []Param, offsets []float64) (*Prescription, error) {
	if len(params) != len(offsets) {
		return nil, fmt.Errorf("lens: %d parameters for %d offsets", len(params), len(offsets))
	}
	q := p.Clone()
	dn := make([]float64, len(p.Cards))
	anyIndex := false
	for i, pr := range params {
		v0, err := pr.nominal(p)
		if err != nil {
			return nil, err
		}
		x := offsets[i]
		c := &q.Cards[pr.Card]
		switch pr.Kind {
		case ParamRadius:
			c.RadiusMM = v0 + x
		case ParamCurvature:
			if cc := v0 + x; cc == 0 {
				c.RadiusMM = math.Inf(1)
			} else {
				c.RadiusMM = 1 / cc
			}
		case ParamThickness:
			c.ThicknessMM = v0 + x
		case ParamConic:
			c.ConicK = v0 + x
		case ParamCoeff:
			c.Coeffs[pr.Order] = v0 + x
		case ParamIndex:
			dn[pr.Card] += x
			anyIndex = true
		}
	}
	if anyIndex {
		prev := p.mediumWrap
		shift := dn
		q.mediumWrap = func(card int, m geo.Medium) geo.Medium {
			if prev != nil {
				m = prev(card, m)
			}
			if card >= 0 && card < len(shift) && shift[card] != 0 {
				m = meltMedium{base: m, dn: shift[card]}
			}
			return m
		}
	}
	return q, nil
}

// field identifies the prescription storage a parameter writes. Radius and
// curvature are two coordinates of the same RadiusMM field, so they alias: an
// offset of one would be overwritten by the other within one perturbation.
func (q Param) field() (kind ParamKind, card, order int) {
	kind = q.Kind
	if kind == ParamCurvature {
		kind = ParamRadius
	}
	order = 0
	if kind == ParamCoeff {
		order = q.Order
	}
	return kind, q.Card, order
}

// checkParams rejects an empty list, invalid references and parameters that
// write the same prescription field (for example ParamRadius and
// ParamCurvature of one card).
func checkParams(p *Prescription, params []Param) error {
	if len(params) == 0 {
		return fmt.Errorf("lens: no parameters")
	}
	for i, a := range params {
		if _, err := a.nominal(p); err != nil {
			return err
		}
		ka, ca, oa := a.field()
		for _, b := range params[:i] {
			kb, cb, ob := b.field()
			if ka == kb && ca == cb && oa == ob {
				return fmt.Errorf("lens: parameters %s and %s address the same prescription field", b, a)
			}
		}
	}
	return nil
}
