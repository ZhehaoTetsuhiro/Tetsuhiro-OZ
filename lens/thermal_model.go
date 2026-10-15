package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// Thermal drift of a prescription (G9).
//
// Model (see geo/thermo_model.go for the equations): the temperature change
// dT = T − 20 °C is applied to the card list as follows.
//
//   - Each glass card i is an element of glass g_i. Its surfaces scale by
//     s = 1 + α ΔT: the radius, the aspheric coefficients A_k·s^-(2k+3), and the
//     clear semi-aperture. A surface between two glasses (cemented) takes the
//     mean of the two CTEs (an assumption: a cemented pair cannot follow two
//     different expansions; the mean is the minimal choice). A surface that
//     touches no glass (a mirror in air) takes MountCTE.
//   - The gap after card i scales by 1 + α_gap ΔT, with α_gap the CTE of g_i
//     (glass gap) or MountCTE (air gap, a mechanical spacer).
//   - The index of each glass card is the catalogue index plus Δn_rel(λ, ΔT)
//     (relative to air), through the mediumWrap hook.
//
// The object distance and the stop are not changed. The prescription's
// Validate, JSON form and the optimisation variables are not affected: the
// result is an ordinary Prescription.

// ThermalOptions controls the temperature model. The zero value is a rigid
// mount: air gaps and air-facing mirror surfaces do not change.
type ThermalOptions struct {
	// MountCTE is the linear expansion coefficient (1/K) of the spacers that
	// set the air gaps and of any surface that touches no glass. 0 means a
	// rigid mount.
	MountCTE float64
}

// isAirName reports whether a glass name denotes air (the same set as
// resolveGlass).
func isAirName(name string) bool {
	switch name {
	case "", "air", "Air", "AIR", "vacuum":
		return true
	}
	return false
}

// Thermal returns a copy of the prescription at a temperature change dT (K)
// from 20 °C. dT = 0 returns a prescription whose system is bit-identical to
// the original one (every scale factor is exactly 1 and the index shift is
// not applied). An error is returned for a glass without thermo-optic data
// (every catalogue glass has data, see geo.ThermalGlassNames).
func (p *Prescription) Thermal(dT float64, opts *ThermalOptions) (*Prescription, error) {
	if math.IsNaN(dT) || math.IsInf(dT, 0) {
		return nil, fmt.Errorf("lens: thermal: temperature change %v is not finite", dT)
	}
	var o ThermalOptions
	if opts != nil {
		o = *opts
	}
	if math.IsNaN(o.MountCTE) || math.IsInf(o.MountCTE, 0) {
		return nil, fmt.Errorf("lens: thermal: mount CTE %v is not finite", o.MountCTE)
	}
	if err := p.Validate(); err != nil {
		return nil, err
	}
	n := len(p.Cards)

	// Thermo-optic data per card (nil for air).
	glass := make([]*geo.ThermalGlass, n)
	for i, c := range p.Cards {
		if isAirName(c.Glass) {
			continue
		}
		g, ok := geo.ThermalGlassByName(c.Glass)
		if !ok {
			return nil, fmt.Errorf("lens: thermal: card %d: glass %q has no thermo-optic data", i, c.Glass)
		}
		glass[i] = &g
	}
	glassCTE := func(i int) (float64, bool) {
		if i < 0 || i >= n || glass[i] == nil {
			return 0, false
		}
		return glass[i].CTE, true
	}

	q := p.Clone()
	for i := 0; i < n; i++ {
		// Surface i sits between card i-1 (before) and card i (after).
		before, hasBefore := glassCTE(i - 1)
		after, hasAfter := glassCTE(i)
		var surfCTE float64
		switch {
		case hasBefore && hasAfter:
			surfCTE = 0.5 * (before + after)
		case hasAfter:
			surfCTE = after
		case hasBefore:
			surfCTE = before
		default:
			surfCTE = o.MountCTE
		}
		s := 1 + surfCTE*dT
		c := &q.Cards[i]
		if !math.IsInf(c.RadiusMM, 0) && c.RadiusMM != 0 {
			c.RadiusMM *= s
		}
		if len(c.Coeffs) > 0 && s != 1 {
			for k := range c.Coeffs {
				c.Coeffs[k] *= math.Pow(s, -float64(2*k+3))
			}
		}
		c.SemiAperMM *= s
		if i < n-1 {
			gapCTE := o.MountCTE
			if a, ok := glassCTE(i); ok {
				gapCTE = a
			}
			c.ThicknessMM *= 1 + gapCTE*dT
		}
	}

	prev := p.mediumWrap
	q.mediumWrap = func(card int, m geo.Medium) geo.Medium {
		if prev != nil {
			m = prev(card, m)
		}
		if card >= 0 && card < n && glass[card] != nil {
			return geo.ThermalMedium{Base: m, Glass: *glass[card], DT: dT}
		}
		return m
	}
	if dT != 0 && q.Name != "" {
		q.Name = fmt.Sprintf("%s @%+.4gK", q.Name, dT)
	}
	return q, nil
}

// ThermalDrift is the first-order drift of a prescription between 20 °C and
// 20 °C + DT at one wavelength. The per-kelvin values are secant slopes
// (value difference divided by DT), not derivatives.
type ThermalDrift struct {
	DT         float64 // temperature change (K)
	EFL0, EFL1 float64 // effective focal length (mm) at 20 °C and at 20 °C+DT
	BFL0, BFL1 float64 // back focal length (mm) at 20 °C and at 20 °C+DT
	EFLPerK    float64 // (EFL1 − EFL0)/DT (mm/K)
	BFLPerK    float64 // (BFL1 − BFL0)/DT (mm/K)
}

// ThermalDriftOf computes the paraxial EFL/BFL at 20 °C and at 20 °C + dT with
// the given options. dT must be non-zero for the per-kelvin values.
func (p *Prescription) ThermalDriftOf(lambdaNm, dT float64, opts *ThermalOptions) (ThermalDrift, error) {
	if dT == 0 {
		return ThermalDrift{}, fmt.Errorf("lens: thermal drift needs dT != 0")
	}
	par0, err := p.Paraxial(lambdaNm)
	if err != nil {
		return ThermalDrift{}, err
	}
	hot, err := p.Thermal(dT, opts)
	if err != nil {
		return ThermalDrift{}, err
	}
	par1, err := hot.Paraxial(lambdaNm)
	if err != nil {
		return ThermalDrift{}, err
	}
	d := ThermalDrift{DT: dT, EFL0: par0.EFL, EFL1: par1.EFL, BFL0: par0.BFL, BFL1: par1.BFL}
	d.EFLPerK = (d.EFL1 - d.EFL0) / dT
	d.BFLPerK = (d.BFL1 - d.BFL0) / dT
	return d, nil
}
