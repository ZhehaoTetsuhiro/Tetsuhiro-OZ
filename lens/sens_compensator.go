package lens

import (
	"fmt"
	"math"
)

// Perturbation is a tolerance draw: Params[i] is moved by Offsets[i] from its
// nominal value.
type Perturbation struct {
	Params  []Param
	Offsets []float64
}

// CompensatorSpec selects the compensating variables of Compensate.
type CompensatorSpec struct {
	// Vars are prescription quantities re-optimised after the perturbation
	// (radii, air spacings, ...). Each starts at its value and may move by at
	// most its half-range in VarRange (same length as Vars). A nil VarRange
	// gives max(1, 0.1·|v|) for radii and thicknesses and max(0.01, 0.1·|v|)
	// for the other kinds, with v the nominal value.
	Vars     []VarSpec
	VarRange []float64
	// ImageShift adds the axial position of the image plane as a compensator
	// (defocus compensation). It needs a fixed image plane in the merit.
	ImageShift bool
	ShiftRange float64 // half-range of the image shift (mm); 0 = 1 mm
	MaxIter    int     // DLS iterations; 0 = 100
}

// CompensationResult is the outcome of Compensate.
type CompensationResult struct {
	// Nominal is the merit of the unperturbed prescription, uncompensated (mm).
	Nominal float64
	// CompensatedNominal is the merit of the unperturbed prescription after the
	// same compensators have been optimised on it (mm).
	CompensatedNominal float64
	// Perturbed is the merit after the perturbation, before compensation (mm).
	Perturbed float64
	// Compensated is the merit after re-optimising the compensators (mm).
	Compensated float64
	// Recovery is (Perturbed − Compensated) / (Perturbed − CompensatedNominal):
	// the share of the tolerance-induced degradation, measured against the
	// compensated nominal, that the compensators remove. 1 means the perturbed
	// system is restored to the compensated nominal. NaN when the perturbation
	// does not degrade the compensated merit.
	Recovery float64
	// ImageShift is the compensated image-plane shift (mm); 0 without ImageShift.
	ImageShift float64
	// Prescription is the compensated prescription (the image shift is not part of it).
	Prescription *Prescription
	// Opt is the optimizer trace of the perturbed run.
	Opt *Result
}

// Compensate applies the perturbation to p, then re-optimises the compensators
// by damped least squares on the merit residuals. The optimizer starts from
// the perturbed state and only accepts decreasing steps, so Compensated ≤
// Perturbed holds by construction. The unperturbed prescription is compensated
// with the same variables and bounds, which gives the reference of Recovery.
//
// The merit should have a fixed image plane (see MeritSpec.Fixed); with a
// floating focus the image shift has no effect and ImageShift is rejected.
func Compensate(p *Prescription, pert Perturbation, merit MeritSpec, comp CompensatorSpec) (*CompensationResult, error) {
	if merit.FloatingFocus && comp.ImageShift {
		return nil, fmt.Errorf("lens: image-shift compensation needs a fixed image plane in the merit")
	}
	nv := len(comp.Vars)
	if nv == 0 && !comp.ImageShift {
		return nil, fmt.Errorf("lens: no compensator selected")
	}
	if comp.VarRange != nil && len(comp.VarRange) != nv {
		return nil, fmt.Errorf("lens: %d ranges for %d compensator variables", len(comp.VarRange), nv)
	}
	m := merit.withDefaults()
	nominal, err := m.Value(p)
	if err != nil {
		return nil, fmt.Errorf("lens: compensation: nominal merit: %w", err)
	}
	ranges := make([]float64, nv)
	for i, s := range comp.Vars {
		v, err := p.GetVar(s)
		if err != nil {
			return nil, err
		}
		ranges[i] = defaultCompRange(s, v)
		if comp.VarRange != nil {
			ranges[i] = comp.VarRange[i]
		}
		if !(ranges[i] > 0) {
			return nil, fmt.Errorf("lens: compensator %s range must be positive", s)
		}
	}
	pq, err := applyParams(p, pert.Params, pert.Offsets)
	if err != nil {
		return nil, err
	}
	perturbed, err := m.Value(pq)
	if err != nil {
		return nil, fmt.Errorf("lens: compensation: perturbed merit: %w", err)
	}
	ref, err := runCompensator(p, m, comp, ranges)
	if err != nil {
		return nil, fmt.Errorf("lens: compensation of the nominal: %w", err)
	}
	run, err := runCompensator(pq, m, comp, ranges)
	if err != nil {
		return nil, fmt.Errorf("lens: compensation of the perturbed system: %w", err)
	}
	out := &CompensationResult{
		Nominal:            nominal,
		CompensatedNominal: ref.merit,
		Perturbed:          perturbed,
		Compensated:        run.merit,
		Recovery:           math.NaN(),
		ImageShift:         run.shift,
		Prescription:       run.p,
		Opt:                run.opt,
	}
	if perturbed-ref.merit > 0 {
		out.Recovery = (perturbed - run.merit) / (perturbed - ref.merit)
	}
	return out, nil
}

// compRun is the outcome of one compensator optimisation.
type compRun struct {
	p     *Prescription
	shift float64
	merit float64
	opt   *Result
}

// runCompensator optimises the compensators of comp on start by damped least
// squares, each within ±ranges[i] of its starting value and the image shift
// within ±ShiftRange. The merit of the result is recomputed from the final
// prescription.
func runCompensator(start *Prescription, m MeritSpec, comp CompensatorSpec, ranges []float64) (*compRun, error) {
	nv := len(comp.Vars)
	nVars := nv
	if comp.ImageShift {
		nVars++
	}
	vars := make([]Variable, nVars)
	for i, s := range comp.Vars {
		v, err := start.GetVar(s)
		if err != nil {
			return nil, err
		}
		vars[i] = Variable{Name: s.String(), Value: v, Lo: v - ranges[i], Hi: v + ranges[i]}
	}
	if comp.ImageShift {
		r := comp.ShiftRange
		if r <= 0 {
			r = 1
		}
		vars[nv] = Variable{Name: "image", Value: 0, Lo: -r, Hi: r}
	}
	build := func(x []float64) (*Prescription, MeritSpec, error) {
		q := start.Clone()
		if err := q.SetVars(comp.Vars, x[:nv]); err != nil {
			return nil, m, err
		}
		mm := m
		if comp.ImageShift {
			mm.ImageZ = m.ImageZ + x[nv]
		}
		return q, mm, nil
	}
	residual := func(x []float64) ([]float64, error) {
		q, mm, err := build(x)
		if err != nil {
			return nil, err
		}
		r, err := mm.Residuals(q)
		if err != nil {
			// A failed trace keeps the residual length fixed and is penalised.
			r = make([]float64, m.Len())
			for i := range r {
				r[i] = failPenalty
			}
		}
		return r, nil
	}
	opt, err := Minimize(Problem{Variables: vars, Residual: residual, MaxIter: comp.MaxIter})
	if err != nil {
		return nil, err
	}
	q, mm, err := build(opt.X)
	if err != nil {
		return nil, err
	}
	merit, err := mm.Value(q)
	if err != nil {
		return nil, err
	}
	out := &compRun{p: q, merit: merit, opt: opt}
	if comp.ImageShift {
		out.shift = opt.X[nv]
	}
	return out, nil
}

// defaultCompRange is the default half-range of a compensator variable.
func defaultCompRange(s VarSpec, v float64) float64 {
	switch s.Kind {
	case VarRadius, VarThickness:
		return math.Max(1, 0.1*math.Abs(v))
	}
	return math.Max(0.01, 0.1*math.Abs(v))
}
