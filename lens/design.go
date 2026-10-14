package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// DesignOptions controls the automatic design drivers. The zero value is
// valid: every field falls back to a documented default (see applyDefaults).
type DesignOptions struct {
	Fields       []Field  // fields to correct (default: on-axis + 0.5°)
	Samples      []Sample // pupil sampling (default UnitDiskSamples(3, 6))
	SpotTargetMM float64  // residual scaling target (default max(2µm, efl/20000))
	FocalTolMM   float64  // focal-length tolerance (default max(1µm, efl/20000))
	ColorTolMM   float64  // axial colour tolerance in mm (default 0.01)
	MaxIter      int      // accepted DLS iterations (default 150)
	SkipPolish   bool     // return the analytic solve without DLS polishing
}

func (o *DesignOptions) applyDefaults(efl float64) {
	if o.Samples == nil {
		o.Samples = UnitDiskSamples(3, 6)
	}
	if o.Fields == nil {
		o.Fields = []Field{{0, 0}, {0, 0.0087}} // on-axis + 0.5°
	}
	if o.SpotTargetMM <= 0 {
		o.SpotTargetMM = math.Max(2e-3, efl/20000)
	}
	if o.FocalTolMM <= 0 {
		o.FocalTolMM = math.Max(1e-3, efl/20000)
	}
	if o.ColorTolMM <= 0 {
		o.ColorTolMM = 0.01
	}
	if o.MaxIter <= 0 {
		o.MaxIter = 150
	}
}

// DesignAchromat designs a cemented achromatic doublet (crown + flint) with
// effective focal length eflMM at the d-line. It is the "自主设计" path of
// the package: a closed-form thin-lens achromat solve (power split by Abbe
// numbers, element shapes chosen to share the crown power) followed — unless
// SkipPolish is set — by a damped-least-squares polish against the spot
// diagrams of the configured fields, the focal length and the F–C focal
// shift.
//
// semiAperMM is the clear semi-diameter of every card (use efl/(2·F#)); 0
// selects efl/10. The returned Result is the optimizer trace (nil when the
// polish is skipped). The design is deterministic.
func DesignAchromat(eflMM float64, crownGlass, flintGlass string, semiAperMM float64, opts *DesignOptions) (*Prescription, *Result, error) {
	if !(eflMM > 0) {
		return nil, nil, fmt.Errorf("lens: efl must be positive")
	}
	g1, ok := geo.Glass(crownGlass)
	if !ok {
		return nil, nil, fmt.Errorf("lens: unknown crown glass %q (known: %v)", crownGlass, geo.GlassNames())
	}
	g2, ok := geo.Glass(flintGlass)
	if !ok {
		return nil, nil, fmt.Errorf("lens: unknown flint glass %q (known: %v)", flintGlass, geo.GlassNames())
	}
	n1 := g1.N(geo.Wavelengthd)
	n2 := g2.N(geo.Wavelengthd)
	V1 := g1.AbbeNumber()
	V2 := g2.AbbeNumber()
	if !(V1 > V2) {
		return nil, nil, fmt.Errorf("lens: crown %q (V=%.1f) must have a higher Abbe number than flint %q (V=%.1f)",
			crownGlass, V1, flintGlass, V2)
	}

	// Thin-lens achromat: phi1 + phi2 = phi with phi1/phi2 = V1/V2.
	phi := 1 / eflMM
	phi1 := phi * V1 / (V1 + V2)
	phi2 := phi - phi1
	// Crown split equally between its two surfaces (cB = -cA); flint takes
	// the remainder on the cemented surface and the rear.
	cA := phi1 / (2 * (n1 - 1))
	cB := -cA
	cC := cB - phi2/(n2-1)

	t1 := math.Max(1.5, eflMM*0.03)
	t2 := math.Max(0.8, eflMM*0.012)
	semi := semiAperMM
	if !(semi > 0) {
		semi = eflMM / 10
	}
	presc := &Prescription{
		Name:      fmt.Sprintf("%s/%s achromat f=%.1fmm", crownGlass, flintGlass, eflMM),
		ObjectZMM: math.Inf(-1),
		Stop:      0,
		Cards: []Card{
			{Name: "crown front", RadiusMM: 1 / cA, ThicknessMM: t1, Glass: crownGlass, SemiAperMM: semi},
			{Name: "cemented surface", RadiusMM: 1 / cB, ThicknessMM: t2, Glass: flintGlass, SemiAperMM: semi},
			{Name: "flint rear", RadiusMM: 1 / cC, ThicknessMM: eflMM * 0.95, Glass: "", SemiAperMM: semi},
		},
	}
	if err := presc.Validate(); err != nil {
		return nil, nil, err
	}
	// Refine the three curvatures on the *thick-lens* paraxial engine so
	// the start already satisfies EFL(d) = eflMM and EFL(F) = EFL(C)
	// exactly; the thin-lens solve above is only the Newton seed.
	if err := achromatStart(presc, eflMM); err != nil {
		return nil, nil, err
	}

	o := DesignOptions{}
	if opts != nil {
		o = *opts
	}
	o.applyDefaults(eflMM)
	if o.SkipPolish {
		return presc, nil, nil
	}

	// DLS polish: variables are the three radii and the two element
	// thicknesses. Bounds keep the model in a physically sane neighbourhood
	// of the analytic solve; in particular the element thicknesses are
	// capped at 8% of the focal length, because a very thick element adds a
	// chromatic shift of the principal plane (~ t·Δn/n²) that no thin-lens
	// power balance can cancel and would otherwise trap the optimizer in a
	// badly-coloured local minimum.
	specs := []VarSpec{
		{0, VarRadius, 0}, {1, VarRadius, 0}, {2, VarRadius, 0},
		{0, VarThickness, 0}, {1, VarThickness, 0},
	}
	rLim := 30 * eflMM
	tLo := 0.5
	tHi := eflMM * 0.08
	work := presc.Clone()
	cur, err := work.Vars(specs)
	if err != nil {
		return nil, nil, err
	}
	vars := make([]Variable, len(specs))
	for i, s := range specs {
		lo, hi := -rLim, rLim
		switch s.Kind {
		case VarThickness:
			lo, hi = tLo, tHi
		}
		vars[i] = Variable{Name: s.String(), Value: cur[i], Lo: lo, Hi: hi}
	}
	prob := Problem{
		Variables: vars,
		Residual: func(x []float64) ([]float64, error) {
			p := work.Clone()
			if err := p.SetVars(specs, x); err != nil {
				return nil, err
			}
			return designResiduals(p, o, eflMM)
		},
		MaxIter: o.MaxIter,
	}
	res, err := Minimize(prob)
	if err != nil {
		return nil, nil, err
	}
	final := presc.Clone()
	if err := final.SetVars(specs, res.X); err != nil {
		return nil, nil, err
	}
	if err := final.Validate(); err != nil {
		return nil, nil, err
	}
	final.Name = fmt.Sprintf("%s/%s achromat f=%.1fmm (designed)", crownGlass, flintGlass, eflMM)
	return final, res, nil
}

// achromatStart refines the curvatures of a freshly created cemented
// doublet so that the *thick-lens* first order satisfies
//
//	EFL(d)  = eflMM          and          EFL(F) = EFL(C)
//
// holding the crown's equal-split bending (the cemented surface keeps
// cB = -cA, the thin-lens shape rule). Two residuals, two unknowns
// (cA and cC), solved with a deterministic damped Newton iteration that
// evaluates the actual geo paraxial engine — so thickness effects on the
// focal length and on the chromatic principal-plane shift, which a
// thin-lens solve cannot see, are corrected exactly at the start.
func achromatStart(p *Prescription, eflMM float64) error {
	cA := 1 / p.Cards[0].RadiusMM
	cC := 1 / p.Cards[2].RadiusMM
	if cA == 0 || cC == 0 || math.IsInf(cA, 0) || math.IsInf(cC, 0) {
		return fmt.Errorf("lens: achromat start solve needs finite front/rear radii")
	}
	resid := func(cA, cC float64) ([2]float64, error) {
		p.Cards[0].RadiusMM = 1 / cA
		p.Cards[1].RadiusMM = -1 / cA
		p.Cards[2].RadiusMM = 1 / cC
		pd, err := p.Paraxial(geo.Wavelengthd)
		if err != nil {
			return [2]float64{}, err
		}
		pf, err := p.Paraxial(geo.WavelengthF)
		if err != nil {
			return [2]float64{}, err
		}
		pc, err := p.Paraxial(geo.WavelengthC)
		if err != nil {
			return [2]float64{}, err
		}
		return [2]float64{
			(pd.EFL - eflMM) / eflMM,
			(pf.EFL - pc.EFL) / eflMM,
		}, nil
	}
	const tol = 1e-10
	for it := 0; it < 60; it++ {
		r0, err := resid(cA, cC)
		if err != nil {
			return err
		}
		if math.Abs(r0[0]) < tol && math.Abs(r0[1]) < tol {
			return nil
		}
		hA := math.Abs(cA)*1e-6 + 1e-12
		hC := math.Abs(cC)*1e-6 + 1e-12
		ra, err := resid(cA+hA, cC)
		if err != nil {
			return err
		}
		rc, err := resid(cA, cC+hC)
		if err != nil {
			return err
		}
		j11 := (ra[0] - r0[0]) / hA
		j12 := (rc[0] - r0[0]) / hC
		j21 := (ra[1] - r0[1]) / hA
		j22 := (rc[1] - r0[1]) / hC
		det := j11*j22 - j12*j21
		if math.Abs(det) < 1e-30 {
			return fmt.Errorf("lens: achromat start solve is singular at iteration %d", it)
		}
		dcA := (-j22*r0[0] + j12*r0[1]) / det
		dcC := (j21*r0[0] - j11*r0[1]) / det
		// Damped step: never move a curvature by more than 25% per
		// iteration, so the solve cannot jump out of the physical branch.
		if m := 0.25 * math.Abs(cA); math.Abs(dcA) > m {
			dcA = math.Copysign(m, dcA)
		}
		if m := 0.25 * math.Abs(cC); math.Abs(dcC) > m {
			dcC = math.Copysign(m, dcC)
		}
		cA += dcA
		cC += dcC
	}
	return fmt.Errorf("lens: achromat start solve did not converge")
}

// designResiduals assembles the merit vector of the design drivers: spot
// residuals for every field, the focal length and the axial colour. Fields
// that cannot be characterised contribute failPenalty of the right length,
// so the optimizer keeps a usable, fixed-size residual vector.
func designResiduals(p *Prescription, o DesignOptions, eflMM float64) ([]float64, error) {
	var out []float64
	perField := 2 * len(o.Samples)
	for _, f := range o.Fields {
		r, err := SpotResiduals(p, geo.Wavelengthd, []Field{f}, o.Samples, o.SpotTargetMM)
		if err != nil {
			for i := 0; i < perField; i++ {
				out = append(out, failPenalty)
			}
			continue
		}
		out = append(out, r...)
	}
	if e, err := EFLResidual(p, geo.Wavelengthd, eflMM, o.FocalTolMM); err == nil {
		out = append(out, e, 0)
	} else {
		out = append(out, failPenalty, failPenalty)
	}
	if e, err := ChromaticResidual(p, o.ColorTolMM); err == nil {
		out = append(out, e, 0)
	} else {
		out = append(out, failPenalty, failPenalty)
	}
	return out, nil
}

// Describe returns a compact human-readable report of the prescription's
// first-order properties and spot quality at the d-line (used by the CLI
// and by the test suite as a review aid).
func Describe(p *Prescription) (string, error) {
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		return "", err
	}
	samples := UnitDiskSamples(3, 6)
	rmsOn, geoOn, nOn, err := SpotRMS(p, geo.Wavelengthd, Field{0, 0}, samples)
	if err != nil {
		return "", err
	}
	rmsOff, geoOff, nOff, err := SpotRMS(p, geo.Wavelengthd, Field{0, 0.0087}, samples)
	if err != nil {
		return "", err
	}
	chrom, err := ChromaticResidual(p, 1.0)
	if err != nil {
		return "", err
	}
	return fmt.Sprintf(
		"%s\n  λd %.1fnm: EFL=%.4fmm  BFL=%.4fmm  F/#=%.2f  ENP r=%.3fmm @z=%.3fmm\n"+
			"  F–C focal shift = %.4f mm\n"+
			"  spot (d): on-axis RMS=%.4fmm GEO=%.4fmm (%d rays) | 0.5° RMS=%.4fmm GEO=%.4fmm (%d rays)",
		p.Name, geo.Wavelengthd, par.EFL, par.BFL, par.FNumber,
		par.EntrancePupilRadius, par.EntrancePupilZ, chrom,
		rmsOn, geoOn, nOn, rmsOff, geoOff, nOff), nil
}
