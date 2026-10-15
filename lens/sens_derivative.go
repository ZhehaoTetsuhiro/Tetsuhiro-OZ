package lens

import (
	"fmt"
	"math"
)

// ScalarFunc evaluates a scalar quantity of a prescription (a merit, a focal
// length, ...).
type ScalarFunc func(p *Prescription) (float64, error)

// SensOptions tunes Sensitivity. The zero value is valid.
//
// Step selection. A central difference has two error sources: truncation (large
// steps; the function is not yet polynomial on the scale of the step) and
// roundoff (small steps; the merit is a sum of ray coordinates with absolute
// error ~1e-14 mm). For a merit on a fixed image plane the truncation scale can
// be small: a parameter step moves the best focus across the plane (see
// TestSensStepConsistency for measured values). Sensitivity therefore evaluates
// a ladder of relative steps RelStep, RelStep/10, ..., MinRelStep for each
// parameter, extrapolates each with Richardson, and reports the stage whose
// neighbourhood agrees best. The uncertainty u of a stage is the largest of its
// own Richardson estimate and its disagreement with the next coarser and next
// finer stage. Converged means u ≤ Tol·|∂f/∂x|. Fixed evaluates RelStep only.
type SensOptions struct {
	// RelStep is the coarsest relative step, as a fraction of each parameter's
	// scale (Param.scale). Default 1e-3.
	RelStep float64
	// MinRelStep is the finest relative step of the ladder. Default 1e-7.
	MinRelStep float64
	// Levels is the number of Richardson levels K >= 1 (default 3): central
	// differences are formed at h, h/2, ..., h/2^(K-1) and extrapolated.
	Levels int
	// Tol is the relative uncertainty required for Converged (default 1e-5).
	Tol float64
	// Fixed evaluates RelStep only, without the ladder; Converged is then
	// reported as true and carries no information.
	Fixed bool
}

// SensResult holds the first derivatives of a scalar with respect to a set of
// parameters, evaluated at the nominal prescription.
type SensResult struct {
	Params    []Param
	Nominal   float64   // f at the unperturbed prescription
	Deriv     []float64 // ∂f/∂x_i (Richardson-extrapolated at the selected stage)
	Err       []float64 // uncertainty u of Deriv (absolute, see SensOptions)
	Naive     []float64 // plain central difference at the selected base step
	RelSteps  []float64 // selected relative step per parameter
	Steps     []float64 // selected base step h_i in parameter units
	Converged []bool    // u ≤ Tol·|Deriv| (always true with Fixed)
	Evals     int       // number of evaluations of f
}

// Sensitivity returns ∂f/∂x_i for every parameter at p, where f is any scalar
// function of the prescription (a merit built with MeritSpec.Func, an
// EFLScalar, ...). Each derivative is a central difference extrapolated with
// Richardson's rule (Romberg table on the halving sequence h, h/2, ...). The
// error of a central difference is an even power series in h, so each
// extrapolation level removes the next power: with K levels the truncation
// error is O(h^(2K)) before roundoff. The step selection is described on
// SensOptions.
//
// The function is evaluated on perturbed clones of p; p itself is not
// modified. The result is deterministic.
func Sensitivity(p *Prescription, params []Param, f ScalarFunc, opts SensOptions) (*SensResult, error) {
	if f == nil {
		return nil, fmt.Errorf("lens: nil scalar function")
	}
	if err := checkParams(p, params); err != nil {
		return nil, err
	}
	o := opts
	if !(o.RelStep > 0) {
		o.RelStep = 1e-3
	}
	if !(o.MinRelStep > 0) || o.MinRelStep > o.RelStep {
		o.MinRelStep = 1e-7
		if o.MinRelStep > o.RelStep {
			o.MinRelStep = o.RelStep
		}
	}
	if o.Levels <= 0 {
		o.Levels = 3
	}
	if !(o.Tol > 0) {
		o.Tol = 1e-5
	}
	ladder := []float64{o.RelStep}
	if !o.Fixed {
		ladder = ladder[:0]
		for rs := o.RelStep; rs >= o.MinRelStep*(1-1e-9); rs /= 10 {
			ladder = append(ladder, rs)
		}
	}
	nominal, err := f(p)
	if err != nil {
		return nil, fmt.Errorf("lens: sensitivity: nominal evaluation failed: %w", err)
	}
	n := len(params)
	res := &SensResult{
		Params:    append([]Param(nil), params...),
		Nominal:   nominal,
		Deriv:     make([]float64, n),
		Err:       make([]float64, n),
		Naive:     make([]float64, n),
		RelSteps:  make([]float64, n),
		Steps:     make([]float64, n),
		Converged: make([]bool, n),
		Evals:     1,
	}
	x := make([]float64, n)
	for i, pr := range params {
		v0, err := pr.nominal(p)
		if err != nil {
			return nil, err
		}
		scale := pr.scale(v0)
		D := make([]float64, len(ladder))
		E := make([]float64, len(ladder))
		N := make([]float64, len(ladder))
		for s, rs := range ladder {
			h := rs * scale
			d := make([]float64, o.Levels)
			for k := range d {
				hk := h / math.Pow(2, float64(k))
				fp, err := evalAt(p, params, f, x, i, +hk)
				if err != nil {
					return nil, err
				}
				fm, err := evalAt(p, params, f, x, i, -hk)
				if err != nil {
					return nil, err
				}
				res.Evals += 2
				d[k] = (fp - fm) / (2 * hk)
			}
			D[s], E[s] = romberg(d)
			N[s] = d[0]
		}
		best, bestU := 0, math.Inf(1)
		for s := range ladder {
			u := E[s]
			if s > 0 {
				u = math.Max(u, math.Abs(D[s]-D[s-1]))
			}
			if s+1 < len(ladder) {
				u = math.Max(u, math.Abs(D[s]-D[s+1]))
			}
			if u < bestU {
				best, bestU = s, u
			}
		}
		res.Deriv[i], res.Err[i], res.Naive[i] = D[best], bestU, N[best]
		res.RelSteps[i], res.Steps[i] = ladder[best], ladder[best]*scale
		res.Converged[i] = o.Fixed || bestU <= o.Tol*math.Abs(D[best])
	}
	return res, nil
}

// evalAt evaluates f with parameter i moved by step s and every other
// parameter at its nominal value. x is scratch space of length len(params).
func evalAt(p *Prescription, params []Param, f ScalarFunc, x []float64, i int, s float64) (float64, error) {
	for k := range x {
		x[k] = 0
	}
	x[i] = s
	q, err := applyParams(p, params, x)
	if err != nil {
		return 0, err
	}
	v, err := f(q)
	if err != nil {
		return 0, fmt.Errorf("lens: sensitivity: %s%+.3g: %w", params[i], s, err)
	}
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return 0, fmt.Errorf("lens: sensitivity: non-finite value at %s%+.3g", params[i], s)
	}
	return v, nil
}

// romberg extrapolates a sequence of central differences d[k] computed at the
// halving steps h/2^k. It returns the most-extrapolated value R[K-1][K-1] and
// the estimate |R[K-1][K-1] − R[K-1][K-2]| of its truncation error.
func romberg(d []float64) (ext, errEst float64) {
	K := len(d)
	R := make([][]float64, K)
	for k := range R {
		R[k] = make([]float64, k+1)
		R[k][0] = d[k]
		for j := 1; j <= k; j++ {
			w := math.Pow(4, float64(j)) - 1
			R[k][j] = R[k][j-1] + (R[k][j-1]-R[k-1][j-1])/w
		}
	}
	ext = R[K-1][K-1]
	if K >= 2 {
		errEst = math.Abs(R[K-1][K-1] - R[K-1][K-2])
	}
	return ext, errEst
}
