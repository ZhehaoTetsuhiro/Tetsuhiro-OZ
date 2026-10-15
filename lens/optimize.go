// Package lens implements lens-system design on top of the geo core of
// Tetsuhiro-OZ: a sequential prescription model, first-order solves, merit
// functions built from real-ray analysis, and a damped-least-squares
// optimizer that turns design goals ("spot below X mrad, EFL = Y mm,
// chromatic focal shift below Z") into concrete radii/thicknesses.
//
// The package is pure Go standard library, float64 throughout, and follows
// the conventions of package geo (mm, nm, x right / y up / z along the
// optical axis, R > 0 puts the centre of curvature on the +z side).
package lens

import (
	"fmt"
	"math"
)

// Variable is one optimizable quantity (a radius, a thickness, ...).
//
// Lo and Hi are hard bounds: the optimizer clamps every trial value into
// [Lo, Hi] and never evaluates outside. Lo == Hi freezes the variable. Step
// is the finite-difference step used for the Jacobian; Step <= 0 asks for
// the default, 1e-6 * max(1, |Value|).
type Variable struct {
	Name  string
	Value float64
	Lo    float64
	Hi    float64
	Step  float64
}

// clampToBounds returns v clamped into [Lo, Hi].
func (v Variable) clamp(x float64) float64 {
	if x < v.Lo {
		return v.Lo
	}
	if x > v.Hi {
		return v.Hi
	}
	return x
}

// frozen reports whether the variable has no freedom (Lo == Hi).
func (v Variable) frozen() bool { return v.Lo == v.Hi }

// fmtVars formats a variable vector for messages.
func fmtVars(x []float64) string {
	s := "["
	for i, v := range x {
		if i > 0 {
			s += " "
		}
		s += fmt.Sprintf("%.6g", v)
	}
	return s + "]"
}

// Problem is a least-squares problem: a set of bounded variables and a
// residual function. Minimize looks for the variable vector that minimizes
// the sum of squared residuals
//
//	S(x) = Σ_k r_k(x)² .
//
// Residuals carry their own scaling and weights: callers divide by target
// tolerances before returning them, so that equal residual values mean equal
// design importance (a standard lens-design practice).
type Problem struct {
	Variables []Variable
	// Residual evaluates the residual vector at x. len(x) == len(Variables).
	// Returning an error aborts the optimization with that error.
	Residual func(x []float64) ([]float64, error)

	// MaxIter is the number of accepted-state iterations (default 100).
	MaxIter int
	// Lam0 is the initial Marquardt damping (default 1e-3).
	Lam0 float64
	// TolRel stops the run when an accepted iteration reduces the merit by
	// less than TolRel relative (default 1e-12).
	TolRel float64
	// Central selects central differences (2n evaluations per Jacobian)
	// instead of forward differences (n+1).
	Central bool
	// Trace, when non-nil, is called after every accepted iteration with
	// (iteration, merit, x-copy) — handy for CLI progress output.
	Trace func(iter int, merit float64, x []float64)
}

// Result describes the outcome of Minimize.
type Result struct {
	X         []float64 // best variable vector found (a copy)
	Merit     float64   // S(X)
	RMS       float64   // sqrt(S(X) / m) over the last residual vector
	Iter      int       // accepted iterations performed
	Evals     int       // residual-function evaluations
	Converged bool      // TolRel reached within MaxIter
	History   []float64 // merit after each accepted iteration (incl. start)
}

// PrescriptionProblem builds a least-squares Problem over design variables of
// a prescription. specs names the quantities to vary — any VarKind, including
// the polynomial-asphere term coefficients (VarPolyTerm) and the Forbes Q
// coefficients (VarQCoef) — and bounds returns the [lo, hi] of each variable
// from its spec and current value (a deterministic, kind-aware rule belongs to
// the caller: radii, thicknesses and coefficients have different natural
// ranges). merit is the residual function of a *prescription* clone, e.g.
//
//	func(q *lens.Prescription) ([]float64, error) {
//		return lens.SpotResiduals(q, geo.Wavelengthd, fields, samples, 0.01)
//	}
//
// The residual wrapper clones the prescription for every evaluation, so a run
// never mutates p and the residuals are a pure function of the variable
// vector: the same x always gives the same value, and Minimize is therefore
// bit-for-bit reproducible. All specs are validated (via Vars) before the
// optimizer starts, so a mis-specified variable is an error instead of a
// silently constant residual.
//
// The variable vector follows specs order, and Variable.Name is the VarSpec
// label ("s1.poly0", "s1.q1", …), so traces and reports identify the asphere
// coefficients.
func PrescriptionProblem(p *Prescription, specs []VarSpec, bounds func(VarSpec, float64) (lo, hi float64), merit func(*Prescription) ([]float64, error)) (Problem, error) {
	if p == nil {
		return Problem{}, fmt.Errorf("lens: nil prescription")
	}
	if len(specs) == 0 {
		return Problem{}, fmt.Errorf("lens: no variables")
	}
	if merit == nil {
		return Problem{}, fmt.Errorf("lens: nil merit function")
	}
	cur, err := p.Vars(specs)
	if err != nil {
		return Problem{}, err
	}
	vars := make([]Variable, len(specs))
	for i, s := range specs {
		lo, hi := bounds(s, cur[i])
		if lo > hi {
			return Problem{}, fmt.Errorf("lens: variable %s: empty bounds [%g, %g]", s, lo, hi)
		}
		vars[i] = Variable{Name: s.String(), Value: cur[i], Lo: lo, Hi: hi}
	}
	return Problem{
		Variables: vars,
		Residual: func(x []float64) ([]float64, error) {
			q := p.Clone()
			if err := q.SetVars(specs, x); err != nil {
				return nil, err
			}
			return merit(q)
		},
	}, nil
}

// Minimize solves p by damped least squares (Levenberg–Marquardt with
// Marquardt's diagonal scaling). It is deterministic: identical inputs give
// bit-identical outputs (no randomness, fixed tie-breaks).
//
// Algorithm: starting from the current variable values,
//
//  1. evaluate the residual vector r and merit S = Σr²;
//  2. build the Jacobian J by finite differences inside the bounds;
//  3. solve (JᵀJ + λ·diag(JᵀJ)) δ = −Jᵀr  (λ → diag entries replaced by λ when
//     a diagonal entry is zero, e.g. frozen variables);
//  4. trial x' = clamp(x + δ); if S(x') < S, accept and set λ /= 3, else
//     reject and set λ *= 5, retrying the solve (at most 20 attempts per
//     iteration, λ capped at 1e12);
//  5. stop when an accepted step reduces S by less than TolRel·S or a step
//     length of zero coincides with the rejection limit.
func Minimize(p Problem) (*Result, error) {
	if len(p.Variables) == 0 {
		return nil, fmt.Errorf("lens: no variables")
	}
	if p.Residual == nil {
		return nil, fmt.Errorf("lens: nil residual function")
	}
	maxIter := p.MaxIter
	if maxIter <= 0 {
		maxIter = 100
	}
	lam := p.Lam0
	if !(lam > 0) {
		lam = 1e-3
	}
	tolRel := p.TolRel
	if !(tolRel > 0) {
		tolRel = 1e-12
	}

	n := len(p.Variables)
	x := make([]float64, n)
	for i, v := range p.Variables {
		x[i] = v.clamp(v.Value)
	}
	eval := 0
	res, err := p.Residual(x)
	if err != nil {
		return nil, err
	}
	eval++
	if len(res) == 0 {
		return nil, fmt.Errorf("lens: residual vector is empty")
	}
	m := len(res)

	merit := 0.0
	for _, r := range res {
		merit += r * r
	}

	out := &Result{Iter: 0, Evals: 0, History: []float64{merit}}
	if p.Trace != nil {
		p.Trace(0, merit, append([]float64(nil), x...))
	}

	// Scratch space.
	j := make([][]float64, m)
	for k := range j {
		j[k] = make([]float64, n)
	}
	jtj := make([]float64, n*n)
	jtjDamp := make([]float64, n*n)
	jtr := make([]float64, n)
	delta := make([]float64, n)
	xtrial := make([]float64, n)
	rTrial := make([]float64, 0, m)

	for iter := 1; iter <= maxIter; iter++ {
		// ---- Jacobian (finite differences, respecting bounds) ----
		for i, v := range p.Variables {
			if v.frozen() {
				for k := 0; k < m; k++ {
					j[k][i] = 0
				}
				continue
			}
			h := v.Step
			if h <= 0 {
				h = 1e-6 * math.Max(1, math.Abs(x[i]))
			}
			// Keep both probe points inside the bounds, with a tiny inset.
			lo, hi := v.Lo, v.Hi
			xi := x[i]
			var xp, xm float64
			if p.Central {
				xp, xm = xi+h, xi-h
				if xp > hi || xm < lo {
					// Shift the stencil to fit.
					if xi-h < lo {
						xm = lo
						xp = math.Min(hi, lo+2*h)
					} else {
						xp = hi
						xm = math.Max(lo, hi-2*h)
					}
				}
				if xp == xm {
					xp = math.Min(hi, xi+h)
					xm = math.Max(lo, xi-h)
				}
			} else {
				xp = math.Min(hi, xi+h)
				if xp == xi { // at the upper bound: probe backwards
					xp = math.Max(lo, xi-h)
				}
				xm = xi
			}
			// Evaluate at xp.
			copy(xtrial, x)
			xtrial[i] = xp
			rp, err := p.Residual(xtrial)
			if err != nil {
				return nil, err
			}
			eval++
			if len(rp) != m {
				return nil, fmt.Errorf("lens: residual length changed (%d -> %d)", m, len(rp))
			}
			if p.Central {
				copy(xtrial, x)
				xtrial[i] = xm
				rm, err := p.Residual(xtrial)
				if err != nil {
					return nil, err
				}
				eval++
				if len(rm) != m {
					return nil, fmt.Errorf("lens: residual length changed (%d -> %d)", m, len(rm))
				}
				hEff := xp - xm
				if hEff == 0 {
					for k := 0; k < m; k++ {
						j[k][i] = 0
					}
					continue
				}
				for k := 0; k < m; k++ {
					j[k][i] = (rp[k] - rm[k]) / hEff
				}
			} else {
				hEff := xp - xm
				if hEff == 0 {
					for k := 0; k < m; k++ {
						j[k][i] = 0
					}
					continue
				}
				for k := 0; k < m; k++ {
					j[k][i] = (rp[k] - res[k]) / hEff
				}
			}
		}

		// ---- Normal equations ----
		for a := 0; a < n; a++ {
			for b := a; b < n; b++ {
				s := 0.0
				for k := 0; k < m; k++ {
					s += j[k][a] * j[k][b]
				}
				jtj[a*n+b] = s
				jtj[b*n+a] = s
			}
			s := 0.0
			for k := 0; k < m; k++ {
				s += j[k][a] * res[k]
			}
			jtr[a] = s
		}

		// ---- Damped solve with rejection loop ----
		accepted := false
		var meritNew float64
		for attempt := 0; attempt < 20; attempt++ {
			copy(jtjDamp, jtj)
			for a := 0; a < n; a++ {
				d := jtj[a*n+a]
				if d <= 0 {
					d = 1
				}
				jtjDamp[a*n+a] = d + lam*d
			}
			// Solve (JᵀJ + λD) δ = −Jᵀr via Cholesky.
			copy(delta, jtr)
			for a := range delta {
				delta[a] = -delta[a]
			}
			if err := cholSolve(jtjDamp, delta, n); err != nil {
				lam *= 5
				if lam > 1e12 {
					break
				}
				continue
			}
			moved := false
			for a := 0; a < n; a++ {
				xtrial[a] = p.Variables[a].clamp(x[a] + delta[a])
				if xtrial[a] != x[a] {
					moved = true
				}
			}
			if !moved {
				lam *= 5
				if lam > 1e12 {
					break
				}
				continue
			}
			rTrial = rTrial[:0]
			rTrial, err = p.Residual(xtrial)
			if err != nil {
				return nil, err
			}
			eval++
			if len(rTrial) != m {
				return nil, fmt.Errorf("lens: residual length changed (%d -> %d)", m, len(rTrial))
			}
			meritNew = 0
			for _, r := range rTrial {
				meritNew += r * r
			}
			if meritNew < merit {
				accepted = true
				break
			}
			lam *= 5
			if lam > 1e12 {
				break
			}
		}

		if !accepted {
			// No improvement: the damped model finds no better point, the
			// classic local-minimum stop.
			out.Converged = true
			break
		}

		// Accept: commit state.
		copy(x, xtrial)
		copy(res, rTrial)
		old := merit
		merit = meritNew
		lam /= 3
		if lam < 1e-12 {
			lam = 1e-12
		}
		out.Iter = iter
		out.History = append(out.History, merit)
		if p.Trace != nil {
			p.Trace(iter, merit, append([]float64(nil), x...))
		}
		if old-merit <= tolRel*math.Max(old, 1e-300) && iter > 2 {
			out.Converged = true
			out.X = append([]float64(nil), x...)
			out.Merit = merit
			out.RMS = math.Sqrt(merit / float64(m))
			out.Evals = eval
			return out, nil
		}
	}

	out.X = append([]float64(nil), x...)
	out.Merit = merit
	out.RMS = math.Sqrt(merit / float64(m))
	out.Evals = eval
	return out, nil
}

// cholSolve solves A·x = b in place for a symmetric positive-definite A
// stored row-major (n×n); b is overwritten with the solution. Returns an
// error when A is not positive definite (the caller then increases λ).
func cholSolve(a []float64, b []float64, n int) error {
	// In-place Cholesky: a = L·Lᵀ (lower triangle).
	for i := 0; i < n; i++ {
		for jj := 0; jj <= i; jj++ {
			s := a[i*n+jj]
			for k := 0; k < jj; k++ {
				s -= a[i*n+k] * a[jj*n+k]
			}
			if i == jj {
				if !(s > 0) {
					return fmt.Errorf("lens: matrix not positive definite")
				}
				a[i*n+i] = math.Sqrt(s)
			} else {
				a[i*n+jj] = s / a[jj*n+jj]
			}
		}
	}
	// Forward substitution: L·y = b.
	for i := 0; i < n; i++ {
		s := b[i]
		for k := 0; k < i; k++ {
			s -= a[i*n+k] * b[k]
		}
		b[i] = s / a[i*n+i]
	}
	// Back substitution: Lᵀ·x = y.
	for i := n - 1; i >= 0; i-- {
		s := b[i]
		for k := i + 1; k < n; k++ {
			s -= a[k*n+i] * b[k]
		}
		b[i] = s / a[i*n+i]
	}
	return nil
}
