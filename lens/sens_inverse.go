package lens

import (
	"fmt"
	"math"
)

// AllocMode selects how the individual tolerance contributions combine into
// the merit increment of the linear model.
type AllocMode int

const (
	// AllocLinear is the worst-case (arithmetic) budget: ΔM = Σ |s_i|·t_i, all
	// parameters moved to their limits in the unfavourable direction.
	AllocLinear AllocMode = iota
	// AllocRSS is the statistical (root-sum-square) budget for independent
	// random tolerances: ΔM = sqrt(Σ (s_i·t_i)²).
	AllocRSS
)

// String implements fmt.Stringer.
func (m AllocMode) String() string {
	switch m {
	case AllocLinear:
		return "linear"
	case AllocRSS:
		return "rss"
	}
	return fmt.Sprintf("AllocMode(%d)", int(m))
}

// Allocation is the result of AllocateTolerances.
//
// Equal contribution: every parameter that is not at its limit contributes the
// same share of the budget, |s_i|·t_i = B/N (linear) or (s_i·t_i)² = B²/N
// (RSS). A parameter whose tolerance would exceed its limit is held at the
// limit and the remaining budget is shared among the others (water filling),
// so the budget is always fully used.
type Allocation struct {
	Mode   AllocMode
	Budget float64 // merit increment budget ΔM (mm)
	Params []Param
	Sens   []float64 // ∂M/∂x_i (input, mm per parameter unit)
	Tol    []float64 // allocated tolerance t_i (half-range in parameter units)
	Share  []float64 // budget share: |s_i|·t_i (linear) or (s_i·t_i)² (RSS)
	Capped []bool    // true when t_i sits at its limit (or has zero sensitivity)
	// Predicted is the linear-model merit increment of the allocated
	// tolerances: Σ|s_i|·t_i (linear) or sqrt(Σ(s_i·t_i)²) (RSS), recomputed
	// from Tol.
	Predicted float64
	// RelErr is |Predicted − Budget| / Budget.
	RelErr float64
}

// AllocateTolerances distributes a merit budget over the parameters in
// proportion to their sensitivities (inverse sensitivity). sens[i] is
// ∂M/∂x_i and params[i] names the parameter; maxTol (optional, same length)
// caps each tolerance, and a non-positive or +Inf entry means no cap.
//
// A parameter with zero sensitivity needs a cap: it consumes no budget and
// gets its cap as tolerance. The function errors when the budget cannot be
// placed, i.e. when every parameter is capped before the budget is used up.
func AllocateTolerances(params []Param, sens []float64, budget float64, mode AllocMode, maxTol []float64) (*Allocation, error) {
	n := len(params)
	if n == 0 {
		return nil, fmt.Errorf("lens: no parameters to allocate")
	}
	if len(sens) != n {
		return nil, fmt.Errorf("lens: %d sensitivities for %d parameters", len(sens), n)
	}
	if maxTol != nil && len(maxTol) != n {
		return nil, fmt.Errorf("lens: %d tolerance limits for %d parameters", len(maxTol), n)
	}
	if mode != AllocLinear && mode != AllocRSS {
		return nil, fmt.Errorf("lens: unknown allocation mode %d", int(mode))
	}
	if !(budget > 0) || math.IsInf(budget, 0) {
		return nil, fmt.Errorf("lens: budget must be positive and finite (%v)", budget)
	}

	a := &Allocation{
		Mode:   mode,
		Budget: budget,
		Params: append([]Param(nil), params...),
		Sens:   append([]float64(nil), sens...),
		Tol:    make([]float64, n),
		Share:  make([]float64, n),
		Capped: make([]bool, n),
	}
	// Work in share units: q_i = |s_i|·t_i (linear) or (s_i·t_i)² (RSS); the
	// budget is Q = B or B².
	Q := budget
	if mode == AllocRSS {
		Q = budget * budget
	}
	tmax := make([]float64, n)
	qcap := make([]float64, n)
	active := make([]bool, n)
	nActive := 0
	for i := range params {
		if math.IsNaN(sens[i]) || math.IsInf(sens[i], 0) {
			return nil, fmt.Errorf("lens: non-finite sensitivity for %s", params[i])
		}
		tmax[i] = math.Inf(1)
		if maxTol != nil && !math.IsInf(maxTol[i], 1) {
			if !(maxTol[i] > 0) {
				return nil, fmt.Errorf("lens: tolerance limit for %s must be positive", params[i])
			}
			tmax[i] = maxTol[i]
		}
		s := math.Abs(sens[i])
		if s == 0 {
			if math.IsInf(tmax[i], 1) {
				return nil, fmt.Errorf("lens: %s has zero sensitivity and no tolerance limit", params[i])
			}
			a.Tol[i] = tmax[i]
			a.Capped[i] = true
			continue
		}
		if mode == AllocLinear {
			qcap[i] = s * tmax[i]
		} else {
			qcap[i] = (s * tmax[i]) * (s * tmax[i])
		}
		active[i] = true
		nActive++
	}

	// Water filling: a parameter whose cap lies below the equal share is held at
	// its cap; the rest of the budget is re-shared among the others until no
	// further cap binds.
	q := make([]float64, n)
	rem := Q
	for nActive > 0 {
		share := rem / float64(nActive)
		capped := false
		for i := range params {
			if active[i] && qcap[i] <= share {
				q[i] = qcap[i]
				rem -= qcap[i]
				active[i] = false
				nActive--
				a.Capped[i] = true
				capped = true
			}
		}
		if !capped {
			for i := range params {
				if active[i] {
					q[i] = share
					active[i] = false
				}
			}
			nActive = 0
			rem = 0
		}
	}
	if rem > 1e-12*Q {
		return nil, fmt.Errorf("lens: budget %g cannot be placed: all parameters reach their limits", budget)
	}

	for i := range params {
		s := math.Abs(sens[i])
		if s == 0 {
			continue
		}
		if a.Capped[i] {
			a.Tol[i] = tmax[i]
		} else if mode == AllocLinear {
			a.Tol[i] = q[i] / s
		} else {
			a.Tol[i] = math.Sqrt(q[i]) / s
		}
		a.Share[i] = q[i]
	}
	var lin, sq float64
	for i := range params {
		c := math.Abs(sens[i]) * a.Tol[i]
		lin += c
		sq += c * c
	}
	if mode == AllocLinear {
		a.Predicted = lin
	} else {
		a.Predicted = math.Sqrt(sq)
	}
	a.RelErr = math.Abs(a.Predicted-budget) / budget
	return a, nil
}
