package coating

import (
	"errors"
	"fmt"
	"math"
	"math/rand"
	"strings"
)

// Quantity selects which radiometric quantity of a stack a Target constrains.
type Quantity uint8

const (
	// Reflectance constrains R (the reflected power fraction).
	Reflectance Quantity = iota
	// Transmittance constrains T.
	Transmittance
	// Absorptance constrains A = 1 - R - T.
	Absorptance
)

// String implements fmt.Stringer.
func (q Quantity) String() string {
	switch q {
	case Reflectance:
		return "R"
	case Transmittance:
		return "T"
	case Absorptance:
		return "A"
	default:
		return fmt.Sprintf("Quantity(%d)", uint8(q))
	}
}

// Target is one design goal: the quantity Qty of the stack at
// (LambdaNm, AngleDeg, Pol) should be as close as possible to Value, weighted
// by Weight in the merit function.
//
// Log switches to a logarithmic (relative) error,
//
//	residual = sqrt(Weight)·(log10(max(calc, 1e-12)) - log10(max(Value, 1e-12))),
//
// which is the right choice for low-reflectance or low-transmittance goals:
// it asks for a fixed *relative* error (R = 1e-4 and 2e-4 are one octave apart)
// instead of letting a single large value dominate the merit. Linear targets
// use residual = sqrt(Weight)·(calc - Value). The conventions for Weight are
// therefore "relative importance in the residual vector" (the merit function
// is the sum of squared residuals).
type Target struct {
	LambdaNm float64
	AngleDeg float64
	Pol      Polarization
	Qty      Quantity
	Value    float64
	Weight   float64
	Log      bool
}

// RTarget builds a reflectance target at normal incidence for unpolarized
// light, turning the logarithmic scale on automatically for low-reflectance
// goals (Value <= 0.05).
func RTarget(lambdaNm, value, weight float64) Target {
	return Target{LambdaNm: lambdaNm, Pol: PolUnpolarized, Qty: Reflectance, Value: value, Weight: weight, Log: value <= 0.05}
}

// TTarget builds a transmittance target at normal incidence for unpolarized
// light, with the logarithmic scale for low-transmittance goals (Value <= 0.05).
func TTarget(lambdaNm, value, weight float64) Target {
	return Target{LambdaNm: lambdaNm, Pol: PolUnpolarized, Qty: Transmittance, Value: value, Weight: weight, Log: value <= 0.05}
}

// ATarget builds an absorptance target at normal incidence for unpolarized
// light.
func ATarget(lambdaNm, value, weight float64) Target {
	return Target{LambdaNm: lambdaNm, Pol: PolUnpolarized, Qty: Absorptance, Value: value, Weight: weight}
}

// At returns a copy of the target evaluated at another angle and polarization.
func (t Target) At(angleDeg float64, pol Polarization) Target {
	t.AngleDeg, t.Pol = angleDeg, pol
	return t
}

// LogScale returns a copy of the target that uses the logarithmic residual.
func (t Target) LogScale() Target { t.Log = true; return t }

// LinearScale returns a copy of the target that uses the linear residual.
func (t Target) LinearScale() Target { t.Log = false; return t }

// String implements fmt.Stringer.
func (t Target) String() string {
	scale := "linear"
	if t.Log {
		scale = "log10"
	}
	return fmt.Sprintf("%s=%.6g at %.6g nm, %.4g deg, %s (weight %.4g, %s)",
		t.Qty, t.Value, t.LambdaNm, t.AngleDeg, t.Pol, t.Weight, scale)
}

// quantity returns the constrained quantity of the stack at the target's
// operating point.
func (t Target) quantity(s Stack) (float64, error) {
	smp, err := s.Evaluate(t.LambdaNm, t.AngleDeg, t.Pol)
	if err != nil {
		return 0, err
	}
	switch t.Qty {
	case Reflectance:
		return smp.R, nil
	case Transmittance:
		return smp.T, nil
	case Absorptance:
		return smp.A, nil
	default:
		return 0, fmt.Errorf("coating: invalid quantity %d in target", uint8(t.Qty))
	}
}

// logFloor keeps the logarithmic residual finite for a vanishing quantity.
const logFloor = 1e-12

// residual returns the weighted residual of one target for a stack.
func (t Target) residual(s Stack) (float64, error) {
	v, err := t.quantity(s)
	if err != nil {
		return 0, err
	}
	w := t.Weight
	if w == 0 {
		w = 1
	}
	if t.Log {
		d := math.Log10(math.Max(v, logFloor)) - math.Log10(math.Max(t.Value, logFloor))
		return math.Sqrt(w) * d, nil
	}
	return math.Sqrt(w) * (v - t.Value), nil
}

// Merit returns the weighted sum of squared target residuals of a stack: the
// value minimised by Optimize (a merit of zero is a perfect fit).
func Merit(stack Stack, targets []Target) (float64, error) {
	m := 0.0
	for _, t := range targets {
		r, err := t.residual(stack)
		if err != nil {
			return 0, err
		}
		m += r * r
	}
	return m, nil
}

// TargetValues returns the achieved value of every target for a stack, in the
// same order as the targets. Together with Merit it is what a design report
// prints.
func TargetValues(stack Stack, targets []Target) ([]float64, error) {
	out := make([]float64, len(targets))
	for i, t := range targets {
		v, err := t.quantity(stack)
		if err != nil {
			return nil, err
		}
		out[i] = v
	}
	return out, nil
}

// BandTargets builds one target per wavelength of the grid from fromNm to toNm
// (step stepNm, inclusive) for the given quantity, angle and polarization.
func BandTargets(fromNm, toNm, stepNm, angleDeg float64, pol Polarization, qty Quantity, value, weight float64, logScale bool) ([]Target, error) {
	lambdas, err := grid(fromNm, toNm, stepNm)
	if err != nil {
		return nil, err
	}
	out := make([]Target, len(lambdas))
	for i, l := range lambdas {
		out[i] = Target{LambdaNm: l, AngleDeg: angleDeg, Pol: pol, Qty: qty, Value: value, Weight: weight, Log: logScale}
	}
	return out, nil
}

// ---------------------------------------------------------------------------
// Options and results
// ---------------------------------------------------------------------------

// Options configures Optimize and GrowLayers. The zero value is usable: every
// field falls back to the documented default (a negative/zero value means
// "unset", except MinThicknessNm == 0, which really does allow layers to be
// driven to zero thickness).
type Options struct {
	// MaxIter is the number of damped least squares iterations per start
	// (default 200).
	MaxIter int
	// Tolerance stops the iteration when the relative merit improvement of an
	// accepted step drops below it (default 1e-9).
	Tolerance float64
	// Restarts is the number of additional random starts beyond the given
	// initial design (default 0). The best merit over all starts wins.
	Restarts int
	// Seed makes the restart sequence reproducible (default 1).
	Seed int64
	// RestartJitter is the log-normal sigma of the restart perturbation
	// (default 0.6: a typical thickness is perturbed by about +-60%).
	RestartJitter float64
	// MinThicknessNm and MaxThicknessNm bound every free thickness (defaults
	// 1 nm and 2000 nm). MinThicknessNm may be set to exactly 0 to let the
	// optimiser annihilate layers; a negative value selects the default.
	MinThicknessNm, MaxThicknessNm float64
	// Free is an optional per-layer mask: Free[i] == false freezes layer i.
	// A nil or short mask leaves the missing layers free.
	Free []bool

	// GrowthRounds is the number of needle-like layer-insertion rounds
	// performed by GrowLayers (default 0 = none).
	GrowthRounds int
	// GrowthThicknessNm is the probe thickness inserted by GrowLayers
	// (default 10 nm).
	GrowthThicknessNm float64
}

// DefaultOptions returns the standard optimiser configuration: 200 iterations,
// 1 nm to 2000 nm thickness bounds, no restarts and no layer growth.
func DefaultOptions() Options {
	return Options{
		MaxIter:           200,
		Tolerance:         1e-9,
		Restarts:          0,
		Seed:              1,
		RestartJitter:     0.6,
		MinThicknessNm:    1,
		MaxThicknessNm:    2000,
		GrowthRounds:      0,
		GrowthThicknessNm: 10,
	}
}

func (o Options) withDefaults() Options {
	d := DefaultOptions()
	if o.MaxIter <= 0 {
		o.MaxIter = d.MaxIter
	}
	if o.Tolerance <= 0 {
		o.Tolerance = d.Tolerance
	}
	if o.Restarts < 0 {
		o.Restarts = d.Restarts
	}
	if o.Seed == 0 {
		o.Seed = d.Seed
	}
	if o.RestartJitter <= 0 {
		o.RestartJitter = d.RestartJitter
	}
	if o.MinThicknessNm < 0 {
		o.MinThicknessNm = d.MinThicknessNm
	}
	if o.MaxThicknessNm <= 0 {
		o.MaxThicknessNm = d.MaxThicknessNm
	}
	if o.MaxThicknessNm < o.MinThicknessNm {
		o.MaxThicknessNm = o.MinThicknessNm
	}
	if o.GrowthRounds < 0 {
		o.GrowthRounds = d.GrowthRounds
	}
	if o.GrowthThicknessNm <= 0 {
		o.GrowthThicknessNm = d.GrowthThicknessNm
	}
	return o
}

// Result is the outcome of an optimisation run.
type Result struct {
	// Stack is the optimised design (thicknesses rounded only by arithmetic).
	Stack Stack
	// Targets are the goals the design was fitted to.
	Targets []Target
	// Merit is the final weighted sum of squared residuals.
	Merit float64
	// RmsResidual is sqrt(Merit/len(Targets)).
	RmsResidual float64
	// InitialMerit is the merit of the starting design of the winning start.
	InitialMerit float64
	// Iterations counts all damped least squares iterations over all starts.
	Iterations int
	// Starts is the number of designs tried (1 + Options.Restarts).
	Starts int
	// History is the merit after every accepted step of the winning start
	// (History[0] is the merit of its starting design).
	History []float64
}

// Achieved returns the value of every target for the optimised stack.
func (r *Result) Achieved() ([]float64, error) { return TargetValues(r.Stack, r.Targets) }

// Improvement returns how much better the result is than its starting point,
// as InitialMerit/Merit (1 means "no improvement").
func (r *Result) Improvement() float64 {
	if r.Merit <= 0 {
		return math.Inf(1)
	}
	return r.InitialMerit / r.Merit
}

// Summary renders the design and its target residuals as a multi-line report.
func (r *Result) Summary() string {
	var b strings.Builder
	fmt.Fprintf(&b, "coating design: %d layers, %.1f nm total, merit %.6g (%d starts, %d iterations)\n",
		len(r.Stack.Layers), r.Stack.TotalThicknessNm(), r.Merit, r.Starts, r.Iterations)
	for i, l := range r.Stack.Layers {
		fmt.Fprintf(&b, "  layer %2d: %-10s %8.2f nm\n", i, l.Material.Name, l.ThicknessNm)
	}
	values, err := r.Achieved()
	for i, t := range r.Targets {
		if err != nil {
			fmt.Fprintf(&b, "  target %d: %s (unevaluated: %v)\n", i, t, err)
			continue
		}
		fmt.Fprintf(&b, "  target %2d: %-42s achieved %.6g\n", i, t.String(), values[i])
	}
	return b.String()
}

// ---------------------------------------------------------------------------
// Damped least squares core
// ---------------------------------------------------------------------------

// residualFunc fills out with the residual vector of the design x.
type residualFunc func(x, out []float64) error

// jacobianStepNm and jacobianRelStep control the finite-difference step of the
// numerically evaluated Jacobian: h_j = max(jacobianStepNm, jacobianRelStep·x_j).
// The thickness dependence of a thin-film spectrum is smooth and the smallest
// phase step introduced is 2*pi*h/lambda ~ 1e-5 rad, far above double rounding.
const (
	jacobianStepNm  = 1e-3
	jacobianRelStep = 1e-4
)

// sumSquares returns the sum of squares of v.
func sumSquares(v []float64) float64 {
	s := 0.0
	for _, x := range v {
		s += x * x
	}
	return s
}

// solveLinear solves M·x = b in place (Gaussian elimination with partial
// pivoting) and reports whether the matrix was numerically invertible. Both M
// and b are destroyed.
func solveLinear(M [][]float64, b []float64) bool {
	n := len(b)
	for k := 0; k < n; k++ {
		p, best := k, math.Abs(M[k][k])
		for i := k + 1; i < n; i++ {
			if v := math.Abs(M[i][k]); v > best {
				best, p = v, i
			}
		}
		if best == 0 || math.IsNaN(best) {
			return false
		}
		if p != k {
			M[k], M[p] = M[p], M[k]
			b[k], b[p] = b[p], b[k]
		}
		piv := M[k][k]
		for i := k + 1; i < n; i++ {
			f := M[i][k] / piv
			if f == 0 {
				continue
			}
			for j := k; j < n; j++ {
				M[i][j] -= f * M[k][j]
			}
			b[i] -= f * b[k]
		}
	}
	for i := n - 1; i >= 0; i-- {
		s := b[i]
		for j := i + 1; j < n; j++ {
			s -= M[i][j] * b[j]
		}
		b[i] = s / M[i][i]
	}
	for _, v := range b {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return false
		}
	}
	return true
}

// dls minimises the sum of squared residuals with a Levenberg-Marquardt style
// damped least squares iteration:
//
//  1. the residual vector r (m entries) and its Jacobian J (central finite
//     differences, one-sided where a variable sits on its bound) are evaluated
//     at x;
//  2. the step solves (JᵗJ + mu·diag(JᵗJ))·delta = -Jᵗr;
//  3. a backtracking line search accepts the step only if it lowers the merit,
//     otherwise mu is increased (more damping) and the step is recomputed.
//
// Thickness bounds are enforced by projecting every trial point onto the box
// [lo, hi]; the iteration ends when no step improves the merit, when the
// relative improvement falls below tol, or after maxIter iterations.
func dls(x0, lo, hi []float64, m int, residue residualFunc, maxIter int, tol float64) (x []float64, merit float64, iters int, history []float64, err error) {
	n := len(x0)
	if n == 0 {
		return nil, 0, 0, nil, errors.New("coating: no free variables")
	}
	if m <= 0 {
		return nil, 0, 0, nil, errors.New("coating: empty residual vector")
	}
	probe := make([]float64, n)

	x = append([]float64(nil), x0...)
	r := make([]float64, m)
	if err := residue(x, r); err != nil {
		return nil, 0, 0, nil, err
	}
	merit = sumSquares(r)
	history = append(history, merit)

	rn := make([]float64, m)
	rp := make([]float64, m)
	rm := make([]float64, m)
	J := make([][]float64, m)
	for i := range J {
		J[i] = make([]float64, n)
	}
	A := make([][]float64, n)
	M := make([][]float64, n)
	for i := range A {
		A[i] = make([]float64, n)
		M[i] = make([]float64, n)
	}
	g := make([]float64, n)
	delta := make([]float64, n)
	xNew := make([]float64, n)

	mu := 1e-3
	converged := false
	for iter := 0; iter < maxIter && !converged; iter++ {
		iters = iter + 1

		// --- Jacobian by finite differences ------------------------------
		for j := 0; j < n; j++ {
			h := math.Max(jacobianStepNm, jacobianRelStep*math.Abs(x[j]))
			hp, hm := h, h
			if x[j]+hp > hi[j] {
				hp = 0
			}
			if x[j]-hm < lo[j] {
				hm = 0
			}
			if hp == 0 && hm == 0 {
				for i := range J {
					J[i][j] = 0
				}
				continue
			}
			copy(probe, x)
			if hp > 0 {
				probe[j] = x[j] + hp
				if err := residue(probe, rp); err != nil {
					return nil, 0, 0, nil, err
				}
			}
			if hm > 0 {
				probe[j] = x[j] - hm
				if err := residue(probe, rm); err != nil {
					return nil, 0, 0, nil, err
				}
			}
			for i := 0; i < m; i++ {
				switch {
				case hp == 0:
					J[i][j] = (r[i] - rm[i]) / hm
				case hm == 0:
					J[i][j] = (rp[i] - r[i]) / hp
				default:
					J[i][j] = (rp[i] - rm[i]) / (hp + hm)
				}
			}
		}

		// --- normal equations --------------------------------------------
		for i := range A {
			for j := range A[i] {
				A[i][j] = 0
			}
		}
		for j := range g {
			g[j] = 0
		}
		for i := 0; i < m; i++ {
			Ji := J[i]
			ri := r[i]
			for a := 0; a < n; a++ {
				if Ji[a] == 0 {
					continue
				}
				g[a] += Ji[a] * ri
				for b := a; b < n; b++ {
					A[a][b] += Ji[a] * Ji[b]
				}
			}
		}
		for a := 0; a < n; a++ {
			for b := 0; b < a; b++ {
				A[a][b] = A[b][a]
			}
		}
		grad := 0.0
		for _, v := range g {
			grad = math.Max(grad, math.Abs(v))
		}
		if grad < 1e-14 {
			break // at a stationary point: nothing left to gain
		}

		// --- damped step with backtracking --------------------------------
		stepped := false
		for attempt := 0; attempt < 24 && !stepped; attempt++ {
			for a := 0; a < n; a++ {
				d := A[a][a]
				if d < 1e-14 {
					d = 1e-14
				}
				for b := 0; b < n; b++ {
					M[a][b] = A[a][b]
					if a == b {
						M[a][b] += mu * d
					}
				}
				delta[a] = -g[a]
			}
			if !solveLinear(M, delta) {
				mu *= 10
				continue
			}
			alpha := 1.0
			for ls := 0; ls < 12; ls++ {
				for a := 0; a < n; a++ {
					v := x[a] + alpha*delta[a]
					if v < lo[a] {
						v = lo[a]
					}
					if v > hi[a] {
						v = hi[a]
					}
					xNew[a] = v
				}
				if err := residue(xNew, rn); err != nil {
					return nil, 0, 0, nil, err
				}
				mn := sumSquares(rn)
				if mn < merit {
					improvement := merit - mn
					copy(x, xNew)
					copy(r, rn)
					merit = mn
					history = append(history, merit)
					if improvement <= tol*merit+1e-300 {
						converged = true
					}
					mu = math.Max(mu/3, 1e-12)
					stepped = true
					break
				}
				alpha *= 0.5
			}
			if !stepped {
				mu *= 10
				if mu > 1e10 {
					break
				}
			}
		}
		if !stepped {
			break // no damped step improves the merit: local minimum
		}
		if merit == 0 {
			break
		}
	}
	return x, merit, iters, history, nil
}

// problem binds a stack, its free layer indices and the targets into the
// residual function used by the solver. The scratch layer slice makes the hot
// path allocation free.
type problem struct {
	stack   Stack
	layers  []Layer
	free    []int
	targets []Target
}

func newProblem(base Stack, targets []Target, free []int) *problem {
	p := &problem{stack: base, targets: targets, free: free}
	p.layers = append([]Layer(nil), base.Layers...)
	p.stack.Layers = p.layers
	return p
}

// residue implements residualFunc. out must be exactly len(targets) long.
func (p *problem) residue(x, out []float64) error {
	if len(out) != len(p.targets) {
		return fmt.Errorf("coating: internal residual length mismatch: %d vs %d", len(out), len(p.targets))
	}
	for k, idx := range p.free {
		p.layers[idx].ThicknessNm = x[k]
	}
	for i, t := range p.targets {
		r, err := t.residual(p.stack)
		if err != nil {
			return err
		}
		out[i] = r
	}
	return nil
}

// freeIndices expands the Free mask into the list of free layer indices.
func freeIndices(stack Stack, mask []bool) []int {
	var idx []int
	for i := range stack.Layers {
		if i < len(mask) && !mask[i] {
			continue
		}
		idx = append(idx, i)
	}
	return idx
}

// Optimize refines the layer thicknesses of a stack so that the computed
// quantities approach the targets, using damped least squares with thickness
// bounds, an optional set of random restarts and a deterministic RNG.
//
// Only layer thicknesses are varied; the materials and the layer count are
// fixed (use GrowLayers to let the layer count grow). The returned Result
// carries the optimised stack, the merit and the iteration history. Note that
// the residual vector is built from finite-difference sensitivities, so the
// returned thicknesses are optimal to within the numerical accuracy of that
// Jacobian, not to machine precision.
func Optimize(stack Stack, targets []Target, opts Options) (*Result, error) {
	opts = opts.withDefaults()
	if err := stack.Validate(); err != nil {
		return nil, err
	}
	if len(stack.Layers) == 0 {
		return nil, errors.New("coating: Optimize needs at least one layer")
	}
	if len(targets) == 0 {
		return nil, errors.New("coating: Optimize needs at least one target")
	}
	free := freeIndices(stack, opts.Free)
	if len(free) == 0 {
		return nil, errors.New("coating: no free layers to optimise")
	}
	base := stack.ThicknessesNm()
	n := len(free)
	lo := make([]float64, n)
	hi := make([]float64, n)
	x0 := make([]float64, n)
	for k, idx := range free {
		lo[k] = opts.MinThicknessNm
		hi[k] = opts.MaxThicknessNm
		x0[k] = math.Min(math.Max(base[idx], lo[k]), hi[k])
	}
	prob := newProblem(stack, targets, free)

	res := &Result{Targets: targets, Starts: 1}
	xBest, meritBest, iters, hist, err := dls(x0, lo, hi, len(targets), prob.residue, opts.MaxIter, opts.Tolerance)
	if err != nil {
		return nil, err
	}
	initial := math.NaN()
	if len(hist) > 0 {
		initial = hist[0]
	}
	res.Iterations = iters
	res.History = hist
	res.Merit = meritBest
	res.InitialMerit = initial

	rng := rand.New(rand.NewSource(opts.Seed))
	for r := 1; r <= opts.Restarts; r++ {
		xs := perturbStart(rng, x0, lo, hi, opts.RestartJitter)
		xr, merr, it, h, err := dls(xs, lo, hi, len(targets), prob.residue, opts.MaxIter, opts.Tolerance)
		if err != nil {
			return nil, err
		}
		res.Starts++
		res.Iterations += it
		if merr < res.Merit {
			xBest, meritBest = xr, merr
			res.Merit, res.History = merr, h
			if len(h) > 0 {
				res.InitialMerit = h[0]
			}
		}
	}
	design := stack.WithThicknesses(base)
	for k, idx := range free {
		design.Layers[idx].ThicknessNm = xBest[k]
	}
	res.Stack = design
	if nt := len(targets); nt > 0 {
		res.RmsResidual = math.Sqrt(math.Max(res.Merit, 0) / float64(nt))
	}
	return res, nil
}

// perturbStart draws a random restart point: each variable is either replaced
// by a uniform draw over its bounds (25% of the cases, to explore a different
// basin) or multiplied by a log-normal factor, and the result is projected
// onto [lo, hi].
func perturbStart(rng *rand.Rand, x0, lo, hi []float64, jitter float64) []float64 {
	out := make([]float64, len(x0))
	for i := range x0 {
		if rng.Float64() < 0.25 {
			out[i] = lo[i] + rng.Float64()*(hi[i]-lo[i])
			continue
		}
		v := x0[i] * math.Exp(jitter*rng.NormFloat64())
		if v < lo[i] {
			v = lo[i]
		}
		if v > hi[i] {
			v = hi[i]
		}
		out[i] = v
	}
	return out
}

// insertLayer returns a copy of the stack with one layer inserted at position
// pos (0 = front surface, away from the substrate).
func insertLayer(s Stack, pos int, l Layer) Stack {
	out := s.Clone()
	out.Layers = append(out.Layers, Layer{})
	copy(out.Layers[pos+1:], out.Layers[pos:])
	out.Layers[pos] = l
	return out
}

// defaultGrowthMaterials are the materials tried by GrowLayers when the caller
// does not name any: a low and a high index oxide, the classic pair for
// dielectric multilayer design.
var defaultGrowthMaterials = []Material{SiO2, Ta2O5, TiO2, MgF2}

// growthEpsilon is the relative merit improvement a growth round must reach to
// be accepted.
const growthEpsilon = 1e-6

// GrowLayers optimises a stack and then lets the layer count grow: after every
// round it probes the insertion of a thin layer of each candidate material at
// every interface, refines the thickness of the inserted layer alone (a one
// variable damped least squares fit), keeps the insertion with the best merit
// and finally re-optimises all thicknesses together. This is the finite
// difference surrogate of the classical needle (adjoint) synthesis: the exact
// adjoint sensitivity is not implemented, the probe is explicit instead.
//
// Rounds stop when no insertion improves the merit by more than a relative
// 1e-6 or when the joint re-optimisation does not beat the previous design, so
// the result never has more layers than the targets can pay for.
func GrowLayers(stack Stack, targets []Target, materials []Material, opts Options) (*Result, error) {
	opts = opts.withDefaults()
	if len(materials) == 0 {
		materials = defaultGrowthMaterials
	}
	for _, m := range materials {
		if !m.Valid() {
			return nil, fmt.Errorf("coating: growth material %q is not valid", m.Name)
		}
	}
	best, err := Optimize(stack, targets, opts)
	if err != nil {
		return nil, err
	}
	for round := 0; round < opts.GrowthRounds; round++ {
		var (
			candStack Stack
			candMerit = math.Inf(1)
			candThick float64
			found     bool
		)
		nl := len(best.Stack.Layers)
		for pos := 0; pos <= nl; pos++ {
			for _, m := range materials {
				probe := insertLayer(best.Stack, pos, Layer{Material: m, ThicknessNm: opts.GrowthThicknessNm})
				sub := opts
				sub.MaxIter = 60
				sub.Restarts = 0
				sub.MinThicknessNm = 0 // the probe layer may vanish again
				sub.Free = make([]bool, len(probe.Layers))
				for i := range sub.Free {
					sub.Free[i] = i == pos
				}
				res, err := Optimize(probe, targets, sub)
				if err != nil {
					continue
				}
				if res.Merit < candMerit {
					candMerit = res.Merit
					candStack = res.Stack
					candThick = res.Stack.Layers[pos].ThicknessNm
					found = true
				}
			}
		}
		if !found || candMerit >= best.Merit*(1-growthEpsilon) || candThick <= 0 {
			break
		}
		// Joint re-optimisation of the grown stack.
		joint, err := Optimize(candStack, targets, opts)
		if err != nil {
			return nil, err
		}
		joint.Starts += best.Starts
		joint.Iterations += best.Iterations
		if joint.Merit >= best.Merit {
			break
		}
		best = joint
	}
	return best, nil
}
