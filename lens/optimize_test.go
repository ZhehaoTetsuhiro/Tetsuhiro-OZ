package lens

import (
	"math"
	"testing"
)

// A pure quadratic residual: S(x,y) = (x-3)² + 2(y+1)², minimum at (3,-1).
func parabolaResidual(x []float64) ([]float64, error) {
	return []float64{x[0] - 3, math.Sqrt2 * (x[1] + 1)}, nil
}

func TestMinimizeParabola(t *testing.T) {
	p := Problem{
		Variables: []Variable{
			{Name: "x", Value: 0, Lo: -100, Hi: 100},
			{Name: "y", Value: 0, Lo: -100, Hi: 100},
		},
		Residual: parabolaResidual,
	}
	res, err := Minimize(p)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(res.X[0]-3) > 1e-8 || math.Abs(res.X[1]+1) > 1e-8 {
		t.Fatalf("parabola: got %v, want [3 -1]", res.X)
	}
	if res.Merit > 1e-16 {
		t.Fatalf("parabola: merit %g, want ~0", res.Merit)
	}
	if res.Iter < 1 || res.Evals < 1 {
		t.Fatalf("parabola: degenerate counters %+v", res)
	}
}

func TestMinimizeExponential(t *testing.T) {
	const aStar, bStar = 2.5, 0.7
	ts := make([]float64, 40)
	ys := make([]float64, 40)
	for i := range ts {
		ts[i] = 0.05 * float64(i+1)
		ys[i] = aStar * math.Exp(bStar*ts[i])
	}
	resid := func(x []float64) ([]float64, error) {
		out := make([]float64, len(ts))
		for i, tv := range ts {
			out[i] = x[0]*math.Exp(x[1]*tv) - ys[i]
		}
		return out, nil
	}
	p := Problem{
		Variables: []Variable{
			{Name: "a", Value: 1, Lo: 1e-9, Hi: 100},
			{Name: "b", Value: 1, Lo: -5, Hi: 5},
		},
		Residual: resid,
		MaxIter:  200,
	}
	got, err := Minimize(p)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(got.X[0]-aStar) > 1e-6 || math.Abs(got.X[1]-bStar) > 1e-6 {
		t.Fatalf("exponential: got a=%g b=%g, want %g %g", got.X[0], got.X[1], aStar, bStar)
	}
	if !got.Converged {
		t.Fatalf("exponential: not converged after %d iterations (merit %g)", got.Iter, got.Merit)
	}
	if got.RMS > 1e-10 {
		t.Fatalf("exponential: RMS %g too large", got.RMS)
	}
}

func TestMinimizeRespectsBounds(t *testing.T) {
	p := Problem{
		Variables: []Variable{
			{Name: "x", Value: 0, Lo: -100, Hi: 100},
			{Name: "y", Value: 0, Lo: -0.5, Hi: 0.5}, // optimum (-1) is out of bounds
		},
		Residual: parabolaResidual,
	}
	res, err := Minimize(p)
	if err != nil {
		t.Fatal(err)
	}
	if res.X[1] != -0.5 {
		t.Fatalf("bounds: y = %g, want clamped to -0.5", res.X[1])
	}
	if math.Abs(res.X[0]-3) > 1e-8 {
		t.Fatalf("bounds: x = %g, want 3", res.X[0])
	}
}

func TestMinimizeFrozenVariable(t *testing.T) {
	p := Problem{
		Variables: []Variable{
			{Name: "x", Value: 0, Lo: -100, Hi: 100},
			{Name: "frozen", Value: 7, Lo: 7, Hi: 7},
		},
		Residual: parabolaResidual,
	}
	res, err := Minimize(p)
	if err != nil {
		t.Fatal(err)
	}
	if res.X[1] != 7 {
		t.Fatalf("frozen: y moved to %g", res.X[1])
	}
	if math.Abs(res.X[0]-3) > 1e-8 {
		t.Fatalf("frozen: x = %g, want 3", res.X[0])
	}
}

func TestMinimizeDeterministic(t *testing.T) {
	mk := func() Problem {
		return Problem{
			Variables: []Variable{
				{Name: "x", Value: -2, Lo: -100, Hi: 100},
				{Name: "y", Value: 4, Lo: -100, Hi: 100},
			},
			Residual: parabolaResidual,
		}
	}
	r1, err := Minimize(mk())
	if err != nil {
		t.Fatal(err)
	}
	r2, err := Minimize(mk())
	if err != nil {
		t.Fatal(err)
	}
	if r1.X[0] != r2.X[0] || r1.X[1] != r2.X[1] || r1.Merit != r2.Merit {
		t.Fatalf("not deterministic: %v/%g vs %v/%g", r1.X, r1.Merit, r2.X, r2.Merit)
	}
}
