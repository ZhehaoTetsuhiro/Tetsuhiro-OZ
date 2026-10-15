// Numeric helpers of the library generator: rounding, the bounded root
// searches of the paraxial solves, the object-distance solve of the finite
// conjugate examples, and the real-ray solves of the mirror systems. The
// measurement itself lives in tools/libcheck/measure, shared with libcheck.
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// ---------- rounding ----------

func r2(x float64) float64 { return math.Round(x*1e2) / 1e2 }
func r3(x float64) float64 { return math.Round(x*1e3) / 1e3 }
func r4(x float64) float64 { return math.Round(x*1e4) / 1e4 }
func r6(x float64) float64 { return math.Round(x*1e6) / 1e6 }
func r7(x float64) float64 { return math.Round(x*1e7) / 1e7 }

// ---------- surface geometry ----------

// sagConic is the axial sag (mm) of a conic profile of vertex radius R and
// conic constant K at radius r.
func sagConic(r, R, K float64) float64 {
	if math.IsInf(R, 0) || R == 0 {
		return 0
	}
	c := 1 / R
	rad := 1 - (1+K)*c*c*r*r
	if rad < 0 {
		return math.NaN()
	}
	return c * r * r / (1 + math.Sqrt(rad))
}

// ---------- generic 1D root finding ----------

// solve1D finds a root of f on [lo,hi] by scanning for a sign change and
// bisecting. f may return NaN (skipped). Returns an error when no bracket
// exists.
func solve1D(f func(float64) float64, lo, hi float64, scan, iters int) (float64, error) {
	if scan < 8 {
		scan = 8
	}
	prevX, prevY := lo, f(lo)
	for i := 1; i <= scan; i++ {
		x := lo + (hi-lo)*float64(i)/float64(scan)
		y := f(x)
		if math.IsNaN(prevY) {
			prevX, prevY = x, y
			continue
		}
		if math.IsNaN(y) {
			continue
		}
		if y == 0 {
			return x, nil
		}
		if (prevY < 0) != (y < 0) {
			a, b := prevX, x
			fa := prevY
			for k := 0; k < iters; k++ {
				m := 0.5 * (a + b)
				fm := f(m)
				if math.IsNaN(fm) {
					break
				}
				if fm == 0 {
					return m, nil
				}
				if (fa < 0) != (fm < 0) {
					b = m
				} else {
					a, fa = m, fm
				}
			}
			return 0.5 * (a + b), nil
		}
		prevX, prevY = x, y
	}
	return 0, fmt.Errorf("no sign change on [%g,%g]", lo, hi)
}

// ---------- paraxial solves ----------

// scaleCards sets the radii of the cards idx to orig/s (curvatures × s).
func scaleCards(p *lens.Prescription, orig []float64, idx []int, s float64) {
	for k, i := range idx {
		if math.IsInf(orig[k], 0) {
			p.Cards[i].RadiusMM = orig[k]
			continue
		}
		p.Cards[i].RadiusMM = orig[k] / s
	}
}

// solveEFL scales the radii of the cards idx by a common factor (curvature ×
// s; thicknesses untouched) so that the d-line paraxial EFL equals target.
//
// The root is bracketed, not guessed. The residual 1/EFL(s) − 1/target is
// sampled on a logarithmic grid of s spanning three decades either side of
// the starting point s0 = |EFL(1)/target|; the first sign change of the
// finite samples is refined by bisection in log s (200 steps). The 1/EFL form
// is continuous through EFL = ∞, so a pole cannot be mistaken for a root.
// The refined solution must reproduce the target to a relative residual of
// 1e-6. When there is no sign change, or the solution misses the target, the
// radii are restored and an error is returned: the caller skips the
// candidate and records it in FAILURES.md.
func solveEFL(p *lens.Prescription, idx []int, target float64) error {
	if !(target != 0) || math.IsInf(target, 0) || math.IsNaN(target) {
		return fmt.Errorf("solveEFL: invalid target %v", target)
	}
	orig := make([]float64, len(idx))
	for k, i := range idx {
		orig[k] = p.Cards[i].RadiusMM
	}
	restore := func() { scaleCards(p, orig, idx, 1) }
	// g(x) is the power residual at s = 10^x; NaN when the system cannot be
	// evaluated (invalid geometry, no paraxial image, zero EFL).
	g := func(x float64) float64 {
		s := math.Pow(10, x)
		scaleCards(p, orig, idx, s)
		par, err := p.Paraxial(geo.Wavelengthd)
		if err != nil || math.IsNaN(par.EFL) || par.EFL == 0 {
			return math.NaN()
		}
		return 1/par.EFL - 1/target
	}
	seed := 1.0
	if e, err := p.Paraxial(geo.Wavelengthd); err == nil && e.EFL != 0 && !math.IsInf(e.EFL, 0) && !math.IsNaN(e.EFL) {
		seed = math.Abs(e.EFL / target)
	}
	if !(seed > 0) || math.IsInf(seed, 0) {
		seed = 1
	}
	center := math.Log10(seed)
	lo, hi := center-3, center+3
	const grid = 600
	xPrev, gPrev, have := 0.0, 0.0, false
	xa, xb, found := 0.0, 0.0, false
	for k := 0; k <= grid; k++ {
		x := lo + (hi-lo)*float64(k)/float64(grid)
		gv := g(x)
		if math.IsNaN(gv) {
			continue
		}
		if gv == 0 {
			xa, xb, found = x, x, true
			break
		}
		if have && (gPrev < 0) != (gv < 0) {
			xa, xb, found = xPrev, x, true
			break
		}
		xPrev, gPrev, have = x, gv, true
	}
	if !found {
		restore()
		return fmt.Errorf("solveEFL: no sign change of the power residual for s in [1e%.2f, 1e%.2f] (target %g mm)", lo, hi, target)
	}
	if xa != xb {
		ga := g(xa)
		for it := 0; it < 200 && xb-xa > 1e-15; it++ {
			xm := 0.5 * (xa + xb)
			gm := g(xm)
			if math.IsNaN(gm) {
				break
			}
			if gm == 0 {
				xa, xb = xm, xm
				break
			}
			if (ga < 0) != (gm < 0) {
				xb = xm
			} else {
				xa, ga = xm, gm
			}
		}
	}
	s := math.Pow(10, 0.5*(xa+xb))
	scaleCards(p, orig, idx, s)
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		restore()
		return fmt.Errorf("solveEFL: %v", err)
	}
	if rel := math.Abs(par.EFL-target) / math.Abs(target); !(rel <= 1e-6) {
		restore()
		return fmt.Errorf("solveEFL: residual %.3g (target %g, got %g)", rel, target, par.EFL)
	}
	return nil
}

// solveObjectDistance finds the finite object position z (mm) at which the
// paraxial transverse magnification of p equals mTarget (negative for a real
// inverted image).
//
// The search is clipped to the region in front of the optical surfaces: the
// object must lie at z < VertexZ(0) and z < VertexZ(1) (surface 0 is the stop
// plane, surface 1 the first lens surface). The distance d = clip − z is
// scanned on a logarithmic grid over [1e-3, 1e3]·scale mm; the first sign
// change of the residual 1/m(z) − 1/mTarget (continuous through the pole at
// m = ∞) is refined by bisection in log d. The solution must lie inside the
// clipped region and reproduce mTarget to a relative 1e-6. An out-of-range
// or unsolvable target returns an error.
func solveObjectDistance(p *lens.Prescription, mTarget, scale float64) (float64, error) {
	if !(mTarget != 0) || math.IsInf(mTarget, 0) || math.IsNaN(mTarget) || !(scale > 0) {
		return 0, fmt.Errorf("object-distance solve: invalid target %v or scale %v", mTarget, scale)
	}
	sys, err := p.System()
	if err != nil {
		return 0, err
	}
	clip := sys.VertexZ(0)
	if sys.NumSurfaces() > 1 {
		clip = math.Min(clip, sys.VertexZ(1))
	}
	// residual at the object position z; NaN when the paraxial solve fails
	resid := func(z float64) float64 {
		p.ObjectZMM = z
		par, err := p.Paraxial(geo.Wavelengthd)
		if err != nil || math.IsNaN(par.Magnification) || math.IsInf(par.Magnification, 0) || par.Magnification == 0 {
			return math.NaN()
		}
		return 1/par.Magnification - 1/mTarget
	}
	at := func(lx float64) float64 { return resid(clip - math.Pow(10, lx)) }
	lo, hi := math.Log10(1e-3*scale), math.Log10(1e3*scale)
	const grid = 800
	xPrev, gPrev, have := 0.0, 0.0, false
	xa, xb, found := 0.0, 0.0, false
	for k := 0; k <= grid; k++ {
		x := lo + (hi-lo)*float64(k)/float64(grid)
		gv := at(x)
		if math.IsNaN(gv) {
			continue
		}
		if gv == 0 {
			xa, xb, found = x, x, true
			break
		}
		if have && (gPrev < 0) != (gv < 0) {
			xa, xb, found = xPrev, x, true
			break
		}
		xPrev, gPrev, have = x, gv, true
	}
	if !found {
		p.ObjectZMM = clip - 1
		return 0, fmt.Errorf("object-distance solve: no sign change for m=%g with the object in front of z=%g (distance scan [1e-3, 1e3]×%g mm)", mTarget, clip, scale)
	}
	if xa != xb {
		ga := at(xa)
		for it := 0; it < 200 && xb-xa > 1e-15; it++ {
			xm := 0.5 * (xa + xb)
			gm := at(xm)
			if math.IsNaN(gm) {
				break
			}
			if gm == 0 {
				xa, xb = xm, xm
				break
			}
			if (ga < 0) != (gm < 0) {
				xb = xm
			} else {
				xa, ga = xm, gm
			}
		}
	}
	z := clip - math.Pow(10, 0.5*(xa+xb))
	if !(z < clip) {
		return 0, fmt.Errorf("object-distance solve: solution z=%g is not in front of surface 0 (z=%g)", z, clip)
	}
	p.ObjectZMM = z
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		return 0, fmt.Errorf("object-distance solve: %v", err)
	}
	if rel := math.Abs(par.Magnification-mTarget) / math.Abs(mTarget); !(rel <= 1e-6) {
		return 0, fmt.Errorf("object-distance solve: magnification residual %.3g (target %g, got %g)", rel, mTarget, par.Magnification)
	}
	return z, nil
}

// solveAfocal adjusts the thickness of card gapIdx until the system is
// afocal, i.e. its paraxial EFL diverges (u' = 0 for a parallel input).
func solveAfocal(p *lens.Prescription, gapIdx int, lo, hi float64) error {
	f := func(t float64) float64 {
		p.Cards[gapIdx].ThicknessMM = t
		par, err := p.Paraxial(geo.Wavelengthd)
		if err != nil {
			return math.NaN()
		}
		e := par.EFL
		if math.IsInf(e, 0) || math.IsNaN(e) {
			return 0
		}
		return 1 / e
	}
	t, err := solve1D(f, lo, hi, 400, 120)
	if err != nil {
		return fmt.Errorf("solveAfocal: %v", err)
	}
	p.Cards[gapIdx].ThicknessMM = t
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		return err
	}
	if e := math.Abs(par.EFL); !math.IsInf(e, 0) && e < 1e5 {
		return fmt.Errorf("solveAfocal: residual EFL %g", par.EFL)
	}
	return nil
}

func spanOf(sys *geo.System) float64 {
	var s float64
	for _, t := range sys.Thickness {
		s += math.Abs(t)
	}
	return s
}
