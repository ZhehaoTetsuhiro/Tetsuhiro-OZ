package optics

import (
	"fmt"
	"math"
)

// ---------------------------------------------------------------------------
// GRIN plate (grin_plate) — a plane-parallel plate of thickness d whose
// refractive index varies across the aperture.
//
// Model. The plate is a thin element (sub-millimetre transverse variation of a
// slab whose thickness cannot be resolved by the wave grid is exactly the
// regime of a phase mask): the field picks up the phase the plate adds
// relative to the reference plane,
//
//	t(x,y) = exp( i·k0·d·(n0 − 1 + Δn(x,y)) ),   k0 = 2π/λ,
//
// i.e. the transverse phase k0·Δn(x,y)·d (the physics the element is for: the
// local optical path difference of the index perturbation) plus the uniform
// piston of a slab of index n0 in air, k0·(n0−1)·d, which uses the same (n−1)
// convention as the axicon and wedge elements; set n0 = 1 to drop it. Use a
// `medium` (split-step) or `biaxial` element instead when the *bulk*
// propagation inside the glass matters (a long GRIN rod, index gradients along
// z); this element is the thin/paraxial picture, like every element in the
// library.
//
// Δn(x,y) is either the even radial polynomial
//
//	Δn(r) = c2·r² + c4·r⁴ + c6·r⁶        (r in metres)
//
// or a user expression (the scripted-element language, see docs/KERNEL.md
// §4.1) in the coordinates x, y, r, th, wl, k and the element's own numeric
// parameters. A purely quadratic profile is a thin lens:
//
//	φ = k0·c2·d·r² = −k0·r²/(2f)  →  f = −1/(2·c2·d),
//
// so c2 < 0 (index highest on the axis) converges, like a real GRIN rod; the
// c4/c6 terms are the spherical-aberration dial of the same profile. The
// expression is checked for finiteness at construction (probeElement), since a
// NaN in a thin element spreads over the whole field after the next
// propagation.
// ---------------------------------------------------------------------------

// GRINParaxialFocalLength returns the focal length of the quadratic term of a
// GRIN plate: f = −1/(2·c2·d) (m); c2 ≤ 0 converges, and c2 = 0 has no power
// (infinite focal length).
func GRINParaxialFocalLength(c2, d float64) float64 {
	if c2 == 0 || d == 0 {
		return math.Inf(1)
	}
	return -1 / (2 * c2 * d)
}

// grinPlateEl implements the "grin_plate" element.
type grinPlateEl struct {
	thickness float64 // d, m
	n0        float64 // base index of the plate
	c2, c4    float64 // Δn polynomial coefficients (1/m², 1/m⁴, ...)
	c6        float64
	expr      *CompiledExpr // optional Δn(x,y) expression
	aperture  float64       // semi-aperture radius, m (0 = unlimited)
	block     bool          // true (default): opaque outside the semi-aperture
	edgeSigma float64
	x0, y0    float64
}

func newGRINPlate(p map[string]any) (Element, error) {
	e := &grinPlateEl{
		thickness: pfd(p, "d", 2e-3),
		n0:        pfd(p, "n0", 1.5),
		c2:        pfd(p, "c2", -50),
		c4:        pfd(p, "c4", 0),
		c6:        pfd(p, "c6", 0),
		aperture:  pfd(p, "aperture", 0),
		block:     pbool(p, "block", true),
		edgeSigma: pfd(p, "edge_sigma", 0),
		x0:        pfd(p, "x", 0),
		y0:        pfd(p, "y", 0),
	}
	if !(e.thickness > 0) || math.IsInf(e.thickness, 0) {
		return nil, fmt.Errorf("grin_plate: thickness d must be > 0")
	}
	if !(e.n0 > 0) {
		return nil, fmt.Errorf("grin_plate: base index n0 must be > 0")
	}
	for _, c := range []float64{e.c2, e.c4, e.c6} {
		if math.IsNaN(c) || math.IsInf(c, 0) {
			return nil, fmt.Errorf("grin_plate: Δn coefficients must be finite")
		}
	}
	if e.aperture < 0 {
		return nil, fmt.Errorf("grin_plate: aperture must be >= 0")
	}
	if e.edgeSigma < 0 {
		return nil, fmt.Errorf("grin_plate: edge_sigma must be >= 0")
	}
	if src := ps(p, "expr", ""); src != "" {
		params := map[string]float64{
			"n0": e.n0, "d": e.thickness, "c2": e.c2, "c4": e.c4, "c6": e.c6,
		}
		ce, err := CompileExpr(src, params)
		if err != nil {
			return nil, fmt.Errorf("grin_plate: expr: %v", err)
		}
		e.expr = ce
		if err := probeElement(e); err != nil {
			return nil, fmt.Errorf("grin_plate: expr: %v", err)
		}
	}
	return e, nil
}

// deltaN returns Δn(x,y) at one point. env is nil for the polynomial path and
// is filled in for the expression path (r and th are derived only when the
// expression actually reads them).
func (e *grinPlateEl) deltaN(x, y, r2 float64, env *exprEnv) float64 {
	if e.expr != nil {
		env.x, env.y = x, y
		if e.expr.NeedsRadius() || e.expr.NeedsTheta() {
			env.r = math.Sqrt(r2)
		}
		if e.expr.NeedsTheta() {
			env.th = math.Atan2(y, x)
		}
		return e.expr.EvalEnv(env)
	}
	return r2 * (e.c2 + r2*(e.c4+r2*e.c6))
}

// Apply multiplies the field by exp(i·k0·d·(n0 − 1 + Δn(x,y))).
func (e *grinPlateEl) Apply(f *Field, ctx *Context) error {
	wl := ctx.Wavelength
	if !(wl > 0) {
		return fmt.Errorf("grin_plate: wavelength must be > 0")
	}
	k0 := 2 * math.Pi / wl
	piston := k0 * e.thickness * (e.n0 - 1)
	scale := k0 * e.thickness
	ap, sig := e.aperture, e.edgeSigma
	nx := f.N
	var env *exprEnv
	if e.expr != nil {
		env = &exprEnv{wl: wl, k: k0}
	}
	for j := 0; j < nx; j++ {
		dy := f.Y(j) - e.y0
		if env != nil {
			env.y = dy
		}
		for i := 0; i < nx; i++ {
			dx := f.X(i) - e.x0
			r2 := dx*dx + dy*dy
			w := 1.0
			if ap > 0 {
				d := ap - math.Sqrt(r2)
				if sig > 0 {
					if d < -8*sig {
						w = 0
					} else {
						w = real(smoothStep(d, sig))
					}
				} else if d < 0 {
					w = 0
				}
			}
			if w == 0 {
				if e.block {
					idx := j*nx + i
					f.Ex[idx] = 0
					if f.Polarized {
						f.Ey[idx] = 0
					}
				}
				continue
			}
			phase := piston + scale*e.deltaN(dx, dy, r2, env)
			t := cexpI(phase)
			if w != 1 {
				if e.block {
					t *= complex(w, 0)
				} else {
					t = complex(1-w, 0) + complex(w, 0)*t
				}
			}
			idx := j*nx + i
			f.Ex[idx] *= t
			if f.Polarized {
				f.Ey[idx] *= t
			}
		}
	}
	return nil
}

func init() {
	RegisterElement("grin_plate", newGRINPlate)
}
