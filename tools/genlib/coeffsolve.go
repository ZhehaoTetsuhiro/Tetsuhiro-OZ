// Coefficient solving for the aspheric example families: the bounded
// coordinate search that fits polynomial and Forbes coefficients to the
// engine's own spot RMS, and the small engine wrappers the families of
// catalog_asphere.go and catalog_posed.go are built from.
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/tools/libcheck/measure"
)

// ---------- engine wrappers ----------

// spotOnAxis returns the on-axis spot RMS (mm) of p at the d line, measured
// exactly the way the library measurement does (lens.SpotRMS at the paraxial
// image plane over the unit-disk sampling of the stop).
func spotOnAxis(p *lens.Prescription) (float64, error) {
	rms, _, _, err := lens.SpotRMS(p, geo.Wavelengthd, lens.Field{}, lens.UnitDiskSamples(3, 6))
	return rms, err
}

// marginalHeightAt traces the axial marginal ray (the ray through 99.9% of
// the stop aperture of an object at infinity) and returns the radius of its
// hit point at the given card. It is the beam radius that scales the
// aspheric coefficients of that surface.
func marginalHeightAt(p *lens.Prescription, card int) (float64, error) {
	if !math.IsInf(p.ObjectZMM, -1) {
		return 0, fmt.Errorf("marginal height: only an object at infinity is supported")
	}
	semi := measure.StopAperture(p)
	if !(semi > 0) {
		return 0, fmt.Errorf("marginal height: no usable stop aperture")
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		return 0, err
	}
	sys, err := p.System()
	if err != nil {
		return 0, err
	}
	res := sys.Trace(geo.NewRay(geo.V(0, 0.999*semi, par.EntrancePupilZ), geo.V(0, 0, 1), geo.Wavelengthd))
	if !res.OK() {
		return 0, fmt.Errorf("marginal height: %v", res.Error())
	}
	for _, st := range res.Steps {
		if st.Surface == card {
			return math.Hypot(st.Point.X, st.Point.Y), nil
		}
	}
	return 0, fmt.Errorf("marginal height: ray does not reach card %d", card)
}

// exitAngleDeg returns the angle (degrees) between the mean direction of the
// parallel on-axis bundle after the last surface and the +z axis: the fold of
// a decentred or tilted system, measured on the real tracer.
func exitAngleDeg(p *lens.Prescription) (float64, error) {
	semi := measure.StopAperture(p)
	sys, err := p.System()
	if err != nil {
		return 0, err
	}
	pts, dirs, _ := fanOf(sys, semi)
	if len(pts) < 3 {
		return 0, fmt.Errorf("exit angle: only %d rays reach the last surface", len(pts))
	}
	var mean geo.Vec3
	for _, d := range dirs {
		mean = mean.Add(d)
	}
	u, ok := mean.Unit()
	if !ok {
		return 0, fmt.Errorf("exit angle: no mean exit direction")
	}
	return math.Atan2(math.Hypot(u.X, u.Y), u.Z) * 180 / math.Pi, nil
}

// fanOf traces the same 21-ray fan of a parallel on-axis bundle that the
// library measurement uses (heights 0.02 … 0.98 of the stop aperture).
func fanOf(sys *geo.System, semi float64) (pts, dirs []geo.Vec3, hs []float64) {
	var span float64
	for _, t := range sys.Thickness {
		span += math.Abs(t)
	}
	start := sys.VertexZ(0) - 2*span - 200
	for i := 0; i < 21; i++ {
		h := semi * (0.02 + 0.96*float64(i)/20)
		res := sys.Trace(geo.NewRay(geo.V(0, h, start), geo.V(0, 0, 1), geo.Wavelengthd))
		if !res.OK() {
			continue
		}
		last := res.Steps[len(res.Steps)-1]
		pts = append(pts, last.Point)
		dirs = append(dirs, last.DirOut)
		hs = append(hs, h)
	}
	return pts, dirs, hs
}

// foldAtCardDeg returns the angle (degrees) between the incoming and the
// outgoing direction of the traced fan at one surface: the fold a mirror
// imposes on the beam, read off the trace itself. Only the steps of rays that
// reach the surface count.
func foldAtCardDeg(p *lens.Prescription, card int) (float64, error) {
	semi := measure.StopAperture(p)
	sys, err := p.System()
	if err != nil {
		return 0, err
	}
	var span float64
	for _, t := range sys.Thickness {
		span += math.Abs(t)
	}
	start := sys.VertexZ(0) - 2*span - 200
	var in, out geo.Vec3
	n := 0
	for i := 0; i < 21; i++ {
		h := semi * (0.02 + 0.96*float64(i)/20)
		res := sys.Trace(geo.NewRay(geo.V(0, h, start), geo.V(0, 0, 1), geo.Wavelengthd))
		for _, st := range res.Steps {
			if st.Surface != card {
				continue
			}
			if st.Status == geo.StatusOK {
				in, out = in.Add(st.DirIn), out.Add(st.DirOut)
				n++
			}
			break
		}
	}
	if n < 3 {
		return 0, fmt.Errorf("fold angle: only %d rays interact with surface %d", n, card)
	}
	u, ok1 := in.Unit()
	v, ok2 := out.Unit()
	if !ok1 || !ok2 {
		return 0, fmt.Errorf("fold angle: degenerate directions at surface %d", card)
	}
	d := math.Max(-1, math.Min(1, u.Dot(v)))
	return math.Acos(d) * 180 / math.Pi, nil
}

// ---------- coefficient search ----------

// coeffVar is one optimisable coefficient: a named pair of accessors into a
// prescription (a polynomial term, a Forbes coefficient, …).
type coeffVar struct {
	name string
	get  func(*lens.Prescription) float64
	set  func(*lens.Prescription, float64)
}

// bracket is the search interval of one coefficient.
type bracket struct{ lo, hi float64 }

// min1D minimises g on [lo,hi]: a coarse scan of `scan` points locates the
// basin, a golden-section search refines it, and the best point of the scan
// is kept whenever the refinement is worse (so a non-unimodal g can only lose
// to the scan, never to the refinement). Deterministic: fixed grid and
// iteration count, no randomness.
func min1D(g func(float64) float64, lo, hi float64, scan, iters int) (float64, float64) {
	if !(hi > lo) {
		return lo, g(lo)
	}
	best, bestV := lo, math.Inf(1)
	for i := 0; i <= scan; i++ {
		x := lo + (hi-lo)*float64(i)/float64(scan)
		v := g(x)
		if !math.IsNaN(v) && v < bestV {
			best, bestV = x, v
		}
	}
	if math.IsInf(bestV, 1) {
		return best, bestV
	}
	d := (hi - lo) / float64(scan)
	a, b := math.Max(lo, best-d), math.Min(hi, best+d)
	const phi = 0.6180339887498949
	c, e := b-phi*(b-a), a+phi*(b-a)
	fc, fe := g(c), g(e)
	for i := 0; i < iters; i++ {
		if fc < fe {
			b, e, fe = e, c, fc
			c = b - phi*(b-a)
			fc = g(c)
		} else {
			a, c, fc = c, e, fe
			e = a + phi*(b-a)
			fe = g(e)
		}
	}
	x := 0.5 * (a + b)
	v := g(x)
	if math.IsNaN(v) || v > bestV {
		return best, bestV
	}
	return x, v
}

// solveCoeffs fits the coefficient slots vs of p to the on-axis spot RMS by
// a bounded coordinate search: `passes` sweeps over the slots in order, each
// slot scanned on a 121-point grid of its bracket and refined by golden
// section (60 iterations). It leaves the best coefficients in p and returns
// the measured RMS. The search evaluates nothing but the public engine
// (lens.SpotRMS), so every number the families record is an engine number.
func solveCoeffs(p *lens.Prescription, vs []coeffVar, br []bracket, passes int) (float64, error) {
	if len(vs) == 0 || len(vs) != len(br) {
		return 0, fmt.Errorf("solveCoeffs: %d slots, %d brackets", len(vs), len(br))
	}
	best := math.Inf(1)
	for pass := 0; pass < passes; pass++ {
		for i, v := range vs {
			g := func(x float64) float64 {
				v.set(p, x)
				r, err := spotOnAxis(p)
				if err != nil {
					return math.NaN()
				}
				return r
			}
			x, val := min1D(g, br[i].lo, br[i].hi, 121, 60)
			if math.IsNaN(val) {
				return 0, fmt.Errorf("solveCoeffs: %s has no usable value on [%g,%g]", v.name, br[i].lo, br[i].hi)
			}
			v.set(p, x)
			if val < best {
				best = val
			}
		}
	}
	rms, err := spotOnAxis(p)
	if err != nil {
		return 0, fmt.Errorf("solveCoeffs: final measurement: %w", err)
	}
	return rms, nil
}

// sig6 rounds x to six significant digits (0, NaN and ±Inf pass through):
// the frozen coefficient form written into the library.
func sig6(x float64) float64 {
	if x == 0 || math.IsNaN(x) || math.IsInf(x, 0) {
		return x
	}
	e := math.Pow(10, 6-math.Ceil(math.Log10(math.Abs(x))))
	return math.Round(x*e) / e
}

// polyVar returns the coefficient slot of one polynomial term of a card. The
// card must already carry the term (see polyTermsOf).
func polyVar(p *lens.Prescription, card int, power int) coeffVar {
	idx := -1
	for i, t := range p.Cards[card].PolyTerms {
		if t.Power == power {
			idx = i
		}
	}
	return coeffVar{
		name: fmt.Sprintf("card %d a%d", card, power),
		get: func(q *lens.Prescription) float64 {
			if idx < 0 {
				return 0
			}
			return q.Cards[card].PolyTerms[idx].Coeff
		},
		set: func(q *lens.Prescription, v float64) {
			if idx < 0 {
				return
			}
			q.Cards[card].PolyTerms[idx].Coeff = v
		},
	}
}

// polyTermsOf installs the polynomial terms of the powers on the card, all
// zero; the slots returned by polyVar then address them by power.
func polyTermsOf(p *lens.Prescription, card int, powers []int) {
	terms := make([]lens.PolyTermCard, len(powers))
	for i, pw := range powers {
		terms[i] = lens.PolyTermCard{Power: pw, Coeff: 0}
	}
	p.Cards[card].PolyTerms = terms
}

// qconVar returns the coefficient slot of one Forbes Qcon coefficient.
func qconVar(p *lens.Prescription, card, m int) coeffVar {
	return coeffVar{
		name: fmt.Sprintf("card %d s%d", card, m),
		get:  func(q *lens.Prescription) float64 { return q.Cards[card].QAsphere.Coeffs[m] },
		set:  func(q *lens.Prescription, v float64) { q.Cards[card].QAsphere.Coeffs[m] = v },
	}
}

// qbfsVar returns the coefficient slot of one Forbes Qbfs coefficient. The
// base radius is re-derived on every write so that the vertex curvature (and
// with it the first-order design) stays at c0: the Qbfs departure contributes
// 2·S(0)/ρ_max² to the vertex curvature (geo.QAsphere.Curvature), so
// R = 1/(c0 − 2·S(0)/ρ²) is exact.
func qbfsVar(p *lens.Prescription, card, m int, c0 float64) coeffVar {
	return coeffVar{
		name: fmt.Sprintf("card %d a%d", card, m),
		get:  func(q *lens.Prescription) float64 { return q.Cards[card].QAsphere.Coeffs[m] },
		set: func(q *lens.Prescription, v float64) {
			qc := q.Cards[card].QAsphere
			qc.Coeffs[m] = v
			q.Cards[card].RadiusMM = qbfsRadius(qc, c0)
		},
	}
}

// qbfsRadius returns the base radius that keeps the vertex curvature of the
// Qbfs card at c0 (see qbfsVar). Only the public geo.QAsphere.Curvature is
// used, so the value is the engine's own curvature convention.
func qbfsRadius(qc *lens.QCard, c0 float64) float64 {
	ref := geo.QAsphere{Form: geo.QBfs, R: 1, NormR: qc.NormRadiusMM, Coeffs: qc.Coeffs}
	s0 := (ref.Curvature() - 1) * qc.NormRadiusMM * qc.NormRadiusMM / 2
	den := c0 - 2*s0/(qc.NormRadiusMM*qc.NormRadiusMM)
	if den == 0 || math.IsNaN(den) || math.IsInf(den, 0) {
		return math.NaN()
	}
	return 1 / den
}
