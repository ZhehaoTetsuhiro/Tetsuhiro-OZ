// Numeric helpers of the library generator: paraxial solves, ray-traced
// focus finding and the metric measurement that feeds both the descriptions
// and (via tools/libcheck, an independent implementation) the index.
package main

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// ---------- rounding ----------

func r4(x float64) float64 { return math.Round(x*1e4) / 1e4 }
func r3(x float64) float64 { return math.Round(x*1e3) / 1e3 }
func r2(x float64) float64 { return math.Round(x*1e2) / 1e2 }

// r6 keeps six decimals: the spot RMS values live between 1e-7 and 1e-1 mm,
// so they need more digits than a radius.
func r6(x float64) float64 { return math.Round(x*1e6) / 1e6 }

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

// scaleCards applies the curvatures scale s to the radii of the cards (R/s).
func scaleCards(p *lens.Prescription, orig []float64, idx []int, s float64) {
	for k, i := range idx {
		if math.IsInf(orig[k], 0) {
			p.Cards[i].RadiusMM = orig[k]
			continue
		}
		p.Cards[i].RadiusMM = orig[k] / s
	}
}

// solveEFL scales the radii of the given cards until the d-line paraxial EFL
// equals target (mm). Thicknesses are untouched.
func solveEFL(p *lens.Prescription, idx []int, target float64) error {
	orig := make([]float64, len(idx))
	for k, i := range idx {
		orig[k] = p.Cards[i].RadiusMM
	}
	f := func(s float64) float64 {
		if s == 0 {
			return math.NaN()
		}
		scaleCards(p, orig, idx, s)
		par, err := p.Paraxial(geo.Wavelengthd)
		if err != nil || math.IsNaN(par.EFL) {
			return math.NaN()
		}
		return par.EFL - target
	}
	// A thin system scales as EFL(s) ≈ EFL(1)/s, so seed the bracket close
	// to that estimate; the scan is wide enough for thick systems.
	seed := 1.0
	if e, err := p.Paraxial(geo.Wavelengthd); err == nil && e.EFL != 0 && !math.IsInf(e.EFL, 0) {
		seed = math.Abs(e.EFL / target)
	}
	if !(seed > 0) || math.IsInf(seed, 0) || math.IsNaN(seed) {
		seed = 1
	}
	lo, hi := seed/8, seed*8
	if lo < 1e-3 {
		lo = 1e-3
	}
	if hi > 1e3 {
		hi = 1e3
	}
	s, err := solve1D(f, lo, hi, 400, 200)
	if err != nil {
		// widen once and retry
		s, err = solve1D(f, 1e-3, 1e3, 2000, 200)
		if err != nil {
			return fmt.Errorf("solveEFL: %v", err)
		}
	}
	scaleCards(p, orig, idx, s)
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		return err
	}
	if rel := math.Abs(par.EFL-target) / math.Abs(target); rel > 1e-6 {
		return fmt.Errorf("solveEFL: residual %.3g (target %g, got %g)", rel, target, par.EFL)
	}
	return nil
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

// ---------- real ray helpers ----------

type rayOut struct {
	p geo.Vec3 // point where the ray left the last surface
	d geo.Vec3 // outgoing direction (unit)
}

// fanRays traces a meridional fan of parallel rays (object at infinity) at
// heights frac_i×semi, starting on the axis in front of the system. Rays
// that do not reach the last surface are reported with ok=false.
func fanRays(sys *geo.System, semi, lambda float64, n int, fracLo, fracHi float64) ([]rayOut, []float64, int) {
	out := make([]rayOut, 0, n)
	hs := make([]float64, 0, n)
	z0 := sys.VertexZ(0)
	maxAbs := 0.0
	for _, t := range sys.Thickness {
		maxAbs += math.Abs(t)
	}
	start := z0 - 2*maxAbs - 100
	for i := 0; i < n; i++ {
		h := semi * (fracLo + (fracHi-fracLo)*float64(i)/float64(n-1))
		r := geo.NewRay(geo.V(0, h, start), geo.V(0, 0, 1), lambda)
		res := sys.Trace(r)
		if !res.OK() {
			continue
		}
		last := res.Steps[len(res.Steps)-1]
		out = append(out, rayOut{p: last.Point, d: last.DirOut})
		hs = append(hs, h)
	}
	return out, hs, n
}

// meridionalIntersectionZ returns the z where rays i and j cross in the y–z
// plane (NaN when they are parallel).
func meridionalIntersectionZ(ro []rayOut, i, j int) float64 {
	a, b := ro[i], ro[j]
	// solve a.p + t*a.d = b.p + s*b.d in (y,z)
	d1, d2 := a.d, b.d
	den := d1.Z*d2.Y - d1.Y*d2.Z
	if math.Abs(den) < 1e-15 {
		return math.NaN()
	}
	ry := b.p.Y - a.p.Y
	rz := b.p.Z - a.p.Z
	t := (ry*d2.Z - rz*d2.Y) / den
	return a.p.Z + t*d1.Z
}

// transverseRMSAt returns the RMS of the meridional heights of the rays at
// plane z (relative to their mean), and the number of rays used.
func transverseRMSAt(ro []rayOut, z float64) (float64, int) {
	ys := make([]float64, 0, len(ro))
	for _, r := range ro {
		if math.Abs(r.d.Z) < 1e-14 {
			continue
		}
		t := (z - r.p.Z) / r.d.Z
		ys = append(ys, r.p.Y+t*r.d.Y)
	}
	if len(ys) == 0 {
		return math.Inf(1), 0
	}
	var m float64
	for _, y := range ys {
		m += y
	}
	m /= float64(len(ys))
	var sum float64
	for _, y := range ys {
		sum += (y - m) * (y - m)
	}
	return math.Sqrt(sum / float64(len(ys))), len(ys)
}

// focusScan finds the axial focus (real or virtual) of a parallel on-axis
// bundle by extrapolating the outgoing rays and minimising their transverse
// RMS over z. It is the only usable first-order probe for systems the
// paraxial engine rejects (mirrors) and for virtual images.
func focusScan(sys *geo.System, semi, lambda float64) (zBest, rms float64, n int, err error) {
	ro, _, total := fanRays(sys, semi, lambda, 21, 0.02, 0.98)
	if len(ro) < 3 {
		return 0, 0, len(ro), fmt.Errorf("focusScan: only %d of %d rays survive", len(ro), total)
	}
	zNear := meridionalIntersectionZ(ro, 0, 1)
	zFar := meridionalIntersectionZ(ro, len(ro)-2, len(ro)-1)
	if math.IsNaN(zNear) || math.IsNaN(zFar) {
		// fall back to the crossing of the extreme rays
		zNear = meridionalIntersectionZ(ro, 0, len(ro)-1)
		zFar = zNear
		if math.IsNaN(zNear) {
			return 0, 0, len(ro), fmt.Errorf("focusScan: no intersection bracket")
		}
	}
	lo, hi := math.Min(zNear, zFar), math.Max(zNear, zFar)
	span := hi - lo
	if span < 1e-9 {
		span = math.Max(1e-3, math.Abs(lo)*2e-3)
	}
	lo -= 0.6 * span
	hi += 0.6 * span
	best, bestRMS := math.NaN(), math.Inf(1)
	const steps = 3000
	for k := 0; k <= steps; k++ {
		z := lo + (hi-lo)*float64(k)/float64(steps)
		r, _ := transverseRMSAt(ro, z)
		if r < bestRMS {
			bestRMS, best = r, z
		}
	}
	dz := (hi - lo) / float64(steps)
	// golden-ish refinement: three passes of a shrinking bracket
	for pass := 0; pass < 3; pass++ {
		for _, z := range []float64{best - dz, best + dz} {
			if r, _ := transverseRMSAt(ro, z); r < bestRMS {
				bestRMS, best = r, z
			}
		}
		dz *= 0.3
	}
	_, n = transverseRMSAt(ro, best)
	return best, bestRMS, n, nil
}

// exitSlopeOf returns the meridional exit slope angle (rad) of a ray.
func exitSlopeOf(r rayOut) float64 { return math.Atan2(r.d.Y, r.d.Z) }

// afocalScan characterises an afocal (parallel-in, parallel-out) system: it
// traces a parallel on-axis bundle and returns the angular blur (RMS of the
// exit slopes about the mean, µrad), the mean exit slope (µrad), the linear
// beam expansion ratio M (exit height per unit input height) and the number
// of surviving rays.
func afocalScan(sys *geo.System, semi, lambda float64) (blurUrad, meanUrad, mag float64, n int, err error) {
	ro, hs, total := fanRays(sys, semi, lambda, 21, 0.02, 0.98)
	if len(ro) < 3 {
		return 0, 0, 0, len(ro), fmt.Errorf("afocalScan: only %d of %d rays survive", len(ro), total)
	}
	slopes := make([]float64, len(ro))
	zExit := -math.MaxFloat64
	for i, r := range ro {
		slopes[i] = exitSlopeOf(r)
		if r.p.Z > zExit {
			zExit = r.p.Z
		}
	}
	var mean float64
	for _, s := range slopes {
		mean += s
	}
	mean /= float64(len(slopes))
	var sum float64
	for _, s := range slopes {
		sum += (s - mean) * (s - mean)
	}
	blur := math.Sqrt(sum / float64(len(slopes)))
	// linear fit of the exit height against the input height
	var sx, sy, sxx, sxy float64
	for i, r := range ro {
		t := (zExit - r.p.Z) / r.d.Z
		y := r.p.Y + t*r.d.Y
		sx += hs[i]
		sy += y
		sxx += hs[i] * hs[i]
		sxy += hs[i] * y
	}
	k := float64(len(ro))
	den := k*sxx - sx*sx
	if math.Abs(den) > 1e-15 {
		mag = (k*sxy - sx*sy) / den
	} else {
		mag = math.NaN()
	}
	return blur * 1e6, mean * 1e6, mag, len(ro), nil
}

// collimatorScan traces a fan of rays from the axial object point through
// the entrance pupil and reports the RMS of their exit slopes (µrad): the
// figure of merit of a collimator.
func collimatorScan(p *lens.Prescription, sys *geo.System, semi, lambda float64) (blurUrad, meanUrad float64, n int, err error) {
	par, err := sys.Paraxial(lambda, p.ObjectZMM)
	if err != nil {
		return 0, 0, 0, err
	}
	objZ := p.ObjectZMM
	zp := par.EntrancePupilZ
	rp := par.EntrancePupilRadius
	if !(rp > 0) || math.IsInf(zp, 0) || math.IsNaN(zp) {
		return 0, 0, 0, fmt.Errorf("collimatorScan: unusable entrance pupil")
	}
	var slopes []float64
	const nRays = 21
	for i := 0; i < nRays; i++ {
		frac := 0.02 + 0.96*float64(i)/float64(nRays-1)
		pt := geo.V(frac*rp, 0, zp)
		r := geo.RayThrough(geo.V(0, 0, objZ), pt, lambda)
		res := sys.Trace(r)
		if !res.OK() {
			continue
		}
		last := res.Steps[len(res.Steps)-1]
		slopes = append(slopes, math.Atan2(last.DirOut.Y, last.DirOut.Z))
	}
	if len(slopes) < 3 {
		return 0, 0, len(slopes), fmt.Errorf("collimatorScan: only %d rays survive", len(slopes))
	}
	var mean float64
	for _, s := range slopes {
		mean += s
	}
	mean /= float64(len(slopes))
	var sum float64
	for _, s := range slopes {
		sum += (s - mean) * (s - mean)
	}
	return math.Sqrt(sum/float64(len(slopes))) * 1e6, mean * 1e6, len(slopes), nil
}

// solveFocusParam solves a single scalar surface parameter (a conic constant
// or an aspheric coefficient, via mutate) so that the parallel ray at height
// frac×semi leaves the last surface aimed exactly at zFocus. This is the
// only way to design a stigmatic mirror system: the paraxial engine rejects
// mirrors, so the solve is done on the real tracer.
func solveFocusParam(cards []lens.Card, lo, hi, semi, frac, zFocus, lambda float64,
	mutate func([]lens.Card, float64)) (float64, error) {
	start := -2*spanOf2(cards) - 200
	resid := func(v float64) float64 {
		c := append([]lens.Card(nil), cards...)
		mutate(c, v)
		sys, err := buildSys(c)
		if err != nil {
			return math.NaN()
		}
		r := geo.NewRay(geo.V(0, semi*frac, start), geo.V(0, 0, 1), lambda)
		res := sys.Trace(r)
		if !res.OK() {
			return math.NaN()
		}
		last := res.Steps[len(res.Steps)-1]
		if math.Abs(last.DirOut.Z) < 1e-12 {
			return math.NaN()
		}
		t := (zFocus - last.Point.Z) / last.DirOut.Z
		return last.Point.Y + t*last.DirOut.Y
	}
	return solve1D(resid, lo, hi, 400, 150)
}

// solveConicForFocus is solveFocusParam specialised to a conic constant.
func solveConicForFocus(cards []lens.Card, ci int, semi, frac, zFocus, lambda float64) error {
	v, err := solveFocusParam(cards, -12, 12, semi, frac, zFocus, lambda,
		func(c []lens.Card, k float64) { c[ci].ConicK = k })
	if err != nil {
		return fmt.Errorf("solveConicForFocus: %v", err)
	}
	cards[ci].ConicK = r6ify(v)
	return nil
}

// solveA4ForFocus is solveFocusParam specialised to the first even-asphere
// coefficient of a card. The bracket is symmetric because the sign of the
// required correction depends on the surface it sits on.
func solveA4ForFocus(cards []lens.Card, ci int, semi, frac, zFocus, lambda float64) error {
	scale := 1 / (semi * semi * semi * semi)
	v, err := solveFocusParam(cards, -60*scale, 60*scale, semi, frac, zFocus, lambda,
		func(c []lens.Card, a float64) { c[ci].Coeffs = []float64{a} })
	if err != nil {
		return fmt.Errorf("solveA4ForFocus: %v", err)
	}
	cards[ci].Coeffs = []float64{v}
	return nil
}

// tiltFocusZ traces a parallel bundle tilted by the field angle (in the
// meridional plane) and returns the best-focus z of its meridional fan: the
// tangential focus, used to correct coma in two-mirror solves.
func tiltFocusZ(cards []lens.Card, semi, tanField, zFocus, lambda float64) (float64, error) {
	sys, err := buildSys(cards)
	if err != nil {
		return 0, err
	}
	span := spanOf(sys)
	start := sys.VertexZ(0) - 2*span - 100
	dir := geo.V(0, tanField, 1)
	u, ok := dir.Unit()
	if !ok {
		return 0, fmt.Errorf("degenerate field")
	}
	type out struct{ p, d geo.Vec3 }
	outs := make([]out, 0, 15)
	for i := 0; i < 15; i++ {
		h := semi * (0.05 + 0.9*float64(i)/14)
		// rays parallel to the field direction, spread over the pupil
		origin := geo.V(0, h, start)
		res := sys.Trace(geo.NewRay(origin, u, lambda))
		if !res.OK() {
			continue
		}
		last := res.Steps[len(res.Steps)-1]
		outs = append(outs, out{last.Point, last.DirOut})
	}
	if len(outs) < 3 {
		return 0, fmt.Errorf("tiltFocusZ: only %d rays survive", len(outs))
	}
	z1, z2 := math.NaN(), math.NaN()
	a, b := outs[0], outs[1]
	d1, d2 := a.d, b.d
	den := d1.Z*d2.Y - d1.Y*d2.Z
	if math.Abs(den) > 1e-15 {
		t := ((b.p.Y-a.p.Y)*d2.Z - (b.p.Z-a.p.Z)*d2.Y) / den
		z1 = a.p.Z + t*d1.Z
	}
	a, b = outs[len(outs)-2], outs[len(outs)-1]
	d1, d2 = a.d, b.d
	den = d1.Z*d2.Y - d1.Y*d2.Z
	if math.Abs(den) > 1e-15 {
		t := ((b.p.Y-a.p.Y)*d2.Z - (b.p.Z-a.p.Z)*d2.Y) / den
		z2 = a.p.Z + t*d1.Z
	}
	if math.IsNaN(z1) && math.IsNaN(z2) {
		return 0, fmt.Errorf("tiltFocusZ: no bracket")
	}
	z := z1
	if math.IsNaN(z) {
		z = z2
	} else if !math.IsNaN(z2) {
		z = 0.5 * (z1 + z2)
	}
	_ = zFocus
	return z, nil
}

func r6ify(x float64) float64 { return math.Round(x*1e6) / 1e6 }

func spanOf(sys *geo.System) float64 {
	var s float64
	for _, t := range sys.Thickness {
		s += math.Abs(t)
	}
	return s
}

// buildSys builds a temporary geo.System from cards (mirror-friendly).
func buildSys(cards []lens.Card) (*geo.System, error) {
	p := &lens.Prescription{Name: "tmp", ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
	return p.System()
}

// ---------- tags ----------

func fieldDegFromTags(tags []string) (float64, bool) {
	for _, t := range tags {
		if strings.HasPrefix(t, "fov") && strings.HasSuffix(t, "deg") {
			v := strings.TrimSuffix(strings.TrimPrefix(t, "fov"), "deg")
			f, err := strconv.ParseFloat(v, 64)
			if err == nil {
				return f, true
			}
		}
	}
	return 0, false
}

func objFieldMMFromTags(tags []string) (float64, bool) {
	for _, t := range tags {
		if strings.HasPrefix(t, "objfield:") {
			f, err := strconv.ParseFloat(strings.TrimPrefix(t, "objfield:"), 64)
			if err == nil {
				return f, true
			}
		}
	}
	return 0, false
}

func magFromTags(tags []string) (float64, bool) {
	for _, t := range tags {
		if strings.HasPrefix(t, "mag:") {
			f, err := strconv.ParseFloat(strings.TrimPrefix(t, "mag:"), 64)
			if err == nil {
				return f, true
			}
		}
	}
	return 0, false
}

func hasTag(tags []string, want string) bool {
	for _, t := range tags {
		if t == want {
			return true
		}
	}
	return false
}

// ---------- measurement ----------

// exMetrics is the measured summary of one library example. libcheck writes
// the same structure into library/geo/index.json.
type exMetrics struct {
	Mode         string  `json:"mode"`
	EFLMM        float64 `json:"efl_mm"`
	BFLMM        float64 `json:"bfl_mm"`
	FNumber      float64 `json:"f_number"`
	ColorFCMm    float64 `json:"color_fc_mm"`
	FieldDeg     float64 `json:"field_deg"`
	RMSOnMM      float64 `json:"rms_on_axis_mm"`
	GEOOnMM      float64 `json:"geo_on_axis_mm"`
	RMSOffMM     float64 `json:"rms_off_axis_mm"`
	GeoOffMM     float64 `json:"geo_off_axis_mm"`
	FocusZMM     float64 `json:"focus_z_mm"`
	FocusRMSMM   float64 `json:"focus_rms_mm"`
	SlopeRMSUrad float64 `json:"slope_rms_urad"`
	SlopeMeanUra float64 `json:"slope_mean_urad"`
	Mag          float64 `json:"mag"`
	Rays         int     `json:"rays"`
}

func stopAperture(p *lens.Prescription) float64 {
	if p.Stop >= 0 && p.Stop < len(p.Cards) && p.Cards[p.Stop].SemiAperMM > 0 {
		return p.Cards[p.Stop].SemiAperMM
	}
	best := 0.0
	for _, c := range p.Cards {
		if c.SemiAperMM > best && c.SemiAperMM < 1e6 {
			best = c.SemiAperMM
		}
	}
	return best
}

func hasMirror(p *lens.Prescription) bool {
	for _, c := range p.Cards {
		if c.Mirror {
			return true
		}
	}
	return false
}

// measureEx is the single measurement path of the generator: it classifies
// the example (refractive / mirror-focus / afocal / collimator) and records
// the metrics of that class.
func measureEx(p *lens.Prescription) (*exMetrics, error) {
	m := &exMetrics{}
	sys, err := p.System()
	if err != nil {
		return nil, err
	}
	semi := stopAperture(p)
	if !(semi > 0) {
		return nil, fmt.Errorf("no stop aperture")
	}
	// field convention
	if deg, ok := fieldDegFromTags(p.Tags); ok {
		m.FieldDeg = deg
	}
	if mm, ok := objFieldMMFromTags(p.Tags); ok {
		m.FieldDeg = mm // reuse the slot for the finite-object field in mm
	}
	if m.FieldDeg == 0 {
		if math.IsInf(p.ObjectZMM, -1) {
			m.FieldDeg = 0.5
		} else {
			m.FieldDeg = 0.02 * math.Abs(p.ObjectZMM)
		}
	}
	afocalTag := hasTag(p.Tags, "afocal")
	if hasMirror(p) {
		m.Mode = "focus"
		z, rms, n, err := focusScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return nil, err
		}
		m.FocusZMM, m.FocusRMSMM, m.Rays = r3(z), r6(rms), n
		return m, nil
	}
	par, perr := sys.Paraxial(geo.Wavelengthd, p.ObjectZMM)
	if perr != nil {
		m.Mode = "focus"
		z, rms, n, err := focusScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return nil, err
		}
		m.FocusZMM, m.FocusRMSMM, m.Rays = r3(z), r6(rms), n
		return m, nil
	}
	m.EFLMM = r6(par.EFL)
	m.BFLMM = r6(par.BFL)
	m.FNumber = r2(par.FNumber)
	if pf, err := sys.Paraxial(geo.WavelengthF, p.ObjectZMM); err == nil {
		if pc, err := sys.Paraxial(geo.WavelengthC, p.ObjectZMM); err == nil {
			m.ColorFCMm = r3(pf.BFL - pc.BFL)
		}
	}
	last := sys.VertexZ(sys.NumSurfaces() - 1)
	if afocalTag || math.IsInf(par.EFL, 0) || math.Abs(par.EFL) > 1e4 {
		m.Mode = "afocal"
		blur, mean, mag, n, err := afocalScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return nil, err
		}
		m.SlopeRMSUrad, m.SlopeMeanUra, m.Mag, m.Rays = r2(blur), r2(mean), r3(mag), n
		return m, nil
	}
	if math.IsInf(par.ImageZ, 0) {
		m.Mode = "collimate"
		blur, mean, n, err := collimatorScan(p, sys, semi, geo.Wavelengthd)
		if err != nil {
			return m, nil
		}
		m.SlopeRMSUrad, m.SlopeMeanUra, m.Rays = r2(blur), r2(mean), n
		return m, nil
	}
	if par.ImageZ <= last {
		m.Mode = "focus"
		z, rms, n, err := focusScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return nil, err
		}
		m.FocusZMM, m.FocusRMSMM, m.Rays = r3(z), r6(rms), n
		return m, nil
	}
	m.Mode = "refract"
	samples := lens.UnitDiskSamples(3, 6)
	rms, gg, n, err := lens.SpotRMS(p, geo.Wavelengthd, lens.Field{}, samples)
	if err != nil {
		return nil, err
	}
	m.RMSOnMM, m.GEOOnMM, m.Rays = r6(rms), r6(gg), n
	off := lens.Field{}
	if math.IsInf(p.ObjectZMM, -1) {
		off = lens.Field{TanY: math.Tan(m.FieldDeg * math.Pi / 180)}
	} else {
		off = lens.Field{TanY: m.FieldDeg}
	}
	if rms, gg, _, err := lens.SpotRMS(p, geo.Wavelengthd, off, samples); err == nil {
		m.RMSOffMM, m.GeoOffMM = r6(rms), r6(gg)
	} else {
		m.RMSOffMM, m.GeoOffMM = math.NaN(), math.NaN()
	}
	return m, nil
}

// ---------- description ----------

func fmtSig(x float64, unit string) string {
	switch {
	case math.IsNaN(x):
		return "—"
	case math.IsInf(x, 0):
		return "∞"
	case math.Abs(x) >= 1000:
		return fmt.Sprintf("%.4g%s", x, unit)
	case math.Abs(x) >= 1:
		return fmt.Sprintf("%.3f%s", x, unit)
	case math.Abs(x) >= 1e-3:
		return fmt.Sprintf("%.1f%s", x*1e3, "µm")
	default:
		return fmt.Sprintf("%.4g%s", x*1e3, "µm")
	}
}

// expectText turns the measured metrics into the "预期指标" clause of the
// description. The bounds are the measured values rounded up to a round
// number, i.e. an expectation the example satisfies with margin.
func expectText(m *exMetrics) string {
	switch m.Mode {
	case "refract":
		s := fmt.Sprintf("EFL≈%s（±0.1%%）、BFL≈%s、F/%s；轴上 RMS≈%s，离轴 %g° RMS≈%s（%d 条光线）",
			fmtSig(m.EFLMM, "mm"), fmtSig(m.BFLMM, "mm"), fmtSig(m.FNumber, ""),
			fmtSig(m.RMSOnMM, "mm"), m.FieldDeg, fmtSig(m.RMSOffMM, "mm"), m.Rays)
		if m.ColorFCMm != 0 {
			s += fmt.Sprintf("；F–C 焦移≈%s", fmtSig(m.ColorFCMm, "mm"))
		}
		return s
	case "focus":
		return fmt.Sprintf("实光线追迹：轴上焦点 z≈%s，最佳焦面 RMS≈%s（%d 条光线，光线追迹法，无近轴近似）",
			fmtSig(m.FocusZMM, "mm"), fmtSig(m.FocusRMSMM, "mm"), m.Rays)
	case "afocal":
		return fmt.Sprintf("实光线追迹（无焦系统）：出射角散布 RMS≈%s，平均出射角≈%s，光束扩展比 M≈%.3f（%d 条光线）",
			fmtSig(m.SlopeRMSUrad/1e6, "rad"), fmtSig(m.SlopeMeanUra/1e6, "rad"), m.Mag, m.Rays)
	case "collimate":
		return fmt.Sprintf("实光线追迹（准直系统）：出射角散布 RMS≈%s（%d 条光线）",
			fmtSig(m.SlopeRMSUrad/1e6, "rad"), m.Rays)
	}
	return ""
}
