// Package measure is the single definition of what the example library
// records about a prescription: how it is classified (refractive, real or
// virtual focus, afocal, collimator), how it is measured with the OZ engine
// (package geo for ray tracing and first-order optics, package lens for the
// spot diagrams), which physical-sanity assertions it must pass, the
// parameter fingerprint that identifies it, and the machine-readable
// expectation clause written into every description.
//
// tools/genlib measures each candidate and writes the clause into its
// description; tools/libcheck reloads the saved JSON, measures it again with
// this same code and compares the numbers of the clause with the fresh
// measurement.
package measure

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"unicode"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// Measurement modes.
const (
	ModeRefract   = "refract"   // paraxial EFL/BFL/F-number and spot diagrams
	ModeFocus     = "focus"     // real-ray best focus of a parallel bundle (mirrors, virtual images)
	ModeAfocal    = "afocal"    // exit-angle spread of a parallel bundle (telescopes, expanders)
	ModeCollimate = "collimate" // exit-angle spread of a point source on the front focal plane
	ModePosed     = "posed"     // real-ray best focus of a decentred or tilted system
)

// The description clause is compared with a fresh measurement with the
// tolerance |got-want| <= CompareTol*|want| + CompareAbs. The clause prints
// nine significant digits, so rounding alone stays far below CompareTol.
const (
	CompareTol = 1e-6
	CompareAbs = 1e-12
)

// Metrics is the measured summary of one prescription. Only the fields of
// its Mode are meaningful; the others are zero.
type Metrics struct {
	Mode          string  `json:"mode"`
	EFLMM         float64 `json:"efl_mm,omitempty"`
	BFLMM         float64 `json:"bfl_mm,omitempty"`
	FNumber       float64 `json:"f_number,omitempty"`
	ColorFCMM     float64 `json:"color_fc_mm,omitempty"`
	Field         float64 `json:"field,omitempty"`
	FieldUnit     string  `json:"field_unit,omitempty"`
	RMSAxisMM     float64 `json:"rms_on_axis_mm,omitempty"`
	GEOAxisMM     float64 `json:"geo_on_axis_mm,omitempty"`
	RMSFieldMM    float64 `json:"rms_field_mm,omitempty"`
	GEOFieldMM    float64 `json:"geo_field_mm,omitempty"`
	FocusZMM      float64 `json:"focus_z_mm,omitempty"`
	FocusRMSMM    float64 `json:"focus_rms_mm,omitempty"`
	SlopeRMSUrad  float64 `json:"slope_rms_urad,omitempty"`
	SlopeMeanUrad float64 `json:"slope_mean_urad,omitempty"`
	Mag           float64 `json:"mag,omitempty"`
	Rays          int     `json:"rays"`
	FieldRays     int     `json:"field_rays,omitempty"`
	StopSemiMM    float64 `json:"stop_semi_aperture_mm"`
}

// Measure classifies p and measures it. Every number it returns is finite:
// quantities that do not apply to the mode are zero.
func Measure(p *lens.Prescription) (*Metrics, error) {
	sys, err := p.System()
	if err != nil {
		return nil, err
	}
	semi := StopAperture(p)
	if !(semi > 0) {
		return nil, fmt.Errorf("no usable stop aperture (F-number undefined)")
	}
	m := &Metrics{StopSemiMM: semi, Field: fieldOf(p), FieldUnit: fieldUnit(p)}
	if posed(p) {
		// A decentred or tilted surface makes the first-order engine refuse the
		// system (geo: first- and third-order optics are coaxial only), so the
		// example is measured on the real tracer alone.
		return posedMetrics(m, sys, semi)
	}
	if HasMirror(p) {
		// the paraxial engine rejects reflecting surfaces: real-ray focus
		return focusMetrics(m, sys, semi)
	}
	par, perr := sys.Paraxial(geo.Wavelengthd, p.ObjectZMM)
	if perr != nil {
		return focusMetrics(m, sys, semi)
	}
	m.EFLMM, m.BFLMM, m.FNumber = par.EFL, par.BFL, par.FNumber
	if pf, err := sys.Paraxial(geo.WavelengthF, p.ObjectZMM); err == nil {
		if pc, err := sys.Paraxial(geo.WavelengthC, p.ObjectZMM); err == nil {
			m.ColorFCMM = pf.BFL - pc.BFL
		}
	}
	if hasTag(p.Tags, "afocal") || math.IsInf(par.EFL, 0) || math.Abs(par.EFL) > 1e4 {
		return afocalMetrics(m, sys, semi)
	}
	if math.IsInf(par.ImageZ, 0) || math.Abs(par.ImageZ) > 1e6 {
		return collimateMetrics(m, p, sys)
	}
	if par.ImageZ <= sys.VertexZ(sys.NumSurfaces()-1) {
		return focusMetrics(m, sys, semi)
	}
	return refractMetrics(m, p)
}

func focusMetrics(m *Metrics, sys *geo.System, semi float64) (*Metrics, error) {
	z, rms, n, err := focusScan(sys, semi, geo.Wavelengthd)
	if err != nil {
		return nil, err
	}
	m.Mode = ModeFocus
	m.FocusZMM, m.FocusRMSMM, m.Rays = z, rms, n
	m.clearParaxial()
	return m, nil
}

// posedMetrics measures a system whose first-order engine is unavailable
// because a surface is decentred or tilted (geo's paraxial engine is defined
// for coaxial systems only). The measurement is real-ray only: the best-focus
// plane of a parallel on-axis bundle in the plane perpendicular to the mean
// exit direction (see posedScan).
func posedMetrics(m *Metrics, sys *geo.System, semi float64) (*Metrics, error) {
	z, rms, n, err := posedScan(sys, semi, geo.Wavelengthd)
	if err != nil {
		return nil, err
	}
	m.Mode = ModePosed
	m.FocusZMM, m.FocusRMSMM, m.Rays = z, rms, n
	m.clearParaxial()
	return m, nil
}

func afocalMetrics(m *Metrics, sys *geo.System, semi float64) (*Metrics, error) {
	blur, mean, mag, n, err := afocalScan(sys, semi, geo.Wavelengthd)
	if err != nil {
		return nil, err
	}
	m.Mode = ModeAfocal
	m.SlopeRMSUrad, m.SlopeMeanUrad, m.Mag, m.Rays = blur, mean, mag, n
	m.clearParaxial()
	return m, nil
}

func collimateMetrics(m *Metrics, p *lens.Prescription, sys *geo.System) (*Metrics, error) {
	blur, mean, n, err := collimatorScan(p, sys, geo.Wavelengthd)
	if err != nil {
		return nil, err
	}
	m.Mode = ModeCollimate
	m.SlopeRMSUrad, m.SlopeMeanUrad, m.Rays = blur, mean, n
	m.clearParaxial()
	return m, nil
}

func refractMetrics(m *Metrics, p *lens.Prescription) (*Metrics, error) {
	samples := lens.UnitDiskSamples(3, 6)
	rmsOn, geoOn, nOn, err := lens.SpotRMS(p, geo.Wavelengthd, lens.Field{}, samples)
	if err != nil {
		return nil, err
	}
	rmsOff, geoOff, nOff, err := lens.SpotRMS(p, geo.Wavelengthd, offField(p, m.Field), samples)
	if err != nil {
		return nil, err
	}
	m.Mode = ModeRefract
	m.RMSAxisMM, m.GEOAxisMM, m.Rays = rmsOn, geoOn, nOn
	m.RMSFieldMM, m.GEOFieldMM, m.FieldRays = rmsOff, geoOff, nOff
	return m, nil
}

// clearParaxial zeroes the first-order quantities, which only refractive
// systems carry (afocal systems have infinite EFL, collimators and mirrors
// have no finite paraxial image).
func (m *Metrics) clearParaxial() {
	m.EFLMM, m.BFLMM, m.FNumber, m.ColorFCMM = 0, 0, 0, 0
	m.RMSAxisMM, m.GEOAxisMM, m.RMSFieldMM, m.GEOFieldMM, m.FieldRays = 0, 0, 0, 0, 0
}

func offField(p *lens.Prescription, deg float64) lens.Field {
	if math.IsInf(p.ObjectZMM, -1) {
		return lens.Field{TanY: math.Tan(deg * math.Pi / 180)}
	}
	return lens.Field{TanY: deg}
}

// fieldUnit is "deg" for an object at infinity (field angle) and "mm" for a
// finite object (object height).
func fieldUnit(p *lens.Prescription) string {
	if math.IsInf(p.ObjectZMM, -1) {
		return "deg"
	}
	return "mm"
}

// fieldOf returns the field of the example: the tagged field angle (deg)
// for an object at infinity, the tagged object height (mm) for a finite
// object, else a default.
func fieldOf(p *lens.Prescription) float64 {
	for _, t := range p.Tags {
		if strings.HasPrefix(t, "fov") && strings.HasSuffix(t, "deg") {
			if v, err := strconv.ParseFloat(strings.TrimSuffix(strings.TrimPrefix(t, "fov"), "deg"), 64); err == nil {
				return v
			}
		}
	}
	for _, t := range p.Tags {
		if strings.HasPrefix(t, "objfield:") {
			if v, err := strconv.ParseFloat(strings.TrimPrefix(t, "objfield:"), 64); err == nil {
				return v
			}
		}
	}
	if math.IsInf(p.ObjectZMM, -1) {
		return 0.5
	}
	return 0.02 * math.Abs(p.ObjectZMM)
}

// StopAperture is the clear semi-aperture of the stop card (or the largest
// finite aperture when the stop card has none).
func StopAperture(p *lens.Prescription) float64 {
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

// HasMirror reports whether any surface of p reflects.
func HasMirror(p *lens.Prescription) bool {
	for _, c := range p.Cards {
		if c.Mirror {
			return true
		}
	}
	return false
}

// posed reports whether any card of p is decentred or tilted (the
// coordinate-break fields of lens.Card). Such a system has no first-order
// optics: geo's paraxial engine refuses it, so it is measured on the real
// tracer alone (ModePosed).
func posed(p *lens.Prescription) bool {
	for _, c := range p.Cards {
		if c.DecenterXMM != 0 || c.DecenterYMM != 0 || c.TiltXDeg != 0 || c.TiltYDeg != 0 || c.TiltZDeg != 0 {
			return true
		}
	}
	return false
}

func hasTag(tags []string, want string) bool {
	for _, t := range tags {
		if t == want {
			return true
		}
	}
	return false
}

func magTag(tags []string) (float64, bool) {
	for _, t := range tags {
		if strings.HasPrefix(t, "mag:") {
			if v, err := strconv.ParseFloat(strings.TrimPrefix(t, "mag:"), 64); err == nil {
				return v, true
			}
		}
	}
	return 0, false
}

// ---------------------------------------------------------------- scans

// fanRays traces parallel rays (object at infinity, travelling along +z)
// at heights semi·frac, frac ∈ [fracLo, fracHi], starting in front of the
// system. Rays that do not reach the last surface are dropped.
func fanRays(sys *geo.System, semi, lambda float64, n int, fracLo, fracHi float64) (pts, dirs []geo.Vec3, hs []float64) {
	var span float64
	for _, t := range sys.Thickness {
		span += math.Abs(t)
	}
	start := sys.VertexZ(0) - 2*span - 200
	for i := 0; i < n; i++ {
		h := semi * (fracLo + (fracHi-fracLo)*float64(i)/float64(n-1))
		res := sys.Trace(geo.NewRay(geo.V(0, h, start), geo.V(0, 0, 1), lambda))
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

func crossZ(p1, d1, p2, d2 geo.Vec3) float64 {
	den := d1.Z*d2.Y - d1.Y*d2.Z
	if math.Abs(den) < 1e-15 {
		return math.NaN()
	}
	t := ((p2.Y-p1.Y)*d2.Z - (p2.Z-p1.Z)*d2.Y) / den
	return p1.Z + t*d1.Z
}

// rmsAt is the RMS (about the mean) of the meridional heights of the rays
// on the plane z.
func rmsAt(pts, dirs []geo.Vec3, z float64) float64 {
	ys := make([]float64, 0, len(pts))
	for i := range pts {
		if math.Abs(dirs[i].Z) < 1e-14 {
			continue
		}
		ys = append(ys, pts[i].Y+(z-pts[i].Z)*dirs[i].Y/dirs[i].Z)
	}
	if len(ys) == 0 {
		return math.Inf(1)
	}
	var m, sum float64
	for _, y := range ys {
		m += y
	}
	m /= float64(len(ys))
	for _, y := range ys {
		sum += (y - m) * (y - m)
	}
	return math.Sqrt(sum / float64(len(ys)))
}

// focusScan locates the axial focus (real or virtual) of a parallel
// on-axis bundle: the plane that minimises the transverse RMS of the traced
// rays, searched on a grid bracketed by the crossings of the extreme ray
// pairs and refined three times.
func focusScan(sys *geo.System, semi, lambda float64) (z, rms float64, n int, err error) {
	pts, dirs, _ := fanRays(sys, semi, lambda, 21, 0.02, 0.98)
	if len(pts) < 3 {
		return 0, 0, len(pts), fmt.Errorf("only %d of 21 rays reach the last surface", len(pts))
	}
	zA := crossZ(pts[0], dirs[0], pts[1], dirs[1])
	zB := crossZ(pts[len(pts)-2], dirs[len(pts)-2], pts[len(pts)-1], dirs[len(pts)-1])
	if math.IsNaN(zA) {
		zA = crossZ(pts[0], dirs[0], pts[len(pts)-1], dirs[len(pts)-1])
	}
	if math.IsNaN(zB) {
		zB = zA
	}
	if math.IsNaN(zA) && math.IsNaN(zB) {
		return 0, 0, len(pts), fmt.Errorf("no ray-pair crossing brackets the focus")
	}
	if math.IsNaN(zA) {
		zA = zB
	}
	lo, hi := math.Min(zA, zB), math.Max(zA, zB)
	span := hi - lo
	if span < 1e-9 {
		span = math.Max(1e-3, math.Abs(lo)*2e-3)
	}
	lo -= 0.6 * span
	hi += 0.6 * span
	best, bestRMS := math.NaN(), math.Inf(1)
	const steps = 3000
	for k := 0; k <= steps; k++ {
		zz := lo + (hi-lo)*float64(k)/float64(steps)
		if r := rmsAt(pts, dirs, zz); r < bestRMS {
			bestRMS, best = r, zz
		}
	}
	dz := (hi - lo) / float64(steps)
	for pass := 0; pass < 3; pass++ {
		for _, zz := range []float64{best - dz, best + dz} {
			if r := rmsAt(pts, dirs, zz); r < bestRMS {
				bestRMS, best = r, zz
			}
		}
		dz *= 0.3
	}
	return best, bestRMS, len(pts), nil
}

// posedScan locates the real-ray best focus of a parallel on-axis bundle
// through a system the first-order engine refuses (decentred or tilted
// surfaces). The traced exit rays are projected onto planes perpendicular to
// their mean exit direction, and the plane of least transverse RMS about its
// centroid is searched. The search window is bracketed by the crossings of
// adjacent traced ray pairs projected on the mean direction, so an oblique
// exit beam does not move the focus out of the window. The returned z is the
// z coordinate of the centroid of the rays on that plane (mm).
func posedScan(sys *geo.System, semi, lambda float64) (focusZ, rms float64, n int, err error) {
	pts, dirs, _ := fanRays(sys, semi, lambda, 21, 0.02, 0.98)
	if len(pts) < 3 {
		return 0, 0, len(pts), fmt.Errorf("only %d of 21 rays reach the last surface", len(pts))
	}
	var sum geo.Vec3
	for _, d := range dirs {
		sum = sum.Add(d)
	}
	axis, ok := sum.Unit()
	if !ok {
		return 0, 0, len(pts), fmt.Errorf("the exit bundle has no mean direction")
	}
	ref := pts[0]
	// at returns the point of ray i on the plane {axis·(x − ref) = u}.
	at := func(i int, u float64) (geo.Vec3, bool) {
		den := axis.Dot(dirs[i])
		if math.Abs(den) < 1e-9 {
			return geo.Zero, false
		}
		t := (u - axis.Dot(pts[i].Sub(ref))) / den
		return pts[i].Add(dirs[i].Mul(t)), true
	}
	// rmsPlane is the transverse RMS of the ray points on that plane.
	rmsPlane := func(u float64) (float64, bool) {
		var cx, cy, cz float64
		xs := make([]geo.Vec3, 0, len(pts))
		for i := range pts {
			p, ok := at(i, u)
			if !ok {
				continue
			}
			xs = append(xs, p)
			cx, cy, cz = cx+p.X, cy+p.Y, cz+p.Z
		}
		if len(xs) < 3 {
			return 0, false
		}
		k := float64(len(xs))
		cx, cy, cz = cx/k, cy/k, cz/k
		var s float64
		for _, p := range xs {
			s += (p.X-cx)*(p.X-cx) + (p.Y-cy)*(p.Y-cy) + (p.Z-cz)*(p.Z-cz)
		}
		return math.Sqrt(s / k), true
	}
	var us []float64
	for i := 0; i+1 < len(pts); i++ {
		if u, ok := pairCrossU(pts[i], dirs[i], pts[i+1], dirs[i+1], axis, ref); ok {
			us = append(us, u)
		}
	}
	if len(us) < 2 {
		return 0, 0, len(pts), fmt.Errorf("no ray-pair crossing brackets the focus")
	}
	sort.Float64s(us)
	lo, hi := us[0], us[len(us)-1]
	span := hi - lo
	if span < 1e-9 {
		span = math.Max(1e-3, math.Abs(lo)*2e-3)
	}
	lo -= 0.5*span + 1
	hi += 0.5*span + 1
	best, bestRMS := math.NaN(), math.Inf(1)
	const steps = 2000
	for k := 0; k <= steps; k++ {
		u := lo + (hi-lo)*float64(k)/float64(steps)
		if r, ok := rmsPlane(u); ok && r < bestRMS {
			bestRMS, best = r, u
		}
	}
	if math.IsNaN(best) {
		return 0, 0, len(pts), fmt.Errorf("no usable focus plane on [%g,%g]", lo, hi)
	}
	du := (hi - lo) / float64(steps)
	for pass := 0; pass < 3; pass++ {
		for _, u := range []float64{best - du, best + du} {
			if r, ok := rmsPlane(u); ok && r < bestRMS {
				bestRMS, best = r, u
			}
		}
		du *= 0.3
	}
	var cz float64
	cnt := 0
	for i := range pts {
		if p, ok := at(i, best); ok {
			cz += p.Z
			cnt++
		}
	}
	if cnt == 0 {
		return 0, 0, len(pts), fmt.Errorf("no ray lands on the best-focus plane")
	}
	return cz / float64(cnt), bestRMS, len(pts), nil
}

// pairCrossU returns the mean-direction coordinate of the midpoint of the
// shortest segment between two traced rays, or ok = false for an (almost)
// parallel pair, whose crossing is not defined.
func pairCrossU(p1, d1, p2, d2, axis, ref geo.Vec3) (float64, bool) {
	f := d1.Dot(p2.Sub(p1))
	g := d2.Dot(p2.Sub(p1))
	b := d1.Dot(d2)
	den := 1 - b*b
	if math.Abs(den) < 1e-12 {
		return 0, false
	}
	s := (f - b*g) / den
	u := (b*f - g) / den
	m1 := p1.Add(d1.Mul(s))
	m2 := p2.Add(d2.Mul(u))
	mid := geo.V(m1.X+0.5*(m2.X-m1.X), m1.Y+0.5*(m2.Y-m1.Y), m1.Z+0.5*(m2.Z-m1.Z))
	return axis.Dot(mid.Sub(ref)), true
}

// afocalScan characterises a parallel-in, parallel-out system: the RMS
// (about the mean) of the exit slopes in µrad, the mean exit slope in µrad,
// and the linear beam expansion (exit height per input height) fitted over
// the bundle.
func afocalScan(sys *geo.System, semi, lambda float64) (blurUrad, meanUrad, mag float64, n int, err error) {
	pts, dirs, hs := fanRays(sys, semi, lambda, 21, -0.98, 0.98)
	if len(pts) < 5 {
		return 0, 0, 0, len(pts), fmt.Errorf("only %d of 21 rays reach the last surface", len(pts))
	}
	slopes := make([]float64, len(pts))
	zExit := -math.MaxFloat64
	for i := range pts {
		slopes[i] = math.Atan2(dirs[i].Y, dirs[i].Z)
		if pts[i].Z > zExit {
			zExit = pts[i].Z
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
	var sx, sy, sxx, sxy float64
	for i := range pts {
		t := (zExit - pts[i].Z) / dirs[i].Z
		y := pts[i].Y + t*dirs[i].Y
		sx += hs[i]
		sy += y
		sxx += hs[i] * hs[i]
		sxy += hs[i] * y
	}
	k := float64(len(pts))
	if den := k*sxx - sx*sx; math.Abs(den) > 1e-15 {
		mag = (k*sxy - sx*sy) / den
	} else {
		mag = math.NaN()
	}
	return blur * 1e6, mean * 1e6, mag, len(pts), nil
}

// collimatorScan traces a fan of rays from the axial object point through
// the entrance pupil and reports the RMS and the mean of the exit slopes
// (µrad): the figure of merit of a collimator.
func collimatorScan(p *lens.Prescription, sys *geo.System, lambda float64) (blurUrad, meanUrad float64, n int, err error) {
	par, err := sys.Paraxial(lambda, p.ObjectZMM)
	if err != nil {
		return 0, 0, 0, err
	}
	if !(par.EntrancePupilRadius > 0) || math.IsInf(par.EntrancePupilZ, 0) || math.IsNaN(par.EntrancePupilZ) {
		return 0, 0, 0, fmt.Errorf("unusable entrance pupil")
	}
	var slopes []float64
	const nRays = 21
	for i := 0; i < nRays; i++ {
		frac := -0.98 + 1.96*float64(i)/float64(nRays-1)
		// the fan lies in the y–z plane, like every other scan, so the y slope is the one measured
		pt := geo.V(0, frac*par.EntrancePupilRadius, par.EntrancePupilZ)
		res := sys.Trace(geo.RayThrough(geo.V(0, 0, p.ObjectZMM), pt, lambda))
		if !res.OK() {
			continue
		}
		last := res.Steps[len(res.Steps)-1]
		slopes = append(slopes, math.Atan2(last.DirOut.Y, last.DirOut.Z))
	}
	if len(slopes) < 5 {
		return 0, 0, len(slopes), fmt.Errorf("only %d of 21 rays reach the last surface", len(slopes))
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

// ---------------------------------------------------------------- checks

// Check returns the physical-sanity problems of a measurement; an empty
// result means the example passes. These are the only assertions of the
// library: the generator and libcheck both apply them.
func Check(p *lens.Prescription, m *Metrics) []string {
	var probs []string
	fail := func(format string, a ...interface{}) { probs = append(probs, fmt.Sprintf(format, a...)) }
	finite := func(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }
	if m.Rays < 16 && m.Mode != ModeRefract {
		fail("only %d of 21 traced rays survive", m.Rays)
	}
	switch m.Mode {
	case ModeFocus:
		if !finite(m.FocusZMM) || math.Abs(m.FocusZMM) > 1e6 {
			fail("focus z = %v is not a usable finite position", m.FocusZMM)
		}
		if !(m.FocusRMSMM > 0) || !finite(m.FocusRMSMM) {
			fail("focus RMS = %v is not positive and finite", m.FocusRMSMM)
		}
		if HasMirror(p) {
			if m.FocusRMSMM > m.StopSemiMM {
				fail("focus RMS %.6gmm exceeds the aperture %.6gmm", m.FocusRMSMM, m.StopSemiMM)
			}
		} else if m.FocusRMSMM > 0.5*m.StopSemiMM {
			fail("virtual-focus RMS %.6gmm exceeds half the aperture %.6gmm", m.FocusRMSMM, m.StopSemiMM)
		}
	case ModeAfocal:
		if !(m.SlopeRMSUrad > 0) || !finite(m.SlopeRMSUrad) {
			fail("exit-angle RMS = %v µrad is not positive and finite", m.SlopeRMSUrad)
		}
		if m.SlopeRMSUrad > 2e5 {
			fail("exit-angle RMS %.1f mrad exceeds 200 mrad", m.SlopeRMSUrad/1e3)
		}
		if math.Abs(m.SlopeMeanUrad) > 1e5 {
			fail("mean exit angle %.1f mrad exceeds 100 mrad", m.SlopeMeanUrad/1e3)
		}
		if !finite(m.Mag) || math.Abs(m.Mag) < 1e-4 || math.Abs(m.Mag) > 1e3 {
			fail("beam ratio %v is out of range", m.Mag)
		}
		// An expander's tag is its expansion ratio |M|. An inverting (Keplerian)
		// beam has a negative measured ratio, so the magnitudes are compared.
		if want, ok := magTag(p.Tags); ok {
			if rel := math.Abs(math.Abs(m.Mag)-want) / math.Abs(want); !(rel <= 0.2) {
				fail("measured magnification %.4g differs from the tagged %.4g by %.1f%%", m.Mag, want, rel*100)
			}
		}
	case ModePosed:
		if !finite(m.FocusZMM) || math.Abs(m.FocusZMM) > 1e6 {
			fail("focus z = %v is not a usable finite position", m.FocusZMM)
		}
		if !(m.FocusRMSMM > 0) || !finite(m.FocusRMSMM) {
			fail("focus RMS = %v is not positive and finite", m.FocusRMSMM)
		}
		if m.FocusRMSMM > m.StopSemiMM {
			fail("focus RMS %.6gmm exceeds the aperture %.6gmm", m.FocusRMSMM, m.StopSemiMM)
		}
	case ModeCollimate:
		if !(m.SlopeRMSUrad > 0) || !finite(m.SlopeRMSUrad) {
			fail("exit-angle RMS = %v µrad is not positive and finite", m.SlopeRMSUrad)
		}
		if m.SlopeRMSUrad > 2e5 {
			fail("exit-angle RMS %.1f mrad exceeds 200 mrad", m.SlopeRMSUrad/1e3)
		}
	case ModeRefract:
		samples := len(lens.UnitDiskSamples(3, 6))
		if !finite(m.EFLMM) || m.EFLMM == 0 {
			fail("EFL = %v is degenerate", m.EFLMM)
			break
		}
		if math.Abs(m.EFLMM) < 0.5 || math.Abs(m.EFLMM) > 1e5 {
			fail("EFL %.6gmm is outside [0.5, 1e5] mm", m.EFLMM)
		}
		if !(m.FNumber > 0) || m.FNumber > 1e3 {
			fail("F-number %v is outside (0, 1e3]", m.FNumber)
		}
		if !finite(m.ColorFCMM) {
			fail("F-C focal shift = %v is not finite", m.ColorFCMM)
		}
		if !(m.RMSAxisMM > 0) || !finite(m.RMSAxisMM) {
			fail("on-axis spot RMS = %v is not positive and finite", m.RMSAxisMM)
		}
		if m.RMSAxisMM > 0.5*math.Abs(m.EFLMM) {
			fail("on-axis spot RMS %.6gmm exceeds half the focal length", m.RMSAxisMM)
		}
		if m.Rays < samples*3/5 {
			fail("only %d of %d sampled rays reach the image plane", m.Rays, samples)
		}
		if !(m.RMSFieldMM > 0) || !finite(m.RMSFieldMM) {
			fail("field spot RMS = %v is not positive and finite", m.RMSFieldMM)
		}
		if m.RMSFieldMM > math.Abs(m.EFLMM) {
			fail("field spot RMS %.6gmm exceeds the focal length", m.RMSFieldMM)
		}
		if m.FieldRays < samples*3/5 {
			fail("only %d of %d sampled off-axis rays reach the image plane", m.FieldRays, samples)
		}
	default:
		fail("unknown mode %q", m.Mode)
	}
	return probs
}

// ---------------------------------------------------------------- clause

// num formats a number of the clause: nine significant digits.
func num(v float64) string { return strconv.FormatFloat(v, 'g', 9, 64) }

type field struct {
	key  string
	val  float64
	unit string
}

// fields lists the numbers the clause of a mode carries, in clause order.
func fields(m *Metrics) []field {
	var fs []field
	switch m.Mode {
	case ModeRefract:
		fs = []field{
			{"EFL", m.EFLMM, "mm"}, {"BFL", m.BFLMM, "mm"}, {"FNO", m.FNumber, ""},
			{"FC", m.ColorFCMM, "mm"}, {"RMSAX", m.RMSAxisMM, "mm"}, {"FIELD", m.Field, m.FieldUnit},
			{"RMSFLD", m.RMSFieldMM, "mm"}, {"RAYS", float64(m.Rays), ""}, {"FRAYS", float64(m.FieldRays), ""},
		}
	case ModeFocus:
		fs = []field{{"FOCUSZ", m.FocusZMM, "mm"}, {"FOCUSRMS", m.FocusRMSMM, "mm"}, {"RAYS", float64(m.Rays), ""}}
	case ModePosed:
		fs = []field{{"FOCUSZ", m.FocusZMM, "mm"}, {"FOCUSRMS", m.FocusRMSMM, "mm"}, {"RAYS", float64(m.Rays), ""}}
	case ModeAfocal:
		fs = []field{{"SLOPERMS", m.SlopeRMSUrad, "µrad"}, {"SLOPEMEAN", m.SlopeMeanUrad, "µrad"},
			{"MAG", m.Mag, ""}, {"RAYS", float64(m.Rays), ""}}
	case ModeCollimate:
		fs = []field{{"SLOPERMS", m.SlopeRMSUrad, "µrad"}, {"SLOPEMEAN", m.SlopeMeanUrad, "µrad"},
			{"RAYS", float64(m.Rays), ""}}
	}
	return fs
}

// ExpectClause renders the expectation clause written at the end of every
// description. It is machine-readable (KEY=value pairs separated by "；"),
// so that libcheck can compare it with a fresh measurement.
func ExpectClause(m *Metrics) string {
	parts := []string{"MODE=" + m.Mode}
	for _, f := range fields(m) {
		parts = append(parts, f.key+"="+num(f.val)+f.unit)
	}
	return "预期指标（引擎计算；单位 mm、µrad、deg；libcheck 复核容差：相对 1e-6）：" + strings.Join(parts, "；")
}

var pairRe = regexp.MustCompile(`([A-Z][A-Z0-9]*)=([^；。\s]*)`)

// ParseClause reads the expectation clause of a description. The number of
// a pair is its leading numeric part; a trailing unit is ignored.
func ParseClause(desc string) (mode string, vals map[string]float64, err error) {
	i := strings.Index(desc, "MODE=")
	if i < 0 {
		return "", nil, fmt.Errorf("description has no expectation clause (MODE=)")
	}
	clause := desc[i:]
	if j := strings.Index(clause, "。"); j >= 0 {
		clause = clause[:j]
	}
	vals = map[string]float64{}
	for _, mt := range pairRe.FindAllStringSubmatch(clause, -1) {
		key, raw := mt[1], mt[2]
		if key == "MODE" {
			mode = raw
			continue
		}
		if _, dup := vals[key]; dup {
			return "", nil, fmt.Errorf("clause repeats key %s", key)
		}
		numPart := strings.TrimRightFunc(raw, func(r rune) bool { return !unicode.IsDigit(r) && r != '.' })
		v, perr := strconv.ParseFloat(numPart, 64)
		if perr != nil {
			return "", nil, fmt.Errorf("clause value %s=%q is not a number", key, raw)
		}
		vals[key] = v
	}
	if mode == "" {
		return "", nil, fmt.Errorf("clause has no MODE")
	}
	return mode, vals, nil
}

// Compare checks a parsed clause against a fresh measurement and returns
// the mismatches (empty = the description matches the recomputed numbers).
func Compare(m *Metrics, mode string, vals map[string]float64) []string {
	var probs []string
	if mode != m.Mode {
		probs = append(probs, fmt.Sprintf("description mode %q, measured mode %q", mode, m.Mode))
	}
	seen := map[string]bool{}
	for _, f := range fields(m) {
		seen[f.key] = true
		got, ok := vals[f.key]
		if !ok {
			probs = append(probs, "clause misses "+f.key)
			continue
		}
		if math.Abs(got-f.val) > CompareTol*math.Abs(f.val)+CompareAbs {
			probs = append(probs, fmt.Sprintf("%s: description %s, recomputed %s", f.key, num(got), num(f.val)))
		}
	}
	var extra []string
	for k := range vals {
		if !seen[k] {
			extra = append(extra, k)
		}
	}
	sort.Strings(extra)
	for _, k := range extra {
		probs = append(probs, "clause has unexpected key "+k)
	}
	return probs
}

// ---------------------------------------------------------------- identity

// Canonical rewrites the prescription into its JSON-safe form: flat
// surfaces (RadiusMM ±Inf) are written as 0, which the lens package reads
// identically (its documented convention). It does not change the optics.
func Canonical(p *lens.Prescription) {
	for i := range p.Cards {
		if math.IsInf(p.Cards[i].RadiusMM, 0) {
			p.Cards[i].RadiusMM = 0
		}
	}
}

type fpCard struct {
	RadiusMM    float64   `json:"radius_mm"`
	ConicK      float64   `json:"conic_k"`
	Coeffs      []float64 `json:"coeffs"`
	ThicknessMM float64   `json:"thickness_mm"`
	Glass       string    `json:"glass"`
	SemiAperMM  float64   `json:"semi_aperture_mm"`
	Mirror      bool      `json:"mirror"`
	// The aspheric descriptions and the coordinate break of the card are part
	// of the design: two examples that differ only in their polynomial terms,
	// their Forbes coefficients or their decentre/tilt are different optics.
	// Every field is omitempty so that prescriptions written before they
	// existed keep their fingerprints byte for byte.
	PolyTerms   []lens.PolyTermCard `json:"poly_terms,omitempty"`
	QAsphere    *lens.QCard         `json:"q_asphere,omitempty"`
	DecenterXMM float64             `json:"decenter_x_mm,omitempty"`
	DecenterYMM float64             `json:"decenter_y_mm,omitempty"`
	TiltXDeg    float64             `json:"tilt_x_deg,omitempty"`
	TiltYDeg    float64             `json:"tilt_y_deg,omitempty"`
	TiltZDeg    float64             `json:"tilt_z_deg,omitempty"`
}

type fpDoc struct {
	ObjectZMM interface{} `json:"object_z_mm"`
	Stop      int         `json:"stop"`
	Cards     []fpCard    `json:"cards"`
}

// Fingerprint is the parameter fingerprint of a prescription: a SHA-256 of
// its canonical parameter list (object position, stop, and per card the
// radius, conic, coefficients, thickness, glass, aperture, mirror flag, the
// polynomial or Forbes aspheric description and the decentre/tilt).
// Names, categories, tags and descriptions are excluded, so two entries
// with the same fingerprint are the same optical design.
func Fingerprint(p *lens.Prescription) (string, error) {
	d := fpDoc{Stop: p.Stop}
	if math.IsInf(p.ObjectZMM, -1) {
		d.ObjectZMM = "inf"
	} else {
		d.ObjectZMM = p.ObjectZMM
	}
	for _, c := range p.Cards {
		r := c.RadiusMM
		if math.IsInf(r, 0) {
			r = 0
		}
		d.Cards = append(d.Cards, fpCard{RadiusMM: r, ConicK: c.ConicK, Coeffs: c.Coeffs,
			ThicknessMM: c.ThicknessMM, Glass: c.Glass, SemiAperMM: c.SemiAperMM, Mirror: c.Mirror,
			PolyTerms: c.PolyTerms, QAsphere: c.QAsphere,
			DecenterXMM: c.DecenterXMM, DecenterYMM: c.DecenterYMM,
			TiltXDeg: c.TiltXDeg, TiltYDeg: c.TiltYDeg, TiltZDeg: c.TiltZDeg})
	}
	buf, err := json.Marshal(d)
	if err != nil {
		return "", err
	}
	sum := sha256.Sum256(buf)
	return hex.EncodeToString(sum[:]), nil
}
