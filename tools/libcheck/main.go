// Command libcheck validates the built-in geometric-optics example library.
//
// For every library/geo/*.json it
//
//  1. loads the file with lens.LoadJSON (so the JSON contract, the glass
//     names and the card geometry are all re-checked),
//  2. builds the geo.System,
//  3. classifies the example and measures it with the OZ engine:
//     refractive systems get a paraxial pass (EFL/BFL/F-number/F–C colour
//     shift) plus spot RMS at the axial and the off-axis field; systems the
//     paraxial engine rejects (mirrors) and systems with a virtual image are
//     measured by real ray tracing (parallel bundle → best axial focus);
//     afocal systems (telescopes, beam expanders) by the exit-angle spread
//     of a parallel bundle; collimators by the exit-angle spread of a point
//     source on the front focal plane,
//  4. asserts that every number is finite and physically sane, and
//  5. writes library/geo/index.json (file → name/category/difficulty/
//     description/tags + the measured metrics) and a summary table.
//
// It exits non-zero if any example fails, if two examples are optically
// identical, or if a file misses its metadata.
//
//	go run ./tools/libcheck [-dir library/geo] [-index library/geo/index.json]
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// ---------------------------------------------------------------- measurement

// metrics is the measured summary of one example, written verbatim into the
// index. Not every field applies to every mode; the zero value means
// "not applicable".
type metrics struct {
	Mode         string  `json:"mode"`
	EFLMM        float64 `json:"efl_mm,omitempty"`
	BFLMM        float64 `json:"bfl_mm,omitempty"`
	FNumber      float64 `json:"f_number,omitempty"`
	ColorFCMm    float64 `json:"color_fc_mm,omitempty"`
	FieldDeg     float64 `json:"field_deg,omitempty"`
	RMSOnMM      float64 `json:"rms_on_axis_mm,omitempty"`
	GEOOnMM      float64 `json:"geo_on_axis_mm,omitempty"`
	RMSOffMM     float64 `json:"rms_off_axis_mm,omitempty"`
	GEOOffMM     float64 `json:"geo_off_axis_mm,omitempty"`
	FocusZMM     float64 `json:"focus_z_mm,omitempty"`
	FocusRMSMM   float64 `json:"focus_rms_mm,omitempty"`
	SlopeRMSUrad float64 `json:"slope_rms_urad,omitempty"`
	SlopeMeanRad float64 `json:"slope_mean_rad,omitempty"`
	Mag          float64 `json:"mag,omitempty"`
	Rays         int     `json:"rays"`
}

// entry is one record of library/geo/index.json.
type entry struct {
	File        string   `json:"file"`
	Name        string   `json:"name"`
	Category    string   `json:"category"`
	Difficulty  int      `json:"difficulty"`
	Description string   `json:"description"`
	Tags        []string `json:"tags"`
	Quality     string   `json:"quality,omitempty"`
	Mode        string   `json:"validation"`
	Cards       int      `json:"cards"`
	Mirror      bool     `json:"has_mirror,omitempty"`
	Metrics     metrics  `json:"metrics"`
}

type indexFile struct {
	Tool       string           `json:"tool"`
	Dir        string           `json:"dir"`
	Count      int              `json:"count"`
	Categories map[string]int   `json:"categories"`
	Modes      map[string]int   `json:"validation_modes"`
	Quality    map[string]int   `json:"quality"`
	Sweeps     int              `json:"sweep_examples"`
	Entries    []entry          `json:"entries"`
	Summary    map[string]stats `json:"category_stats"`
}

type stats struct {
	Count      int     `json:"count"`
	RMSMedian  float64 `json:"rms_median_mm,omitempty"`
	RMSMin     float64 `json:"rms_min_mm,omitempty"`
	RMSMax     float64 `json:"rms_max_mm,omitempty"`
	BestRMS    float64 `json:"rms_best_mm,omitempty"`
	BestFile   string  `json:"rms_best_file,omitempty"`
	WorstFile  string  `json:"rms_worst_file,omitempty"`
	QualityOpt int     `json:"optimized"`
}

// fanRays traces a meridional fan of parallel rays through a system.
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

func rmsAt(pts, dirs []geo.Vec3, z float64) float64 {
	var sum, m float64
	n := 0
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
	for _, y := range ys {
		m += y
	}
	m /= float64(len(ys))
	for _, y := range ys {
		sum += (y - m) * (y - m)
		n++
	}
	return math.Sqrt(sum / float64(n))
}

// focusScan finds the axial focus (real or virtual) of a parallel bundle by
// minimising the transverse RMS of the traced rays over the image plane.
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
		return 0, 0, len(pts), fmt.Errorf("no ray-pair crossing to bracket the focus")
	}
	if math.IsNaN(zA) {
		zA = zB
	}
	if math.IsNaN(zB) {
		zB = zA
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
		r := rmsAt(pts, dirs, zz)
		if r < bestRMS {
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

// afocalScan characterises a parallel-in/parallel-out system.
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
// the entrance pupil and reports the spread of the exit slopes.
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
		pt := geo.V(frac*par.EntrancePupilRadius, 0, par.EntrancePupilZ)
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

// ---------------------------------------------------------------- validation

func fieldAngle(p *lens.Prescription) (float64, bool) {
	for _, t := range p.Tags {
		if strings.HasPrefix(t, "fov") && strings.HasSuffix(t, "deg") {
			var f float64
			if _, err := fmt.Sscanf(strings.TrimSuffix(strings.TrimPrefix(t, "fov"), "deg"), "%g", &f); err == nil {
				return f, true
			}
		}
	}
	return 0, false
}

func objFieldMM(p *lens.Prescription) (float64, bool) {
	for _, t := range p.Tags {
		if strings.HasPrefix(t, "objfield:") {
			var f float64
			if _, err := fmt.Sscanf(strings.TrimPrefix(t, "objfield:"), "%g", &f); err == nil {
				return f, true
			}
		}
	}
	return 0, false
}

func magTag(p *lens.Prescription) (float64, bool) {
	for _, t := range p.Tags {
		if strings.HasPrefix(t, "mag:") {
			var f float64
			if _, err := fmt.Sscanf(strings.TrimPrefix(t, "mag:"), "%g", &f); err == nil {
				return f, true
			}
		}
	}
	return 0, false
}

func hasTag(p *lens.Prescription, want string) bool {
	for _, t := range p.Tags {
		if t == want {
			return true
		}
	}
	return false
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

func mirrorIn(p *lens.Prescription) bool {
	for _, c := range p.Cards {
		if c.Mirror {
			return true
		}
	}
	return false
}

// validate measures one prescription and returns the metrics, the validation
// mode and the list of failed assertions.
func validate(p *lens.Prescription) (*metrics, string, []string, error) {
	var problems []string
	fail := func(format string, a ...interface{}) { problems = append(problems, fmt.Sprintf(format, a...)) }
	finite := func(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }

	sys, err := p.System()
	if err != nil {
		return nil, "", nil, fmt.Errorf("System(): %w", err)
	}
	semi := stopAperture(p)
	if !(semi > 0) {
		return nil, "", nil, fmt.Errorf("no usable stop aperture (F-number undefined)")
	}
	m := &metrics{}
	if deg, ok := fieldAngle(p); ok {
		m.FieldDeg = deg
	} else if mm, ok := objFieldMM(p); ok {
		m.FieldDeg = mm
	} else if math.IsInf(p.ObjectZMM, -1) {
		m.FieldDeg = 0.5
	} else {
		m.FieldDeg = 0.02 * math.Abs(p.ObjectZMM)
	}

	// ---- mirror systems: the paraxial engine rejects them, real trace only
	if mirrorIn(p) {
		z, rms, n, err := focusScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return nil, "focus", nil, err
		}
		m.Mode = "focus"
		m.FocusZMM, m.FocusRMSMM, m.Rays = z, rms, n
		if !finite(z) || math.Abs(z) > 1e6 {
			fail("focus z = %v is not a usable finite position", z)
		}
		if !(rms > 0) || !finite(rms) {
			fail("best-focus RMS = %v is not positive and finite", rms)
		}
		if rms > semi {
			fail("best-focus RMS %.6gmm exceeds the aperture %.6gmm (focus scan failed?)", rms, semi)
		}
		if n < 16 {
			fail("only %d of 21 traced rays survive", n)
		}
		return m, "focus", problems, nil
	}

	par, perr := sys.Paraxial(geo.Wavelengthd, p.ObjectZMM)
	if perr != nil {
		// a non-mirror system the paraxial engine rejects: real trace
		z, rms, n, err := focusScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return nil, "focus", nil, err
		}
		m.Mode, m.FocusZMM, m.FocusRMSMM, m.Rays = "focus", z, rms, n
		if !(rms > 0) || !finite(rms) {
			fail("best-focus RMS = %v is not positive and finite", rms)
		}
		return m, "focus", problems, nil
	}
	m.EFLMM, m.BFLMM, m.FNumber = par.EFL, par.BFL, par.FNumber
	if pf, err := sys.Paraxial(geo.WavelengthF, p.ObjectZMM); err == nil {
		if pc, err := sys.Paraxial(geo.WavelengthC, p.ObjectZMM); err == nil {
			m.ColorFCMm = pf.BFL - pc.BFL
		}
	}
	if !finite(m.ColorFCMm) {
		fail("F–C focal shift = %v is not finite", m.ColorFCMm)
	}

	// ---- afocal: telescopes and beam expanders
	if hasTag(p, "afocal") || math.IsInf(par.EFL, 0) || math.Abs(par.EFL) > 1e4 {
		blur, mean, mag, n, err := afocalScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return nil, "afocal", nil, err
		}
		m.Mode = "afocal"
		m.SlopeRMSUrad, m.SlopeMeanRad, m.Mag, m.Rays = blur, mean*1e-6, mag, n
		if !(blur > 0) || !finite(blur) {
			fail("exit-angle RMS = %v µrad is not positive and finite", blur)
		}
		if blur > 2e5 {
			fail("exit-angle RMS %.1f mrad exceeds 200 mrad", blur/1e3)
		}
		if math.Abs(mean) > 1e5 {
			fail("mean exit angle %.1f mrad exceeds 100 mrad", mean/1e3)
		}
		if !finite(mag) || math.Abs(mag) < 0.05 || math.Abs(mag) > 200 {
			fail("beam/magnification ratio %v is out of range", mag)
		}
		if want, ok := magTag(p); ok {
			if rel := math.Abs(mag-want) / math.Abs(want); rel > 0.2 {
				fail("measured magnification %.4g differs from the tagged %.4g by %.1f%%", mag, want, rel*100)
			}
		}
		if n < 16 {
			fail("only %d of 21 traced rays survive", n)
		}
		return m, "afocal", problems, nil
	}

	if !finite(par.EFL) || par.EFL == 0 {
		fail("EFL = %v is degenerate", par.EFL)
		return m, "refract", problems, nil
	}
	if math.Abs(par.EFL) < 0.5 || math.Abs(par.EFL) > 1e5 {
		fail("EFL %.6gmm is outside [0.5, 1e5] mm", par.EFL)
	}
	if !(par.FNumber > 0) || par.FNumber > 1e3 {
		fail("F-number %v is outside (0, 1e3]", par.FNumber)
	}

	// ---- collimator: a point source on the front focal plane
	if math.IsInf(par.ImageZ, 0) || math.Abs(par.ImageZ) > 1e6 {
		blur, mean, n, err := collimatorScan(p, sys, geo.Wavelengthd)
		if err != nil {
			return m, "collimate", problems, err
		}
		m.Mode = "collimate"
		m.SlopeRMSUrad, m.SlopeMeanRad, m.Rays = blur, mean*1e-6, n
		if !(blur > 0) || !finite(blur) {
			fail("exit-angle RMS = %v µrad is not positive and finite", blur)
		}
		if blur > 2e5 {
			fail("exit-angle RMS %.1f mrad exceeds 200 mrad", blur/1e3)
		}
		if n < 16 {
			fail("only %d of 21 traced rays survive", n)
		}
		return m, "collimate", problems, nil
	}

	// ---- virtual image: real trace of the focal region
	last := sys.VertexZ(sys.NumSurfaces() - 1)
	if par.ImageZ <= last {
		z, rms, n, err := focusScan(sys, semi, geo.Wavelengthd)
		if err != nil {
			return m, "focus", problems, err
		}
		m.Mode, m.FocusZMM, m.FocusRMSMM, m.Rays = "focus", z, rms, n
		if !(rms > 0) || !finite(rms) {
			fail("virtual-focus RMS = %v is not positive and finite", rms)
		}
		if rms > 0.5*semi {
			fail("virtual-focus RMS %.6gmm exceeds half the aperture %.6gmm", rms, semi)
		}
		return m, "focus", problems, nil
	}

	// ---- the normal case: paraxial + spot diagrams
	m.Mode = "refract"
	samples := lens.UnitDiskSamples(3, 6)
	rmsOn, geoOn, nOn, err := lens.SpotRMS(p, geo.Wavelengthd, lens.Field{}, samples)
	if err != nil {
		return m, "refract", problems, err
	}
	m.RMSOnMM, m.GEOOnMM, m.Rays = rmsOn, geoOn, nOn
	if !(rmsOn > 0) || !finite(rmsOn) {
		fail("on-axis spot RMS = %v is not positive and finite", rmsOn)
	}
	if rmsOn > 0.5*math.Abs(par.EFL) {
		fail("on-axis spot RMS %.6gmm exceeds half the focal length %.6gmm", rmsOn, par.EFL)
	}
	if nOn < len(samples)*3/5 {
		fail("only %d of %d sampled rays reach the image plane", nOn, len(samples))
	}
	off := lens.Field{}
	if math.IsInf(p.ObjectZMM, -1) {
		off = lens.Field{TanY: math.Tan(m.FieldDeg * math.Pi / 180)}
	} else {
		off = lens.Field{TanY: m.FieldDeg}
	}
	rmsOff, geoOff, nOff, err := lens.SpotRMS(p, geo.Wavelengthd, off, samples)
	if err != nil {
		fail("off-axis spot: %v", err)
	} else {
		m.RMSOffMM, m.GEOOffMM = rmsOff, geoOff
		if !(rmsOff > 0) || !finite(rmsOff) {
			fail("off-axis spot RMS = %v is not positive and finite", rmsOff)
		}
		if rmsOff > math.Abs(par.EFL) {
			fail("off-axis spot RMS %.6gmm exceeds the focal length %.6gmm", rmsOff, par.EFL)
		}
		if nOff < len(samples)*3/5 {
			fail("only %d of %d sampled off-axis rays reach the image plane", nOff, len(samples))
		}
	}
	return m, "refract", problems, nil
}

// ---------------------------------------------------------------- reporting

func round(x float64, dec int) float64 {
	p := math.Pow(10, float64(dec))
	return math.Round(x*p) / p
}

func fmtMetric(x float64) string {
	switch {
	case math.IsNaN(x):
		return "—"
	case math.Abs(x) >= 1e4:
		return fmt.Sprintf("%.3g", x)
	case math.Abs(x) >= 1:
		return fmt.Sprintf("%.3f", x)
	case math.Abs(x) >= 1e-3:
		return fmt.Sprintf("%.4f", x)
	default:
		return fmt.Sprintf("%.3g", x)
	}
}

func main() {
	dir := flag.String("dir", "library/geo", "library directory to validate")
	indexPath := flag.String("index", "library/geo/index.json", "index file to write (\"\" = do not write)")
	quiet := flag.Bool("quiet", false, "only print failures and the summary")
	flag.Parse()

	files, err := filepath.Glob(filepath.Join(*dir, "*.json"))
	if err != nil {
		fmt.Fprintln(os.Stderr, "libcheck:", err)
		os.Exit(1)
	}
	sort.Strings(files)
	if len(files) == 0 {
		fmt.Fprintf(os.Stderr, "libcheck: no JSON prescriptions in %s\n", *dir)
		os.Exit(1)
	}

	idx := indexFile{
		Tool:       "oz libcheck",
		Dir:        *dir,
		Categories: map[string]int{},
		Modes:      map[string]int{},
		Quality:    map[string]int{},
		Summary:    map[string]stats{},
	}
	sigs := map[string]string{}
	failed := 0
	bad := []string{}
	rmsByCat := map[string][]float64{}

	for _, f := range files {
		base := filepath.Base(f)
		p, err := lens.LoadJSON(f)
		if err != nil {
			fmt.Printf("FAIL %-52s load: %v\n", base, err)
			failed++
			bad = append(bad, base)
			continue
		}
		// metadata + naming invariants
		var mprob []string
		if p.Name == "" {
			mprob = append(mprob, "missing name")
		}
		if p.Category == "" {
			mprob = append(mprob, "missing category")
		}
		if p.Description == "" {
			mprob = append(mprob, "missing description")
		}
		if p.Difficulty < 1 || p.Difficulty > 5 {
			mprob = append(mprob, fmt.Sprintf("difficulty %d out of 1..5", p.Difficulty))
		}
		if len(p.Tags) < 3 {
			mprob = append(mprob, "fewer than 3 tags")
		}
		for _, r := range base {
			if r > 127 {
				mprob = append(mprob, "non-ASCII file name")
				break
			}
		}
		// duplicate detection on the optical content
		sig := signature(p)
		if prev, ok := sigs[sig]; ok {
			mprob = append(mprob, "duplicate of "+prev)
		} else {
			sigs[sig] = base
		}

		mts, mode, probs, err := validate(p)
		if err != nil {
			fmt.Printf("FAIL %-52s %s: %v\n", base, mode, err)
			failed++
			bad = append(bad, base)
			continue
		}
		probs = append(mprob, probs...)

		qual := ""
		for _, t := range p.Tags {
			if strings.HasPrefix(t, "quality:") {
				qual = strings.TrimPrefix(t, "quality:")
			}
		}
		e := entry{
			File: base, Name: p.Name, Category: p.Category, Difficulty: p.Difficulty,
			Description: p.Description, Tags: p.Tags, Quality: qual, Mode: mode,
			Cards: len(p.Cards), Mirror: mirrorIn(p),
		}
		if mts != nil {
			mtt := *mts
			mtt.EFLMM = round(mtt.EFLMM, 6)
			mtt.BFLMM = round(mtt.BFLMM, 6)
			mtt.FNumber = round(mtt.FNumber, 4)
			mtt.ColorFCMm = round(mtt.ColorFCMm, 6)
			mtt.RMSOnMM = round(mtt.RMSOnMM, 7)
			mtt.GEOOnMM = round(mtt.GEOOnMM, 7)
			mtt.RMSOffMM = round(mtt.RMSOffMM, 7)
			mtt.GEOOffMM = round(mtt.GEOOffMM, 7)
			mtt.FocusZMM = round(mtt.FocusZMM, 6)
			mtt.FocusRMSMM = round(mtt.FocusRMSMM, 7)
			mtt.SlopeRMSUrad = round(mtt.SlopeRMSUrad, 3)
			mtt.SlopeMeanRad = round(mtt.SlopeMeanRad, 9)
			mtt.Mag = round(mtt.Mag, 6)
			e.Metrics = mtt
			if v := keyRMS(mtt); !math.IsNaN(v) {
				rmsByCat[p.Category] = append(rmsByCat[p.Category], v)
			}
		}
		idx.Entries = append(idx.Entries, e)
		idx.Categories[p.Category]++
		idx.Modes[mode]++
		if qual != "" {
			idx.Quality[qual]++
		}
		if hasTag(p, "sweep:efl") || hasPrefixTag(p, "sweep:") {
			idx.Sweeps++
		}
		if len(probs) > 0 {
			failed++
			bad = append(bad, base)
			fmt.Printf("FAIL %-52s %s\n", base, strings.Join(probs, "; "))
			continue
		}
		if !*quiet {
			fmt.Printf("ok   %-52s %-9s %s\n", base, mode, metricLine(e))
		}
	}

	idx.Count = len(idx.Entries)
	// per-category metric summary
	for cat, xs := range rmsByCat {
		sort.Float64s(xs)
		st := stats{Count: len(xs), BestRMS: xs[0], RMSMin: xs[0], RMSMax: xs[len(xs)-1]}
		st.RMSMedian = xs[len(xs)/2]
		idx.Summary[cat] = st
	}
	for _, e := range idx.Entries {
		st, ok := idx.Summary[e.Category]
		if !ok {
			continue
		}
		if e.Quality == "optimized" {
			st.QualityOpt++
		}
		v := keyRMS(e.Metrics)
		if !math.IsNaN(v) {
			if st.BestFile == "" || v < st.BestRMS {
				st.BestRMS, st.BestFile = v, e.File
			}
			if v >= st.RMSMax {
				st.WorstFile = e.File
			}
		}
		idx.Summary[e.Category] = st
	}

	if *indexPath != "" {
		buf, err := json.MarshalIndent(idx, "", " ")
		if err != nil {
			fmt.Fprintln(os.Stderr, "libcheck: index:", err)
			os.Exit(1)
		}
		if err := os.WriteFile(*indexPath, append(buf, '\n'), 0o644); err != nil {
			fmt.Fprintln(os.Stderr, "libcheck: index:", err)
			os.Exit(1)
		}
	}

	// ---- summary table
	fmt.Printf("\n=== libcheck 汇总 ===\n")
	fmt.Printf("文件 %d 个，类别 %d 个，失败 %d 个\n", idx.Count, len(idx.Categories), failed)
	cats := make([]string, 0, len(idx.Categories))
	for c := range idx.Categories {
		cats = append(cats, c)
	}
	sort.Strings(cats)
	fmt.Printf("%-16s %6s %6s %-12s %-12s %-12s %-9s\n", "类别", "数量", "优化", "RMS最小值", "RMS中位数", "RMS最大值", "最差文件")
	for _, c := range cats {
		st := idx.Summary[c]
		fmt.Printf("%-16s %6d %6d %-12s %-12s %-12s %s\n", c, idx.Categories[c], st.QualityOpt,
			fmtMetric(st.RMSMin), fmtMetric(st.RMSMedian), fmtMetric(st.RMSMax), st.WorstFile)
	}
	modes := make([]string, 0, len(idx.Modes))
	for m := range idx.Modes {
		modes = append(modes, m)
	}
	sort.Strings(modes)
	for _, m := range modes {
		fmt.Printf("校验模式 %-10s %4d\n", m, idx.Modes[m])
	}
	qs := make([]string, 0, len(idx.Quality))
	for q := range idx.Quality {
		qs = append(qs, q)
	}
	sort.Strings(qs)
	for _, q := range qs {
		fmt.Printf("设计质量 %-10s %4d\n", q, idx.Quality[q])
	}
	fmt.Printf("参数扫掠示例 %d\n", idx.Sweeps)
	if *indexPath != "" {
		fmt.Printf("索引已写入 %s\n", *indexPath)
	}
	if failed > 0 {
		fmt.Printf("\n失败文件：%s\n", strings.Join(bad, ", "))
		os.Exit(1)
	}
}

func hasPrefixTag(p *lens.Prescription, pre string) bool {
	for _, t := range p.Tags {
		if strings.HasPrefix(t, pre) {
			return true
		}
	}
	return false
}

// keyRMS is the single number a category summary uses for an example.
func keyRMS(m metrics) float64 {
	switch {
	case m.RMSOnMM > 0:
		return m.RMSOnMM
	case m.FocusRMSMM > 0:
		return m.FocusRMSMM
	case m.SlopeRMSUrad > 0:
		return m.SlopeRMSUrad * 1e-6 // radians
	}
	return math.NaN()
}

func metricLine(e entry) string {
	m := e.Metrics
	switch m.Mode {
	case "refract":
		return fmt.Sprintf("EFL=%-9s F/%-6s RMS(轴向)=%-10s RMS(%g°)=%s",
			fmtMetric(m.EFLMM), fmtMetric(m.FNumber), fmtMetric(m.RMSOnMM), m.FieldDeg, fmtMetric(m.RMSOffMM))
	case "focus":
		return fmt.Sprintf("实光线追迹 焦点 z=%-10s RMS=%s (%d rays)", fmtMetric(m.FocusZMM), fmtMetric(m.FocusRMSMM), m.Rays)
	case "afocal":
		return fmt.Sprintf("无焦 出射角 RMS=%s µrad M=%s (%d rays)", fmtMetric(m.SlopeRMSUrad), fmtMetric(m.Mag), m.Rays)
	case "collimate":
		return fmt.Sprintf("准直 出射角 RMS=%s µrad (%d rays)", fmtMetric(m.SlopeRMSUrad), m.Rays)
	}
	return ""
}

// signature is the optical identity of a prescription (name excluded): two
// library entries with the same signature are the same lens.
func signature(p *lens.Prescription) string {
	obj := "inf"
	if !math.IsInf(p.ObjectZMM, 0) {
		obj = fmt.Sprintf("%.9g", p.ObjectZMM)
	}
	var b strings.Builder
	fmt.Fprintf(&b, "obj=%s stop=%d |", obj, p.Stop)
	for _, c := range p.Cards {
		fmt.Fprintf(&b, " R=%.9g K=%.9g t=%.9g g=%s a=%.9g m=%v c=%v |",
			c.RadiusMM, c.ConicK, c.ThicknessMM, c.Glass, c.SemiAperMM, c.Mirror, c.Coeffs)
	}
	return b.String()
}
