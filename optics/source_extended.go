package optics

import (
	"fmt"
	"math"
)

// ExtendedSourceParams describes a spatially incoherent extended emitter: a
// rectangular patch of radiating surface whose points are mutually incoherent.
//
// The kernel models it the only way that is exact for incoherent light: the
// patch is sampled into point sources, each simulated as its own coherent
// unit, and the recorded intensity is the sum over the points. Two points of
// the patch are not coherent with each other, so no single complex field can
// represent the emitter.
//
// Width/Height are the full extents (m) along the two axes spanning the
// emitter plane (u horizontal, v vertical). A zero extent collapses that axis
// to a single on-axis sample, which degenerates to a plain point source — the
// model's continuity check.
//
// SamplesX/SamplesY are the sample counts per axis (>= 1). The total power is
// split between the samples: evenly for profile "uniform", weighted by
// exp(-2 r^2/w^2) for "gaussian".
//
// Standoff is the distance (m) from the emitter plane at which each point's
// spherical wavefront is sampled. The kernel builds a source field on the
// source's own plane, so a point emitter is represented by a spherical source
// placed R in front of it with curvature radius R: the virtual point then sits
// exactly on the emitter whatever R is, and the propagated result is
// independent of R. R is therefore a pure sampling choice — a larger R samples
// a flatter wavefront (kinder to the grid), and it must stay below the distance
// to the first element in front of the emitter, which is what the automatic
// default (0.95 of that distance) enforces.
//
// Note on absolute power: the kernel normalises a source's power over the grid
// window at the source plane, so the power an emitter delivers into a small
// downstream aperture converges as the window grows. The relative weights
// between the emitter's points, and the incoherent sum, are exact.
type ExtendedSourceParams struct {
	Width    float64 // full extent along u (m); 0 = one on-axis sample
	Height   float64 // full extent along v (m); 0 = one on-axis sample
	SamplesX int     // samples along u (>= 1)
	SamplesY int     // samples along v (>= 1)
	Power    float64 // total emitted power (W), split over the samples
	Profile  string  // "uniform" (default) or "gaussian"
	Standoff float64 // sampling distance in front of the emitter (m); 0 = auto
}

// maxExtendedSamples bounds the point sources one extended emitter may expand
// to. Each point is a full simulation of the whole light path, so the cost is
// linear in this number; the bound keeps one mis-typed JSON field from turning
// into hours of computation.
const maxExtendedSamples = 4096

// extendedDefaultStandoff is the sampling distance used when nothing sits in
// front of the emitter (the light escapes anyway, so only the sampling plane
// has to be somewhere sane).
const extendedDefaultStandoff = 0.05

// extendedSampleCentres returns the coordinates of n samples evenly covering
// the extent [−ext/2, +ext/2]: the cell centres. n <= 1 (or a zero extent)
// collapses to a single on-axis sample.
func extendedSampleCentres(ext float64, n int) []float64 {
	if n <= 1 || ext <= 0 {
		return []float64{0}
	}
	cell := ext / float64(n)
	out := make([]float64, n)
	for i := range out {
		out[i] = -ext/2 + (float64(i)+0.5)*cell
	}
	return out
}

// extendedWeights returns the per-sample power fractions for a profile,
// normalised so they sum to one. The gaussian profile is
// exp(-2 r^2 / w^2) with w the 1/e^2 *half*-extent of the axis (width/2), the
// same convention as the gaussian source's waist. A collapsed axis contributes
// nothing to the exponent (its sample sits on the axis of that direction).
func extendedWeights(us, vs []float64, profile string, width, height float64) [][]float64 {
	w := make([][]float64, len(vs))
	var sum float64
	// sig2 is the squared 1/e^2 half-extent per axis; a collapsed axis (zero
	// width/height) keeps sig2 = 0 and contributes nothing to the exponent,
	// because its single sample sits on the axis of that direction.
	sig2u, sig2v := 0.0, 0.0
	if width > 0 {
		sig2u = width * width / 4
	}
	if height > 0 {
		sig2v = height * height / 4
	}
	for j, v := range vs {
		row := make([]float64, len(us))
		for i, u := range us {
			switch profile {
			case "gaussian":
				var r2 float64
				if sig2u > 0 {
					r2 += u * u / sig2u
				}
				if sig2v > 0 {
					r2 += v * v / sig2v
				}
				row[i] = math.Exp(-2 * r2)
			default:
				row[i] = 1
			}
			sum += row[i]
		}
		w[j] = row
	}
	if sum > 0 {
		for j := range w {
			for i := range w[j] {
				w[j][i] /= sum
			}
		}
	}
	return w
}

// extendedSourceSpec is the parsed, validated form of an extended source.
type extendedSourceSpec struct {
	p  ExtendedSourceParams
	us []float64
	vs []float64
	w  [][]float64
}

// parseExtendedSource validates the extended-source parameters. A nil map is a
// valid (empty) parameter set, which degenerates to a single on-axis point.
func parseExtendedSource(params map[string]any) (*extendedSourceSpec, error) {
	s := &extendedSourceSpec{p: ExtendedSourceParams{
		Width:    pfd(params, "width", 0),
		Height:   pfd(params, "height", 0),
		SamplesX: pi_(params, "samples_x", 1),
		SamplesY: pi_(params, "samples_y", 1),
		Power:    pfd(params, "power", 1e-3),
		Profile:  ps(params, "profile", "uniform"),
		Standoff: pfd(params, "standoff", 0),
	}}
	if s.p.Width < 0 || s.p.Height < 0 {
		return nil, fmt.Errorf("extended: width and height must be >= 0 (0 collapses that axis to one on-axis sample)")
	}
	for _, axis := range []struct {
		key string
		n   int
	}{{"samples_x", s.p.SamplesX}, {"samples_y", s.p.SamplesY}} {
		v, ok := params[axis.key]
		if !ok || v == nil {
			continue
		}
		f, err := asFloat(v)
		if err != nil {
			return nil, fmt.Errorf("extended: %s must be a whole number, got %v", axis.key, v)
		}
		if f != math.Trunc(f) {
			return nil, fmt.Errorf("extended: %s must be a whole number, got %v", axis.key, v)
		}
	}
	if s.p.SamplesX < 1 || s.p.SamplesY < 1 {
		return nil, fmt.Errorf("extended: samples_x and samples_y must be >= 1 (got %d, %d)",
			s.p.SamplesX, s.p.SamplesY)
	}
	if s.p.SamplesX*s.p.SamplesY > maxExtendedSamples {
		return nil, fmt.Errorf("extended: samples_x=%d, samples_y=%d => %d point sources exceeds the limit of %d "+
			"(each point is a full simulation of the light path)",
			s.p.SamplesX, s.p.SamplesY, s.p.SamplesX*s.p.SamplesY, maxExtendedSamples)
	}
	switch s.p.Profile {
	case "", "uniform", "gaussian":
		if s.p.Profile == "" {
			s.p.Profile = "uniform"
		}
	default:
		return nil, fmt.Errorf("extended: unknown profile %q (have: uniform, gaussian)", s.p.Profile)
	}
	if s.p.Standoff < 0 {
		return nil, fmt.Errorf("extended: standoff must be >= 0 (0 = automatic)")
	}
	for _, f := range []float64{s.p.Width, s.p.Height, s.p.Standoff} {
		if math.IsNaN(f) || math.IsInf(f, 0) {
			return nil, fmt.Errorf("extended: width/height/standoff must be finite")
		}
	}
	if !(s.p.Power > 0) {
		return nil, fmt.Errorf("extended: power must be > 0")
	}
	s.us = extendedSampleCentres(s.p.Width, s.p.SamplesX)
	s.vs = extendedSampleCentres(s.p.Height, s.p.SamplesY)
	s.w = extendedWeights(s.us, s.vs, s.p.Profile, s.p.Width, s.p.Height)
	return s, nil
}

// ExpandExtendedSources rewrites every extended source of a configuration
// into the point sources it stands for. The result carries one entry per
// emitted point, in row-major order (v outer, u inner), each as its own
// incoherent unit (an empty group key, which the unit splitter turns into a
// per-source unit). Warnings describe anything the expansion had to adjust;
// the caller adds them to the run so nothing happens silently.
//
// The expansion is a pure function of its inputs: identical input gives
// identical output, which keeps the simulation deterministic.
func ExpandExtendedSources(sources []SourceSpec, scene *SceneSpec) ([]SourceSpec, []string) {
	var out []SourceSpec
	var warns []string
	seen := 0
	for i := range sources {
		s := &sources[i]
		if s.Type != "extended" {
			out = append(out, *s)
			continue
		}
		spec, err := parseExtendedSource(s.Params)
		if err != nil {
			// An invalid emitter is reported, not silently dropped: its light
			// would otherwise vanish without a trace.
			warns = append(warns, fmt.Sprintf("扩展光源 %s 参数无效（%v），已忽略该光源", s.ResolvedLabel(i), err))
			continue
		}
		if len(out)+len(spec.us)*len(spec.vs) > maxExpandedSources {
			warns = append(warns, fmt.Sprintf("扩展光源 %s 展开后光源总数超过 %d，已忽略该光源",
				s.ResolvedLabel(i), maxExpandedSources))
			continue
		}
		base := s.ResolvedID(i)
		pos, dir := s.Position(), s.Direction()
		dist := emitterStandoff(scene, pos, dir)
		uAxis, vAxis := emitterAxes(dir)
		for j, v := range spec.vs {
			for k, u := range spec.us {
				seen++
				// The point itself sits on the emitter plane; the source is
				// placed R in front of it with curvature radius R, so the
				// sampled wavefront is centred exactly on the point whatever
				// R is (see ExtendedSourceParams.Standoff).
				r := extendedStandoff(spec.p.Standoff, dist)
				pt := pos.Add(uAxis.Scale(u)).Add(vAxis.Scale(v))
				src := pt.Add(dir.Scale(r))
				params := map[string]any{
					"power":  spec.p.Power * spec.w[j][k],
					"radius": r,
				}
				if pbool(s.Params, "converging", false) {
					params["converging"] = true
				}
				out = append(out, SourceSpec{
					ID:         fmt.Sprintf("%s_p%d", base, seen),
					Label:      fmt.Sprintf("%s·点%d", s.ResolvedLabel(i), seen),
					Type:       "spherical",
					Params:     params,
					Pos:        &src,
					Dir:        &dir,
					Wavelength: s.EffectiveWavelength(0),
					// Group stays empty: every point is its own incoherent
					// unit, which is the definition of an incoherent emitter.
				})
			}
		}
	}
	return out, warns
}

// maxExpandedSources bounds the total number of point sources one scene may
// end up with after the expansion (across all extended emitters).
const maxExpandedSources = 8192

// hasExtendedSource reports whether any source needs the incoherent-emitter
// expansion, so scenes without one skip the copy entirely.
func hasExtendedSource(sources []SourceSpec) bool {
	for i := range sources {
		if sources[i].Type == "extended" {
			return true
		}
	}
	return false
}

// extendedStandoff picks the sampling distance R: the requested one when
// given, otherwise 0.95 of the distance to the first element in front of the
// emitter. The cap keeps the virtual point behind the sampling plane (and the
// source in front of the first element it illuminates), so the routing is
// unchanged. R never affects the propagated field: the wavefront stays centred
// on the emitter point for every R > 0, only the sampling changes.
func extendedStandoff(requested, dist float64) float64 {
	r := requested
	if r <= 0 {
		if dist > 0 {
			r = 0.95 * dist
		} else {
			r = extendedDefaultStandoff
		}
	}
	if dist > 0 && r > 0.95*dist {
		r = 0.95 * dist
	}
	return r
}

// emitterStandoff returns the distance from the emitter plane to the first
// component in front of it along the emission direction (0 when none).
func emitterStandoff(scene *SceneSpec, from, dir Vec3) float64 {
	if scene == nil {
		return 0
	}
	best := math.Inf(1)
	for i := range scene.Components {
		g := resolveGeom(&scene.Components[i])
		denom := dir.Dot(g.n)
		if math.Abs(denom) < sceneParallelCos {
			continue
		}
		if t := g.pos.Sub(from).Dot(g.n) / denom; t > sceneEpsilon && t < best {
			best = t
		}
	}
	if math.IsInf(best, 1) {
		return 0
	}
	return best
}

// emitterAxes returns an orthonormal pair (u, v) spanning the plane
// perpendicular to the emission direction d: u horizontal, v vertical.
func emitterAxes(d Vec3) (Vec3, Vec3) {
	up := v3(0, 1, 0)
	if math.Abs(d.Dot(up)) > 0.999 {
		up = v3(1, 0, 0)
	}
	u := up.Cross(d).Unit()
	v := d.Cross(u).Unit()
	return u, v
}
