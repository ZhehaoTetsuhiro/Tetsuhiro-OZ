package geo

import (
	"fmt"
	"math"
	"sort"
)

// ---------------------------------------------------------------------------
// Ghost analysis: every double-reflection path of a system
// ---------------------------------------------------------------------------
//
// GhostAnalysis enumerates the double-reflection paths of a system (see
// GhostPairs), estimates the relative energy R_I·R_J of each of them, keeps
// the paths above the requested energy threshold, sorts them by energy
// (descending, ties by surface index pair) and reports for each path
//
//   - the paraxial landing point on the image plane, from the reversed ABCD
//     chain of the path (ParaxialGhostLanding),
//   - the real-ray spot of a small-aperture sampling about the reference
//     (pupil-centre) ray: centroid, RMS and GEO radius (ReplayGhost).
//
// Everything is deterministic: equal inputs give bit-identical reports (no map
// iteration, no randomness, no time).
//
// The analysis inherits the limitations of the ghost model (see the GhostPath
// documentation): the relative energy is a normal-incidence Fresnel product
// without coating phase and without polarization, exactly two reflections are
// followed, and mirror or coordinate-break systems are reported as unsupported
// paths rather than analysed.

// GhostOptions configures a ghost analysis. Field, Pupil and Points describe
// the small-aperture sampling exactly as they do for SpotDiagram: the launch
// rays are LaunchRay(Field, Pupil, point, LambdaNm) for every sampling point,
// so an existing pattern (Hexapolar, Rings, Grid) scaled to a small pupil
// radius can be used directly.
type GhostOptions struct {
	// Field is the object point (an infinite field or a finite object).
	Field Field
	// Pupil locates the sampled bundle: Z the plane the rays start on, Radius
	// the sampling semi-diameter (mm, positive) and Offset the transverse
	// centre of the sampling. For a coaxial system the entrance pupil
	// (Paraxial.EntrancePupilZ) with a small Radius is the natural choice; the
	// reference ray is the ray through the centre of the sampling.
	Pupil Pupil
	// Points are the normalized sampling points; nil means Hexapolar(1), the
	// centre plus its six-point ring.
	Points []PupilPoint
	// ImageZ is the axial position (mm) of the image plane the ghost spots are
	// evaluated on.
	ImageZ float64
	// LambdaNm is the wavelength (nm).
	LambdaNm float64
	// MinRelative drops every path whose relative energy R_I·R_J is smaller
	// than this value (0 keeps all of them, including the unsupported ones).
	MinRelative float64
	// MaxPaths keeps at most this many paths after sorting (0 keeps all).
	MaxPaths int
	// Reflectance optionally supplies the intensity reflectance of every
	// surface (len = NumSurfaces; NaN entries fall back to the normal-incidence
	// Fresnel value). It is the "given per-surface R" variant of the energy
	// estimate, e.g. for coated surfaces.
	Reflectance []float64
}

// GhostSpot is the real-ray ghost spot on the image plane.
type GhostSpot struct {
	N        int     // rays that reached the image plane
	NTotal   int     // rays launched
	Centroid Vec3    // centroid of the landings (zero when N == 0)
	RMS      float64 // RMS radius about the centroid (mm)
	GEO      float64 // largest distance from the centroid (mm)
}

// String implements fmt.Stringer.
func (g GhostSpot) String() string {
	return fmt.Sprintf("spot[n=%d/%d centroid=%s rms=%.6g geo=%.6g]",
		g.N, g.NTotal, g.Centroid, g.RMS, g.GEO)
}

// GhostResult is the outcome of one ghost path.
type GhostResult struct {
	Pair     GhostPair
	Energy   float64 // relative energy R_I·R_J
	ReflectI float64 // intensity reflectance of the first reflection site
	ReflectJ float64 // intensity reflectance of the second reflection site
	// Paraxial is the paraxial landing point (mm, global frame, z = ImageZ);
	// its transverse coordinates are NaN when the first-order engine does not
	// apply (a decentred or tilted system) or when the path was skipped.
	Paraxial Vec3
	// Chief is the real-ray landing point of the reference (pupil-centre) ray
	// and ChiefStatus its status; the centroid of the whole sampling is
	// Spot.Centroid.
	Chief       Vec3
	ChiefStatus Status
	Spot        GhostSpot
	// Note is non-empty when the path could not be analysed (a mirror or
	// coordinate-break element in the path, a reflectance outside [0,1], a
	// non-finite energy). Such paths carry Energy 0 and are ranked last.
	Note string
}

// Skipped reports whether the path was not analysed.
func (g GhostResult) Skipped() bool { return g.Note != "" }

// String implements fmt.Stringer.
func (g GhostResult) String() string {
	if g.Skipped() {
		return fmt.Sprintf("ghost[%s: %s]", g.Pair, g.Note)
	}
	return fmt.Sprintf("ghost[%s E=%.6g RI=%.6g RJ=%.6g paraxial=%s chief=%s %s]",
		g.Pair, g.Energy, g.ReflectI, g.ReflectJ, g.Paraxial, g.Chief, g.Spot)
}

// GhostReport is the result of a ghost analysis: the enumerated path count and
// the reported paths.
type GhostReport struct {
	Lambda     float64       // wavelength (nm)
	ImageZ     float64       // image plane (mm)
	Enumerated int           // number of double-reflection paths of the system
	Paths      []GhostResult // kept paths, sorted by energy (descending)
}

// String implements fmt.Stringer.
func (r *GhostReport) String() string {
	s := fmt.Sprintf("ghost report[λ=%.6gnm z=%.6gmm paths %d/%d]", r.Lambda, r.ImageZ, len(r.Paths), r.Enumerated)
	for i := range r.Paths {
		s += "\n  " + r.Paths[i].String()
	}
	return s
}

// GhostResultLess orders ghost results: relative energy descending, ties by
// the surface index pair (I ascending, then J ascending). It is the ordering
// of GhostReport.Paths and a strict weak ordering (energies of a report are
// finite, unsupported paths carry 0).
func GhostResultLess(a, b GhostResult) bool {
	if a.Energy != b.Energy {
		return a.Energy > b.Energy
	}
	if a.Pair.I != b.Pair.I {
		return a.Pair.I < b.Pair.I
	}
	return a.Pair.J < b.Pair.J
}

// ghostCandidate is one enumerated path with the energy the ordering uses.
type ghostCandidate struct {
	result GhostResult
	energy float64
}

// GhostAnalysis analyses every double-reflection path of the system (see
// GhostOptions, GhostResult and GhostReport). It returns an error only for a
// system or an option set that cannot be analysed at all; individual paths the
// model does not cover are reported with a Note.
func (s *System) GhostAnalysis(opt GhostOptions) (*GhostReport, error) {
	if err := s.Validate(); err != nil {
		return nil, err
	}
	if !finite(opt.ImageZ) {
		return nil, fmt.Errorf("geo: ghost analysis: image plane z is not finite")
	}
	if !finitePositive(opt.LambdaNm) {
		return nil, fmt.Errorf("geo: ghost analysis: wavelength %g is not positive and finite", opt.LambdaNm)
	}
	if !finitePositive(opt.Pupil.Radius) {
		return nil, fmt.Errorf("geo: ghost analysis: pupil radius %g is not positive and finite", opt.Pupil.Radius)
	}
	if !finite(opt.Pupil.Z) || !finite(opt.Pupil.Offset.X) || !finite(opt.Pupil.Offset.Y) {
		return nil, fmt.Errorf("geo: ghost analysis: pupil plane or offset is not finite")
	}
	if opt.MinRelative < 0 || math.IsNaN(opt.MinRelative) {
		return nil, fmt.Errorf("geo: ghost analysis: minimum relative energy %g is negative", opt.MinRelative)
	}
	if opt.MaxPaths < 0 {
		return nil, fmt.Errorf("geo: ghost analysis: MaxPaths %d is negative", opt.MaxPaths)
	}
	pts := opt.Points
	if len(pts) == 0 {
		pts = Hexapolar(1)
	}
	for _, p := range pts {
		if !finite(p.PX) || !finite(p.PY) {
			return nil, fmt.Errorf("geo: ghost analysis: sampling point %s is not finite", p)
		}
	}
	// The reference ray is the ray through the centre of the sampling: the ray
	// the paraxial prediction belongs to and the one traced as GhostResult.Chief.
	refRay := LaunchRay(opt.Field, opt.Pupil, PupilPoint{}, opt.LambdaNm)
	if !refRay.Valid() {
		return nil, fmt.Errorf("geo: ghost analysis: the reference ray of the sampling is degenerate")
	}
	paraxialRay, ok := GhostParaxialOfRay(refRay, s.VertexZ(0))
	if !ok {
		return nil, fmt.Errorf("geo: ghost analysis: the reference ray does not travel towards +z")
	}

	// The relative energies are algebraic: the ordering, the threshold and the
	// top-N cut do not depend on any trace.
	cands := make([]ghostCandidate, 0, len(GhostPairs(s.NumSurfaces())))
	for _, pr := range GhostPairs(s.NumSurfaces()) {
		res := GhostResult{Pair: pr, Paraxial: V(math.NaN(), math.NaN(), opt.ImageZ)}
		energy, ri, rj, err := s.GhostEnergy(pr, opt.LambdaNm, opt.Reflectance)
		switch {
		case err != nil:
			res.Note = err.Error()
		case !finite(energy):
			res.Note = fmt.Sprintf("geo: ghost energy of %s is not finite", pr)
		default:
			if _, perr := s.GhostPathOf(pr); perr != nil {
				res.Note = perr.Error()
				break
			}
			res.Energy, res.ReflectI, res.ReflectJ = energy, ri, rj
		}
		if res.Skipped() {
			energy = 0
		}
		cands = append(cands, ghostCandidate{result: res, energy: energy})
	}
	// Deterministic order: energy descending, then (I, J) ascending.
	sort.SliceStable(cands, func(a, b int) bool {
		if cands[a].energy != cands[b].energy {
			return cands[a].energy > cands[b].energy
		}
		return GhostResultLess(cands[a].result, cands[b].result)
	})

	report := &GhostReport{Lambda: opt.LambdaNm, ImageZ: opt.ImageZ, Enumerated: len(cands)}
	for _, c := range cands {
		if c.energy < opt.MinRelative {
			continue
		}
		if opt.MaxPaths > 0 && len(report.Paths) >= opt.MaxPaths {
			break
		}
		res := c.result
		if !res.Skipped() {
			if err := s.fillGhostResult(&res, paraxialRay, opt, pts, refRay); err != nil {
				res.Note = err.Error()
			}
		}
		report.Paths = append(report.Paths, res)
	}
	return report, nil
}

// fillGhostResult completes one result: the paraxial landing point of the
// reference ray and the real-ray spot of the sampling plus the reference ray.
func (s *System) fillGhostResult(res *GhostResult, paraxialRay GhostParaxial, opt GhostOptions, pts []PupilPoint, refRay Ray) error {
	path, err := s.GhostPathOf(res.Pair)
	if err != nil {
		return err
	}
	res.Paraxial = V(math.NaN(), math.NaN(), opt.ImageZ)
	if px, perr := s.ParaxialGhostLanding(res.Pair, paraxialRay, opt.ImageZ, opt.LambdaNm); perr == nil {
		res.Paraxial = px
	}
	chief := s.ReplayGhost(path, refRay, opt.ImageZ)
	res.Chief, res.ChiefStatus = chief.Landing, chief.Status

	// One pass over the rays; the reference ray is counted once, whether or
	// not the sampling itself contains the pupil centre.
	landings := make([]Vec3, 0, len(pts)+1)
	refInSet := false
	for _, p := range pts {
		if p.PX == 0 && p.PY == 0 {
			refInSet = true
		}
	}
	res.Spot.NTotal = len(pts)
	if !refInSet {
		res.Spot.NTotal++
	}
	for _, p := range pts {
		hit := s.ReplayGhost(path, LaunchRay(opt.Field, opt.Pupil, p, opt.LambdaNm), opt.ImageZ)
		if hit.OK() {
			landings = append(landings, hit.Landing)
		}
	}
	if !refInSet && chief.OK() {
		landings = append(landings, chief.Landing)
	}
	if len(landings) == 0 {
		return nil
	}
	sum := Zero
	for _, v := range landings {
		sum = sum.Add(v)
	}
	res.Spot.N = len(landings)
	res.Spot.Centroid = sum.Div(float64(len(landings)))
	var rms2, geo float64
	for _, v := range landings {
		d := v.Sub(res.Spot.Centroid).Norm()
		rms2 += d * d
		if d > geo {
			geo = d
		}
	}
	res.Spot.RMS = math.Sqrt(rms2 / float64(len(landings)))
	res.Spot.GEO = geo
	return nil
}
