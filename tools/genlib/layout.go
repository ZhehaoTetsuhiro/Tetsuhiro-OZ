// The example catalogue machinery of library/geo: thin-lens seeds, layout
// expansion, paraxial scaling, automatic apertures and the DLS polish
// wrapper. Every prescription in the library is produced by one of the two
// solvers here (closed form, or the package's own least-squares polisher),
// never by hand-typed numbers.
package main

import (
	"fmt"
	"math"
	"strings"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// ---------- library container ----------

type item struct {
	file string
	p    *lens.Prescription
	note string
}

type lib struct {
	items []item
	seen  map[string]bool
	byCat map[string]int
	// nir marks the examples that are designed for the near infrared: their
	// 1064/1550 nm focal difference is measured as well.
	nir map[string]bool
}

func newLib() *lib {
	return &lib{seen: map[string]bool{}, byCat: map[string]int{}, nir: map[string]bool{}}
}

func (l *lib) add(file string, p *lens.Prescription, cat string, diff int, tg []string, note string) {
	if l.seen[file] {
		panic("duplicate file name " + file)
	}
	l.seen[file] = true
	p.Category = cat
	p.Difficulty = diff
	p.Tags = tg
	if p.Name == "" {
		p.Name = file
	}
	l.items = append(l.items, item{file: file, p: p, note: note})
	l.byCat[cat]++
}

// ---------- thin-lens seeds ----------

func idxOf(glass string) float64 {
	g, ok := geo.Glass(glass)
	if !ok {
		panic("unknown glass " + glass)
	}
	return g.N(geo.Wavelengthd)
}

func abbeOf(glass string) float64 {
	g, _ := geo.Glass(glass)
	return g.AbbeNumber()
}

// seedFlatPCX is the curvature pair of a plano-convex lens whose flat face
// meets the light first.
func seedFlatPCX(glass string, f float64) []float64 {
	n := idxOf(glass)
	return []float64{0, -1 / ((n - 1) * f)}
}

// seedConvPCX is the curvature pair of a plano-convex lens whose convex face
// meets the light first (the best form for an infinite conjugate).
func seedConvPCX(glass string, f float64) []float64 {
	n := idxOf(glass)
	return []float64{1 / ((n - 1) * f), 0}
}

// seedBic is the symmetric biconvex pair.
func seedBic(glass string, f float64) []float64 {
	n := idxOf(glass)
	return []float64{1 / (2 * (n - 1) * f), -1 / (2 * (n - 1) * f)}
}

// seedMen is the thin-lens seed of a meniscus whose rear curvature is 0.55
// times its front curvature: the sign of f picks a positive or a negative
// meniscus. It is a real curvature pair (mm⁻¹), not a ratio, so that a
// meniscus can be mixed with other elements inside one layout.
func seedMen(glass string, f float64) []float64 {
	n := idxOf(glass)
	k := 1 / (0.45 * (n - 1) * f)
	return []float64{k, 0.55 * k}
}

func seedMenP(glass string, f float64) []float64 { return seedMen(glass, f) }
func seedMenN(glass string, f float64) []float64 { return seedMen(glass, f) }

// seedAchromat is the thin-lens curvature triple of a cemented crown+flint
// doublet of focal length f: the crown takes its half of the power split by
// the Abbe numbers, and its two surfaces share the power equally.
func seedAchromat(crown, flint string, f float64) []float64 {
	n1, n2 := idxOf(crown), idxOf(flint)
	V1, V2 := abbeOf(crown), abbeOf(flint)
	phi := 1 / f
	phi1 := phi * V1 / (V1 + V2)
	phi2 := phi - phi1
	cA := phi1 / (2 * (n1 - 1))
	cB := -cA
	cC := cB - phi2/(n2-1)
	return []float64{cA, cB, cC}
}

// ---------- layout machinery ----------

// grp is one element of a layout: curvature *ratios* (from the seed
// functions above — the whole layout is scaled by the paraxial solve, so
// their units do not matter), the glass after every surface, the element
// thickness and the air gap that follows it.
type grp struct {
	name    string
	curv    []float64
	glasses []string // one entry per surface: medium *after* that surface
	thick   float64  // axial thickness of the element
	after   float64  // air gap after the element
}

// cardsOf expands a layout into prescription cards behind the aperture-stop
// plane, which sits at z = 0 and is `lead` in front of the first surface.
func cardsOf(gs []grp, lead float64) []lens.Card {
	cards := []lens.Card{{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: lead, SemiAperMM: 1}}
	for gi, g := range gs {
		n := len(g.curv)
		for i, c := range g.curv {
			r := math.Inf(1)
			if c != 0 {
				r = 1 / c
			}
			t := g.thick / float64(n-1)
			if i == n-1 {
				t = g.after
				if gi == len(gs)-1 {
					t = math.Max(10, g.after) // trailing air: ignored by System()
				}
			}
			cards = append(cards, lens.Card{
				Name:        g.name,
				RadiusMM:    r,
				ThicknessMM: t,
				Glass:       g.glasses[i],
			})
		}
	}
	return cards
}

// radiiIdx lists the card indices that carry a finite curvature.
func radiiIdx(cards []lens.Card) []int {
	var out []int
	for i, c := range cards {
		if !math.IsInf(c.RadiusMM, 0) && c.RadiusMM != 0 {
			out = append(out, i)
		}
	}
	return out
}

// buildLayout assembles a prescription from a layout, scales all curvatures
// so that the paraxial EFL equals efl, rounds the geometry and fills the
// apertures so that no used ray is vignetted.
func buildLayout(name string, gs []grp, lead, efl, pupilSemi, objZ float64, objTags ...string) (*lens.Prescription, error) {
	p := &lens.Prescription{Name: name, ObjectZMM: objZ, Stop: 0, Cards: cardsOf(gs, lead)}
	if err := finishLens(p, pupilSemi, efl); err != nil {
		return nil, err
	}
	if len(objTags) > 0 {
		p.Tags = objTags
	}
	return p, nil
}

// finishLens rounds the geometry, guarantees manufacturing room between
// neighbouring surfaces, gives every non-stop card an aperture that clears
// the used beam and re-solves the focal length after those edits. It is
// idempotent enough to be called again after a polish.
func finishLens(p *lens.Prescription, pupilSemi, targetEFL float64) error {
	if pupilSemi > 0 {
		p.Cards[0].SemiAperMM = r3(pupilSemi)
	}
	// Exact similarity step: scaling every radius *and* every thickness by a
	// multiplies the EFL by exactly a (the paraxial equations are
	// similarity-covariant), so the seed layout lands on the requested focal
	// length in one step. The loops below then only repair rounding and
	// manufacturing clearance.
	if targetEFL != 0 {
		if par, err := p.Paraxial(geo.Wavelengthd); err == nil {
			if e := par.EFL; !math.IsNaN(e) && !math.IsInf(e, 0) && e != 0 {
				a := targetEFL / e
				// Only when the seed is already in the right power class:
				// a near-afocal seed would otherwise collapse to nothing.
				if a > 1e-3 && a < 1e3 {
					for i := range p.Cards {
						if !math.IsInf(p.Cards[i].RadiusMM, 0) {
							p.Cards[i].RadiusMM *= a
						}
						p.Cards[i].ThicknessMM *= a
					}
				}
			}
		}
	}
	for pass := 0; pass < 3; pass++ {
		for i := range p.Cards {
			if !math.IsInf(p.Cards[i].RadiusMM, 0) {
				p.Cards[i].RadiusMM = r4(p.Cards[i].RadiusMM)
			}
		}
		if err := clearance(p, pupilSemi); err != nil {
			return err
		}
		if err := solveEFL(p, radiiIdx(p.Cards), targetEFL); err != nil {
			return err
		}
		if err := autoApertures(p, pupilSemi); err != nil {
			return err
		}
	}
	for i := range p.Cards {
		p.Cards[i].ThicknessMM = r3(p.Cards[i].ThicknessMM)
		if !math.IsInf(p.Cards[i].RadiusMM, 0) {
			p.Cards[i].RadiusMM = r4(p.Cards[i].RadiusMM)
		}
	}
	if err := p.Validate(); err != nil {
		return err
	}
	_, err := p.System()
	return err
}

// clearance grows the axial room of surfaces whose sag would leave no gap
// (and of a concave first surface, whose rim would sit in front of the
// aperture stop).
func clearance(p *lens.Prescription, pupilSemi float64) error {
	for i := 1; i < len(p.Cards); i++ {
		semi := pupilSemi
		if p.Cards[i].SemiAperMM > 0 {
			semi = math.Max(semi, p.Cards[i].SemiAperMM)
		}
		s := sagConic(semi, p.Cards[i].RadiusMM, p.Cards[i].ConicK)
		if math.IsNaN(s) {
			continue
		}
		if i == 1 && s < 0 && -s+0.5 > p.Cards[0].ThicknessMM {
			p.Cards[0].ThicknessMM = r3(math.Max(2, -s*1.2+0.5))
		}
		if i+1 < len(p.Cards) && !p.Cards[i].Mirror {
			s2 := sagConic(semi, p.Cards[i+1].RadiusMM, p.Cards[i+1].ConicK)
			if !math.IsNaN(s2) {
				if edge := p.Cards[i].ThicknessMM + s2 - s; edge < 0.6 {
					p.Cards[i].ThicknessMM = r3(p.Cards[i].ThicknessMM + (0.6 - edge))
				}
			}
		}
	}
	return nil
}

// autoApertures traces on-axis marginal rays through the entrance pupil and
// sets every card's clear semi-aperture to clear the used beam with margin.
// The aperture-stop card keeps its exact value: it defines the F-number.
func autoApertures(p *lens.Prescription, pupilSemi float64) error {
	if !(pupilSemi > 0) {
		return nil
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil || math.IsNaN(par.EntrancePupilZ) || math.IsInf(par.EntrancePupilZ, 0) {
		return nil
	}
	sys, err := p.System()
	if err != nil {
		return err
	}
	maxY := make([]float64, sys.NumSurfaces())
	launch := func(sx, sy float64) {
		var r geo.Ray
		if math.IsInf(p.ObjectZMM, -1) {
			r = geo.NewRay(geo.V(sx*pupilSemi, sy*pupilSemi, par.EntrancePupilZ), geo.V(0, 0, 1), geo.Wavelengthd)
		} else {
			r = geo.RayThrough(geo.V(0, 0, p.ObjectZMM), geo.V(sx*pupilSemi, sy*pupilSemi, par.EntrancePupilZ), geo.Wavelengthd)
		}
		res := sys.Trace(r)
		if !res.OK() {
			return
		}
		for _, st := range res.Steps {
			if y := math.Abs(st.Point.Y); y > maxY[st.Surface] {
				maxY[st.Surface] = y
			}
		}
	}
	const k = 0.999
	for _, s := range [][2]float64{{0, -k}, {0, k}, {-k, 0}, {k, 0}, {0.7 * k, 0.7 * k}, {-0.7 * k, -0.7 * k}} {
		launch(s[0], s[1])
	}
	for i := 1; i < len(p.Cards); i++ {
		need := 1.25*maxY[i] + 0.4
		if need > p.Cards[i].SemiAperMM {
			p.Cards[i].SemiAperMM = r3(need)
		}
	}
	return nil
}

// ---------- quality helpers ----------

// spotOnOff returns the d-line on-axis and off-axis spot RMS of a refractive
// prescription.
func spotOnOff(p *lens.Prescription, fieldDeg float64) (on, off float64, err error) {
	samples := lens.UnitDiskSamples(3, 6)
	on, _, _, err = lens.SpotRMS(p, geo.Wavelengthd, lens.Field{}, samples)
	if err != nil {
		return 0, 0, err
	}
	f := lens.Field{}
	if math.IsInf(p.ObjectZMM, -1) {
		f = lens.Field{TanY: math.Tan(fieldDeg * math.Pi / 180)}
	} else {
		f = lens.Field{TanY: 0.02 * math.Abs(p.ObjectZMM)}
	}
	off, _, _, err = lens.SpotRMS(p, geo.Wavelengthd, f, samples)
	if err != nil {
		return on, 0, err
	}
	return on, off, nil
}

// polish runs the packaged damped-least-squares polisher on the curvatures
// (and the first few thicknesses) against the spot diagrams of two fields,
// the focal length and the axial colour.
func polish(p *lens.Prescription, efl, fieldDeg float64, maxIter int) (*lens.Prescription, *lens.Result, error) {
	specs := make([]lens.VarSpec, 0, 12)
	for _, i := range radiiIdx(p.Cards) {
		specs = append(specs, lens.VarSpec{Card: i, Kind: lens.VarRadius})
	}
	for i := 1; i < len(p.Cards)-1 && i <= 3; i++ {
		specs = append(specs, lens.VarSpec{Card: i, Kind: lens.VarThickness})
	}
	work := p.Clone()
	cur, err := work.Vars(specs)
	if err != nil {
		return nil, nil, err
	}
	rLim := 20 * math.Abs(efl)
	vars := make([]lens.Variable, len(specs))
	for i, s := range specs {
		lo, hi := -rLim, rLim
		if s.Kind == lens.VarThickness {
			lo, hi = 0.4, math.Max(2, 0.12*math.Abs(efl))
		}
		vars[i] = lens.Variable{Name: s.String(), Value: cur[i], Lo: lo, Hi: hi}
	}
	samples := lens.UnitDiskSamples(3, 6)
	fields := []lens.Field{{}, {}}
	if math.IsInf(p.ObjectZMM, -1) {
		fields[1] = lens.Field{TanY: math.Tan(fieldDeg * math.Pi / 180)}
	} else {
		fields[1] = lens.Field{TanY: 0.02 * math.Abs(p.ObjectZMM)}
	}
	target := math.Max(2e-3, math.Abs(efl)/2000)
	tol := math.Max(1e-3, math.Abs(efl)/20000)
	const fail = 10.0
	resid := func(x []float64) ([]float64, error) {
		q := work.Clone()
		if err := q.SetVars(specs, x); err != nil {
			return nil, err
		}
		var out []float64
		for _, f := range fields {
			r, err := lens.SpotResiduals(q, geo.Wavelengthd, []lens.Field{f}, samples, target)
			if err != nil {
				for i := 0; i < 2*len(samples); i++ {
					out = append(out, fail)
				}
				continue
			}
			out = append(out, r...)
		}
		if e, err := lens.EFLResidual(q, geo.Wavelengthd, efl, tol); err == nil {
			out = append(out, e)
		} else {
			out = append(out, fail)
		}
		if e, err := lens.ChromaticResidual(q, 0.05); err == nil {
			out = append(out, e)
		} else {
			out = append(out, fail)
		}
		return out, nil
	}
	if maxIter <= 0 {
		maxIter = 120
	}
	res, err := lens.Minimize(lens.Problem{Variables: vars, Residual: resid, MaxIter: maxIter})
	if err != nil {
		return nil, nil, err
	}
	final := work.Clone()
	if err := final.SetVars(specs, res.X); err != nil {
		return nil, nil, err
	}
	if err := final.Validate(); err != nil {
		return nil, nil, err
	}
	return final, res, nil
}

// polishIfBetter keeps the polished prescription only when it improves both
// the on-axis and the off-axis spot.
func polishIfBetter(p *lens.Prescription, efl, fieldDeg float64, maxIter int) (*lens.Prescription, string, *lens.Result) {
	on0, off0, err := spotOnOff(p, fieldDeg)
	if err != nil {
		return p, "classic", nil
	}
	q, res, err := polish(p, efl, fieldDeg, maxIter)
	if err != nil {
		return p, "classic", nil
	}
	if err := finishLens(q, q.Cards[0].SemiAperMM, efl); err != nil {
		return p, "classic", nil
	}
	on1, off1, err := spotOnOff(q, fieldDeg)
	if err != nil {
		return p, "classic", nil
	}
	if on1 < on0 && off1 <= off0*1.0001 {
		return q, "optimized", res
	}
	return p, "classic", nil
}

// ---------- tag / code helpers ----------

func fnumTag(f float64) string { return "F/" + trimNum(f, 2) }

func eflTag(f float64) string { return "efl:" + trimNum(f, 1) }

func fnumCode(f float64) string { return "f" + strings.ReplaceAll(trimNum(f, 2), ".", "p") }

func trimNum(x float64, dec int) string {
	s := fmt.Sprintf("%.*f", dec, x)
	s = strings.TrimRight(s, "0")
	s = strings.TrimRight(s, ".")
	if s == "" || s == "-" {
		s = "0"
	}
	return s
}

func numCode(x float64) string {
	return strings.NewReplacer(".", "p", "-", "m").Replace(trimNum(x, 2))
}

func glassCode(name string) string {
	switch name {
	case "N-BK7":
		return "nbk7"
	case "FusedSilica", "fused silica":
		return "silica"
	case "N-SF5":
		return "nsf5"
	case "N-SK16":
		return "nsk16"
	case "F2":
		return "f2"
	}
	return strings.ToLower(strings.NewReplacer("-", "", " ", "").Replace(name))
}
