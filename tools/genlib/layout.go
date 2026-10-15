// The example catalogue machinery of library/geo: thin-lens seeds, layout
// expansion, paraxial scaling, automatic apertures, the DLS polish wrapper,
// and the library container with its failure records. Every prescription in
// the library is produced by one of the solvers here (closed form, bounded
// root search, or the package's own least-squares polisher), never by
// hand-typed numbers.
package main

import (
	"fmt"
	"math"
	"strings"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/tools/libcheck/measure"
)

// ---------- generation records ----------

// The records of one generation run: every candidate that was skipped with
// its reason, every candidate merged into an identical earlier design, and
// every file renamed to keep file names unique. main writes them to
// library/geo/FAILURES.md.
var (
	buildFails []string // construction failures: the candidate never reached the library
	candFails  []string // registered candidates skipped at measurement or assertion time
	merges     []string // candidates merged into an identical earlier design
	renames    []string // files renamed to keep names unique
)

func record(dst *[]string, format string, a ...interface{}) {
	msg := strings.TrimRight(fmt.Sprintf(format, a...), "\n")
	*dst = append(*dst, msg)
	fmt.Println("  ! " + msg)
}

// failf records a construction failure: a layout or solve that produced no
// candidate. The message names the family and the parameters.
func failf(format string, a ...interface{}) { record(&buildFails, format, a...) }

// skipf records a registered candidate that failed measurement, the
// library assertions or the description check, and was not written.
func skipf(format string, a ...interface{}) { record(&candFails, format, a...) }

// ---------- library container ----------

type item struct {
	file string
	p    *lens.Prescription
	note string
}

type lib struct {
	items []item
	seen  map[string]bool   // file names in use
	names map[string]bool   // example names in use
	fps   map[string]string // parameter fingerprint -> file that holds it
	byCat map[string]int
	// nir marks the examples designed for the near infrared.
	nir map[string]bool
}

func newLib() *lib {
	return &lib{seen: map[string]bool{}, names: map[string]bool{}, fps: map[string]string{},
		byCat: map[string]int{}, nir: map[string]bool{}}
}

// add registers one candidate. It never panics:
//   - a candidate whose parameters equal an earlier one is merged into it;
//   - a file name already in use gets a readable "-altN" suffix;
//   - a name already in use gets a "（变体N）" suffix.
func (l *lib) add(file string, p *lens.Prescription, cat string, diff int, tg []string, note string) {
	measure.Canonical(p)
	fp, err := measure.Fingerprint(p)
	if err != nil {
		failf("%s: fingerprint: %v", file, err)
		return
	}
	if prev, ok := l.fps[fp]; ok {
		record(&merges, "%s 与 %s 参数指纹相同（同一光学设计），合并", file, prev)
		return
	}
	name := file
	if l.seen[name] {
		base := strings.TrimSuffix(file, ".json")
		for k := 2; ; k++ {
			cand := fmt.Sprintf("%s-alt%d.json", base, k)
			if !l.seen[cand] {
				name = cand
				break
			}
		}
		renames = append(renames, fmt.Sprintf("%s 改名为 %s（同名文件，参数不同）", file, name))
	}
	l.seen[name] = true
	l.fps[fp] = name
	p.Category = cat
	p.Difficulty = diff
	p.Tags = tg
	if p.Name == "" {
		p.Name = name
	}
	if l.names[p.Name] {
		base := p.Name
		for k := 2; ; k++ {
			cand := fmt.Sprintf("%s（变体%d）", base, k)
			if !l.names[cand] {
				p.Name = cand
				break
			}
		}
	}
	l.names[p.Name] = true
	l.items = append(l.items, item{file: name, p: p, note: note})
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
// meniscus.
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

// grp is one element of a layout: curvature ratios (from the seed functions
// above; the whole layout is scaled by the paraxial solve, so their units do
// not matter), the glass after every surface, the element thickness and the
// air gap that follows it.
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
// the used beam, and solves the radii for the target focal length.
//
// The result is frozen in three steps — thicknesses to 0.1 µm, the radii
// re-solved on that geometry, the radii to 0.1 nm — and the frozen
// prescription must reproduce the target EFL to a relative residual of 1e-6.
// Otherwise an error is returned and the candidate is skipped.
func finishLens(p *lens.Prescription, pupilSemi, targetEFL float64) error {
	if pupilSemi > 0 {
		p.Cards[0].SemiAperMM = r3(pupilSemi)
	}
	// Exact similarity step: scaling every radius and every thickness by a
	// multiplies the EFL by exactly a, so the seed lands on the target in one
	// step. Only applied when the seed is already in the right power class.
	if targetEFL != 0 {
		if par, err := p.Paraxial(geo.Wavelengthd); err == nil {
			if e := par.EFL; !math.IsNaN(e) && !math.IsInf(e, 0) && e != 0 {
				if a := targetEFL / e; a > 1e-3 && a < 1e3 {
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
		if targetEFL != 0 {
			if err := solveEFL(p, radiiIdx(p.Cards), targetEFL); err != nil {
				return err
			}
		}
		if err := autoApertures(p, pupilSemi); err != nil {
			return err
		}
	}
	// Freeze the geometry, re-solve the radii on it, freeze the radii.
	for i := range p.Cards {
		p.Cards[i].ThicknessMM = r4(p.Cards[i].ThicknessMM)
	}
	if targetEFL != 0 {
		if err := solveEFL(p, radiiIdx(p.Cards), targetEFL); err != nil {
			return fmt.Errorf("final solve: %w", err)
		}
		for i := range p.Cards {
			if !math.IsInf(p.Cards[i].RadiusMM, 0) {
				p.Cards[i].RadiusMM = r7(p.Cards[i].RadiusMM)
			}
		}
		par, err := p.Paraxial(geo.Wavelengthd)
		if err != nil {
			return fmt.Errorf("frozen prescription: %w", err)
		}
		if rel := math.Abs(par.EFL-targetEFL) / math.Abs(targetEFL); !(rel <= 1e-6) {
			return fmt.Errorf("EFL residual %.3g after freezing (target %g, got %g)", rel, targetEFL, par.EFL)
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

// spotOnOff returns the d-line on-axis and off-axis spot RMS of a refractive
// prescription (the off-axis field is the one of the measurement).
func spotOnOff(p *lens.Prescription, fieldDeg float64) (on, off float64, err error) {
	samples := lens.UnitDiskSamples(3, 6)
	on, _, _, err = lens.SpotRMS(p, geo.Wavelengthd, lens.Field{}, samples)
	if err != nil {
		return 0, 0, err
	}
	f := lens.Field{TanY: 0.02 * math.Abs(p.ObjectZMM)}
	if math.IsInf(p.ObjectZMM, -1) {
		f = lens.Field{TanY: math.Tan(fieldDeg * math.Pi / 180)}
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

// trimNum formats x with dec decimals and drops the trailing zeros of the
// fraction only. (Stripping zeros of the integer part as well turned 100
// into "1" and 1000 into "1": that collision produced the duplicate file
// name mirror-parabola-f1-f4.json.)
func trimNum(x float64, dec int) string {
	s := fmt.Sprintf("%.*f", dec, x)
	if strings.Contains(s, ".") {
		s = strings.TrimRight(s, "0")
		s = strings.TrimSuffix(s, ".")
	}
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
