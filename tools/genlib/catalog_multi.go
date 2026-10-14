// The catalogue: three-element and multi-element camera-type objectives.
package main

import (
	"fmt"
	"math"
)

// ---------- classic layout archetypes ----------
//
// Every archetype is a list of elements with thin-lens seeds; the global
// paraxial scale solve then puts the focal length on target, after which a
// subset is polished by the package's DLS engine. All of them use a front
// aperture stop (the library convention: the stop sits in front of the first
// surface so that every example is sampled through a real entrance pupil and
// so that a concave first surface is always hit from its front side).

func layoutCooke(f float64, crown, flint string) []grp {
	return []grp{
		{name: "L1", curv: seedConvPCX(crown, 0.72*f), glasses: []string{crown, ""}, thick: 0.07 * f, after: 0.12 * f},
		{name: "L2", curv: seedBic(flint, -0.45*f), glasses: []string{flint, ""}, thick: 0.03 * f, after: 0.10 * f},
		{name: "L3", curv: seedBic(crown, 0.80*f), glasses: []string{crown, ""}, thick: 0.05 * f, after: 1.2 * f},
	}
}

func layoutCookeFlattener(f float64, crown, flint string) []grp {
	gs := layoutCooke(f, crown, flint)
	gs[2].after = 0.18 * f
	// weak negative field flattener close to the image
	gs = append(gs, grp{name: "L4", curv: seedBic(flint, -2.5*f), glasses: []string{flint, ""}, thick: 0.025 * f, after: 1.0 * f})
	return gs
}

func layoutTessar(f float64) []grp {
	return []grp{
		{name: "L1", curv: seedConvPCX("N-SK16", 0.85*f), glasses: []string{"N-SK16", ""}, thick: 0.06 * f, after: 0.10 * f},
		{name: "L2/L3", curv: seedAchromat("N-SK16", "N-SF5", 0.9*f),
			glasses: []string{"N-SK16", "N-SF5", ""}, thick: 0.07 * f, after: 0.13 * f},
		{name: "L4", curv: seedMenN("F2", -0.8*f), glasses: []string{"F2", ""}, thick: 0.05 * f, after: 1.1 * f},
	}
}

func layoutDGauss(f float64) []grp {
	return []grp{
		{name: "L1", curv: seedConvPCX("N-SK16", 0.75*f), glasses: []string{"N-SK16", ""}, thick: 0.05 * f, after: 0.02 * f},
		{name: "L2/L3", curv: seedAchromat("N-SK16", "F2", 0.8*f),
			glasses: []string{"N-SK16", "F2", ""}, thick: 0.075 * f, after: 0.07 * f},
		{name: "L4/L5", curv: seedAchromatRev("F2", "N-SK16", 0.8*f),
			glasses: []string{"F2", "N-SK16", ""}, thick: 0.075 * f, after: 0.04 * f},
		{name: "L6", curv: seedBic("N-SK16", 0.95*f), glasses: []string{"N-SK16", ""}, thick: 0.05 * f, after: 0.9 * f},
	}
}

// seedAchromatRev is the doublet solution seen from the flint side: the
// flint takes the front surface and the crown is cemented behind it.
func seedAchromatRev(flintFront, crownRear string, f float64) []float64 {
	s := seedAchromat(crownRear, flintFront, f)
	return []float64{s[2], s[1], s[0]}
}

func layoutSonnar(f float64) []grp {
	return []grp{
		{name: "L1/L2", curv: seedAchromat("N-SK16", "F2", 0.85*f),
			glasses: []string{"N-SK16", "F2", ""}, thick: 0.075 * f, after: 0.08 * f},
		{name: "L3", curv: seedBic("N-SK16", 0.55*f), glasses: []string{"N-SK16", ""}, thick: 0.06 * f, after: 0.10 * f},
		{name: "L4", curv: seedBic("F2", -0.7*f), glasses: []string{"F2", ""}, thick: 0.04 * f, after: 1.0 * f},
	}
}

func layoutPetzval(f float64) []grp {
	return []grp{
		{name: "G1", curv: seedAchromat("N-SK16", "F2", 0.55*f),
			glasses: []string{"N-SK16", "F2", ""}, thick: 0.07 * f, after: 0.35 * f},
		{name: "G2", curv: seedAchromat("N-BK7", "N-SF5", 0.5*f),
			glasses: []string{"N-BK7", "N-SF5", ""}, thick: 0.06 * f, after: 1.1 * f},
	}
}

func layoutSymmetric(f float64) []grp {
	return []grp{
		{name: "F1", curv: seedAchromat("N-SK16", "F2", 0.7*f),
			glasses: []string{"N-SK16", "F2", ""}, thick: 0.08 * f, after: 0.12 * f},
		{name: "F2", curv: seedAchromatRev("F2", "N-SK16", 0.7*f),
			glasses: []string{"F2", "N-SK16", ""}, thick: 0.08 * f, after: 1.0 * f},
	}
}

func layoutAirSpacedTriplet(f float64, crown, flint string) []grp {
	return []grp{
		{name: "L1", curv: seedBic(crown, 1.4*f), glasses: []string{crown, ""}, thick: 0.05 * f, after: 0.25 * f},
		{name: "L2", curv: seedBic(flint, -1.2*f), glasses: []string{flint, ""}, thick: 0.035 * f, after: 0.25 * f},
		{name: "L3", curv: seedBic(crown, 1.1*f), glasses: []string{crown, ""}, thick: 0.05 * f, after: 1.0 * f},
	}
}

func layoutRetrofocus(f float64) []grp {
	return []grp{
		{name: "L1", curv: seedMenN("N-SF5", -0.7*f), glasses: []string{"N-SF5", ""}, thick: 0.05 * f, after: 0.05 * f},
		{name: "L2", curv: seedBic("N-SF5", -1.1*f), glasses: []string{"N-SF5", ""}, thick: 0.035 * f, after: 0.30 * f},
		{name: "L3/L4", curv: seedAchromat("N-SK16", "F2", 0.55*f),
			glasses: []string{"N-SK16", "F2", ""}, thick: 0.075 * f, after: 0.08 * f},
		{name: "L5", curv: seedBic("N-SK16", 0.7*f), glasses: []string{"N-SK16", ""}, thick: 0.05 * f, after: 0.9 * f},
	}
}

func layoutMacro(f float64) []grp {
	return []grp{
		{name: "L1/L2", curv: seedAchromat("N-BK7", "F2", 0.9*f),
			glasses: []string{"N-BK7", "F2", ""}, thick: 0.07 * f, after: 0.25 * f},
		{name: "L3/L4", curv: seedAchromatRev("F2", "N-BK7", 0.85*f),
			glasses: []string{"F2", "N-BK7", ""}, thick: 0.07 * f, after: 0.5 * f},
	}
}

// ---------- family builders ----------

type archetype struct {
	key, zh string
	layout  func(float64) []grp
	field   float64
	diff    int
}

func multiElement(l *lib, a archetype, f, wantFnum float64, polishIt bool) {
	fnum := fnumFor(f, wantFnum)
	semi := f / (2 * fnum)
	p, err := buildLayout(fmt.Sprintf("%s f=%smm F/%s", a.zh, trimNum(f, 0), trimNum(fnum, 2)),
		a.layout(f), 3, f, semi, math.Inf(-1))
	if err != nil {
		fmt.Printf("  ! %s f=%g: %v\n", a.key, f, err)
		return
	}
	quality := "classic"
	stats := ""
	if polishIt {
		q, qual, res := polishIfBetter(p, f, a.field, 150)
		p = q
		quality = qual
		if qual == "optimized" && res != nil {
			stats = fmt.Sprintf("｜DLS 抛光：%d 次接受迭代／%d 次求值，merit=%.4g", res.Iter, res.Evals, res.Merit)
		}
	}
	if quality == "optimized" {
		p.Name = fmt.Sprintf("%s f=%smm F/%s（DLS 抛光）", a.zh, trimNum(f, 0), trimNum(fnum, 2))
	} else {
		p.Name = fmt.Sprintf("%s f=%smm F/%s", a.zh, trimNum(f, 0), trimNum(fnum, 2))
	}
	file := fmt.Sprintf("%s-f%s-%s.json", a.key, trimNum(f, 0), fnumCode(fnum))
	l.add(file, p, "multi", a.diff,
		[]string{a.key, "multi-element", fnumTag(fnum), eflTag(f), "visible",
			"quality:" + quality, fmt.Sprintf("fov%gdeg", a.field), "sweep:efl", "front-stop"},
		fmt.Sprintf("%s｜f=%smm、F/%s、%g° 视场｜经典构造（薄透镜初值 + 近轴标定）%s",
			a.zh, trimNum(f, 0), trimNum(fnum, 2), a.field, stats))
}

func buildMulti(l *lib) {
	cooke := archetype{"triplet-cooke", "三片式 Cooke 型物镜（正-负-正）", func(f float64) []grp { return layoutCooke(f, "N-SK16", "F2") }, 1, 3}
	cookeBK7 := archetype{"triplet-cooke-nbk7", "三片式 Cooke 型物镜（N-BK7/F2 玻璃对）", func(f float64) []grp { return layoutCooke(f, "N-BK7", "N-SF5") }, 1, 3}
	flattener := archetype{"triplet-flattener", "三片式 + 场镜（4 片）", func(f float64) []grp { return layoutCookeFlattener(f, "N-SK16", "F2") }, 2, 4}
	airTri := archetype{"triplet-airspaced", "三片分离式（空气间隔三片）", func(f float64) []grp { return layoutAirSpacedTriplet(f, "N-BK7", "F2") }, 1, 3}
	tessar := archetype{"tessar", "Tessar 型 4 片物镜", layoutTessar, 3, 4}
	dgauss := archetype{"dga型", "Double Gauss 型 6 片物镜", layoutDGauss, 4, 5}
	dgauss.key = "double-gauss"
	dgauss.zh = "Double Gauss 型 6 片物镜"
	sonnar := archetype{"sonnar", "Sonnar 型 5 片物镜", layoutSonnar, 3, 4}
	petzval := archetype{"petzval", "Petzval 型 4 片（两组分离）", layoutPetzval, 2, 4}
	symm := archetype{"symmetric", "对称式 4 片（双胶合对置）", layoutSymmetric, 3, 4}
	retro := archetype{"retrofocus", "反远摄（广角）5 片", layoutRetrofocus, 4, 5}
	macro := archetype{"macro", "微距 4 片（有限共轭优化版）", layoutMacro, 2, 4}

	for _, f := range []float64{35, 50, 75, 100, 135, 150, 200, 300} {
		multiElement(l, cooke, f, 4, true)
	}
	for _, f := range []float64{50, 100, 200} {
		multiElement(l, cooke, f, 5.6, false)
		multiElement(l, cookeBK7, f, 5.6, false)
	}
	for _, f := range []float64{28, 40, 60} {
		multiElement(l, cooke, f, 2.8, false)
	}
	for _, f := range []float64{50, 100, 200} {
		multiElement(l, flattener, f, 4, true)
		multiElement(l, airTri, f, 4, true)
	}
	for _, f := range []float64{35, 50, 85, 100, 135, 200, 300} {
		multiElement(l, tessar, f, 2.8, true)
	}
	for _, f := range []float64{50, 85, 100, 135, 200} {
		multiElement(l, dgauss, f, 1.8, true)
		multiElement(l, dgauss, f, 2.8, false)
	}
	for _, f := range []float64{85, 135, 200} {
		multiElement(l, sonnar, f, 2, false)
	}
	for _, f := range []float64{100, 150, 200, 300, 400, 500} {
		multiElement(l, petzval, f, 5, false)
	}
	for _, f := range []float64{50, 100, 150, 200, 300} {
		multiElement(l, symm, f, 4, false)
	}
	for _, f := range []float64{28, 35} {
		multiElement(l, retro, f, 3.5, false)
	}
	for _, f := range []float64{50, 100} {
		multiElement(l, macro, f, 4, false)
	}
}
