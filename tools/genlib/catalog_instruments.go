// The catalogue: beam expanders, relays, condensers, collimators and the
// rest of the instrument family.
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// buildAfocal assembles a layout whose gap sepCard is solved so that the
// system is afocal: parallel in, parallel out. No focal-length solve is
// applied (an afocal system has none).
func buildAfocal(name string, gs []grp, lead, pupilSemi float64, sepCard int, sepLo, sepHi float64) (*lens.Prescription, error) {
	p := &lens.Prescription{Name: name, ObjectZMM: math.Inf(-1), Stop: 0, Cards: cardsOf(gs, lead)}
	p.Cards[0].SemiAperMM = r3(pupilSemi)
	if err := solveAfocal(p, sepCard, sepLo, sepHi); err != nil {
		return nil, err
	}
	p.Cards[sepCard].ThicknessMM = r3(p.Cards[sepCard].ThicknessMM)
	for i := range p.Cards {
		if !math.IsInf(p.Cards[i].RadiusMM, 0) {
			p.Cards[i].RadiusMM = r4(p.Cards[i].RadiusMM)
		}
	}
	if err := autoApertures(p, pupilSemi); err != nil {
		return nil, err
	}
	if err := p.Validate(); err != nil {
		return nil, err
	}
	if _, err := p.System(); err != nil {
		return nil, err
	}
	return p, nil
}

// expander adds a Galilean (negative eye lens first) or Keplerian (positive
// first) beam expander of expansion ratio mag and input beam radius semiIn.
func expander(l *lib, kepler bool, mag, semiIn, baseF float64, diff int) {
	var gs []grp
	key, zh := "galilean", "伽利略扩束器"
	if kepler {
		key, zh = "kepler", "开普勒扩束器"
		gs = []grp{
			{name: "L1", curv: seedConvPCX("N-BK7", baseF), glasses: []string{"N-BK7", ""}, thick: 0.06 * baseF, after: baseF + mag*baseF},
			{name: "L2", curv: seedConvPCX("N-BK7", mag*baseF), glasses: []string{"N-BK7", ""}, thick: 0.06 * mag * baseF, after: 0.5*mag*baseF + 20},
		}
	} else {
		gs = []grp{
			{name: "L1", curv: seedBic("N-BK7", -baseF), glasses: []string{"N-BK7", ""}, thick: 0.06 * baseF, after: baseF * (mag - 1)},
			{name: "L2", curv: seedConvPCX("N-BK7", mag*baseF), glasses: []string{"N-BK7", ""}, thick: 0.06 * mag * baseF, after: 0.5*mag*baseF + 20},
		}
	}
	p, err := buildAfocal(fmt.Sprintf("%s %.4gx（入射光束半径 %.3gmm）", zh, mag, semiIn), gs, 4, semiIn, 2, baseF*(mag-1)*0.6, baseF*(mag+1)*1.4)
	if err != nil {
		fmt.Printf("  ! expander %s %gx s=%g: %v\n", key, mag, semiIn, err)
		return
	}
	p.Name = fmt.Sprintf("%s %.4gx（入射光束半径 %.3gmm）", zh, mag, semiIn)
	l.add(fmt.Sprintf("expander-%s-m%s-s%s.json", key, numCode(mag), numCode(semiIn)), p, "instrument", 3,
		[]string{"beam-expander", key, "afocal", "mag:" + trimNum(mag, 3), fmt.Sprintf("fov%gdeg", 0.5),
			"visible", "classic", "quality:classic"},
		fmt.Sprintf("%s｜扩展比 %.4gx、入射光束半径 %.3gmm｜经典构造：两片 N-BK7 单透镜按 f1 = %.4gmm、f2 = %.4gmm 布置，间隔由近轴引擎解出使系统无焦；出射角散布与扩展比由实光线追迹实测",
			zh, mag, semiIn, baseF, mag*baseF))
}

// beamReducer is an expander used backwards (mag < 1).
func beamReducer(l *lib, kepler bool, mag, semiIn float64) {
	var gs []grp
	key, zh := "galilean", "伽利略缩束器"
	if kepler {
		key, zh = "kepler", "开普勒缩束器"
		gs = []grp{
			{name: "L1", curv: seedConvPCX("N-BK7", 100), glasses: []string{"N-BK7", ""}, thick: 6, after: 100 + mag*100},
			{name: "L2", curv: seedConvPCX("N-BK7", mag*100), glasses: []string{"N-BK7", ""}, thick: math.Max(3, 6*mag), after: 60},
		}
	} else {
		gs = []grp{
			{name: "L1", curv: seedBic("N-BK7", -100), glasses: []string{"N-BK7", ""}, thick: 6, after: 100 * (mag - 1)},
			{name: "L2", curv: seedConvPCX("N-BK7", mag*100), glasses: []string{"N-BK7", ""}, thick: math.Max(3, 6*mag), after: 60},
		}
	}
	p, err := buildAfocal(fmt.Sprintf("%s %.3gx（入射光束半径 %.4gmm）", zh, mag, semiIn), gs, 4, semiIn, 2, 100*(mag-1)*0.6, 100*(mag+1)*1.4)
	if err != nil {
		fmt.Printf("  ! reducer %s %gx: %v\n", key, mag, err)
		return
	}
	l.add(fmt.Sprintf("reducer-%s-m%s-s%s.json", key, numCode(mag), numCode(semiIn)), p, "instrument", 3,
		[]string{"beam-reducer", key, "afocal", "mag:" + trimNum(mag, 3), "fov0.5deg", "visible",
			"classic", "quality:classic"},
		fmt.Sprintf("%s｜缩小比 %.3gx、入射光束半径 %.4gmm｜经典构造：扩束器的反向使用，间隔由近轴引擎解出使系统无焦",
			zh, mag, semiIn))
}

// relay1to1 is a field-lens relay: an image is carried 1:1 to a new plane,
// with a field lens at the pupil to keep the beam on the following aperture.
func relay1to1(l *lib, f, fieldF float64, withField bool) {
	gs := []grp{
		{name: "L1", curv: seedConvPCX("N-BK7", f), glasses: []string{"N-BK7", ""}, thick: 0.05 * f, after: f},
	}
	if withField {
		gs = append(gs, grp{name: "field lens", curv: seedBic("N-BK7", fieldF), glasses: []string{"N-BK7", ""}, thick: 0.04 * fieldF, after: f})
	} else {
		gs[0].after = 2 * f
	}
	gs = append(gs, grp{name: "L2", curv: seedConvPCX("N-BK7", f), glasses: []string{"N-BK7", ""}, thick: 0.05 * f, after: f})
	semi := math.Min(maxSemi, f/6)
	p, err := buildLayout(fmt.Sprintf("1:1 转像中继 f=%smm", trimNum(f, 0)), gs, 3, f, semi, -2.05*f)
	if err != nil {
		fmt.Printf("  ! relay1to1 f=%g: %v\n", f, err)
		return
	}
	par, _ := p.Paraxial(geo.Wavelengthd)
	file := fmt.Sprintf("relay-1to1-f%s.json", trimNum(f, 0))
	zh := "1:1 转像中继（两片 + 中间场镜）"
	if !withField {
		file = fmt.Sprintf("relay-1to1-plain-f%s.json", trimNum(f, 0))
		zh = "1:1 转像中继（两片，无场镜）"
	}
	l.add(file, p, "instrument", 3,
		[]string{"relay", "1to1", "finite-conjugate", eflTag(f), "visible", "classic", "quality:classic",
			"objfield:" + trimNum(0.02*2.05*f, 2)},
		fmt.Sprintf("%s｜f=%smm、物距 %smm、像距 %smm、近轴倍率 %s｜经典构造：两片单透镜构成 1:1 中继%s，指标在 d 线实测",
			zh, trimNum(f, 0), trimNum(-2.05*f, 1), trimNum(par.ImageZ, 1), trimNum(par.Magnification, 3),
			map[bool]string{true: "，中间场镜把光瞳成像到第二片附近", false: ""}[withField]))
}

// relay4f is a telecentric 4f relay between two focal planes.
func relay4f(l *lib, f1, f2 float64) {
	gs := []grp{
		{name: "L1", curv: seedAchromat("N-BK7", "F2", f1), glasses: []string{"N-BK7", "F2", ""}, thick: 0.05 * f1, after: f1 + f2},
		{name: "L2", curv: seedAchromatRev("F2", "N-BK7", f2), glasses: []string{"F2", "N-BK7", ""}, thick: 0.05 * f2, after: f2},
	}
	semi := math.Min(maxSemi, f1/8)
	p, err := buildLayout(fmt.Sprintf("4f 中继（f1=%smm, f2=%smm）", trimNum(f1, 0), trimNum(f2, 0)), gs, 3, f1, semi, -f1)
	if err != nil {
		fmt.Printf("  ! relay4f %g/%g: %v\n", f1, f2, err)
		return
	}
	par, _ := p.Paraxial(geo.Wavelengthd)
	l.add(fmt.Sprintf("relay-4f-f%s-f%s.json", trimNum(f1, 0), trimNum(f2, 0)), p, "instrument", 4,
		[]string{"relay", "4f", "telecentric", "finite-conjugate", "visible", "classic", "quality:classic",
			"mag:" + trimNum(par.Magnification, 3), "objfield:" + trimNum(0.02*f1, 2)},
		fmt.Sprintf("4f 中继系统｜f1=%smm、f2=%smm、物距 %smm｜经典构造：物在 L1 前焦面、中间像在公共焦面、像在 L2 后焦面（近轴倍率 %s×），两片均为胶合组以压色差",
			trimNum(f1, 0), trimNum(f2, 0), trimNum(-f1, 0), trimNum(par.Magnification, 3)))
}

// condenser is a two-element projector/illuminator condenser working at a
// finite conjugate.
func condenser(l *lib, f float64, ratio float64) {
	gs := []grp{
		{name: "L1", curv: seedConvPCX("N-BK7", 1.6*f), glasses: []string{"N-BK7", ""}, thick: math.Max(3, 0.12*f), after: 0.4 * f},
		{name: "L2", curv: seedConvPCX("N-BK7", f), glasses: []string{"N-BK7", ""}, thick: math.Max(3, 0.12*f), after: 1.2 * f},
	}
	semi := math.Min(maxSemi, f/3)
	objZ := -(1 + 1/ratio) * f * 1.05
	p, err := buildLayout(fmt.Sprintf("聚光镜（f=%smm，%s:1 成像）", trimNum(f, 0), trimNum(ratio, 2)), gs, 3, f, semi, objZ)
	if err != nil {
		fmt.Printf("  ! condenser %g: %v\n", f, err)
		return
	}
	par, _ := p.Paraxial(geo.Wavelengthd)
	l.add(fmt.Sprintf("condenser-f%s-r%s.json", trimNum(f, 0), numCode(ratio)), p, "instrument", 3,
		[]string{"condenser", "illumination", "finite-conjugate", eflTag(f), "visible", "classic",
			"quality:classic", "mag:" + trimNum(par.Magnification, 3), "objfield:" + trimNum(0.02*math.Abs(objZ), 2)},
		fmt.Sprintf("聚光镜｜两片平凸 f=%smm，物距 %smm，%s:1 成像（近轴倍率 %s×）｜经典构造：照明系统常用的两片式聚光镜，第一片承大孔径、第二片收光",
			trimNum(f, 0), trimNum(objZ, 1), trimNum(ratio, 3), trimNum(par.Magnification, 3)))
}

// collimator places a point source exactly on the front focal plane, so the
// exit pencil is parallel: the library's only class whose figure of merit is
// the exit-slope spread.
func collimator(l *lib, crown, flint string, f, fnum float64) {
	seed := seedAchromat(crown, flint, f)
	gs := []grp{{name: "L1", curv: seed, glasses: []string{crown, flint, ""}, thick: math.Max(2, 0.04*f), after: 1.2 * f}}
	semi := math.Min(maxSemi, f/(2*fnum))
	p, err := buildLayout(fmt.Sprintf("准直镜 %s/%s f=%smm", crown, flint, trimNum(f, 0)), gs, 3, f, semi, math.Inf(-1))
	if err != nil {
		fmt.Printf("  ! collimator %g: %v\n", f, err)
		return
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		fmt.Printf("  ! collimator paraxial %g: %v\n", f, err)
		return
	}
	p.ObjectZMM = r3(par.FrontFocalZ)
	if err := autoApertures(p, semi); err != nil {
		fmt.Printf("  ! collimator apertures %g: %v\n", f, err)
		return
	}
	p.Name = fmt.Sprintf("准直镜 %s/%s f=%smm（点光源置于前焦点）", crown, flint, trimNum(f, 0))
	l.add(fmt.Sprintf("collimator-f%s-%s-%s.json", trimNum(f, 0), glassCode(crown), glassCode(flint)), p, "instrument", 3,
		[]string{"collimator", "point-source", "infinite-image", eflTag(f), fnumTag(fnum), "visible",
			"classic", "quality:classic"},
		fmt.Sprintf("准直镜｜%s/%s 双胶合 f=%smm、F/%s｜经典构造：点光源置于前焦点（物距 %smm），出射光束应为平行光；无实像面，指标取实光线追迹的出射角散布",
			crown, flint, trimNum(f, 0), trimNum(fnum, 2), trimNum(p.ObjectZMM, 2)))
}

// telephoto is a long-focus objective with a short physical length: a
// positive front group and a negative rear group.
func telephoto(l *lib, f, telephotoRatio float64) {
	f1 := f * telephotoRatio
	gs := []grp{
		{name: "front group", curv: seedAchromat("N-BK7", "F2", f1), glasses: []string{"N-BK7", "F2", ""}, thick: 0.05 * f1, after: 0.35 * f1},
		{name: "rear group", curv: seedBic("F2", -0.7*f1), glasses: []string{"F2", ""}, thick: 0.03 * f1, after: 0.9 * f},
	}
	semi := math.Min(maxSemi, f/8)
	p, err := buildLayout(fmt.Sprintf("远摄物镜 f=%smm（物理长度比 %s）", trimNum(f, 0), trimNum(telephotoRatio, 2)), gs, 3, f, semi, math.Inf(-1))
	if err != nil {
		fmt.Printf("  ! telephoto %g: %v\n", f, err)
		return
	}
	par, _ := p.Paraxial(geo.Wavelengthd)
	l.add(fmt.Sprintf("telephoto-f%s-ratio%s.json", trimNum(f, 0), numCode(telephotoRatio)), p, "instrument", 4,
		[]string{"telephoto", "long-focus", eflTag(f), "visible", "classic", "quality:classic", "front-stop"},
		fmt.Sprintf("远摄物镜（Telephoto）｜f=%smm、F/%s、前组 f=%smm + 后组负透镜｜经典构造：后组为负透镜把主面后移，总长短于焦距（后焦距 %smm ＞ 总长 %smm）",
			trimNum(f, 0), trimNum(fnumFor(f, 8), 2), trimNum(f1, 0), trimNum(par.BackFocalZ-par.VertexZ[par.StopIndex], 1), trimNum(par.VertexZ[len(par.VertexZ)-1], 1)))
}

// scanLens is an f-theta type scan lens: a positive group with a negative
// corrector, used to focus a collimated beam onto a flat field.
func scanLens(l *lib, f float64) {
	gs := []grp{
		{name: "L1", curv: seedConvPCX("N-BK7", 1.5*f), glasses: []string{"N-BK7", ""}, thick: math.Max(3, 0.05*f), after: 0.3 * f},
		{name: "L2", curv: seedBic("F2", -1.2*f), glasses: []string{"F2", ""}, thick: math.Max(2, 0.03*f), after: 0.25 * f},
		{name: "L3", curv: seedBic("N-BK7", 0.9*f), glasses: []string{"N-BK7", ""}, thick: math.Max(3, 0.05*f), after: 0.8 * f},
	}
	semi := math.Min(maxSemi, f/12)
	p, err := buildLayout(fmt.Sprintf("f-θ 扫描镜头 f=%smm", trimNum(f, 0)), gs, 3, f, semi, math.Inf(-1))
	if err != nil {
		fmt.Printf("  ! scan lens %g: %v\n", f, err)
		return
	}
	l.add(fmt.Sprintf("scanlens-f%s.json", trimNum(f, 0)), p, "instrument", 4,
		[]string{"scan-lens", "f-theta", "flat-field", eflTag(f), "visible", "classic", "quality:classic",
			"fov2deg", "front-stop"},
		fmt.Sprintf("f-θ 扫描镜头（激光打标/扫描用）｜f=%smm、F/%s、三片式（正-负-正）｜经典构造：负片压平像场，使像高近似 ∝ 扫描角（fθ 条件）；像质由实测记录",
			trimNum(f, 0), trimNum(fnumFor(f, 12), 2)))
}

// erector is the 1:1 image-erecting relay of a prism/relay-tube system: two
// achromatic groups with a field stop between them.
func erector(l *lib, f float64) {
	gs := []grp{
		{name: "L1/L2", curv: seedAchromat("N-BK7", "F2", f), glasses: []string{"N-BK7", "F2", ""}, thick: 0.06 * f, after: f},
		{name: "L3/L4", curv: seedAchromatRev("F2", "N-BK7", f), glasses: []string{"F2", "N-BK7", ""}, thick: 0.06 * f, after: f},
	}
	semi := math.Min(maxSemi, f/8)
	p, err := buildLayout(fmt.Sprintf("1:1 正像中继（转像系统）f=%smm", trimNum(f, 0)), gs, 3, f, semi, -1.9*f)
	if err != nil {
		fmt.Printf("  ! erector %g: %v\n", f, err)
		return
	}
	par, _ := p.Paraxial(geo.Wavelengthd)
	l.add(fmt.Sprintf("erector-1to1-f%s.json", trimNum(f, 0)), p, "instrument", 4,
		[]string{"erector", "relay", "1to1", "finite-conjugate", eflTag(f), "visible", "classic",
			"quality:classic", "mag:" + trimNum(par.Magnification, 3), "objfield:" + trimNum(0.02*1.9*f, 2)},
		fmt.Sprintf("转像系统（1:1 正像中继）｜两组胶合 f=%smm，物距 %smm、像距 %smm、近轴倍率 %s｜经典构造：望远镜/潜望镜里的转像镜组，中间像位于两组之间的公共焦面",
			trimNum(f, 0), trimNum(-1.9*f, 1), trimNum(par.ImageZ, 1), trimNum(par.Magnification, 3)))
}

func buildInstruments(l *lib) {
	for _, m := range []float64{2, 3, 4, 5, 6, 8, 10, 15, 20} {
		expander(l, false, m, 2, 25, 3)
		expander(l, true, m, 2, 25, 3)
	}
	for _, m := range []float64{3, 5, 8} {
		expander(l, false, m, 5, 40, 3)
		expander(l, true, m, 5, 40, 3)
	}
	for _, m := range []float64{4, 10} {
		expander(l, true, m, 1, 15, 3)
		expander(l, false, m, 1, 15, 3)
	}
	for _, m := range []float64{2, 4, 8} {
		beamReducer(l, false, 1/m, 4)
	}
	beamReducer(l, true, 0.5, 6)
	beamReducer(l, true, 0.25, 8)
	for _, f := range []float64{50, 100, 200} {
		relay1to1(l, f, f/2, true)
		relay1to1(l, f, 0, false)
	}
	relay4f(l, 100, 100)
	relay4f(l, 150, 75)
	for _, c := range []struct {
		f, ratio float64
	}{{25, 1}, {40, 2}, {60, 3}, {100, 1}} {
		condenser(l, c.f, c.ratio)
	}
	for _, f := range []float64{50, 100, 150, 200, 300, 500} {
		collimator(l, "N-BK7", "F2", f, 5)
	}
	for _, c := range []struct {
		f, ratio float64
	}{{150, 0.85}, {300, 0.75}, {500, 0.7}} {
		telephoto(l, c.f, c.ratio)
	}
	for _, f := range []float64{100, 250} {
		scanLens(l, f)
	}
	for _, f := range []float64{50, 100, 200} {
		erector(l, f)
	}
}
