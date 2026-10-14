// The catalogue: eyepieces, microscope objectives and telescopes.
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// ---------- eyepiece layouts ----------
//
// Eyepieces are traced the physically meaningful way round: the surfaces are
// listed eye-lens-first, so an object at infinity stands for the collimated
// pencil that leaves the eye, exactly the conjugate the eyepiece is designed
// for. The entrance pupil of the trace is the eye pupil (the aperture stop
// of the reversed system). Their EFL comes out ≈ +f_e, the eyepiece focal
// length quoted on the box.

func layoutHuygens(fe float64) []grp {
	const g = "N-BK7"
	return []grp{
		{name: "eye lens", curv: seedConvPCX(g, fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: fe},
		{name: "field lens", curv: seedConvPCX(g, fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: 0.3*fe + 20},
	}
}

func layoutRamsden(fe float64) []grp {
	const g = "N-BK7"
	return []grp{
		{name: "eye lens", curv: seedFlatPCX(g, fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: 0.70 * fe},
		{name: "field lens", curv: seedConvPCX(g, fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: 0.3*fe + 20},
	}
}

func layoutKellner(fe float64) []grp {
	const g = "N-BK7"
	return []grp{
		{name: "eye doublet", curv: seedAchromat(g, "F2", 0.75*fe),
			glasses: []string{g, "F2", ""}, thick: 0.14 * fe, after: 0.70 * fe},
		{name: "field lens", curv: seedConvPCX(g, fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: 0.3*fe + 20},
	}
}

func layoutPlossl(fe float64) []grp {
	const g = "N-BK7"
	return []grp{
		{name: "eye doublet", curv: seedAchromat(g, "F2", 0.62*fe),
			glasses: []string{g, "F2", ""}, thick: 0.16 * fe, after: 0.55 * fe},
		{name: "field doublet", curv: seedAchromatRev("F2", g, 0.62*fe),
			glasses: []string{"F2", g, ""}, thick: 0.16 * fe, after: 0.3*fe + 20},
	}
}

func layoutAbbe(fe float64) []grp {
	const g = "N-SK16"
	return []grp{
		{name: "eye lens", curv: seedConvPCX(g, 0.55*fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: 0.35 * fe},
		{name: "field triplet", curv: seedAchromat(g, "F2", 0.9*fe),
			glasses: []string{g, "F2", ""}, thick: 0.18 * fe, after: 0.2*fe + 20},
	}
}

func layoutErfle(fe float64) []grp {
	const g = "N-BK7"
	return []grp{
		{name: "eye doublet", curv: seedAchromat(g, "F2", 0.9*fe),
			glasses: []string{g, "F2", ""}, thick: 0.16 * fe, after: 0.30 * fe},
		{name: "center lens", curv: seedBic(g, 0.75*fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: 0.30 * fe},
		{name: "field doublet", curv: seedAchromatRev("F2", g, 0.9*fe),
			glasses: []string{"F2", g, ""}, thick: 0.16 * fe, after: 0.3*fe + 20},
	}
}

func layoutWideField(fe float64) []grp {
	const g = "N-SK16"
	return []grp{
		{name: "eye doublet", curv: seedAchromat(g, "F2", 1.1*fe),
			glasses: []string{g, "F2", ""}, thick: 0.18 * fe, after: 0.22 * fe},
		{name: "L3", curv: seedBic(g, 1.2*fe), glasses: []string{g, ""}, thick: 0.10 * fe, after: 0.22 * fe},
		{name: "L4", curv: seedBic("F2", -1.0*fe), glasses: []string{"F2", ""}, thick: 0.09 * fe, after: 0.22 * fe},
		{name: "field doublet", curv: seedAchromatRev("F2", g, 1.1*fe),
			glasses: []string{"F2", g, ""}, thick: 0.18 * fe, after: 0.2*fe + 20},
	}
}

// eyepiece builds one eyepiece example: the eye pupil is the aperture stop.
func eyepiece(l *lib, key, zh string, layout func(float64) []grp, fe, eyePupilMM, fieldDeg float64, diff int, polishIt bool) {
	p, err := buildLayout(fmt.Sprintf("%s f_e=%smm", zh, trimNum(fe, 0)), layout(fe), 12, fe, eyePupilMM/2, math.Inf(-1))
	if err != nil {
		fmt.Printf("  ! eyepiece %s fe=%g: %v\n", key, fe, err)
		return
	}
	quality := "classic"
	stats := ""
	if polishIt {
		q, qual, res := polishIfBetter(p, fe, fieldDeg, 120)
		p = q
		quality = qual
		if qual == "optimized" && res != nil {
			stats = fmt.Sprintf("｜DLS 抛光：%d 次接受迭代／%d 次求值，merit=%.4g", res.Iter, res.Evals, res.Merit)
		}
	}
	p.Name = fmt.Sprintf("%s 目镜 f_e=%smm", zh, trimNum(fe, 0))
	file := fmt.Sprintf("eyepiece-%s-fe%s.json", key, trimNum(fe, 2))
	if quality == "optimized" {
		file = fmt.Sprintf("eyepiece-%s-fe%s-dls.json", key, trimNum(fe, 2))
		p.Name = fmt.Sprintf("%s 目镜 f_e=%smm（DLS 抛光）", zh, trimNum(fe, 2))
	}
	l.add(file, p, "eyepiece", diff,
		[]string{"eyepiece", key, "reverse-traced", "f_e:" + trimNum(fe, 2),
			fmt.Sprintf("fov%gdeg", fieldDeg), "visible", "quality:" + quality, "eye-pupil-stop"},
		fmt.Sprintf("%s｜f_e=%smm（眼端反向追迹：物体在无限远、孔径光阑＝眼瞳直径 %smm）｜经典构造%s；指标在 d 线实测",
			zh, trimNum(fe, 2), trimNum(eyePupilMM, 1), stats))
}

func buildEyepieces(l *lib) {
	for _, fe := range []float64{20, 25, 32, 40, 50} {
		eyepiece(l, "huygens", "惠更斯（Huygens）两片式", layoutHuygens, fe, 3, 5, 2, false)
	}
	for _, fe := range []float64{15, 20, 25, 32, 40} {
		eyepiece(l, "ramsden", "冉斯登（Ramsden）两片式", layoutRamsden, fe, 3, 5, 2, false)
	}
	for _, fe := range []float64{10, 15, 20, 25, 32} {
		eyepiece(l, "kellner", "凯尔纳（Kellner）三片式", layoutKellner, fe, 3, 5, 3, true)
	}
	for _, fe := range []float64{10, 12.5, 15, 20, 25, 32, 40} {
		eyepiece(l, "plossl", "普洛塞尔（Plössl）四片式", layoutPlossl, fe, 4, 5, 3, true)
	}
	for _, fe := range []float64{6, 10, 12.5, 15, 20, 25} {
		eyepiece(l, "abbe", "阿贝无畸变（Abbe orthoscopic）四片式", layoutAbbe, fe, 3, 7, 4, false)
	}
	for _, fe := range []float64{20, 25, 32} {
		eyepiece(l, "erfle", "埃尔弗（Erfle）五片广角", layoutErfle, fe, 5, 8, 4, true)
	}
	for _, fe := range []float64{15, 25} {
		eyepiece(l, "wide6", "六片超广角（Wide-field 型）", layoutWideField, fe, 6, 10, 5, false)
	}
	// different eye pupils and glass pairs at the two most common focal
	// lengths: eye relief and pupil diameter are the parameters that matter
	// in practice.
	for _, fe := range []float64{25, 40} {
		eyepiece(l, "plossl-pup6", "普洛塞尔（6mm 眼瞳）", layoutPlossl, fe, 6, 5, 3, false)
	}
	for _, fe := range []float64{20, 32} {
		eyepiece(l, "huygens-pup5", "惠更斯（5mm 眼瞳）", layoutHuygens, fe, 5, 5, 2, false)
	}
	for _, fe := range []float64{18, 24} {
		eyepiece(l, "kellner-wide", "凯尔纳（大视场版）", layoutKellner, fe, 3.5, 8, 3, false)
	}
	for _, fe := range []float64{9, 16} {
		eyepiece(l, "abbe-short", "阿贝无畸变（短焦版）", layoutAbbe, fe, 2.5, 7, 4, false)
	}
}

// ---------- microscope objectives ----------

// layoutLister is the classic two-group (Lister) objective: two cemented
// doublets well separated, the front one strongly bent.
func layoutLister(front, rear float64) []grp {
	return []grp{
		{name: "front doublet", curv: seedAchromat("N-SK16", "N-SF5", front),
			glasses: []string{"N-SK16", "N-SF5", ""}, thick: 0.06 * front, after: 0.55 * front},
		{name: "rear doublet", curv: seedAchromat("N-BK7", "F2", rear),
			glasses: []string{"N-BK7", "F2", ""}, thick: 0.07 * rear, after: 0.6 * rear},
	}
}

// microObjective builds a finite-conjugate objective: the object distance is
// solved with the paraxial engine so that the transverse magnification is
// the quoted one, and the image distance is reported as the tube length.
func microObjective(l *lib, key, zh string, layout func(float64) []grp, efl, mag, objNA, fieldMM float64, diff int, polishIt bool) {
	f := efl
	semi := math.Abs(f) * objNA / math.Sqrt(1-objNA*objNA)
	if semi > maxSemi {
		semi = maxSemi
	}
	objZ0 := -math.Abs(f) * (1 + 1/mag) * 1.3
	p, err := buildLayout(fmt.Sprintf("%s %gx", zh, mag), layout(f), 2.5, f, semi, objZ0)
	if err != nil {
		fmt.Printf("  ! micro %s %gx: %v\n", key, mag, err)
		return
	}
	// solve the object distance for the quoted magnification; the object
	// must sit outside the front focal point for a real image, so the
	// bracket spans from just inside the focal point out to a few focal
	// lengths.
	target := -math.Abs(mag)
	fmag := func(z float64) float64 {
		p.ObjectZMM = z
		par, err := p.Paraxial(geo.Wavelengthd)
		if err != nil {
			return math.NaN()
		}
		return par.Magnification - target
	}
	ffl := math.NaN()
	if par0, err := p.Paraxial(geo.Wavelengthd); err == nil {
		ffl = par0.FrontFocalZ
	}
	if math.IsNaN(ffl) || math.IsInf(ffl, 0) {
		ffl = -math.Abs(f)
	}
	// Newton's form: the object sits f/mag in front of the front focal point.
	u := math.Abs(f) / math.Abs(mag)
	z, err := solve1D(fmag, ffl-5*u, ffl-u/5, 400, 100)
	if err != nil {
		fmt.Printf("  ! micro %s %gx objective-distance solve: %v\n", key, mag, err)
		p.ObjectZMM = r3(ffl - u)
	} else {
		p.ObjectZMM = r3(z)
	}
	if err := autoApertures(p, semi); err != nil {
		fmt.Printf("  ! micro apertures %s %gx: %v\n", key, mag, err)
	}
	quality := "classic"
	stats := ""
	if polishIt {
		q, qual, res := polishIfBetter(p, f, 0, 150)
		q.ObjectZMM = p.ObjectZMM
		p = q
		quality = qual
		if qual == "optimized" && res != nil {
			stats = fmt.Sprintf("｜DLS 抛光：%d 次接受迭代／%d 次求值，merit=%.4g", res.Iter, res.Evals, res.Merit)
		}
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		fmt.Printf("  ! micro %s %gx paraxial: %v\n", key, mag, err)
		return
	}
	p.Name = fmt.Sprintf("%s 显微物镜 %gx／%g 筒长（物距 %smm）", zh, mag, 150.0, trimNum(p.ObjectZMM, 2))
	file := fmt.Sprintf("micro-%s-%gx.json", key, mag)
	if quality == "optimized" {
		file = fmt.Sprintf("micro-%s-%gx-dls.json", key, mag)
	}
	l.add(file, p, "micro-telescope", diff,
		[]string{"microscope", key, "finite-conjugate", fmt.Sprintf("mag:%s", trimNum(par.Magnification, 4)),
			"NA:" + trimNum(objNA, 2), "visible", "quality:" + quality,
			"objfield:" + trimNum(fieldMM, 3)},
		fmt.Sprintf("%s｜%gx、物方 NA≈%s、筒长 %smm（像距 %smm）｜有限共轭：物距由近轴引擎解出以满足倍率%s",
			zh, mag, trimNum(objNA, 2), trimNum(150.0, 0), trimNum(par.ImageZ-par.VertexZ[par.StopIndex], 2), stats))
}

func buildMicro(l *lib) {
	type mcase struct {
		key, zh string
		efl     float64
		mag     float64
		na      float64
	}
	for _, c := range []mcase{
		{"achr", "单组胶合有限共轭物镜", 40, 4, 0.10},
		{"achr", "单组胶合有限共轭物镜", 16, 10, 0.15},
		{"achr", "单组胶合有限共轭物镜", 8, 20, 0.20},
		{"achr", "单组胶合有限共轭物镜", 4, 40, 0.25},
		{"achr", "单组胶合有限共轭物镜", 2.8, 60, 0.30},
	} {
		field := 0.35 * (10 / c.mag)
		microObjective(l, c.key, c.zh, func(f float64) []grp {
			return []grp{{name: "doublet", curv: seedAchromat("N-BK7", "F2", f), glasses: []string{"N-BK7", "F2", ""}, thick: 0.08 * f, after: 1.2 * f}}
		}, c.efl, c.mag, c.na, field, 3, c.mag <= 20)
	}
	// Lister-type two-group objectives
	for _, c := range []struct {
		mag, na, front, rear float64
	}{{10, 0.18, 12, 20}, {20, 0.25, 8, 14}, {40, 0.30, 4, 8}} {
		field := 0.35 * (10 / c.mag)
		microObjective(l, "lister", "Lister 式双组物镜", func(f float64) []grp {
			return layoutLister(f*0.75, f)
		}, c.front, c.mag, c.na, field, 4, true)
	}
	// long working distance variants (weaker front group, larger object gap)
	for _, c := range []struct {
		mag, na, efl float64
	}{{10, 0.12, 16}, {20, 0.15, 8}} {
		microObjective(l, "lwd", "长工作距离物镜（LWD）", func(f float64) []grp {
			return []grp{
				{name: "front doublet", curv: seedAchromat("N-BK7", "F2", f), glasses: []string{"N-BK7", "F2", ""}, thick: 0.08 * f, after: 1.6 * f},
			}
		}, c.efl, c.mag, c.na, 0.35, 4, false)
	}
}

// ---------- telescopes ----------

// telescope builds an afocal telescope: the separation is solved on the
// paraxial engine so that the exit pencil is collimated (EFL → ∞).
func telescope(l *lib, key, zh string, gs []grp, sepCard int, fobj, feye, fieldDeg float64, diff int) {
	p := &lens.Prescription{Name: fmt.Sprintf("%s f_obj=%smm f_eye=%smm", zh, trimNum(fobj, 0), trimNum(feye, 0)),
		ObjectZMM: math.Inf(-1), Stop: 0, Cards: cardsOf(gs, 5)}
	semi := math.Min(maxSemi, fobj/12)
	p.Cards[0].SemiAperMM = r3(semi)
	if err := solveAfocal(p, sepCard, (fobj+feye)*0.5, (fobj+feye)*1.5); err != nil {
		fmt.Printf("  ! telescope %s: %v\n", key, err)
		return
	}
	p.Cards[sepCard].ThicknessMM = r3(p.Cards[sepCard].ThicknessMM)
	if err := autoApertures(p, semi); err != nil {
		fmt.Printf("  ! telescope apertures %s: %v\n", key, err)
		return
	}
	mag := fobj / math.Abs(feye)
	p.Name = fmt.Sprintf("%s 倍率 %.3gx", zh, mag)
	l.add(fmt.Sprintf("telescope-%s-%s-%s.json", key, trimNum(fobj, 0), trimNum(math.Abs(feye), 0)), p, "micro-telescope", diff,
		[]string{"telescope", key, "afocal", "mag:" + trimNum(mag, 3), "visible", "classic",
			"quality:classic", fmt.Sprintf("fov%gdeg", fieldDeg), "front-stop"},
		fmt.Sprintf("%s｜物镜 f=%smm、目镜 f=%smm、倍率≈%.3gx｜经典构造：分离量由近轴引擎解出使系统无焦（出射平行光），倍率与出射角散布由实光线追迹实测",
			zh, trimNum(fobj, 0), trimNum(feye, 0), mag))
}

func buildTelescopes(l *lib) {
	keplerObjective := func(fobj float64) grp {
		return grp{name: "objective", curv: seedAchromat("N-BK7", "F2", fobj),
			glasses: []string{"N-BK7", "F2", ""}, thick: 0.025 * fobj, after: fobj}
	}
	keplerEye := func(feye float64) grp {
		return grp{name: "eye lens", curv: seedConvPCX("N-BK7", feye), glasses: []string{"N-BK7", ""}, thick: 0.1 * feye, after: 0.5*feye + 20}
	}
	galileanEye := func(feye float64) grp {
		return grp{name: "eye lens", curv: seedBic("N-BK7", feye), glasses: []string{"N-BK7", ""}, thick: 0.1 * math.Abs(feye), after: 0.5*math.Abs(feye) + 20}
	}
	for _, c := range []struct{ fobj, feye float64 }{
		{400, 25}, {600, 25}, {800, 20}, {1000, 25}, {1200, 20}, {300, 25}, {500, 15}, {750, 10}, {1500, 25},
	} {
		telescope(l, "kepler", "开普勒式（Keplerian）望远镜", []grp{keplerObjective(c.fobj), keplerEye(c.feye)}, 3, c.fobj, c.feye, 0.5, 4)
	}
	for _, c := range []struct{ fobj, feye float64 }{
		{400, -25}, {600, -20}, {800, -25}, {1000, -20}, {300, -15}, {500, -20}, {800, -10}, {1200, -25}, {240, -12},
	} {
		telescope(l, "galilei", "伽利略式（Galilean）望远镜（目镜为负透镜）", []grp{keplerObjective(c.fobj), galileanEye(c.feye)}, 3, c.fobj, c.feye, 0.5, 4)
	}
	// Keplerian finder/spotting scopes with an achromatic eye lens
	for _, c := range []struct{ fobj, feye float64 }{
		{180, 20}, {240, 24}, {300, 20}, {360, 18},
	} {
		telescope(l, "finder", "寻星镜／瞄准镜（开普勒式，目镜为胶合组）",
			[]grp{
				keplerObjective(c.fobj),
				{name: "eye doublet", curv: seedAchromat("N-BK7", "F2", c.feye), glasses: []string{"N-BK7", "F2", ""}, thick: 0.12 * c.feye, after: 0.5*c.feye + 20},
			}, 3, c.fobj, c.feye, 0.5, 4)
	}
	// compound microscope: objective + tube lens (finite -> infinity corrected)
	for _, c := range []struct {
		fobj, ftube, mag float64
	}{{18, 180, 10}, {9, 180, 20}, {4.5, 180, 40}} {
		f := c.fobj
		gs := []grp{
			{name: "objective", curv: seedAchromat("N-BK7", "F2", f), glasses: []string{"N-BK7", "F2", ""}, thick: 0.08 * f, after: 1.4 * f},
			{name: "tube lens", curv: seedAchromat("N-SK16", "N-SF5", c.ftube), glasses: []string{"N-SK16", "N-SF5", ""}, thick: 0.03 * c.ftube, after: 0.6 * c.ftube},
		}
		p, err := buildLayout(fmt.Sprintf("无限远矫正显微物镜组 %gx（物镜 %smm + 镜筒透镜 %smm）", c.mag, trimNum(f, 0), trimNum(c.ftube, 0)),
			gs, 2.5, c.ftube, 4, -1.35*f)
		if err != nil {
			fmt.Printf("  ! compound micro %gx: %v\n", c.mag, err)
			continue
		}
		// distance from the objective to the tube lens = the parfocal length
		p.Cards[4].ThicknessMM = r3(0.55 * c.ftube)
		if err := fillInfRun(p, 4); err != nil {
			fmt.Printf("  ! compound micro fill %gx: %v\n", c.mag, err)
			continue
		}
		par, perr := p.Paraxial(geo.Wavelengthd)
		if perr != nil {
			fmt.Printf("  ! compound micro paraxial %gx: %v\n", c.mag, perr)
			continue
		}
		l.add(fmt.Sprintf("micro-compound-%gx.json", c.mag), p, "micro-telescope", 5,
			[]string{"microscope", "infinity-corrected", "compound", "mag:" + trimNum(par.Magnification, 3),
				"visible", "classic", "quality:classic", "objfield:0.2"},
			fmt.Sprintf("无限远矫正显微镜组｜%gx（物镜 f=%smm + 镜筒透镜 f=%smm）｜经典构造：物镜把物点送到无限远、镜筒透镜成像于 %smm（倍率 %s×，近轴实测）",
				c.mag, trimNum(f, 1), trimNum(c.ftube, 0), trimNum(par.ImageZ, 1), trimNum(par.Magnification, 2)))
	}
}

// fillInfRun re-runs the aperture fill after a geometry edit that does not
// change the focal length.
func fillInfRun(p *lens.Prescription, _ int) error {
	if err := p.Validate(); err != nil {
		return err
	}
	_, err := p.System()
	if err != nil {
		return err
	}
	return autoApertures(p, p.Cards[0].SemiAperMM)
}
