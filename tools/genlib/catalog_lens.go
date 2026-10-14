// The catalogue: single-lens and cemented-doublet families.
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// maxSemi caps the clear semi-aperture of the library examples: a very
// large aperture buys nothing for a validation library and makes the models
// unwieldy, so the F-number of a long lens is relaxed instead.
const maxSemi = 60.0

// fnumFor returns the F-number actually used for a requested one, relaxed so
// that the clear semi-aperture never exceeds maxSemi.
func fnumFor(efl, want float64) float64 {
	min := math.Abs(efl) / (2 * maxSemi)
	if want < min {
		return math.Round(min*100) / 100
	}
	return want
}

func thickFor(efl float64) float64 {
	t := 0.05 * math.Abs(efl)
	switch {
	case t < 1.2:
		return 1.2
	case t > 20:
		return 20
	}
	return t
}

// ---------- singlets ----------

type shapeT struct {
	key, zh string
	seed    func(string, float64) []float64
}

func singlet(l *lib, sh shapeT, glass string, fl, wantFnum float64, fov float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	f := fl
	g := grp{
		name:    "L1",
		curv:    sh.seed(glass, f),
		glasses: []string{glass, ""},
		thick:   thickFor(efl),
		after:   0.2*efl + 10,
	}
	file := fmt.Sprintf("singlet-%s-f%s-%s.json", sh.key, trimNum(efl, 2), glassCode(glass))
	name := fmt.Sprintf("%s单透镜 %s f=%.4gmm F/%s", sh.zh, glass, efl, trimNum(fnum, 2))
	p, err := buildLayout(name, []grp{g}, 2.5, f, semi, math.Inf(-1))
	if err != nil {
		fmt.Printf("  ! singlet %s %s f=%g: %v\n", sh.key, glass, f, err)
		return
	}
	note := fmt.Sprintf("单透镜（%s）｜%s｜设计焦距 %smm，口径 %smm，F/%s",
		sh.zh, glass, trimNum(efl, 3), trimNum(semi, 3), trimNum(fnum, 2))
	if sh.key == "pcx" {
		note += "｜凸面朝物为平行光最佳形式（球差最小）"
	}
	if sh.key == "pcxflat" {
		note += "｜平面朝物：球差约为最佳形式的 4 倍，作为对照"
	}
	if fl < 0 {
		note += "｜负透镜：平行光入射为虚像，用实光线追迹外推虚焦面"
	}
	l.add(file, p, "singlet", 1, []string{
		"singlet", sh.key, glass, fnumTag(fnum), eflTag(efl),
		"visible", "classic", "quality:analytic", fmt.Sprintf("fov%gdeg", fov), "sweep:f",
	}, note)
}

func buildSinglets(l *lib) {
	shapes := map[string]shapeT{
		"pcx":     {"pcx", "平凸（凸面朝物）", seedConvPCX},
		"pcxflat": {"pcxflat", "平凸（平面朝物）", seedFlatPCX},
		"dcx":     {"dcx", "双凸", seedBic},
	}
	pos := []float64{25, 30, 35, 40, 50, 60, 75, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000}
	for _, glass := range []string{"N-BK7", "FusedSilica", "N-SK16"} {
		for _, f := range pos {
			singlet(l, shapes["pcx"], glass, f, 5, 0.5)
		}
	}
	for _, glass := range []string{"N-BK7", "FusedSilica"} {
		for _, f := range []float64{25, 35, 50, 75, 100, 150, 200, 300, 400, 500, 750, 1000} {
			singlet(l, shapes["dcx"], glass, f, 5, 0.5)
		}
	}
	for _, f := range []float64{50, 100, 200, 400} {
		singlet(l, shapes["pcxflat"], "N-BK7", f, 5, 0.5)
	}
	// plano-concave and negative meniscus: virtual image, traced by the
	// real-ray focus scan.
	neg := shapeT{"pcm", "平凹（平面朝物）", seedFlatPCX}
	negC := shapeT{"pcmconc", "平凹（凹面朝物）", seedConvPCX}
	menN := shapeT{"menn", "负弯月", seedMenN}
	for _, glass := range []string{"N-BK7", "N-SF5"} {
		for _, f := range []float64{-25, -50, -100, -150, -200, -300, -500, -1000} {
			singlet(l, neg, glass, f, 5, 0.5)
		}
	}
	for _, f := range []float64{-50, -100, -200, -400} {
		singlet(l, negC, "N-BK7", f, 5, 0.5)
	}
	for _, glass := range []string{"N-BK7", "F2"} {
		for _, f := range []float64{-50, -100, -200, -400} {
			singlet(l, menN, glass, f, 5, 0.5)
		}
	}
	// positive meniscus
	menP := shapeT{"menp", "正弯月", seedMenP}
	for _, glass := range []string{"N-BK7", "N-SK16"} {
		for _, f := range []float64{50, 75, 100, 150, 200, 300, 400, 500} {
			singlet(l, menP, glass, f, 5, 0.5)
		}
	}
	// F-number sweep at f = 100 (N-BK7, best form)
	for _, fn := range []float64{1.4, 2, 2.4, 2.8, 3.5, 4, 5.6, 8, 11, 16} {
		efl, semi := 100.0, 100/(2*fn)
		g := grp{name: "L1", curv: seedConvPCX("N-BK7", efl), glasses: []string{"N-BK7", ""}, thick: thickFor(efl), after: 30}
		p, err := buildLayout(fmt.Sprintf("平凸单透镜 N-BK7 f=100mm F/%s（口径扫描）", trimNum(fn, 2)),
			[]grp{g}, 2.5, efl, semi, math.Inf(-1))
		if err != nil {
			fmt.Printf("  ! singlet fnum %g: %v\n", fn, err)
			continue
		}
		l.add(fmt.Sprintf("singlet-pcx-f100-nbk7-%s.json", fnumCode(fn)), p, "singlet", 1,
			[]string{"singlet", "pcx", "N-BK7", fnumTag(fn), "efl:100", "visible", "classic",
				"quality:analytic", "fov0.5deg", "sweep:Fnum"},
			fmt.Sprintf("平凸单透镜｜N-BK7｜f=100mm 固定，口径随 F/%s 变化：观察球差随口径的三次方增长", trimNum(fn, 2)))
	}
	// bending family: curvature ratio c2/c1 swept at fixed f and F-number
	for _, k := range []float64{-2, -1.5, -1, -0.5, 0, 0.5, 0.9, 1.3, 2} {
		target := 100.0
		sign := 1.0
		if k > 1 {
			// a meniscus with |c2/c1| > 1 has negative power: build it from
			// negative curvatures so the seed already has the right sign
			target, sign = -100, -1
		}
		g := grp{name: "L1", curv: []float64{sign, sign * k}, glasses: []string{"N-BK7", ""}, thick: thickFor(100), after: 30}
		p, err := buildLayout(fmt.Sprintf("弯曲族 N-BK7 f=%smm（c2/c1=%.2f）", trimNum(target, 0), k),
			[]grp{g}, 2.5, target, 10, math.Inf(-1))
		if err != nil {
			fmt.Printf("  ! bend k=%g: %v\n", k, err)
			continue
		}
		l.add(fmt.Sprintf("singlet-bend-k%s-f%s-nbk7.json", numCode(k), trimNum(math.Abs(target), 0)), p, "singlet", 2,
			[]string{"singlet", "bending", "N-BK7", "F/5", "visible", "classic", "quality:analytic",
				"fov0.5deg", "sweep:bending", "efl:" + trimNum(target, 0)},
			fmt.Sprintf("单透镜弯曲族｜N-BK7｜f=%smm、F/5 固定，第二个面的曲率比为 %s：弯月形（比值→1）球差显著下降", trimNum(target, 0), trimNum(k, 2)))
	}
	// aspheric singlets: the front surface carries A4
	for _, f := range []float64{25, 50, 100, 200} {
		semi := f / (2 * fnumFor(f, 3))
		g := grp{name: "L1", curv: seedConvPCX("N-BK7", f), glasses: []string{"N-BK7", ""}, thick: thickFor(f), after: 30}
		p, err := buildPrescriptionAsphere(fmt.Sprintf("非球面单透镜 N-BK7 f=%smm F/3", trimNum(f, 0)), []grp{g}, 2.5, f, semi)
		if err != nil {
			fmt.Printf("  ! asphere f=%g: %v\n", f, err)
			continue
		}
		l.add(fmt.Sprintf("singlet-asphere-f%s-nbk7.json", trimNum(f, 0)), p, "singlet", 2,
			[]string{"singlet", "asphere", "N-BK7", "F/3", "visible", "classic", "quality:analytic",
				"fov0.5deg", "efl:" + trimNum(f, 0), "sweep:shape"},
			fmt.Sprintf("非球面单透镜｜N-BK7｜f=%smm、F/%s：球面项 A4 由实光线追迹求解，使边缘球差归零（施密特式修正）", trimNum(f, 0), trimNum(fnumFor(f, 3), 2)))
	}
}

// buildPrescriptionAsphere builds a biconvex-free aspheric singlet whose
// second surface carries an A4 term solved on the real trace so that the
// marginal focus coincides with the paraxial one.
func buildPrescriptionAsphere(name string, gs []grp, lead, efl, semi float64) (*lens.Prescription, error) {
	p, err := buildLayout(name, gs, lead, efl, semi, math.Inf(-1))
	if err != nil {
		return nil, err
	}
	// Put the asphere on the second surface (the convex one), then solve A4.
	last := len(p.Cards) - 1
	cards := append([]lens.Card(nil), p.Cards...)
	cards[last].Coeffs = []float64{0}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		return nil, err
	}
	zFocus := par.BackFocalZ
	h := semi * 0.98
	z0 := -(spanOf2(cards) * 2) - 100
	resid := func(a4 float64) float64 {
		c := append([]lens.Card(nil), cards...)
		c[last].Coeffs = []float64{a4}
		sys, err := buildSysFromCards(c)
		if err != nil {
			return math.NaN()
		}
		r := geo.NewRay(geo.V(0, h, z0), geo.V(0, 0, 1), geo.Wavelengthd)
		res := sys.Trace(r)
		if !res.OK() {
			return math.NaN()
		}
		st := res.Steps[len(res.Steps)-1]
		if math.Abs(st.DirOut.Z) < 1e-12 {
			return math.NaN()
		}
		t := (zFocus - st.Point.Z) / st.DirOut.Z
		return st.Point.Y + t*st.DirOut.Y
	}
	a4, err := solve1D(resid, -1e-3/(semi*semi), 1e-3/(semi*semi), 200, 120)
	if err != nil {
		return nil, fmt.Errorf("asphere solve: %v", err)
	}
	p.Cards = cards
	p.Cards[last].Coeffs = []float64{a4}
	if err := p.Validate(); err != nil {
		return nil, err
	}
	if _, err := p.System(); err != nil {
		return nil, err
	}
	return p, nil
}

func spanOf2(cards []lens.Card) float64 {
	var s float64
	for _, c := range cards {
		s += math.Abs(c.ThicknessMM)
	}
	return s
}

func buildSysFromCards(cards []lens.Card) (*geo.System, error) {
	return buildSys(cards)
}

// ---------- cemented achromats ----------

func achromatFilename(kind, crown, flint string, efl float64) string {
	return fmt.Sprintf("achromat-%s-f%s-%s-%s.json", kind, trimNum(math.Abs(efl), 3), glassCode(crown), glassCode(flint))
}

// achromatDesigned runs the package's automatic achromat designer.
func achromatDesigned(l *lib, crown, flint string, efl, wantFnum float64, fieldDeg float64) {
	fnum := fnumFor(efl, wantFnum)
	semi := efl / (2 * fnum)
	opts := &lens.DesignOptions{}
	if fieldDeg > 0 {
		opts.Fields = []lens.Field{{}, {TanY: math.Tan(fieldDeg * math.Pi / 180)}}
	}
	p, res, err := lens.DesignAchromat(efl, crown, flint, semi, opts)
	if err != nil {
		fmt.Printf("  ! achromat dls %s/%s f=%g: %v\n", crown, flint, efl, err)
		return
	}
	stats := ""
	if res != nil {
		stats = fmt.Sprintf("｜DLS 抛光：%d 次接受迭代／%d 次求值，merit=%.4g", res.Iter, res.Evals, res.Merit)
	}
	name := fmt.Sprintf("消色差双胶合 %s/%s f=%.4gmm F/%s（自主设计）", crown, flint, efl, trimNum(fnum, 2))
	p.Name = name
	l.add(achromatFilename("dls", crown, flint, efl), p, "achromat", 2,
		[]string{"achromat", "doublet", crown, flint, fnumTag(fnum), eflTag(efl), "visible",
			"optimized", "quality:optimized", fmt.Sprintf("fov%gdeg", fieldDeg), "sweep:efl"},
		fmt.Sprintf("双胶合消色差｜%s+%s｜f=%smm、F/%s，%g° 视场｜自主设计：薄透镜解→厚透镜牛顿消色差→DLS 抛光%s",
			crown, flint, trimNum(efl, 0), trimNum(fnum, 2), fieldDeg, stats))
}

// achromatSolved uses the same designer without the DLS polish: the
// closed-form thick-lens achromat, useful as a known-quality reference.
func achromatSolved(l *lib, crown, flint string, efl, wantFnum float64) {
	fnum := fnumFor(efl, wantFnum)
	semi := efl / (2 * fnum)
	p, _, err := lens.DesignAchromat(efl, crown, flint, semi, &lens.DesignOptions{SkipPolish: true})
	if err != nil {
		fmt.Printf("  ! achromat solve %s/%s f=%g f/%g: %v\n", crown, flint, efl, fnum, err)
		return
	}
	p.Name = fmt.Sprintf("消色差双胶合 %s/%s f=%.4gmm F/%s（解析解）", crown, flint, efl, trimNum(fnum, 2))
	file := achromatFilename("solve", crown, flint, efl)
	if math.Abs(fnum-5) > 1e-9 {
		file = fmt.Sprintf("achromat-solve-f%s-%s-%s-%s.json", trimNum(efl, 0), fnumCode(fnum), glassCode(crown), glassCode(flint))
	}
	l.add(file, p, "achromat", 2,
		[]string{"achromat", "doublet", crown, flint, fnumTag(fnum), eflTag(efl), "visible",
			"classic", "quality:analytic", "fov0.5deg", "sweep:efl"},
		fmt.Sprintf("双胶合消色差｜%s+%s｜f=%smm、F/%s｜解析构造（薄透镜消色差解 + 厚透镜牛顿修正，未经 DLS 抛光）",
			crown, flint, trimNum(efl, 0), trimNum(fnum, 2)))
}

// achromatReversed reverses a designed doublet (flint first). Reversing a
// prescription means mirroring it in z — the surfaces come in the opposite
// order and every radius changes sign — which keeps the focal length and the
// first-order colour correction and changes the higher-order aberrations.
func achromatReversed(l *lib, crown, flint string, efl, wantFnum float64) {
	fnum := fnumFor(efl, wantFnum)
	semi := efl / (2 * fnum)
	base, _, err := lens.DesignAchromat(efl, crown, flint, semi, &lens.DesignOptions{SkipPolish: true})
	if err != nil {
		fmt.Printf("  ! achromat rev %s/%s f=%g: %v\n", crown, flint, efl, err)
		return
	}
	n := len(base.Cards)
	cards := make([]lens.Card, 0, n)
	for j := 0; j < n; j++ {
		src := base.Cards[n-1-j]
		t := base.Cards[0].ThicknessMM
		g := ""
		if j < n-1 {
			t = base.Cards[n-2-j].ThicknessMM
			g = base.Cards[n-2-j].Glass
		}
		cards = append(cards, lens.Card{
			Name:        src.Name,
			RadiusMM:    -src.RadiusMM,
			ConicK:      src.ConicK,
			Coeffs:      src.Coeffs,
			ThicknessMM: t,
			Glass:       g,
			SemiAperMM:  r3(semi * 1.2),
		})
	}
	p := &lens.Prescription{
		Name:      fmt.Sprintf("反向消色差双胶合 %s/%s f=%.4gmm（%s 在前）", crown, flint, efl, flint),
		ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards,
	}
	if err := finishLens(p, semi, efl); err != nil {
		fmt.Printf("  ! achromat rev finish %s/%s f=%g: %v\n", crown, flint, efl, err)
		return
	}
	l.add(achromatFilename("rev", crown, flint, efl), p, "achromat", 3,
		[]string{"achromat", "doublet", "reversed", flint, crown, fnumTag(fnum), eflTag(efl), "visible",
			"classic", "quality:classic", "fov0.5deg", "sweep:orientation"},
		fmt.Sprintf("反向双胶合（火石在前）｜%s/%s｜f=%smm、F/%s｜经典构造：把 %s/%s 消色差对调方向（处方按 z 镜像反排、半径变号），焦距与一阶色差不变、高级像差改变",
			flint, crown, trimNum(efl, 0), trimNum(fnum, 2), crown, flint))
}

// nirAchromat builds a genuine two-colour achromat for the near infrared: the
// thin-lens power split uses Abbe numbers evaluated between 1064 nm and
// 1550 nm (from the same Sellmeier data), and the whole prescription is then
// scaled so that EFL(1064) is exact. The residual EFL(1550) difference is
// measured by both this generator and libcheck.
func nirAchromat(l *lib, crown, flint string, efl float64, fnum float64) {
	const l1, l2 = 1064.0, 1550.0
	semi := efl / (2 * fnum)
	g1, ok1 := geo.Glass(crown)
	g2, ok2 := geo.Glass(flint)
	if !ok1 || !ok2 {
		return
	}
	n1, n2 := g1.N(l1), g2.N(l1)
	v1 := (g1.N(l1) - 1) / (g1.N(l1) - g1.N(l2))
	v2 := (g2.N(l1) - 1) / (g2.N(l1) - g2.N(l2))
	if !(v1 > v2) {
		fmt.Printf("  ! nir achromat %s/%s: V(%g)=%.2f is not above V(%g)=%.2f\n", crown, flint, l1, v1, l1, v2)
		return
	}
	phi := 1 / efl
	phi1 := phi * v1 / (v1 + v2)
	phi2 := phi - phi1
	cA := phi1 / (2 * (n1 - 1))
	cB := -cA
	cC := cB - phi2/(n2-1)
	cards := []lens.Card{
		{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 3, SemiAperMM: r3(semi)},
		{Name: "crown front", RadiusMM: 1 / cA, ThicknessMM: math.Max(2, 0.03*efl), Glass: crown, SemiAperMM: r3(semi * 1.2)},
		{Name: "cemented", RadiusMM: 1 / cB, ThicknessMM: math.Max(1, 0.012*efl), Glass: flint, SemiAperMM: r3(semi * 1.2)},
		{Name: "flint rear", RadiusMM: 1 / cC, ThicknessMM: 0.5 * efl, Glass: "", SemiAperMM: r3(semi * 1.2)},
	}
	p := &lens.Prescription{Name: fmt.Sprintf("近红外消色差双胶合 %s/%s f=%smm", crown, flint, trimNum(efl, 0)),
		ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
	// Scale all curvatures together (which preserves the power split) so
	// that EFL(1064) lands exactly on the target.
	orig := []float64{1 / cA, 1 / cB, 1 / cC}
	eval := func(s float64) float64 {
		for k, i := range []int{1, 2, 3} {
			p.Cards[i].RadiusMM = orig[k] / s
		}
		par, err := p.Paraxial(l1)
		if err != nil || math.IsNaN(par.EFL) {
			return math.NaN()
		}
		return par.EFL - efl
	}
	s, err := solve1D(eval, 0.2, 5, 200, 120)
	if err != nil {
		fmt.Printf("  ! nir achromat %s/%s f=%g: %v\n", crown, flint, efl, err)
		return
	}
	for k, i := range []int{1, 2, 3} {
		p.Cards[i].RadiusMM = r4(orig[k] / s)
	}
	if err := autoApertures(p, semi); err != nil {
		fmt.Printf("  ! nir achromat apertures %s/%s f=%g: %v\n", crown, flint, efl, err)
		return
	}
	if err := p.Validate(); err != nil {
		fmt.Printf("  ! nir achromat validate %s/%s f=%g: %v\n", crown, flint, efl, err)
		return
	}
	pa, _ := p.Paraxial(l1)
	pb, _ := p.Paraxial(l2)
	delta := math.NaN()
	if pa != nil && pb != nil {
		delta = pa.EFL - pb.EFL
	}
	p.Name = fmt.Sprintf("近红外消色差双胶合 %s/%s f=%smm（1064/1550nm）", crown, flint, trimNum(efl, 0))
	file := fmt.Sprintf("achromat-nir-f%s-%s-%s.json", trimNum(efl, 0), glassCode(crown), glassCode(flint))
	l.add(file, p, "achromat", 4,
		[]string{"achromat", "doublet", "nir", crown, flint, fnumTag(fnum), eflTag(efl),
			"classic", "quality:analytic", "fov0.5deg", "sweep:band"},
		fmt.Sprintf("近红外消色差双胶合｜%s+%s｜f=%smm、F/%s｜经典构造：功率分配用 1064nm/1550nm 之间的阿贝数（同一 Sellmeier 数据）解析求出、再由近轴引擎按 EFL(1064)=目标 定标；1064 与 1550nm 的焦距差实测 %s（由 libcheck 复核）",
			crown, flint, trimNum(efl, 0), trimNum(fnum, 2), fmtSig(delta, "mm")))
	l.nir[file] = true
}

// doubletFinite builds a cemented doublet used at a finite conjugate (1:1
// imaging), where the chromatic and spherical balance differs from the
// infinite-conjugate design.
func doubletFinite(l *lib, crown, flint string, efl float64, ratio float64) {
	semi := math.Min(30, efl/10)
	g := grp{name: "L1", curv: seedAchromat(crown, flint, efl), glasses: []string{crown, flint, ""},
		thick: math.Max(1.5, 0.03*efl), after: 0.3 * efl}
	objZ := -2 * efl * ratio
	p, err := buildLayout(fmt.Sprintf("%s/%s 双胶合 1:%.3g 成像", crown, flint, ratio),
		[]grp{g}, 3, efl, semi, objZ)
	if err != nil {
		fmt.Printf("  ! doublet finite %s/%s: %v\n", crown, flint, err)
		return
	}
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		return
	}
	l.add(fmt.Sprintf("doublet-finite-m%s-%s-%s.json", trimNum(math.Abs(par.Magnification), 2), glassCode(crown), glassCode(flint)),
		p, "achromat", 3,
		[]string{"achromat", "doublet", "finite-conjugate", crown, flint, eflTag(efl), "visible",
			"classic", "quality:analytic", fmt.Sprintf("mag:%s", trimNum(par.Magnification, 3)),
			"objfield:" + trimNum(0.03*math.Abs(objZ), 2), "sweep:conjugate"},
		fmt.Sprintf("有限共轭双胶合（%s/%s）｜物距 %smm（≈1:%.3g 成像）、f=%smm｜经典构造：无限共轭消色差解的有限共轭使用，剩余像差由实测记录",
			crown, flint, trimNum(objZ, 1), ratio, trimNum(efl, 0)))
}

func buildAchromats(l *lib) {
	dlsEfl := []float64{40, 50, 60, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 750, 1000}
	pairs := [][2]string{{"N-BK7", "F2"}, {"N-SK16", "N-SF5"}}
	for _, pr := range pairs {
		for _, f := range dlsEfl {
			achromatDesigned(l, pr[0], pr[1], f, 5, 0.5)
		}
	}
	solveEfl := []float64{45, 55, 70, 90, 110, 140, 160, 220, 260, 320, 450, 550, 650, 800, 900, 1200}
	for _, pr := range [][2]string{{"N-BK7", "N-SF5"}, {"FusedSilica", "F2"}, {"FusedSilica", "N-SF5"}} {
		for _, f := range solveEfl {
			achromatSolved(l, pr[0], pr[1], f, 5)
		}
	}
	// F-number variants at f = 100
	for _, fn := range []float64{2.8, 4, 5.6, 8, 11, 16} {
		achromatSolved(l, "N-BK7", "F2", 100, fn)
	}
	// reversed (flint first)
	for _, f := range []float64{75, 100, 150, 200, 300, 500} {
		achromatReversed(l, "N-BK7", "F2", f, 5)
	}
	// near infrared two-colour achromats
	for _, f := range []float64{100, 150, 200, 250, 300, 400, 500} {
		nirAchromat(l, "FusedSilica", "N-SF5", f, 5)
	}
	// In the near infrared F2 disperses more gently than N-BK7, so the
	// crown/flint roles of the pair swap: F2 is the "crown" of the NIR pair.
	for _, f := range []float64{120, 350} {
		nirAchromat(l, "F2", "N-BK7", f, 5)
	}
	// finite-conjugate doublets
	for _, c := range []struct {
		f     float64
		ratio float64
	}{{50, 1}, {100, 2}, {100, 4}, {150, 6}} {
		doubletFinite(l, "N-BK7", "F2", c.f, c.ratio)
	}
}
