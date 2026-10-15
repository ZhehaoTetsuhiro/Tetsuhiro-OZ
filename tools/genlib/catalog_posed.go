// The catalogue: decentred and tilted systems (the coordinate-break fields of
// lens.Card). geo's first-order and third-order engines are defined for
// coaxial systems only and refuse a posed prescription, so every example here
// is measured exclusively on the real tracer: the library measurement falls
// back to the "posed" mode (best-focus plane of a parallel on-axis bundle).
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/lens"
)

// posedMinRays is the smallest number of the 21 traced fan rays that must
// reach the last surface for a posed example to be registered.
const posedMinRays = 16

// posedCheck measures p and applies the generator's own quality gate for the
// posed families: the measurement must be a posed/focus real-ray focus with
// enough rays and an RMS small against the aperture.
func posedCheck(file string, p *lens.Prescription, maxRMSSemi float64) (float64, bool) {
	m, err := measureOf(p)
	if err != nil {
		failf("%s: %v", file, err)
		return 0, false
	}
	if m.Mode != "posed" && m.Mode != "focus" {
		failf("%s: measured mode %s, want posed/focus (a posed system has no first-order optics)", file, m.Mode)
		return 0, false
	}
	if m.Rays < posedMinRays {
		failf("%s: only %d of 21 rays reach the last surface", file, m.Rays)
		return 0, false
	}
	if !(m.FocusRMSMM > 0) || m.FocusRMSMM > maxRMSSemi*m.StopSemiMM {
		failf("%s: focus RMS %.6gmm exceeds %.3g x the aperture %.6gmm", file, m.FocusRMSMM, maxRMSSemi, m.StopSemiMM)
		return 0, false
	}
	return m.FocusRMSMM, true
}

// tiltedMirror builds one mirror folded by a tilt about x (the meridional
// plane of the traced fan). The mirror radius is R = -2f for the paraxial
// focus, K selects the conic (0 sphere, -1 paraboloid). The description
// states the engine-measured fold of the exit bundle (180° − 2θ for a mirror
// tilted by θ) and the measured best-focus position of the folded beam.
func tiltedMirror(l *lib, key, zh string, f, wantFnum, tiltDeg, K float64) {
	fnum := fnumFor(f, wantFnum)
	semi := math.Min(120, f/(2*fnum))
	p := &lens.Prescription{Name: fmt.Sprintf("%s（倾斜 %s°）f=%smm F/%s", zh, trimNum(tiltDeg, 2), trimNum(f, 0), trimNum(fnum, 2)),
		ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
			{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 60, SemiAperMM: r3(semi)},
			{Name: "primary", RadiusMM: -2 * f, ConicK: K, Mirror: true, ThicknessMM: 200, SemiAperMM: r3(semi * 1.02), TiltXDeg: tiltDeg},
		}}
	if err := p.Validate(); err != nil {
		failf("posed mirror %s: %v", key, err)
		return
	}
	file := fmt.Sprintf("posed-mirror-%s-tilt%s-f%s-%s.json", key, numCode(tiltDeg), trimNum(f, 0), fnumCode(fnum))
	rms, ok := posedCheck(file, p, 0.25)
	if !ok {
		return
	}
	fold, err := exitAngleDeg(p)
	if err != nil {
		failf("posed mirror %s: %v", file, err)
		return
	}
	m, _ := measureOf(p)
	note := fmt.Sprintf("%s｜倾斜 θ=%s°（tilt_x_deg，坐标断点）｜f=%smm、F/%s、口径 %smm｜实光线追迹：平行光入射经倾斜镜折转，实测出射束轴与入射光轴夹角 %s°（薄镜几何折转 180°−2θ=%s°；镜面曲率与球差使实光束轴略偏）｜焦点由平行光束的实光线最佳焦面定位：RMS=%smm（一阶引擎对 pose 系统报错，指标全部来自追迹）",
		zh, trimNum(tiltDeg, 2), trimNum(f, 0), trimNum(fnum, 2), trimNum(semi, 3), trimNum(fold, 4),
		trimNum(180-2*tiltDeg, 4), trimNum(rms, 4))
	l.add(file, p, "tilt", 4, []string{
		"posed", "tilt", "mirror", key, "single-mirror", fnumTag(fnum), eflTag(f), "visible",
		"trace-validated", "quality:classic", "sweep:tilt", fmt.Sprintf("tilt:%sdeg", trimNum(tiltDeg, 2)),
		"metrics:" + m.Mode,
	}, note)
}

// ---------- the families of this catalogue ----------

func buildPosed(l *lib) {
	// tilted single mirrors: the fold grows with the tilt, the meridional fan
	// stays in the plane of the fan
	for _, t := range []float64{2, 5, 7, 10, 12, 15, 20, 25, 30, 35, 40, 45} {
		tiltedMirror(l, "sphere", "球面倾斜反射镜", 200, 4, t, 0)
	}
	for _, t := range []float64{5, 10, 15, 20, 30, 45} {
		tiltedMirror(l, "sphere", "球面倾斜反射镜", 100, 4, t, 0)
	}
	for _, t := range []float64{10, 20, 30} {
		tiltedMirror(l, "sphere", "球面倾斜反射镜", 400, 4, t, 0)
	}
	for _, t := range []float64{10, 20, 30} {
		tiltedMirror(l, "sphere", "球面倾斜反射镜", 200, 8, t, 0)
	}
	for _, t := range []float64{15, 30} {
		tiltedMirror(l, "sphere", "球面倾斜反射镜", 300, 5.6, t, 0)
	}
	tiltedMirror(l, "sphere", "球面倾斜反射镜", 100, 8, 20, 0)
	for _, t := range []float64{10, 20, 30, 45} {
		tiltedMirror(l, "parabola", "抛物面倾斜反射镜", 200, 4, t, -1)
	}
	for _, t := range []float64{15, 30} {
		tiltedMirror(l, "parabola", "抛物面倾斜反射镜", 100, 5, t, -1)
	}
	tiltedMirror(l, "parabola", "抛物面倾斜反射镜", 400, 5, 20, -1)
	tiltedMirror(l, "parabola", "抛物面倾斜反射镜", 300, 8, 15, -1)
	// decentred singlets: the lateral shift of the whole element
	for _, dy := range []float64{0.05, 0.1, 0.25, 0.5, 1, 2} {
		posedSinglet(l, "N-BK7", 100, 4, 0, dy, 0, 0)
	}
	for _, dx := range []float64{0.1, 0.25, 0.5, 1} {
		posedSinglet(l, "N-BK7", 100, 4, dx, 0, 0, 0)
	}
	for _, dy := range []float64{0.25, 0.5, 1} {
		posedSinglet(l, "N-BK7", 200, 5.6, 0, dy, 0, 0)
	}
	for _, dy := range []float64{0.1, 0.25} {
		posedSinglet(l, "N-BK7", 50, 2.8, 0, dy, 0, 0)
	}
	posedSinglet(l, "FusedSilica", 100, 4, 0, 0.5, 0, 0)
	posedSinglet(l, "N-SK16", 100, 4, 0.25, 0, 0, 0)
	// tilted singlets
	for _, tx := range []float64{0.25, 0.5, 1, 2, 5, 10} {
		posedSinglet(l, "N-BK7", 100, 4, 0, 0, tx, 0)
	}
	for _, ty := range []float64{0.5, 1, 3} {
		posedSinglet(l, "N-BK7", 100, 4, 0, 0, 0, ty)
	}
	for _, tx := range []float64{1, 2, 5} {
		posedSinglet(l, "N-BK7", 200, 5.6, 0, 0, tx, 0)
	}
	posedSinglet(l, "N-BK7", 50, 2.8, 0, 0, 1, 0)
	posedSinglet(l, "N-SK16", 100, 4, 0, 0, 0, 2)
	// plane-parallel plate with a pose, in the parallel and in the converging beam
	for _, t := range []float64{2, 5, 10, 20, 30} {
		posedPlate(l, "N-BK7", 100, 4, t, 0, "front")
	}
	for _, t := range []float64{5, 10, 20, 30, 45} {
		posedPlate(l, "N-BK7", 100, 4, t, 0, "conv")
	}
	for _, dx := range []float64{0.5, 1, 2} {
		posedPlate(l, "N-BK7", 100, 4, 0, dx, "front")
	}
	posedPlate(l, "N-BK7", 100, 4, 0, 1, "conv")
	posedPlate(l, "N-BK7", 100, 4, 10, 0.5, "front")
	posedPlate(l, "N-BK7", 100, 4, 20, 1, "front")
	for _, t := range []float64{10, 20} {
		posedPlate(l, "FusedSilica", 200, 5.6, t, 0, "conv")
	}
	posedPlate(l, "N-SK16", 100, 4, 10, 0, "front")
	// folded arm: flat mirror + the singlet rotated onto the folded beam
	for _, t := range []float64{20, 30, 45} {
		tiltedFlatMirrorAndLens(l, "N-BK7", 100, 4, t)
	}
	tiltedFlatMirrorAndLens(l, "N-BK7", 200, 5, 45)
	tiltedFlatMirrorAndLens(l, "N-BK7", 50, 2.8, 45)
}

// posedSinglet poses both surfaces of a plano-convex singlet with the same
// decentre or tilt (the element moves as a whole; the gap between its two
// surfaces stays on the unposed axis, which the description states).
func posedSinglet(l *lib, glass string, fl, wantFnum, decentreX, decentreY, tiltX, tiltY float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	p, err := buildLayout(fmt.Sprintf("偏心/倾斜单透镜 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2)),
		[]grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		failf("posed singlet %s f=%g: %v", glass, fl, err)
		return
	}
	for i := 1; i < len(p.Cards); i++ {
		p.Cards[i].DecenterXMM = decentreX
		p.Cards[i].DecenterYMM = decentreY
		p.Cards[i].TiltXDeg = tiltX
		p.Cards[i].TiltYDeg = tiltY
	}
	what := ""
	var tag string
	switch {
	case decentreX != 0 && decentreY == 0:
		what, tag = fmt.Sprintf("横向偏心 x=%smm", trimNum(decentreX, 3)), "decenter-x"
	case decentreY != 0 && decentreX == 0:
		what, tag = fmt.Sprintf("横向偏心 y=%smm", trimNum(decentreY, 3)), "decenter-y"
	case decentreX != 0 && decentreY != 0:
		what, tag = fmt.Sprintf("横向偏心 (%s, %s)mm", trimNum(decentreX, 3), trimNum(decentreY, 3)), "decenter-xy"
	case tiltX != 0:
		what, tag = fmt.Sprintf("倾斜 x 轴 %s°", trimNum(tiltX, 3)), "tilt-x"
	default:
		what, tag = fmt.Sprintf("倾斜 y 轴 %s°", trimNum(tiltY, 3)), "tilt-y"
	}
	code := fmt.Sprintf("dcx%s-dcy%s-tx%s-ty%s", numCode(decentreX), numCode(decentreY), numCode(tiltX), numCode(tiltY))
	file := fmt.Sprintf("posed-singlet-%s-f%s-%s-%s.json", code, trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	rms, ok := posedCheck(file, p, 0.1)
	if !ok {
		return
	}
	m, err := measureOf(p)
	if err != nil {
		failf("posed singlet %s: %v", file, err)
		return
	}
	p.Name = fmt.Sprintf("偏心/倾斜单透镜 %s f=%.4gmm F/%s（%s）", glass, efl, trimNum(fnum, 2), what)
	note := fmt.Sprintf("偏心/倾斜单透镜｜%s｜f=%smm、F/%s｜%s（decenter_*_mm / tilt_*_deg，坐标断点；两面同 pose，面间距仍沿原光轴）｜pose 系统无一阶量：焦点由平行光束的实光线最佳焦面定位，实测最佳焦面 z=%smm、RMS=%smm｜像差（彗差/像散）随偏置或倾角增长，可与同一 f 数字的轴对称例对照",
		glass, trimNum(efl, 0), trimNum(fnum, 2), what, trimNum(m.FocusZMM, 5), trimNum(rms, 4))
	l.add(file, p, "tilt", 4, []string{
		"posed", tag, "singlet", glass, fnumTag(fnum), eflTag(efl), "visible",
		"trace-validated", "quality:classic", "sweep:" + tag, "mode-witness:posed-real-rays",
	}, note)
}

// posedPlate inserts a plane-parallel plate into the beam of the same
// plano-convex singlet. pos = "front" puts the plate into the parallel beam
// in front of the lens, pos = "conv" puts it into the converging beam behind
// it. Plate and lens carry the pose; the plate is a two-surface element with
// the same pose on both surfaces.
func posedPlate(l *lib, glass string, fl, wantFnum, tiltDeg, decenter float64, pos string) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	p, err := buildLayout(fmt.Sprintf("倾斜平板 + 单透镜 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2)),
		[]grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		failf("posed plate %s f=%g: %v", glass, fl, err)
		return
	}
	plate := []lens.Card{
		{Name: "plate front", RadiusMM: math.Inf(1), ThicknessMM: 3, Glass: glass, SemiAperMM: r3(1.2 * semi), TiltXDeg: tiltDeg, DecenterXMM: decenter},
		{Name: "plate rear", RadiusMM: math.Inf(1), ThicknessMM: 5, SemiAperMM: r3(1.2 * semi), TiltXDeg: tiltDeg, DecenterXMM: decenter},
	}
	stop := p.Cards[0]
	lensCards := p.Cards[1:]
	var cards []lens.Card
	switch pos {
	case "front":
		cards = append([]lens.Card{stop}, plate[0], plate[1])
		cards = append(cards, lensCards...)
	case "conv":
		// the plate sits 30 mm behind the lens (or the mid-gap, whichever is
		// smaller) so that the beam radius at the plate stays useful
		cards = append([]lens.Card{stop}, lensCards[0])
		rear := lensCards[1]
		gap := math.Min(0.3*efl, rear.ThicknessMM)
		rear.ThicknessMM = gap
		cards = append(cards, rear, plate[0], plate[1])
	default:
		failf("posed plate %s: unknown position %q", glass, pos)
		return
	}
	q := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
	file := fmt.Sprintf("posed-plate-%s-tilt%s-dcx%s-f%s-%s-%s.json", pos, numCode(tiltDeg), numCode(decenter), trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	rms, ok := posedCheck(file, q, 0.1)
	if !ok {
		return
	}
	m, err := measureOf(q)
	if err != nil {
		failf("posed plate %s: %v", file, err)
		return
	}
	zhn := "平行光中"
	if pos == "conv" {
		zhn = "会聚光中"
	}
	q.Name = fmt.Sprintf("倾斜平板（%s）+ 单透镜 %s f=%.4gmm F/%s（tilt %s°）", zhn, glass, efl, trimNum(fnum, 2), trimNum(tiltDeg, 2))
	note := fmt.Sprintf("倾斜平板 + 球面单透镜｜%s｜f=%smm、F/%s｜%s毫米厚平板置于%s（tilt_x_deg=%s、decenter_x_mm=%s，坐标断点；平板两面同 pose）｜平板使光束平移（会聚光中另引入像散）：实测最佳焦面 z=%smm、RMS=%smm｜一阶引擎对 pose 系统报错，指标全部来自实光线追迹",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(3, 0), zhn, trimNum(tiltDeg, 2), trimNum(decenter, 3),
		trimNum(m.FocusZMM, 5), trimNum(rms, 4))
	l.add(file, q, "tilt", 3, []string{
		"posed", "tilt", "plate", "singlet", glass, fnumTag(fnum), eflTag(efl), "visible",
		"trace-validated", "quality:classic", "sweep:tilt", "plane-parallel",
	}, note)
}

// tiltedFlatMirrorAndLens folds a collimated beam by 45° with a plane mirror
// and then focuses it with the plano-convex singlet placed on the folded arm:
// the lens card carries the coordinate break that rotates its axis onto the
// folded beam. This is the smallest complete folded layout of the sequential
// engine (the mirror folds, the next surface is placed with its own pose).
func tiltedFlatMirrorAndLens(l *lib, glass string, fl, wantFnum, tiltDeg float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	p, err := buildLayout(fmt.Sprintf("折叠光路 + 单透镜 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2)),
		[]grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		failf("fold+lens %s f=%g: %v", glass, fl, err)
		return
	}
	const gap = 40.0 // mirror to lens along the folded arm
	thick := thickFor(efl)
	front, rear := p.Cards[1], p.Cards[2]
	// A mirror tilted by θ about x sends the incoming +z beam along
	// d = (0, sin 2θ, −cos 2θ). The folded arm is placed by the poses of its
	// surfaces: the decentre carries the y part of the arm, the (possibly
	// negative) chain thickness of the preceding card carries its z part
	// (lens.Card has no axial pose parameter), and the lens axis is rotated
	// onto the beam by tilt_x = 2θ − 180°.
	rad := 2 * tiltDeg * math.Pi / 180
	sin, cos := math.Sin(rad), math.Cos(rad)
	axis := 2*tiltDeg - 180
	q := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: []lens.Card{
		p.Cards[0],
		{Name: "fold mirror", RadiusMM: 0, Mirror: true, ThicknessMM: r6(-gap * cos), SemiAperMM: r3(1.5 * semi), TiltXDeg: tiltDeg},
		{Name: "L1 (folded arm)", RadiusMM: front.RadiusMM, ThicknessMM: r6(-thick * cos), Glass: front.Glass,
			SemiAperMM: r3(1.5 * semi), DecenterYMM: r6(gap * sin), TiltXDeg: axis},
		{Name: "L1 rear (folded arm)", RadiusMM: rear.RadiusMM, ThicknessMM: 0,
			SemiAperMM: r3(1.5 * semi), DecenterYMM: r6((gap + thick) * sin), TiltXDeg: axis},
	}}
	file := fmt.Sprintf("posed-fold-tilt%s-lens-f%s-%s-%s.json", numCode(tiltDeg), trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	rms, ok := posedCheck(file, q, 0.1)
	if !ok {
		return
	}
	fold, err := foldAtCardDeg(q, 1)
	if err != nil {
		failf("fold+lens %s: %v", file, err)
		return
	}
	// the fold of a plane mirror tilted by θ is exactly 180° − 2θ: the angle
	// read off the traced steps at the mirror is checked against the closed
	// form
	wantFold := 180 - 2*tiltDeg
	if math.Abs(fold-wantFold) > 1e-9 {
		failf("fold+lens %s: measured fold %.9g° differs from 180°-2θ = %.9g°", file, fold, wantFold)
		return
	}
	m, _ := measureOf(q)
	q.Name = fmt.Sprintf("折转光路（%s° 折轴镜）+ 折叠臂单透镜 %s f=%.4gmm F/%s",
		trimNum(wantFold, 0), glass, efl, trimNum(fnum, 2))
	note := fmt.Sprintf("折叠光路：平面折轴镜 + 折叠臂单透镜｜%s｜f=%smm、F/%s｜折轴镜倾斜 θ=%s°（tilt_x_deg），折叠臂长 %smm，透镜两面用 pose 放到折叠臂上（decenter_y_mm=%s、%s，tilt_x_deg=%s；链厚承接折叠臂的 z 分量）｜实光线在折轴镜上的折转 %s°（解析 180°−2θ=%s°，实测与解析一致）｜最佳焦面 z=%smm、RMS=%smm｜全部指标来自追迹（一阶引擎对 pose 系统报错）",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(tiltDeg, 2), trimNum(gap, 0),
		trimNum(gap*sin, 5), trimNum((gap+thick)*sin, 5), trimNum(axis, 2),
		trimNum(fold, 4), trimNum(wantFold, 4),
		trimNum(m.FocusZMM, 5), trimNum(rms, 4))
	l.add(file, q, "tilt", 5, []string{
		"posed", "tilt", "fold", "flat-mirror", "singlet", glass, fnumTag(fnum), eflTag(efl),
		"visible", "trace-validated", "quality:classic",
	}, note)
}
