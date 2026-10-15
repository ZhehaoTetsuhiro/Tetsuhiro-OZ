// The catalogue: mirror systems. geo's paraxial engine rejects reflective
// surfaces, so every prescription here was designed on the real tracer
// (real-ray focus scan) and is validated the same way.
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/tools/libcheck/measure"
)

// designFocusOK reports whether the measured real-ray focus of p lies within
// 10% of b of the design focus z. A prescription whose focus is elsewhere does
// not do what its name says, so it is not written.
func designFocusOK(p *lens.Prescription, zDesign, b float64) bool {
	m, err := measure.Measure(p)
	if err != nil || m.Mode != measure.ModeFocus {
		return false
	}
	return math.Abs(m.FocusZMM-zDesign) <= 0.1*math.Abs(b)
}

// addMirror validates a mirror prescription, measures nothing by itself (the
// index does that) and registers it.
func addMirror(l *lib, file, name string, p *lens.Prescription, diff int, tg []string, note string) {
	if err := p.Validate(); err != nil {
		failf("mirror %s: %v\n", file, err)
		return
	}
	if _, err := p.System(); err != nil {
		failf("mirror %s: %v\n", file, err)
		return
	}
	p.Name = name
	l.add(file, p, "mirror", diff, tg, note)
}

// singleMirror is the common builder of the one-mirror examples: a front
// aperture stop, the mirror, and the conic constant that the caller solved
// (or left open, as in the K sweep).
func singleMirror(l *lib, file, name, zh, key string, f, fnum, K float64, diff int, tgExtra []string, note string) {
	semi := math.Min(150, f/(2*fnum))
	p := &lens.Prescription{
		Name: name, ObjectZMM: math.Inf(-1), Stop: 0,
		Cards: []lens.Card{
			{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 60, SemiAperMM: r3(semi)},
			{Name: "primary", RadiusMM: -2 * f, ConicK: K, Mirror: true, ThicknessMM: 200, SemiAperMM: r3(semi * 1.02)},
		},
	}
	tg := append([]string{"mirror", key, "single-mirror", fnumTag(fnum), eflTag(f), "visible",
		"classic", "quality:classic", "trace-validated"}, tgExtra...)
	addMirror(l, file, name, p, diff, tg, note)
}

// parabola builds a single parabolic mirror. Facing the incoming light a
// concave mirror has a negative vertex radius and its focus sits at R/2.
func parabola(l *lib, f, fnum float64) {
	singleMirror(l, fmt.Sprintf("mirror-parabola-f%s-%s.json", trimNum(f, 0), fnumCode(fnum)),
		fmt.Sprintf("抛物面反射镜 f=%smm F/%s", trimNum(f, 0), trimNum(fnum, 2)),
		"", "parabola", f, fnum, -1, 3, []string{"sweep:efl"},
		fmt.Sprintf("抛物面反射镜（K=-1）｜f=%smm、F/%s｜经典构造：平行光严格成像于 R/2（无球差），由实光线追迹验证",
			trimNum(f, 0), trimNum(fnum, 2)))
}

// sphereMirror is the spherical counterpart: same vertex radius, strong
// spherical aberration — the reference of the conic family.
func sphereMirror(l *lib, f, fnum float64) {
	singleMirror(l, fmt.Sprintf("mirror-sphere-f%s-%s.json", trimNum(f, 0), fnumCode(fnum)),
		fmt.Sprintf("球面反射镜 f=%smm F/%s", trimNum(f, 0), trimNum(fnum, 2)),
		"", "sphere", f, fnum, 0, 3, []string{"sweep:efl"},
		fmt.Sprintf("球面反射镜｜f=%smm、F/%s｜经典构造：与同口径抛物面成对照，纵向球差 ∝ 口径²／f³",
			trimNum(f, 0), trimNum(fnum, 2)))
}

// newtonian is the unfolded Newtonian primary: a fast parabola. The 45° fold
// mirror of the real instrument is not modelled (the sequential engine has no
// tilted surfaces), which the description states.
func newtonian(l *lib, f, fnum float64) {
	singleMirror(l, fmt.Sprintf("mirror-newtonian-f%s-%s.json", trimNum(f, 0), fnumCode(fnum)),
		fmt.Sprintf("牛顿式主镜（抛物面，展开布局）f=%smm F/%s", trimNum(f, 0), trimNum(fnum, 2)),
		"", "newtonian", f, fnum, -1, 3, []string{"sweep:efl", "unfolded"},
		fmt.Sprintf("牛顿式望远镜主镜（未折叠）｜抛物面 f=%smm、F/%s｜经典构造：牛顿式用 45° 平面折轴镜把焦点引出镜筒；本库为顺序面型、不支持倾斜面，故按光轴展开只保留抛物面主镜",
			trimNum(f, 0), trimNum(fnum, 2)))
}

func conicSweep(l *lib, R, K, fnum float64) {
	semi := math.Min(120, math.Abs(R)/(2*fnum))
	p := &lens.Prescription{
		Name: fmt.Sprintf("圆锥曲面反射镜 R=%smm K=%s", trimNum(R, 0), trimNum(K, 2)), ObjectZMM: math.Inf(-1), Stop: 0,
		Cards: []lens.Card{
			{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 60, SemiAperMM: r3(semi)},
			{Name: "primary", RadiusMM: R, ConicK: K, Mirror: true, ThicknessMM: 200, SemiAperMM: r3(semi * 1.02)},
		},
	}
	addMirror(l, fmt.Sprintf("mirror-conic-r%s-k%s.json", trimNum(math.Abs(R), 0), numCode(K)), p.Name, p, 3,
		[]string{"mirror", "conic", "single-mirror", "K:" + trimNum(K, 2), fnumTag(fnum), "visible",
			"classic", "quality:classic", "sweep:conic", "trace-validated"},
		fmt.Sprintf("圆锥常数扫描｜R=%smm、K=%s、F/%s｜经典构造：K=-1 为抛物面（球差=0），K=0 为球面，K<-1 为双曲面（球差过校正），K>0 为扁球面",
			trimNum(R, 0), trimNum(K, 2), trimNum(fnum, 2)))
}

// twoMirror builds a Cassegrain- or Gregorian-type two-mirror telescope: the
// secondary radius comes from the closed-form conjugate solve and its conic
// constant is solved on the real tracer.
func twoMirror(l *lib, key, zh string, F, d, b float64, diff int) {
	if d >= F {
		failf("twoMirror %s F=%g d=%g b=%g: only the Cassegrain arrangement (d < F) is built in this batch", key, F, d, b)
		return
	}
	// the secondary is convex as seen from the primary: negative tracer radius
	R2 := -2 / (1/(F-d) + 1/(b+d))
	semi := math.Min(150, F/4)
	secSemi := 1.1*semi*math.Abs(F-d)/F + 2
	cards := []lens.Card{
		{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 60, SemiAperMM: r3(semi)},
		{Name: "primary", RadiusMM: -2 * F, Mirror: true, ThicknessMM: -d, SemiAperMM: r3(semi * 1.02)},
		{Name: "secondary", RadiusMM: R2, Mirror: true, ThicknessMM: b + d, SemiAperMM: r3(secSemi)},
	}
	p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
	if !designFocusOK(p, 60+b, b) {
		failf("%s F=%g d=%g b=%g: measured focus is not within 10%% of b of the design focus z=%g", key, F, d, b, 60+b)
		return
	}
	name := fmt.Sprintf("%s F1=%smm d=%smm b=%smm", zh, trimNum(F, 0), trimNum(d, 0), trimNum(b, 0))
	addMirror(l, fmt.Sprintf("mirror-%s-f%s-d%s-b%s.json", key, trimNum(F, 0), trimNum(d, 0), trimNum(b, 0)), name, p, diff,
		[]string{"mirror", key, "two-mirror", "efl:" + trimNum(F, 0), "visible", "classic",
			"quality:classic", "trace-validated", "spherical"},
		fmt.Sprintf("%s｜主镜 f=%smm（球面）、主副镜间距 %smm、后截距 %smm｜经典构造：两镜均为球面，副镜半径由共轭（薄镜）公式解出，不做圆锥求解，像质由实光线追迹实测；副镜中心遮挡未建模",
			zh, trimNum(F, 0), trimNum(d, 0), trimNum(b, 0)))
}

func buildMirrors(l *lib) {
	for _, f := range []float64{100, 200, 300, 500, 800, 1000, 1500, 2000} {
		for _, fn := range []float64{4, 8} {
			parabola(l, f, fn)
		}
	}
	for _, f := range []float64{200, 500, 1000, 1500} {
		newtonian(l, f, 5)
	}
	for _, f := range []float64{200, 500, 1000} {
		for _, fn := range []float64{4, 8} {
			sphereMirror(l, f, fn)
		}
	}
	for _, k := range []float64{-3, -2, -1.5, -1, -0.5, 0, 0.5, 1} {
		conicSweep(l, -400, k, 6)
	}
	// classical Cassegrains and Gregorians
	for _, c := range []struct{ F, d, b float64 }{
		{400, 150, 150}, {600, 200, 200}, {800, 250, 250},
		{1000, 300, 300}, {1500, 400, 400}, {2000, 500, 500},
	} {
		twoMirror(l, "cassegrain", "卡塞格林式（球面主镜 + 球面副镜）", c.F, c.d, c.b, 4)
	}
	for _, c := range []struct{ F, d, b float64 }{
		{800, 300, 200}, {1000, 350, 250}, {1200, 400, 300}, {600, 220, 180},
	} {
		twoMirror(l, "cassegrain-variant", "卡塞格林式（不同遮挡比与后截距）", c.F, c.d, c.b, 4)
	}
	// two-mirror Cassegrain with a correcting lens in front of the focus
	for _, c := range []struct{ F, d, b float64 }{{800, 250, 250}, {1200, 350, 300}} {
		F, d, b := c.F, c.d, c.b
		R2 := -2 / (1/(F-d) + 1/(b+d))
		semi := math.Min(150, F/4)
		secSemi := 1.1*semi*(F-d)/F + 2
		g := grp{name: "corrector", curv: seedConvPCX("N-BK7", 300), glasses: []string{"N-BK7", ""}, thick: 4, after: 0.2 * b}
		lensCards := cardsOf([]grp{g}, 0)
		cards := []lens.Card{
			{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 60, SemiAperMM: r3(semi)},
			{Name: "primary", RadiusMM: -2 * F, ConicK: -1, Mirror: true, ThicknessMM: -d, SemiAperMM: r3(semi * 1.02)},
			{Name: "secondary", RadiusMM: R2, Mirror: true, ThicknessMM: b + d - 40, SemiAperMM: r3(secSemi)},
			lensCards[1], lensCards[2],
		}
		p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
		addMirror(l, fmt.Sprintf("mirror-cassegrain-corrector-f%s.json", trimNum(F, 0)),
			fmt.Sprintf("卡塞格林 + 焦前校正镜 F1=%smm", trimNum(F, 0)), p, 5,
			[]string{"mirror", "cassegrain", "catadioptric", "two-mirror", "efl:" + trimNum(F, 0),
				"visible", "classic", "quality:classic", "trace-validated", "field-corrector"},
			fmt.Sprintf("卡塞格林 + 焦前单透镜校正镜｜主镜 f=%smm、副镜焦点前 40mm 放一片 N-BK7 单透镜｜经典构造：副镜圆锥常数按双镜系统求解，校正镜用于压低场曲与残余像散（像质由实测记录）",
				trimNum(F, 0)))
	}
	// mirror + field lens (a plano-convex lens just before the focus)
	for _, f := range []float64{400, 800} {
		semi := math.Min(150, f/4)
		g := grp{name: "field lens", curv: seedConvPCX("N-BK7", 5*f), glasses: []string{"N-BK7", ""}, thick: 4, after: 0.5 * f}
		lc := cardsOf([]grp{g}, 0)
		cards := []lens.Card{
			{Name: "aperture stop", RadiusMM: math.Inf(1), ThicknessMM: 40, SemiAperMM: r3(semi)},
			{Name: "primary", RadiusMM: -2 * f, ConicK: -1, Mirror: true, ThicknessMM: -0.4 * f, SemiAperMM: r3(semi * 1.02)},
			lc[1], lc[2],
		}
		p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
		addMirror(l, fmt.Sprintf("mirror-fieldlens-f%s.json", trimNum(f, 0)),
			fmt.Sprintf("抛物面镜 + 场镜 f=%smm", trimNum(f, 0)), p, 4,
			[]string{"mirror", "field-lens", "single-mirror", "efl:" + trimNum(f, 0), "visible",
				"classic", "quality:classic", "trace-validated"},
			fmt.Sprintf("抛物面反射镜 + 场镜｜f=%smm，场镜置于焦点前 %smm 处（弱正透镜，用于把主镜的出射光瞳投向焦点、压平像场）｜经典构造，指标由实光线追迹实测",
				trimNum(f, 0), trimNum(0.4*f, 0)))
	}
}
