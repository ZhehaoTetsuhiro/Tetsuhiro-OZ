// The catalogue: aspheric families. Every design here starts from one of the
// spherical builders of catalog_lens.go (or the package's achromat designer)
// and adds a polynomial (poly_terms) or Forbes Q (q_asphere) departure whose
// coefficients are fitted by the bounded coordinate search of coeffsolve.go
// to the engine's own on-axis spot RMS. A candidate is registered only when
// the fitted surface really lowers the measured RMS of its spherical base.
package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/tools/libcheck/measure"
)

// asphereFamily is one polynomial family: the powers of its terms, the number
// of coordinate-search passes, and the largest RMS ratio (asphere/sphere)
// that still counts as a real improvement.
type asphereFamily struct {
	key      string
	zh       string
	powers   []int
	passes   int
	maxRatio float64
}

// solvedAsphereSinglet builds a plano-convex (convex first) singlet with the
// same layout idiom as the spherical singlet family, measures its on-axis
// spot RMS, fits the given polynomial powers on its convex front surface and
// returns the finished prescription together with the spherical RMS, the
// aspheric RMS after freezing the coefficients to six significant digits, and
// the marginal height at the aspherised surface.
func solvedAsphereSinglet(glass string, fl, wantFnum float64, powers []int, passes int) (*lens.Prescription, float64, float64, float64, error) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	name := fmt.Sprintf("多项式非球面单透镜 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2))
	p, err := buildLayout(name, []grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		return nil, 0, 0, 0, err
	}
	r0, err := spotOnAxis(p)
	if err != nil {
		return nil, 0, 0, 0, fmt.Errorf("spherical measurement: %w", err)
	}
	const card = 1
	mh, err := marginalHeightAt(p, card)
	if err != nil {
		return nil, 0, 0, 0, err
	}
	polyTermsOf(p, card, powers)
	var vs []coeffVar
	var br []bracket
	for _, pw := range powers {
		vs = append(vs, polyVar(p, card, pw))
		sc := 0.05 / math.Pow(mh, float64(pw))
		br = append(br, bracket{lo: -sc, hi: sc})
	}
	if _, err := solveCoeffs(p, vs, br, passes); err != nil {
		return nil, 0, 0, 0, err
	}
	// freeze the coefficients to six significant digits
	for i := range p.Cards[card].PolyTerms {
		p.Cards[card].PolyTerms[i].Coeff = sig6(p.Cards[card].PolyTerms[i].Coeff)
	}
	if err := autoApertures(p, semi); err != nil {
		return nil, 0, 0, 0, err
	}
	if err := p.Validate(); err != nil {
		return nil, 0, 0, 0, err
	}
	r1, err := spotOnAxis(p)
	if err != nil {
		return nil, 0, 0, 0, fmt.Errorf("aspheric measurement: %w", err)
	}
	return p, r0, r1, mh, nil
}

// termsZh renders the powers as "r⁴+r⁶" style text for the descriptions.
func termsZh(powers []int) string {
	s := ""
	for i, pw := range powers {
		if i > 0 {
			s += "+"
		}
		s += fmt.Sprintf("r%s", supDigits(pw))
	}
	return s
}

// supDigits renders an integer with Unicode superscript digits (readability
// of the Chinese descriptions only; the file names stay ASCII).
func supDigits(n int) string {
	sup := []rune("⁰¹²³⁴⁵⁶⁷⁸⁹")
	if n < 0 {
		return fmt.Sprint(n)
	}
	out := ""
	for _, r := range fmt.Sprint(n) {
		out += string(sup[r-'0'])
	}
	return out
}

// asphereSinglet registers one polynomial-asphere singlet: the spherical and
// the aspheric RMS are both engine numbers and the comparison is stated in
// the description, which is the point of the family (spherical aberration
// drops by the measured factor).
func asphereSinglet(l *lib, glass string, fl, wantFnum float64, ap asphereFamily) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	file := fmt.Sprintf("asphere-poly-%s-f%s-%s-%s.json", ap.key, trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	p, r0, r1, mh, err := solvedAsphereSinglet(glass, fl, wantFnum, ap.powers, ap.passes)
	if err != nil {
		failf("asphere %s: %v", file, err)
		return
	}
	if !(r1 > 0) || !(r1 <= ap.maxRatio*r0) {
		failf("%s: aspheric RMS %.6gmm does not improve on the spherical %.6gmm (ratio %.4g > %.4g)",
			file, r1, r0, r1/r0, ap.maxRatio)
		return
	}
	p.Name = fmt.Sprintf("多项式非球面单透镜 %s f=%.4gmm F/%s（%s）", glass, efl, trimNum(fnum, 2), ap.zh)
	note := fmt.Sprintf("多项式非球面单透镜（poly_terms）｜%s｜f=%smm、F/%s、口径 %smm｜球差修正：球面基底实测轴上 RMS=%smm，前表面加 %s 项（系数由引擎点列 RMS 坐标搜索拟合、六位有效数字冻结）后 RMS=%smm，降至 %.1f%%｜该面边缘光线高度 %smm（非球面项的量纲参考）",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(efl/fnum, 3), trimNum(r0, 4), termsZh(ap.powers),
		trimNum(r1, 4), 100*r1/r0, trimNum(mh, 4))
	l.add(file, p, "asphere", 2, []string{
		"asphere", "poly-terms", "singlet", glass, fnumTag(fnum), eflTag(efl), "visible",
		"optimized", "quality:optimized", "compare:sphere", "sweep:efl", termsTag(ap.powers),
	}, note)
}

// termsTag renders the powers as an ASCII tag ("a4", "a4a6a8").
func termsTag(powers []int) string {
	s := ""
	for _, pw := range powers {
		s += fmt.Sprintf("a%d", pw)
	}
	return s
}

// achromatAsphere adds a polynomial asphere to one surface of the package's
// DLS-designed achromat and registers the corrected doublet. The base is the
// same design the achromat family records, so the description can state the
// measured gain of the aspheric correction on a colour-corrected base.
func achromatAsphere(l *lib, crown, flint string, efl, wantFnum float64, surf int, ap asphereFamily) {
	fnum := fnumFor(efl, wantFnum)
	semi := efl / (2 * fnum)
	opts := &lens.DesignOptions{Fields: []lens.Field{{}, {TanY: math.Tan(0.5 * math.Pi / 180)}}}
	p, _, err := lens.DesignAchromat(efl, crown, flint, semi, opts)
	if err != nil {
		failf("achromat asphere %s/%s f=%g F/%s s%d: %v", crown, flint, efl, trimNum(fnum, 2), surf, err)
		return
	}
	if surf < 0 || surf >= len(p.Cards) {
		failf("achromat asphere %s/%s f=%g F/%s: surface %d out of range", crown, flint, efl, trimNum(fnum, 2), surf)
		return
	}
	r0, err := spotOnAxis(p)
	if err != nil {
		failf("achromat asphere %s/%s f=%g F/%s: %v", crown, flint, efl, trimNum(fnum, 2), err)
		return
	}
	mh, err := marginalHeightAt(p, surf)
	if err != nil {
		failf("achromat asphere %s/%s f=%g F/%s: %v", crown, flint, efl, trimNum(fnum, 2), err)
		return
	}
	polyTermsOf(p, surf, ap.powers)
	var vs []coeffVar
	var br []bracket
	for _, pw := range ap.powers {
		vs = append(vs, polyVar(p, surf, pw))
		sc := 0.05 / math.Pow(mh, float64(pw))
		br = append(br, bracket{lo: -sc, hi: sc})
	}
	if _, err := solveCoeffs(p, vs, br, ap.passes); err != nil {
		failf("achromat asphere %s/%s f=%g F/%s s%d: %v", crown, flint, efl, trimNum(fnum, 2), surf, err)
		return
	}
	for i := range p.Cards[surf].PolyTerms {
		p.Cards[surf].PolyTerms[i].Coeff = sig6(p.Cards[surf].PolyTerms[i].Coeff)
	}
	if err := autoApertures(p, semi); err != nil {
		failf("achromat asphere %s/%s f=%g F/%s: %v", crown, flint, efl, trimNum(fnum, 2), err)
		return
	}
	if err := p.Validate(); err != nil {
		failf("achromat asphere %s/%s f=%g F/%s: %v", crown, flint, efl, trimNum(fnum, 2), err)
		return
	}
	r1, err := spotOnAxis(p)
	if err != nil {
		failf("achromat asphere %s/%s f=%g F/%s: %v", crown, flint, efl, trimNum(fnum, 2), err)
		return
	}
	if !(r1 > 0) || !(r1 <= ap.maxRatio*r0) {
		failf("achromat-asphere s%d f=%g %s/%s: aspheric RMS %.6gmm does not improve on the base %.6gmm (ratio %.4g > %.4g)",
			surf, efl, crown, flint, r1, r0, r1/r0, ap.maxRatio)
		return
	}
	file := fmt.Sprintf("achromat-asphere-s%d-%s-f%s-%s-%s.json", surf, ap.key, trimNum(efl, 0), glassCode(crown), glassCode(flint))
	p.Name = fmt.Sprintf("非球面修正双胶合 %s/%s f=%.4gmm F/%s（第 %d 面 %s）", crown, flint, efl, trimNum(fnum, 2), surf, ap.zh)
	note := fmt.Sprintf("双胶合消色差 + 非球面修正｜%s+%s｜f=%smm、F/%s｜基底为 DLS 自主设计消色差（同 achromat-dls 家族），在第 %d 面追加 %s 项（系数由引擎点列 RMS 坐标搜索拟合、六位有效数字冻结）｜基底实测轴上 RMS=%smm，修正后 RMS=%smm，降至 %.1f%%",
		crown, flint, trimNum(efl, 0), trimNum(fnum, 2), surf, termsZh(ap.powers), trimNum(r0, 4), trimNum(r1, 4), 100*r1/r0)
	l.add(file, p, "asphere", 3, []string{
		"asphere", "poly-terms", "achromat", "doublet", crown, flint, fnumTag(fnum), eflTag(efl),
		"visible", "optimized", "quality:optimized", "compare:sphere", termsTag(ap.powers),
	}, note)
}

// qconEquivalent registers the Forbes Qcon form of an asphere that was fitted
// as a polynomial r⁴ term: for m = 0 the Qcon departure is s0·u⁴ with
// u = r/ρ_max, so s0 = a4·ρ_max⁴ is the identical surface, and the measured
// RMS of the two forms must agree (the description states both).
func qconEquivalent(l *lib, glass string, fl, wantFnum float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	p, r0, r1, mh, err := solvedAsphereSinglet(glass, fl, wantFnum, []int{4}, 1)
	if err != nil {
		failf("qcon equivalent %s f=%g: %v", glass, fl, err)
		return
	}
	if !(r1 > 0) || !(r1 <= 0.5*r0) {
		failf("qcon-equivalent f=%g %s: the polynomial base does not improve (%.6g vs %.6g)", efl, glass, r1, r0)
		return
	}
	const card = 1
	rho := sig6(mh)
	a4 := p.Cards[card].PolyTerms[0].Coeff
	p.Cards[card].PolyTerms = nil
	p.Cards[card].QAsphere = &lens.QCard{Form: "qcon", NormRadiusMM: rho, Coeffs: []float64{sig6(a4 * rho * rho * rho * rho)}}
	if err := autoApertures(p, efl/(2*fnum)); err != nil {
		failf("qcon equivalent %s f=%g: %v", glass, fl, err)
		return
	}
	if err := p.Validate(); err != nil {
		failf("qcon equivalent %s f=%g: %v", glass, fl, err)
		return
	}
	rq, err := spotOnAxis(p)
	if err != nil {
		failf("qcon equivalent %s f=%g: %v", glass, fl, err)
		return
	}
	if rel := math.Abs(rq-r1) / r1; !(rel <= 1e-4) {
		failf("qcon-equivalent f=%g %s: qcon RMS %.9g differs from the polynomial %.9g (relative %.3g)", efl, glass, rq, r1, rel)
		return
	}
	file := fmt.Sprintf("qasphere-qcon-equiv-f%s-%s-%s.json", trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	p.Name = fmt.Sprintf("Forbes Qcon 非球面单透镜 %s f=%.4gmm F/%s（与多项式面等价）", glass, efl, trimNum(fnum, 2))
	note := fmt.Sprintf("Forbes Qcon 非球面（q_asphere.form=qcon）｜%s｜f=%smm、F/%s｜ρ_max=%smm｜s0=a4·ρ_max⁴（m=0 项 u⁴ 即 r⁴/ρ_max⁴，与多项式 a4 面逐点等价，故实测 RMS 与多项式面一致）｜球面基底 RMS=%smm，Qcon 修正后 RMS=%smm（多项式面 %.6gmm，相对差 %s）",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(rho, 4), trimNum(r0, 4), trimNum(rq, 4), r1, trimNum(math.Abs(rq-r1)/r1, 4))
	l.add(file, p, "qasphere", 3, []string{
		"qasphere", "qcon", "forbes", "singlet", glass, fnumTag(fnum), eflTag(efl), "visible",
		"optimized", "quality:optimized", "compare:poly", "sweep:efl",
	}, note)
}

// qconFit registers a Qcon asphere whose coefficients (s0, s1, …) are fitted
// directly on the engine: the base is the spherical singlet and the departure
// is optimised as a Forbes coefficient vector, not converted from a
// polynomial.
func qconFit(l *lib, glass string, fl, wantFnum float64, terms int, maxRatio float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	p, err := buildLayout(fmt.Sprintf("Forbes Qcon 非球面单透镜 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2)),
		[]grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		failf("qcon fit %s f=%g: %v", glass, fl, err)
		return
	}
	r0, err := spotOnAxis(p)
	if err != nil {
		failf("qcon fit %s f=%g: %v", glass, fl, err)
		return
	}
	const card = 1
	mh, err := marginalHeightAt(p, card)
	if err != nil {
		failf("qcon fit %s f=%g: %v", glass, fl, err)
		return
	}
	rho := sig6(mh)
	p.Cards[card].QAsphere = &lens.QCard{Form: "qcon", NormRadiusMM: rho, Coeffs: make([]float64, terms)}
	var vs []coeffVar
	var br []bracket
	for m := 0; m < terms; m++ {
		vs = append(vs, qconVar(p, card, m))
		sc := 0.03
		if m == 1 {
			sc = 0.02
		} else if m >= 2 {
			sc = 0.01
		}
		br = append(br, bracket{lo: -sc, hi: sc})
	}
	if _, err := solveCoeffs(p, vs, br, 2); err != nil {
		failf("qcon fit %s f=%g: %v", glass, fl, err)
		return
	}
	for i := range p.Cards[card].QAsphere.Coeffs {
		p.Cards[card].QAsphere.Coeffs[i] = sig6(p.Cards[card].QAsphere.Coeffs[i])
	}
	if err := autoApertures(p, semi); err != nil {
		failf("qcon fit %s f=%g: %v", glass, fl, err)
		return
	}
	if err := p.Validate(); err != nil {
		failf("qcon fit %s f=%g: %v", glass, fl, err)
		return
	}
	r1, err := spotOnAxis(p)
	if err != nil {
		failf("qcon fit %s f=%g: %v", glass, fl, err)
		return
	}
	if !(r1 > 0) || !(r1 <= maxRatio*r0) {
		failf("qcon-fit m=%d f=%g %s: RMS %.6gmm does not improve on the spherical %.6gmm (ratio %.4g > %.4g)",
			terms, efl, glass, r1, r0, r1/r0, maxRatio)
		return
	}
	file := fmt.Sprintf("qasphere-qcon-m%d-f%s-%s-%s.json", terms, trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	p.Name = fmt.Sprintf("Forbes Qcon 非球面单透镜 %s f=%.4gmm F/%s（%d 项直接拟合）", glass, efl, trimNum(fnum, 2), terms)
	note := fmt.Sprintf("Forbes Qcon 非球面（q_asphere.form=qcon，直接拟合）｜%s｜f=%smm、F/%s｜ρ_max=%smm、%d 个系数 s0…s%d｜球面基底实测轴上 RMS=%smm，Qcon 多项式修正后 RMS=%smm，降至 %.1f%%（系数由引擎点列 RMS 坐标搜索拟合）",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(rho, 4), terms, terms-1, trimNum(r0, 4), trimNum(r1, 4), 100*r1/r0)
	l.add(file, p, "qasphere", 3, []string{
		"qasphere", "qcon", "forbes", "singlet", glass, fnumTag(fnum), eflTag(efl), "visible",
		"optimized", "quality:optimized", "compare:sphere", fmt.Sprintf("coeffs:%d", terms),
	}, note)
}

// qbfsFit registers a Forbes Qbfs asphere (best-fit sphere base, departure
// along the normal) whose coefficients are fitted on the engine. The base
// radius is re-derived on every coefficient change so that the vertex
// curvature — and with it the first-order design — stays exactly at the
// spherical base value (see qbfsVar).
func qbfsFit(l *lib, glass string, fl, wantFnum float64, terms int, maxRatio float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	p, err := buildLayout(fmt.Sprintf("Forbes Qbfs 非球面单透镜 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2)),
		[]grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		failf("qbfs fit %s f=%g: %v", glass, fl, err)
		return
	}
	r0, err := spotOnAxis(p)
	if err != nil {
		failf("qbfs fit %s f=%g: %v", glass, fl, err)
		return
	}
	const card = 1
	mh, err := marginalHeightAt(p, card)
	if err != nil {
		failf("qbfs fit %s f=%g: %v", glass, fl, err)
		return
	}
	rho := sig6(mh)
	c0 := 1 / p.Cards[card].RadiusMM
	p.Cards[card].QAsphere = &lens.QCard{Form: "qbfs", NormRadiusMM: rho, Coeffs: make([]float64, terms)}
	var vs []coeffVar
	var br []bracket
	for m := 0; m < terms; m++ {
		vs = append(vs, qbfsVar(p, card, m, c0))
		sc := 0.02
		if m == 1 {
			sc = 0.01
		} else if m >= 2 {
			sc = 0.005
		}
		br = append(br, bracket{lo: -sc, hi: sc})
	}
	if _, err := solveCoeffs(p, vs, br, 2); err != nil {
		failf("qbfs fit %s f=%g: %v", glass, fl, err)
		return
	}
	for i := range p.Cards[card].QAsphere.Coeffs {
		p.Cards[card].QAsphere.Coeffs[i] = sig6(p.Cards[card].QAsphere.Coeffs[i])
	}
	p.Cards[card].RadiusMM = qbfsRadius(p.Cards[card].QAsphere, c0)
	if err := autoApertures(p, semi); err != nil {
		failf("qbfs fit %s f=%g: %v", glass, fl, err)
		return
	}
	if err := p.Validate(); err != nil {
		failf("qbfs fit %s f=%g: %v", glass, fl, err)
		return
	}
	r1, err := spotOnAxis(p)
	if err != nil {
		failf("qbfs fit %s f=%g: %v", glass, fl, err)
		return
	}
	if !(r1 > 0) || !(r1 <= maxRatio*r0) {
		failf("qbfs-fit m=%d f=%g %s: RMS %.6gmm does not improve on the spherical %.6gmm (ratio %.4g > %.4g)",
			terms, efl, glass, r1, r0, r1/r0, maxRatio)
		return
	}
	file := fmt.Sprintf("qasphere-qbfs-m%d-f%s-%s-%s.json", terms, trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	p.Name = fmt.Sprintf("Forbes Qbfs 非球面单透镜 %s f=%.4gmm F/%s（%d 项直接拟合）", glass, efl, trimNum(fnum, 2), terms)
	note := fmt.Sprintf("Forbes Qbfs 非球面（q_asphere.form=qbfs，最佳拟合球基底）｜%s｜f=%smm、F/%s｜ρ_max=%smm、%d 个系数 a0…a%d｜基底半径随系数重解以保持顶点曲率（=一阶焦距不变；Qbfs 面 u²(1−u²) 项含 r² 分量，曲率 c+2S(0)/ρ_max² 由 geo.QAsphere.Curvature 定义）｜球面基底实测轴上 RMS=%smm，Qbfs 修正后 RMS=%smm，降至 %.1f%%",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(rho, 4), terms, terms-1, trimNum(r0, 4), trimNum(r1, 4), 100*r1/r0)
	l.add(file, p, "qasphere", 3, []string{
		"qasphere", "qbfs", "forbes", "singlet", glass, fnumTag(fnum), eflTag(efl), "visible",
		"optimized", "quality:optimized", "compare:sphere", fmt.Sprintf("coeffs:%d", terms),
	}, note)
}

// qconZero registers the degenerate zero-coefficient Qcon: a q_asphere with
// the conic base (R, K) and no departure at all (an empty coefficient list).
// The description states the equivalence with the bare conic and the check
// compares the two systems' measured metrics in the generator.
func qconZero(l *lib, glass string, fl, wantFnum, K float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	p, err := buildLayout(fmt.Sprintf("Qcon 零系数退化 %s f=%.4gmm F/%s K=%s", glass, efl, trimNum(fnum, 2), trimNum(K, 2)),
		[]grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		failf("qcon zero %s f=%g K=%g: %v", glass, fl, K, err)
		return
	}
	const card = 1
	mh, err := marginalHeightAt(p, card)
	if err != nil {
		failf("qcon zero %s f=%g K=%g: %v", glass, fl, K, err)
		return
	}
	ref := p.Clone()
	ref.Cards[card].ConicK = K
	p.Cards[card].ConicK = K
	p.Cards[card].QAsphere = &lens.QCard{Form: "qcon", NormRadiusMM: sig6(mh)}
	if err := autoApertures(p, semi); err != nil {
		failf("qcon zero %s f=%g K=%g: %v", glass, fl, K, err)
		return
	}
	if err := autoApertures(ref, semi); err != nil {
		failf("qcon zero ref %s f=%g K=%g: %v", glass, fl, K, err)
		return
	}
	m0, err := measureOf(ref)
	if err != nil {
		failf("qcon zero ref %s f=%g K=%g: %v", glass, fl, K, err)
		return
	}
	m1, err := measureOf(p)
	if err != nil {
		failf("qcon zero %s f=%g K=%g: %v", glass, fl, K, err)
		return
	}
	if math.Abs(m1.RMSAxisMM-m0.RMSAxisMM) > 1e-9*m0.RMSAxisMM {
		failf("qcon-zero f=%g %s K=%g: the zero-coefficient Qcon does not reproduce the conic (RMS %.9g vs %.9g)",
			efl, glass, K, m1.RMSAxisMM, m0.RMSAxisMM)
		return
	}
	file := fmt.Sprintf("qasphere-qcon-zero-f%s-%s-k%s-%s.json", trimNum(efl, 2), fnumCode(fnum), numCode(K), glassCode(glass))
	p.Name = fmt.Sprintf("Forbes Qcon 零系数退化单透镜 %s f=%.4gmm F/%s K=%s", glass, efl, trimNum(fnum, 2), trimNum(K, 2))
	note := fmt.Sprintf("Forbes Qcon 零系数退化例｜%s｜f=%smm、F/%s、基底 K=%s｜q_asphere.form=qcon 且系数表为空（或全零）：u⁴ 级数项全部为零，面型严格退化为基底圆锥面 z=conic(R,K)，与非球面描述缺省时的 conic_k 系统逐点相同（实测轴上 RMS=%smm，与直接写 conic_k=%s 的同一系统一致，由生成器比对）",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(K, 2), trimNum(m1.RMSAxisMM, 4), trimNum(K, 2))
	l.add(file, p, "qasphere", 2, []string{
		"qasphere", "qcon", "forbes", "degenerate", "zero-coefficients", "singlet", glass,
		fnumTag(fnum), eflTag(efl), "visible", "analytic", "quality:analytic", "K:" + trimNum(K, 2),
	}, note)
}

// qbfsZero registers the degenerate zero-coefficient Qbfs: the departure
// coefficients vanish, so the surface is exactly the best-fit sphere (K = 0)
// and the first-order design is unchanged.
func qbfsZero(l *lib, glass string, fl, wantFnum float64) {
	efl := math.Abs(fl)
	fnum := fnumFor(fl, wantFnum)
	semi := efl / (2 * fnum)
	g := grp{name: "L1", curv: seedConvPCX(glass, fl), glasses: []string{glass, ""}, thick: thickFor(efl), after: 0.2*efl + 10}
	p, err := buildLayout(fmt.Sprintf("Qbfs 零系数退化 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2)),
		[]grp{g}, 2.5, fl, semi, math.Inf(-1))
	if err != nil {
		failf("qbfs zero %s f=%g: %v", glass, fl, err)
		return
	}
	const card = 1
	mh, err := marginalHeightAt(p, card)
	if err != nil {
		failf("qbfs zero %s f=%g: %v", glass, fl, err)
		return
	}
	ref := p.Clone()
	p.Cards[card].QAsphere = &lens.QCard{Form: "qbfs", NormRadiusMM: sig6(mh), Coeffs: []float64{0, 0}}
	if err := autoApertures(p, semi); err != nil {
		failf("qbfs zero %s f=%g: %v", glass, fl, err)
		return
	}
	if err := autoApertures(ref, semi); err != nil {
		failf("qbfs zero ref %s f=%g: %v", glass, fl, err)
		return
	}
	m0, err := measureOf(ref)
	if err != nil {
		failf("qbfs zero ref %s f=%g: %v", glass, fl, err)
		return
	}
	m1, err := measureOf(p)
	if err != nil {
		failf("qbfs zero %s f=%g: %v", glass, fl, err)
		return
	}
	if math.Abs(m1.RMSAxisMM-m0.RMSAxisMM) > 1e-9*m0.RMSAxisMM {
		failf("qbfs-zero f=%g %s: the zero-coefficient Qbfs does not reproduce the sphere (RMS %.9g vs %.9g)",
			efl, glass, m1.RMSAxisMM, m0.RMSAxisMM)
		return
	}
	file := fmt.Sprintf("qasphere-qbfs-zero-f%s-%s-%s.json", trimNum(efl, 2), fnumCode(fnum), glassCode(glass))
	p.Name = fmt.Sprintf("Forbes Qbfs 零系数退化单透镜 %s f=%.4gmm F/%s", glass, efl, trimNum(fnum, 2))
	note := fmt.Sprintf("Forbes Qbfs 零系数退化例｜%s｜f=%smm、F/%s｜q_asphere.form=qbfs 且系数表为全零：法向级数项全部为零，面型退化为最佳拟合球（K=0 的球面），与普通球面系统逐点相同（实测轴上 RMS=%smm，与无 q_asphere 描述的同一系统一致，由生成器比对）",
		glass, trimNum(efl, 0), trimNum(fnum, 2), trimNum(m1.RMSAxisMM, 4))
	l.add(file, p, "qasphere", 2, []string{
		"qasphere", "qbfs", "forbes", "degenerate", "zero-coefficients", "singlet", glass,
		fnumTag(fnum), eflTag(efl), "visible", "analytic", "quality:analytic",
	}, note)
}

// measureOf runs the library measurement of one prescription (the same code
// path libcheck applies), used for the generator's own equivalence checks.
func measureOf(p *lens.Prescription) (*measure.Metrics, error) {
	return measure.Measure(p)
}

// ---------- the families of this catalogue ----------

// The polynomial families: powers, coordinate-search passes and the largest
// asphere/sphere RMS ratio that still counts as a real improvement.
var (
	famA4   = asphereFamily{key: "a4", zh: "单项 r⁴", powers: []int{4}, passes: 1, maxRatio: 0.5}
	famA46  = asphereFamily{key: "a46", zh: "r⁴+r⁶", powers: []int{4, 6}, passes: 2, maxRatio: 0.35}
	famA468 = asphereFamily{key: "a468", zh: "r⁴+r⁶+r⁸", powers: []int{4, 6, 8}, passes: 2, maxRatio: 0.35}
	famOdd  = asphereFamily{key: "o35", zh: "奇次 r³+r⁵", powers: []int{3, 5}, passes: 2, maxRatio: 0.85}
)

func buildAspheres(l *lib) {
	// single r⁴ correction: focal length sweep, F-number sweep, glasses
	for _, f := range []float64{25, 50, 100, 200, 400} {
		asphereSinglet(l, "N-BK7", f, 4, famA4)
	}
	for _, fn := range []float64{2.8, 3.5, 5.6, 8, 11} {
		asphereSinglet(l, "N-BK7", 100, fn, famA4)
	}
	// the fast singlets need higher-order terms: a single r⁴ cannot hold the
	// spherical aberration of F/2…F/2.8
	for _, fn := range []float64{2, 2.4} {
		asphereSinglet(l, "N-BK7", 100, fn, famA468)
	}
	asphereSinglet(l, "N-BK7", 100, 2.8, famA46)
	for _, c := range []struct {
		glass string
		f, fn float64
	}{{"FusedSilica", 50, 2.8}, {"FusedSilica", 50, 4}, {"FusedSilica", 200, 5.6},
		{"N-SK16", 50, 2.8}, {"N-SK16", 50, 4}, {"N-SK16", 200, 5.6},
		{"F2", 100, 4}, {"F2", 200, 4}} {
		asphereSinglet(l, c.glass, c.f, c.fn, famA4)
	}
	// two- and three-term corrections
	for _, c := range []struct {
		f, fn float64
	}{{50, 4}, {100, 4}, {100, 5.6}, {200, 4}, {400, 4}, {50, 2.8}, {200, 2.8}, {100, 8}} {
		asphereSinglet(l, "N-BK7", c.f, c.fn, famA46)
	}
	for _, c := range []struct {
		f, fn float64
	}{{50, 4}, {100, 4}, {200, 4}, {400, 4}, {50, 2}, {150, 2.8}, {200, 2.8}, {400, 5.6}} {
		asphereSinglet(l, "N-BK7", c.f, c.fn, famA468)
	}
	// odd powers (r³, r⁵)
	for _, c := range []struct {
		glass string
		f, fn float64
	}{{"N-BK7", 50, 4}, {"N-BK7", 100, 4}, {"N-BK7", 100, 5.6}, {"N-BK7", 200, 4},
		{"N-BK7", 400, 4}, {"FusedSilica", 200, 5.6}, {"N-BK7", 75, 3.5}} {
		asphereSinglet(l, c.glass, c.f, c.fn, famOdd)
	}
	// cemented achromat + one polynomial correction
	for _, pr := range [][2]string{{"N-BK7", "F2"}, {"N-SK16", "N-SF5"}} {
		for _, f := range []float64{50, 100, 150, 200, 300, 400} {
			achromatAsphere(l, pr[0], pr[1], f, 5, 2, famA4)
		}
		achromatAsphere(l, pr[0], pr[1], 100, 5, 0, famA4)
		achromatAsphere(l, pr[0], pr[1], 100, 5, 1, famA4)
		achromatAsphere(l, pr[0], pr[1], 100, 5, 2, famA46)
		achromatAsphere(l, pr[0], pr[1], 200, 5, 2, famA46)
		achromatAsphere(l, pr[0], pr[1], 100, 2.8, 2, famA4)
	}
	// Forbes Qcon: exact conversion of the r⁴ fit, direct fits, degenerate zeros
	for _, c := range []struct {
		glass string
		f, fn float64
	}{{"N-BK7", 50, 4}, {"N-BK7", 100, 2.8}, {"N-BK7", 100, 4}, {"N-BK7", 100, 5.6},
		{"N-BK7", 200, 4}, {"N-BK7", 400, 4}, {"N-BK7", 300, 4}, {"FusedSilica", 100, 4}} {
		qconEquivalent(l, c.glass, c.f, c.fn)
	}
	for _, c := range []struct {
		glass string
		f, fn float64
	}{{"N-BK7", 50, 4}, {"N-BK7", 100, 2.8}, {"N-BK7", 100, 4}, {"N-BK7", 100, 5.6},
		{"N-BK7", 200, 4}, {"N-BK7", 400, 4}, {"N-BK7", 300, 4}, {"N-BK7", 500, 4}} {
		qconFit(l, c.glass, c.f, c.fn, 2, 0.35)
	}
	for _, f := range []float64{100, 200, 400} {
		qconFit(l, "N-BK7", f, 4, 3, 0.35)
	}
	for _, c := range []struct {
		glass string
		f, K  float64
	}{{"N-BK7", 100, 0}, {"N-BK7", 100, -0.5}, {"N-BK7", 100, -1}, {"N-BK7", 100, -2},
		{"N-BK7", 200, -1}, {"FusedSilica", 100, -1}, {"N-BK7", 400, -1.5}} {
		qconZero(l, c.glass, c.f, 4, c.K)
	}
	// Forbes Qbfs: direct fits and degenerate zeros
	for _, c := range []struct {
		glass string
		f, fn float64
	}{{"N-BK7", 50, 4}, {"N-BK7", 100, 2.8}, {"N-BK7", 100, 4}, {"N-BK7", 100, 5.6},
		{"N-BK7", 200, 4}, {"N-BK7", 400, 4}, {"N-BK7", 300, 4}, {"N-BK7", 500, 4}} {
		qbfsFit(l, c.glass, c.f, c.fn, 1, 0.7)
	}
	for _, c := range []struct {
		f, fn float64
	}{{100, 2.8}, {100, 4}, {200, 4}, {400, 4}, {100, 5.6}} {
		qbfsFit(l, "N-BK7", c.f, c.fn, 2, 0.45)
	}
	qbfsFit(l, "N-BK7", 100, 4, 3, 0.45)
	qbfsFit(l, "FusedSilica", 100, 4, 2, 0.45)
	qbfsZero(l, "N-BK7", 100, 4)
	qbfsZero(l, "N-BK7", 200, 4)
	qbfsZero(l, "FusedSilica", 100, 4)
}
