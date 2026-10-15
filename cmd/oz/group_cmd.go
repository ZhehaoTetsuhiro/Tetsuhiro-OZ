package main

// `oz lens group ...`: the editable element-level lens group (lens.Group) on
// the command line, and `oz lens parax`: the first-order report of a
// prescription with mirror support (geo.Paraxial).
//
// Every subcommand calls the shared core in group_api.go, so `oz lens group
// image --json out.json` writes the response body of POST
// /api/lens/group/image and the flag names match the request fields.

import (
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// multiFlag collects a flag that may be repeated (-p a -p b).
type multiFlag []string

func (m *multiFlag) String() string     { return strings.Join(*m, ",") }
func (m *multiFlag) Set(v string) error { *m = append(*m, v); return nil }

// flagSet reports whether a flag was given on the command line, so that "not
// given" is distinguishable from the zero value.
func flagSet(fs *flag.FlagSet, name string) bool {
	found := false
	fs.Visit(func(f *flag.Flag) {
		if f.Name == name {
			found = true
		}
	})
	return found
}

// groupBundleFlags are the bundle flags shared by image / afocal / dispersion:
// the object, the field, the pupil and the sampling.
type groupBundleFlags struct {
	objectZ     string
	fieldDeg    float64
	fieldTan    float64
	pupil       float64
	pupilZ      float64
	rings       int
	lambda      float64
	hasFieldTan bool
}

// register declares the bundle flags on a flag set.
func (f *groupBundleFlags) register(fs *flag.FlagSet) {
	fs.StringVar(&f.objectZ, "object-z", "", "物面 z mm（\"inf\" = 无穷远；缺省用 group 自己的 object_z_mm）")
	fs.Float64Var(&f.fieldDeg, "field-deg", 0, "无穷远物的视场角（度，子午 y 方向）")
	fs.Float64Var(&f.fieldTan, "field-tan", 0, "无穷远物=场角切线；有限物=物高 mm（给定时优先于 --field-deg）")
	fs.Float64Var(&f.pupil, "pupil", 0, "采样光瞳半径 mm（0 = 自动：光阑元件的口径）")
	fs.Float64Var(&f.pupilZ, "pupil-z", 0, "显式光瞳平面的 z mm（配合 --pupil 使用）")
	fs.IntVar(&f.rings, "rings", 3, "Hexapolar 采样环数（1–8）")
	fs.Float64Var(&f.lambda, "lambda", dLineNm, "波长 nm")
}

// finish records what only the parse can tell (which field flag was given).
func (f *groupBundleFlags) finish(fs *flag.FlagSet) { f.hasFieldTan = flagSet(fs, "field-tan") }

// apply fills the bundle fields of a request.
func (f *groupBundleFlags) apply(req *groupRequest) error {
	raw, err := objectZFlag(f.objectZ)
	if err != nil {
		return err
	}
	req.ObjectZMM = raw
	if f.hasFieldTan {
		req.FieldTan = &f.fieldTan
	} else {
		req.FieldDeg = f.fieldDeg
	}
	req.PupilRadiusMM = f.pupil
	req.PupilZMM = f.pupilZ
	req.Rings = f.rings
	req.LambdaNm = f.lambda
	return nil
}

// cmdLensGroup dispatches the group subcommands.
func cmdLensGroup(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("lens group: 需要子命令 show | layout | image | afocal | dispersion | set | to-prescription")
	}
	sub, rest := args[0], args[1:]
	switch sub {
	case "show":
		return cmdGroupShow(rest)
	case "layout":
		return cmdGroupLayout(rest)
	case "image":
		return cmdGroupImage(rest)
	case "afocal":
		return cmdGroupAfocal(rest)
	case "dispersion":
		return cmdGroupDispersion(rest)
	case "set":
		return cmdGroupSet(rest)
	case "to-prescription":
		return cmdGroupToPrescription(rest)
	}
	return fmt.Errorf("lens group: 未知子命令 %q（可用: show, layout, image, afocal, dispersion, set, to-prescription）", sub)
}

// groupFileArg parses the flags of a group subcommand and returns the path of
// the single <group.json> argument (it may sit anywhere among the flags).
func groupFileArg(sub string, fs *flag.FlagSet, args []string) (string, error) {
	pos, err := parseFlagsAll(fs, args)
	if err != nil {
		return "", err
	}
	if len(pos) != 1 {
		return "", fmt.Errorf("lens group %s: 需要一个 <group.json> 参数", sub)
	}
	return pos[0], nil
}

// loadGroupFile reads and validates a group file.
func loadGroupFile(path string) (*lens.Group, error) {
	g, err := lens.LoadGroupJSON(path)
	if err != nil {
		return nil, fmt.Errorf("lens group: %w", err)
	}
	return g, nil
}

// objectZFlag turns the --object-z value into a request field: "" leaves the
// group's (or the prescription's) own object, "inf" is an object at infinity,
// anything else is a finite object plane.
func objectZFlag(v string) (raw []byte, err error) {
	s := strings.TrimSpace(v)
	switch strings.ToLower(s) {
	case "":
		return nil, nil
	case "inf", "-inf", "+inf", "infinity":
		return []byte(`"inf"`), nil
	}
	z, err := strconv.ParseFloat(s, 64)
	if err != nil || math.IsNaN(z) || math.IsInf(z, 0) {
		return nil, fmt.Errorf("--object-z %q 不是有限数字或 \"inf\"", v)
	}
	return []byte(strconv.FormatFloat(z, 'g', -1, 64)), nil
}

// cmdGroupShow prints the one-line-per-element summary of a group.
func cmdGroupShow(args []string) error {
	fs := flag.NewFlagSet("lens group show", flag.ExitOnError)
	path, err := groupFileArg("show", fs, args)
	if err != nil {
		return err
	}
	g, err := loadGroupFile(path)
	if err != nil {
		return err
	}
	fmt.Println(g.String())
	return nil
}

// cmdGroupLayout prints the surface table of the layout (the global vertex,
// axis, normal and fold count of every surface) and the exit frame.
func cmdGroupLayout(args []string) error {
	fs := flag.NewFlagSet("lens group layout", flag.ExitOnError)
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	path, err := groupFileArg("layout", fs, args)
	if err != nil {
		return err
	}
	g, err := loadGroupFile(path)
	if err != nil {
		return err
	}
	res, err := runGroupLayout(&groupRequest{Group: g})
	if err != nil {
		return err
	}
	printGroupLayout(path, res)
	return writeJSONFile(*jsonOut, res)
}

// printGroupLayout prints the layout table.
func printGroupLayout(path string, res *groupLayoutResponse) {
	fmt.Printf("布局 %s: %d 面（面号 = Group.System() 里的下标，含坐标断点），折转 %d 次\n",
		path, res.NumSurfaces, res.Folds)
	fmt.Printf("  %4s %-5s %6s  %-30s %-24s %-24s %4s\n",
		"元件", "面", "面号", "顶点 (x, y, z) mm", "轴向", "法向", "折转")
	for _, s := range res.Surfaces {
		fmt.Printf("  %4d %-5s %6d  (%+10.6g, %+10.6g, %+10.6g)  (%+.4g, %+.4g, %+.4g)  (%+.4g, %+.4g, %+.4g)  %4d\n",
			s.Element, s.Face, s.Surface,
			s.Vertex.X, s.Vertex.Y, s.Vertex.Z, s.Axis.X, s.Axis.Y, s.Axis.Z, s.Normal.X, s.Normal.Y, s.Normal.Z, s.Folds)
	}
	fmt.Printf("  出射：顶点 (%+.6g, %+.6g, %+.6g) mm，光束方向 (%+.6g, %+.6g, %+.6g)\n",
		res.ExitVertex.X, res.ExitVertex.Y, res.ExitVertex.Z, res.ExitAxis.X, res.ExitAxis.Y, res.ExitAxis.Z)
}

// cmdGroupImage prints the best-focus image of a bundle through the group.
func cmdGroupImage(args []string) error {
	fs := flag.NewFlagSet("lens group image", flag.ExitOnError)
	var bf groupBundleFlags
	bf.register(fs)
	scan := fs.Int("scan", 0, "最佳焦面搜索的粗采样点数（0 = 65）")
	rangeMM := fs.Float64("range", 0, "焦面搜索窗半宽 mm（0 = 自动）")
	tolMM := fs.Float64("tol", 0, "焦面搜索的轴向容差 mm（0 = 1e-9）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	path, err := groupFileArg("image", fs, args)
	if err != nil {
		return err
	}
	bf.finish(fs)
	g, err := loadGroupFile(path)
	if err != nil {
		return err
	}
	req := &groupRequest{Group: g, Scan: *scan, RangeMM: *rangeMM, TolMM: *tolMM}
	if err := bf.apply(req); err != nil {
		return err
	}
	res, err := runGroupImage(req)
	if err != nil {
		return err
	}
	printGroupImage(path, res)
	return writeJSONFile(*jsonOut, res)
}

// printGroupImage prints one image search.
func printGroupImage(path string, res *groupImageResponse) {
	obj := "∞"
	if res.ObjectZMM != nil {
		obj = fmt.Sprintf("%.6g mm", *res.ObjectZMM)
	}
	fmt.Printf("成像位置 %s: λ=%.4f nm，物 %s，场值 %v，%d 条光线（%d 条被拦）\n",
		path, res.LambdaNm, obj, res.FieldTanY, res.Rays+res.Failed, res.Failed)
	fmt.Printf("  像点 (%+.9g, %+.9g, %+.9g) mm；RMS %.6g mm，最大半径 %.6g mm\n",
		res.Point.X, res.Point.Y, res.Point.Z, res.RMSMM, res.GeoRadiusMM)
	fmt.Printf("  焦面：过 (%+.6g, %+.6g, %+.6g) 法向 (%+.6g, %+.6g, %+.6g)，焦面距离 %+.9g mm\n",
		res.PlaneO.X, res.PlaneO.Y, res.PlaneO.Z, res.Axis.X, res.Axis.Y, res.Axis.Z, res.FocusMM)
	fmt.Printf("  主光线落点 (%+.6g, %+.6g, %+.6g)，方向 (%+.6g, %+.6g, %+.6g)\n",
		res.ChiefHit.X, res.ChiefHit.Y, res.ChiefHit.Z, res.ChiefDir.X, res.ChiefDir.Y, res.ChiefDir.Z)
}

// cmdGroupAfocal prints the angular magnification of an afocal group.
func cmdGroupAfocal(args []string) error {
	fs := flag.NewFlagSet("lens group afocal", flag.ExitOnError)
	var bf groupBundleFlags
	bf.register(fs)
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	path, err := groupFileArg("afocal", fs, args)
	if err != nil {
		return err
	}
	bf.finish(fs)
	g, err := loadGroupFile(path)
	if err != nil {
		return err
	}
	req := &groupRequest{Group: g}
	if err := bf.apply(req); err != nil {
		return err
	}
	res, err := runGroupAfocal(req)
	if err != nil {
		return err
	}
	printGroupAfocal(path, res)
	return writeJSONFile(*jsonOut, res)
}

// printGroupAfocal prints the angular magnification.
func printGroupAfocal(path string, res *groupAfocalResponse) {
	fmt.Printf("角放大率 %s: λ=%.4f nm，%d 条光线，采样 %d 环\n", path, res.LambdaNm, res.Rays, res.Samples)
	fmt.Printf("  Γ = %+.12g；入/出射束直径比 %.9g（|Γ| = %.9g）；出射准直残差 %.3g rad\n",
		res.Gamma, res.BeamRatio, math.Abs(res.Gamma), res.SpreadRad)
	fmt.Printf("  入射轴 (%+.4g, %+.4g, %+.4g) → 出射轴 (%+.4g, %+.4g, %+.4g)\n",
		res.AxisIn.X, res.AxisIn.Y, res.AxisIn.Z, res.AxisOut.X, res.AxisOut.Y, res.AxisOut.Z)
	fmt.Printf("  探测束：入射 (%+.6g, %+.6g, %+.6g) → 出射 (%+.6g, %+.6g, %+.6g)\n",
		res.EntryDir.X, res.EntryDir.Y, res.EntryDir.Z, res.ExitDir.X, res.ExitDir.Y, res.ExitDir.Z)
}

// cmdGroupDispersion prints the axial focus shift between two wavelengths.
func cmdGroupDispersion(args []string) error {
	fs := flag.NewFlagSet("lens group dispersion", flag.ExitOnError)
	var bf groupBundleFlags
	bf.register(fs)
	lambdaA := fs.Float64("lambda-a", geo.WavelengthF, "波长 A（蓝，默认 F 线 486.1327 nm）")
	lambdaB := fs.Float64("lambda-b", geo.WavelengthC, "波长 B（红，默认 C 线 656.2725 nm）")
	scan := fs.Int("scan", 0, "最佳焦面搜索的粗采样点数（0 = 65）")
	rangeMM := fs.Float64("range", 0, "焦面搜索窗半宽 mm（0 = 自动）")
	tolMM := fs.Float64("tol", 0, "焦面搜索的轴向容差 mm（0 = 1e-9）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	path, err := groupFileArg("dispersion", fs, args)
	if err != nil {
		return err
	}
	bf.finish(fs)
	g, err := loadGroupFile(path)
	if err != nil {
		return err
	}
	req := &groupRequest{Group: g, Scan: *scan, RangeMM: *rangeMM, TolMM: *tolMM, LambdaANm: *lambdaA, LambdaBNm: *lambdaB}
	if err := bf.apply(req); err != nil {
		return err
	}
	res, err := runGroupDispersion(req)
	if err != nil {
		return err
	}
	printGroupDispersion(path, res)
	return writeJSONFile(*jsonOut, res)
}

// printGroupDispersion prints the F–C style focus shift.
func printGroupDispersion(path string, res *groupDispersionResponse) {
	fmt.Printf("色散 %s: λa=%.4f nm → 焦面 %+.9g mm；λb=%.4f nm → 焦面 %+.9g mm（沿出射轴，从末面顶点量起）\n",
		path, res.LambdaANm, res.A.FocusMM, res.LambdaBNm, res.B.FocusMM)
	fmt.Printf("  焦移 (a − b) = %+.9g mm\n", res.ShiftMM)
}

// cmdGroupSet applies the -p edits to a group and writes the result.
func cmdGroupSet(args []string) error {
	fs := flag.NewFlagSet("lens group set", flag.ExitOnError)
	var params multiFlag
	fs.Var(&params, "p", "参数编辑，形如 e0.front.R=120、e2.tilt_x_deg=45（可重复，按顺序应用）")
	out := fs.String("o", "", "把更新后的 group 写成 JSON 到该路径")
	path, err := groupFileArg("set", fs, args)
	if err != nil {
		return err
	}
	if len(params) == 0 {
		return fmt.Errorf("lens group set: 至少给一个 -p 参数")
	}
	g, err := loadGroupFile(path)
	if err != nil {
		return err
	}
	sets := make([]groupSetItem, 0, len(params))
	for _, p := range params {
		name, value, ok := strings.Cut(p, "=")
		if !ok {
			return fmt.Errorf("lens group set: -p %q 缺少 \"=\"（形如 e0.front.R=120）", p)
		}
		v, err := strconv.ParseFloat(strings.TrimSpace(value), 64)
		if err != nil {
			return fmt.Errorf("lens group set: -p %q 的值 %q 不是数字", p, value)
		}
		sets = append(sets, groupSetItem{Param: strings.TrimSpace(name), Value: v})
	}
	res, err := runGroupSet(&groupSetRequest{Group: g, Sets: sets})
	if err != nil {
		return err // lens 的参数错误原样回显
	}
	for _, s := range sets {
		if v, err := res.Group.GetParam(mustParam(s.Param)); err == nil {
			fmt.Printf("  %s = %+.10g\n", s.Param, v)
		}
	}
	fmt.Println(res.Group.String())
	if *out != "" {
		if err := res.Group.SaveJSON(*out); err != nil {
			return err
		}
		fmt.Fprintln(os.Stderr, "已写入", *out)
	}
	return nil
}

// mustParam parses a parameter that runGroupSet already accepted.
func mustParam(s string) lens.GroupParam {
	p, err := parseGroupParam(s)
	if err != nil {
		return lens.GroupParam{}
	}
	return p
}

// cmdGroupToPrescription converts a group into the card form.
func cmdGroupToPrescription(args []string) error {
	fs := flag.NewFlagSet("lens group to-prescription", flag.ExitOnError)
	out := fs.String("o", "", "处方 JSON 输出路径（必填）")
	path, err := groupFileArg("to-prescription", fs, args)
	if err != nil {
		return err
	}
	if *out == "" {
		return fmt.Errorf("lens group to-prescription: 需要 -o 指定处方 JSON 输出路径")
	}
	g, err := loadGroupFile(path)
	if err != nil {
		return err
	}
	p, err := runGroupPrescription(&groupRequest{Group: g})
	if err != nil {
		return err
	}
	if err := p.SaveJSON(*out); err != nil {
		return err
	}
	fmt.Printf("处方 ← %s: %d 张卡片，停位 %d，物距 %s\n", path, len(p.Cards), p.Stop, describeObject(p.ObjectZMM))
	fmt.Fprintln(os.Stderr, "已写入", *out)
	return nil
}

// describeObject formats an object position.
func describeObject(z float64) string {
	if math.IsInf(z, -1) {
		return "∞"
	}
	return fmt.Sprintf("%.6g mm", z)
}

// cmdLensParax prints the first-order report of a prescription.
func cmdLensParax(args []string) error {
	fs := flag.NewFlagSet("lens parax", flag.ExitOnError)
	objectZ := fs.String("object-z", "", "物面 z mm（缺省用处方自己的 object_z_mm；\"inf\" = 无穷远）")
	fieldTan := fs.Float64("field-tan", 0, "参考光线的场切线（afocal 的 Γ 需要 ≠0）")
	lambda := fs.Float64("lambda", dLineNm, "波长 nm")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	pos, err := parseFlagsAll(fs, args)
	if err != nil {
		return err
	}
	if len(pos) != 1 {
		return fmt.Errorf("lens parax: 需要一个 <处方.json> 参数")
	}
	rxPath := pos[0]
	p, err := lens.LoadJSON(rxPath)
	if err != nil {
		return err
	}
	raw, err := objectZFlag(*objectZ)
	if err != nil {
		return err
	}
	req := &paraxRequest{Prescription: *p, ObjectZMM: raw, LambdaNm: *lambda}
	if flagSet(fs, "field-tan") {
		req.FieldTan = fieldTan
	}
	res, err := runParax(req)
	if err != nil {
		return err
	}
	printParax(rxPath, res)
	return writeJSONFile(*jsonOut, res)
}

// printParax prints the first-order report.
func printParax(path string, res *paraxResponse) {
	obj := "∞"
	if res.ObjectZMM != nil {
		obj = fmt.Sprintf("z=%.6g mm", *res.ObjectZMM)
	}
	fmt.Printf("一阶 %s: λ=%.4f nm，物 %s\n", path, res.LambdaNm, obj)
	fmt.Printf("  EFL %s mm（习惯焦距 Parity·EFL %s）；FFL %s；BFL %s\n",
		f2s(res.EFLMM), f2s(res.FocalLenMM), f2s(res.FFLMM), f2s(res.BFLMM))
	fmt.Printf("  像面 z %s mm；横向放大率 m %s\n", f2s(res.ImageZMM), f2s(res.Magnifi))
	fmt.Printf("  入瞳 z %s mm r %s mm；出瞳 z %s mm r %s mm\n",
		f2s(res.EntrPupilZ), f2s(res.EntrPupilR), f2s(res.ExitPupilZ), f2s(res.ExitPupilR))
	af := "否"
	if res.Afocal {
		af = "是"
	}
	fmt.Printf("  镜面 %d 面；奇偶性 Parity %+.0f；无焦（afocal）%s\n", res.Mirrors, res.Parity, af)
	if res.Gamma != nil {
		fmt.Printf("  角放大率 Γ = %+.12g\n", *res.Gamma)
	}
	if res.GammaNote != "" {
		fmt.Printf("  Γ：%s\n", res.GammaNote)
	}
	if res.FCShiftMM != nil {
		fmt.Printf("  F–C 焦移 = %+.9g mm\n", *res.FCShiftMM)
	}
	if res.FCNote != "" {
		fmt.Printf("  F–C：%s\n", res.FCNote)
	}
}

// f2s formats an optional number ("—" for a quantity that does not exist).
func f2s(v *float64) string {
	if v == nil {
		return "—"
	}
	return fmt.Sprintf("%+.9g", *v)
}
