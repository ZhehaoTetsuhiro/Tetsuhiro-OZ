// Command oz is the Tetsuhiro-OZ command line: autonomous lens design,
// prescription evaluation and deterministic Monte-Carlo tolerance analysis
// on top of the OZ kernel.
//
// Subcommands:
//
//	oz design achromat --efl 100 [--crown N-BK7 --flint F2 ...] [--out lens.json]
//	oz eval lens.json
//	oz tol lens.json --radius 0.05 --thickness 0.05 --index 0.0005 --runs 2000 \
//	    --max-rms-on 0.02 --max-focal 0.1 --target-efl 100
//	oz coating design|spectrum ...
//	oz serve [-addr :1130] [WOS options]   the one HTTP service: API + GUI (docs/API.md)
//	oz version
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"strings"
	"time"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/plot"
)

// version of the OZ preview; keep in sync with /VERSION at the repo root.
const version = "0.1.0-Preview"

func main() {
	if len(os.Args) < 2 {
		usage()
		os.Exit(2)
	}
	var err error
	switch os.Args[1] {
	case "design":
		err = cmdDesign(os.Args[2:])
	case "eval":
		err = cmdEval(os.Args[2:])
	case "tol", "tolerance":
		err = cmdTol(os.Args[2:])
	case "lens":
		err = cmdLens(os.Args[2:])
	case "geo":
		err = cmdGeo(os.Args[2:])
	case "doe":
		err = cmdDoe(os.Args[2:])
	case "coating":
		err = cmdCoating(os.Args[2:])
	case "serve":
		err = cmdServe(os.Args[2:])
	case "version", "-v", "--version":
		fmt.Println("oz " + version)
	case "help", "-h", "--help":
		usage()
	default:
		fmt.Fprintf(os.Stderr, "oz: unknown command %q\n\n", os.Args[1])
		usage()
		os.Exit(2)
	}
	if err != nil {
		fmt.Fprintln(os.Stderr, "oz:", err)
		os.Exit(1)
	}
}

func usage() {
	fmt.Print(`oz — Tetsuhiro-OZ 光学设计命令行 ` + version + `

用法:（Go 惯例：选项在前，文件参数在最后；lens/geo 分析子命令的文件放前放后都可以）
  oz design achromat [选项] --out lens.json
      自主设计消色差双胶合（薄透镜解 → 厚透镜牛顿消色差 → DLS 抛光）。
      选项: --efl 焦距mm(必) --crown --flint --semi --field 视场角(度) --iter --no-polish --out
  oz eval <处方.json>
      打印一阶性质 / F–C 色差 / 点列图报告。
  oz tol [选项] <处方.json>
      选项: --radius R --thickness T --index dN [--conic K] [--coef A]
            --runs N --seed S --target-efl E
            --max-focal M --max-color M --max-rms-on M --max-rms-off M [--json out.json]
      Monte-Carlo 公差分析（确定性；至少给一个 max-* 判据）。
  oz coating design --template <t> --out stack.json [--lambda 550 --from 420 --to 680
      --pairs 12 --cutoff 550 --inc Air --sub BK7]
      模板: single-ar | vcoat | broadband-ar | hr | lwp | swp | splitter（自主设计膜系）
  oz coating spectrum <膜系.json> [--from 400 --to 700 --step 10 --angle 0 --pol s|p|u]
      打印 R/T/A 光谱表。
  oz lens thermal [--delta-t 10 --lambda 587.5618 --mount-cte 0] <处方.json>
      热漂移：ΔT 下 EFL/BFL 的变化与每开尔文斜率（G9）。
  oz lens sens [--metric efl|merit --inverse --budget ΔM --compensate ...] <处方.json>
      灵敏度表 ∂f/∂x；--inverse 逆灵敏度公差分配；--compensate 补偿器验算（G6）。
  oz lens group show|layout|image|afocal|dispersion|set|to-prescription <group.json> ...
      元件级透镜组（含镜、3D 折转，见 docs/API.md §4）：
      show      打印元件表；layout 逐面顶点/轴向/法向/折转（--json）；
      image     实光线成像位置（--field-deg/--field-tan/--lambda/--pupil/--scan/--range/--json）；
      afocal    角放大率 Γ 与准直残差（--lambda/--pupil/--json）；
      dispersion 两波长焦移（--lambda-a F --lambda-b C，--json）；
      set       逐参数编辑（-p e0.front.R=120 -p e2.tilt_x_deg=45，-o out.json）；
      to-prescription  转成卡片处方（-o rx.json）。
  oz lens parax [--object-z mm --field-tan t --lambda nm] <处方.json>   （含镜系统的一阶）
      EFL/FFL/BFL/像面/放大率/入出瞳/镜数/Parity/Afocal/Γ/F–C 焦移（--json）。
  oz geo mtf [--fields 0,1,2 --cycles-per-mm 10,20,30,40,50 --rays N --bins N] <处方.json>
      几何 MTF（子午/弧矢）与场曲（T/S/Petzval 离焦）摘要（G11）。
  oz geo ghost [--min-relative R --max-paths N --rings N] <处方.json>   （G10，见 docs/API.md）
      双反射鬼像路径表：面对 (I,J)、相对能量 R_I·R_J、近轴落点与杂散光斑。
  oz doe focus --focal-mm F --levels N [--lambda 632.8 --grid-size 1024 ...]   （G8）
      衍射光学元件一级效率 η 与解析 η_N 对比、焦点偏差。
  oz serve [选项]
      启动唯一的 HTTP 服务：GUI 在 /，API 在 /api/。WOS 的全部参数与环境变量同样适用：
      -addr（默认 :1130；WOS 旧端口用 -addr :1120）  -max-run-mb 512  -elements 目录1,目录2
      -no-elements-watch  -gpu（或环境变量 WOS_GPU=1）  -check-elements  -gen-go 文件 [-gen-name 名]
      环境变量 WOS_MUSA_HELPER / WOS_MUSA_SOCKET 指定 MUSA helper（见 docs/GPU.md）。
      几何/膜系：/api/version、/api/lens/{design,eval,tol}、/api/coating/{design,spectrum}
      波动/量子：/api/catalog、/api/health、/api/elements、/api/validate、/api/simulate、
                 /api/runs/...、/api/quantum、/api/convert（完整说明见 docs/API.md）
  oz version | oz help
`)
}

func cmdDesign(args []string) error {
	fs := flag.NewFlagSet("design", flag.ExitOnError)
	sub := ""
	if len(args) > 0 && !strings.HasPrefix(args[0], "-") {
		sub, args = args[0], args[1:]
	}
	switch sub {
	case "", "achromat":
	default:
		return fmt.Errorf("design: unknown lens type %q (have: achromat)", sub)
	}
	efl := fs.Float64("efl", 0, "有效焦距 mm（必填）")
	crown := fs.String("crown", "N-BK7", "冕玻璃名 (geo 玻璃目录)")
	flint := fs.String("flint", "F2", "火石玻璃名")
	semi := fs.Float64("semi", 0, "半口径 mm（0 = efl/10）")
	fieldDeg := fs.Float64("field", 0.5, "离轴视场角（度）")
	iter := fs.Int("iter", 0, "DLS 迭代上限（0 = 默认）")
	noPolish := fs.Bool("no-polish", false, "只做解析解 + 厚透镜消色差，不抛光")
	out := fs.String("out", "", "处方 JSON 输出路径")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *efl <= 0 {
		return fmt.Errorf("design achromat: --efl 必填且为正")
	}
	opts := &lens.DesignOptions{SkipPolish: *noPolish, MaxIter: *iter}
	if *fieldDeg > 0 {
		opts.Fields = []lens.Field{{TanX: 0, TanY: 0}, {TanX: 0, TanY: math.Tan(*fieldDeg * math.Pi / 180)}}
	}
	t0 := time.Now()
	p, res, err := lens.DesignAchromat(*efl, *crown, *flint, *semi, opts)
	if err != nil {
		return err
	}
	if res != nil {
		fmt.Fprintf(os.Stderr, "optimizer: %d 次接受迭代 / %d 次求值, merit=%.6g, %.2fs\n",
			res.Iter, res.Evals, res.Merit, time.Since(t0).Seconds())
	}
	if rep, err := lens.Describe(p); err == nil {
		fmt.Println(rep)
	}
	if *out != "" {
		if err := p.SaveJSON(*out); err != nil {
			return err
		}
		fmt.Fprintln(os.Stderr, "已写入", *out)
	}
	return nil
}

func cmdEval(args []string) error {
	fs := flag.NewFlagSet("eval", flag.ExitOnError)
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return fmt.Errorf("eval: 需要一个 <处方.json> 参数")
	}
	p, err := lens.LoadJSON(fs.Arg(0))
	if err != nil {
		return err
	}
	rep, err := lens.Describe(p)
	if err != nil {
		return err
	}
	fmt.Println(rep)
	return nil
}

func cmdTol(args []string) error {
	fs := flag.NewFlagSet("tol", flag.ExitOnError)
	radius := fs.Float64("radius", 0, "半径公差 mm（±，均匀）")
	thick := fs.Float64("thickness", 0, "厚度公差 mm（±）")
	conic := fs.Float64("conic", 0, "圆锥常数公差（±）")
	coef := fs.Float64("coef", 0, "非球面系数公差（±）")
	index := fs.Float64("index", 0, "玻璃熔炼 Δn（±，整批共享）")
	runs := fs.Int("runs", 2000, "Monte-Carlo 次数")
	seed := fs.Int64("seed", 1, "随机种子（确定性）")
	maxOn := fs.Float64("max-rms-on", 0, "轴上点列 RMS 判据 mm")
	maxOff := fs.Float64("max-rms-off", 0, "离轴点列 RMS 判据 mm")
	maxColor := fs.Float64("max-color", 0, "F–C 色差判据 mm")
	maxFocal := fs.Float64("max-focal", 0, "焦距误差判据 mm")
	target := fs.Float64("target-efl", 0, "目标焦距 mm（配合 --max-focal）")
	jsonOut := fs.String("json", "", "报告 JSON 输出路径")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return fmt.Errorf("tol: 需要一个 <处方.json> 参数")
	}
	p, err := lens.LoadJSON(fs.Arg(0))
	if err != nil {
		return err
	}
	spec := lens.ToleranceSpec{
		RadiusMM: *radius, ThicknessMM: *thick,
		ConicK: *conic, CoefAbs: *coef, IndexDelta: *index,
	}
	crit := lens.ToleranceCriteria{
		TargetEFLMM: *target, MaxFocalErrMM: *maxFocal,
		MaxColorShiftMM: *maxColor, MaxRMSOnAxisMM: *maxOn, MaxRMSOffAxisMM: *maxOff,
	}
	t0 := time.Now()
	rep, err := lens.RunTolerances(p, spec, crit, &lens.ToleranceOptions{
		Runs: *runs, Seed: *seed,
	})
	if err != nil {
		return err
	}
	fmt.Printf("%s\n  （%d 次运行，%.2fs）\n", rep.String(), rep.Runs, time.Since(t0).Seconds())
	if *jsonOut != "" {
		buf, err := json.MarshalIndent(rep, "", "  ")
		if err != nil {
			return err
		}
		if err := os.WriteFile(*jsonOut, append(buf, '\n'), 0o644); err != nil {
			return err
		}
		fmt.Fprintln(os.Stderr, "已写入", *jsonOut)
	}
	return nil
}

// cmdCoating dispatches the coating subcommands.
func cmdCoating(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("coating: 需要子命令 design | spectrum")
	}
	sub, rest := args[0], args[1:]
	switch sub {
	case "design":
		return cmdCoatingDesign(rest)
	case "spectrum":
		return cmdCoatingSpectrum(rest)
	}
	return fmt.Errorf("coating: 未知子命令 %q（可用: design, spectrum）", sub)
}

// cmdCoatingDesign runs one of the built-in coating templates and saves the
// designed stack.
func cmdCoatingDesign(args []string) error {
	fs := flag.NewFlagSet("coating design", flag.ExitOnError)
	template := fs.String("template", "", "模板: single-ar|vcoat|broadband-ar|hr|lwp|swp|splitter")
	lambda := fs.Float64("lambda", 550, "设计波长 nm")
	from := fs.Float64("from", 420, "宽带 AR 起始波长 nm")
	to := fs.Float64("to", 680, "宽带 AR 结束波长 nm")
	pairs := fs.Int("pairs", 12, "高反堆 (H L) 对数")
	cutoff := fs.Float64("cutoff", 0, "边沿滤光片截止波长 nm")
	inc := fs.String("inc", "Air", "入射介质（材料库名）")
	sub := fs.String("sub", "BK7", "衬底（材料库名）")
	out := fs.String("out", "", "膜系 JSON 输出路径（必填）")
	if err := fs.Parse(args); err != nil {
		return err
	}
	incM, ok := coating.LookupMaterial(*inc)
	if !ok {
		return fmt.Errorf("coating design: 未知入射介质 %q", *inc)
	}
	subM, ok := coating.LookupMaterial(*sub)
	if !ok {
		return fmt.Errorf("coating design: 未知衬底 %q", *sub)
	}
	var (
		stack coating.Stack
		res   *coating.Result
		err   error
	)
	switch *template {
	case "single-ar":
		stack, res, err = coating.SingleLayerAR(incM, subM, *lambda, coating.Options{})
	case "vcoat":
		stack, res, err = coating.VCoat(incM, subM, *lambda, coating.Options{})
	case "broadband-ar":
		stack, res, err = coating.BroadbandAR(incM, subM, *from, *to, coating.Options{})
	case "hr":
		stack, res, err = coating.HighReflector(incM, subM, *lambda, *pairs, coating.Options{})
	case "lwp":
		if !(*cutoff > 0) {
			return fmt.Errorf("coating design: lwp 需要 --cutoff")
		}
		stack, res, err = coating.LongWavePass(incM, subM, *cutoff, coating.Options{})
	case "swp":
		if !(*cutoff > 0) {
			return fmt.Errorf("coating design: swp 需要 --cutoff")
		}
		stack, res, err = coating.ShortWavePass(incM, subM, *cutoff, coating.Options{})
	case "splitter":
		stack, res, err = coating.Splitter50(incM, subM, *lambda, coating.Options{})
	default:
		return fmt.Errorf("coating design: 未知模板 %q（可用: single-ar, vcoat, broadband-ar, hr, lwp, swp, splitter）", *template)
	}
	if err != nil {
		return err
	}
	if res != nil {
		fmt.Print(res.Summary())
	}
	if *out == "" {
		return fmt.Errorf("coating design: 需要 --out 指定膜系 JSON 输出路径")
	}
	if err := stack.SaveJSON(*out); err != nil {
		return err
	}
	fmt.Fprintln(os.Stderr, "已写入", *out)
	return nil
}

// cmdCoatingSpectrum prints the R/T/A spectrum of a saved stack.
func cmdCoatingSpectrum(args []string) error {
	fs := flag.NewFlagSet("coating spectrum", flag.ExitOnError)
	from := fs.Float64("from", 400, "起始波长 nm")
	to := fs.Float64("to", 700, "结束波长 nm")
	step := fs.Float64("step", 10, "波长步长 nm")
	angle := fs.Float64("angle", 0, "入射角（度）")
	pol := fs.String("pol", "u", "偏振 s|p|u")
	svgOut := fs.String("svg", "", "把 R/T 光谱写成 SVG 图到该路径（可嵌入报告/网页）")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 1 {
		return fmt.Errorf("coating spectrum: 需要一个 <膜系.json> 参数")
	}
	stack, err := coating.LoadStackJSON(fs.Arg(0))
	if err != nil {
		return err
	}
	p, err := coating.ParsePolarization(*pol)
	if err != nil {
		return err
	}
	rows, err := stack.SpectrumRange(*from, *to, *step, *angle, p)
	if err != nil {
		return err
	}
	fmt.Printf("# %d layers, %.1f nm total, incidence %s, %.1f deg, pol %s\n",
		len(stack.Layers), stack.TotalThicknessNm(), stack.Incident.Name, *angle, p)
	fmt.Printf("# lambda_nm	R	T	A\n")
	for _, r := range rows {
		fmt.Printf("%.0f	%.6f	%.6f	%.6f\n", r.LambdaNm, r.R, r.T, r.A)
	}
	if *svgOut != "" {
		ptsR := make([][2]float64, 0, len(rows))
		ptsT := make([][2]float64, 0, len(rows))
		for _, r := range rows {
			ptsR = append(ptsR, [2]float64{r.LambdaNm, r.R})
			ptsT = append(ptsT, [2]float64{r.LambdaNm, r.T})
		}
		ch := plot.NewChart(720, 420)
		ch.SetTitle(fmt.Sprintf("R/T spectrum: %d layers, %.1f deg, pol %s", len(stack.Layers), *angle, p))
		ch.SetLabels("wavelength (nm)", "R / T")
		ch.SetRange(0, 0, 0, 1.02) // x 自动、y 固定 0..1.02
		ch.AddLine("R", "#ff6b6b", ptsR)
		ch.AddLine("T", "#4dabf7", ptsT)
		if err := ch.WriteFile(*svgOut); err != nil {
			return err
		}
		fmt.Printf("# SVG 图已写入 %s\n", *svgOut)
	}
	return nil
}
