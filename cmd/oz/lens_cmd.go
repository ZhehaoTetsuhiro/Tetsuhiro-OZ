package main

// `oz lens` subcommands: thermal drift and sensitivity analysis. Both call the
// shared core in analysis.go, so `oz lens sens` and POST /api/lens/sens run the
// same computation (the same holds for thermal).

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"

	"tetsuhiro-oz/lens"
)

// cmdLens dispatches the lens analysis subcommands.
func cmdLens(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("lens: 需要子命令 thermal | sens | group | parax")
	}
	sub, rest := args[0], args[1:]
	switch sub {
	case "thermal":
		return cmdLensThermal(rest)
	case "sens":
		return cmdLensSens(rest)
	case "group":
		return cmdLensGroup(rest)
	case "parax":
		return cmdLensParax(rest)
	}
	return fmt.Errorf("lens: 未知子命令 %q（可用: thermal, sens, group, parax）", sub)
}

// parseFlagsAll parses a command line whose positional argument (the
// prescription path) may sit anywhere among the flags, so all of
//
//	oz lens thermal <处方.json> --delta-t 10
//	oz lens thermal --delta-t 10 <处方.json>
//	oz lens thermal --delta-t 10 <处方.json> --json out.json
//
// work. flag.FlagSet.Parse stops at the first non-flag argument, so the tail
// behind each positional is parsed again and the positionals are collected in
// order.
func parseFlagsAll(fs *flag.FlagSet, args []string) ([]string, error) {
	var pos []string
	rest := args
	for {
		if err := fs.Parse(rest); err != nil {
			return nil, err
		}
		left := fs.Args()
		if len(left) == 0 {
			return pos, nil
		}
		pos = append(pos, left[0])
		rest = left[1:]
	}
}

// cmdLensThermal prints the first-order thermal drift of a prescription.
func cmdLensThermal(args []string) error {
	fs := flag.NewFlagSet("lens thermal", flag.ExitOnError)
	dt := fs.Float64("delta-t", 10, "温度变化 ΔT（K，相对 20 °C，非 0）")
	lambda := fs.Float64("lambda", dLineNm, "波长 nm")
	mount := fs.Float64("mount-cte", 0, "隔圈/镜筒的线膨胀系数（1/K，默认刚性镜筒）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	pos, err := parseFlagsAll(fs, args)
	if err != nil {
		return err
	}
	if len(pos) != 1 {
		return fmt.Errorf("lens thermal: 需要一个 <处方.json> 参数")
	}
	rxPath := pos[0]
	p, err := lens.LoadJSON(rxPath)
	if err != nil {
		return err
	}
	t0 := time.Now()
	res, err := runLensThermal(lensThermalRequest{Prescription: *p, DeltaTK: *dt, LambdaNm: *lambda, MountCTE: *mount})
	if err != nil {
		return err
	}
	fmt.Printf("热漂移 ΔT=%+g K（λ=%.4f nm，镜筒 α=%g 1/K）: %s\n", res.DeltaTK, res.LambdaNm, res.MountCTE, rxPath)
	fmt.Printf("  EFL  %+.6f mm → %+.6f mm   （%+.4g mm/K）\n", res.EFL0MM, res.EFL1MM, res.EFLPerKMM)
	fmt.Printf("  BFL  %+.6f mm → %+.6f mm   （%+.4g mm/K）\n", res.BFL0MM, res.BFL1MM, res.BFLPerKMM)
	fmt.Fprintf(os.Stderr, "（%.2fs）\n", time.Since(t0).Seconds())
	return writeJSONFile(*jsonOut, res)
}

// cmdLensSens prints the sensitivity table of a prescription and, on request,
// the inverse-sensitivity allocation and/or a compensator verification.
func cmdLensSens(args []string) error {
	fs := flag.NewFlagSet("lens sens", flag.ExitOnError)
	lambda := fs.Float64("lambda", dLineNm, "波长 nm")
	metric := fs.String("metric", metricEFL, "标量函数: efl（∂EFL/∂x）| merit（固定像面点列 merit）")
	relStep := fs.Float64("rel-step", 0, "最粗相对步长（0 = 1e-3）")
	fixed := fs.Bool("fixed-step", false, "只用固定步长求导，不跑 Richardson 阶梯（快，但不给收敛判据）")
	inverse := fs.Bool("inverse", false, "逆灵敏度：把 budget 的 merit 增量按灵敏度分配到各参数（自动用 merit 灵敏度）")
	budget := fs.Float64("budget", 0, "逆灵敏度的 merit 预算 ΔM（mm，必填 > 0）")
	mode := fs.String("mode", "rss", "分配模式: rss（统计）| linear（最坏情况）")
	maxCurv := fs.Float64("max-curv", 0, "曲率公差上限 1/mm（0 = 2e-5，≈R=50 mm 处的 0.05 mm 半径公差）")
	maxThk := fs.Float64("max-thk", 0, "厚度公差上限 mm（0 = 0.05）")
	maxDN := fs.Float64("max-dn", 0, "玻璃折射率公差上限 Δn（0 = 5e-4）")
	maxK := fs.Float64("max-k", 0, "圆锥常数公差上限 K（0 = 0.01）")
	maxCoef := fs.Float64("max-coef", 0, "非球面系数公差上限（0 = 1e-7）")
	compensate := fs.Bool("compensate", false, "补偿器验算：按 --perturb-* 扰动后重新优化补偿器，报告恢复率")
	pRadius := fs.Float64("perturb-radius", 0, "扰动量：每个有限半径 +dR（mm）")
	pThk := fs.Float64("perturb-thickness", 0, "扰动量：每个间隔 +dt（mm）")
	pIndex := fs.Float64("perturb-index", 0, "扰动量：每种玻璃 +Δn")
	compVars := fs.String("comp-vars", "", "补偿器变量（逗号分隔，如 s1.thk,s0.R；留空=只用像面补偿）")
	compRange := fs.String("comp-range", "", "补偿器半行程（mm，逗号分隔，与 --comp-vars 一一对应；0=默认）")
	imageShift := fs.Bool("image-shift", false, "把像面位置加入补偿器（配合 --comp-vars 使用）")
	shiftRange := fs.Float64("shift-range", 0, "像面补偿半行程 mm（0 = 1 mm）")
	compIter := fs.Int("comp-iter", 0, "补偿器 DLS 迭代上限（0 = 100）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	pos, err := parseFlagsAll(fs, args)
	if err != nil {
		return err
	}
	if len(pos) != 1 {
		return fmt.Errorf("lens sens: 需要一个 <处方.json> 参数")
	}
	rxPath := pos[0]
	if *inverse && !(*budget > 0) {
		return fmt.Errorf("lens sens: --inverse 需要 --budget（merit 增量预算 mm）")
	}
	if *inverse || *compensate {
		if *metric != metricMerit {
			*metric = metricMerit // --inverse/--compensate 使用 merit 灵敏度（runLensSens 也会校验）
		}
	}
	p, err := lens.LoadJSON(rxPath)
	if err != nil {
		return err
	}
	req := lensSensRequest{
		Prescription: *p, LambdaNm: *lambda, Metric: *metric, RelStep: *relStep, FixedStep: *fixed,
		Inverse: *inverse, BudgetMM: *budget, Mode: *mode,
		MaxCurvature: *maxCurv, MaxThickness: *maxThk, MaxIndex: *maxDN, MaxConic: *maxK, MaxCoeff: *maxCoef,
		Compensate: *compensate, PerturbRadiusMM: *pRadius, PerturbThicknessMM: *pThk, PerturbIndex: *pIndex,
		ImageShift: *imageShift, ShiftRangeMM: *shiftRange, MaxIter: *compIter,
	}
	if *compVars != "" {
		req.CompVars = strings.Split(*compVars, ",")
	}
	if *compRange != "" {
		for _, s := range strings.Split(*compRange, ",") {
			v, err := strconv.ParseFloat(strings.TrimSpace(s), 64)
			if err != nil {
				return fmt.Errorf("lens sens: --comp-range 项 %q 不是数字", s)
			}
			req.CompRangesMM = append(req.CompRangesMM, v)
		}
	}
	t0 := time.Now()
	res, err := runLensSens(req)
	if err != nil {
		return err
	}
	printSensTable(rxPath, res)
	if res.Alloc != nil {
		printAllocTable(res.Alloc)
	}
	if res.Comp != nil {
		printCompTable(req, res.Comp)
	}
	fmt.Fprintf(os.Stderr, "（%.2fs）\n", time.Since(t0).Seconds())
	return writeJSONFile(*jsonOut, res)
}

// printSensTable prints the sensitivity rows.
func printSensTable(path string, res *lensSensResponse) {
	unit := "mm / 参数单位"
	if res.Metric == metricEFL {
		unit = "mm EFL / 参数单位"
	}
	fmt.Printf("灵敏度表 %s: metric=%s，λ=%.4f nm，f(名义)=%.6g mm，%d 次求值（单位: %s）\n",
		path, res.Metric, res.LambdaNm, res.Nominal, res.Evals, unit)
	fmt.Printf("  %-10s %14s %12s %14s %12s %12s  %s\n", "参数", "∂f/∂x", "不确定度", "naive", "相对步长", "步长", "收敛")
	for _, r := range res.Rows {
		cv := "no"
		if r.Converged {
			cv = "yes"
		}
		fmt.Printf("  %-10s %14.6g %12.3g %14.6g %12.3g %12.3g  %s\n",
			r.Param, r.Deriv, r.Err, r.Naive, r.RelStep, r.Step, cv)
	}
}

// printAllocTable prints the inverse-sensitivity allocation.
func printAllocTable(a *allocTable) {
	fmt.Printf("逆灵敏度分配（%s，预算 ΔM=%.6g mm，线性模型预测 %.6g mm，相对误差 %.2e）:\n",
		a.Mode, a.BudgetMM, a.Predicted, a.RelErr)
	fmt.Printf("  %-10s %14s %14s %14s  %s\n", "参数", "∂M/∂x", "分配公差", "预算份额", "封顶")
	for _, r := range a.Rows {
		cap := ""
		if r.Capped {
			cap = "yes"
		}
		fmt.Printf("  %-10s %14.6g %14.6g %14.6g  %s\n", r.Param, r.Sens, r.Tol, r.Share, cap)
	}
}

// printCompTable prints the compensator verification.
func printCompTable(req lensSensRequest, c *compTable) {
	fmt.Printf("补偿器验算（扰动: dR=%g mm, dt=%g mm, dn=%g"+
		"; 补偿器: %s）:\n", req.PerturbRadiusMM, req.PerturbThicknessMM, req.PerturbIndex, compVarLabel(req))
	fmt.Printf("  名义 merit            %.6g mm（补偿后 %.6g mm）\n", c.NominalMM, c.CompensatedNominalMM)
	fmt.Printf("  扰动后 merit          %.6g mm\n", c.PerturbedMM)
	fmt.Printf("  补偿后 merit          %.6g mm（%d 次迭代）\n", c.CompensatedMM, c.Iterations)
	if c.Recovery != nil {
		fmt.Printf("  恢复率 recovery       %.4f\n", *c.Recovery)
	} else {
		fmt.Printf("  恢复率 recovery       n/a（扰动未使 merit 恶化）\n")
	}
	fmt.Printf("  像面补偿量            %+.6g mm\n", c.ImageShiftMM)
}

// compVarLabel names the compensators of a request for the report.
func compVarLabel(req lensSensRequest) string {
	parts := append([]string(nil), req.CompVars...)
	if len(parts) == 0 || req.ImageShift {
		parts = append(parts, "像面(shift)")
	}
	return strings.Join(parts, ",")
}

// writeJSONFile writes the response as indented JSON when a path is given.
func writeJSONFile(path string, v any) error {
	if path == "" {
		return nil
	}
	buf, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return err
	}
	if err := os.WriteFile(path, append(buf, '\n'), 0o644); err != nil {
		return err
	}
	fmt.Fprintln(os.Stderr, "已写入", path)
	return nil
}

// parseFloatList parses a comma-separated list of numbers (CLI list options).
func parseFloatList(flagName, s string) ([]float64, error) {
	if strings.TrimSpace(s) == "" {
		return nil, nil
	}
	parts := strings.Split(s, ",")
	out := make([]float64, len(parts))
	for i, p := range parts {
		v, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
		if err != nil || math.IsNaN(v) || math.IsInf(v, 0) {
			return nil, fmt.Errorf("%s: 项 %q 不是有限数字", flagName, p)
		}
		out[i] = v
	}
	return out, nil
}
