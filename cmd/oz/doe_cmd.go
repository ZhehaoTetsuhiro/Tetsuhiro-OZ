package main

// `oz doe focus`: first-order efficiency and focus deviation of an N-level
// focusing DOE. Calls the shared core in analysis.go (runDOEFocus), which the
// HTTP endpoint POST /api/doe/focus uses as well.

import (
	"flag"
	"fmt"
	"os"
	"time"
)

// cmdDoe dispatches the DOE analysis subcommands.
func cmdDoe(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("doe: 需要子命令 focus")
	}
	sub, rest := args[0], args[1:]
	switch sub {
	case "focus":
		return cmdDoeFocus(rest)
	}
	return fmt.Errorf("doe: 未知子命令 %q（可用: focus）", sub)
}

// cmdDoeFocus prints the first-order efficiency and the focus deviation of a
// focusing DOE.
func cmdDoeFocus(args []string) error {
	fs := flag.NewFlagSet("doe focus", flag.ExitOnError)
	focal := fs.Float64("focal-mm", 0, "设计焦距 mm（必填，> 0）")
	levels := fs.Int("levels", 0, "相位级数 N（必填，≥ 2）")
	lambda := fs.Float64("lambda", 0, "波长 nm（0 = 632.8）")
	waist := fs.Float64("waist-um", 0, "DOE 上高斯照明 1/e 振幅半径 µm（0 = 由 Fresnel 数给出）")
	nf := fs.Float64("fresnel-number", 0, "Fresnel 数 N_F = w0²/(λf)（0 = 32）")
	grid := fs.Int("grid-size", 0, "网格 M（2 的幂，0 = 1024；2048 给参考精度，≈3 s）")
	pix := fs.Float64("pixel-um", 0, "DOE 像素间距 µm（0 = 8·w0/M，与 doe 测试同配置）")
	win := fs.Float64("window-radii", 0, "效率窗口半径（焦斑半径的倍数，0 = 2.5）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if fs.NArg() != 0 {
		return fmt.Errorf("doe focus: 多余的参数 %q（本子命令不需要处方文件）", fs.Arg(0))
	}
	req := doeFocusRequest{
		FocalMM: *focal, Levels: *levels, LambdaNm: *lambda, WaistUM: *waist,
		FresnelNumber: *nf, GridSize: *grid, PixelUM: *pix, WindowRadii: *win,
	}
	t0 := time.Now()
	res, err := runDOEFocus(req)
	if err != nil {
		return err
	}
	printDOEFocus(res)
	fmt.Fprintf(os.Stderr, "（%.2fs）\n", time.Since(t0).Seconds())
	return writeJSONFile(*jsonOut, res)
}

// printDOEFocus prints the efficiency and focus table.
func printDOEFocus(res *doeFocusResponse) {
	fmt.Printf("量化相位 DOE: f=%.6g mm（%.6g µm），N=%d 级，λ=%.4g nm，M=%d，像素 %.6g µm，w0=%.6g µm（N_F=%.4g）\n",
		res.FocalMM, res.FocalUM, res.Levels, res.LambdaNm, res.GridSize, res.PixelUM, res.WaistUM, res.FresnelNumber)
	fmt.Printf("  一级效率 η        %.6f（解析 η_N=%.6f，相对差 %.3e）\n", res.Eta, res.EtaAnalytic, res.EtaRelErr)
	fmt.Printf("  焦点偏差          %+.6g µm（%+.4g 像素；轴对称峰值 %.6g µm，设计焦面 %.6g µm）\n",
		res.FocusDeviationUM, res.FocusDeviationPX, res.AxialPeakUM, res.FocalUM)
	fmt.Printf("  焦面峰值横向偏移  %.4g 像素\n", res.PeakOffsetPX)
	fmt.Printf("  焦斑半径 w_f      %.6g µm；窗口份额（理想高斯）%.6f\n", res.SpotRadiusUM, res.WindowFraction)
}
