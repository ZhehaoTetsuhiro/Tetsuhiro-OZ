package main

// `oz geo` subcommands: geometric MTF with field curvature (`geo mtf`) and the
// ghost-path table (`geo ghost`). Both call the shared core in analysis.go.

import (
	"flag"
	"fmt"
	"os"
	"time"

	"tetsuhiro-oz/lens"
)

// cmdGeo dispatches the geometry analysis subcommands.
func cmdGeo(args []string) error {
	if len(args) < 1 {
		return fmt.Errorf("geo: 需要子命令 mtf | ghost")
	}
	sub, rest := args[0], args[1:]
	switch sub {
	case "mtf":
		return cmdGeoMTF(rest)
	case "ghost":
		return cmdGeoGhost(rest)
	}
	return fmt.Errorf("geo: 未知子命令 %q（可用: mtf, ghost）", sub)
}

// cmdGeoGhost prints the double-reflection (ghost) path table of a
// prescription.
func cmdGeoGhost(args []string) error {
	fs := flag.NewFlagSet("geo ghost", flag.ExitOnError)
	lambda := fs.Float64("lambda", dLineNm, "波长 nm")
	field := fs.Float64("field", 0, "子午视场角（度）")
	imageZ := fs.Float64("image-z", 0, "像面轴向位置 mm（0 = 近轴像面）")
	pupilR := fs.Float64("pupil-radius", 0, "采样光瞳半径 mm（0 = 入瞳半径的 1e-3）")
	rings := fs.Int("rings", 0, "采样环数（Hexapolar，0 = 1）")
	minRel := fs.Float64("min-relative", 0, "相对能量阈值（低于它的路径不列出）")
	maxPaths := fs.Int("max-paths", 0, "最多列出前 N 条路径（0 = 全部）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	pos, err := parseFlagsAll(fs, args)
	if err != nil {
		return err
	}
	if len(pos) != 1 {
		return fmt.Errorf("geo ghost: 需要一个 <处方.json> 参数")
	}
	rxPath := pos[0]
	p, err := lens.LoadJSON(rxPath)
	if err != nil {
		return err
	}
	req := geoGhostRequest{
		Prescription: *p, LambdaNm: *lambda, FieldDeg: *field, ImageZMM: *imageZ,
		PupilRadiusMM: *pupilR, Rings: *rings, MinRelative: *minRel, MaxPaths: *maxPaths,
	}
	t0 := time.Now()
	res, err := runGeoGhost(req)
	if err != nil {
		return err
	}
	printGeoGhost(rxPath, res)
	fmt.Fprintf(os.Stderr, "（%.2fs）\n", time.Since(t0).Seconds())
	return writeJSONFile(*jsonOut, res)
}

// printGeoGhost prints the ghost table.
func printGeoGhost(path string, res *geoGhostResponse) {
	fmt.Printf("鬼像路径 %s: λ=%.4f nm，像面 z=%.6g mm，视场 θ=%g°，光瞳 z=%.6g mm r=%.6g mm，%d 条路径取 %d 条\n",
		path, res.LambdaNm, res.ImageZMM, res.FieldDeg, res.PupilZMM, res.PupilRMM, res.Enumerated, len(res.Paths))
	fmt.Printf("  %-18s %13s %10s %10s  %-22s %-22s %s\n",
		"面对 (I,J)", "相对能量 R_I·R_J", "R_I", "R_J", "近轴落点 (x,y) mm", "杂散光斑质心 (x,y) mm", "RMS mm")
	for _, r := range res.Paths {
		if r.Note != "" {
			fmt.Printf("  %-18s  %s\n", r.Pair, r.Note)
			continue
		}
		par := "n/a"
		if r.ParaxialXMM != nil {
			par = fmt.Sprintf("(%+.6g, %+.6g)", *r.ParaxialXMM, *r.ParaxialYMM)
		}
		fmt.Printf("  %-18s %13.6g %10.4g %10.4g  %-22s (%+.6g, %+.6g)    %.6g\n",
			r.Pair, r.Energy, r.ReflectI, r.ReflectJ, par,
			r.SpotCentroidXMM, r.SpotCentroidYMM, r.SpotRMSMM)
	}
}

// cmdGeoMTF prints the geometric MTF and the focal surfaces of a prescription.
func cmdGeoMTF(args []string) error {
	fs := flag.NewFlagSet("geo mtf", flag.ExitOnError)
	lambda := fs.Float64("lambda", dLineNm, "波长 nm")
	fieldsStr := fs.String("fields", "0", "子午视场角（度，逗号分隔，如 0,1,2）")
	freqStr := fs.String("cycles-per-mm", "10,20,30,40,50", "MTF 频率（cycles/mm，逗号分隔）")
	rays := fs.Int("rays", 0, "每个视场的光线数（0 = 2000）")
	bins := fs.Int("bins", 0, "LSF 分箱数（0 = 201）")
	imageZ := fs.Float64("image-z", 0, "像面轴向位置 mm（0 = 近轴像面）")
	rimFrac := fs.Float64("rim-frac", 0, "场曲/焦点对的入瞳半径比例（0 = 1e-3）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	pos, err := parseFlagsAll(fs, args)
	if err != nil {
		return err
	}
	if len(pos) != 1 {
		return fmt.Errorf("geo mtf: 需要一个 <处方.json> 参数")
	}
	rxPath := pos[0]
	fields, err := parseFloatList("--fields", *fieldsStr)
	if err != nil {
		return err
	}
	freqs, err := parseFloatList("--cycles-per-mm", *freqStr)
	if err != nil {
		return err
	}
	p, err := lens.LoadJSON(rxPath)
	if err != nil {
		return err
	}
	req := geoMTFRequest{
		Prescription: *p, LambdaNm: *lambda, FieldsDeg: fields, CyclesPerMM: freqs,
		Rays: *rays, Bins: *bins, ImageZMM: *imageZ, RimFrac: *rimFrac,
	}
	t0 := time.Now()
	res, err := runGeoMTF(req)
	if err != nil {
		return err
	}
	printGeoMTF(rxPath, res)
	fmt.Fprintf(os.Stderr, "（%.2fs）\n", time.Since(t0).Seconds())
	return writeJSONFile(*jsonOut, res)
}

// printGeoMTF prints the MTF table and the focal surfaces.
func printGeoMTF(path string, res *geoMTFResponse) {
	fmt.Printf("几何 MTF %s: λ=%.4f nm，像面 z=%.6g mm，每视场 %d 条光线，%d 箱（每视场各自分箱）\n",
		path, res.LambdaNm, res.ImageZMM, res.Rays, res.Bins)
	fmt.Printf("  %12s", "f (cyc/mm)")
	for _, f := range res.Frequencies {
		fmt.Printf(" %10.4g", f)
	}
	fmt.Println()
	for _, fld := range res.Fields {
		fmt.Printf("  视场 θ=%.6g°（tan y=%.6g，%d/%d 条光线到像面，RMS=%.4g mm，质心 (%.4g, %.4g) mm）\n",
			fld.ThetaDeg, fld.TanY, fld.RaysHit, fld.RaysTotal, fld.RMSRadiusMM, fld.CentroidXMM, fld.CentroidYMM)
		fmt.Printf("     子午 Δ=%.4g mm Nyquist=%.6g；弧矢 Δ=%.4g mm Nyquist=%.6g（cycles/mm）\n",
			fld.BinWidthTMM, fld.NyquistTCPM, fld.BinWidthSMM, fld.NyquistSCPM)
		fmt.Printf("  %12s", "MTF_t 子午")
		for _, v := range fld.MTFTangential {
			fmt.Printf(" %10.4f", v)
		}
		fmt.Println()
		fmt.Printf("  %12s", "MTF_s 弧矢")
		for _, v := range fld.MTFSagittal {
			fmt.Printf(" %10.4f", v)
		}
		fmt.Println()
	}
	if res.FociNote != "" {
		fmt.Printf("  场曲：%s\n", res.FociNote)
	}
	if len(res.Foci) > 0 {
		fmt.Printf("  场曲（焦点位置相对像面，mm；T=子午 S=弧矢 P=Petzval 面 (3S−T)/2）:\n")
		fmt.Printf("    %10s %14s %14s %14s\n", "θ (deg)", "T 离焦", "S 离焦", "P 离焦")
		for _, fp := range res.Foci {
			fmt.Printf("    %10.6g %14.6g %14.6g %14.6g\n", fp.ThetaDeg, fp.SagTMM, fp.SagSMM, fp.SagPMM)
		}
	}
}
