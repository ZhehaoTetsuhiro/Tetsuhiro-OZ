package main

// `oz scan`: sweep one parameter of a wave-optics configuration and report
// one metric per value. It calls the shared core in the optics package
// (optics.RunScan), which the HTTP endpoint POST /api/scan uses as well, so
// the command line and the service cannot drift apart.

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"

	"tetsuhiro-oz/optics"
)

// cmdScan sweeps one configuration parameter and prints the resulting curve.
func cmdScan(args []string) error {
	fs := flag.NewFlagSet("scan", flag.ExitOnError)
	param := fs.String("param", "", "参数路径（必填）：scene.components[0].roll | scene.components[0].params.f | scene.components[0].shape.params.width | sources[0].params.tilt_x | wavelength ...")
	values := fs.String("values", "", "取值列表（必填，逗号分隔；支持 a-b 等步长范围）")
	metric := fs.String("metric", "power", "指标: power|peak|centroid_x|centroid_y|rms_x|rms_y|strehl|intensity_max")
	plane := fs.String("plane", "", "读取哪个记录平面（id 或 label；默认第一个）")
	jsonOut := fs.String("json", "", "把结果写成 JSON 到该路径")
	if err := fs.Parse(args); err != nil {
		return err
	}
	if *param == "" {
		return fmt.Errorf("scan: --param 必填（如 scene.components[0].roll）")
	}
	if *values == "" {
		return fmt.Errorf("scan: --values 必填（如 0,15,30,45 或 0-90-5）")
	}
	vals, err := parseScanValues(*values)
	if err != nil {
		return err
	}
	rest := fs.Args()
	if len(rest) != 1 {
		return fmt.Errorf("scan: 需要一个配置 JSON 文件（用法: oz scan --param ... --values ... config.json）")
	}
	buf, err := os.ReadFile(rest[0])
	if err != nil {
		return fmt.Errorf("scan: %w", err)
	}
	var cfg optics.Config
	if err := json.Unmarshal(buf, &cfg); err != nil {
		return fmt.Errorf("scan: %s: %w", rest[0], err)
	}
	spec := optics.ScanSpec{Param: *param, Values: vals, Metric: *metric, Plane: *plane}
	t0 := time.Now()
	rep, err := optics.RunScan(cfg, spec)
	if err != nil {
		return err
	}
	printScan(rep, *param, *metric)
	fmt.Fprintf(os.Stderr, "（%d 个点，%.2fs）\n", len(rep.Points), time.Since(t0).Seconds())
	return writeJSONFile(*jsonOut, rep)
}

// parseScanValues reads a comma-separated list of numbers, with an optional
// "a-b" or "a-b-step" range: "0,15,30", "0-90", "0-90-5", "0,1.5e-3,3e-3".
//
// The separator is a '-' that is *not* part of a number: a '-' immediately
// after another '-' or at the start of a token is a sign, not a range, and
// scientific notation ("5e-7") keeps its '-' too. The split therefore runs on
// '-' characters that follow a digit.
func parseScanValues(s string) ([]float64, error) {
	var out []float64
	for _, part := range strings.Split(s, ",") {
		part = strings.TrimSpace(part)
		if part == "" {
			continue
		}
		bits := splitRanges(part)
		if len(bits) == 1 {
			v, err := strconv.ParseFloat(bits[0], 64)
			if err != nil {
				return nil, fmt.Errorf("scan: 无法解析取值 %q", part)
			}
			out = append(out, v)
			continue
		}
		if len(bits) > 3 {
			return nil, fmt.Errorf("scan: 范围 %q 最多是 a-b-step 三段", part)
		}
		lo, err1 := strconv.ParseFloat(bits[0], 64)
		hi, err2 := strconv.ParseFloat(bits[1], 64)
		step := 1.0
		var err3 error
		if len(bits) == 3 {
			step, err3 = strconv.ParseFloat(bits[2], 64)
		}
		if err1 != nil || err2 != nil || err3 != nil {
			return nil, fmt.Errorf("scan: 无法解析范围 %q", part)
		}
		if step == 0 {
			return nil, fmt.Errorf("scan: 范围 %q 的步长不能为 0", part)
		}
		if (hi > lo) != (step > 0) {
			return nil, fmt.Errorf("scan: 范围 %q 的步长方向与起止相反", part)
		}
		n := int(math.Floor((hi-lo)/step + 1e-9))
		if n < 0 {
			n = 0
		}
		if n+1 > 4096 {
			return nil, fmt.Errorf("scan: 范围 %q 展开成 %d 个点，超过 4096 上限", part, n+1)
		}
		for k := 0; k <= n; k++ {
			out = append(out, lo+step*float64(k))
		}
	}
	if len(out) == 0 {
		return nil, fmt.Errorf("scan: --values 没有给出任何取值")
	}
	return out, nil
}

// splitRanges splits "a-b-step" on the range separators, keeping signs and
// scientific-notation exponents intact.
func splitRanges(s string) []string {
	var out []string
	start := 0
	for i := 1; i < len(s); i++ {
		if s[i] != '-' {
			continue
		}
		// A '-' after a digit (or a closing '.') separates two numbers;
		// anything else is a sign or an exponent.
		prev := s[i-1]
		if prev >= '0' && prev <= '9' || prev == '.' {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	return append(out, s[start:])
}

// printScan prints the sweep as an aligned table: value, metric, and the
// per-point error when a value could not be simulated.
func printScan(rep *optics.ScanReport, param, metric string) {
	fmt.Printf("参数扫描 %s（指标 %s）\n", param, metricOr(rep.Metric, metric))
	for i, p := range rep.Points {
		if p.Error != "" {
			fmt.Printf("  [%3d] %-14.6g  错误: %s\n", i, p.Value, firstLine(p.Error))
			continue
		}
		fmt.Printf("  [%3d] %-14.6g  %.6g\n", i, p.Value, p.Metric)
	}
	if len(rep.Warnings) > 0 {
		fmt.Printf("警告（去重）:\n")
		for _, w := range rep.Warnings {
			fmt.Printf("  %s\n", w)
		}
	}
}

// metricOr returns the report's metric, or the requested one when the report
// carries the default.
func metricOr(got, want string) string {
	if got != "" {
		return got
	}
	return want
}

// firstLine keeps only the first line of a (possibly multi-line) error text.
func firstLine(s string) string {
	if i := strings.IndexByte(s, '\n'); i >= 0 {
		return strings.TrimSpace(s[:i])
	}
	return s
}
