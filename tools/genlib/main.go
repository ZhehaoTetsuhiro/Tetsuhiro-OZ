// Command genlib generates the built-in geometric-optics example library in
// library/geo. Every prescription is built by an optical solve, measured with
// the OZ engine (tools/libcheck/measure), checked by the library's
// assertions and written as JSON whose Chinese description states the
// measured expectation clause. Candidates that cannot be solved or fail the
// assertions are skipped and listed in library/geo/FAILURES.md; the library
// is then re-validated by tools/libcheck, which writes library/geo/index.json.
//
//	go run ./tools/genlib            # writes library/geo/*.json and FAILURES.md
//	go run ./tools/genlib -dry       # measures and reports only
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"tetsuhiro-oz/tools/libcheck/measure"
)

const outDir = "library/geo"

func main() {
	dry := flag.Bool("dry", false, "measure and report, do not write files")
	flag.Parse()

	l := newLib()
	buildSinglets(l)
	buildAchromats(l)
	buildMulti(l)
	buildEyepieces(l)
	buildMicro(l)
	buildTelescopes(l)
	buildMirrors(l)
	buildInstruments(l)
	buildAspheres(l)
	buildPosed(l)

	if !*dry {
		if err := os.MkdirAll(outDir, 0o755); err != nil {
			fmt.Fprintln(os.Stderr, "genlib:", err)
			os.Exit(1)
		}
		// the library is regenerated as a whole: stale examples must not survive
		old, _ := filepath.Glob(filepath.Join(outDir, "*.json"))
		for _, f := range old {
			if filepath.Base(f) != "index.json" {
				if err := os.Remove(f); err != nil {
					fmt.Fprintln(os.Stderr, "genlib:", err)
					os.Exit(1)
				}
			}
		}
	}

	byMode := map[string]int{}
	byQuality := map[string]int{}
	byCat := map[string]int{}
	written := 0
	for _, it := range l.items {
		m, err := measure.Measure(it.p)
		if err != nil {
			skipf("%s: measure: %v", it.file, err)
			continue
		}
		if probs := measure.Check(it.p, m); len(probs) > 0 {
			skipf("%s: assertion: %s", it.file, strings.Join(probs, "; "))
			continue
		}
		it.p.Description = it.note + "。" + measure.ExpectClause(m)
		// the clause must read back to exactly the measured numbers
		mode, vals, err := measure.ParseClause(it.p.Description)
		if err != nil {
			skipf("%s: clause: %v", it.file, err)
			continue
		}
		if probs := measure.Compare(m, mode, vals); len(probs) > 0 {
			skipf("%s: clause mismatch: %s", it.file, strings.Join(probs, "; "))
			continue
		}
		it.p.Tags = append(it.p.Tags, "mode:"+m.Mode)
		if err := checkMetadata(it); err != nil {
			skipf("%s: %v", it.file, err)
			continue
		}
		if !*dry {
			if err := it.p.SaveJSON(filepath.Join(outDir, it.file)); err != nil {
				fmt.Fprintln(os.Stderr, "genlib: write:", err)
				os.Exit(1)
			}
		}
		written++
		byCat[it.p.Category]++
		byMode[m.Mode]++
		for _, t := range it.p.Tags {
			if strings.HasPrefix(t, "quality:") {
				byQuality[strings.TrimPrefix(t, "quality:")]++
			}
		}
	}

	registered := len(l.items)
	fmt.Printf("\n== 生成完成：登记候选 %d = 写入 %d + 跳过 %d；构造失败 %d（未登记）；合并 %d；改名 %d\n",
		registered, written, len(candFails), len(buildFails), len(merges), len(renames))
	printCounts("类别（写入）", byCat)
	printCounts("模式", byMode)
	printCounts("质量", byQuality)

	if !*dry {
		path := filepath.Join(outDir, "FAILURES.md")
		if err := writeFailures(path, registered, written); err != nil {
			fmt.Fprintln(os.Stderr, "genlib:", err)
			os.Exit(1)
		}
		fmt.Printf("失败清单已写入 %s：构造失败 %d 条，候选跳过 %d 条\n", path, len(buildFails), len(candFails))
	}
}

func printCounts(title string, m map[string]int) {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		fmt.Printf("  %s %-16s %4d\n", title, k, m[k])
	}
}

// checkMetadata verifies the descriptive invariants of one candidate.
func checkMetadata(it item) error {
	p := it.p
	if p.Category == "" {
		return fmt.Errorf("missing category")
	}
	if p.Difficulty < 1 || p.Difficulty > 5 {
		return fmt.Errorf("difficulty %d out of range", p.Difficulty)
	}
	if len(p.Tags) < 3 {
		return fmt.Errorf("too few tags")
	}
	if p.Name == "" || p.Description == "" {
		return fmt.Errorf("missing name/description")
	}
	for _, r := range it.file {
		if r > 127 {
			return fmt.Errorf("file name is not ASCII")
		}
	}
	return nil
}

// uncovered lists the families that are outside this batch.
var uncovered = []string{
	"自由曲面（XY 多项式离轴面，缺口 G3）与泽尼克面（G4）：依赖尚未合并的面型能力，本批未生成（若某例需要该能力则跳过并记录于此）",
	"施密特校正板与施密特-卡塞格林的解析 A4 校正板：需要校正板设计解，本次未做",
	"里奇-克雷蒂安（双曲面，两镜圆锥常数联立求解）与格利高里式：RC 需要圆锥求解；格利高里式的副镜共轭符号未核实",
	"DOE / 二元光学、热光学（dn/dT、热膨胀）、偏振（琼斯矩阵）、杂散光（二次反射）：本批未涉及",
	"卡塞格林的中心遮挡、牛顿式的完整 45° 折轴镜筒：本批只做单折轴镜折转与折叠臂透镜的折转光路",
	"公差进阶（灵敏度、逆灵敏度、补偿器）、膜系联合评价：本批未涉及",
}

// writeFailures writes the generation record. The two stages are reported
// separately so that the counts reconcile:
//
//	registered candidates = written + skipped (measurement/assertion/description)
//	construction failures = layouts or solves that produced no candidate
func writeFailures(path string, registered, written int) error {
	var b strings.Builder
	b.WriteString("# library/geo 生成失败清单\n\n")
	b.WriteString("由 `go run ./tools/genlib` 生成（确定性输出，无时间戳）；校验与索引由 `go run ./tools/libcheck` 完成。\n\n")
	b.WriteString("## 计数口径\n\n")
	fmt.Fprintf(&b, "- 登记候选（进入库的候选条目）：%d = 写入 %d + 候选跳过 %d\n", registered, written, len(candFails))
	fmt.Fprintf(&b, "- 构造失败（求解、布局或质量门槛失败，未产生候选，未登记）：%d\n", len(buildFails))
	fmt.Fprintf(&b, "- 合并（参数指纹与已登记条目相同，不计入候选）：%d\n", len(merges))
	fmt.Fprintf(&b, "- 同名文件改名（参数不同，不计入失败）：%d\n", len(renames))
	fmt.Fprintf(&b, "- 失败记录合计：%d = 构造失败 %d + 候选跳过 %d（逐条不重复）\n\n", len(buildFails)+len(candFails), len(buildFails), len(candFails))

	b.WriteString("## 候选跳过（已登记，测量 / 断言 / 描述复核未通过，未写入 JSON）\n\n")
	if len(candFails) == 0 {
		b.WriteString("（无）\n\n")
	}
	for i, f := range candFails {
		fmt.Fprintf(&b, "%d. %s\n", i+1, f)
	}
	b.WriteString("\n## 构造失败（求解/布局/质量门槛失败，未产生候选）\n\n")
	if len(buildFails) == 0 {
		b.WriteString("（无）\n\n")
	}
	for i, f := range buildFails {
		fmt.Fprintf(&b, "%d. %s\n", i+1, f)
	}
	b.WriteString("\n## 合并（参数完全相同，保留先出现的条目）\n\n")
	if len(merges) == 0 {
		b.WriteString("（无）\n\n")
	}
	for i, m := range merges {
		fmt.Fprintf(&b, "%d. %s\n", i+1, m)
	}
	b.WriteString("\n## 同名改名（参数不同，文件名加 -altN 后缀）\n\n")
	if len(renames) == 0 {
		b.WriteString("（无）\n\n")
	}
	for i, r := range renames {
		fmt.Fprintf(&b, "%d. %s\n", i+1, r)
	}
	b.WriteString("\n## 本批未覆盖（留给后续批次）\n\n")
	for _, u := range uncovered {
		fmt.Fprintf(&b, "- %s\n", u)
	}
	return os.WriteFile(path, []byte(b.String()), 0o644)
}
