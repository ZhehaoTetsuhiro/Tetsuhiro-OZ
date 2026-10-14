// Command genlib generates the built-in geometric-optics example library in
// library/geo: every prescription is built by an optical solve, measured with
// the OZ engine, and written as JSON with a Chinese description that states
// the expected metrics. tools/libcheck then re-loads and re-validates every
// file independently and writes library/geo/index.json.
//
//	go run ./tools/genlib            # writes library/geo/*.json
//	go run ./tools/genlib -dry       # summary only
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"tetsuhiro-oz/lens"
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

	if !*dry {
		if err := os.MkdirAll(outDir, 0o755); err != nil {
			fmt.Fprintln(os.Stderr, "genlib:", err)
			os.Exit(1)
		}
	}

	type stat struct {
		mode    string
		quality string
	}
	byMode := map[string]int{}
	byQuality := map[string]int{}
	sweeps := 0
	written, failed := 0, 0

	for _, it := range l.items {
		m, err := measureEx(it.p)
		if err != nil {
			fmt.Printf("MEASURE-FAIL %-52s %v\n", it.file, err)
			failed++
			continue
		}
		it.p.Description = it.note + "。预期指标：" + expectText(m)
		it.p.Tags = append(it.p.Tags, "mode:"+m.Mode)
		if !*dry {
			if err := it.p.SaveJSON(filepath.Join(outDir, it.file)); err != nil {
				fmt.Printf("WRITE-FAIL   %-52s %v\n", it.file, err)
				failed++
				continue
			}
		}
		written++
		byMode[m.Mode]++
		for _, t := range it.p.Tags {
			if strings.HasPrefix(t, "quality:") {
				byQuality[strings.TrimPrefix(t, "quality:")]++
			}
			if strings.HasPrefix(t, "sweep:") {
				sweeps++
				break
			}
		}
	}

	// sanity: no two examples may share the same optical content
	sig := map[string]string{}
	dupes := 0
	for _, it := range l.items {
		s, err := signature(it.p)
		if err != nil {
			continue
		}
		if prev, ok := sig[s]; ok {
			fmt.Printf("DUPLICATE    %-52s == %s\n", it.file, prev)
			dupes++
			continue
		}
		sig[s] = it.file
	}
	if err := checkMetadata(l); err != nil {
		fmt.Println("META-FAIL:", err)
		failed++
	}

	fmt.Printf("\n== 生成完成：%d 个示例（写入 %d，失败 %d，重复 %d）\n", len(l.items), written, failed, dupes)
	cats := make([]string, 0, len(l.byCat))
	for c := range l.byCat {
		cats = append(cats, c)
	}
	sort.Strings(cats)
	for _, c := range cats {
		fmt.Printf("  %-16s %4d\n", c, l.byCat[c])
	}
	modes := make([]string, 0, len(byMode))
	for m := range byMode {
		modes = append(modes, m)
	}
	sort.Strings(modes)
	for _, m := range modes {
		fmt.Printf("  模式 %-12s %4d\n", m, byMode[m])
	}
	qs := make([]string, 0, len(byQuality))
	for q := range byQuality {
		qs = append(qs, q)
	}
	sort.Strings(qs)
	for _, q := range qs {
		fmt.Printf("  质量 %-12s %4d\n", q, byQuality[q])
	}
	fmt.Printf("  参数扫掠标记 %d\n", sweeps)
	if failed > 0 || dupes > 0 {
		os.Exit(1)
	}
}

// signature is the optical identity of a prescription: name, object
// position, stop and the exact card list. Two examples with the same
// signature would be the same lens under two file names.
func signature(p *lens.Prescription) (string, error) {
	type sig struct {
		Object float64      `json:"object"`
		Stop   int          `json:"stop"`
		Cards  []lens.Card  `json:"cards"`
	}
	buf, err := json.Marshal(sig{Object: p.ObjectZMM, Stop: p.Stop, Cards: p.Cards})
	if err != nil {
		return "", err
	}
	if math.IsInf(p.ObjectZMM, 0) {
		return "inf|" + string(buf), nil
	}
	return string(buf), nil
}

// checkMetadata verifies the library invariants that the index relies on.
func checkMetadata(l *lib) error {
	for _, it := range l.items {
		if it.p.Category == "" {
			return fmt.Errorf("%s: missing category", it.file)
		}
		if it.p.Difficulty < 1 || it.p.Difficulty > 5 {
			return fmt.Errorf("%s: difficulty %d out of range", it.file, it.p.Difficulty)
		}
		if len(it.p.Tags) < 3 {
			return fmt.Errorf("%s: too few tags", it.file)
		}
		if it.p.Name == "" || it.p.Description == "" {
			return fmt.Errorf("%s: missing name/description", it.file)
		}
		for _, r := range it.file {
			if r > 127 {
				return fmt.Errorf("%s: file name is not ASCII", it.file)
			}
		}
	}
	return nil
}
