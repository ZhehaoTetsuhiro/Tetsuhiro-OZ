// Command libcheck validates the built-in geometric-optics example library
// (library/geo). For every *.json file (index.json excluded) it
//
//  1. reloads the file with lens.LoadJSON, which re-checks the JSON contract,
//     the glass names and the card geometry;
//  2. checks the metadata: name, category, description, difficulty 1..5,
//     at least three tags, an ASCII file name;
//  3. measures the example with tools/libcheck/measure (the OZ engine) and
//     applies that package's physical-sanity assertions;
//  4. parses the expectation clause of the description and compares every
//     number with the fresh measurement (relative tolerance 1e-6, see
//     measure.CompareTol);
//  5. checks that names and parameter fingerprints are unique across the
//     library.
//
// It then writes library/geo/index.json and exits 1 if any example fails.
//
//	go run ./tools/libcheck [-dir library/geo] [-index library/geo/index.json]
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/tools/libcheck/measure"
)

// entry is one record of index.json.
type entry struct {
	File        string           `json:"file"`
	Name        string           `json:"name"`
	Category    string           `json:"category"`
	Difficulty  int              `json:"difficulty"`
	Description string           `json:"description"`
	Tags        []string         `json:"tags"`
	Quality     string           `json:"quality,omitempty"`
	Validation  string           `json:"validation"`
	Cards       int              `json:"cards"`
	HasMirror   bool             `json:"has_mirror,omitempty"`
	Fingerprint string           `json:"fingerprint"`
	Metrics     *measure.Metrics `json:"metrics"`
}

// indexFile is library/geo/index.json.
type indexFile struct {
	Tool       string         `json:"tool"`
	Dir        string         `json:"dir"`
	Tolerance  string         `json:"description_tolerance"`
	Count      int            `json:"count"`
	Categories map[string]int `json:"categories"`
	Validation map[string]int `json:"validation_modes"`
	Quality    map[string]int `json:"quality"`
	SweepCount int            `json:"sweep_examples"`
	Entries    []entry        `json:"entries"`
}

// tagValue returns the value of the first tag with the given prefix.
func tagValue(tags []string, prefix string) (string, bool) {
	for _, t := range tags {
		if strings.HasPrefix(t, prefix) {
			return strings.TrimPrefix(t, prefix), true
		}
	}
	return "", false
}

func asciiOnly(s string) bool {
	for _, r := range s {
		if r > 127 {
			return false
		}
	}
	return true
}

// metadataProblems checks the descriptive fields that the index relies on.
func metadataProblems(base string, p *lens.Prescription) []string {
	var probs []string
	if p.Name == "" {
		probs = append(probs, "missing name")
	}
	if p.Category == "" {
		probs = append(probs, "missing category")
	}
	if p.Description == "" {
		probs = append(probs, "missing description")
	}
	if p.Difficulty < 1 || p.Difficulty > 5 {
		probs = append(probs, fmt.Sprintf("difficulty %d out of 1..5", p.Difficulty))
	}
	if len(p.Tags) < 3 {
		probs = append(probs, "fewer than 3 tags")
	}
	if !asciiOnly(base) {
		probs = append(probs, "non-ASCII file name")
	}
	return probs
}

func main() {
	dir := flag.String("dir", "library/geo", "library directory to validate")
	indexPath := flag.String("index", "library/geo/index.json", "index file to write (\"\" = do not write)")
	quiet := flag.Bool("quiet", false, "print only failures and the summary")
	flag.Parse()

	files, err := filepath.Glob(filepath.Join(*dir, "*.json"))
	if err != nil {
		fmt.Fprintln(os.Stderr, "libcheck:", err)
		os.Exit(1)
	}
	sort.Strings(files)
	var targets []string
	for _, f := range files {
		if filepath.Base(f) != "index.json" {
			targets = append(targets, f)
		}
	}
	if len(targets) == 0 {
		fmt.Fprintf(os.Stderr, "libcheck: no JSON prescriptions in %s\n", *dir)
		os.Exit(1)
	}

	idx := indexFile{
		Tool:       "oz libcheck",
		Dir:        *dir,
		Tolerance:  fmt.Sprintf("description numbers vs recomputed: |got-want| <= %g*|want| + %g", measure.CompareTol, measure.CompareAbs),
		Categories: map[string]int{},
		Validation: map[string]int{},
		Quality:    map[string]int{},
	}
	names := map[string]string{} // name -> file
	fps := map[string]string{}   // parameter fingerprint -> file
	failed := 0
	var bad []string
	reject := func(base string, probs []string) {
		failed++
		bad = append(bad, base)
		fmt.Printf("FAIL %-56s %s\n", base, strings.Join(probs, "; "))
	}

	for _, f := range targets {
		base := filepath.Base(f)
		p, err := lens.LoadJSON(f)
		if err != nil {
			reject(base, []string{"load: " + err.Error()})
			continue
		}
		probs := metadataProblems(base, p)
		if prev, ok := names[p.Name]; ok {
			probs = append(probs, "name duplicates "+prev)
		} else {
			names[p.Name] = base
		}
		fp, err := measure.Fingerprint(p)
		if err != nil {
			probs = append(probs, "fingerprint: "+err.Error())
		} else if prev, ok := fps[fp]; ok {
			probs = append(probs, "parameter fingerprint duplicates "+prev)
		} else {
			fps[fp] = base
		}

		m, err := measure.Measure(p)
		if err != nil {
			reject(base, append(probs, "measure: "+err.Error()))
			continue
		}
		probs = append(probs, measure.Check(p, m)...)

		mode, vals, err := measure.ParseClause(p.Description)
		if err != nil {
			probs = append(probs, err.Error())
		} else {
			probs = append(probs, measure.Compare(m, mode, vals)...)
		}
		if tagMode, ok := tagValue(p.Tags, "mode:"); !ok || tagMode != m.Mode {
			probs = append(probs, fmt.Sprintf("tag mode:%s does not match measured mode %s", tagMode, m.Mode))
		}

		if len(probs) > 0 {
			reject(base, probs)
			continue
		}
		quality, _ := tagValue(p.Tags, "quality:")
		e := entry{
			File: base, Name: p.Name, Category: p.Category, Difficulty: p.Difficulty,
			Description: p.Description, Tags: p.Tags, Quality: quality, Validation: m.Mode,
			Cards: len(p.Cards), HasMirror: measure.HasMirror(p), Fingerprint: fp, Metrics: m,
		}
		idx.Entries = append(idx.Entries, e)
		idx.Categories[p.Category]++
		idx.Validation[m.Mode]++
		if quality != "" {
			idx.Quality[quality]++
		}
		if _, ok := tagValue(p.Tags, "sweep:"); ok {
			idx.SweepCount++
		}
		if !*quiet {
			fmt.Printf("ok   %-56s %-9s %s\n", base, m.Mode, summaryLine(m))
		}
	}

	idx.Count = len(idx.Entries)
	fmt.Printf("\n=== libcheck 汇总 ===\n")
	fmt.Printf("文件 %d 个，通过 %d 个，失败 %d 个，类别 %d 个\n", len(targets), idx.Count, failed, len(idx.Categories))
	cats := make([]string, 0, len(idx.Categories))
	for c := range idx.Categories {
		cats = append(cats, c)
	}
	sort.Strings(cats)
	for _, c := range cats {
		fmt.Printf("  %-18s %4d\n", c, idx.Categories[c])
	}
	modes := make([]string, 0, len(idx.Validation))
	for m := range idx.Validation {
		modes = append(modes, m)
	}
	sort.Strings(modes)
	for _, m := range modes {
		fmt.Printf("校验模式 %-10s %4d\n", m, idx.Validation[m])
	}
	fmt.Printf("参数扫掠示例 %d\n", idx.SweepCount)

	if failed > 0 {
		if *indexPath != "" {
			_ = os.Remove(*indexPath) // a partial index would not match the files
		}
		fmt.Printf("\n失败文件：%s\n", strings.Join(bad, ", "))
		os.Exit(1)
	}
	if *indexPath != "" {
		buf, err := json.MarshalIndent(idx, "", "  ")
		if err != nil {
			fmt.Fprintln(os.Stderr, "libcheck: index:", err)
			os.Exit(1)
		}
		if err := os.WriteFile(*indexPath, append(buf, '\n'), 0o644); err != nil {
			fmt.Fprintln(os.Stderr, "libcheck: index:", err)
			os.Exit(1)
		}
		if err := checkIndexMatches(*indexPath, targets); err != nil {
			fmt.Fprintln(os.Stderr, "libcheck: index:", err)
			os.Exit(1)
		}
		fmt.Printf("索引已写入 %s：%d 个 JSON 文件 ↔ %d 条索引记录，一一对应\n", *indexPath, len(targets), idx.Count)
	}
}

// checkIndexMatches re-reads the index and requires that its file names are
// exactly the JSON files of the library (no missing, extra or duplicate entry).
func checkIndexMatches(indexPath string, targets []string) error {
	buf, err := os.ReadFile(indexPath)
	if err != nil {
		return err
	}
	var back indexFile
	if err := json.Unmarshal(buf, &back); err != nil {
		return err
	}
	want := map[string]bool{}
	for _, f := range targets {
		want[filepath.Base(f)] = true
	}
	seen := map[string]bool{}
	for _, e := range back.Entries {
		if seen[e.File] {
			return fmt.Errorf("index lists %s twice", e.File)
		}
		seen[e.File] = true
		if !want[e.File] {
			return fmt.Errorf("index lists %s, which is not a library file", e.File)
		}
	}
	for f := range want {
		if !seen[f] {
			return fmt.Errorf("library file %s is missing from the index", f)
		}
	}
	if back.Count != len(back.Entries) || len(back.Entries) != len(want) {
		return fmt.Errorf("index count %d, entries %d, files %d", back.Count, len(back.Entries), len(want))
	}
	return nil
}

// summaryLine is the one-line metric summary of a passing example.
func summaryLine(m *measure.Metrics) string {
	switch m.Mode {
	case measure.ModeRefract:
		return fmt.Sprintf("EFL=%.6g F/%.3g RMS(轴)=%.4g RMS(%g%s)=%.4g", m.EFLMM, m.FNumber, m.RMSAxisMM, m.Field, m.FieldUnit, m.RMSFieldMM)
	case measure.ModeFocus:
		return fmt.Sprintf("焦点 z=%.6g RMS=%.4g (%d rays)", m.FocusZMM, m.FocusRMSMM, m.Rays)
	case measure.ModeAfocal:
		return fmt.Sprintf("无焦 出射角 RMS=%.4g µrad M=%.4g (%d rays)", m.SlopeRMSUrad, m.Mag, m.Rays)
	case measure.ModeCollimate:
		return fmt.Sprintf("准直 出射角 RMS=%.4g µrad (%d rays)", m.SlopeRMSUrad, m.Rays)
	}
	return ""
}
