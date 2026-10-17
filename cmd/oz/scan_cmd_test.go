package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// scanFixture writes the 1016-style layout (extended source, one slit, one
// screen) to a temp file and returns its path.
func scanFixture(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "scene.json")
	cfg := `{"grid":{"size":64,"width":0.01},"wavelength":6.328e-7,"polarized":false,"method":"asm",
 "sources":[{"id":"src","type":"extended","pos":{"x":0,"y":0,"z":-0.1},
   "params":{"power":0.001,"samples_x":2,"samples_y":1,"standoff":0.04}}],
 "scene":{"components":[
   {"id":"slit","type":"aperture","label":"单缝","pos":{"x":0,"y":0,"z":-0.05},
    "shape":{"kind":"rectangle","params":{"width":0.0004,"height":0.01}}},
   {"id":"det","type":"sensor","label":"屏","pos":{"x":0,"y":0,"z":0.05},"yaw":3.141592653589793,
    "shape":{"kind":"circle","params":{"radius":0.007}}}]}}`
	if err := os.WriteFile(path, []byte(cfg), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// TestScanCLIRunsAndPrints drives `oz scan` end to end and checks the table.
func TestScanCLIRunsAndPrints(t *testing.T) {
	cfgPath := scanFixture(t)
	out := captureStdout(t, func() {
		if err := cmdScan([]string{
			"--param", "scene.components[0].roll",
			"--values", "0,0.5,1.0",
			"--metric", "power",
			cfgPath,
		}); err != nil {
			t.Fatal(err)
		}
	})
	if !strings.Contains(out, "scene.components[0].roll") {
		t.Errorf("output does not name the swept parameter:\n%s", out)
	}
	// Three points, three table rows, and no error lines.
	if n := strings.Count(out, "["); n < 3 {
		t.Errorf("expected 3 table rows, got %d:\n%s", n, out)
	}
	if strings.Contains(out, "错误") {
		t.Errorf("a valid sweep reported errors:\n%s", out)
	}
}

// TestScanCLIWritesJSON checks the --json export.
func TestScanCLIWritesJSON(t *testing.T) {
	cfgPath := scanFixture(t)
	outPath := filepath.Join(t.TempDir(), "scan.json")
	if err := cmdScan([]string{
		"--param", "wavelength", "--values", "5e-7-8e-7-1.5e-7",
		"--metric", "power", "--json", outPath, cfgPath,
	}); err != nil {
		t.Fatal(err)
	}
	buf, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{`"param": "wavelength"`, `"metric": "power"`, `"points"`} {
		if !strings.Contains(string(buf), want) {
			t.Errorf("exported JSON lacks %s:\n%s", want, buf)
		}
	}
	if !strings.Contains(string(buf), "5e-07") && !strings.Contains(string(buf), "5e-7") {
		t.Errorf("exported JSON lacks the swept values:\n%s", buf)
	}
}

// TestScanCLIRejectsBadInput covers the guard rails: missing flags, an
// unusable value list, a missing file and a broken JSON body.
func TestScanCLIRejectsBadInput(t *testing.T) {
	cfgPath := scanFixture(t)
	bad := filepath.Join(t.TempDir(), "bad.json")
	if err := os.WriteFile(bad, []byte("{not json"), 0o644); err != nil {
		t.Fatal(err)
	}
	cases := []struct {
		name string
		args []string
	}{
		{"no param", []string{"--values", "0,1", cfgPath}},
		{"no values", []string{"--param", "wavelength", cfgPath}},
		{"bad values", []string{"--param", "wavelength", "--values", "abc", cfgPath}},
		{"empty values", []string{"--param", "wavelength", "--values", "", cfgPath}},
		{"missing file", []string{"--param", "wavelength", "--values", "6e-7", filepath.Join(t.TempDir(), "nope.json")}},
		{"broken json", []string{"--param", "wavelength", "--values", "6e-7", bad}},
		{"no config", []string{"--param", "wavelength", "--values", "6e-7"}},
	}
	for _, tc := range cases {
		if err := cmdScan(tc.args); err == nil {
			t.Errorf("%s: cmdScan accepted it", tc.name)
		}
	}
}

// TestParseScanValues covers the value-list grammar, including ranges.
func TestParseScanValues(t *testing.T) {
	got, err := parseScanValues("0,15,30")
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 3 || got[0] != 0 || got[1] != 15 || got[2] != 30 {
		t.Errorf("list = %v", got)
	}
	got, err = parseScanValues("0-90-30")
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 4 || got[3] != 90 {
		t.Errorf("range = %v", got)
	}
	got, err = parseScanValues("0-4")
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 5 || got[4] != 4 {
		t.Errorf("implicit step 1 = %v", got)
	}
	got, err = parseScanValues("-1.5e-3,0,1.5e-3")
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 3 || got[0] != -1.5e-3 {
		t.Errorf("signed list = %v", got)
	}
	// A descending range needs an explicit negative step.
	got, err = parseScanValues("90-0--30")
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 4 || got[0] != 90 || got[3] != 0 {
		t.Errorf("descending range = %v", got)
	}
	// Scientific notation must survive the split.
	got, err = parseScanValues("5e-7-8e-7-1.5e-7")
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != 3 || got[0] != 5e-7 || got[2] != 8e-7 {
		t.Errorf("scientific range = %v", got)
	}
	for _, bad := range []string{"", "  ", "abc", "0-90-0", "0-90--5"} {
		if _, err := parseScanValues(bad); err == nil {
			t.Errorf("parseScanValues(%q) accepted it", bad)
		}
	}
}

// captureStdout runs fn with os.Stdout redirected and returns what it wrote.
func captureStdout(t *testing.T, fn func()) string {
	t.Helper()
	r, w, err := os.Pipe()
	if err != nil {
		t.Fatal(err)
	}
	old := os.Stdout
	os.Stdout = w
	done := make(chan string, 1)
	go func() {
		var b strings.Builder
		buf := make([]byte, 4096)
		for {
			n, err := r.Read(buf)
			if n > 0 {
				b.Write(buf[:n])
			}
			if err != nil {
				break
			}
		}
		done <- b.String()
	}()
	fn()
	os.Stdout = old
	_ = w.Close()
	return <-done
}
