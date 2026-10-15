package main

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"

	"tetsuhiro-oz/lens"
)

// cliDesign writes a 60 mm f/5 achromat prescription into a temp directory and
// returns its path. It is the design every CLI analysis test starts from.
func cliDesign(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "lens.json")
	if err := cmdDesign([]string{"achromat", "--efl", "60", "--semi", "6", "--out", path}); err != nil {
		t.Fatalf("design: %v", err)
	}
	return path
}

// readJSONFile decodes a JSON file written by a --json option.
func readJSONFile(t *testing.T, path string, v any) []byte {
	t.Helper()
	buf, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	if err := json.Unmarshal(buf, v); err != nil {
		t.Fatalf("decode %s: %v", path, err)
	}
	return buf
}

// TestLensThermalCLI drives `oz lens thermal`: the JSON report carries the
// designed focal length and the secant-slope identity, and two runs are
// byte-identical (determinism).
func TestLensThermalCLI(t *testing.T) {
	rx := cliDesign(t)
	dir := t.TempDir()
	run := func(out string) []byte {
		t.Helper()
		if err := cmdLens([]string{"thermal", "--delta-t", "10", "--mount-cte", "1.2e-5", "--json", out, rx}); err != nil {
			t.Fatalf("lens thermal: %v", err)
		}
		var res lensThermalResponse
		return readJSONFile(t, out, &res)
	}
	var res lensThermalResponse
	buf1 := run(filepath.Join(dir, "t1.json"))
	if err := json.Unmarshal(buf1, &res); err != nil {
		t.Fatal(err)
	}
	if math.Abs(res.EFL0MM-60)/60 > 0.01 {
		t.Fatalf("efl0_mm = %g, want the designed 60 mm", res.EFL0MM)
	}
	if res.DeltaTK != 10 || res.MountCTE != 1.2e-5 {
		t.Errorf("echoed ΔT/α = %g/%g, want 10/1.2e-5", res.DeltaTK, res.MountCTE)
	}
	if d := math.Abs(res.EFLPerKMM - (res.EFL1MM-res.EFL0MM)/10); d > 1e-12 {
		t.Errorf("efl_per_k_mm %.17g ≠ (efl1−efl0)/10 %.17g (Δ %.3g)", res.EFLPerKMM, (res.EFL1MM-res.EFL0MM)/10, d)
	}
	if d := math.Abs(res.BFLPerKMM - (res.BFL1MM-res.BFL0MM)/10); d > 1e-12 {
		t.Errorf("bfl_per_k_mm %.17g ≠ (bfl1−bfl0)/10 %.17g (Δ %.3g)", res.BFLPerKMM, (res.BFL1MM-res.BFL0MM)/10, d)
	}
	if !(math.Abs(res.EFLPerKMM) > 0) {
		t.Errorf("efl_per_k_mm = 0, want a non-zero drift")
	}
	if buf2 := run(filepath.Join(dir, "t2.json")); string(buf1) != string(buf2) {
		t.Errorf("thermal CLI is not deterministic:\n%s\n%s", buf1, buf2)
	}
	// Missing ΔT and an unreadable file are errors.
	if err := cmdLens([]string{"thermal", "--delta-t", "0", rx}); err == nil {
		t.Error("thermal with delta-t 0 accepted")
	}
	if err := cmdLens([]string{"thermal", "--delta-t", "10", filepath.Join(dir, "nope.json")}); err == nil {
		t.Error("thermal with a missing prescription accepted")
	}
	if err := cmdLens([]string{"nope"}); err == nil {
		t.Error("unknown lens subcommand accepted")
	}
}

// TestLensSensCLI drives `oz lens sens` in its three modes: the EFL table (with
// a finite-difference cross-check), the inverse allocation and the compensator
// verification.
func TestLensSensCLI(t *testing.T) {
	rx := cliDesign(t)
	dir := t.TempDir()

	// 1. Sensitivity table.
	table := filepath.Join(dir, "sens.json")
	if err := cmdLens([]string{"sens", "--lambda", "587.5618", "--json", table, rx}); err != nil {
		t.Fatalf("lens sens: %v", err)
	}
	var res lensSensResponse
	readJSONFile(t, table, &res)
	if res.Metric != metricEFL || len(res.Rows) != 7 {
		t.Fatalf("table: metric %q with %d rows, want efl with 7", res.Metric, len(res.Rows))
	}
	byName := map[string]sensRow{}
	for _, r := range res.Rows {
		byName[r.Param] = r
		if !r.Converged {
			t.Errorf("%s did not converge: %+v", r.Param, r)
		}
	}
	first := byName["s0.c"]
	if !(first.Deriv < 0) {
		t.Errorf("∂EFL/∂c(s0) = %g, want negative for the front surface of a positive lens", first.Deriv)
	}
	if !(byName["s2.c"].Deriv > 0) {
		t.Errorf("∂EFL/∂c(s2) = %g, want positive for the rear surface", byName["s2.c"].Deriv)
	}
	// Independent central difference of the paraxial EFL.
	p, err := lens.LoadJSON(rx)
	if err != nil {
		t.Fatal(err)
	}
	c0 := 1 / p.Cards[0].RadiusMM
	h := math.Abs(c0) * 1e-5
	p.Cards[0].RadiusMM = 1 / (c0 + h)
	par, err := p.Paraxial(dLineNm)
	if err != nil {
		t.Fatal(err)
	}
	plus := par.EFL
	p.Cards[0].RadiusMM = 1 / (c0 - h)
	if par, err = p.Paraxial(dLineNm); err != nil {
		t.Fatal(err)
	}
	fd := (plus - par.EFL) / (2 * h)
	if math.Abs(first.Deriv-fd) > 1e-4*math.Abs(fd) {
		t.Errorf("CLI ∂EFL/∂c(s0) %.6g vs central difference %.6g (rel %.3g)", first.Deriv, fd, math.Abs(first.Deriv-fd)/math.Abs(fd))
	}

	// 2. Inverse sensitivity.
	inv := filepath.Join(dir, "inv.json")
	if err := cmdLens([]string{"sens", "--inverse", "--budget", "0.005", "--json", inv, rx}); err != nil {
		t.Fatalf("lens sens --inverse: %v", err)
	}
	var ires lensSensResponse
	readJSONFile(t, inv, &ires)
	if ires.Alloc == nil || len(ires.Alloc.Rows) != 7 {
		t.Fatalf("inverse answer incomplete: %+v", ires.Alloc)
	}
	if ires.Alloc.RelErr > 1e-9 {
		t.Errorf("allocation relative error %g, want ≤ 1e-9", ires.Alloc.RelErr)
	}
	if err := cmdLens([]string{"sens", "--inverse", rx}); err == nil {
		t.Error("--inverse without --budget accepted")
	}

	// 3. Compensator verification: the perturbation degrades the merit and the
	// refocus removes at least part of the degradation.
	comp := filepath.Join(dir, "comp.json")
	if err := cmdLens([]string{"sens", "--compensate", "--perturb-radius", "0.05", "--perturb-index", "0.0005", "--json", comp, rx}); err != nil {
		t.Fatalf("lens sens --compensate: %v", err)
	}
	var cres lensSensResponse
	readJSONFile(t, comp, &cres)
	c := cres.Comp
	if c == nil {
		t.Fatal("compensation result missing")
	}
	if !(c.PerturbedMM > c.NominalMM) {
		t.Errorf("perturbed merit %.6g ≤ nominal %.6g: perturbations should degrade the spot", c.PerturbedMM, c.NominalMM)
	}
	if !(c.CompensatedMM <= c.PerturbedMM) {
		t.Errorf("compensated merit %.6g > perturbed %.6g", c.CompensatedMM, c.PerturbedMM)
	}
	if c.Recovery == nil || *c.Recovery <= 0.5 {
		t.Errorf("recovery %v, want a refocus that removes more than half of the degradation", c.Recovery)
	}
	if c.ImageShiftMM == 0 {
		t.Error("image shift compensator did nothing")
	}
	// Determinism of the whole table.
	again := filepath.Join(dir, "sens2.json")
	if err := cmdLens([]string{"sens", "--lambda", "587.5618", "--json", again, rx}); err != nil {
		t.Fatal(err)
	}
	b1, _ := os.ReadFile(table)
	b2, _ := os.ReadFile(again)
	if string(b1) != string(b2) {
		t.Errorf("sens CLI is not deterministic:\n%s\n%s", b1, b2)
	}
}

// TestGeoMTFCLI drives `oz geo mtf`: on-axis tangential and sagittal MTF agree,
// the MTF does not rise with frequency, the field curvature obeys the
// third-order tangential/sagittal/Petzval relation, and two runs agree.
func TestGeoMTFCLI(t *testing.T) {
	rx := cliDesign(t)
	dir := t.TempDir()
	out := filepath.Join(dir, "mtf.json")
	if err := cmdGeo([]string{"mtf", "--fields", "0,1.5", "--cycles-per-mm", "10,20,30,60", "--json", out, rx}); err != nil {
		t.Fatalf("geo mtf: %v", err)
	}
	var res geoMTFResponse
	buf1 := readJSONFile(t, out, &res)
	if len(res.Fields) != 2 || len(res.Foci) != 2 {
		t.Fatalf("mtf answer: %d fields, %d focal rows, want 2 each (%s)", len(res.Fields), len(res.Foci), res.FociNote)
	}
	for i, fld := range res.Fields {
		if len(fld.MTFTangential) != 4 || len(fld.MTFSagittal) != 4 {
			t.Fatalf("field %d: %d/%d MTF values, want 4", i, len(fld.MTFTangential), len(fld.MTFSagittal))
		}
		for k, v := range fld.MTFTangential {
			if v < 0 || v > 1 {
				t.Errorf("field %d: MTF_t[%d] = %g outside [0,1]", i, k, v)
			}
			if k > 0 && v > fld.MTFTangential[k-1]+1e-12 {
				t.Errorf("field %d: MTF_t rises with frequency: %v", i, fld.MTFTangential)
			}
		}
		if fld.RaysHit < fld.RaysTotal*99/100 {
			t.Errorf("field %d: only %d of %d rays reached the image plane", i, fld.RaysHit, fld.RaysTotal)
		}
		if fld.NyquistTCPM < 60 {
			t.Errorf("field %d: Nyquist %.6g below the highest requested 60 cycles/mm", i, fld.NyquistTCPM)
		}
		if i == 0 {
			for k := range fld.MTFTangential {
				if d := math.Abs(fld.MTFTangential[k] - fld.MTFSagittal[k]); d > 1e-3 {
					t.Errorf("on-axis: MTF_t %.6g vs MTF_s %.6g (|Δ| %.3g > 1e-3)", fld.MTFTangential[k], fld.MTFSagittal[k], d)
				}
			}
		}
	}
	for _, fp := range res.Foci {
		if d := math.Abs((fp.SagTMM - fp.SagPMM) - 3*(fp.SagSMM-fp.SagPMM)); d > 1e-6*math.Max(1, math.Abs(fp.SagTMM)) {
			t.Errorf("θ=%g: (T−P) = %.6g vs 3(S−P) = %.6g", fp.ThetaDeg, fp.SagTMM-fp.SagPMM, 3*(fp.SagSMM-fp.SagPMM))
		}
	}
	again := filepath.Join(dir, "mtf2.json")
	if err := cmdGeo([]string{"mtf", "--fields", "0,1.5", "--cycles-per-mm", "10,20,30,60", "--json", again, rx}); err != nil {
		t.Fatal(err)
	}
	buf2, _ := os.ReadFile(again)
	if string(buf1) != string(buf2) {
		t.Errorf("geo mtf CLI is not deterministic:\n%s\n%s", buf1, buf2)
	}
	if err := cmdGeo([]string{"mtf", "--fields", "0,91", rx}); err == nil {
		t.Error("field angle 91° accepted")
	}
	if err := cmdGeo([]string{"nope"}); err == nil {
		t.Error("unknown geo subcommand accepted")
	}
}

// TestGeoGhostCLI drives `oz geo ghost`: the path table is sorted by energy,
// every energy is the product of its two Fresnel reflectances, the classical
// inner-surface pair of the cemented doublet is present, and a coated variant
// replaces the reflectances.
func TestGeoGhostCLI(t *testing.T) {
	rx := cliDesign(t)
	dir := t.TempDir()
	out := filepath.Join(dir, "ghost.json")
	if err := cmdGeo([]string{"ghost", "--json", out, rx}); err != nil {
		t.Fatalf("geo ghost: %v", err)
	}
	var res geoGhostResponse
	buf1 := readJSONFile(t, out, &res)
	if res.Enumerated != 3 || len(res.Paths) != 3 {
		t.Fatalf("ghost answer: %d paths of %d enumerated, want 3/3", len(res.Paths), res.Enumerated)
	}
	prev := math.Inf(1)
	inner := false
	for _, p := range res.Paths {
		if p.Energy > prev {
			t.Errorf("paths are not sorted by energy: %+v", res.Paths)
		}
		prev = p.Energy
		if d := math.Abs(p.Energy - p.ReflectI*p.ReflectJ); d > 1e-15 {
			t.Errorf("%s: energy %.17g ≠ R_I·R_J %.17g (Δ %.3g)", p.Pair, p.Energy, p.ReflectI*p.ReflectJ, d)
		}
		if p.ParaxialXMM == nil || p.ParaxialYMM == nil {
			t.Errorf("%s: no paraxial landing point for a coaxial system", p.Pair)
		}
		if p.Note != "" {
			t.Errorf("%s: unexpected note %q", p.Pair, p.Note)
		}
		if p.I == 2 && p.J == 1 {
			inner = true
		}
	}
	if !inner {
		t.Errorf("the classical (I=2, J=1) inner-surface ghost is missing: %+v", res.Paths)
	}
	again := filepath.Join(dir, "ghost2.json")
	if err := cmdGeo([]string{"ghost", "--json", again, rx}); err != nil {
		t.Fatal(err)
	}
	buf2, _ := os.ReadFile(again)
	if string(buf1) != string(buf2) {
		t.Errorf("geo ghost CLI is not deterministic:\n%s\n%s", buf1, buf2)
	}
	if err := cmdGeo([]string{"ghost"}); err == nil {
		t.Error("ghost without a prescription accepted")
	}
	if err := cmdGeo([]string{"ghost", "--rings", "99", rx}); err == nil {
		t.Error("ghost with 99 sampling rings accepted")
	}
}

// TestDOEFocusCLI drives `oz doe focus`: the first-order efficiency matches the
// analytic η_N = [sin(π/N)/(π/N)]² to 1e-3 at the default grid, the axial peak
// lies within a pixel of the design plane, and two runs agree bit-for-bit.
func TestDOEFocusCLI(t *testing.T) {
	dir := t.TempDir()
	out := filepath.Join(dir, "focus.json")
	if err := cmdDoe([]string{"focus", "--focal-mm", "10", "--levels", "4", "--json", out}); err != nil {
		t.Fatalf("doe focus: %v", err)
	}
	var res doeFocusResponse
	readJSONFile(t, out, &res)
	wantEta := 8 / (math.Pi * math.Pi) // η_4 = [sin(π/4)/(π/4)]²
	if d := math.Abs(res.EtaAnalytic - wantEta); d > 1e-12 {
		t.Errorf("eta_analytic %.17g, want 8/π² = %.17g", res.EtaAnalytic, wantEta)
	}
	if d := math.Abs(res.Eta - res.EtaAnalytic); d > 1e-3 {
		t.Errorf("η %.6f vs η_4 %.6f (|Δ| %.3g > 1e-3)", res.Eta, res.EtaAnalytic, d)
	}
	if math.Abs(res.FocusDeviationPX) > 1 {
		t.Errorf("focus deviation %.4g px, want ≤ 1 px", res.FocusDeviationPX)
	}
	if res.Levels != 4 || res.GridSize != 1024 || res.LambdaNm != 632.8 {
		t.Errorf("echoed parameters %+v, want levels 4, grid 1024, λ 632.8 nm", res)
	}
	// Determinism at a cheap grid.
	small := func(name string) []byte {
		t.Helper()
		p := filepath.Join(dir, name)
		if err := cmdDoe([]string{"focus", "--focal-mm", "10", "--levels", "4", "--grid-size", "128", "--json", p}); err != nil {
			t.Fatal(err)
		}
		b, err := os.ReadFile(p)
		if err != nil {
			t.Fatal(err)
		}
		return b
	}
	if a, b := small("a.json"), small("b.json"); string(a) != string(b) {
		t.Errorf("doe focus CLI is not deterministic:\n%s\n%s", a, b)
	}
	if err := cmdDoe([]string{"focus", "--levels", "4"}); err == nil {
		t.Error("doe focus without --focal-mm accepted")
	}
	if err := cmdDoe([]string{"focus", "--focal-mm", "10", "--levels", "1"}); err == nil {
		t.Error("doe focus with one level accepted")
	}
	if err := cmdDoe([]string{"focus", "--focal-mm", "10", "--levels", "4", "stray.json"}); err == nil {
		t.Error("doe focus with a stray argument accepted")
	}
	if err := cmdDoe([]string{"nope"}); err == nil {
		t.Error("unknown doe subcommand accepted")
	}
}
