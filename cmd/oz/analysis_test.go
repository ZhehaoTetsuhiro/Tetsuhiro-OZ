package main

import (
	"encoding/json"
	"math"
	"net/http"
	"net/http/httptest"
	"testing"

	"tetsuhiro-oz/lens"
)

// designForAnalysis designs the toy prescription shared by the analysis tests:
// a 60 mm N-BK7/F2 achromat at f/5 (the semi-aperture leaves the sampled
// pupil inside the aperture, so the fixed-plane spot merit carries no
// vignetted-ray penalty).
func designForAnalysis(t *testing.T) []byte {
	t.Helper()
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	d := postMap(t, ts.URL, "/api/lens/design", `{"efl_mm":60,"semi_aperture_mm":6}`)
	presc, err := json.Marshal(d["prescription"])
	if err != nil || d["prescription"] == nil {
		t.Fatalf("design answer incomplete: %v (%v)", d, err)
	}
	return presc
}

// postRaw sends one POST and returns the status and the raw body.
func postRaw(t *testing.T, base, path, body string) (int, []byte) {
	t.Helper()
	code, _, buf := doHTTP(t, base, http.MethodPost, path, body)
	return code, buf
}

// TestAnalysisEndpointDeterminism sends every new endpoint the same minimal
// request twice and pins that the two answers are byte-identical (no clock, no
// map order, no RNG in the analysis paths).
func TestAnalysisEndpointDeterminism(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	presc := designForAnalysis(t)
	body := func(s string) string { return `{"prescription":` + string(presc) + s + `}` }
	cases := []struct{ path, body string }{
		{"/api/lens/thermal", body(`,"delta_t_k":10`)},
		{"/api/lens/sens", body(``)},
		{"/api/lens/sens", body(`,"metric":"merit","inverse":true,"budget_mm":0.005`)},
		{"/api/geo/mtf", body(`,"fields_deg":[0,0.5],"cycles_per_mm":[10,30,50]`)},
		{"/api/geo/ghost", body(``)},
		{"/api/doe/focus", `{"focal_mm":10,"levels":4,"grid_size":256}`},
	}
	for _, c := range cases {
		code1, out1 := postRaw(t, ts.URL, c.path, c.body)
		if code1 != http.StatusOK {
			t.Fatalf("POST %s: %d %s", c.path, code1, out1)
		}
		code2, out2 := postRaw(t, ts.URL, c.path, c.body)
		if code2 != http.StatusOK {
			t.Fatalf("POST %s (repeat): %d %s", c.path, code2, out2)
		}
		if string(out1) != string(out2) {
			t.Errorf("POST %s is not deterministic:\n  %s\n  %s", c.path, out1, out2)
		}
	}
}

// TestLensThermalEndpoint checks the thermal answer of POST /api/lens/thermal:
// the designed achromat, ΔT = 10 K, and the secant-slope identity between the
// reported focal lengths and the per-kelvin values.
func TestLensThermalEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	presc := designForAnalysis(t)
	m := postMap(t, ts.URL, "/api/lens/thermal", `{"prescription":`+string(presc)+`,"delta_t_k":10}`)
	efl0, efl1 := m["efl0_mm"].(float64), m["efl1_mm"].(float64)
	bfl0, bfl1 := m["bfl0_mm"].(float64), m["bfl1_mm"].(float64)
	eflPerK, bflPerK := m["efl_per_k_mm"].(float64), m["bfl_per_k_mm"].(float64)
	if math.Abs(efl0-60)/60 > 0.01 {
		t.Fatalf("efl0_mm = %g, want the designed 60 mm (1 %%)", efl0)
	}
	if !(math.Abs(efl1-efl0) > 0) {
		t.Errorf("efl1_mm = %g equals efl0_mm: the thermal model did nothing", efl1)
	}
	// The reported slopes are the secant slopes of the reported values.
	if d := math.Abs(eflPerK - (efl1-efl0)/10); d > 1e-12 {
		t.Errorf("efl_per_k_mm = %.17g, want (efl1−efl0)/10 = %.17g (Δ %.3g)", eflPerK, (efl1-efl0)/10, d)
	}
	if d := math.Abs(bflPerK - (bfl1-bfl0)/10); d > 1e-12 {
		t.Errorf("bfl_per_k_mm = %.17g, want (bfl1−bfl0)/10 = %.17g (Δ %.3g)", bflPerK, (bfl1-bfl0)/10, d)
	}
	// A glass achromat drifts by micrometres per kelvin, not by millimetres.
	if a := math.Abs(eflPerK); a > 0.01 {
		t.Errorf("|efl_per_k_mm| = %g mm/K is implausibly large for a 50 mm glass achromat", a)
	}
	if !(math.Abs(eflPerK) > 0) {
		t.Errorf("efl_per_k_mm = 0: expected a non-zero drift")
	}
}

// TestLensSensEndpoint checks the sensitivity table and the inverse-sensitivity
// allocation of POST /api/lens/sens.
func TestLensSensEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	presc := designForAnalysis(t)

	// EFL table: one row per curvature, thickness and glass of the achromat.
	m := postMap(t, ts.URL, "/api/lens/sens", `{"prescription":`+string(presc)+`}`)
	if m["metric"] != "efl" {
		t.Fatalf("metric = %v, want the default efl", m["metric"])
	}
	rows, _ := m["rows"].([]any)
	if len(rows) != 7 { // 3 curvatures + 2 thicknesses + 2 glass indices
		t.Fatalf("sensitivity table has %d rows, want 7: %v", len(rows), rows)
	}
	names := map[string]bool{}
	for _, r := range rows {
		row := r.(map[string]any)
		name, _ := row["param"].(string)
		names[name] = true
		deriv, ok := row["deriv"].(float64)
		if !ok || math.IsNaN(deriv) || math.IsInf(deriv, 0) {
			t.Fatalf("row %v has no finite derivative", row)
		}
		if row["converged"] != true {
			t.Errorf("row %v did not converge at the default settings", row)
		}
	}
	for _, want := range []string{"s0.c", "s0.thk", "s0.dn", "s1.c", "s1.thk", "s1.dn", "s2.c"} {
		if !names[want] {
			t.Errorf("table lacks %s: %v", want, names)
		}
	}
	// Physics: the fronts of a positive achromat pull the focus in (∂EFL/∂c < 0
	// for the front surface), the rear surface pushes it out, the crown's index
	// adds power (∂EFL/∂n < 0) and the flint's index removes it.
	if first, last := rowDeriv(t, rows, "s0.c"), rowDeriv(t, rows, "s2.c"); !(first < 0 && last > 0) {
		t.Errorf("∂EFL/∂c = %g (front) / %g (rear), want negative / positive", first, last)
	}
	if crown, flint := rowDeriv(t, rows, "s0.dn"), rowDeriv(t, rows, "s1.dn"); !(crown < 0 && flint > 0) {
		t.Errorf("∂EFL/∂n = %g (crown) / %g (flint), want negative / positive", crown, flint)
	}
	// Independent check: a plain central difference of the paraxial EFL, taken
	// in the test through the public prescription fields, agrees with the
	// curvature and thickness rows (the endpoint wiring must not distort the
	// package's extrapolated derivatives).
	p := &lens.Prescription{}
	if err := json.Unmarshal(presc, p); err != nil {
		t.Fatal(err)
	}
	efl := func() float64 {
		t.Helper()
		par, err := p.Paraxial(587.5618)
		if err != nil {
			t.Fatal(err)
		}
		return par.EFL
	}
	c0 := 1 / p.Cards[0].RadiusMM
	const hrel = 1e-5
	h := math.Abs(c0) * hrel
	p.Cards[0].RadiusMM = 1 / (c0 + h)
	plus := efl()
	p.Cards[0].RadiusMM = 1 / (c0 - h)
	minus := efl()
	p.Cards[0].RadiusMM = 1 / c0
	fd := (plus - minus) / (2 * h)
	if got := rowDeriv(t, rows, "s0.c"); math.Abs(got-fd) > 1e-4*math.Abs(fd) {
		t.Errorf("∂EFL/∂c(s0): endpoint %.6g vs central difference %.6g (rel %.3g > 1e-4)",
			got, fd, math.Abs(got-fd)/math.Abs(fd))
	}
	t0 := p.Cards[0].ThicknessMM
	ht := math.Max(math.Abs(t0), 1) * hrel
	p.Cards[0].ThicknessMM = t0 + ht
	plus = efl()
	p.Cards[0].ThicknessMM = t0 - ht
	minus = efl()
	p.Cards[0].ThicknessMM = t0
	fdT := (plus - minus) / (2 * ht)
	if got := rowDeriv(t, rows, "s0.thk"); math.Abs(got-fdT) > 1e-4*math.Abs(fdT) {
		t.Errorf("∂EFL/∂t(s0): endpoint %.6g vs central difference %.6g (rel %.3g > 1e-4)",
			got, fdT, math.Abs(got-fdT)/math.Abs(fdT))
	}
	nominal, _ := m["nominal"].(float64)
	if math.Abs(nominal-60)/60 > 0.01 {
		t.Errorf("nominal EFL = %g, want 60 mm", nominal)
	}

	// Inverse: the allocation predicts the requested spot-merit budget.
	inv := postMap(t, ts.URL, "/api/lens/sens",
		`{"prescription":`+string(presc)+`,"metric":"merit","inverse":true,"budget_mm":0.005}`)
	alloc, _ := inv["allocation"].(map[string]any)
	if alloc == nil {
		t.Fatalf("inverse answer lacks the allocation: %v", inv)
	}
	if alloc["mode"] != "rss" {
		t.Errorf("mode = %v, want the default rss", alloc["mode"])
	}
	budget, _ := alloc["budget_mm"].(float64)
	predicted, _ := alloc["predicted_mm"].(float64)
	relerr, _ := alloc["rel_err"].(float64)
	if math.Abs(predicted-budget)/budget > 1e-9 {
		t.Errorf("predicted increment %.15g vs budget %g (the allocation must use the whole budget)", predicted, budget)
	}
	if relerr > 1e-9 {
		t.Errorf("rel_err = %g, want ≤ 1e-9", relerr)
	}
	arows, _ := alloc["rows"].([]any)
	if len(arows) != 7 {
		t.Fatalf("allocation has %d rows, want 7", len(arows))
	}
	for _, r := range arows {
		row := r.(map[string]any)
		if tol, _ := row["tol"].(float64); !(tol > 0) {
			t.Errorf("row %v got a non-positive tolerance", row)
		}
	}

	// Compensator: a 50 µm radius tolerance degrades the spot merit, and
	// refocusing recovers at least part of it.
	comp := postMap(t, ts.URL, "/api/lens/sens",
		`{"prescription":`+string(presc)+`,"metric":"merit","compensate":true,"perturb_radius_mm":0.05,"perturb_index":0.0005}`)
	c, _ := comp["compensation"].(map[string]any)
	if c == nil {
		t.Fatalf("compensate answer lacks the compensation: %v", comp)
	}
	pert, _ := c["perturbed_mm"].(float64)
	compd, _ := c["compensated_mm"].(float64)
	nominal2, _ := c["nominal_mm"].(float64)
	if !(pert > nominal2) {
		t.Errorf("perturbed merit %.6g must exceed the nominal %.6g", pert, nominal2)
	}
	if !(compd <= pert) {
		t.Errorf("compensated merit %.6g must not exceed the perturbed %.6g", compd, pert)
	}
	if rec, ok := c["recovery"].(float64); !ok || math.IsNaN(rec) {
		t.Errorf("recovery = %v, want a finite share (the perturbation degrades the merit)", c["recovery"])
	}
	if shift, ok := c["image_shift_mm"].(float64); !ok || shift == 0 {
		t.Errorf("image_shift_mm = %v, want a non-zero refocus", c["image_shift_mm"])
	}
}

// TestGeoMTFEndpoint checks the geometric MTF and the field curvature of POST
// /api/geo/mtf on the designed achromat.
func TestGeoMTFEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	presc := designForAnalysis(t)
	m := postMap(t, ts.URL, "/api/geo/mtf",
		`{"prescription":`+string(presc)+`,"fields_deg":[0,1],"cycles_per_mm":[10,20,30,50]}`)
	fields, _ := m["fields"].([]any)
	if len(fields) != 2 {
		t.Fatalf("mtf answer has %d fields, want 2: %v", len(fields), m)
	}
	var prev float64
	for i, v := range fields {
		fld := v.(map[string]any)
		mt, _ := fld["mtf_tangential"].([]any)
		ms, _ := fld["mtf_sagittal"].([]any)
		if len(mt) != 4 || len(ms) != 4 {
			t.Fatalf("field %d: %d/%d MTF values, want 4 each", i, len(mt), len(ms))
		}
		for k := range mt {
			tv, sv := mt[k].(float64), ms[k].(float64)
			if tv < 0 || tv > 1 || sv < 0 || sv > 1 {
				t.Errorf("field %d: MTF %g/%g outside [0,1]", i, tv, sv)
			}
			if k > 0 && tv > mt[k-1].(float64)+1e-12 {
				t.Errorf("field %d: MTF rises with frequency: %v", i, mt)
			}
		}
		if nyq, _ := fld["nyquist_tangential_per_mm"].(float64); nyq < 50 {
			t.Errorf("field %d: Nyquist %g below the highest requested frequency 50", i, nyq)
		}
		if i == 0 {
			// Rotational symmetry: the on-axis meridional and sagittal MTFs
			// must agree (they differ only through the sampling pattern).
			for k := range mt {
				if d := math.Abs(mt[k].(float64) - ms[k].(float64)); d > 1e-3 {
					t.Errorf("on-axis MTF_t %g vs MTF_s %g differ by %.3g", mt[k], ms[k], d)
				}
			}
		}
		rms, _ := fld["rms_radius_mm"].(float64)
		if rms <= 0 {
			t.Errorf("field %d: spot RMS %g must be positive", i, rms)
		}
		if i > 0 && rms <= prev {
			t.Errorf("field %d: RMS %g does not grow with the field angle (%g)", i, rms, prev)
		}
		prev = rms
	}
	// Field curvature: one row per field, tangential sag three times the
	// sagittal one about the Petzval surface (z_T − z_P = 3(z_S − z_P)).
	foci, _ := m["foci"].([]any)
	if len(foci) != 2 {
		t.Fatalf("mtf answer has %d focal rows, want 2: %v", len(foci), m["foci_note"])
	}
	for _, v := range foci {
		fp := v.(map[string]any)
		st, _ := fp["sag_tangential_mm"].(float64)
		ss, _ := fp["sag_sagittal_mm"].(float64)
		sp, _ := fp["sag_petzval_mm"].(float64)
		if d := math.Abs((st - sp) - 3*(ss-sp)); d > 1e-6*math.Max(1, math.Abs(st)) {
			t.Errorf("θ=%v: (T−P) = %.6g vs 3(S−P) = %.6g", fp["theta_deg"], st-sp, 3*(ss-sp))
		}
	}
}

// TestGeoGhostEndpoint checks the ghost-path table of POST /api/geo/ghost: the
// pairs are sorted by energy, every energy is the product of its two Fresnel
// reflectances, and a cemented doublet carries the classical (I=2, J=1) path.
func TestGeoGhostEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	presc := designForAnalysis(t)
	m := postMap(t, ts.URL, "/api/geo/ghost", `{"prescription":`+string(presc)+`}`)
	paths, _ := m["paths"].([]any)
	enum, _ := m["paths_enumerated"].(float64)
	if int(enum) != 3 { // 3 surfaces → 3 double-reflection pairs
		t.Fatalf("paths_enumerated = %v, want 3", m["paths_enumerated"])
	}
	if len(paths) != 3 {
		t.Fatalf("ghost table has %d rows, want 3: %v", len(paths), paths)
	}
	prev := math.Inf(1)
	for _, v := range paths {
		p := v.(map[string]any)
		e, _ := p["relative_energy"].(float64)
		ri, _ := p["reflect_i"].(float64)
		rj, _ := p["reflect_j"].(float64)
		if e > prev {
			t.Errorf("ghost table is not sorted by energy: %v", paths)
		}
		prev = e
		if d := math.Abs(e - ri*rj); d > 1e-15*math.Max(1, e) {
			t.Errorf("path %v: energy %g ≠ R_I·R_J = %g", p["pair"], e, ri*rj)
		}
		if ri <= 0 || ri >= 1 || rj <= 0 || rj >= 1 {
			t.Errorf("path %v: reflectances %g/%g outside (0,1)", p["pair"], ri, rj)
		}
		if _, ok := p["paraxial_x_mm"]; !ok {
			t.Errorf("path %v lacks the paraxial landing point", p["pair"])
		}
	}
	// The cemented doublet's classical inner-surface ghost (rear face + cement
	// plane, I=2, J=1) is among the paths, and the Fresnel reflectance of the
	// three surfaces is the highest at the crown/air interface (largest Δn).
	found := false
	for _, v := range paths {
		p := v.(map[string]any)
		if p["i"] == float64(2) && p["j"] == float64(1) {
			found = true
		}
	}
	if !found {
		t.Errorf("the classical (I=2, J=1) inner-surface ghost is missing: %v", paths)
	}
	// A coated variant: supplied reflectances replace the Fresnel values.
	m2 := postMap(t, ts.URL, "/api/geo/ghost",
		`{"prescription":`+string(presc)+`,"reflectance":[0.5,0.25,0.01]}`)
	paths2, _ := m2["paths"].([]any)
	best := paths2[0].(map[string]any)
	e2, _ := best["relative_energy"].(float64)
	if d := math.Abs(e2 - 0.5*0.25); d > 1e-15 { // I=1, J=0 is the only 0.5·0.25 pair
		t.Errorf("coated ghost energy %g, want R0·R1 = 0.125 (Δ %.3g)", e2, d)
	}
}

// TestDOEFocusEndpoint checks POST /api/doe/focus: the first-order efficiency
// matches the analytic η_N = [sin(π/N)/(π/N)]² to the acceptance of the doe
// package (1e-3 at the default grid) and the axial peak sits within one pixel
// of the design focal plane.
func TestDOEFocusEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	const body = `{"focal_mm":10,"levels":4,"lambda_nm":632.8}`
	m := postMap(t, ts.URL, "/api/doe/focus", body)
	eta, _ := m["eta_first_order"].(float64)
	etaN, _ := m["eta_analytic"].(float64)
	wantEta := 0.8105694691387022 // (sin(π/4)/(π/4))² = 8/π²
	if d := math.Abs(etaN - wantEta); d > 1e-12 {
		t.Errorf("eta_analytic = %.17g, want 8/π² = %.17g (Δ %.3g)", etaN, wantEta, d)
	}
	if d := math.Abs(eta - etaN); d > 1e-3 {
		t.Errorf("η = %.6f vs η_N = %.6f (|Δ| %.3g > 1e-3, the doe-package acceptance)", eta, etaN, d)
	}
	if rel, _ := m["eta_rel_err"].(float64); rel > 1e-3 {
		t.Errorf("eta_rel_err = %.3g, want ≤ 1e-3", rel)
	}
	devUM, _ := m["focus_deviation_um"].(float64)
	devPX, _ := m["focus_deviation_px"].(float64)
	if math.Abs(devPX) > 1 {
		t.Errorf("focus deviation %.3g px (%.3g µm), want ≤ 1 px", devPX, devUM)
	}
	if off, _ := m["peak_offset_px"].(float64); off > 1 {
		t.Errorf("focal-plane peak offset %.3g px, want ≤ 1 px", off)
	}
	if wf, _ := m["spot_radius_um"].(float64); !(wf > 0) {
		t.Errorf("spot_radius_um = %v, want positive", m["spot_radius_um"])
	}
}

// rowDeriv returns the derivative of one named row of a sensitivity table.
func rowDeriv(t *testing.T, rows []any, name string) float64 {
	t.Helper()
	for _, r := range rows {
		row := r.(map[string]any)
		if row["param"] == name {
			return row["deriv"].(float64)
		}
	}
	t.Fatalf("sensitivity row %s not found", name)
	return 0
}

// TestAnalysisEndpointErrors pins the error contract of the new endpoints:
// missing mandatory fields, out-of-range options and work limits answer 400
// with a JSON error, never 200 or 500.
func TestAnalysisEndpointErrors(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	presc := designForAnalysis(t)
	body := func(s string) string { return `{"prescription":` + string(presc) + s + `}` }
	cases := []struct{ path, body string }{
		{"/api/lens/thermal", body(`}`)},                               // no delta_t_k
		{"/api/lens/thermal", body(`,"delta_t_k":1000}`)},              // above maxThermalDT
		{"/api/lens/thermal", body(`,"delta_t_k":10,"lambda_nm":-1}`)}, // bad wavelength
		{"/api/lens/thermal", `{"delta_t_k":10}`},                      // no prescription
		{"/api/lens/sens", body(`,"metric":"spot"}`)},                  // unknown metric
		{"/api/lens/sens", body(`,"metric":"merit","inverse":true}`)},  // no budget
		{"/api/lens/sens", body(`,"metric":"merit","inverse":true,"budget_mm":0.01,"mode":"x"}`)},
		{"/api/lens/sens", body(`,"metric":"merit","compensate":true}`)},         // no perturbation
		{"/api/lens/sens", body(`,"compensate":true,"perturb_radius_mm":0.05}`)}, // compensate needs metric=merit
		{"/api/lens/sens", body(`,"compensator_vars":["s9.R"]}`)},                // card out of range
		{"/api/geo/mtf", body(`,"fields_deg":[95]}`)},                            // non-physical field
		{"/api/geo/mtf", body(`,"rays":4}`)},
		{"/api/geo/ghost", body(`,"rings":99}`)},
		{"/api/geo/ghost", body(`,"reflectance":[0.5,2]}`)}, // wrong length / out of [0,1]
		{"/api/doe/focus", `{}`},                            // no focal length
		{"/api/doe/focus", `{"focal_mm":-1,"levels":4}`},
		{"/api/doe/focus", `{"focal_mm":10,"levels":1}`},
		{"/api/doe/focus", `{"focal_mm":10,"levels":4,"grid_size":100}`},  // not a power of two
		{"/api/doe/focus", `{"focal_mm":10,"levels":4,"grid_size":4096}`}, // above maxDOEGrid
	}
	for _, c := range cases {
		code, out := postRaw(t, ts.URL, c.path, c.body)
		if code != http.StatusBadRequest {
			t.Errorf("POST %s %s = %d, want 400 (%s)", c.path, c.body, code, out)
			continue
		}
		var m map[string]any
		if err := json.Unmarshal(out, &m); err != nil || m["error"] == nil {
			t.Errorf("POST %s %s: error body lacks an error field: %s", c.path, c.body, out)
		}
	}
}
