package main

import (
	"bytes"
	"encoding/json"
	"errors"
	"flag"
	"io"
	"math"
	"net/http"
	"net/http/httptest"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"tetsuhiro-oz/optics"
)

// serviceSimConfig is a small sequential simulation used to drive POST
// /api/simulate, /api/validate and /api/convert through the unified mux.
const serviceSimConfig = `{"grid":{"size":64,"width":0.004},"wavelength":6.328e-7,"method":"asm",
 "source":{"type":"gaussian","params":{"waist":0.0005}},
 "elements":[{"type":"propagate","params":{"distance":0.1}},
             {"type":"lens","params":{"f":0.2}},
             {"type":"sensor","params":{"label":"focus"}}]}`

// serviceQuantumConfig is a Hong-Ou-Mandel setup (|1,1⟩ through a 50:50 beam
// splitter).
const serviceQuantumConfig = `{"modes":2,"cutoff":4,"state":{"type":"fock","params":{"occupation":[1,1]}},` +
	`"gates":[{"type":"beam_splitter","params":{"mode0":0,"mode1":1,"reflectivity":0.5}}]}`

// testElementDef is a scripted element definition (a metalens phase profile).
const testElementDef = `{
  "name": "metalens",
  "label": "超表面透镜",
  "behavior": "transmit",
  "class": "lens",
  "params": [{"key": "f", "label": "焦距", "unit": "m", "kind": "float", "default": 0.3}],
  "phase": "-k*(sqrt(r*r+f*f) - f)"
}`

// doHTTP performs one request and returns status, headers and body.
func doHTTP(t *testing.T, base, method, path, body string) (int, http.Header, []byte) {
	t.Helper()
	var rd io.Reader
	if body != "" {
		rd = strings.NewReader(body)
	}
	req, err := http.NewRequest(method, base+path, rd)
	if err != nil {
		t.Fatal(err)
	}
	if body != "" {
		req.Header.Set("Content-Type", "application/json")
	}
	// Redirects are answers in their own right here: report the 3xx itself.
	noFollow := &http.Client{CheckRedirect: func(*http.Request, []*http.Request) error {
		return http.ErrUseLastResponse
	}}
	res, err := noFollow.Do(req)
	if err != nil {
		t.Fatal(err)
	}
	defer res.Body.Close()
	buf, err := io.ReadAll(res.Body)
	if err != nil {
		t.Fatal(err)
	}
	return res.StatusCode, res.Header, buf
}

// postMap posts a JSON body and decodes a 200 JSON object answer.
func postMap(t *testing.T, base, path, body string) map[string]any {
	t.Helper()
	code, _, buf := doHTTP(t, base, http.MethodPost, path, body)
	if code != http.StatusOK {
		t.Fatalf("POST %s: %d %s", path, code, buf)
	}
	var m map[string]any
	if err := json.Unmarshal(buf, &m); err != nil {
		t.Fatalf("POST %s: not JSON: %v", path, err)
	}
	return m
}

// TestServiceMountsEveryEndpoint sends the smallest request each route accepts
// through the one mux and checks the answer's status. Every path must be
// routed: a 404 here would mean a mount is missing.
func TestServiceMountsEveryEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	cases := []struct {
		method, path, body string
		want               int
	}{
		// WOS simulation API.
		{http.MethodGet, "/api/catalog", "", http.StatusOK},
		{http.MethodGet, "/api/health", "", http.StatusOK},
		{http.MethodGet, "/api/elements", "", http.StatusOK},
		{http.MethodPost, "/api/elements/reload", "", http.StatusOK},
		{http.MethodPost, "/api/validate", serviceSimConfig, http.StatusOK},
		{http.MethodPost, "/api/quantum", serviceQuantumConfig, http.StatusOK},
		{http.MethodPost, "/api/convert", serviceSimConfig, http.StatusOK},
		{http.MethodPost, "/api/simulate", serviceSimConfig, http.StatusAccepted},
		{http.MethodGet, "/api/runs/no-such-run", "", http.StatusNotFound},
		// OZ geometry and coating API. eval, tol and the spectrum take a real
		// prescription or stack, so they are exercised in TestOZEndpointChain.
		// The analysis endpoints below take one too; an empty body is routed
		// (400, not 404) and proves the mount.
		{http.MethodGet, "/api/version", "", http.StatusOK},
		{http.MethodPost, "/api/lens/design", `{"efl_mm":60,"semi_aperture_mm":6}`, http.StatusOK},
		{http.MethodPost, "/api/lens/thermal", `{}`, http.StatusBadRequest},
		{http.MethodPost, "/api/lens/sens", `{}`, http.StatusBadRequest},
		{http.MethodPost, "/api/geo/mtf", `{}`, http.StatusBadRequest},
		{http.MethodPost, "/api/geo/ghost", `{}`, http.StatusBadRequest},
		{http.MethodPost, "/api/doe/focus", `{"focal_mm":10,"levels":4,"grid_size":128}`, http.StatusOK},
		{http.MethodPost, "/api/coating/design", `{"template":"single-ar","lambda_nm":550}`, http.StatusOK},
		// GUI and the pre-merge address.
		{http.MethodGet, "/", "", http.StatusOK},
		{http.MethodGet, "/oz/", "", http.StatusMovedPermanently},
		{http.MethodGet, "/app.js", "", http.StatusOK},
		{http.MethodGet, "/style.css", "", http.StatusOK},
		{http.MethodGet, "/favicon.svg", "", http.StatusOK},
		{http.MethodGet, "/fonts/rajdhani-700.woff2", "", http.StatusOK},
	}
	for _, c := range cases {
		code, _, body := doHTTP(t, ts.URL, c.method, c.path, c.body)
		if code != c.want {
			t.Errorf("%s %s = %d, want %d (%s)", c.method, c.path, code, c.want, bytes.TrimSpace(body))
		}
	}
}

// TestServiceMethodContract checks that each endpoint answers 405 to a method
// it does not serve, in both the OZ and the WOS half of the mux.
func TestServiceMethodContract(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	for _, c := range []struct{ method, path string }{
		{http.MethodPost, "/api/version"},
		{http.MethodGet, "/api/lens/design"},
		{http.MethodGet, "/api/lens/eval"},
		{http.MethodGet, "/api/lens/tol"},
		{http.MethodGet, "/api/lens/thermal"},
		{http.MethodGet, "/api/lens/sens"},
		{http.MethodGet, "/api/geo/mtf"},
		{http.MethodGet, "/api/geo/ghost"},
		{http.MethodGet, "/api/doe/focus"},
		{http.MethodGet, "/api/coating/design"},
		{http.MethodGet, "/api/coating/spectrum"},
		{http.MethodGet, "/api/simulate"},
		{http.MethodPost, "/api/catalog"},
		{http.MethodGet, "/api/elements/reload"},
		{http.MethodGet, "/api/quantum"},
	} {
		code, _, _ := doHTTP(t, ts.URL, c.method, c.path, "")
		if code != http.StatusMethodNotAllowed {
			t.Errorf("%s %s = %d, want 405", c.method, c.path, code)
		}
	}
}

// TestServiceOwnership checks that the exact OZ paths and the /api/ prefix do
// not shadow each other: /api/version is answered by OZ, /api/health by WOS.
func TestServiceOwnership(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	_, _, v := doHTTP(t, ts.URL, http.MethodGet, "/api/version", "")
	if !strings.Contains(string(v), `"name":"tetsuhiro-oz"`) {
		t.Errorf("/api/version not answered by OZ: %s", v)
	}
	_, _, h := doHTTP(t, ts.URL, http.MethodGet, "/api/health", "")
	if !strings.Contains(string(h), `"gpu"`) {
		t.Errorf("/api/health not answered by WOS: %s", h)
	}
}

// TestOZEndpointChain drives the geometry and coating endpoints end to end:
// design a lens, evaluate and tolerance-analyse the result, design a coating,
// scan its spectrum. Each step feeds the next.
func TestOZEndpointChain(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	base := ts.URL

	d := postMap(t, base, "/api/lens/design", `{"efl_mm":60,"semi_aperture_mm":6,"field_deg":0.5}`)
	presc, err := json.Marshal(d["prescription"])
	if err != nil || d["prescription"] == nil || d["optimizer"] == nil {
		t.Fatalf("design answer incomplete: %v", d)
	}

	// Wiring check: the designed focal length is the requested one.
	e := postMap(t, base, "/api/lens/eval", string(presc))
	metrics, _ := e["metrics"].(map[string]any)
	efl, ok := metrics["efl_mm"].(float64)
	if !ok {
		t.Fatalf("eval lacks efl_mm: %v", e)
	}
	if math.Abs(efl-60)/60 > 0.01 {
		t.Errorf("designed EFL = %g mm, requested 60 mm (more than 1 %% off)", efl)
	}

	tolBody, _ := json.Marshal(map[string]any{
		"prescription": d["prescription"],
		"spec":         map[string]any{"radius_mm": 0.05, "thickness_mm": 0.05, "index_delta": 0.0005},
		"criteria":     map[string]any{"target_efl_mm": 60, "max_focal_err_mm": 0.1},
		"runs":         20,
		"seed":         3,
	})
	tol := postMap(t, base, "/api/lens/tol", string(tolBody))
	if len(tol) == 0 {
		t.Fatal("tolerance report is empty")
	}

	// The analysis endpoints accept the designed prescription as well: design →
	// thermal / sensitivity / MTF / ghost in one chain.
	for _, c := range []struct{ path, tail string }{
		{"/api/lens/thermal", `,"delta_t_k":10`},
		{"/api/lens/sens", `,"metric":"merit","inverse":true,"budget_mm":0.005`},
		{"/api/geo/mtf", `,"fields_deg":[0,0.5],"cycles_per_mm":[10,30]`},
		{"/api/geo/ghost", ``},
	} {
		m := postMap(t, base, c.path, `{"prescription":`+string(presc)+c.tail+`}`)
		if len(m) == 0 {
			t.Errorf("%s answered an empty object", c.path)
		}
	}
	if m := postMap(t, base, "/api/doe/focus", `{"focal_mm":10,"levels":4,"grid_size":128}`); m["eta_first_order"] == nil {
		t.Errorf("doe focus answer incomplete: %v", m)
	}

	c := postMap(t, base, "/api/coating/design", `{"template":"single-ar","lambda_nm":550,"substrate":"BK7"}`)
	stack, err := json.Marshal(c["stack"])
	if err != nil || c["stack"] == nil {
		t.Fatalf("coating design lacks stack: %v", c)
	}
	s := postMap(t, base, "/api/coating/spectrum",
		`{"stack":`+string(stack)+`,"from_nm":500,"to_nm":600,"step_nm":10,"pol":"u"}`)
	rows, _ := s["rows"].([]any)
	if len(rows) != 11 {
		t.Fatalf("spectrum has %d rows, want 11 (500..600 nm step 10)", len(rows))
	}
	// Physics check: a single-layer AR must cut the bare-glass reflectance of
	// BK7 at 550 nm (Fresnel, n = 1.5185: R0 = ((n-1)/(n+1))² ≈ 4.24 %) by more
	// than half.
	var r550 float64
	for _, r := range rows {
		m := r.(map[string]any)
		if m["lambda_nm"] == float64(550) {
			r550 = m["R"].(float64)
		}
	}
	if !(r550 < 0.0212) {
		t.Errorf("R(550 nm) = %.4f, want < 0.0212 (half of the bare-glass Fresnel value)", r550)
	}

	// Defaults: an empty spectrum request scans 400–700 nm in 5 nm steps.
	s2 := postMap(t, base, "/api/coating/spectrum", `{"stack":`+string(stack)+`}`)
	if rows2, _ := s2["rows"].([]any); len(rows2) != 61 {
		t.Errorf("default spectrum has %d rows, want 61", len(rows2))
	}
}

// TestOZEndpointErrors checks the OZ error answers: missing or invalid input
// gets a 400 with a JSON error, never a 200 or a crash.
func TestOZEndpointErrors(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	base := ts.URL
	for _, c := range []struct {
		method, path, body string
		want               int
	}{
		{http.MethodPost, "/api/lens/design", `{}`, http.StatusBadRequest},
		{http.MethodPost, "/api/lens/design", `{"efl_mm":-5}`, http.StatusBadRequest},
		{http.MethodPost, "/api/lens/design", `{"efl_mm":"x"}`, http.StatusBadRequest},
		{http.MethodPost, "/api/lens/eval", `{"cards":`, http.StatusBadRequest},
		{http.MethodPost, "/api/coating/design", `{"template":"nope"}`, http.StatusBadRequest},
		{http.MethodPost, "/api/coating/design", `{"template":"single-ar","incident":"Unobtainium"}`, http.StatusBadRequest},
		{http.MethodPost, "/api/coating/design", `{"template":"lwp"}`, http.StatusBadRequest},
		{http.MethodPost, "/api/coating/spectrum", `{"stack":{}, "pol":"x"}`, http.StatusBadRequest},
		// Work limits: refused before any computation (see the maxSpectrumPoints,
		// maxTolRuns and maxDesignIter constants).
		{http.MethodPost, "/api/coating/spectrum", `{"stack":{}, "from_nm":400, "to_nm":700, "step_nm":0.001}`, http.StatusBadRequest},
		{http.MethodPost, "/api/lens/tol", `{"prescription":{}, "runs":20001}`, http.StatusBadRequest},
		{http.MethodPost, "/api/lens/design", `{"efl_mm":100, "max_iter":10001}`, http.StatusBadRequest},
	} {
		code, hdr, body := doHTTP(t, base, c.method, c.path, c.body)
		if code != c.want {
			t.Errorf("%s %s %s = %d, want %d (%s)", c.method, c.path, c.body, code, c.want, body)
			continue
		}
		if !strings.HasPrefix(hdr.Get("Content-Type"), "application/json") {
			t.Errorf("%s %s error is not JSON: %q", c.method, c.path, hdr.Get("Content-Type"))
		}
		var m map[string]any
		if err := json.Unmarshal(body, &m); err != nil || m["error"] == nil {
			t.Errorf("%s %s error body lacks an error field: %s", c.method, c.path, body)
		}
	}
}

// TestGUIServedAtRoot checks the embedded GUI: the page, its assets and the
// MIME types the browser needs (font/woff2 for the vendored font).
func TestGUIServedAtRoot(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	code, hdr, page := doHTTP(t, ts.URL, http.MethodGet, "/", "")
	if code != http.StatusOK || !strings.HasPrefix(hdr.Get("Content-Type"), "text/html") {
		t.Fatalf("GET / = %d %s", code, hdr.Get("Content-Type"))
	}
	if !strings.Contains(string(page), `<meta name="viewport" content="width=device-width, initial-scale=1">`) {
		t.Error("GET / lacks the viewport meta tag")
	}
	if cc := hdr.Get("Cache-Control"); cc != "no-cache" {
		t.Errorf("GET / Cache-Control = %q, want no-cache (GUI must revalidate after an upgrade)", cc)
	}
	for path, ctype := range map[string]string{
		"/app.js":                   "javascript",
		"/oz.js":                    "javascript",
		"/style.css":                "text/css",
		"/favicon.svg":              "image/svg+xml",
		"/fonts/rajdhani-700.woff2": "font/woff2",
	} {
		code, hdr, _ := doHTTP(t, ts.URL, http.MethodGet, path, "")
		if code != http.StatusOK || !strings.Contains(hdr.Get("Content-Type"), ctype) {
			t.Errorf("GET %s = %d %q, want %s", path, code, hdr.Get("Content-Type"), ctype)
		}
	}
	if code, _, _ := doHTTP(t, ts.URL, http.MethodGet, "/no-such-file", ""); code != http.StatusNotFound {
		t.Errorf("GET /no-such-file = %d, want 404", code)
	}
}

// TestParseServeFlagsAcceptsWOSOptions checks that `oz serve` takes every
// option the retired wos binary took, with the same defaults.
func TestParseServeFlagsAcceptsWOSOptions(t *testing.T) {
	o, err := parseServeFlags(nil)
	if err != nil {
		t.Fatal(err)
	}
	if o.Addr != ":1130" || o.MaxRunMB != 512 || o.GPU || o.Check || o.NoWatch || o.GenGo != "" {
		t.Fatalf("defaults = %+v", o)
	}
	o, err = parseServeFlags([]string{
		"-addr", ":7007", "-max-run-mb", "64", "-elements", "a,b",
		"-no-elements-watch", "-gpu", "-check-elements",
		"-gen-go", "x.json", "-gen-name", "Y",
	})
	if err != nil {
		t.Fatal(err)
	}
	if o.Addr != ":7007" || o.MaxRunMB != 64 || o.ExtraDirs != "a,b" || !o.NoWatch || !o.GPU || !o.Check ||
		o.GenGo != "x.json" || o.GenName != "Y" {
		t.Fatalf("parsed = %+v", o)
	}
	if _, err := parseServeFlags([]string{"-no-such-flag"}); err == nil {
		t.Error("unknown flag accepted")
	}
	if _, err := parseServeFlags([]string{"stray"}); err == nil {
		t.Error("positional argument accepted")
	}
	if _, err := parseServeFlags([]string{"-h"}); !errors.Is(err, flag.ErrHelp) {
		t.Errorf("-h returned %v, want flag.ErrHelp", err)
	}
}

// restoreElements returns the scripted-element registry to the defaults after a
// test that changed the definition directories.
func restoreElements(t *testing.T) {
	t.Helper()
	t.Cleanup(func() {
		optics.SetElementDirs(optics.DefaultElementDirs())
		optics.ReloadElementDefinitions(nil)
	})
}

// TestServeCheckElements runs the -check-elements path of `oz serve` on a good
// and a broken definition directory.
func TestServeCheckElements(t *testing.T) {
	restoreElements(t)
	good := t.TempDir()
	if err := os.WriteFile(filepath.Join(good, "metalens.json"), []byte(testElementDef), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := cmdServe([]string{"-check-elements", "-no-elements-watch", "-elements", good}); err != nil {
		t.Fatalf("valid definitions rejected: %v", err)
	}
	bad := t.TempDir()
	if err := os.WriteFile(filepath.Join(bad, "broken.json"), []byte(`{"name": "broken", "phase": "1 +"}`), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := cmdServe([]string{"-check-elements", "-no-elements-watch", "-elements", bad}); err == nil {
		t.Error("broken definition passed -check-elements")
	}
}

// TestRunGenGoWritesNativeElement checks the -gen-go exporter: the output is a
// Go source file for the optics package that names the requested type.
func TestRunGenGoWritesNativeElement(t *testing.T) {
	dir := t.TempDir()
	def := filepath.Join(dir, "metalens.json")
	if err := os.WriteFile(def, []byte(testElementDef), 0o644); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(dir, "out.go.txt")
	f, err := os.Create(out)
	if err != nil {
		t.Fatal(err)
	}
	saved := os.Stdout
	os.Stdout = f
	err = runGenGo(def, "metalens_native")
	os.Stdout = saved
	f.Close()
	if err != nil {
		t.Fatal(err)
	}
	src, err := os.ReadFile(out)
	if err != nil {
		t.Fatal(err)
	}
	for _, want := range []string{"package optics", "metalens_native"} {
		if !strings.Contains(string(src), want) {
			t.Errorf("generated source lacks %q", want)
		}
	}
}

// TestOZRedirectsToGUITab checks that the pre-merge /oz/ address lands on the
// geometry/coating tab of the GUI.
func TestOZRedirectsToGUITab(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	code, hdr, _ := doHTTP(t, ts.URL, http.MethodGet, "/oz/", "")
	if code != http.StatusMovedPermanently || hdr.Get("Location") != "/#oz" {
		t.Errorf("GET /oz/ = %d Location %q, want 301 /#oz", code, hdr.Get("Location"))
	}
}

// TestGUIScriptsParse runs `node --check` on the two GUI scripts. The binary is
// taken from $OZ_NODE, else from PATH; without node the test is skipped and the
// check has to be run by hand.
func TestGUIScriptsParse(t *testing.T) {
	node := os.Getenv("OZ_NODE")
	if node == "" {
		found, err := exec.LookPath("node")
		if err != nil {
			t.Skip("node not installed (set OZ_NODE); syntax check skipped")
		}
		node = found
	}
	for _, f := range []string{"app.js", "oz.js"} {
		out, err := exec.Command(node, "--check", filepath.Join("..", "..", "server", "web", f)).CombinedOutput()
		if err != nil {
			t.Errorf("node --check %s: %v\n%s", f, err, out)
		}
	}
}

// TestSpectrumPointLimit checks the boundary of the spectrum work limit: exactly
// maxSpectrumPoints points are computed, one more is refused. The step is a
// binary fraction so that the point count is exact in floating point, and the
// range stays inside the validity of the BK7 dispersion model (400–2900 nm).
func TestSpectrumPointLimit(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	c := postMap(t, ts.URL, "/api/coating/design", `{"template":"single-ar","lambda_nm":550}`)
	stack, err := json.Marshal(c["stack"])
	if err != nil {
		t.Fatal(err)
	}
	const step = 0.25
	spectrum := func(points int) (int, []byte) {
		to := 400 + float64(points-1)*step
		body := `{"stack":` + string(stack) + `,"from_nm":400,"to_nm":` + strconv.FormatFloat(to, 'f', -1, 64) +
			`,"step_nm":` + strconv.FormatFloat(step, 'f', -1, 64) + `}`
		code, _, out := doHTTP(t, ts.URL, http.MethodPost, "/api/coating/spectrum", body)
		return code, out
	}
	code, out := spectrum(maxSpectrumPoints)
	if code != http.StatusOK {
		t.Fatalf("at the limit: %d %s", code, out)
	}
	var m struct {
		Rows []json.RawMessage `json:"rows"`
	}
	if err := json.Unmarshal(out, &m); err != nil || len(m.Rows) != maxSpectrumPoints {
		t.Fatalf("at the limit: %d rows, want %d (%v)", len(m.Rows), maxSpectrumPoints, err)
	}
	if code, out := spectrum(maxSpectrumPoints + 1); code != http.StatusBadRequest {
		t.Fatalf("one point over the limit: %d %s, want 400", code, out)
	}
}
