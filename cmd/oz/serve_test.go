package main

import (
	"bytes"
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
)

// TestServeEndpoints drives the HTTP API end to end: design a lens, design a
// coating, scan its spectrum — all over real HTTP against ozMux().
func TestServeEndpoints(t *testing.T) {
	srv := httptest.NewServer(ozMux())
	defer srv.Close()

	get := func(path string) string {
		t.Helper()
		resp, err := http.Get(srv.URL + path)
		if err != nil {
			t.Fatal(err)
		}
		defer resp.Body.Close()
		if resp.StatusCode != 200 {
			t.Fatalf("GET %s: %d", path, resp.StatusCode)
		}
		var buf bytes.Buffer
		if _, err := buf.ReadFrom(resp.Body); err != nil {
			t.Fatal(err)
		}
		return buf.String()
	}
	post := func(path, body string) map[string]any {
		t.Helper()
		resp, err := http.Post(srv.URL+path, "application/json", strings.NewReader(body))
		if err != nil {
			t.Fatal(err)
		}
		defer resp.Body.Close()
		if resp.StatusCode != 200 {
			var buf bytes.Buffer
			buf.ReadFrom(resp.Body)
			t.Fatalf("POST %s: %d %s", path, resp.StatusCode, buf.String())
		}
		var out map[string]any
		if err := json.NewDecoder(resp.Body).Decode(&out); err != nil {
			t.Fatal(err)
		}
		return out
	}

	// Version, then the GUI. The geometry/coating tab is part of the one page
	// at /; the page, its stylesheet and its script must stay mobile-ready and
	// must reach every OZ endpoint.
	if !strings.Contains(get("/api/version"), "0.1.0-Preview") {
		t.Fatal("version endpoint missing version")
	}
	page := get("/")
	if !strings.Contains(page, "Tetsuhiro OZ") {
		t.Fatal("GUI page missing")
	}
	for _, want := range []string{
		`<meta name="viewport" content="width=device-width, initial-scale=1">`,
		`id="ozModeBtn"`, `id="ozEditor"`, `id="ozOutput"`, `id="ozPlot"`,
		`<script src="oz.js"></script>`,
	} {
		if !strings.Contains(page, want) {
			t.Errorf("GUI page missing %q", want)
		}
	}
	css := get("/style.css")
	for _, want := range []string{
		"height: clamp(",            // plot height follows the viewport width, capped
		"@media (max-width: 700px)", // single-column layout on phones
		"@media (pointer: coarse)",  // touch targets on touch screens
		"#ozGroupPlot",              // the 3D group canvas keeps its own responsive height
		"touch-action: none",        // ... and receives the drag gestures on touch screens
	} {
		if !strings.Contains(css, want) {
			t.Errorf("style.css missing %q", want)
		}
	}
	js := get("/oz.js")
	for _, want := range []string{
		"devicePixelRatio", `addEventListener("resize"`, "orientationchange",
		"/api/lens/design", "/api/lens/eval", "/api/lens/tol",
		"/api/lens/thermal", "/api/lens/sens",
		"/api/lens/group/layout", "/api/lens/group/image", "/api/lens/group/afocal",
		"/api/lens/group/dispersion", "/api/lens/group/set", "/api/lens/group/prescription",
		"ozGroup", "ozGroupPlot", "pointerdown",
		"/api/coating/design", "/api/coating/spectrum", "/api/version",
	} {
		if !strings.Contains(js, want) {
			t.Errorf("oz.js missing %q", want)
		}
	}

	// Lens design.
	d := post("/api/lens/design", `{"efl_mm":60,"semi_aperture_mm":6,"field_deg":0.5}`)
	if d["prescription"] == nil || d["optimizer"] == nil {
		t.Fatalf("design response incomplete: %v", d)
	}
	presc, _ := json.Marshal(d["prescription"])

	// Lens eval of the designed prescription.
	e := post("/api/lens/eval", string(presc))
	if e["report"] == nil {
		t.Fatalf("eval response incomplete: %v", e)
	}

	// Coating design + spectrum.
	c := post("/api/coating/design", `{"template":"single-ar","lambda_nm":550,"substrate":"BK7"}`)
	stack, ok := c["stack"]
	if !ok {
		t.Fatalf("coating design response incomplete: %v", c)
	}
	specBody, _ := json.Marshal(map[string]any{
		"stack": stack, "from_nm": 500, "to_nm": 600, "step_nm": 10, "pol": "u",
	})
	s := post("/api/coating/spectrum", string(specBody))
	rows, ok := s["rows"].([]any)
	if !ok || len(rows) == 0 {
		t.Fatalf("spectrum response incomplete: %v", s)
	}
	first := rows[0].(map[string]any)
	if _, ok := first["R"]; !ok {
		t.Fatalf("spectrum row without R: %v", first)
	}
}
