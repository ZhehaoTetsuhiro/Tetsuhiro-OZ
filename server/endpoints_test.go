package server

import (
	"bytes"
	"encoding/binary"
	"encoding/json"
	"image/png"
	"io"
	"math"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
	"time"

	"tetsuhiro-oz/optics"
)

// simConfigJSON is the smallest sequential system that exercises a full run:
// a Gaussian beam, a propagation, a lens and a sensor plane (128×128 grid).
const simConfigJSON = `{"grid":{"size":128,"width":0.004},"wavelength":6.328e-7,"method":"asm",
 "source":{"type":"gaussian","params":{"waist":0.0005}},
 "elements":[{"type":"propagate","params":{"distance":0.1}},
             {"type":"lens","params":{"f":0.2}},
             {"type":"sensor","params":{"label":"focus"}}]}`

// httpResult is one completed request.
type httpResult struct {
	code   int
	header http.Header
	body   []byte
}

func doReq(t *testing.T, method, url, body string) httpResult {
	t.Helper()
	var rd io.Reader
	if body != "" {
		rd = strings.NewReader(body)
	}
	req, err := http.NewRequest(method, url, rd)
	if err != nil {
		t.Fatal(err)
	}
	if body != "" {
		req.Header.Set("Content-Type", "application/json")
	}
	res, err := http.DefaultClient.Do(req)
	if err != nil {
		t.Fatal(err)
	}
	defer res.Body.Close()
	buf, err := io.ReadAll(res.Body)
	if err != nil {
		t.Fatal(err)
	}
	return httpResult{code: res.StatusCode, header: res.Header, body: buf}
}

func decodeMap(t *testing.T, b []byte) map[string]any {
	t.Helper()
	var m map[string]any
	if err := json.Unmarshal(b, &m); err != nil {
		t.Fatalf("response is not a JSON object: %v\n%s", err, b)
	}
	return m
}

// submitAndWait posts a config to /api/simulate and polls the run until it is
// done, returning its status document.
func submitAndWait(t *testing.T, base, cfg string) map[string]any {
	t.Helper()
	sub := doReq(t, http.MethodPost, base+"/api/simulate", cfg)
	if sub.code != http.StatusAccepted {
		t.Fatalf("POST /api/simulate: %d %s", sub.code, sub.body)
	}
	id, _ := decodeMap(t, sub.body)["run_id"].(string)
	if id == "" {
		t.Fatalf("no run_id in %s", sub.body)
	}
	deadline := time.Now().Add(120 * time.Second)
	for {
		r := doReq(t, http.MethodGet, base+"/api/runs/"+id, "")
		if r.code != http.StatusOK {
			t.Fatalf("GET run: %d %s", r.code, r.body)
		}
		meta := decodeMap(t, r.body)
		switch meta["status"] {
		case StatusDone:
			return meta
		case StatusError:
			t.Fatalf("run failed: %v", meta["error"])
		}
		if time.Now().After(deadline) {
			t.Fatal("run did not finish in time")
		}
		time.Sleep(20 * time.Millisecond)
	}
}

// TestEndpointCatalog covers GET /api/catalog.
func TestEndpointCatalog(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	r := doReq(t, http.MethodGet, ts.URL+"/api/catalog", "")
	if r.code != http.StatusOK {
		t.Fatalf("status %d", r.code)
	}
	m := decodeMap(t, r.body)
	for _, k := range []string{"elements", "sources", "quantum", "methods"} {
		if _, ok := m[k]; !ok {
			t.Errorf("catalog lacks %q", k)
		}
	}
	if els, _ := m["elements"].([]any); len(els) == 0 {
		t.Error("catalog lists no elements")
	}
	if r := doReq(t, http.MethodPost, ts.URL+"/api/catalog", "{}"); r.code != http.StatusMethodNotAllowed {
		t.Errorf("POST /api/catalog = %d, want 405", r.code)
	}
}

// TestEndpointHealth covers GET /api/health, including the gpu object.
func TestEndpointHealth(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	r := doReq(t, http.MethodGet, ts.URL+"/api/health", "")
	if r.code != http.StatusOK {
		t.Fatalf("status %d", r.code)
	}
	m := decodeMap(t, r.body)
	if m["status"] != "ok" {
		t.Errorf("status = %v", m["status"])
	}
	gpu, ok := m["gpu"].(map[string]any)
	if !ok {
		t.Fatalf("health lacks gpu object: %s", r.body)
	}
	for _, k := range []string{"available", "enabled", "backend", "attempts", "failures", "degraded"} {
		if _, ok := gpu[k]; !ok {
			t.Errorf("gpu object lacks %q", k)
		}
	}
}

// TestEndpointValidate covers POST /api/validate for a valid and an invalid
// configuration and the method and body guards.
func TestEndpointValidate(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	url := ts.URL + "/api/validate"

	ok := decodeMap(t, doReq(t, http.MethodPost, url, simConfigJSON).body)
	if ok["ok"] != true {
		t.Fatalf("valid config reported issues: %v", ok)
	}
	bad := strings.Replace(simConfigJSON, `"size":128`, `"size":3`, 1)
	r := doReq(t, http.MethodPost, url, bad)
	if r.code != http.StatusOK {
		t.Fatalf("invalid config status %d", r.code)
	}
	m := decodeMap(t, r.body)
	if m["ok"] != false {
		t.Fatalf("grid size 3 accepted: %s", r.body)
	}
	if iss, _ := m["issues"].([]any); len(iss) == 0 {
		t.Errorf("no issues listed: %s", r.body)
	}
	if r := doReq(t, http.MethodGet, url, ""); r.code != http.StatusMethodNotAllowed {
		t.Errorf("GET /api/validate = %d, want 405", r.code)
	}
	if r := doReq(t, http.MethodPost, url, "{not json"); r.code != http.StatusBadRequest {
		t.Errorf("malformed body = %d, want 400", r.code)
	}
}

// TestEndpointSimulateAndRuns covers the run lifecycle: POST /api/simulate
// (202), GET /api/runs/{id}, and every sub-resource of a finished run (planes
// as binary, PNG, colour; profiles; inspect; scene) plus the error answers.
func TestEndpointSimulateAndRuns(t *testing.T) {
	ts := httptest.NewServer(New(256 << 20).Handler())
	defer ts.Close()
	base := ts.URL

	// Invalid configs are refused before a run is created.
	if r := doReq(t, http.MethodPost, base+"/api/simulate", strings.Replace(simConfigJSON, `"size":128`, `"size":3`, 1)); r.code != http.StatusBadRequest {
		t.Fatalf("invalid simulate = %d, want 400", r.code)
	}
	if r := doReq(t, http.MethodGet, base+"/api/simulate", ""); r.code != http.StatusMethodNotAllowed {
		t.Fatalf("GET /api/simulate = %d, want 405", r.code)
	}

	meta := submitAndWait(t, base, simConfigJSON)
	if g, _ := meta["grid"].(map[string]any); g == nil || g["size"] != float64(128) {
		t.Fatalf("grid metadata = %v", meta["grid"])
	}
	planes, _ := meta["planes"].([]any)
	if len(planes) == 0 {
		t.Fatal("finished run lists no planes")
	}
	pid, _ := planes[len(planes)-1].(map[string]any)["id"].(string)
	if pid == "" {
		t.Fatalf("plane has no id: %v", planes[len(planes)-1])
	}

	// Binary plane: N×N float32 little-endian, grid size in a header.
	bin := doReq(t, http.MethodGet, base+"/api/runs/"+meta["run_id"].(string)+"/planes/"+pid+"?field=total&fmt=bin", "")
	if bin.code != http.StatusOK {
		t.Fatalf("plane bin: %d %s", bin.code, bin.body)
	}
	if got := bin.header.Get("X-Grid-Size"); got != "128" {
		t.Errorf("X-Grid-Size = %q", got)
	}
	if len(bin.body) != 128*128*4 {
		t.Fatalf("plane bin is %d bytes, want %d", len(bin.body), 128*128*4)
	}
	var sum float64
	for i := 0; i < len(bin.body); i += 4 {
		v := math.Float32frombits(binary.LittleEndian.Uint32(bin.body[i:]))
		if v < 0 || math.IsNaN(float64(v)) {
			t.Fatalf("total intensity %v at %d", v, i/4)
		}
		sum += float64(v)
	}
	if sum <= 0 {
		t.Error("total intensity is zero everywhere")
	}

	// PNG views: intensity (log) and real-wavelength colour.
	for _, q := range []string{"field=total&fmt=png&scale=log&cmap=inferno", "field=color&fmt=png"} {
		r := doReq(t, http.MethodGet, base+"/api/runs/"+meta["run_id"].(string)+"/planes/"+pid+"?"+q, "")
		if r.code != http.StatusOK || r.header.Get("Content-Type") != "image/png" {
			t.Fatalf("png (%s): %d %s", q, r.code, r.header.Get("Content-Type"))
		}
		if _, err := png.Decode(bytes.NewReader(r.body)); err != nil {
			t.Fatalf("png (%s) does not decode: %v", q, err)
		}
	}
	if r := doReq(t, http.MethodGet, base+"/api/runs/"+meta["run_id"].(string)+"/planes/"+pid+"?field=color&fmt=bin", ""); r.code != http.StatusBadRequest {
		t.Errorf("field=color with fmt=bin = %d, want 400", r.code)
	}

	// Profile: 1-D cut, one value per grid column.
	prof := doReq(t, http.MethodGet, base+"/api/runs/"+meta["run_id"].(string)+"/profiles/"+pid+"?axis=x&field=total", "")
	if prof.code != http.StatusOK {
		t.Fatalf("profile: %d %s", prof.code, prof.body)
	}
	pm := decodeMap(t, prof.body)
	if v, _ := pm["v"].([]any); len(v) != 128 {
		t.Errorf("profile has %d samples, want 128", len(v))
	}

	// Inspect: local Jones/Stokes readout at a pixel, also via the alias path.
	for _, u := range []string{
		"/inspect/" + pid + "?x=64&y=64",
		"/planes/" + pid + "/inspect?x=64&y=64",
	} {
		r := doReq(t, http.MethodGet, base+"/api/runs/"+meta["run_id"].(string)+u, "")
		if r.code != http.StatusOK {
			t.Fatalf("inspect %s: %d %s", u, r.code, r.body)
		}
		im := decodeMap(t, r.body)
		if _, ok := im["stokes"]; !ok {
			t.Errorf("inspect lacks stokes: %s", r.body)
		}
	}

	// Scene: present exactly when the run metadata carries one.
	_, hasScene := meta["scene"]
	sc := doReq(t, http.MethodGet, base+"/api/runs/"+meta["run_id"].(string)+"/scene", "")
	if hasScene && sc.code != http.StatusOK {
		t.Errorf("scene metadata present but /scene = %d", sc.code)
	}
	if !hasScene && sc.code != http.StatusNotFound {
		t.Errorf("no scene in metadata but /scene = %d", sc.code)
	}

	// Error answers.
	rid := meta["run_id"].(string)
	for _, c := range []struct {
		method, path string
		want         int
	}{
		{http.MethodGet, "/api/runs/does-not-exist", http.StatusNotFound},
		{http.MethodGet, "/api/runs/" + rid + "/planes/no-such-plane", http.StatusNotFound},
		{http.MethodGet, "/api/runs/" + rid + "/nope", http.StatusNotFound},
		{http.MethodPost, "/api/runs/" + rid, http.StatusMethodNotAllowed},
		{http.MethodGet, "/api/runs/" + rid + "/planes/" + pid + "?fmt=tiff", http.StatusBadRequest},
	} {
		if r := doReq(t, c.method, base+c.path, ""); r.code != c.want {
			t.Errorf("%s %s = %d, want %d", c.method, c.path, r.code, c.want)
		}
	}
}

// TestEndpointConvert covers POST /api/convert: an old element sequence is
// turned into a positioned scene plus its sources.
func TestEndpointConvert(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	r := doReq(t, http.MethodPost, ts.URL+"/api/convert", simConfigJSON)
	if r.code != http.StatusOK {
		t.Fatalf("convert: %d %s", r.code, r.body)
	}
	m := decodeMap(t, r.body)
	if _, ok := m["scene"]; !ok {
		t.Errorf("convert lacks scene: %s", r.body)
	}
	if s, _ := m["sources"].([]any); len(s) == 0 {
		t.Errorf("convert lacks sources: %s", r.body)
	}
	if r := doReq(t, http.MethodGet, ts.URL+"/api/convert", ""); r.code != http.StatusMethodNotAllowed {
		t.Errorf("GET /api/convert = %d, want 405", r.code)
	}
}

// TestEndpointQuantum covers POST /api/quantum in JSON, SVG and PNG form. The
// Hong-Ou-Mandel check is a physics assertion: a 50:50 beam splitter acting on
// |1,1⟩ gives no coincidence, so P(1,1) = 0 and the joint distribution is
// normalised.
func TestEndpointQuantum(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	body := `{"modes":2,"cutoff":4,"state":{"type":"fock","params":{"occupation":[1,1]}},` +
		`"gates":[{"type":"beam_splitter","params":{"mode0":0,"mode1":1,"reflectivity":0.5}}]}`
	url := ts.URL + "/api/quantum"

	r := doReq(t, http.MethodPost, url, body)
	if r.code != http.StatusOK {
		t.Fatalf("quantum json: %d %s", r.code, r.body)
	}
	var res optics.QuantumResult
	if err := json.Unmarshal(r.body, &res); err != nil {
		t.Fatal(err)
	}
	if math.Abs(res.Norm-1) > 1e-9 {
		t.Errorf("norm = %g, want 1", res.Norm)
	}
	base := res.Cutoff + 1
	idx11 := 1 + base*1 // little-endian index of |1,1⟩
	if len(res.JointFull) != base*base {
		t.Fatalf("joint_full has %d entries, want %d", len(res.JointFull), base*base)
	}
	if p := res.JointFull[idx11]; math.Abs(p) > 1e-12 {
		t.Errorf("HOM P(1,1) = %g, want 0", p)
	}
	if len(res.MeanN) != 2 {
		t.Errorf("mean_photons has %d modes, want 2", len(res.MeanN))
	}

	svg := doReq(t, http.MethodPost, url+"?fmt=svg", body)
	if svg.code != http.StatusOK || !strings.HasPrefix(svg.header.Get("Content-Type"), "image/svg+xml") {
		t.Errorf("quantum svg: %d %s", svg.code, svg.header.Get("Content-Type"))
	}
	pngR := doReq(t, http.MethodPost, url+"?fmt=png", body)
	if pngR.code != http.StatusOK {
		t.Fatalf("quantum png: %d", pngR.code)
	}
	if _, err := png.Decode(bytes.NewReader(pngR.body)); err != nil {
		t.Errorf("quantum png does not decode: %v", err)
	}
	if r := doReq(t, http.MethodPost, url, "{"); r.code != http.StatusBadRequest {
		t.Errorf("malformed quantum body = %d, want 400", r.code)
	}
	if r := doReq(t, http.MethodGet, url, ""); r.code != http.StatusMethodNotAllowed {
		t.Errorf("GET /api/quantum = %d, want 405", r.code)
	}
}
