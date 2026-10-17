package server

import (
	"encoding/json"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
	"time"
)

// scanConfig is a small scene with a rotatable slit and an extended source,
// matching the GUI's scan panel.
const scanConfig = `{"grid":{"size":64,"width":0.01},"wavelength":6.328e-7,"polarized":false,"method":"asm",
 "sources":[{"id":"src","type":"extended","pos":{"x":0,"y":0,"z":-0.1},
   "params":{"power":0.001,"samples_x":2,"samples_y":1,"standoff":0.04}}],
 "scene":{"components":[
   {"id":"slit","type":"aperture","label":"单缝","pos":{"x":0,"y":0,"z":-0.05},
    "shape":{"kind":"rectangle","params":{"width":0.0004,"height":0.01}}},
   {"id":"det","type":"sensor","label":"屏","pos":{"x":0,"y":0,"z":0.05},"yaw":3.141592653589793,
    "shape":{"kind":"circle","params":{"radius":0.007}}}]}}`

// scanRun posts a sweep and polls until the run is finished.
func scanRun(t *testing.T, base, body string) map[string]any {
	t.Helper()
	sub := doReq(t, http.MethodPost, base+"/api/scan", body)
	if sub.code != http.StatusAccepted {
		t.Fatalf("POST /api/scan: %d %s", sub.code, sub.body)
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
			t.Fatal("scan did not finish in time")
		}
		time.Sleep(20 * time.Millisecond)
	}
}

// scanPoints extracts the points of a finished scan run. wantParam, when not
// empty, is asserted to be the swept parameter.
func scanPoints(t *testing.T, meta map[string]any, wantParam string) []struct {
	Value  float64 `json:"value"`
	Metric float64 `json:"metric"`
	Error  string  `json:"error"`
} {
	t.Helper()
	raw, err := json.Marshal(meta["scan"])
	if err != nil {
		t.Fatal(err)
	}
	var rep struct {
		Param  string `json:"param"`
		Metric string `json:"metric"`
		Points []struct {
			Value  float64 `json:"value"`
			Metric float64 `json:"metric"`
			Error  string  `json:"error"`
		} `json:"points"`
	}
	if err := json.Unmarshal(raw, &rep); err != nil {
		t.Fatal(err)
	}
	if wantParam != "" && rep.Param != wantParam {
		t.Errorf("param = %q, want %q", rep.Param, wantParam)
	}
	return rep.Points
}

// TestScanEndpointRunsASweep drives POST /api/scan end to end: the sweep is
// accepted, runs, and every point reports the swept value with a metric.
func TestScanEndpointRunsASweep(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	meta := scanRun(t, ts.URL, `{"config":`+scanConfig+`,"scan":{"param":"scene.components[0].roll","values":[0,0.5,1.0],"metric":"power"}}`)
	if meta["kind"] != "scan" {
		t.Errorf("run kind = %v, want scan", meta["kind"])
	}
	pts := scanPoints(t, meta, "scene.components[0].roll")
	if len(pts) != 3 {
		t.Fatalf("got %d points, want 3", len(pts))
	}
	for i, want := range []float64{0, 0.5, 1.0} {
		if pts[i].Value != want {
			t.Errorf("point %d: value %g, want %g", i, pts[i].Value, want)
		}
		if pts[i].Error != "" {
			t.Errorf("point %d: error %q", i, pts[i].Error)
		}
		if pts[i].Metric <= 0 {
			t.Errorf("point %d: metric %g, want > 0", i, pts[i].Metric)
		}
	}
}

// TestScanEndpointRejectsBadRequest covers the request-level validation.
func TestScanEndpointRejectsBadRequest(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	cases := []struct {
		name, body string
		want       int
	}{
		{"no param", `{"config":` + scanConfig + `,"scan":{"values":[0],"metric":"power"}}`, http.StatusBadRequest},
		{"no values", `{"config":` + scanConfig + `,"scan":{"param":"scene.components[0].roll","metric":"power"}}`, http.StatusBadRequest},
		{"bad json", `{`, http.StatusBadRequest},
		{"invalid config", `{"config":{"grid":{"size":0,"width":0},"sources":[]},"scan":{"param":"wavelength","values":[6e-7]}}`, http.StatusBadRequest},
	}
	for _, tc := range cases {
		if r := doReq(t, http.MethodPost, ts.URL+"/api/scan", tc.body); r.code != tc.want {
			t.Errorf("%s: got %d, want %d (%s)", tc.name, r.code, tc.want, r.body)
		}
	}
	if r := doReq(t, http.MethodGet, ts.URL+"/api/scan", ""); r.code != http.StatusMethodNotAllowed {
		t.Errorf("GET /api/scan = %d, want 405", r.code)
	}
}

// TestScanEndpointSweepsWavelength is the end-to-end physics check over HTTP:
// sweeping the wavelength must report every requested value, and the metric
// must change between them.
func TestScanEndpointSweepsWavelength(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	meta := scanRun(t, ts.URL, `{"config":`+scanConfig+`,"scan":{"param":"wavelength","values":[5e-7,8e-7],"metric":"power"}}`)
	pts := scanPoints(t, meta, "wavelength")
	if len(pts) != 2 {
		t.Fatalf("got %d points, want 2", len(pts))
	}
	if pts[0].Value != 5e-7 || pts[1].Value != 8e-7 {
		t.Errorf("values = %g, %g; want 5e-07, 8e-07", pts[0].Value, pts[1].Value)
	}
	if pts[0].Error != "" || pts[1].Error != "" {
		t.Errorf("errors: %q / %q", pts[0].Error, pts[1].Error)
	}
	if pts[0].Metric == pts[1].Metric {
		t.Errorf("the metric did not change with wavelength (%g == %g)", pts[0].Metric, pts[1].Metric)
	}
}

// TestScanEndpointReportsBadValue: a value that makes the configuration
// invalid is reported per point, and the sweep still finishes.
func TestScanEndpointReportsBadValue(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	meta := scanRun(t, ts.URL, `{"config":`+scanConfig+`,"scan":{"param":"scene.components[0].shape.params.width","values":[0.0004,-1,0.0008],"metric":"power"}}`)
	pts := scanPoints(t, meta, "scene.components[0].shape.params.width")
	if len(pts) != 3 {
		t.Fatalf("got %d points, want 3", len(pts))
	}
	if pts[0].Error != "" || pts[2].Error != "" {
		t.Errorf("valid values must succeed: %+v", pts)
	}
	if pts[1].Error == "" {
		t.Error("a negative slit width must be reported as an error for that point")
	}
}

// TestScanEndpointUnknownRun: a scan run's id is a normal run id, so an
// unknown id still 404s and a scan run has no planes to serve.
func TestScanEndpointUnknownRun(t *testing.T) {
	ts := httptest.NewServer(New(64 << 20).Handler())
	defer ts.Close()
	meta := scanRun(t, ts.URL, `{"config":`+scanConfig+`,"scan":{"param":"wavelength","values":[6e-7],"metric":"power"}}`)
	rid, _ := meta["run_id"].(string)
	if r := doReq(t, http.MethodGet, ts.URL+"/api/runs/"+rid+"/planes/det", ""); r.code != http.StatusNotFound {
		t.Errorf("planes of a scan run = %d, want 404 (a sweep has no planes)", r.code)
	}
	if r := doReq(t, http.MethodGet, ts.URL+"/api/runs/does-not-exist", ""); r.code != http.StatusNotFound {
		t.Errorf("unknown run = %d, want 404", r.code)
	}
	// The report must name the metric it reports.
	raw, _ := json.Marshal(meta["scan"])
	if !strings.Contains(string(raw), `"metric":"power"`) {
		t.Errorf("report lacks the metric name: %s", raw)
	}
}
