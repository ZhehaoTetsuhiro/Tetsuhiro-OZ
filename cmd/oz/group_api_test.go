package main

// HTTP tests of the lens-group endpoints (POST /api/lens/group/*) and of
// POST /api/lens/parax (task wiring2). The fixtures are the ones of
// group_cli_test.go, so the CLI and the service are checked against the same
// closed forms; every endpoint is called twice and must answer byte-identically
// (the determinism clause of docs/API.md §3).

import (
	"encoding/json"
	"math"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// groupPost posts a body and returns the raw answer (any status).
func groupPost(t *testing.T, base, path, body string) (int, []byte) {
	t.Helper()
	code, _, buf := doHTTP(t, base, http.MethodPost, path, body)
	return code, buf
}

// groupCall posts a body, requires 200 and decodes the answer into v.
func groupCall(t *testing.T, base, path, body string, v any) []byte {
	t.Helper()
	code, buf := groupPost(t, base, path, body)
	if code != http.StatusOK {
		t.Fatalf("POST %s: %d %s", path, code, buf)
	}
	if v != nil {
		if err := json.Unmarshal(buf, v); err != nil {
			t.Fatalf("POST %s: %v (%s)", path, err, buf)
		}
	}
	return buf
}

// groupJSON marshals a group as the wire form used at the top level of the
// requests (the "flat" form) or nested under "group".
func groupJSON(t *testing.T, g *lens.Group) string {
	t.Helper()
	buf, err := json.Marshal(g)
	if err != nil {
		t.Fatal(err)
	}
	return string(buf)
}

// TestGroupEndpointDeterminism posts the same body twice to every new endpoint
// and requires byte-identical answers.
func TestGroupEndpointDeterminism(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	cass := groupJSON(t, gCassegrainFixture())
	dbl := groupJSON(t, gDoubletFixture(t, 100, 12, 6, 3))
	kep := groupJSON(t, gKeplerFixture(t, false))
	cases := []struct{ path, body string }{
		{"/api/lens/group/layout", cass},
		{"/api/lens/group/image", cass},
		{"/api/lens/group/prescription", cass},
		{"/api/lens/group/afocal", kep},
		{"/api/lens/group/dispersion", dbl},
		{"/api/lens/group/set", `{"group":` + cass + `,"sets":[{"param":"e1.front.R","value":-1800}]}`},
	}
	for _, c := range cases {
		code1, b1 := groupPost(t, ts.URL, c.path, c.body)
		code2, b2 := groupPost(t, ts.URL, c.path, c.body)
		if code1 != http.StatusOK || code2 != http.StatusOK {
			t.Errorf("%s answered %d/%d: %s", c.path, code1, code2, b1)
			continue
		}
		if string(b1) != string(b2) {
			t.Errorf("%s is not deterministic:\n%s\n%s", c.path, b1, b2)
		}
	}
	// parax: the library Cassegrain prescription.
	rx, err := os.ReadFile(filepath.Join("..", "..", "library", "geo", "mirror-cassegrain-f1000-d300-b300.json"))
	if err != nil {
		t.Skipf("library file not readable: %v", err)
	}
	body := `{"prescription":` + string(rx) + `}`
	_, b1 := groupPost(t, ts.URL, "/api/lens/parax", body)
	_, b2 := groupPost(t, ts.URL, "/api/lens/parax", body)
	if string(b1) != string(b2) {
		t.Errorf("parax is not deterministic:\n%s\n%s", b1, b2)
	}
}

// TestGroupLayoutEndpoint checks the surface table of the Cassegrain over HTTP
// (vertices 0 / 60 / −240, one fold in the folded Kepler) and that the flat and
// the nested request forms are the same request.
func TestGroupLayoutEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	cass := groupJSON(t, gCassegrainFixture())
	var res groupLayoutResponse
	flat := groupCall(t, ts.URL, "/api/lens/group/layout", cass, &res)
	if res.NumSurfaces != 3 || res.Folds != 0 {
		t.Fatalf("%d surfaces, %d folds, want 3 and 0", res.NumSurfaces, res.Folds)
	}
	for i, wantZ := range []float64{0, 60, -240} {
		if d := math.Abs(res.Surfaces[i].Vertex.Z - wantZ); d > 1e-12 {
			t.Errorf("surface %d vertex z = %g, want %g", i, res.Surfaces[i].Vertex.Z, wantZ)
		}
	}
	if d := math.Abs(res.ExitVertex.Z + 240); d > 1e-12 {
		t.Errorf("exit vertex z = %g, want −240", res.ExitVertex.Z)
	}
	if res.ExitAxis.Z != 1 || res.Surfaces[0].Face != "front" || res.Surfaces[0].Axis.Z != 1 {
		t.Errorf("exit axis/face/axis: %+v", res)
	}
	// Nested form = flat form.
	nested := groupCall(t, ts.URL, "/api/lens/group/layout", `{"group":`+cass+`}`, nil)
	if string(nested) != string(flat) {
		t.Errorf("the flat and the nested request forms differ:\n%s\n%s", flat, nested)
	}
	// Folded Kepler: one fold, the chain runs along +y after the mirror.
	var fk groupLayoutResponse
	groupCall(t, ts.URL, "/api/lens/group/layout", groupJSON(t, gKeplerFixture(t, true)), &fk)
	if fk.NumSurfaces != 6 || fk.Folds != 1 {
		t.Fatalf("folded: %d surfaces, %d folds, want 6 and 1", fk.NumSurfaces, fk.Folds)
	}
	if d := math.Abs(fk.ExitAxis.Y - 1); d > 1e-12 {
		t.Errorf("folded exit axis y = %g, want 1", fk.ExitAxis.Y)
	}
	if fk.Surfaces[2].Folds != 1 || fk.Surfaces[2].Vertex.Z != 15 {
		t.Errorf("post-fold surface: %+v", fk.Surfaces[2])
	}
}

// TestGroupImageEndpoint checks the image measurement of the parabola (R/2) and
// of the Cassegrain (a virtual focus 601.7 mm past the secondary, matching the
// paraxial −840 mm within the spherical aberration) plus the finite-object path.
func TestGroupImageEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	var res groupImageResponse
	groupCall(t, ts.URL, "/api/lens/group/image", groupJSON(t, gParabolaFixture()), &res)
	if d := math.Abs(res.Point.Z + 200); d > 1e-9 {
		t.Errorf("image z = %g, want −200 (R/2)", res.Point.Z)
	}
	if d := math.Abs(res.FocusMM - 200); d > 1e-9 {
		t.Errorf("focus distance = %g, want 200", res.FocusMM)
	}
	if res.RMSMM > 1e-10 || res.Rays != 37 || res.Failed != 0 || res.Samples != 3 {
		t.Errorf("parabola image: %+v", res)
	}
	if res.ObjectZMM != nil || res.LambdaNm != geo.Wavelengthd {
		t.Errorf("parabola bundle: object %v, lambda %g", res.ObjectZMM, res.LambdaNm)
	}

	var cass groupImageResponse
	groupCall(t, ts.URL, "/api/lens/group/image", groupJSON(t, gCassegrainFixture()), &cass)
	t.Logf("cassegrain image (%g, %g, %g), focus %g mm, RMS %g mm, rays %d",
		cass.Point.X, cass.Point.Y, cass.Point.Z, cass.FocusMM, cass.RMSMM, cass.Rays)
	// The real best focus sits 1.7 mm behind the paraxial image (the two
	// spherical mirrors); the threshold is the aberration budget of that pair.
	if d := math.Abs(cass.Point.Z + 841.709370442); d > 0.05 {
		t.Errorf("cassegrain image z = %g, want −841.709 (Δ %g)", cass.Point.Z, d)
	}
	if d := math.Abs(cass.Point.Z + 840); d > 2 {
		t.Errorf("the real image is %g mm away from the paraxial −840 mm", d)
	}
	if cass.RMSMM > 0.2 || cass.Rays != 37 || cass.Failed != 0 {
		t.Errorf("cassegrain image: %+v", cass)
	}

	// Finite object: the doublet at 1000 mm still focuses (and the answer is
	// finite and deterministic).
	var fin groupImageResponse
	groupCall(t, ts.URL, "/api/lens/group/image", `{"group":`+groupJSON(t, gDoubletFixture(t, 100, 12, 6, 3))+
		`,"object_z_mm":-1000,"field_tan":0.5}`, &fin)
	if fin.ObjectZMM == nil || *fin.ObjectZMM != -1000 {
		t.Fatalf("finite object not echoed: %+v", fin.ObjectZMM)
	}
	if fin.FieldTanY != 0.5 {
		t.Errorf("field echo = %g, want the object height 0.5 mm", fin.FieldTanY)
	}
	for _, v := range []float64{fin.Point.X, fin.Point.Y, fin.Point.Z, fin.RMSMM, fin.FocusMM} {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			t.Fatalf("finite-object image is not finite: %+v", fin)
		}
	}
}

// TestGroupAfocalEndpoint checks Γ of the Kepler telescope over HTTP and that a
// focal system is refused with the explanation of the lens package.
func TestGroupAfocalEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	var res groupAfocalResponse
	groupCall(t, ts.URL, "/api/lens/group/afocal", groupJSON(t, gKeplerFixture(t, false)), &res)
	if d := math.Abs(res.Gamma + 1); d > 1e-9 {
		t.Errorf("Γ = %g, want −1", res.Gamma)
	}
	if d := math.Abs(res.BeamRatio - 1); d > 0.02 {
		t.Errorf("D_in/D_out = %g, want ≈ 1", res.BeamRatio)
	}
	if res.AxisOut.Z != 1 || res.Rays == 0 {
		t.Errorf("kepler afocal: %+v", res)
	}
	code, buf := groupPost(t, ts.URL, "/api/lens/group/afocal", groupJSON(t, gCassegrainFixture()))
	if code != http.StatusBadRequest || !strings.Contains(string(buf), "not afocal") {
		t.Errorf("a focal system must be refused: %d %s", code, buf)
	}
}

// TestGroupDispersionEndpoint checks the zero dispersion of a mirror group and
// the residual F–C shift of the cemented achromat.
func TestGroupDispersionEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	var mirror groupDispersionResponse
	groupCall(t, ts.URL, "/api/lens/group/dispersion", groupJSON(t, gParabolaFixture()), &mirror)
	if mirror.ShiftMM != 0 || mirror.A.FocusMM != mirror.B.FocusMM {
		t.Errorf("a pure mirror system must not disperse: %+v", mirror)
	}
	if mirror.LambdaANm != geo.WavelengthF || mirror.LambdaBNm != geo.WavelengthC {
		t.Errorf("default wavelengths %g/%g, want F/C", mirror.LambdaANm, mirror.LambdaBNm)
	}
	var dbl groupDispersionResponse
	groupCall(t, ts.URL, "/api/lens/group/dispersion", groupJSON(t, gDoubletFixture(t, 100, 12, 6, 3)), &dbl)
	if !(dbl.ShiftMM < 0) {
		t.Errorf("normal dispersion: shift %g, want < 0", dbl.ShiftMM)
	}
	if d := math.Abs(dbl.ShiftMM + 0.4065509316); d > 0.01 {
		t.Errorf("doublet F–C shift = %g (Δ %g)", dbl.ShiftMM, d)
	}
	if d := math.Abs(dbl.A.FocusMM - dbl.B.FocusMM); d < 1e-3 {
		t.Errorf("the achromat shows no measurable colour: %g mm", d)
	}
}

// TestGroupSetEndpoint edits a group over HTTP: the updated group is returned,
// a bad parameter and a rejected edit are 400 with the error of the API, and
// the work limit on the number of edits is enforced.
func TestGroupSetEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	cass := groupJSON(t, gCassegrainFixture())
	var res groupSetResponse
	groupCall(t, ts.URL, "/api/lens/group/set",
		`{"group":`+cass+`,"sets":[{"param":"e1.front.R","value":-1600},{"param":"e1.tilt_x_deg","value":0}]}`, &res)
	if res.Applied != 2 || res.Group == nil || len(res.Group.Elements) != 3 {
		t.Fatalf("set answer: %+v", res)
	}
	if got := res.Group.Elements[1].Front.RadiusMM; got != -1600 {
		t.Errorf("e1.front.R = %g, want −1600", got)
	}
	// The edit is visible in the next measurement: R = −1600 moves the focus.
	var img groupImageResponse
	groupCall(t, ts.URL, "/api/lens/group/image", `{"group":`+groupJSON(t, res.Group)+`}`, &img)
	if math.Abs(img.FocusMM) < 1 {
		t.Errorf("the edited group still focuses at %g mm", img.FocusMM)
	}
	updated, err := json.Marshal(res.Group)
	if err != nil {
		t.Fatal(err)
	}
	var img2 groupImageResponse
	groupCall(t, ts.URL, "/api/lens/group/image", `{"group":`+string(updated)+`}`, &img2)
	if img2.FocusMM != img.FocusMM {
		t.Errorf("the returned group is not the group that was measured: %g vs %g", img2.FocusMM, img.FocusMM)
	}

	for _, c := range []struct{ name, body, want string }{
		{"missing group", `{"sets":[{"param":"e1.front.R","value":1}]}`, "缺少 group"},
		{"missing sets", `{"group":` + cass + `}`, "缺少 sets"},
		{"unknown parameter", `{"group":` + cass + `,"sets":[{"param":"e1.front.bogus","value":1}]}`, "无法识别"},
		{"rejected edit", `{"group":` + cass + `,"sets":[{"param":"e0.rear.R","value":10}]}`, "has no rear surface"},
		{"out of range", `{"group":` + cass + `,"sets":[{"param":"e9.front.R","value":10}]}`, "out of range"},
	} {
		code, buf := groupPost(t, ts.URL, "/api/lens/group/set", c.body)
		if code != http.StatusBadRequest {
			t.Errorf("%s: code %d, want 400 (%s)", c.name, code, buf)
			continue
		}
		if !strings.Contains(string(buf), c.want) {
			t.Errorf("%s: error %s lacks %q", c.name, buf, c.want)
		}
	}
	// The edit limit is refused before any work happens.
	sets := make([]string, maxGroupSets+1)
	for i := range sets {
		sets[i] = `{"param":"e1.front.R","value":-2000}`
	}
	code, buf := groupPost(t, ts.URL, "/api/lens/group/set",
		`{"group":`+cass+`,"sets":[`+strings.Join(sets, ",")+`]}`)
	if code != http.StatusBadRequest || !strings.Contains(string(buf), "超过上限") {
		t.Errorf("the edit limit is not enforced: %d %s", code, buf)
	}
}

// TestGroupPrescriptionEndpoint converts the Cassegrain group over HTTP and
// compares the cards with the library file.
func TestGroupPrescriptionEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	var out struct {
		Prescription *lens.Prescription `json:"prescription"`
	}
	groupCall(t, ts.URL, "/api/lens/group/prescription", groupJSON(t, gCassegrainFixture()), &out)
	if out.Prescription == nil || len(out.Prescription.Cards) != 3 {
		t.Fatalf("prescription answer: %+v", out)
	}
	want, err := lens.LoadJSON(filepath.Join("..", "..", "library", "geo", "mirror-cassegrain-f1000-d300-b300.json"))
	if err != nil {
		t.Skipf("library file not readable: %v", err)
	}
	for i := range want.Cards {
		a, b := out.Prescription.Cards[i], want.Cards[i]
		if a.RadiusMM != b.RadiusMM || a.ThicknessMM != b.ThicknessMM || a.Mirror != b.Mirror {
			t.Errorf("card %d: %+v vs library %+v", i, a, b)
		}
	}
	// A 3D fold in the middle of the chain has no card form: 400 with the
	// explanation of the lens package.
	code, buf := groupPost(t, ts.URL, "/api/lens/group/prescription", groupJSON(t, gKeplerFixture(t, true)))
	if code != http.StatusBadRequest || !strings.Contains(string(buf), "out-of-plane 3D fold") {
		t.Errorf("a mid-chain fold must be refused: %d %s", code, buf)
	}
}

// TestParaxEndpoint checks POST /api/lens/parax: the library Cassegrain (two
// mirrors, EFL = −857.142857, image −840, parity +1, F–C exactly zero), the
// finite-object override and the afocal Kepler with its field tangent.
func TestParaxEndpoint(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	rx, err := os.ReadFile(filepath.Join("..", "..", "library", "geo", "mirror-cassegrain-f1000-d300-b300.json"))
	if err != nil {
		t.Skipf("library file not readable: %v", err)
	}
	var res paraxResponse
	groupCall(t, ts.URL, "/api/lens/parax", `{"prescription":`+string(rx)+`}`, &res)
	if res.Mirrors != 2 || res.Parity != 1 || res.Afocal {
		t.Fatalf("mirrors %d, parity %g, afocal %v", res.Mirrors, res.Parity, res.Afocal)
	}
	if res.EFLMM == nil || math.Abs(*res.EFLMM+857.1428571428571) > 1e-9 {
		t.Errorf("EFL = %v", res.EFLMM)
	}
	if res.ImageZMM == nil || math.Abs(*res.ImageZMM+840) > 1e-9 {
		t.Errorf("ImageZ = %v", res.ImageZMM)
	}
	if res.FCShiftMM == nil || *res.FCShiftMM != 0 {
		t.Errorf("a pure mirror system has no F–C shift: %v", res.FCShiftMM)
	}
	if res.Gamma != nil || res.GammaNote == "" {
		t.Errorf("Γ of an infinite object: %v (%q)", res.Gamma, res.GammaNote)
	}

	// A single mirror with a finite object: ImageZ = s·R/(R+2s) = −1500.
	one := &lens.Prescription{
		Name: "one-mirror", ObjectZMM: math.Inf(-1),
		Cards: []lens.Card{{Name: "mirror", RadiusMM: -2000, Mirror: true, SemiAperMM: 80}},
	}
	buf, err := json.Marshal(one)
	if err != nil {
		t.Fatal(err)
	}
	var r2 paraxResponse
	groupCall(t, ts.URL, "/api/lens/parax", `{"prescription":`+string(buf)+`,"object_z_mm":-3000}`, &r2)
	if r2.ImageZMM == nil || math.Abs(*r2.ImageZMM+1500) > 1e-9 {
		t.Errorf("finite-object image z = %v, want −1500", r2.ImageZMM)
	}
	if r2.Gamma == nil || math.Abs(*r2.Gamma+2) > 1e-9 {
		t.Errorf("finite-object γ = %v, want −2", r2.Gamma)
	}

	// The afocal Kepler prescription: Γ needs the field tangent.
	keplerPath := filepath.Join(t.TempDir(), "kepler.json")
	if err := cmdLens([]string{"group", "to-prescription", "-o", keplerPath,
		writeGroup(t, "kepler", gKeplerFixture(t, false))}); err != nil {
		t.Fatal(err)
	}
	krx, err := os.ReadFile(keplerPath)
	if err != nil {
		t.Fatal(err)
	}
	var r3 paraxResponse
	groupCall(t, ts.URL, "/api/lens/parax", `{"prescription":`+string(krx)+`,"field_tan":0.001}`, &r3)
	if !r3.Afocal || r3.EFLMM != nil {
		t.Errorf("the Kepler prescription: afocal %v, EFL %v", r3.Afocal, r3.EFLMM)
	}
	if r3.Gamma == nil || math.Abs(*r3.Gamma+1) > 1e-9 {
		t.Errorf("kepler Γ = %v, want −1", r3.Gamma)
	}
	if r3.FCNote == "" {
		t.Error("the F–C shift of an afocal system is not explained")
	}
}

// TestGroupEndpointErrors checks the error contract: every new endpoint is
// POST-only (405 otherwise, over the whole mux), an invalid request is a 400
// with a JSON error, and the work limits are refused before any computation.
func TestGroupEndpointErrors(t *testing.T) {
	ts := httptest.NewServer(ozMux())
	defer ts.Close()
	cass := groupJSON(t, gCassegrainFixture())
	paths := []string{
		"/api/lens/group/layout", "/api/lens/group/image", "/api/lens/group/afocal",
		"/api/lens/group/dispersion", "/api/lens/group/set", "/api/lens/group/prescription",
		"/api/lens/parax",
	}
	for _, p := range paths {
		code, _, buf := doHTTP(t, ts.URL, http.MethodGet, p, "")
		if code != http.StatusMethodNotAllowed {
			t.Errorf("GET %s = %d, want 405 (%s)", p, code, buf)
		}
	}
	cases := []struct {
		name, path, body string
	}{
		{"empty body", "/api/lens/group/layout", `{}`},
		{"elements null", "/api/lens/group/layout", `{"elements":[]}`},
		{"unknown kind", "/api/lens/group/layout", `{"elements":[{"kind":"banana","front":{}}]}`},
		{"stop with glass", "/api/lens/group/layout",
			`{"elements":[{"kind":"stop","front":{},"glass":"N-BK7","thickness_mm":0}]}`},
		{"bad object", "/api/lens/group/layout", `{"elements":[{"kind":"stop","front":{"semi_aperture_mm":5}}],"object_z_mm":"near"}`},
		{"malformed body", "/api/lens/group/image", `{"elements":`},
		{"rings over the limit", "/api/lens/group/image", `{"group":` + cass + `,"samples":99}`},
		{"scan over the limit", "/api/lens/group/image", `{"group":` + cass + `,"scan":99999}`},
		{"negative range", "/api/lens/group/image", `{"group":` + cass + `,"range_mm":-1}`},
		{"negative lambda", "/api/lens/group/image", `{"group":` + cass + `,"lambda_nm":-5}`},
		{"finite object with field_deg", "/api/lens/group/image", `{"group":` + cass + `,"object_z_mm":-5000,"field_deg":0.5}`},
		{"afocal of a mirror system", "/api/lens/group/afocal", `{"group":` + cass + `}`},
		{"set without sets", "/api/lens/group/set", `{"group":` + cass + `,"sets":[]}`},
		{"prescription of a 3D fold", "/api/lens/group/prescription", groupJSON(t, gKeplerFixture(t, true))},
		{"parax without a prescription", "/api/lens/parax", `{}`},
		{"parax bad object", "/api/lens/parax", `{"prescription":{"cards":[{"radius_mm":-100}]},"object_z_mm":"x"}`},
	}
	for _, c := range cases {
		code, buf := groupPost(t, ts.URL, c.path, c.body)
		if code != http.StatusBadRequest {
			t.Errorf("%s: %d, want 400 (%s)", c.name, code, buf)
			continue
		}
		var m map[string]any
		if err := json.Unmarshal(buf, &m); err != nil || m["error"] == nil {
			t.Errorf("%s: the body is not a JSON error: %s", c.name, buf)
		}
	}
}

// TestGroupExamplesInGUI decodes the example groups embedded in
// server/web/oz.js (marked with /*oz:example*/ … /*:oz:example*/): the GUI
// ships real, valid lens groups, and the numbers of the Cassegrain example are
// the library ones.
func TestGroupExamplesInGUI(t *testing.T) {
	src, err := os.ReadFile(filepath.Join("..", "..", "server", "web", "oz.js"))
	if err != nil {
		t.Fatalf("oz.js: %v", err)
	}
	text := string(src)
	const open, close = "/*oz:example*/", "/*:oz:example*/"
	var groups []*lens.Group
	rest := text
	for {
		i := strings.Index(rest, open)
		if i < 0 {
			break
		}
		rest = rest[i+len(open):]
		j := strings.Index(rest, close)
		if j < 0 {
			t.Fatalf("unterminated example block in oz.js")
		}
		var g lens.Group
		if err := json.Unmarshal([]byte(rest[:j]), &g); err != nil {
			t.Fatalf("example %d is not a valid group JSON: %v", len(groups), err)
		}
		if err := g.Validate(); err != nil {
			t.Fatalf("example %d does not validate: %v", len(groups), err)
		}
		if _, err := g.System(); err != nil {
			t.Fatalf("example %d has no system: %v", len(groups), err)
		}
		groups = append(groups, &g)
		rest = rest[j:]
	}
	if len(groups) < 2 {
		t.Fatalf("%d example groups in oz.js, want at least 2", len(groups))
	}
	// The first example is the library Cassegrain: three elements, the two
	// mirror radii of library/geo/mirror-cassegrain-f1000-d300-b300.json.
	cass := groups[0]
	if len(cass.Elements) != 3 || cass.Elements[1].Kind != lens.ElementMirror {
		t.Fatalf("example 0: %+v", cass.Elements)
	}
	if cass.Elements[1].Front.RadiusMM != -2000 || cass.Elements[2].Front.RadiusMM != -646.1538461538461 {
		t.Errorf("example 0 radii: %g, %g", cass.Elements[1].Front.RadiusMM, cass.Elements[2].Front.RadiusMM)
	}
	// A lens example (doublet) and a folded example must be present.
	var lenses, folded int
	for _, g := range groups {
		for _, e := range g.Elements {
			if e.Kind == lens.ElementLens && e.Rear != nil {
				lenses++
			}
			if e.Kind == lens.ElementMirror && (e.TiltXDeg != 0 || e.TiltYDeg != 0) {
				folded++
			}
		}
	}
	if lenses == 0 || folded == 0 {
		t.Errorf("the examples lack a cemented lens (%d) or a folding mirror (%d)", lenses, folded)
	}
}
