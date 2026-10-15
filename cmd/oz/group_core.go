package main

// Shared core of `oz lens group ...` / `oz lens parax` and of the
// POST /api/lens/group/* and POST /api/lens/parax endpoints (the handlers are
// in group_http.go).
//
// The command line and the service call the same run* functions below, so the
// two cannot drift apart: `oz lens group image --json out.json` writes exactly
// the response body of POST /api/lens/group/image for the same group and
// bundle. Every function is a pure function of its request (no clock, no map
// iteration, no randomness), so equal requests give byte-identical answers.

import (
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"strconv"
	"strings"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// Work limits of the group endpoints. A request above a limit is refused with
// 400 before any computation starts (the policy of maxSpectrumPoints and
// friends in serve.go).
const (
	maxGroupRings    = 8    // /api/lens/group/*: Hexapolar sampling rings
	maxGroupScan     = 4096 // /api/lens/group/image: focus-search samples
	maxGroupSets     = 512  // /api/lens/group/set: parameter edits
	maxGroupElements = 512  // group size
)

// ---------------------------------------------------------------------------
// The numeric wire form
// ---------------------------------------------------------------------------

// vec3Out is the numeric JSON form of a geo.Vec3 (mm, global frame).
type vec3Out struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
	Z float64 `json:"z"`
}

// vec3OutOf converts a geo.Vec3.
func vec3OutOf(v geo.Vec3) vec3Out { return vec3Out{X: v.X, Y: v.Y, Z: v.Z} }

// finPtr returns a pointer to v, or nil when v is not finite (NaN or ±Inf), so
// that a quantity which does not exist (an afocal EFL, the magnification of an
// infinite object) is serialised as null instead of failing the marshaller.
func finPtr(v float64) *float64 {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return nil
	}
	return &v
}

// ---------------------------------------------------------------------------
// The group request (the group wire form plus the bundle parameters)
// ---------------------------------------------------------------------------

// groupRequest is the body of the group endpoints. The group may be given in
// either of two ways:
//
//   - nested: {"group": <lens.Group wire form>, ...} — the form /set uses;
//   - flat: the wire form of lens.Group at the top level of the body, i.e.
//     {"name":…, "elements":[…], "stop":…, "origin_z_mm":…}, plus the bundle
//     parameters (object_z_mm, field_deg, field_tan, pupil_radius_mm, samples,
//     lambda_nm, ...).
//
// The two forms are mutually exclusive. object_z_mm belongs to the group (its
// design object) and is what the bundle is launched from: "inf" (or a missing
// field) is an object at infinity, a number is a finite object plane.
type groupRequest struct {
	// Flat group wire form.
	Name      string         `json:"name,omitempty"`
	Elements  []lens.Element `json:"elements,omitempty"`
	Stop      *int           `json:"stop,omitempty"`
	OriginZMM *float64       `json:"origin_z_mm,omitempty"`
	Group     *lens.Group    `json:"group,omitempty"`

	// object_z_mm is shared: a number or the string "inf".
	ObjectZMM json.RawMessage `json:"object_z_mm,omitempty"`

	// Bundle parameters of image / afocal / dispersion.
	FieldDeg      float64  `json:"field_deg,omitempty"` // 视场角（度，无穷远物，子午 y 方向）
	FieldTan      *float64 `json:"field_tan,omitempty"` // 无穷远物=场角切线；有限物=物高 mm
	PupilRadiusMM float64  `json:"pupil_radius_mm,omitempty"`
	PupilZMM      float64  `json:"pupil_z_mm,omitempty"`
	Rings         int      `json:"samples,omitempty"` // Hexapolar 环数（0 = 3）
	LambdaNm      float64  `json:"lambda_nm,omitempty"`

	// Image-search options (image / dispersion).
	RangeMM float64 `json:"range_mm,omitempty"`
	Scan    int     `json:"scan,omitempty"`
	TolMM   float64 `json:"tol_mm,omitempty"`

	// Dispersion wavelengths (defaults: F and C).
	LambdaANm float64 `json:"lambda_a_nm,omitempty"`
	LambdaBNm float64 `json:"lambda_b_nm,omitempty"`
}

// parseObjectZ accepts a number or the string "inf" (the convention of
// lens.Group). found is false when the field is absent.
func parseObjectZ(raw json.RawMessage) (z float64, found bool, err error) {
	if len(raw) == 0 || string(raw) == "null" {
		return 0, false, nil
	}
	var f float64
	if err := json.Unmarshal(raw, &f); err == nil {
		return f, true, nil
	}
	var s string
	if err := json.Unmarshal(raw, &s); err != nil {
		return 0, false, errors.New(`object_z_mm 必须是数字或 "inf"`)
	}
	switch s {
	case "inf", "-inf":
		return math.Inf(-1), true, nil
	}
	return 0, false, fmt.Errorf(`object_z_mm %q 无法识别（用数字或 "inf"）`, s)
}

// group resolves the request to a validated group.
func (r *groupRequest) group() (*lens.Group, error) {
	z, hasZ, err := parseObjectZ(r.ObjectZMM)
	if err != nil {
		return nil, err
	}
	var g *lens.Group
	if r.Group != nil {
		if len(r.Elements) > 0 {
			return nil, errors.New("请求里 group 与 elements 只能给一个")
		}
		g = r.Group.Clone()
		if hasZ {
			g.ObjectZMM = z
		}
	} else {
		if len(r.Elements) == 0 {
			return nil, errors.New("请求缺少 elements（或嵌套的 group）")
		}
		if len(r.Elements) > maxGroupElements {
			return nil, fmt.Errorf("elements %d 个，超过上限 %d", len(r.Elements), maxGroupElements)
		}
		g = &lens.Group{Name: r.Name, Elements: r.Elements, ObjectZMM: math.Inf(-1)}
		if r.Stop != nil {
			g.Stop = *r.Stop
		}
		if r.OriginZMM != nil {
			g.OriginZ = *r.OriginZMM
		}
		if hasZ {
			g.ObjectZMM = z
		}
	}
	if err := g.Validate(); err != nil {
		return nil, err
	}
	return g, nil
}

// objectRaw encodes an object position for the ObjectZMM raw field of the
// requests: a number, or the string "inf" for an object at infinity.
func objectRaw(z float64) json.RawMessage {
	if math.IsInf(z, -1) {
		return json.RawMessage(`"inf"`)
	}
	buf, err := json.Marshal(z)
	if err != nil {
		return json.RawMessage(`"inf"`)
	}
	return buf
}

// bundleFieldEcho reports the field of a bundle for the response: the tangent
// of the field angle for an object at infinity, the object height in mm for a
// finite one.
func bundleFieldEcho(b lens.Bundle) float64 {
	if b.Object.Finite {
		return b.Object.HY
	}
	return b.Object.TanY
}

// bundle builds the real-ray bundle of the request for the resolved group.
// The object comes from the group (object_z_mm); the field is the meridional
// tangent of the field angle for an object at infinity, and the object height
// in mm for a finite one.
func (r *groupRequest) bundle(g *lens.Group) (lens.Bundle, error) {
	lambda := r.LambdaNm
	if lambda == 0 {
		lambda = dLineNm
	}
	if math.IsNaN(lambda) || math.IsInf(lambda, 0) || !(lambda > 0) {
		return lens.Bundle{}, fmt.Errorf("lambda_nm %v 必须为正且有限", r.LambdaNm)
	}
	rings := defI(r.Rings, 3)
	if rings < 1 || rings > maxGroupRings {
		return lens.Bundle{}, fmt.Errorf("samples %d 超出 [1, %d]（Hexapolar 环数）", r.Rings, maxGroupRings)
	}
	if math.IsNaN(r.PupilRadiusMM) || math.IsInf(r.PupilRadiusMM, 0) || r.PupilRadiusMM < 0 {
		return lens.Bundle{}, fmt.Errorf("pupil_radius_mm %v 必须是非负有限数", r.PupilRadiusMM)
	}
	if math.IsNaN(r.PupilZMM) || math.IsInf(r.PupilZMM, 0) {
		return lens.Bundle{}, fmt.Errorf("pupil_z_mm %v 必须是有限数", r.PupilZMM)
	}
	b := lens.Bundle{
		Samples: geo.Hexapolar(rings), LambdaNm: lambda,
		Pupil: geo.Pupil{Z: r.PupilZMM, Radius: r.PupilRadiusMM},
	}
	if math.IsInf(g.ObjectZMM, -1) {
		tan := 0.0
		switch {
		case r.FieldTan != nil:
			tan = *r.FieldTan
		default:
			if math.IsNaN(r.FieldDeg) || math.IsInf(r.FieldDeg, 0) {
				return lens.Bundle{}, fmt.Errorf("field_deg %v 必须是有限数", r.FieldDeg)
			}
			if math.Abs(r.FieldDeg) >= 90 {
				return lens.Bundle{}, fmt.Errorf("field_deg %g 超出 (−90, 90)", r.FieldDeg)
			}
			tan = math.Tan(r.FieldDeg * math.Pi / 180)
		}
		if math.IsNaN(tan) || math.IsInf(tan, 0) {
			return lens.Bundle{}, fmt.Errorf("field_tan %v 必须是有限数", tan)
		}
		b.Object = geo.FieldTan(0, tan)
		return b, nil
	}
	// Finite object: the field is the object height in mm (field_tan).
	if r.FieldTan == nil {
		if r.FieldDeg != 0 {
			return lens.Bundle{}, errors.New("有限物距的视场请用 field_tan 给出物高 mm（field_deg 只用于无穷远物）")
		}
		b.Object = geo.FieldHeight(g.ObjectZMM, 0, 0)
		return b, nil
	}
	if math.IsNaN(*r.FieldTan) || math.IsInf(*r.FieldTan, 0) {
		return lens.Bundle{}, errors.New("field_tan 必须是有限数")
	}
	b.Object = geo.FieldHeight(g.ObjectZMM, 0, *r.FieldTan)
	return b, nil
}

// imageOptions validates and returns the focus-search options.
func (r *groupRequest) imageOptions() (lens.ImageOptions, error) {
	if r.Scan < 0 || r.Scan > maxGroupScan {
		return lens.ImageOptions{}, fmt.Errorf("scan %d 超出 [0, %d]", r.Scan, maxGroupScan)
	}
	if math.IsNaN(r.RangeMM) || math.IsInf(r.RangeMM, 0) || r.RangeMM < 0 {
		return lens.ImageOptions{}, fmt.Errorf("range_mm %v 必须是非负有限数", r.RangeMM)
	}
	if math.IsNaN(r.TolMM) || math.IsInf(r.TolMM, 0) || r.TolMM < 0 {
		return lens.ImageOptions{}, fmt.Errorf("tol_mm %v 必须是非负有限数", r.TolMM)
	}
	return lens.ImageOptions{RangeMM: r.RangeMM, Scan: r.Scan, TolMM: r.TolMM}, nil
}

// ---------------------------------------------------------------------------
// Responses
// ---------------------------------------------------------------------------

// groupLayoutSurface is one surface of the layout table.
type groupLayoutSurface struct {
	Element int     `json:"element"`
	Face    string  `json:"face"`    // "front" | "rear"
	Surface int     `json:"surface"` // index in Group.System() (breaks included)
	Vertex  vec3Out `json:"vertex"`
	Axis    vec3Out `json:"axis"`   // local +z (the optical axis of the surface)
	Normal  vec3Out `json:"normal"` // surface normal at the vertex
	Folds   int     `json:"folds"`  // chain folds in front of the surface
}

// groupLayoutResponse is the answer of POST /api/lens/group/layout and of
// `oz lens group layout`: the global placement of every surface plus the exit
// frame (the vertex of the last surface and the outgoing beam direction).
type groupLayoutResponse struct {
	Surfaces    []groupLayoutSurface `json:"surfaces"`
	ExitVertex  vec3Out              `json:"exit_vertex"`
	ExitAxis    vec3Out              `json:"exit_axis"`
	Folds       int                  `json:"folds"` // total number of folds
	NumSurfaces int                  `json:"num_surfaces"`
}

// groupImageResponse is the answer of POST /api/lens/group/image.
type groupImageResponse struct {
	Point       vec3Out  `json:"point"`
	RMSMM       float64  `json:"rms_mm"`
	GeoRadiusMM float64  `json:"geo_radius_mm"`
	Axis        vec3Out  `json:"axis"`
	PlaneO      vec3Out  `json:"plane_o"`
	FocusMM     float64  `json:"focus_mm"`
	Rays        int      `json:"rays"`
	Failed      int      `json:"failed"`
	ChiefHit    vec3Out  `json:"chief_hit"`
	ChiefDir    vec3Out  `json:"chief_dir"`
	LambdaNm    float64  `json:"lambda_nm"`
	FieldTanY   float64  `json:"field_tan_y"`
	ObjectZMM   *float64 `json:"object_z_mm"` // null = object at infinity
	Samples     int      `json:"samples"`     // Hexapolar rings
}

// groupAfocalResponse is the answer of POST /api/lens/group/afocal.
type groupAfocalResponse struct {
	Gamma     float64 `json:"gamma"`
	BeamRatio float64 `json:"beam_ratio"`
	SpreadRad float64 `json:"spread_rad"`
	Rays      int     `json:"rays"`
	EntryDir  vec3Out `json:"entry_dir"`
	ExitDir   vec3Out `json:"exit_dir"`
	AxisIn    vec3Out `json:"axis_in"`
	AxisOut   vec3Out `json:"axis_out"`
	LambdaNm  float64 `json:"lambda_nm"`
	Samples   int     `json:"samples"`
}

// groupDispersionResponse is the answer of POST /api/lens/group/dispersion:
// the axial focus shift FocusMM(a) − FocusMM(b) and both image searches.
type groupDispersionResponse struct {
	ShiftMM   float64            `json:"shift_mm"`
	LambdaANm float64            `json:"lambda_a_nm"`
	LambdaBNm float64            `json:"lambda_b_nm"`
	A         groupImageResponse `json:"a"`
	B         groupImageResponse `json:"b"`
}

// groupSetItem is one parameter edit of /api/lens/group/set.
type groupSetItem struct {
	Param string  `json:"param"` // GroupParam.String() syntax (+ friendly aliases)
	Value float64 `json:"value"`
}

// groupSetRequest is the body of POST /api/lens/group/set: a group, and the
// edits applied to it in order.
type groupSetRequest struct {
	Group *lens.Group    `json:"group"`
	Sets  []groupSetItem `json:"sets"`
}

// groupSetResponse is the updated group.
type groupSetResponse struct {
	Group   *lens.Group `json:"group"`
	Applied int         `json:"applied"`
}

// paraxRequest is the body of POST /api/lens/parax: a prescription plus the
// object position (absent = the prescription's own), the reference field
// tangent and the wavelength.
type paraxRequest struct {
	Prescription lens.Prescription `json:"prescription"`
	ObjectZMM    json.RawMessage   `json:"object_z_mm,omitempty"`
	FieldTan     *float64          `json:"field_tan,omitempty"`
	LambdaNm     float64           `json:"lambda_nm,omitempty"`
}

// paraxResponse is the answer of POST /api/lens/parax and of `oz lens parax`.
// Quantities that do not exist for the system (an afocal EFL, the magnification
// of an infinite object, the F–C shift of an afocal system) are null and the
// note field says why.
type paraxResponse struct {
	LambdaNm  float64  `json:"lambda_nm"`
	ObjectZMM *float64 `json:"object_z_mm"` // null = object at infinity

	Mirrors int     `json:"mirrors"`
	Parity  float64 `json:"parity"`
	Afocal  bool    `json:"afocal"`

	EFLMM      *float64 `json:"efl_mm"`
	FocalLenMM *float64 `json:"focal_length_mm"` // Parity()·EFL, the customary focal length
	FFLMM      *float64 `json:"ffl_mm"`
	BFLMM      *float64 `json:"bfl_mm"`
	ImageZMM   *float64 `json:"image_z_mm"`
	Magnifi    *float64 `json:"magnification"`
	EntrPupilZ *float64 `json:"entrance_pupil_z_mm"`
	EntrPupilR *float64 `json:"entrance_pupil_radius_mm"`
	ExitPupilZ *float64 `json:"exit_pupil_z_mm"`
	ExitPupilR *float64 `json:"exit_pupil_radius_mm"`

	Gamma     *float64 `json:"gamma"`
	GammaNote string   `json:"gamma_note,omitempty"`
	FCShiftMM *float64 `json:"fc_shift_mm"`
	FCNote    string   `json:"fc_note,omitempty"`
}

// ---------------------------------------------------------------------------
// The runs
// ---------------------------------------------------------------------------

// runGroupLayout places every surface of the group in the global frame.
func runGroupLayout(req *groupRequest) (*groupLayoutResponse, error) {
	g, err := req.group()
	if err != nil {
		return nil, err
	}
	lay, err := g.Layout()
	if err != nil {
		return nil, err
	}
	v, axis, err := g.ExitFrame()
	if err != nil {
		return nil, err
	}
	out := &groupLayoutResponse{
		Surfaces:    make([]groupLayoutSurface, 0, len(lay)),
		ExitVertex:  vec3OutOf(v),
		ExitAxis:    vec3OutOf(axis),
		NumSurfaces: len(lay),
	}
	for _, s := range lay {
		if s.Folds > out.Folds {
			out.Folds = s.Folds
		}
		out.Surfaces = append(out.Surfaces, groupLayoutSurface{
			Element: s.Element, Face: s.Face.String(), Surface: s.Surface,
			Vertex: vec3OutOf(s.Vertex), Axis: vec3OutOf(s.Axis), Normal: vec3OutOf(s.Normal),
			Folds: s.Folds,
		})
	}
	return out, nil
}

// groupImageResult maps one lens.ImageResult onto the wire form.
func groupImageResult(res *lens.ImageResult, lambda, fieldTan float64, objectZ float64, rings int) groupImageResponse {
	out := groupImageResponse{
		Point: vec3OutOf(res.Point), RMSMM: res.RMS, GeoRadiusMM: res.GeoRadius,
		Axis: vec3OutOf(res.Axis), PlaneO: vec3OutOf(res.PlaneO), FocusMM: res.FocusMM,
		Rays: res.Rays, Failed: res.Failed,
		ChiefHit: vec3OutOf(res.ChiefHit), ChiefDir: vec3OutOf(res.ChiefDir),
		LambdaNm: lambda, FieldTanY: fieldTan,
		ObjectZMM: finPtr(objectZ), Samples: rings,
	}
	return out
}

// runGroupImage finds the best-focus image of the bundle through the group.
func runGroupImage(req *groupRequest) (*groupImageResponse, error) {
	g, err := req.group()
	if err != nil {
		return nil, err
	}
	b, err := req.bundle(g)
	if err != nil {
		return nil, err
	}
	opt, err := req.imageOptions()
	if err != nil {
		return nil, err
	}
	res, err := g.ImagePosition(b, opt)
	if err != nil {
		return nil, err
	}
	out := groupImageResult(res, b.LambdaNm, bundleFieldEcho(b), g.ObjectZMM, defI(req.Rings, 3))
	for _, v := range []float64{out.RMSMM, out.GeoRadiusMM, out.FocusMM} {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return nil, errors.New("成像结果不是有限数")
		}
	}
	return &out, nil
}

// runGroupAfocal measures the angular magnification of an afocal group.
func runGroupAfocal(req *groupRequest) (*groupAfocalResponse, error) {
	g, err := req.group()
	if err != nil {
		return nil, err
	}
	b, err := req.bundle(g)
	if err != nil {
		return nil, err
	}
	res, err := g.AngularMagnification(b)
	if err != nil {
		return nil, err
	}
	return &groupAfocalResponse{
		Gamma: res.Gamma, BeamRatio: res.BeamRatio, SpreadRad: res.SpreadRad, Rays: res.Rays,
		EntryDir: vec3OutOf(res.EntryDir), ExitDir: vec3OutOf(res.ExitDir),
		AxisIn: vec3OutOf(res.AxisIn), AxisOut: vec3OutOf(res.AxisOut),
		LambdaNm: b.LambdaNm, Samples: defI(req.Rings, 3),
	}, nil
}

// runGroupDispersion measures the axial focus shift between two wavelengths.
func runGroupDispersion(req *groupRequest) (*groupDispersionResponse, error) {
	g, err := req.group()
	if err != nil {
		return nil, err
	}
	b, err := req.bundle(g)
	if err != nil {
		return nil, err
	}
	opt, err := req.imageOptions()
	if err != nil {
		return nil, err
	}
	la := req.LambdaANm
	if la == 0 {
		la = geo.WavelengthF
	}
	lb := req.LambdaBNm
	if lb == 0 {
		lb = geo.WavelengthC
	}
	for _, l := range [...]float64{la, lb} {
		if math.IsNaN(l) || math.IsInf(l, 0) || !(l > 0) {
			return nil, fmt.Errorf("波长 %v 必须为正且有限", l)
		}
	}
	res, err := g.Dispersion(b, la, lb, opt)
	if err != nil {
		return nil, err
	}
	rings := defI(req.Rings, 3)
	return &groupDispersionResponse{
		ShiftMM: res.ShiftMM, LambdaANm: la, LambdaBNm: lb,
		A: groupImageResult(res.A, la, bundleFieldEcho(b), g.ObjectZMM, rings),
		B: groupImageResult(res.B, lb, bundleFieldEcho(b), g.ObjectZMM, rings),
	}, nil
}

// runGroupSet applies the edits to the group in order. The error of
// lens.Group.SetParam is returned as it is (the CLI and the service echo it
// verbatim).
func runGroupSet(req *groupSetRequest) (*groupSetResponse, error) {
	if req.Group == nil {
		return nil, errors.New("请求缺少 group")
	}
	if len(req.Sets) == 0 {
		return nil, errors.New("请求缺少 sets（至少一项 {param, value}）")
	}
	if len(req.Sets) > maxGroupSets {
		return nil, fmt.Errorf("sets %d 项，超过上限 %d", len(req.Sets), maxGroupSets)
	}
	g := req.Group.Clone()
	if err := g.Validate(); err != nil {
		return nil, err
	}
	for i, s := range req.Sets {
		p, err := parseGroupParam(s.Param)
		if err != nil {
			return nil, fmt.Errorf("第 %d 项: %w", i+1, err)
		}
		if math.IsNaN(s.Value) {
			return nil, fmt.Errorf("第 %d 项 %s: 值不是数字", i+1, s.Param)
		}
		if err := g.SetParam(p, s.Value); err != nil {
			return nil, err
		}
	}
	if err := g.Validate(); err != nil {
		return nil, err
	}
	return &groupSetResponse{Group: g, Applied: len(req.Sets)}, nil
}

// runGroupPrescription converts the group into the card (surface by surface)
// form, which fails with a descriptive error for the layouts a card chain
// cannot express (see lens.Group.Prescription).
func runGroupPrescription(req *groupRequest) (*lens.Prescription, error) {
	g, err := req.group()
	if err != nil {
		return nil, err
	}
	return g.Prescription()
}

// runParax computes the first-order quantities of a prescription, mirrors
// included (geo.System.Paraxial and its Parity/Afocal/AngularMagnification
// additions, plus the F–C focus shift of geo.System.AxialChromatic).
func runParax(req *paraxRequest) (*paraxResponse, error) {
	p := &req.Prescription
	if err := p.Validate(); err != nil {
		return nil, err
	}
	lambda := req.LambdaNm
	if lambda == 0 {
		lambda = dLineNm
	}
	if math.IsNaN(lambda) || math.IsInf(lambda, 0) || !(lambda > 0) {
		return nil, fmt.Errorf("lambda_nm %v 必须为正且有限", req.LambdaNm)
	}
	objectZ := p.ObjectZMM
	if z, found, err := parseObjectZ(req.ObjectZMM); err != nil {
		return nil, err
	} else if found {
		objectZ = z
	}
	sys, err := p.System()
	if err != nil {
		return nil, err
	}
	par, err := sys.Paraxial(lambda, objectZ)
	if err != nil {
		return nil, err
	}
	out := &paraxResponse{
		LambdaNm: lambda, ObjectZMM: finPtr(objectZ),
		Mirrors: par.Mirrors, Parity: par.Parity(), Afocal: par.Afocal(),
		EFLMM: finPtr(par.EFL), FocalLenMM: finPtr(par.Parity() * par.EFL),
		FFLMM: finPtr(par.FFL), BFLMM: finPtr(par.BFL), ImageZMM: finPtr(par.ImageZ),
		Magnifi:    finPtr(par.Magnification),
		EntrPupilZ: finPtr(par.EntrancePupilZ), EntrPupilR: finPtr(par.EntrancePupilRadius),
		ExitPupilZ: finPtr(par.ExitPupilZ), ExitPupilR: finPtr(par.ExitPupilRadius),
	}
	// Angular magnification: it exists for an afocal system (any object) and
	// for a finite object of a system with power. A failure is reported in
	// gamma_note instead of failing the report.
	fieldTan := 0.0
	if req.FieldTan != nil {
		fieldTan = *req.FieldTan
	}
	switch {
	case out.Afocal && fieldTan == 0:
		out.GammaNote = "无焦系统：给 --field-tan / field_tan（≠0）才能给出 Γ"
	case out.Afocal || !par.Infinite:
		if math.IsNaN(fieldTan) || math.IsInf(fieldTan, 0) {
			out.GammaNote = "field_tan 必须是有限数"
			break
		}
		if g, err := par.AngularMagnification(fieldTan); err != nil {
			out.GammaNote = err.Error()
		} else {
			out.Gamma = finPtr(g)
		}
	default:
		out.GammaNote = "无穷远物 + 有光焦度的系统没有角放大率（Γ）：给 --object-z 有限物距或换 afocal 系统"
	}
	// F–C focus shift (F line 486.1327 nm minus C line 656.2725 nm). It only
	// exists for a system with an image: an afocal system has no focal plane.
	if par.Afocal() {
		out.FCNote = "无焦系统没有像面，F–C 焦移不适用（用 Γ）"
	} else if shift, err := sys.AxialChromatic(objectZ, geo.WavelengthF, geo.WavelengthC); err != nil {
		out.FCNote = err.Error()
	} else if math.IsNaN(shift) || math.IsInf(shift, 0) {
		out.FCNote = "F–C 焦移非有限（像面在无穷远）"
	} else {
		out.FCShiftMM = finPtr(shift)
	}
	return out, nil
}

// ---------------------------------------------------------------------------
// The GroupParam text syntax
// ---------------------------------------------------------------------------

// parseGroupParam parses one parameter reference of `oz lens group set -p` and
// of the "param" field of POST /api/lens/group/set. The canonical syntax is
// lens.GroupParam.String():
//
//	e0.front.R      e1.rear.c       e0.front.K      e0.front.A4
//	e0.front.poly0  e0.front.q0     e0.front.xy0    e0.front.z1
//	e0.front.semi   e0.rear.semi    e0.thk          e0.gap
//	e0.z            e0.dx           e0.dy           e0.tx
//	e0.ty           e0.tz
//
// The long aliases of the element fields are accepted as well (gap_mm,
// thickness_mm, position_mm, decenter_x_mm, tilt_x_deg, ...), and the faces may
// be written f/n or 1/2 instead of front/rear: e2.tilt_x_deg = e2.tx.
func parseGroupParam(s string) (lens.GroupParam, error) {
	parts := strings.Split(strings.TrimSpace(s), ".")
	if len(parts) < 2 || len(parts) > 3 {
		return lens.GroupParam{}, fmt.Errorf("参数 %q 无法识别（形如 e0.front.R、e2.tilt_x_deg、e1.gap_mm）", s)
	}
	idx, ok := parseElementRef(parts[0])
	if !ok {
		return lens.GroupParam{}, fmt.Errorf("参数 %q: %q 不是元件引用（形如 e0）", s, parts[0])
	}
	if len(parts) == 2 {
		switch norm(parts[1]) {
		case "thk", "thickness", "thicknessmm", "t":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPThickness}, nil
		case "gap", "gapmm", "g":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPGap}, nil
		case "z", "position", "positionmm", "pos":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPPosition}, nil
		case "dx", "decenterx", "decenterxmm":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPDecenterX}, nil
		case "dy", "decentery", "decenterymm":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPDecenterY}, nil
		case "tx", "tiltx", "tiltxdeg":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPTiltX}, nil
		case "ty", "tilty", "tiltydeg":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPTiltY}, nil
		case "tz", "tiltz", "tiltzdeg":
			return lens.GroupParam{Surf: lens.SurfRef{Element: idx}, Kind: lens.GPTiltZ}, nil
		}
		return lens.GroupParam{}, fmt.Errorf("参数 %q: 元件级参数 %q 无法识别（thk/thickness_mm、gap/gap_mm、z/position_mm、dx/dy/decenter_*_mm、tx/ty/tz/tilt_*_deg）", s, parts[1])
	}
	face, ok := parseFaceRef(parts[1])
	if !ok {
		return lens.GroupParam{}, fmt.Errorf("参数 %q: %q 不是面（front|rear）", s, parts[1])
	}
	ref := lens.SurfRef{Element: idx, Face: face}
	key := parts[2]
	switch norm(key) {
	case "r", "radius", "radiusmm":
		return lens.GroupParam{Surf: ref, Kind: lens.GPRadius}, nil
	case "c", "curvature":
		return lens.GroupParam{Surf: ref, Kind: lens.GPCurvature}, nil
	case "k", "conic", "conick":
		return lens.GroupParam{Surf: ref, Kind: lens.GPConic}, nil
	case "semi", "semiaperture", "semiaperturemm":
		return lens.GroupParam{Surf: ref, Kind: lens.GPSemiAperture}, nil
	}
	// The coefficient lists poly<i>, q<i>, xy<i> take the list index directly.
	for _, c := range [...]struct {
		prefix string
		kind   lens.GroupParamKind
	}{
		{"poly", lens.GPPolyTerm}, {"q", lens.GPQCoef}, {"xy", lens.GPXYTerm},
	} {
		rest, ok := strings.CutPrefix(norm(key), c.prefix)
		if !ok {
			continue
		}
		i, err := strconv.Atoi(rest)
		if err != nil {
			return lens.GroupParam{}, fmt.Errorf("参数 %q: %s 列表下标必须是整数", s, c.prefix)
		}
		return lens.GroupParam{Surf: ref, Kind: c.kind, Order: i}, nil
	}
	// A4, A6, ... (the even asphere coefficients) and z1, z2, ... (the Noll
	// index of the Zernike list, so z1 is Coeffs[0]).
	if rest, ok := strings.CutPrefix(norm(key), "a"); ok {
		n, err := strconv.Atoi(rest)
		if err != nil || n < 4 || n%2 != 0 {
			return lens.GroupParam{}, fmt.Errorf("参数 %q: 非球面系数写 A4、A6、…（收到 %q）", s, key)
		}
		return lens.GroupParam{Surf: ref, Kind: lens.GPCoeff, Order: (n - 4) / 2}, nil
	}
	if rest, ok := strings.CutPrefix(norm(key), "z"); ok {
		j, err := strconv.Atoi(rest)
		if err != nil || j < 1 {
			return lens.GroupParam{}, fmt.Errorf("参数 %q: Zernike 系数写 z<j>（Noll 序号 j ≥ 1，收到 %q）", s, key)
		}
		return lens.GroupParam{Surf: ref, Kind: lens.GPZernike, Order: j - 1}, nil
	}
	return lens.GroupParam{}, fmt.Errorf("参数 %q: 面参数 %q 无法识别（R/c/K/A4…/poly0/q0/xy0/z1/semi）", s, key)
}

// norm lowercases a key and drops underscores and hyphens, so that tilt_x_deg,
// tiltX-deg and TX all read as one key.
func norm(s string) string {
	return strings.ToLower(strings.NewReplacer("_", "", "-", "").Replace(s))
}

// parseElementRef parses "e0" / "E12" into an element index.
func parseElementRef(s string) (int, bool) {
	t := strings.ToLower(strings.TrimSpace(s))
	if !strings.HasPrefix(t, "e") {
		return 0, false
	}
	i, err := strconv.Atoi(t[1:])
	if err != nil || i < 0 {
		return 0, false
	}
	return i, true
}

// parseFaceRef parses "front"/"f"/"1" and "rear"/"r"/"2".
func parseFaceRef(s string) (lens.Face, bool) {
	switch norm(s) {
	case "front", "f", "1":
		return lens.Front, true
	case "rear", "r", "2":
		return lens.Rear, true
	}
	return lens.Front, false
}
