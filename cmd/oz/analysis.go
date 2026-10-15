package main

// Shared analysis core of the new CLI subcommands and HTTP endpoints.
//
// `oz lens thermal|sens`, `oz geo mtf|ghost` and `oz doe focus` and the POST
// /api/lens/thermal, /api/lens/sens, /api/geo/mtf, /api/geo/ghost and
// /api/doe/focus endpoints all call the run* functions below: the command line
// and the service cannot drift apart because there is only one implementation.
// Every function is a pure function of its request (no time, no map iteration,
// no randomness), so equal requests give bit-identical answers.
//
// Requests and responses are the JSON forms of the API (docs/API.md §4);
// the CLI maps its flags onto the same structures and prints the response.

import (
	"errors"
	"fmt"
	"math"

	"tetsuhiro-oz/doe"
	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// Work limits of the new endpoints. A request above a limit is refused with 400
// before any computation starts (same policy as maxSpectrumPoints and
// maxTolRuns in serve.go).
const (
	maxThermalDT   = 500.0  // /api/lens/thermal: |delta_t_k|
	maxSensParams  = 128    // /api/lens/sens: sensitivity parameters
	maxMTFFields   = 64     // /api/geo/mtf: field angles
	maxMTFRays     = 200000 // /api/geo/mtf: sampled rays per field
	maxMTFBins     = 4096   // /api/geo/mtf: LSF bins
	maxMTFFreqs    = 512    // /api/geo/mtf: MTF frequencies
	maxDOEGrid     = 2048   // /api/doe/focus: grid size (power of two)
	maxDOELevels   = 4096   // /api/doe/focus: phase levels
	maxGhostPaths  = 512    // /api/geo/ghost: double-reflection paths
	maxGhostRings  = 8      // /api/geo/ghost: sampling rings (Hexapolar)
	maxSensCompIt  = 10000  // /api/lens/sens: compensator iterations
	maxCompRangeMM = 1e6    // /api/lens/sens: compensator half-range
)

// dLineNm mirrors the default wavelength of the commands that take --lambda:
// the helium d-line of the geo package (587.5618 nm).
const dLineNm = 587.5618

// ---------------------------------------------------------------------------
// Thermal drift (G9) — `oz lens thermal`, POST /api/lens/thermal
// ---------------------------------------------------------------------------

// lensThermalRequest is the input of runLensThermal.
type lensThermalRequest struct {
	Prescription lens.Prescription `json:"prescription"`
	DeltaTK      float64           `json:"delta_t_k"`       // temperature change from 20 °C, K (non-zero)
	LambdaNm     float64           `json:"lambda_nm"`       // 0 = d-line
	MountCTE     float64           `json:"mount_cte_per_k"` // spacer CTE, 1/K (0 = rigid mount)
}

// lensThermalResponse is the JSON form of lens.ThermalDrift with the units in
// the field names.
type lensThermalResponse struct {
	DeltaTK   float64 `json:"delta_t_k"`
	LambdaNm  float64 `json:"lambda_nm"`
	MountCTE  float64 `json:"mount_cte_per_k"`
	EFL0MM    float64 `json:"efl0_mm"`      // EFL at 20 °C
	EFL1MM    float64 `json:"efl1_mm"`      // EFL at 20 °C + delta_t_k
	BFL0MM    float64 `json:"bfl0_mm"`      // BFL at 20 °C
	BFL1MM    float64 `json:"bfl1_mm"`      // BFL at 20 °C + delta_t_k
	EFLPerKMM float64 `json:"efl_per_k_mm"` // secant slope (EFL1 − EFL0)/delta_t_k
	BFLPerKMM float64 `json:"bfl_per_k_mm"` // secant slope (BFL1 − BFL0)/delta_t_k
}

// runLensThermal computes the first-order thermal drift (lens.ThermalDriftOf).
func runLensThermal(req lensThermalRequest) (*lensThermalResponse, error) {
	if req.DeltaTK == 0 {
		return nil, errors.New("delta_t_k 必填且非 0（0 K 没有漂移）")
	}
	if math.IsNaN(req.DeltaTK) || math.IsInf(req.DeltaTK, 0) || math.Abs(req.DeltaTK) > maxThermalDT {
		return nil, fmt.Errorf("delta_t_k %v 超出 ±%g K 的范围", req.DeltaTK, maxThermalDT)
	}
	if math.IsNaN(req.MountCTE) || math.IsInf(req.MountCTE, 0) || math.Abs(req.MountCTE) > 1e-2 {
		return nil, fmt.Errorf("mount_cte_per_k %v 不是合理的线膨胀系数（|α| ≤ 0.01 1/K）", req.MountCTE)
	}
	lambda := req.LambdaNm
	if lambda == 0 {
		lambda = dLineNm
	}
	if math.IsNaN(lambda) || math.IsInf(lambda, 0) || lambda <= 0 {
		return nil, fmt.Errorf("lambda_nm %v 必须为正且有限", req.LambdaNm)
	}
	d, err := req.Prescription.ThermalDriftOf(lambda, req.DeltaTK, &lens.ThermalOptions{MountCTE: req.MountCTE})
	if err != nil {
		return nil, err
	}
	out := &lensThermalResponse{
		DeltaTK: d.DT, LambdaNm: lambda, MountCTE: req.MountCTE,
		EFL0MM: d.EFL0, EFL1MM: d.EFL1, BFL0MM: d.BFL0, BFL1MM: d.BFL1,
		EFLPerKMM: d.EFLPerK, BFLPerKMM: d.BFLPerK,
	}
	for _, v := range []float64{out.EFL0MM, out.EFL1MM, out.BFL0MM, out.BFL1MM, out.EFLPerKMM, out.BFLPerKMM} {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return nil, errors.New("thermal drift is not finite for this prescription")
		}
	}
	return out, nil
}

// ---------------------------------------------------------------------------
// Sensitivity, inverse sensitivity and compensator (G6)
// `oz lens sens`, POST /api/lens/sens
// ---------------------------------------------------------------------------

// Sensitivity metric names of lensSensRequest.Metric.
const (
	metricEFL   = "efl"   // ∂EFL/∂x (mm per parameter unit)
	metricMerit = "merit" // ∂ (fixed-plane spot merit)/∂x (merit mm per unit)
)

// lensSensRequest is the input of runLensSens. The zero value of the option
// fields selects the documented defaults.
type lensSensRequest struct {
	Prescription lens.Prescription `json:"prescription"`
	LambdaNm     float64           `json:"lambda_nm"`  // 0 = d-line
	Metric       string            `json:"metric"`     // "", "efl" (default) or "merit"
	RelStep      float64           `json:"rel_step"`   // coarsest relative step (0 = 1e-3)
	FixedStep    bool              `json:"fixed_step"` // single fixed step, no ladder

	// Inverse sensitivity (AllocateTolerances): distribute budget_mm of the
	// spot-merit increment over the parameters.
	Inverse      bool    `json:"inverse"`
	BudgetMM     float64 `json:"budget_mm"`
	Mode         string  `json:"mode"`                 // "rss" (default) or "linear"
	MaxCurvature float64 `json:"max_curvature_per_mm"` // cap, 1/mm (0 = 2e-5)
	MaxThickness float64 `json:"max_thickness_mm"`     // cap, mm  (0 = 0.05)
	MaxIndex     float64 `json:"max_index"`            // cap, Δn  (0 = 5e-4)
	MaxConic     float64 `json:"max_conic"`            // cap, K   (0 = 0.01)
	MaxCoeff     float64 `json:"max_coeff"`            // cap, A_k (0 = 1e-7)

	// Compensator verification (Compensate): perturb, then re-optimise.
	Compensate         bool      `json:"compensate"`
	PerturbRadiusMM    float64   `json:"perturb_radius_mm"`     // +dR on every finite-radius surface
	PerturbThicknessMM float64   `json:"perturb_thickness_mm"`  // +dt on every gap
	PerturbIndex       float64   `json:"perturb_index"`         // +dn of every glass
	CompVars           []string  `json:"compensator_vars"`      // "s1.R", "s1.thk", ...; empty = image shift only
	CompRangesMM       []float64 `json:"compensator_ranges_mm"` // half-ranges, same length (0 = default)
	ImageShift         bool      `json:"image_shift"`           // add the image plane to the compensators
	ShiftRangeMM       float64   `json:"shift_range_mm"`        // half-range of the image shift (0 = 1 mm)
	MaxIter            int       `json:"max_iter"`              // compensator DLS iterations (0 = 100)
}

// sensRow is one line of the sensitivity table.
type sensRow struct {
	Param     string  `json:"param"`
	Deriv     float64 `json:"deriv"`     // ∂f/∂x_i
	Err       float64 `json:"err"`       // uncertainty of Deriv (mm)
	Naive     float64 `json:"naive"`     // plain central difference at the selected step
	RelStep   float64 `json:"rel_step"`  // selected relative step
	Step      float64 `json:"step"`      // selected step in parameter units
	Converged bool    `json:"converged"` // u ≤ tol·|Deriv|
}

// allocRow is one allocated tolerance of the inverse table.
type allocRow struct {
	Param  string  `json:"param"`
	Sens   float64 `json:"sens"`   // ∂M/∂x_i used for the allocation
	Tol    float64 `json:"tol"`    // allocated half-range, parameter units
	Share  float64 `json:"share"`  // budget share |s·t| (linear) or (s·t)² (rss)
	Capped bool    `json:"capped"` // true when the tolerance sits at its limit
}

// allocTable is the inverse-sensitivity result (lens.Allocation).
type allocTable struct {
	Mode      string     `json:"mode"`
	BudgetMM  float64    `json:"budget_mm"`
	Predicted float64    `json:"predicted_mm"`
	RelErr    float64    `json:"rel_err"`
	Rows      []allocRow `json:"rows"`
}

// compTable is the compensator verification result (lens.CompensationResult).
type compTable struct {
	NominalMM            float64  `json:"nominal_mm"`
	CompensatedNominalMM float64  `json:"compensated_nominal_mm"`
	PerturbedMM          float64  `json:"perturbed_mm"`
	CompensatedMM        float64  `json:"compensated_mm"`
	Recovery             *float64 `json:"recovery"` // null when the perturbation does not degrade the merit
	ImageShiftMM         float64  `json:"image_shift_mm"`
	Iterations           int      `json:"iterations"`
	Converged            bool     `json:"converged"`
}

// lensSensResponse is the answer of `oz lens sens` and POST /api/lens/sens.
type lensSensResponse struct {
	Metric   string      `json:"metric"`
	LambdaNm float64     `json:"lambda_nm"`
	Nominal  float64     `json:"nominal"` // f at the unperturbed prescription (mm)
	Evals    int         `json:"evals"`   // evaluations of f used by the sensitivities
	Rows     []sensRow   `json:"rows"`
	Alloc    *allocTable `json:"allocation,omitempty"`
	Comp     *compTable  `json:"compensation,omitempty"`
}

// sensParams returns the perturbation parameters of a prescription in a fixed
// order (card ascending; per card: curvature, thickness, conic, coefficients,
// index), so the table and the allocation are deterministic:
//
//   - every card contributes its vertex curvature c = 1/R (ParamCurvature,
//     also for planes, whose radius is not a usable coordinate);
//   - every card but the last contributes its axial thickness;
//   - a card with a non-zero conic constant contributes it;
//   - every non-zero aspheric coefficient contributes;
//   - every card in front of a glass contributes the index Δn of that glass.
func sensParams(p *lens.Prescription) []lens.Param {
	out := make([]lens.Param, 0, 4*len(p.Cards))
	for i, c := range p.Cards {
		out = append(out, lens.Param{Kind: lens.ParamCurvature, Card: i})
		if i < len(p.Cards)-1 {
			out = append(out, lens.Param{Kind: lens.ParamThickness, Card: i})
		}
		if c.ConicK != 0 {
			out = append(out, lens.Param{Kind: lens.ParamConic, Card: i})
		}
		for k, a := range c.Coeffs {
			if a != 0 {
				out = append(out, lens.Param{Kind: lens.ParamCoeff, Card: i, Order: k})
			}
		}
		if !isAirGlass(c.Glass) {
			out = append(out, lens.Param{Kind: lens.ParamIndex, Card: i})
		}
	}
	return out
}

// isAirGlass reports whether a card glass name denotes air (the set of
// lens.isAirName, which is not exported).
func isAirGlass(name string) bool {
	switch name {
	case "", "air", "Air", "AIR", "vacuum":
		return true
	}
	return false
}

// capOf returns the tolerance cap of a parameter kind (mm, 1/mm, Δn or the
// kind's own unit).
func (req lensSensRequest) capOf(k lens.ParamKind) float64 {
	switch k {
	case lens.ParamCurvature, lens.ParamRadius:
		return defF(req.MaxCurvature, 2e-5)
	case lens.ParamThickness:
		return defF(req.MaxThickness, 0.05)
	case lens.ParamIndex:
		return defF(req.MaxIndex, 5e-4)
	case lens.ParamConic:
		return defF(req.MaxConic, 0.01)
	case lens.ParamCoeff:
		return defF(req.MaxCoeff, 1e-7)
	}
	return math.Inf(1)
}

// runLensSens computes the sensitivity table of a prescription and, on
// request, the inverse-sensitivity allocation or a compensator verification.
func runLensSens(req lensSensRequest) (*lensSensResponse, error) {
	lambda := req.LambdaNm
	if lambda == 0 {
		lambda = dLineNm
	}
	if math.IsNaN(lambda) || math.IsInf(lambda, 0) || lambda <= 0 {
		return nil, fmt.Errorf("lambda_nm %v 必须为正且有限", req.LambdaNm)
	}
	p := &req.Prescription
	if err := p.Validate(); err != nil {
		return nil, err
	}
	metric := req.Metric
	if metric == "" {
		metric = metricEFL
	}
	if metric != metricEFL && metric != metricMerit {
		return nil, fmt.Errorf("未知 metric %q（可用: efl, merit）", req.Metric)
	}
	if req.Inverse || req.Compensate {
		if metric != metricMerit {
			return nil, errors.New("--inverse/--compensate 使用固定像面点列 merit 的灵敏度（metric=merit）")
		}
	}
	params := sensParams(p)
	if len(params) == 0 {
		return nil, errors.New("prescription 没有可扰动的参数")
	}
	if len(params) > maxSensParams {
		return nil, fmt.Errorf("灵敏度参数 %d 个，超过上限 %d（减少卡数/系数）", len(params), maxSensParams)
	}
	var f lens.ScalarFunc
	var merit lens.MeritSpec
	switch metric {
	case metricEFL:
		f = lens.EFLScalar(lambda)
	case metricMerit:
		m, err := lens.MeritSpec{Kind: lens.MeritSpotRMS, LambdaNm: lambda}.Fixed(p)
		if err != nil {
			return nil, fmt.Errorf("固定像面的点列 merit: %w", err)
		}
		merit = m
		f = m.Func()
	}
	sr, err := lens.Sensitivity(p, params, f, lens.SensOptions{RelStep: req.RelStep, Fixed: req.FixedStep})
	if err != nil {
		return nil, err
	}
	out := &lensSensResponse{Metric: metric, LambdaNm: lambda, Nominal: sr.Nominal, Evals: sr.Evals}
	out.Rows = make([]sensRow, len(params))
	for i, q := range params {
		out.Rows[i] = sensRow{
			Param: q.String(), Deriv: sr.Deriv[i], Err: sr.Err[i], Naive: sr.Naive[i],
			RelStep: sr.RelSteps[i], Step: sr.Steps[i], Converged: sr.Converged[i],
		}
	}

	if req.Inverse {
		alloc, err := runSensInverse(p, params, sr, req)
		if err != nil {
			return nil, err
		}
		out.Alloc = alloc
	}
	if req.Compensate {
		comp, err := runSensCompensate(p, merit, req)
		if err != nil {
			return nil, err
		}
		out.Comp = comp
	}
	return out, nil
}

// runSensInverse distributes the merit budget over the parameters in
// proportion to their sensitivities (lens.AllocateTolerances).
func runSensInverse(p *lens.Prescription, params []lens.Param, sr *lens.SensResult, req lensSensRequest) (*allocTable, error) {
	if !(req.BudgetMM > 0) || math.IsInf(req.BudgetMM, 0) {
		return nil, fmt.Errorf("budget_mm %v 必须为正且有限", req.BudgetMM)
	}
	var mode lens.AllocMode
	switch req.Mode {
	case "", "rss":
		mode = lens.AllocRSS
	case "linear":
		mode = lens.AllocLinear
	default:
		return nil, fmt.Errorf("未知 mode %q（可用: rss, linear）", req.Mode)
	}
	caps := make([]float64, len(params))
	for i, q := range params {
		caps[i] = req.capOf(q.Kind)
	}
	alloc, err := lens.AllocateTolerances(params, sr.Deriv, req.BudgetMM, mode, caps)
	if err != nil {
		return nil, err
	}
	out := &allocTable{
		Mode: alloc.Mode.String(), BudgetMM: alloc.Budget,
		Predicted: alloc.Predicted, RelErr: alloc.RelErr,
		Rows: make([]allocRow, len(params)),
	}
	for i, q := range params {
		out.Rows[i] = allocRow{Param: q.String(), Sens: alloc.Sens[i], Tol: alloc.Tol[i], Share: alloc.Share[i], Capped: alloc.Capped[i]}
	}
	return out, nil
}

// runSensCompensate perturbs the prescription and re-optimises the
// compensators (lens.Compensate). The perturbation moves every radius, gap and
// glass index by the given amount; the compensators are the requested
// prescription quantities plus, unless variables are given, the image plane.
func runSensCompensate(p *lens.Prescription, merit lens.MeritSpec, req lensSensRequest) (*compTable, error) {
	pert := lens.Perturbation{}
	for i, c := range p.Cards {
		if req.PerturbRadiusMM != 0 && !math.IsInf(c.RadiusMM, 0) && c.RadiusMM != 0 {
			pert.Params = append(pert.Params, lens.Param{Kind: lens.ParamRadius, Card: i})
			pert.Offsets = append(pert.Offsets, req.PerturbRadiusMM)
		}
		if req.PerturbThicknessMM != 0 && i < len(p.Cards)-1 {
			pert.Params = append(pert.Params, lens.Param{Kind: lens.ParamThickness, Card: i})
			pert.Offsets = append(pert.Offsets, req.PerturbThicknessMM)
		}
		if req.PerturbIndex != 0 && !isAirGlass(c.Glass) {
			pert.Params = append(pert.Params, lens.Param{Kind: lens.ParamIndex, Card: i})
			pert.Offsets = append(pert.Offsets, req.PerturbIndex)
		}
	}
	if len(pert.Params) == 0 {
		return nil, errors.New("--compensate 需要至少一个非零扰动量（--perturb-radius/-thickness/-index）")
	}
	comp := lens.CompensatorSpec{ImageShift: req.ImageShift, ShiftRange: req.ShiftRangeMM, MaxIter: req.MaxIter}
	if req.MaxIter < 0 || req.MaxIter > maxSensCompIt {
		return nil, fmt.Errorf("max_iter %d 超出上限 %d", req.MaxIter, maxSensCompIt)
	}
	if req.ShiftRangeMM < 0 || req.ShiftRangeMM > maxCompRangeMM {
		return nil, fmt.Errorf("shift_range_mm %v 超出上限 %g", req.ShiftRangeMM, maxCompRangeMM)
	}
	for i, s := range req.CompVars {
		v, err := parseVarSpec(s)
		if err != nil {
			return nil, err
		}
		comp.Vars = append(comp.Vars, v)
		if len(req.CompRangesMM) > i {
			comp.VarRange = append(comp.VarRange, req.CompRangesMM[i])
		}
	}
	if len(req.CompRangesMM) > len(req.CompVars) {
		return nil, fmt.Errorf("compensator_ranges_mm 有 %d 项，compensator_vars 只有 %d 项", len(req.CompRangesMM), len(req.CompVars))
	}
	if len(comp.Vars) == 0 {
		comp.ImageShift = true // the default compensator: refocus
	}
	res, err := lens.Compensate(p, pert, merit, comp)
	if err != nil {
		return nil, err
	}
	out := &compTable{
		NominalMM: res.Nominal, CompensatedNominalMM: res.CompensatedNominal,
		PerturbedMM: res.Perturbed, CompensatedMM: res.Compensated,
		ImageShiftMM: res.ImageShift, Converged: res.Opt != nil && res.Opt.Converged,
	}
	if res.Opt != nil {
		out.Iterations = res.Opt.Iter
	}
	if !math.IsNaN(res.Recovery) {
		r := res.Recovery
		out.Recovery = &r
	}
	return out, nil
}

// parseVarSpec parses a compensator variable such as "s1.R", "s1.thk", "s0.K"
// or "s0.A6" into a lens.VarSpec. The labels are the ones VarSpec.String
// prints; the coefficient index is the label's radial power 2k+4 (A4 → 0).
func parseVarSpec(s string) (lens.VarSpec, error) {
	var card int
	var kind string
	if n, err := fmt.Sscanf(s, "s%d.%3s", &card, &kind); n != 2 || err != nil {
		return lens.VarSpec{}, fmt.Errorf("变量 %q 不是 s<卡号>.<R|thk|K|A4|A6|...> 形式", s)
	}
	v := lens.VarSpec{Card: card}
	switch {
	case kind == "R":
		v.Kind = lens.VarRadius
	case kind == "thk":
		v.Kind = lens.VarThickness
	case kind == "K":
		v.Kind = lens.VarConic
	case len(kind) >= 2 && kind[0] == 'A':
		var power int
		if _, err := fmt.Sscanf(kind[1:], "%d", &power); err != nil || power < 4 || power%2 != 0 {
			return lens.VarSpec{}, fmt.Errorf("变量 %q 的系数标号 %q 不是偶数径向幂 ≥ 4（A4, A6, ...）", s, kind)
		}
		v.Kind = lens.VarCoeff
		v.Order = (power - 4) / 2
	default:
		return lens.VarSpec{}, fmt.Errorf("变量 %q 的种类 %q 未知（R, thk, K, A4, A6, ...）", s, kind)
	}
	return v, nil
}

// ---------------------------------------------------------------------------
// Geometric MTF and field curvature (G11) — `oz geo mtf`, POST /api/geo/mtf
// ---------------------------------------------------------------------------

// geoMTFRequest is the input of runGeoMTF.
type geoMTFRequest struct {
	Prescription lens.Prescription `json:"prescription"`
	LambdaNm     float64           `json:"lambda_nm"`     // 0 = d-line
	FieldsDeg    []float64         `json:"fields_deg"`    // meridional field angles (default [0])
	CyclesPerMM  []float64         `json:"cycles_per_mm"` // MTF frequencies (default 10,20,30,40,50)
	Rays         int               `json:"rays"`          // sampled rays per field (default 2000)
	Bins         int               `json:"bins"`          // LSF bins (default 201)
	ImageZMM     float64           `json:"image_z_mm"`    // 0 = paraxial image plane
	RimFrac      float64           `json:"rim_frac"`      // FieldFoci pupil pair radius (default 1e-3)
}

// geoMTFField is the MTF of one field at the requested frequencies.
type geoMTFField struct {
	ThetaDeg      float64   `json:"theta_deg"`
	TanY          float64   `json:"tan_y"`
	RaysHit       int       `json:"rays_hit"`
	RaysTotal     int       `json:"rays_total"`
	CentroidXMM   float64   `json:"centroid_x_mm"`
	CentroidYMM   float64   `json:"centroid_y_mm"`
	RMSRadiusMM   float64   `json:"rms_radius_mm"`
	BinWidthTMM   float64   `json:"bin_width_tangential_mm"`
	BinWidthSMM   float64   `json:"bin_width_sagittal_mm"`
	NyquistTCPM   float64   `json:"nyquist_tangential_per_mm"`
	NyquistSCPM   float64   `json:"nyquist_sagittal_per_mm"`
	MTFTangential []float64 `json:"mtf_tangential"`
	MTFSagittal   []float64 `json:"mtf_sagittal"`
}

// geoFocalRow is the tangential/sagittal focus of one field (lens-independent
// copy of geo.FocalPoint with the field-curvature sags added).
type geoFocalRow struct {
	ThetaDeg float64 `json:"theta_deg"`
	String   string  `json:"string"`
	ZTMM     float64 `json:"z_tangential_mm"`
	ZSMM     float64 `json:"z_sagittal_mm"`
	YTMM     float64 `json:"y_tangential_mm"`
	YSMM     float64 `json:"y_sagittal_mm"`
	SagTMM   float64 `json:"sag_tangential_mm"` // z_T − image plane
	SagSMM   float64 `json:"sag_sagittal_mm"`   // z_S − image plane
	SagPMM   float64 `json:"sag_petzval_mm"`    // (3·S − T)/2, the Petzval surface
}

// geoMTFResponse is the answer of `oz geo mtf` and POST /api/geo/mtf.
type geoMTFResponse struct {
	LambdaNm    float64       `json:"lambda_nm"`
	ImageZMM    float64       `json:"image_z_mm"`
	Rays        int           `json:"rays"`
	Bins        int           `json:"bins"`
	Frequencies []float64     `json:"cycles_per_mm"`
	Fields      []geoMTFField `json:"fields"`
	Foci        []geoFocalRow `json:"foci"`
	FociNote    string        `json:"foci_note,omitempty"`
}

// runGeoMTF computes the geometric MTF of every requested field from real rays
// (rays → LSF → MTF, geo.BinLSF/geo.LSF.MTF) and the tangential/sagittal focal
// surfaces (geo.System.FieldFoci).
func runGeoMTF(req geoMTFRequest) (*geoMTFResponse, error) {
	lambda := req.LambdaNm
	if lambda == 0 {
		lambda = dLineNm
	}
	if math.IsNaN(lambda) || math.IsInf(lambda, 0) || lambda <= 0 {
		return nil, fmt.Errorf("lambda_nm %v 必须为正且有限", req.LambdaNm)
	}
	fields := req.FieldsDeg
	if len(fields) == 0 {
		fields = []float64{0}
	}
	if len(fields) > maxMTFFields {
		return nil, fmt.Errorf("视场 %d 个，超过上限 %d", len(fields), maxMTFFields)
	}
	freqs := req.CyclesPerMM
	if len(freqs) == 0 {
		freqs = []float64{10, 20, 30, 40, 50}
	}
	if len(freqs) > maxMTFFreqs {
		return nil, fmt.Errorf("频率 %d 个，超过上限 %d", len(freqs), maxMTFFreqs)
	}
	for _, f := range freqs {
		if math.IsNaN(f) || math.IsInf(f, 0) || f < 0 {
			return nil, fmt.Errorf("频率 %v 必须是非负有限数（cycles/mm）", f)
		}
	}
	for _, a := range fields {
		if math.IsNaN(a) || math.IsInf(a, 0) || math.Abs(a) >= 90 {
			return nil, fmt.Errorf("视场角 %v 必须是 |θ| < 90 的有限度数", a)
		}
	}
	rays := defI(req.Rays, 2000)
	bins := defI(req.Bins, 201)
	if rays < 4 || rays > maxMTFRays {
		return nil, fmt.Errorf("rays %d 超出 [4, %d]", req.Rays, maxMTFRays)
	}
	if bins < 2 || bins > maxMTFBins {
		return nil, fmt.Errorf("bins %d 超出 [2, %d]", req.Bins, maxMTFBins)
	}
	p := &req.Prescription
	if err := p.Validate(); err != nil {
		return nil, err
	}
	sys, err := p.System()
	if err != nil {
		return nil, err
	}
	par, err := p.Paraxial(lambda)
	if err != nil {
		return nil, fmt.Errorf("一阶分析（对偏心/倾斜面不适用）: %w", err)
	}
	imageZ := req.ImageZMM
	if imageZ == 0 {
		if math.IsInf(par.ImageZ, 0) || math.IsNaN(par.ImageZ) {
			return nil, errors.New("近轴像面不存在，请用 image_z_mm 指定像面")
		}
		imageZ = par.ImageZ
	}
	if math.IsNaN(imageZ) || math.IsInf(imageZ, 0) {
		return nil, fmt.Errorf("image_z_mm %v 不是有限值", req.ImageZMM)
	}

	// Trace every field, then bin each field on its own spot range: the fields
	// sit at different image-plane positions, so a shared range would combine
	// unrelated spot offsets into one coarse LSF.
	traces := make([]fieldTrace, 0, len(fields))
	pts := geo.EqualAreaDisc(rays)
	for _, a := range fields {
		field := geo.FieldTan(0, math.Tan(a*math.Pi/180))
		if !par.Infinite {
			// Finite object: the angle is taken at the entrance pupil, the
			// standard field-angle convention with a finite conjugate.
			field = geo.FieldHeight(par.ObjectZ, 0, math.Tan(a*math.Pi/180)*(par.EntrancePupilZ-par.ObjectZ))
		}
		pup, err := sys.EntrancePupil(field, lambda)
		if err != nil {
			return nil, err
		}
		sd, err := sys.SpotDiagram(field, pup, pts, imageZ, lambda)
		if err != nil {
			return nil, err
		}
		tr := fieldTrace{thetaDeg: a, field: field, sd: sd}
		for _, sp := range sd.Points {
			if sp.Status.OK() {
				tr.xs = append(tr.xs, sp.Hit.X)
				tr.ys = append(tr.ys, sp.Hit.Y)
			}
		}
		if len(tr.xs) == 0 {
			return nil, fmt.Errorf("视场 %g° 没有任何光线到达像面（口径/渐晕？）", a)
		}
		traces = append(traces, tr)
	}

	// Per-field bin ranges: [min, max + pad), pad = a relative sliver of the
	// span (or of the coordinate magnitude for a degenerate spot) so that
	// BinLSF accepts every sample and no ray is dropped.
	out := &geoMTFResponse{
		LambdaNm: lambda, ImageZMM: imageZ, Rays: rays, Bins: bins,
		Frequencies: freqs, Fields: make([]geoMTFField, 0, len(traces)),
	}
	for _, tr := range traces {
		loY, hiY := spotRange(tr.ys)
		loX, hiX := spotRange(tr.xs)
		lsfT, err := geo.BinLSF(tr.ys, loY, hiY, bins)
		if err != nil {
			return nil, err
		}
		lsfS, err := geo.BinLSF(tr.xs, loX, hiX, bins)
		if err != nil {
			return nil, err
		}
		st := tr.sd.Stats(geo.ChiefReference())
		row := geoMTFField{
			ThetaDeg: tr.thetaDeg, TanY: tr.field.TanY,
			RaysHit: len(tr.xs), RaysTotal: len(tr.sd.Points),
			CentroidXMM: tr.sd.Centroid.X, CentroidYMM: tr.sd.Centroid.Y,
			RMSRadiusMM:   st.RMS,
			MTFTangential: lsfT.MTF(freqs), MTFSagittal: lsfS.MTF(freqs),
		}
		row.BinWidthTMM = lsfT.Width
		row.BinWidthSMM = lsfS.Width
		if lsfT.Width > 0 {
			row.NyquistTCPM = 1 / (2 * lsfT.Width)
		}
		if lsfS.Width > 0 {
			row.NyquistSCPM = 1 / (2 * lsfS.Width)
		}
		out.Fields = append(out.Fields, row)
	}

	// Field curvature: FieldFoci needs an object at infinity and a coaxial
	// system; the note says why the rows are missing when they are.
	if !par.Infinite {
		out.FociNote = "视场角取自无限远物距；有限物距的系统没有场曲行（geo.System.FieldFoci 的定义域）"
	} else {
		foci, err := sys.FieldFoci(fields, lambda, defF(req.RimFrac, 1e-3))
		if err != nil {
			out.FociNote = "场曲计算失败: " + err.Error()
		} else {
			for _, fp := range foci {
				pz := 1.5*fp.ZS - 0.5*fp.ZT // (3·S − T)/2
				out.Foci = append(out.Foci, geoFocalRow{
					ThetaDeg: fp.ThetaDeg,
					String:   fmt.Sprintf("θ=%.6g° T=(%.6g, %.6g) S=(%.6g, %.6g)", fp.ThetaDeg, fp.YT, fp.ZT, fp.YS, fp.ZS),
					ZTMM:     fp.ZT, ZSMM: fp.ZS, YTMM: fp.YT, YSMM: fp.YS,
					SagTMM: fp.ZT - imageZ, SagSMM: fp.ZS - imageZ, SagPMM: pz - imageZ,
				})
			}
		}
	}
	return out, nil
}

// fieldTrace is one traced field: the spot diagram and the image-plane hits of
// the successfully traced rays.
type fieldTrace struct {
	thetaDeg float64
	field    geo.Field
	sd       *geo.SpotDiagram
	xs, ys   []float64
}

// spotRange returns the LSF bin range [lo, hi) of one spot along one
// direction. hi is padded by a relative sliver of the span so that BinLSF
// accepts the maximum sample; a degenerate spot (zero extent) gets a small
// finite range around its position.
func spotRange(vals []float64) (lo, hi float64) {
	lo, hi = math.Inf(1), math.Inf(-1)
	for _, v := range vals {
		if v < lo {
			lo = v
		}
		if v > hi {
			hi = v
		}
	}
	if hi <= lo {
		// Every sample at the same coordinate: a finite, tiny range.
		c := lo
		if math.IsInf(c, 0) || math.IsNaN(c) {
			c = 0
		}
		return c - 1e-9, c + 1e-9
	}
	pad := (hi - lo) * 1e-6
	if pad == 0 {
		pad = 1e-12 * math.Max(1, math.Abs(hi))
	}
	return lo, hi + pad
}

// ---------------------------------------------------------------------------
// DOE focus (G8) — `oz doe focus`, POST /api/doe/focus
// ---------------------------------------------------------------------------

// doeFocusRequest is the input of runDOEFocus. Lengths are taken in mm on the
// outside and in µm inside (the unit of the doe package); the optical
// parameters are the ones of doe.FocusSpec with the waist taken from the
// Fresnel number N_F = w0²/(λf) when waist_um is not given.
type doeFocusRequest struct {
	FocalMM       float64 `json:"focal_mm"`       // design focal length (mm, > 0)
	Levels        int     `json:"levels"`         // phase levels N (>= 2)
	LambdaNm      float64 `json:"lambda_nm"`      // wavelength (0 = 632.8 nm)
	WaistUM       float64 `json:"waist_um"`       // 1/e amplitude radius on the DOE (0 = from fresnel_number)
	FresnelNumber float64 `json:"fresnel_number"` // N_F = w0²/(λf) (0 = 32)
	GridSize      int     `json:"grid_size"`      // grid M (power of two, 0 = 1024)
	PixelUM       float64 `json:"pixel_um"`       // DOE pixel pitch (0 = 8·w0/M, the test configuration)
	WindowRadii   float64 `json:"window_radii"`   // efficiency window radius in spot radii (0 = 2.5)
}

// doeFocusResponse is the answer of `oz doe focus` and POST /api/doe/focus.
type doeFocusResponse struct {
	FocalMM          float64 `json:"focal_mm"`
	FocalUM          float64 `json:"focal_um"`
	LambdaNm         float64 `json:"lambda_nm"`
	Levels           int     `json:"levels"`
	GridSize         int     `json:"grid_size"`
	PixelUM          float64 `json:"pixel_um"`
	WaistUM          float64 `json:"waist_um"`
	FresnelNumber    float64 `json:"fresnel_number"`
	SpotRadiusUM     float64 `json:"spot_radius_um"`     // λF/(πw0), 1/e intensity
	Eta              float64 `json:"eta_first_order"`    // simulated first-order efficiency
	EtaAnalytic      float64 `json:"eta_analytic"`       // η_N = [sin(π/N)/(π/N)]²
	EtaRelErr        float64 `json:"eta_rel_err"`        // |η − η_N| / η_N
	AxialPeakUM      float64 `json:"axial_peak_um"`      // axial intensity maximum
	FocusDeviationUM float64 `json:"focus_deviation_um"` // axial peak − design focal length
	FocusDeviationPX float64 `json:"focus_deviation_px"` // in pixels
	PeakOffsetPX     float64 `json:"peak_offset_px"`     // lateral focal-plane peak offset
	WindowFraction   float64 `json:"window_fraction"`    // ideal Gaussian share in the window
	PowerIn          float64 `json:"incident_power"`
}

// runDOEFocus simulates an N-level focusing DOE (doe.SimulateFocus) and
// reports the first-order efficiency and the focus deviation.
func runDOEFocus(req doeFocusRequest) (*doeFocusResponse, error) {
	if !(req.FocalMM > 0) || math.IsInf(req.FocalMM, 0) {
		return nil, fmt.Errorf("focal_mm %v 必须为正且有限", req.FocalMM)
	}
	if req.Levels < 2 || req.Levels > maxDOELevels {
		return nil, fmt.Errorf("levels %d 超出 [2, %d]", req.Levels, maxDOELevels)
	}
	lambdaNm := defF(req.LambdaNm, 632.8)
	if math.IsNaN(lambdaNm) || lambdaNm <= 0 {
		return nil, fmt.Errorf("lambda_nm %v 必须为正", req.LambdaNm)
	}
	grid := defI(req.GridSize, 1024)
	if grid < 16 || grid > maxDOEGrid || grid&(grid-1) != 0 {
		return nil, fmt.Errorf("grid_size %d 必须是 2 的幂且在 [16, %d]", req.GridSize, maxDOEGrid)
	}
	fUM := req.FocalMM * 1000
	lambdaUM := lambdaNm / 1000
	nf := defF(req.FresnelNumber, 32)
	if nf <= 0 || math.IsInf(nf, 0) {
		return nil, fmt.Errorf("fresnel_number %v 必须为正", req.FresnelNumber)
	}
	w0 := req.WaistUM
	if w0 == 0 {
		w0 = math.Sqrt(nf * lambdaUM * fUM)
	}
	dx := req.PixelUM
	if dx == 0 {
		dx = 8 * w0 / float64(grid)
	}
	if !(w0 > 0) || !(dx > 0) {
		return nil, errors.New("waist_um 与 pixel_um 必须为正")
	}
	res, err := doe.SimulateFocus(doe.FocusSpec{
		Lambda: lambdaUM, F: fUM, Levels: req.Levels,
		W0: w0, M: grid, Dx: dx, WinW: defF(req.WindowRadii, 2.5),
	})
	if err != nil {
		return nil, err
	}
	out := &doeFocusResponse{
		FocalMM: req.FocalMM, FocalUM: fUM, LambdaNm: lambdaNm, Levels: req.Levels,
		GridSize: grid, PixelUM: dx, WaistUM: w0,
		FresnelNumber: w0 * w0 / (lambdaUM * fUM),
		SpotRadiusUM:  res.Wf, Eta: res.Eta, EtaAnalytic: res.EtaAnalytic,
		AxialPeakUM: res.AxialZ, FocusDeviationUM: res.AxialZ - fUM,
		FocusDeviationPX: (res.AxialZ - fUM) / dx, PeakOffsetPX: res.PeakOffsetPx,
		WindowFraction: res.WindowFrac, PowerIn: res.PowerIn,
	}
	if res.EtaAnalytic > 0 {
		out.EtaRelErr = math.Abs(res.Eta-res.EtaAnalytic) / res.EtaAnalytic
	}
	for _, v := range []float64{out.Eta, out.EtaAnalytic, out.EtaRelErr, out.SpotRadiusUM,
		out.AxialPeakUM, out.FocusDeviationUM, out.FocusDeviationPX, out.PeakOffsetPX, out.PowerIn} {
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return nil, errors.New("DOE focus result is not finite")
		}
	}
	return out, nil
}
