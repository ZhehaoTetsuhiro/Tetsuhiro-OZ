package main

import (
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"net/http"
	"runtime"
	"time"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
	"tetsuhiro-oz/server"
)

// defaultServeAddr is OZ's listen address. The retired WOS binary listened on
// :1120; `oz serve -addr :1120` keeps that port.
const defaultServeAddr = ":1130"

// defaultMaxRunMB is the WOS default of -max-run-mb.
const defaultMaxRunMB = 512

// Work limits of the OZ endpoints. A request above a limit is refused with 400
// before any computation starts, so one request cannot hold the CPU or memory
// for an unbounded time. The GUI stays well inside all three.
const (
	maxSpectrumPoints = 10000 // /api/coating/spectrum: (to-from)/step+1
	maxTolRuns        = 20000 // /api/lens/tol: Monte-Carlo runs
	maxDesignIter     = 10000 // /api/lens/design: DLS iterations
)

// cmdServe runs the single OZ HTTP service. It accepts every option of the
// retired `wos` binary (same names and defaults) and serves one mux:
//
//   - the WOS simulation API under /api/ (server.Server.Handler),
//   - the OZ geometry and coating API (/api/version, /api/lens/*, /api/coating/*),
//   - the GUI at / (WOS wave/quantum tabs and the OZ geometry/coating tab).
//
// Environment variables read by the optics kernel need no flag: WOS_GPU=1
// enables the GPU backend, WOS_MUSA_HELPER and WOS_MUSA_SOCKET locate the MUSA
// helper process.
func cmdServe(args []string) error {
	o, err := parseServeFlags(args)
	if errors.Is(err, flag.ErrHelp) {
		return nil
	}
	if err != nil {
		return err
	}

	announceBackend(o.GPU)
	if o.GenGo != "" {
		if err := runGenGo(o.GenGo, o.GenName); err != nil {
			return fmt.Errorf("生成失败: %w", err)
		}
		return nil
	}
	setElementDirs(o.ExtraDirs)
	report := reloadElements()
	if o.Check {
		if len(report.Errors) > 0 {
			return fmt.Errorf("%d 个元件定义有错误", len(report.Errors))
		}
		return nil
	}
	if !o.NoWatch {
		go watchElements()
	}

	srv := server.New(o.MaxRunMB << 20)
	// 后台监视 GPU：启用但持续回退 CPU 时限频打印警告（同样的状态也会出现在
	// /api/health 的 gpu 字段里）。
	srv.WatchGPU(0)
	log.Printf("Tetsuhiro-OZ 服务 oz %s http://localhost%s", version, o.Addr)
	log.Printf("  几何/膜系 API: /api/lens/*, /api/coating/*   波动/量子 API: /api/simulate, /api/quantum, ...   GUI: /")
	return http.ListenAndServe(o.Addr, newServiceMux(srv))
}

// serveOptions is the parsed `oz serve` command line. The names and defaults
// are the ones the retired `wos` binary used.
type serveOptions struct {
	Addr      string
	MaxRunMB  int64
	ExtraDirs string
	Check     bool
	GPU       bool
	GenGo     string
	GenName   string
	NoWatch   bool
}

// parseServeFlags parses the serve flags without exiting the process, so the
// accepted set can be tested.
func parseServeFlags(args []string) (serveOptions, error) {
	var o serveOptions
	fs := flag.NewFlagSet("serve", flag.ContinueOnError)
	fs.StringVar(&o.Addr, "addr", defaultServeAddr, "监听地址")
	fs.Int64Var(&o.MaxRunMB, "max-run-mb", defaultMaxRunMB, "仿真结果在内存中的保留上限（MB，WOS 同名参数）")
	fs.StringVar(&o.ExtraDirs, "elements", "", "额外的元件定义目录（逗号分隔，优先级高于默认目录）")
	fs.BoolVar(&o.Check, "check-elements", false, "只校验元件定义文件，打印结果后退出")
	fs.BoolVar(&o.GPU, "gpu", false, "启用 GPU（CUDA/cuFFT 或 MUSA/muFFT helper）加速（需以 -tags cuda 或 -tags musa 构建；等价于环境变量 WOS_GPU=1）")
	fs.StringVar(&o.GenGo, "gen-go", "", "把元件定义文件导出为 Go 原生元件源码（打印到 stdout）")
	fs.StringVar(&o.GenName, "gen-name", "", "导出元件的类型名（默认使用定义文件里的 name）")
	fs.BoolVar(&o.NoWatch, "no-elements-watch", false, "关闭元件定义文件的自动重载（默认每 2 秒检查一次）")
	if err := fs.Parse(args); err != nil {
		return o, err
	}
	if fs.NArg() > 0 {
		return o, fmt.Errorf("serve: 多余的参数 %q", fs.Arg(0))
	}
	return o, nil
}

// newServiceMux assembles the one HTTP mux. Path ownership is fixed here, so
// no two components claim the same path:
//
//	/api/version, /api/lens/{design,eval,tol}, /api/coating/{design,spectrum}  OZ (exact paths)
//	/api/...                                                                  WOS API (server.Server.Handler)
//	/oz/                                                                      301 -> /#oz (pre-merge address)
//	/                                                                         embedded GUI (server.WebHandler)
//
// Exact OZ paths are more specific than the /api/ prefix, so they win; every
// other /api/ path goes to the WOS mux, which answers 404 for unknown routes.
func newServiceMux(srv *server.Server) *http.ServeMux {
	mux := http.NewServeMux()
	mux.Handle("/api/", srv.Handler())
	mux.HandleFunc("/api/version", recoverJSON(handleVersion))
	mux.HandleFunc("/api/lens/design", recoverJSON(handleLensDesign))
	mux.HandleFunc("/api/lens/eval", recoverJSON(handleLensEval))
	mux.HandleFunc("/api/lens/tol", recoverJSON(handleLensTol))
	mux.HandleFunc("/api/coating/design", recoverJSON(handleCoatingDesign))
	mux.HandleFunc("/api/coating/spectrum", recoverJSON(handleCoatingSpectrum))
	mux.HandleFunc("/oz/", handleOZRedirect)
	mux.Handle("/", server.WebHandler())
	return mux
}

// ozMux returns the complete service mux with default run limits. Tests use it
// to drive the whole service over httptest.
func ozMux() *http.ServeMux {
	return newServiceMux(server.New(defaultMaxRunMB << 20))
}

// writeJSON marshals before touching the response (NaN/Inf in the payload
// makes encoding/json fail; a half-written 200 would hide that).
func writeJSON(w http.ResponseWriter, code int, v any) {
	buf, err := json.Marshal(v)
	if err != nil {
		code = http.StatusInternalServerError
		buf = []byte(`{"error":"response contains non-finite numbers"}`)
	}
	w.Header().Set("Content-Type", "application/json; charset=utf-8")
	w.Header().Set("Access-Control-Allow-Origin", "*")
	w.WriteHeader(code)
	_, _ = w.Write(buf)
}

func writeAPIErr(w http.ResponseWriter, code int, err error) {
	writeJSON(w, code, map[string]any{"error": err.Error()})
}

// recoverJSON turns a panic inside an OZ handler into a JSON 500 instead of a
// dropped connection.
func recoverJSON(h http.HandlerFunc) http.HandlerFunc {
	return func(w http.ResponseWriter, r *http.Request) {
		defer func() {
			if rec := recover(); rec != nil {
				log.Printf("panic in %s %s: %v", r.Method, r.URL.Path, rec)
				writeAPIErr(w, http.StatusInternalServerError, fmt.Errorf("internal error: %v", rec))
			}
		}()
		h(w, r)
	}
}

// decodePOST reads a limited JSON body into v for a POST-only endpoint. It
// answers 405 for any other method and 400 for an unreadable or malformed
// body, and reports whether the handler may continue.
func decodePOST(w http.ResponseWriter, r *http.Request, v any, limit int64) bool {
	if r.Method != http.MethodPost {
		writeAPIErr(w, http.StatusMethodNotAllowed, errors.New("POST only"))
		return false
	}
	body, err := io.ReadAll(http.MaxBytesReader(w, r.Body, limit))
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("failed to read body: %w", err))
		return false
	}
	if err := json.Unmarshal(body, v); err != nil {
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("invalid JSON: %w", err))
		return false
	}
	return true
}

// designRequest is the /api/lens/design body.
type designRequest struct {
	EFLmm      float64 `json:"efl_mm"`
	Crown      string  `json:"crown"`
	Flint      string  `json:"flint"`
	SemiAperMM float64 `json:"semi_aperture_mm"`
	FieldDeg   float64 `json:"field_deg"`
	SkipPolish bool    `json:"skip_polish"`
	MaxIter    int     `json:"max_iter"`
}

// tolRequest is the /api/lens/tol body.
type tolRequest struct {
	Prescription lens.Prescription      `json:"prescription"`
	Spec         lens.ToleranceSpec     `json:"spec"`
	Criteria     lens.ToleranceCriteria `json:"criteria"`
	Runs         int                    `json:"runs"`
	Seed         int64                  `json:"seed"`
}

// coatingDesignRequest is the /api/coating/design body.
type coatingDesignRequest struct {
	Template  string  `json:"template"`
	LambdaNm  float64 `json:"lambda_nm"`
	FromNm    float64 `json:"from_nm"`
	ToNm      float64 `json:"to_nm"`
	Pairs     int     `json:"pairs"`
	CutoffNm  float64 `json:"cutoff_nm"`
	Incident  string  `json:"incident"`
	Substrate string  `json:"substrate"`
}

// coatingSpectrumRequest is the /api/coating/spectrum body.
type coatingSpectrumRequest struct {
	Stack    coating.Stack `json:"stack"`
	FromNm   float64       `json:"from_nm"`
	ToNm     float64       `json:"to_nm"`
	StepNm   float64       `json:"step_nm"`
	AngleDeg float64       `json:"angle_deg"`
	Pol      string        `json:"pol"`
}

// handleVersion answers GET /api/version. The version comes from the constant
// in main.go, which a test keeps equal to the repository VERSION file.
func handleVersion(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodGet {
		writeAPIErr(w, http.StatusMethodNotAllowed, errors.New("GET only"))
		return
	}
	writeJSON(w, http.StatusOK, map[string]any{
		"name": "tetsuhiro-oz", "version": version,
		"go": runtime.Version(), "time": time.Now().Format(time.RFC3339),
	})
}

// handleLensDesign answers POST /api/lens/design: the autonomous achromat
// design, the same call as `oz design achromat`.
func handleLensDesign(w http.ResponseWriter, r *http.Request) {
	var req designRequest
	if !decodePOST(w, r, &req, 1<<20) {
		return
	}
	if !(req.EFLmm > 0) {
		writeAPIErr(w, http.StatusBadRequest, errors.New("efl_mm 必填且为正"))
		return
	}
	if req.MaxIter > maxDesignIter {
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("max_iter %d exceeds the limit %d", req.MaxIter, maxDesignIter))
		return
	}
	opts := &lens.DesignOptions{SkipPolish: req.SkipPolish, MaxIter: req.MaxIter}
	if req.FieldDeg > 0 {
		// Same field convention as the CLI: on-axis plus the given half-angle.
		opts.Fields = []lens.Field{{TanX: 0, TanY: 0}, {TanX: 0, TanY: math.Tan(req.FieldDeg * math.Pi / 180)}}
	}
	p, res, err := lens.DesignAchromat(req.EFLmm, defStr(req.Crown, "N-BK7"), defStr(req.Flint, "F2"), req.SemiAperMM, opts)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	out := map[string]any{"prescription": p}
	if rep, err := lens.Describe(p); err == nil {
		out["report"] = rep
	}
	if res != nil {
		out["optimizer"] = map[string]any{
			"iterations": res.Iter, "evals": res.Evals,
			"merit": res.Merit, "converged": res.Converged,
		}
	}
	writeJSON(w, http.StatusOK, out)
}

// handleLensEval answers POST /api/lens/eval with the report and first-order
// metrics of a prescription (the body is the prescription itself).
func handleLensEval(w http.ResponseWriter, r *http.Request) {
	var p lens.Prescription
	if !decodePOST(w, r, &p, 1<<20) {
		return
	}
	rep, err := lens.Describe(&p)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	metrics := map[string]any{}
	if par, err := p.Paraxial(geo.Wavelengthd); err == nil {
		metrics["efl_mm"] = par.EFL
		metrics["bfl_mm"] = par.BFL
		metrics["f_number"] = par.FNumber
	}
	if chrom, err := lens.ChromaticResidual(&p, 1.0); err == nil {
		metrics["color_shift_mm"] = chrom
	}
	writeJSON(w, http.StatusOK, map[string]any{"report": rep, "metrics": metrics})
}

// handleLensTol answers POST /api/lens/tol: deterministic Monte-Carlo tolerance
// analysis of a prescription.
func handleLensTol(w http.ResponseWriter, r *http.Request) {
	var req tolRequest
	if !decodePOST(w, r, &req, 4<<20) {
		return
	}
	if req.Runs > maxTolRuns {
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("runs %d exceeds the limit %d", req.Runs, maxTolRuns))
		return
	}
	rep, err := lens.RunTolerances(&req.Prescription, req.Spec, req.Criteria, &lens.ToleranceOptions{Runs: req.Runs, Seed: req.Seed})
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	writeJSON(w, http.StatusOK, rep)
}

// handleCoatingDesign answers POST /api/coating/design with one of the built-in
// coating templates.
func handleCoatingDesign(w http.ResponseWriter, r *http.Request) {
	var req coatingDesignRequest
	if !decodePOST(w, r, &req, 1<<20) {
		return
	}
	inc, ok := coating.LookupMaterial(defStr(req.Incident, "Air"))
	if !ok {
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("unknown incident material %q", req.Incident))
		return
	}
	sub, ok := coating.LookupMaterial(defStr(req.Substrate, "BK7"))
	if !ok {
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("unknown substrate material %q", req.Substrate))
		return
	}
	lambda := defF(req.LambdaNm, 550)
	var (
		stack coating.Stack
		res   *coating.Result
		err   error
	)
	switch req.Template {
	case "single-ar":
		stack, res, err = coating.SingleLayerAR(inc, sub, lambda, coating.Options{})
	case "vcoat":
		stack, res, err = coating.VCoat(inc, sub, lambda, coating.Options{})
	case "broadband-ar":
		stack, res, err = coating.BroadbandAR(inc, sub, defF(req.FromNm, 420), defF(req.ToNm, 680), coating.Options{})
	case "hr":
		stack, res, err = coating.HighReflector(inc, sub, lambda, defI(req.Pairs, 12), coating.Options{})
	case "lwp":
		stack, res, err = coating.LongWavePass(inc, sub, req.CutoffNm, coating.Options{})
	case "swp":
		stack, res, err = coating.ShortWavePass(inc, sub, req.CutoffNm, coating.Options{})
	case "splitter":
		stack, res, err = coating.Splitter50(inc, sub, lambda, coating.Options{})
	default:
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("unknown template %q (have: single-ar, vcoat, broadband-ar, hr, lwp, swp, splitter)", req.Template))
		return
	}
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	out := map[string]any{"stack": stack}
	if res != nil {
		out["summary"] = res.Summary()
		out["merit"] = res.Merit
	}
	writeJSON(w, http.StatusOK, out)
}

// handleCoatingSpectrum answers POST /api/coating/spectrum with the R/T/A table
// of a stack.
func handleCoatingSpectrum(w http.ResponseWriter, r *http.Request) {
	var req coatingSpectrumRequest
	if !decodePOST(w, r, &req, 8<<20) {
		return
	}
	pol, err := coating.ParsePolarization(req.Pol)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	from, to, step := defF(req.FromNm, 400), defF(req.ToNm, 700), defF(req.StepNm, 5)
	if step > 0 && to >= from && (to-from)/step+1 > maxSpectrumPoints {
		writeAPIErr(w, http.StatusBadRequest, fmt.Errorf("spectrum has more than %d points; raise step_nm or narrow the range", maxSpectrumPoints))
		return
	}
	rows, err := req.Stack.SpectrumRange(from, to, step, req.AngleDeg, pol)
	if err != nil {
		writeAPIErr(w, http.StatusBadRequest, err)
		return
	}
	out := make([]map[string]float64, len(rows))
	for i, s := range rows {
		out[i] = map[string]float64{"lambda_nm": s.LambdaNm, "R": s.R, "T": s.T, "A": s.A}
	}
	writeJSON(w, http.StatusOK, map[string]any{"rows": out})
}

// handleOZRedirect keeps the pre-merge address of the geometry/coating page
// working: /oz/ moves permanently to the geometry/coating tab of the GUI.
func handleOZRedirect(w http.ResponseWriter, r *http.Request) {
	if r.URL.Path != "/oz/" {
		http.NotFound(w, r)
		return
	}
	http.Redirect(w, r, "/#oz", http.StatusMovedPermanently)
}

// defStr / defF / defI fill request defaults (zero means "not given").
func defStr(v, d string) string {
	if v == "" {
		return d
	}
	return v
}

func defF(v, d float64) float64 {
	if v == 0 {
		return d
	}
	return v
}

func defI(v, d int) int {
	if v == 0 {
		return d
	}
	return v
}
