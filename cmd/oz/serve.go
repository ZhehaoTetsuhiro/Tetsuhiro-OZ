package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"net/http"
	"runtime"
	"time"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// cmdServe runs the OZ HTTP service: a small JSON API over the lens and
// coating packages plus a one-page web front end. It is the programmatic
// face of the design suite (the same calls the CLI makes), meant for local
// use and for wiring examples from the libraries.
func cmdServe(args []string) error {
	fs := flag.NewFlagSet("serve", flag.ExitOnError)
	addr := fs.String("addr", ":1130", "监听地址")
	if err := fs.Parse(args); err != nil {
		return err
	}
	log.Printf("Tetsuhiro-OZ 服务 oz %s http://localhost%s", version, *addr)
	return http.ListenAndServe(*addr, ozMux())
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
	w.WriteHeader(code)
	_, _ = w.Write(buf)
}

func writeAPIErr(w http.ResponseWriter, code int, err error) {
	writeJSON(w, code, map[string]any{"error": err.Error()})
}

// readJSONBody reads a limited JSON body into v.
func readJSONBody(w http.ResponseWriter, r *http.Request, v any, limit int64) error {
	if r.Method != http.MethodPost {
		return fmt.Errorf("POST only")
	}
	body, err := io.ReadAll(http.MaxBytesReader(w, r.Body, limit))
	if err != nil {
		return fmt.Errorf("failed to read body: %w", err)
	}
	if err := json.Unmarshal(body, v); err != nil {
		return fmt.Errorf("invalid JSON: %w", err)
	}
	return nil
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
	Prescription lens.Prescription    `json:"prescription"`
	Spec         lens.ToleranceSpec   `json:"spec"`
	Criteria     lens.ToleranceCriteria `json:"criteria"`
	Runs         int                  `json:"runs"`
	Seed         int64                `json:"seed"`
}

// coatingDesignRequest is the /api/coating/design body.
type coatingDesignRequest struct {
	Template   string  `json:"template"`
	LambdaNm   float64 `json:"lambda_nm"`
	FromNm     float64 `json:"from_nm"`
	ToNm       float64 `json:"to_nm"`
	Pairs      int     `json:"pairs"`
	CutoffNm   float64 `json:"cutoff_nm"`
	Incident   string  `json:"incident"`
	Substrate  string  `json:"substrate"`
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

// ozMux builds the service routes (exported to the test suite via this
// package-level function).
func ozMux() *http.ServeMux {
	mux := http.NewServeMux()

	mux.HandleFunc("/api/version", func(w http.ResponseWriter, r *http.Request) {
		writeJSON(w, http.StatusOK, map[string]any{
			"name": "tetsuhiro-oz", "version": version,
			"go": runtime.Version(), "time": time.Now().Format(time.RFC3339),
		})
	})

	mux.HandleFunc("/api/lens/design", func(w http.ResponseWriter, r *http.Request) {
		var req designRequest
		if err := readJSONBody(w, r, &req, 1<<20); err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
			return
		}
		opts := &lens.DesignOptions{SkipPolish: req.SkipPolish, MaxIter: req.MaxIter}
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
	})

	mux.HandleFunc("/api/lens/eval", func(w http.ResponseWriter, r *http.Request) {
		var p lens.Prescription
		if err := readJSONBody(w, r, &p, 1<<20); err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
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
	})

	mux.HandleFunc("/api/lens/tol", func(w http.ResponseWriter, r *http.Request) {
		var req tolRequest
		if err := readJSONBody(w, r, &req, 4<<20); err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
			return
		}
		rep, err := lens.RunTolerances(&req.Prescription, req.Spec, req.Criteria, &lens.ToleranceOptions{Runs: req.Runs, Seed: req.Seed})
		if err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
			return
		}
		writeJSON(w, http.StatusOK, rep)
	})

	mux.HandleFunc("/api/coating/design", func(w http.ResponseWriter, r *http.Request) {
		var req coatingDesignRequest
		if err := readJSONBody(w, r, &req, 1<<20); err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
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
	})

	mux.HandleFunc("/api/coating/spectrum", func(w http.ResponseWriter, r *http.Request) {
		var req coatingSpectrumRequest
		if err := readJSONBody(w, r, &req, 8<<20); err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
			return
		}
		pol, err := coating.ParsePolarization(req.Pol)
		if err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
			return
		}
		rows, err := req.Stack.SpectrumRange(defF(req.FromNm, 400), defF(req.ToNm, 700), defF(req.StepNm, 5), req.AngleDeg, pol)
		if err != nil {
			writeAPIErr(w, http.StatusBadRequest, err)
			return
		}
		out := make([]map[string]float64, len(rows))
		for i, s := range rows {
			out[i] = map[string]float64{"lambda_nm": s.LambdaNm, "R": s.R, "T": s.T, "A": s.A}
		}
		writeJSON(w, http.StatusOK, map[string]any{"rows": out})
	})

	mux.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		if r.URL.Path != "/" {
			http.NotFound(w, r)
			return
		}
		w.Header().Set("Content-Type", "text/html; charset=utf-8")
		_, _ = io.WriteString(w, indexPage)
	})
	return mux
}

// defStr / defF / defI fill request defaults.
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

const indexPage = `<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8"><title>Tetsuhiro-OZ</title>
<style>
 body{font-family:system-ui,sans-serif;background:#101418;color:#dde3ea;margin:2rem auto;max-width:900px;padding:0 1rem}
 h1{font-weight:600} code{background:#1c232b;padding:.1rem .3rem;border-radius:4px}
 textarea,button{font-family:inherit} textarea{width:100%;height:8rem;background:#161c22;color:#dde3ea;border:1px solid #2a3540;border-radius:6px;padding:.6rem}
 button{background:#2f6feb;color:#fff;border:0;border-radius:6px;padding:.45rem .9rem;cursor:pointer;margin:.4rem .3rem}
 pre{background:#161c22;border:1px solid #2a3540;border-radius:6px;padding:.8rem;overflow:auto;max-height:26rem}
 canvas{background:#0c0f13;border:1px solid #2a3540;border-radius:6px;width:100%;height:220px}
 details{margin:.8rem 0} summary{cursor:pointer;color:#8fb7ff}
</style></head><body>
<h1>Tetsuhiro-OZ <small style="color:#7c8a99">/ v` + version + `</small></h1>
<p>几何光学 · 镜头自主设计 · 镀膜自主设计 —— 本地服务。<code>GET /api/version</code> 查看状态。</p>
<details open><summary>消色差双胶合自主设计 → POST /api/lens/design</summary>
<textarea id="designReq">{"efl_mm":100,"crown":"N-BK7","flint":"F2","semi_aperture_mm":10,"field_deg":0.5}</textarea>
<button onclick="post('/api/lens/design','designReq','designOut')">设计</button>
</details>
<pre id="designOut"></pre>
<details><summary>膜系设计 → POST /api/coating/design</summary>
<textarea id="coatReq">{"template":"broadband-ar","from_nm":420,"to_nm":680,"incident":"Air","substrate":"BK7"}</textarea>
<button onclick="post('/api/coating/design','coatReq','coatOut')">设计</button>
<button onclick="spectrum()">设计并画光谱</button>
</details>
<pre id="coatOut"></pre>
<canvas id="plot" width="860" height="220"></canvas>
<script>
async function post(url, inId, outId){const r=await fetch(url,{method:'POST',headers:{'Content-Type':'application/json'},body:document.getElementById(inId).value});const t=await r.text();try{document.getElementById(outId).textContent=JSON.stringify(JSON.parse(t),null,2);}catch(e){document.getElementById(outId).textContent=t;}return t;}
async function spectrum(){const r=await fetch('/api/coating/design',{method:'POST',headers:{'Content-Type':'application/json'},body:document.getElementById('coatReq').value});const d=JSON.parse(await r.text());document.getElementById('coatOut').textContent=JSON.stringify(d,null,2);if(!d.stack)return;
const s=await fetch('/api/coating/spectrum',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({stack:d.stack,from_nm:400,to_nm:750,step_nm:2,angle_deg:0,pol:'u'})});const sd=JSON.parse(await s.text());
const c=document.getElementById('plot').getContext('2d');c.clearRect(0,0,860,220);c.strokeStyle='#2a3540';c.beginPath();c.moveTo(40,10);c.lineTo(40,200);c.lineTo(850,200);c.stroke();
c.strokeStyle='#5ea1ff';c.beginPath();let first=true;for(const p of sd.rows){const x=40+(p.lambda_nm-400)/(750-400)*810;const y=200-p.R*190;first?(c.moveTo(x,y),first=false):c.lineTo(x,y);}c.stroke();
c.fillStyle='#7c8a99';c.fillText('R(λ), 400–750 nm, 0–100%',46,16);}
</script>
</body></html>
`
