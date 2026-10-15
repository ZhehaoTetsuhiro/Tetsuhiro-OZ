// Command crosscheck exports reference results of the OZ kernels for
// independent comparison:
//
//	go run ./tools/crosscheck            > oz_tmm.json   (coating spectra; see crosscheck.py)
//	go run ./tools/crosscheck --lens P.json              (lens first order + marginal ray; see lenscheck.py)
package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

type layerSpec struct {
	Mat string  `json:"mat"`
	D   float64 `json:"d_nm"`
}

type point struct {
	Lam   float64 `json:"lam_nm"`
	Angle float64 `json:"angle_deg"`
	Pol   string  `json:"pol"`
	R     float64 `json:"R"`
	T     float64 `json:"T"`
	A     float64 `json:"A"`
}

type kase struct {
	Name      string      `json:"name"`
	Incident  string      `json:"incident"`
	Substrate string      `json:"substrate"`
	Layers    []layerSpec `json:"layers"`
	Points    []point     `json:"points"`
}

func lookup(n string) coating.Material {
	m, ok := coating.LookupMaterial(n)
	if !ok {
		fmt.Fprintln(os.Stderr, "crosscheck: missing material", n)
		os.Exit(1)
	}
	return m
}

func grid(from, to, step float64) []float64 {
	var out []float64
	for l := from; l <= to+1e-9; l += step {
		out = append(out, l)
	}
	return out
}

func main() {
	if len(os.Args) >= 3 && os.Args[1] == "--lens" {
		runLens(os.Args[2])
		return
	}
	names := []string{"Air", "MgF2", "SiO2", "TiO2", "Ta2O5", "HfO2", "BK7", "FusedSilica"}
	materials := map[string][][]float64{}
	for _, name := range names {
		m := lookup(name)
		for _, l := range grid(400, 1200, 5) {
			idx := m.IndexNm(l)
			materials[name] = append(materials[name], []float64{l, real(idx), imag(idx)})
		}
	}

	mkstack := func(inc, sub string, ls []layerSpec) coating.Stack {
		layers := make([]coating.Layer, len(ls))
		for i, l := range ls {
			layers[i] = coating.Layer{Material: lookup(l.Mat), ThicknessNm: l.D}
		}
		return coating.NewStack(lookup(inc), lookup(sub), layers...)
	}
	var cases []kase
	add := func(name, inc, sub string, ls []layerSpec, lambdas, angles []float64) {
		st := mkstack(inc, sub, ls)
		k := kase{Name: name, Incident: inc, Substrate: sub, Layers: ls}
		for _, l := range lambdas {
			for _, a := range angles {
				for _, p := range []coating.Polarization{coating.PolS, coating.PolP} {
					s, err := st.Evaluate(l, a, p)
					if err != nil {
						fmt.Fprintln(os.Stderr, "crosscheck:", err)
						os.Exit(1)
					}
					k.Points = append(k.Points, point{l, a, p.String(), s.R, s.T, s.A})
				}
			}
		}
		cases = append(cases, k)
	}

	// Case 1: a single MgF2 quarter-wave-ish layer on BK7, multiple angles.
	add("single MgF2 100nm on BK7", "Air", "BK7",
		[]layerSpec{{"MgF2", 100}}, grid(400, 700, 10), []float64{0, 30, 45, 60})

	// Case 2: a 12-layer quarter-wave stack (TiO2/SiO2) for 550 nm on BK7.
	nH := lookup("TiO2").N(550)
	nL := lookup("SiO2").N(550)
	var qw []layerSpec
	for i := 0; i < 12; i++ {
		if i%2 == 0 {
			qw = append(qw, layerSpec{"TiO2", 550 / (4 * nH)})
		} else {
			qw = append(qw, layerSpec{"SiO2", 550 / (4 * nL)})
		}
	}
	add("12-layer QW stack TiO2/SiO2 @550nm on BK7", "Air", "BK7",
		qw, grid(400, 700, 5), []float64{0, 45})

	// Case 3: an arbitrary five-layer stack on fused silica, oblique angles.
	add("5-layer arbitrary on FusedSilica", "Air", "FusedSilica",
		[]layerSpec{{"Ta2O5", 83}, {"SiO2", 121}, {"Ta2O5", 41}, {"SiO2", 200}, {"HfO2", 63}},
		grid(450, 650, 25), []float64{0, 30, 60})

	enc := json.NewEncoder(os.Stdout)
	enc.SetIndent("", " ")
	if err := enc.Encode(map[string]any{"materials": materials, "cases": cases}); err != nil {
		fmt.Fprintln(os.Stderr, "crosscheck:", err)
		os.Exit(1)
	}
}

// runLens prints the first-order data of a prescription plus a marginal-ray
// reference trace, for the rayoptics cross-check (lenscheck.py).
func runLens(path string) {
	p, err := lens.LoadJSON(path)
	if err != nil {
		fmt.Fprintln(os.Stderr, "lenscheck:", err)
		os.Exit(1)
	}
	sys, err := p.System()
	if err != nil {
		fmt.Fprintln(os.Stderr, "lenscheck:", err)
		os.Exit(1)
	}
	// Cross-check configuration: a "pseudo-infinite" finite object at 1e6 mm
	// in front of surface 0. Both sides (this tool and the rayoptics model)
	// use exactly this distance, so every compared quantity is defined by
	// the same geometry; 1e6 mm keeps the infinity approximation at the 1e-6
	// relative level while staying inside the numerically safe range of the
	// rayoptics surface intersection.
	const objDist = 1e6
	pc := p.Clone()
	pc.ObjectZMM = -objDist
	par, err := pc.Paraxial(geo.Wavelengthd)
	if err != nil {
		fmt.Fprintln(os.Stderr, "lenscheck:", err)
		os.Exit(1)
	}
	type surf struct {
		Name      string  `json:"name"`
		Radius    float64 `json:"radius_mm"`
		Thickness float64 `json:"thickness_mm"`
		Nd        float64 `json:"nd"`
		Vd        float64 `json:"vd"`
		Semi      float64 `json:"semi_mm"`
	}
	surfs := make([]surf, len(p.Cards))
	for i, c := range p.Cards {
		nd, vd := 1.0, 0.0
		if c.Glass != "" {
			g, ok := geo.Glass(c.Glass)
			if !ok {
				fmt.Fprintln(os.Stderr, "lenscheck: unknown glass", c.Glass)
				os.Exit(1)
			}
			nd = g.N(geo.Wavelengthd)
			vd = g.AbbeNumber()
		}
		th := 0.0
		if i < len(p.Cards)-1 {
			th = c.ThicknessMM
		}
		surfs[i] = surf{Name: c.Name, Radius: c.RadiusMM, Thickness: th, Nd: nd, Vd: vd, Semi: c.SemiAperMM}
	}
	pup, err := sys.EntrancePupil(geo.Field{}, geo.Wavelengthd)
	if err != nil {
		fmt.Fprintln(os.Stderr, "lenscheck:", err)
		os.Exit(1)
	}
	// Marginal ray, defined exactly as rayoptics defines it: the real ray
	// that starts at the on-axis object point (z = -objDist) and passes
	// through the rim of the ENTRANCE PUPIL. The stop rim is knife-edge at
	// fraction 1.0, so fall back to 0.999 and export the fraction actually
	// used (the python side traces the same pupil fraction).
	traceMarginal := func(frac float64) (float64, string) {
		objPt := geo.V(0, 0, -objDist)
		aimPt := geo.V(0, pup.Radius*frac, pup.Z)
		dir := aimPt.Sub(objPt).Normalized()
		if dir.IsNaN() {
			return math.NaN(), "bad-aim"
		}
		res := sys.Trace(geo.NewRay(objPt, dir, geo.Wavelengthd))
		if !res.OK() {
			return math.NaN(), res.Status.String()
		}
		zf := par.VertexZ[len(par.VertexZ)-1] + par.BFL
		d := res.Ray.Dir
		if d.Z == 0 {
			return math.NaN(), "escape"
		}
		t := (zf - res.Ray.Origin.Z) / d.Z
		return res.Ray.Origin.Y + d.Y*t, "ok"
	}
	fracUsed := 1.0
	marginalY, mStatus := traceMarginal(fracUsed)
	if mStatus != "ok" {
		fracUsed = 0.999
		marginalY, mStatus = traceMarginal(fracUsed)
	}
	hUsed := pup.Radius * fracUsed
	if math.IsNaN(marginalY) {
		marginalY = 0 // keep the JSON encodable; status carries the truth
	}
	objZ := any(p.ObjectZMM)
	if math.IsInf(p.ObjectZMM, -1) {
		objZ = "inf"
	}
	out := map[string]any{
		"surfaces":                surfs,
		"object_z_mm":             objZ,
		"object_distance_mm":      objDist,
		"wavelength_nm":           geo.Wavelengthd,
		"efl_d_mm":                par.EFL,
		"bfl_d_mm":                par.BFL,
		"fno_d":                   par.FNumber,
		"pupil_radius":            pup.Radius,
		"pupil_z_mm":              pup.Z,
		"marginal_radius_mm":      hUsed,
		"marginal_pupil_fraction": fracUsed,
		"marginal_status":         mStatus,
		"marginal_y_at_focus_mm":  marginalY,
	}
	enc := json.NewEncoder(os.Stdout)
	enc.SetIndent("", " ")
	if err := enc.Encode(out); err != nil {
		fmt.Fprintln(os.Stderr, "lenscheck:", err)
		os.Exit(1)
	}
}
