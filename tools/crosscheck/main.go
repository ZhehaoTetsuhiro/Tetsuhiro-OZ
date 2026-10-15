// Command crosscheck exports reference results of the OZ kernels for
// independent comparison:
//
//	go run ./tools/crosscheck            > oz_tmm.json   (coating spectra; see crosscheck.py)
//	go run ./tools/crosscheck --lens P.json              (lens first order + marginal ray; see lenscheck.py)
package main

import (
	"encoding/json"
	"fmt"
	"os"

	"tetsuhiro-oz/coating"
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
		if err := runLens(os.Args[2], os.Stdout); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
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
