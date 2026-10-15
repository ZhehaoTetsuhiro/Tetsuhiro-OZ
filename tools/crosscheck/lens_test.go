package main

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// Spec thresholds of the rayoptics cross-check. They are fixed by the
// cross-check specification and are not relaxed to make a case pass.
const (
	tolRelFirstOrder = 1e-6 // EFL and F/#, relative
	tolBFLAbs        = 1e-6 // mm, back focal length
	tolPPAbs         = 1e-6 // mm, principal-plane positions
	tolEPRAbs        = 1e-6 // mm, entrance-pupil radius
	tolMarginalAbs   = 1e-3 // mm, transverse intercept of the marginal ray
)

// defaultRayopticsPython is the interpreter with rayoptics installed. The
// environment variable OZ_RAYOPTICS_PYTHON overrides it.
const defaultRayopticsPython = "/root/miniconda3/bin/python3"

// lensCase is one self-built prescription (no library/geo input).
type lensCase struct {
	name string
	p    *lens.Prescription
}

// lensCases returns the test prescriptions: a biconvex singlet, a plano-convex
// singlet, a cemented BK7/F2 doublet, and the doublet and a singlet with the
// aperture stop in the interior (the entrance pupil is then a virtual image).
// Every case has an object at infinity.
func lensCases() []lensCase {
	inf := math.Inf(-1)
	return []lensCase{
		{"singlet biconvex BK7, stop at front", &lens.Prescription{
			Name: "singlet", ObjectZMM: inf, Stop: 0,
			Cards: []lens.Card{
				{Name: "S1", RadiusMM: 100, ThicknessMM: 5, Glass: "N-BK7", SemiAperMM: 12.5},
				{Name: "S2", RadiusMM: -100, Glass: "", SemiAperMM: 12.5},
			}}},
		{"singlet plano-convex BK7, stop at front", &lens.Prescription{
			Name: "plano-convex", ObjectZMM: inf, Stop: 0,
			Cards: []lens.Card{
				{Name: "S1", RadiusMM: 51.68, ThicknessMM: 4, Glass: "N-BK7", SemiAperMM: 12.5},
				{Name: "S2", RadiusMM: 0, Glass: "", SemiAperMM: 12.5},
			}}},
		{"cemented doublet BK7/F2, stop at front", &lens.Prescription{
			Name: "doublet", ObjectZMM: inf, Stop: 0,
			Cards: []lens.Card{
				{Name: "S1", RadiusMM: 60, ThicknessMM: 6, Glass: "N-BK7", SemiAperMM: 12.5},
				{Name: "S2", RadiusMM: -40, ThicknessMM: 2.5, Glass: "F2", SemiAperMM: 12.5},
				{Name: "S3", RadiusMM: -120, Glass: "", SemiAperMM: 12.5},
			}}},
		{"cemented doublet BK7/F2, stop at cemented surface", &lens.Prescription{
			Name: "doublet-interior-stop", ObjectZMM: inf, Stop: 1,
			Cards: []lens.Card{
				{Name: "S1", RadiusMM: 60, ThicknessMM: 6, Glass: "N-BK7", SemiAperMM: 25},
				{Name: "S2", RadiusMM: -40, ThicknessMM: 2.5, Glass: "F2", SemiAperMM: 12.5},
				{Name: "S3", RadiusMM: -120, Glass: "", SemiAperMM: 25},
			}}},
		{"singlet biconvex BK7, stop at back surface", &lens.Prescription{
			Name: "singlet-interior-stop", ObjectZMM: inf, Stop: 1,
			Cards: []lens.Card{
				{Name: "S1", RadiusMM: 100, ThicknessMM: 5, Glass: "N-BK7", SemiAperMM: 25},
				{Name: "S2", RadiusMM: -100, Glass: "", SemiAperMM: 12.5},
			}}},
	}
}

// paraxialRef is the closed-form Gaussian reference of an object-at-infinity
// prescription.
type paraxialRef struct {
	EFL, BFL, EPR, FNumber float64
}

// paraxialReference computes EFL, BFL, the entrance-pupil radius and F/# from
// the reduced ray-transfer (ABCD) recurrence of Gaussian optics: in the
// reduced variables (y, n·u) a surface of curvature c between indices n and
// n' acts as nu' = nu − y·(n'−n)·c, and a gap t in index n' as
// y' = y + (t/n')·nu'. A parallel marginal ray of unit height at the first
// vertex gives EFL = −1/u_k and BFL = −y_k/u_k (image space is air). The
// entrance pupil is the object-side image of the stop: a parallel ray of
// height EPR meets the stop edge, so EPR = a_stop/|y_stop| for unit input.
//
// The recurrence shares no code with package geo or with rayoptics.
func paraxialReference(p *lens.Prescription) (paraxialRef, error) {
	n := len(p.Cards)
	if p.Stop < 0 || p.Stop >= n {
		return paraxialRef{}, fmt.Errorf("stop %d out of range", p.Stop)
	}
	y, nu := 1.0, 0.0
	nBefore := 1.0
	yStop := math.NaN()
	for i, c := range p.Cards {
		nAfter := 1.0
		if c.Glass != "" {
			g, ok := geo.Glass(c.Glass)
			if !ok {
				return paraxialRef{}, fmt.Errorf("unknown glass %q", c.Glass)
			}
			nAfter = g.N(geo.Wavelengthd)
		}
		curv := 0.0
		if c.RadiusMM != 0 && !math.IsInf(c.RadiusMM, 0) {
			curv = 1 / c.RadiusMM
		}
		if i == p.Stop {
			yStop = y
		}
		nu -= y * (nAfter - nBefore) * curv
		if i < n-1 {
			y += c.ThicknessMM / nAfter * nu
		}
		nBefore = nAfter
	}
	ref := paraxialRef{EFL: -1 / nu, BFL: -y / nu}
	ref.EPR = p.Cards[p.Stop].SemiAperMM / math.Abs(yStop)
	ref.FNumber = ref.EFL / (2 * ref.EPR)
	return ref, nil
}

// rayopticsResult is the RESULT line printed by lenscheck.py.
type rayopticsResult struct {
	Checks map[string]float64 `json:"checks"`
	OZ     map[string]float64 `json:"oz"`
	RO     map[string]float64 `json:"ro"`
	StopH  float64            `json:"stop_height_ro"`
	EPRRef float64            `json:"epr_ref"`
	OK     bool               `json:"ok"`
}

// parseRayopticsResult extracts the RESULT line from the checker output.
func parseRayopticsResult(out string) (rayopticsResult, error) {
	var res rayopticsResult
	for _, line := range strings.Split(out, "\n") {
		if rest, ok := strings.CutPrefix(line, "RESULT "); ok {
			if err := json.Unmarshal([]byte(rest), &res); err != nil {
				return res, fmt.Errorf("parse RESULT: %w", err)
			}
			return res, nil
		}
	}
	return res, fmt.Errorf("no RESULT line in checker output:\n%s", out)
}

func relErr(a, b float64) float64 { return math.Abs(a-b) / math.Abs(b) }

// TestLensCrossCheckRayoptics runs the object-at-infinity first-order data and
// the marginal ray of each self-built prescription through Tetsuhiro-OZ and
// through rayoptics (lenscheck.py), and asserts the spec thresholds. Each
// quantity is also compared with the closed-form ABCD reference.
func TestLensCrossCheckRayoptics(t *testing.T) {
	python := os.Getenv("OZ_RAYOPTICS_PYTHON")
	if python == "" {
		python = defaultRayopticsPython
	}
	if _, err := os.Stat(python); err != nil {
		t.Skipf("rayoptics interpreter %s not found: %v", python, err)
	}
	if out, err := exec.Command(python, "-c", "import rayoptics").CombinedOutput(); err != nil {
		t.Skipf("rayoptics not importable with %s: %v\n%s", python, err, out)
	}
	script, err := filepath.Abs("lenscheck.py")
	if err != nil {
		t.Fatal(err)
	}
	for _, tc := range lensCases() {
		t.Run(tc.name, func(t *testing.T) {
			exp, err := exportLens(tc.p)
			if err != nil {
				t.Fatalf("exportLens: %v", err)
			}
			if exp.MarginalStatus != "ok" {
				t.Fatalf("marginal ray status %q", exp.MarginalStatus)
			}
			ref, err := paraxialReference(tc.p)
			if err != nil {
				t.Fatalf("paraxialReference: %v", err)
			}
			path := filepath.Join(t.TempDir(), "lens.json")
			buf, err := json.MarshalIndent(exp, "", " ")
			if err != nil {
				t.Fatal(err)
			}
			if err := os.WriteFile(path, append(buf, '\n'), 0o644); err != nil {
				t.Fatal(err)
			}
			out, runErr := exec.Command(python, script, path).CombinedOutput()
			res, err := parseRayopticsResult(string(out))
			if err != nil {
				t.Fatalf("checker failed (%v): %v", runErr, err)
			}
			t.Logf("\n%s", out)

			// Tetsuhiro-OZ against the closed form.
			if d := relErr(exp.EFL, ref.EFL); d > tolRelFirstOrder {
				t.Errorf("OZ EFL %.15g vs ABCD %.15g: rel %.3e > %.0e", exp.EFL, ref.EFL, d, tolRelFirstOrder)
			}
			if d := math.Abs(exp.BFL - ref.BFL); d > tolBFLAbs {
				t.Errorf("OZ BFL %.15g vs ABCD %.15g: |d| %.3e > %.0e mm", exp.BFL, ref.BFL, d, tolBFLAbs)
			}
			if d := relErr(exp.FNumber, ref.FNumber); d > tolRelFirstOrder {
				t.Errorf("OZ F/# %.15g vs ABCD %.15g: rel %.3e > %.0e", exp.FNumber, ref.FNumber, d, tolRelFirstOrder)
			}
			if d := math.Abs(exp.EPRadius - ref.EPR); d > tolEPRAbs {
				t.Errorf("OZ EPR %.15g vs ABCD %.15g: |d| %.3e > %.0e mm", exp.EPRadius, ref.EPR, d, tolEPRAbs)
			}

			// rayoptics against the closed form.
			if d := relErr(res.RO["efl"], ref.EFL); d > tolRelFirstOrder {
				t.Errorf("rayoptics EFL vs ABCD: rel %.3e > %.0e", d, tolRelFirstOrder)
			}
			if d := math.Abs(res.RO["bfl"] - ref.BFL); d > tolBFLAbs {
				t.Errorf("rayoptics BFL vs ABCD: |d| %.3e > %.0e mm", d, tolBFLAbs)
			}
			if d := relErr(res.RO["fno"], ref.FNumber); d > tolRelFirstOrder {
				t.Errorf("rayoptics F/# vs ABCD: rel %.3e > %.0e", d, tolRelFirstOrder)
			}

			// OZ against rayoptics, recomputed here from the printed values.
			vertexLast := exp.VertexLastZ
			if d := relErr(exp.EFL, res.RO["efl"]); d > tolRelFirstOrder {
				t.Errorf("EFL OZ vs rayoptics: rel %.3e > %.0e", d, tolRelFirstOrder)
			}
			if d := relErr(exp.FNumber, res.RO["fno"]); d > tolRelFirstOrder {
				t.Errorf("F/# OZ vs rayoptics: rel %.3e > %.0e", d, tolRelFirstOrder)
			}
			if d := math.Abs(exp.BFL - res.RO["bfl"]); d > tolBFLAbs {
				t.Errorf("BFL OZ vs rayoptics: |d| %.3e > %.0e mm", d, tolBFLAbs)
			}
			pp := math.Max(math.Abs(exp.FrontPrincipal-res.RO["fpp"]),
				math.Abs(exp.BackPrincipal-vertexLast-res.RO["bpp"]))
			if pp > tolPPAbs {
				t.Errorf("principal planes OZ vs rayoptics: |d| %.3e > %.0e mm", pp, tolPPAbs)
			}
			stopSemi := tc.p.Cards[tc.p.Stop].SemiAperMM
			if tc.p.Stop == 0 {
				if d := math.Abs(exp.EPRadius - stopSemi); d > tolEPRAbs {
					t.Errorf("OZ EPR %.12g vs physical stop semi-aperture %.12g: |d| %.3e", exp.EPRadius, stopSemi, d)
				}
			} else {
				if d := math.Abs(res.StopH - stopSemi); d > tolEPRAbs {
					t.Errorf("rayoptics paraxial stop height %.12g vs stop semi-aperture %.12g: |d| %.3e > %.0e mm",
						res.StopH, stopSemi, d, tolEPRAbs)
				}
			}
			if d := math.Abs(exp.MarginalY - res.RO["marg"]); d > tolMarginalAbs {
				t.Errorf("marginal intercept OZ %.12g vs rayoptics %.12g: |d| %.3e > %.0e mm",
					exp.MarginalY, res.RO["marg"], d, tolMarginalAbs)
			}
			if !res.OK {
				t.Errorf("lenscheck.py verdict FAIL")
			}
			t.Logf("EFL %.12g mm, BFL %.12g mm, F/# %.12g, EPR %.12g mm, marginal y %.12g mm (ABCD: EFL %.12g, BFL %.12g, F/# %.12g, EPR %.12g)",
				exp.EFL, exp.BFL, exp.FNumber, exp.EPRadius, exp.MarginalY, ref.EFL, ref.BFL, ref.FNumber, ref.EPR)
		})
	}
}

// TestExportLensRejectsFiniteObject guards the cross-check contract: only an
// object at infinity has a defined EFL/BFL/F# comparison.
func TestExportLensRejectsFiniteObject(t *testing.T) {
	c := lensCases()[0].p
	c.ObjectZMM = -1000
	if _, err := exportLens(c); err == nil {
		t.Fatal("finite object accepted; want an error")
	}
}
