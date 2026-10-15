// Package geo_test drives the mirror examples of the built-in library through
// the public API (lens.LoadJSON → Prescription.System → geo.Paraxial) and
// cross-checks the first-order image against the real ray trace, which is the
// interface a user of the library sees. The FOCUSZ numbers stored in the
// descriptions are the measured best focus of the parallel fan and are frozen
// here as regression pins; they used to be artifacts of a ray-pair crossing
// sign error in tools/libcheck/measure that was fixed on 2026-10-15 (see the
// in-package TestLibraryMirrorFocusNumbers).
package geo_test

import (
	"math"
	"os"
	"strconv"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// libraryMirror loads one library example relative to the geo package directory.
func libraryMirror(t *testing.T, name string) *lens.Prescription {
	t.Helper()
	path := "../library/geo/" + name
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("library example %s: %v", path, err)
	}
	p, err := lens.LoadJSON(path)
	if err != nil {
		t.Fatalf("LoadJSON(%s): %v", name, err)
	}
	return p
}

// storedFocusZ reads the "FOCUSZ=<number>mm" clause of a library description.
func storedFocusZ(t *testing.T, p *lens.Prescription) float64 {
	t.Helper()
	const key = "FOCUSZ="
	i := strings.Index(p.Description, key)
	if i < 0 {
		t.Fatalf("description of %q has no FOCUSZ clause", p.Name)
	}
	rest := p.Description[i+len(key):]
	j := strings.Index(rest, "mm")
	if j < 0 {
		t.Fatalf("description of %q has no unit after FOCUSZ", p.Name)
	}
	v, err := strconv.ParseFloat(strings.TrimSpace(rest[:j]), 64)
	if err != nil {
		t.Fatalf("FOCUSZ of %q is not a number: %v", p.Name, err)
	}
	return v
}

// axisCrossingZ traces a parallel ray of height h and returns its axis crossing.
func axisCrossingZ(sys *geo.System, h, lambda float64) float64 {
	start := sys.VertexZ(0) - 2*1000 - 200
	res := sys.Trace(geo.NewRay(geo.V(0, h, start), geo.V(0, 0, 1), lambda))
	if !res.OK() {
		return math.NaN()
	}
	p, d := res.Ray.Origin, res.Ray.Dir
	if d.Y == 0 {
		return math.Inf(1)
	}
	return p.Z - p.Y*d.Z/d.Y
}

func checkRelT(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	scale := math.Max(1, math.Abs(want))
	if math.Abs(got-want) > tol*scale {
		t.Errorf("%s: got %.12g, want %.12g (tol %.3g relative)", name, got, want, tol)
	}
}

// TestLibraryMirrorExamplesParaxial checks the two folded mirror layouts of the
// library: the parabolic primary of the Newtonian example and the two-mirror
// example named "Cassegrain".
func TestLibraryMirrorExamplesParaxial(t *testing.T) {
	for _, tc := range []struct {
		file       string
		wantZ      float64
		wantFocusZ float64
		mirrors    int
		vertices   []float64
	}{
		{"mirror-newtonian-f1000-f5.json", -940, -940, 1, []float64{0, 60}},
		{"mirror-parabola-f1000-f8.json", -940, -940, 1, []float64{0, 60}},
		{"mirror-cassegrain-f1000-d300-b300.json", 360, 357.097475, 2, []float64{0, 60, -240}},
	} {
		p := libraryMirror(t, tc.file)
		sys, err := p.System()
		if err != nil {
			t.Fatalf("%s: System: %v", tc.file, err)
		}
		if sys.NumSurfaces() != len(tc.vertices) {
			t.Fatalf("%s: %d surfaces, want %d", tc.file, sys.NumSurfaces(), len(tc.vertices))
		}
		for i, z := range tc.vertices {
			if sys.VertexZ(i) != z {
				t.Errorf("%s: vertex %d at z=%g, want %g (the fold is a negative thickness)", tc.file, i, sys.VertexZ(i), z)
			}
		}
		par, err := sys.Paraxial(geo.Wavelengthd, p.ObjectZMM)
		if err != nil {
			t.Fatalf("%s: Paraxial: %v", tc.file, err)
		}
		if par.Mirrors != tc.mirrors {
			t.Errorf("%s: Mirrors = %d, want %d", tc.file, par.Mirrors, tc.mirrors)
		}
		checkRelT(t, tc.file+" paraxial image", par.ImageZ, tc.wantZ, 1e-12)
		// Real rays: the image of a parallel bundle of vanishing aperture.
		got := axisCrossingZ(sys, 1e-4, geo.Wavelengthd)
		checkRelT(t, tc.file+" real ray image", got, tc.wantZ, 1e-9)
		// The FOCUSZ stored in the description is the measured best focus of
		// the parallel fan (tools/libcheck/measure.focusScan). Since the
		// crossZ sign fix it coincides with the physical focus within the
		// design tolerance; the values are frozen here as regression pins.
		stored := storedFocusZ(t, p)
		checkRelT(t, tc.file+" stored FOCUSZ (frozen)", stored, tc.wantFocusZ, 1e-6)
		if math.Abs(stored-par.ImageZ) > 30 {
			t.Errorf("%s: stored FOCUSZ=%g is not within the design tolerance of the physical focus %g", tc.file, stored, par.ImageZ)
		}
		t.Logf("%s: paraxial ImageZ %.9f, real ray %.9f, stored FOCUSZ %.6f (measured best focus)", tc.file, par.ImageZ, got, stored)
	}
}

// TestLibraryMirrorPrescriptionParaxial checks the prescription level entry
// point of a mirror example and its angular magnification (both examples are
// telescopes: the collimated input stays collimated only for the afocal ones,
// so the check is the consistency of the definition instead of a value).
func TestLibraryMirrorPrescriptionParaxial(t *testing.T) {
	p := libraryMirror(t, "mirror-newtonian-f1000-f5.json")
	par, err := p.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatalf("Prescription.Paraxial: %v", err)
	}
	if !math.IsInf(par.EFL, 0) {
		// The Newtonian example is a single mirror: it has power.
		checkRelT(t, "newtonian EFL = R/2", par.EFL, -1000, 1e-12)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	// A finite object: the prescription's ObjectZMM drives the same engine, and
	// the image follows the mirror equation s' = sR/(R+2s) measured from the
	// vertex of the primary (z = 60).
	src := p.Clone()
	src.ObjectZMM = -3000
	par2, err := src.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatalf("finite object: %v", err)
	}
	s := 3000 + sys.VertexZ(1)
	checkRelT(t, "newtonian finite conjugate", par2.ImageZ, sys.VertexZ(1)+s*(-2000)/(-2000+2*s), 1e-12)
}
