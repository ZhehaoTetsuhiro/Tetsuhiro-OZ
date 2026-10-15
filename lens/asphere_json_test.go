package lens

import (
	"encoding/json"
	"math"
	"testing"

	"tetsuhiro-oz/geo"
)

// Round-trip and dispatch tests for the freeform (xy_polynomial) and Zernike
// (zernike) Card descriptions added by the G3/G4 integration.

func TestFreeformCardJSONRoundTrip(t *testing.T) {
	p := biconvexBK7()
	p.Cards[0].XYPoly = &XYCard{
		NormRadiusMM: 12,
		Terms: []XYTermCard{
			{I: 2, J: 0, Coeff: 1e-3},
			{I: 0, J: 2, Coeff: 1.5e-3},
			{I: 1, J: 1, Coeff: 2e-4},
		},
	}
	dir := t.TempDir()
	path := dir + "/ff.json"
	if err := p.SaveJSON(path); err != nil {
		t.Fatal(err)
	}
	q, err := LoadJSON(path)
	if err != nil {
		t.Fatal(err)
	}
	if q.Cards[0].XYPoly == nil {
		t.Fatal("xy_polynomial did not survive the round trip")
	}
	want, err := json.Marshal(p)
	if err != nil {
		t.Fatal(err)
	}
	got, err := json.Marshal(q)
	if err != nil {
		t.Fatal(err)
	}
	if string(want) != string(got) {
		t.Fatalf("round trip is not byte-stable:\n%s\n%s", want, got)
	}

	// The card must build a geo.XYPolynomial and evaluate the same sag.
	// (XYPolynomial stores terms in canonical order, so look the term up by
	// its (i, j) index rather than by position.)
	prof, err := q.Cards[0].profile()
	if err != nil {
		t.Fatal(err)
	}
	ff, ok := prof.(geo.XYPolynomial)
	if !ok {
		t.Fatalf("profile is %T, want geo.XYPolynomial", prof)
	}
	found := false
	for _, term := range ff.Terms {
		if term.I == 2 && term.J == 0 {
			found = true
			if a := p.Cards[0].XYPoly.Terms[0].Coeff; a != term.Coeff {
				t.Fatalf("coefficient round trip: %g != %g", a, term.Coeff)
			}
		}
	}
	if !found {
		t.Fatal("(2,0) term missing from the canonicalized polynomial")
	}
	sys, err := q.System()
	if err != nil {
		t.Fatal(err)
	}
	// A trace through the freeform front surface must produce finite hits.
	res := sys.Trace(geo.NewRay(geo.V(0, 2, -10), geo.V(0, 0, 1), geo.Wavelengthd))
	if !res.OK() {
		t.Fatalf("trace through the freeform card failed: %v", res.Status)
	}

	// Mutual exclusion is enforced.
	bad := biconvexBK7()
	bad.Cards[0].XYPoly = &XYCard{Terms: []XYTermCard{{I: 2, J: 0, Coeff: 1e-3}}}
	bad.Cards[0].PolyTerms = []PolyTermCard{{Power: 4, Coeff: 1e-3}}
	if err := bad.Cards[0].validateAsphere(0); err == nil {
		t.Error("freeform + poly_terms combination was not rejected")
	}
}

func TestZernikeCardJSONRoundTrip(t *testing.T) {
	p := biconvexBK7()
	p.Cards[0].Zernike = &ZernikeCard{NormRadiusMM: 10, Coeffs: []float64{0, 0, 0, 2e-4}}
	dir := t.TempDir()
	path := dir + "/z.json"
	if err := p.SaveJSON(path); err != nil {
		t.Fatal(err)
	}
	q, err := LoadJSON(path)
	if err != nil {
		t.Fatal(err)
	}
	if q.Cards[0].Zernike == nil {
		t.Fatal("zernike did not survive the round trip")
	}
	want, _ := json.Marshal(p)
	got, _ := json.Marshal(q)
	if string(want) != string(got) {
		t.Fatalf("round trip is not byte-stable:\n%s\n%s", want, got)
	}

	prof, err := q.Cards[0].profile()
	if err != nil {
		t.Fatal(err)
	}
	z, ok := prof.(geo.ZernikeSurface)
	if !ok {
		t.Fatalf("profile is %T, want geo.ZernikeSurface", prof)
	}
	if z.Coeffs[3] != 2e-4 {
		t.Fatalf("coefficient round trip: %g", z.Coeffs[3])
	}

	// Orthonormal flag scales the list by the Noll normalisation factors.
	po := biconvexBK7()
	po.Cards[0].Zernike = &ZernikeCard{NormRadiusMM: 10, Coeffs: []float64{0, 0, 0, 1e-4}, Orthonormal: true}
	profo, err := po.Cards[0].profile()
	if err != nil {
		t.Fatal(err)
	}
	zo := profo.(geo.ZernikeSurface)
	if want := 1e-4 * geo.ZernikeNollNorm(4); math.Abs(zo.Coeffs[3]-want) > 1e-18 {
		t.Fatalf("orthonormal scaling: got %g, want %g", zo.Coeffs[3], want)
	}

	// Dispersion-free trace through the Zernike surface must be finite.
	sys, err := q.System()
	if err != nil {
		t.Fatal(err)
	}
	res := sys.Trace(geo.NewRay(geo.V(0, 3, -10), geo.V(0, 0, 1), geo.Wavelengthd))
	if !res.OK() {
		t.Fatalf("trace through the Zernike card failed: %v", res.Status)
	}
}
