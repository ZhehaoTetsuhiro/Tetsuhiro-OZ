package lens

import (
	"encoding/json"
	"math"
	"reflect"
	"strings"
	"testing"

	"tetsuhiro-oz/geo"
)

// legacyPrescriptionJSON is a prescription as written before the decentre and
// tilt fields existed. It must load unchanged, with every new field zero.
const legacyPrescriptionJSON = `{"name":"legacy biconvex","object_z_mm":"inf","stop":0,"cards":[
  {"name":"front","radius_mm":100,"thickness_mm":6,"glass":"N-BK7","semi_aperture_mm":12.7},
  {"name":"rear","radius_mm":-100,"thickness_mm":95,"semi_aperture_mm":12.7}]}`

// TestCardPoseLegacyJSONUnchanged checks backward compatibility: an old file
// loads as an undecentred, untilted prescription with identical first-order
// numbers, and writing it back produces no decentre/tilt keys.
func TestCardPoseLegacyJSONUnchanged(t *testing.T) {
	var old Prescription
	if err := json.Unmarshal([]byte(legacyPrescriptionJSON), &old); err != nil {
		t.Fatal(err)
	}
	for i, c := range old.Cards {
		if c.posed() || c.DecenterXMM != 0 || c.DecenterYMM != 0 || c.TiltXDeg != 0 || c.TiltYDeg != 0 || c.TiltZDeg != 0 {
			t.Errorf("card %d of a legacy file is posed: %+v", i, c)
		}
	}
	sys, err := old.System()
	if err != nil {
		t.Fatal(err)
	}
	for i, sf := range sys.Surfaces {
		if sf.Pose != nil && !sf.Pose.IsIdentity() {
			t.Errorf("surface %d of a legacy file has pose %+v", i, *sf.Pose)
		}
	}
	par, err := old.Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	ref, err := biconvexBK7().Paraxial(geo.Wavelengthd)
	if err != nil {
		t.Fatal(err)
	}
	if par.EFL != ref.EFL || par.BFL != ref.BFL {
		t.Errorf("legacy EFL/BFL %v/%v differ from the Go-built prescription %v/%v", par.EFL, par.BFL, ref.EFL, ref.BFL)
	}
	out, err := json.Marshal(&old)
	if err != nil {
		t.Fatal(err)
	}
	for _, key := range []string{"decenter", "tilt_"} {
		if strings.Contains(string(out), key) {
			t.Errorf("re-written legacy JSON contains %q: %s", key, out)
		}
	}
}

// TestCardPoseJSONRoundTrip checks that the decentre and tilt fields are read,
// converted (degrees to radians) into the geo surface pose, and survive a
// write/read cycle and Clone unchanged.
func TestCardPoseJSONRoundTrip(t *testing.T) {
	const js = `{"name":"decentred doublet","object_z_mm":"inf","stop":0,"cards":[
	  {"radius_mm":80,"thickness_mm":5,"glass":"N-BK7","semi_aperture_mm":12,
	   "decenter_x_mm":0.5,"decenter_y_mm":-0.25,"tilt_x_deg":2.5,"tilt_y_deg":-1,"tilt_z_deg":30},
	  {"radius_mm":-40,"thickness_mm":2,"glass":"F2","semi_aperture_mm":12},
	  {"radius_mm":-200,"thickness_mm":0,"semi_aperture_mm":12}]}`
	var p Prescription
	if err := json.Unmarshal([]byte(js), &p); err != nil {
		t.Fatal(err)
	}
	c := p.Cards[0]
	if c.DecenterXMM != 0.5 || c.DecenterYMM != -0.25 || c.TiltXDeg != 2.5 || c.TiltYDeg != -1 || c.TiltZDeg != 30 {
		t.Fatalf("fields not read: %+v", c)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	want := geo.CoordBreak{Dx: 0.5, Dy: -0.25, TiltX: 2.5 * math.Pi / 180, TiltY: -math.Pi / 180, TiltZ: math.Pi / 6}
	if sys.Surfaces[0].Pose == nil {
		t.Fatal("decentred card has no pose")
	}
	got := *sys.Surfaces[0].Pose
	for name, pair := range map[string][2]float64{
		"Dx": {got.Dx, want.Dx}, "Dy": {got.Dy, want.Dy}, "TiltX": {got.TiltX, want.TiltX},
		"TiltY": {got.TiltY, want.TiltY}, "TiltZ": {got.TiltZ, want.TiltZ},
	} {
		if math.Abs(pair[0]-pair[1]) > 1e-15 {
			t.Errorf("pose %s = %.17g, want %.17g", name, pair[0], pair[1])
		}
	}
	if p := sys.Surfaces[1].Pose; p != nil && !p.IsIdentity() {
		t.Errorf("the pose of card 0 leaked into card 1: %+v", *p)
	}

	buf, err := json.Marshal(&p)
	if err != nil {
		t.Fatal(err)
	}
	var back Prescription
	if err := json.Unmarshal(buf, &back); err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(p.Cards, back.Cards) {
		t.Errorf("cards changed in the write/read cycle:\n%+v\n%+v", p.Cards, back.Cards)
	}
	if !reflect.DeepEqual(p.Cards, p.Clone().Cards) {
		t.Errorf("Clone lost the pose fields")
	}
	if !strings.Contains(p.String(), "decentre") {
		t.Errorf("String() does not report the decentred card: %s", p.String())
	}
}

// TestCardPoseRejectsNonFinite checks that a non-finite decentre or tilt is
// rejected by Validate (JSON cannot carry NaN, so the struct is set directly).
func TestCardPoseRejectsNonFinite(t *testing.T) {
	for _, mutate := range []func(*Card){
		func(c *Card) { c.DecenterXMM = math.NaN() },
		func(c *Card) { c.TiltYDeg = math.Inf(1) },
	} {
		p := biconvexBK7()
		mutate(&p.Cards[1])
		if err := p.Validate(); err == nil {
			t.Errorf("non-finite pose accepted: %+v", p.Cards[1])
		}
	}
}

// TestCardPoseFirstOrderRejectedTraceWorks checks the behaviour of a decentred
// prescription: the first-order and merit functions (Paraxial, SpotRMS) return
// an error, while the geo system traces real rays.
func TestCardPoseFirstOrderRejectedTraceWorks(t *testing.T) {
	p := biconvexBK7()
	p.Cards[0].DecenterXMM = 0.3
	if _, err := p.Paraxial(geo.Wavelengthd); err == nil {
		t.Error("Paraxial of a decentred prescription: want an error")
	}
	if _, _, _, err := SpotRMS(p, geo.Wavelengthd, Field{}, UnitDiskSamples(2, 6)); err == nil {
		t.Error("SpotRMS of a decentred prescription: want an error (paraxial image plane)")
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	res := sys.Trace(geo.NewRay(geo.V(0.1, 0, -20), geo.V(0, 0, 1), geo.Wavelengthd))
	if !res.OK() {
		t.Errorf("real-ray trace of the decentred prescription failed: %v", res.Status)
	}
}

// TestCardPoseTiltedMirrorFromPrescription repeats the axial mirror check of
// package geo through a prescription: a plane mirror card tilted by 4° about x
// deflects an axial ray by 8° (1e-12 rad).
func TestCardPoseTiltedMirrorFromPrescription(t *testing.T) {
	const js = `{"name":"tilted fold mirror","object_z_mm":"inf","stop":0,"cards":[
	  {"radius_mm":0,"thickness_mm":0,"mirror":true,"tilt_x_deg":4}]}`
	var p Prescription
	if err := json.Unmarshal([]byte(js), &p); err != nil {
		t.Fatal(err)
	}
	sys, err := p.System()
	if err != nil {
		t.Fatal(err)
	}
	in := geo.NewRay(geo.V(0, 0, -10), geo.V(0, 0, 1), geo.Wavelengthd)
	res := sys.Trace(in)
	if !res.OK() {
		t.Fatalf("trace failed: %v", res.Status)
	}
	out := res.Ray.Dir
	back := in.Dir.Neg()
	angle := math.Atan2(out.Cross(back).Norm(), out.Dot(back))
	want := 8 * math.Pi / 180
	if math.Abs(angle-want) > 1e-12 {
		t.Errorf("deflection %.15g rad, want 2θ = %.15g rad", angle, want)
	}
}
