package geo

import (
	"fmt"
	"math"
	"reflect"
	"testing"
)

// ---------------------------------------------------------------------------
// Acceptance (a): the ghost spot of a small aperture against the paraxial chain
// ---------------------------------------------------------------------------

// ghostPrimaryHit returns the image-plane hit of the reference (chief) ray of
// the sampling through the nominal single-pass system: the primary image point
// the ghost spot is displaced from.
func ghostPrimaryHit(t *testing.T, sys *System, opt GhostOptions) Vec3 {
	t.Helper()
	ref := LaunchRay(opt.Field, opt.Pupil, PupilPoint{}, opt.LambdaNm)
	h := sys.traceToPlane(ref, opt.ImageZ)
	if !h.Status.OK() {
		t.Fatalf("the nominal reference ray did not reach the image plane: %v", h.Status)
	}
	return h.Hit
}

// TestGhostDoubletSpotVsParaxial is acceptance (a) for the cemented doublet:
// for every double-reflection path the centroid of a small-aperture sampling
// near the axis must reproduce the paraxial landing point of the reversed ABCD
// chain to better than 1 %. The relative difference is measured against the
// ghost displacement |centroid − primary| — the distance between the ghost
// spot and the nominal image point of the same field, the quantity a ghost
// analysis is about (the reference ray itself sits on the axis of the sampling,
// where a relative difference of positions is ill-conditioned).
//
// The classic two-inner-surface pair of the doublet (I = 2 rear surface,
// J = 1 cement) is asserted to be present and analysed.
func TestGhostDoubletSpotVsParaxial(t *testing.T) {
	sys := ghostDoublet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	opt := GhostOptions{
		Field:    FieldTan(0.5*math.Pi/180, 0),
		Pupil:    Pupil{Z: -5, Radius: 0.25},
		ImageZ:   p.BackFocalZ,
		LambdaNm: Wavelengthd,
	}
	primary := ghostPrimaryHit(t, sys, opt)
	rep, err := sys.GhostAnalysis(opt)
	if err != nil {
		t.Fatal(err)
	}
	if rep.Enumerated != 3 || len(rep.Paths) != 3 {
		t.Fatalf("got %d of %d paths, want 3 of 3", len(rep.Paths), rep.Enumerated)
	}
	sawInner := false
	for _, r := range rep.Paths {
		if r.Skipped() {
			t.Fatalf("%s: %s", r.Pair, r.Note)
		}
		if r.Pair == (GhostPair{I: 2, J: 1}) {
			sawInner = true
		}
		if r.Spot.N != r.Spot.NTotal {
			t.Fatalf("%s: %d of %d rays reached the image plane", r.Pair, r.Spot.N, r.Spot.NTotal)
		}
		delta := r.Spot.Centroid.Sub(r.Paraxial).Norm()
		displacement := r.Spot.Centroid.Sub(primary).Norm()
		if displacement < 1e-3 {
			t.Fatalf("%s: ghost displacement %.3g mm is too small for a relative criterion", r.Pair, displacement)
		}
		rel := delta / displacement
		if rel > 1e-2 {
			t.Errorf("%s: |centroid − paraxial| = %.6g mm is %.6g %% of the ghost displacement %.6g mm, want <= 1 %%",
				r.Pair, delta, 100*rel, displacement)
		}
		t.Logf("%s: centroid %s, paraxial %s, |Δ| = %.3g mm, ghost displacement %.3g mm, relative %.4g %%",
			r.Pair, r.Spot.Centroid, r.Paraxial, delta, displacement, 100*rel)
	}
	if !sawInner {
		t.Error("the classic two-inner-surface pair (I=2, J=1) is missing from the report")
	}
}

// TestGhostSingletSpotVsParaxial is acceptance (a) for a single lens: the only
// pair (I = 1, J = 0) is the reflection off the rear surface followed by the
// reflection off the front surface, the classical single-lens ghost.
func TestGhostSingletSpotVsParaxial(t *testing.T) {
	sys := ghostSinglet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	opt := GhostOptions{
		Field:    FieldTan(0.5*math.Pi/180, 0),
		Pupil:    Pupil{Z: -5, Radius: 0.2},
		ImageZ:   p.BackFocalZ,
		LambdaNm: Wavelengthd,
	}
	primary := ghostPrimaryHit(t, sys, opt)
	rep, err := sys.GhostAnalysis(opt)
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Paths) != 1 || rep.Paths[0].Pair != (GhostPair{I: 1, J: 0}) {
		t.Fatalf("got %v, want the single pair (I=1, J=0)", rep.Paths)
	}
	r := rep.Paths[0]
	delta := r.Spot.Centroid.Sub(r.Paraxial).Norm()
	displacement := r.Spot.Centroid.Sub(primary).Norm()
	if rel := delta / displacement; rel > 1e-2 {
		t.Errorf("%s: |centroid − paraxial| = %.6g mm is %.6g %% of the ghost displacement %.6g mm, want <= 1 %%",
			r.Pair, delta, 100*rel, displacement)
	}
	t.Logf("singlet %s: centroid %s, paraxial %s, |Δ| = %.3g mm, ghost displacement %.3g mm",
		r.Pair, r.Spot.Centroid, r.Paraxial, delta, displacement)
}

// ---------------------------------------------------------------------------
// Acceptance (b): the relative energy is R_I·R_J
// ---------------------------------------------------------------------------

// TestGhostEnergyProductAnalytic is acceptance (b): the reported relative
// energy of every path equals the product of the normal-incidence intensity
// reflectances of its two reflection sites, R = ((n₁−n₂)/(n₁+n₂))², computed
// here from the media of the system with the formula written out (not through
// the package function), to 1e-12. The per-surface override is checked too.
func TestGhostEnergyProductAnalytic(t *testing.T) {
	sys := ghostDoublet(t)
	lam := Wavelengthd
	rOf := func(i int) float64 {
		n1 := sys.Mediums[i].N(lam)
		n2 := sys.Mediums[i+1].N(lam)
		d := (n1 - n2) / (n1 + n2)
		return d * d
	}
	p, err := sys.Paraxial(lam, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	rep, err := sys.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.5*math.Pi/180, 0), Pupil: Pupil{Z: -5, Radius: 0.25},
		ImageZ: p.BackFocalZ, LambdaNm: lam,
	})
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range rep.Paths {
		wantI, wantJ := rOf(r.Pair.I), rOf(r.Pair.J)
		checkClose(t, "reflectance of the first site", r.ReflectI, wantI, 1e-15)
		checkClose(t, "reflectance of the second site", r.ReflectJ, wantJ, 1e-15)
		checkClose(t, "relative energy R_I·R_J", r.Energy, wantI*wantJ, 1e-12)
	}
	// The chain of the doublet: air→BK7, BK7→SF5, SF5→air at 587.562 nm.
	checkClose(t, "R(air, BK7)", rOf(0), math.Pow((1-N_BK7.N(lam))/(1+N_BK7.N(lam)), 2), 1e-18)
	checkClose(t, "R(SF5, air)", rOf(2), math.Pow((N_SF5.N(lam)-1)/(N_SF5.N(lam)+1), 2), 1e-18)

	// A supplied per-surface R replaces the Fresnel value exactly.
	rep, err = sys.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.5*math.Pi/180, 0), Pupil: Pupil{Z: -5, Radius: 0.25},
		ImageZ: p.BackFocalZ, LambdaNm: lam,
		Reflectance: []float64{math.NaN(), 0.5, 0.25},
	})
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range rep.Paths {
		if r.Pair == (GhostPair{I: 2, J: 1}) {
			checkClose(t, "overridden energy of the inner ghost", r.Energy, 0.125, 0)
			checkClose(t, "overridden R_I (surface 2)", r.ReflectI, 0.25, 0)
			checkClose(t, "overridden R_J (surface 1)", r.ReflectJ, 0.5, 0)
		}
	}
}

// ---------------------------------------------------------------------------
// Acceptance (c): energy threshold, top-N and deterministic ordering
// ---------------------------------------------------------------------------

// TestGhostThresholdAndTopN is acceptance (c): the energy threshold drops the
// paths below it (a path exactly at the threshold is kept) and MaxPaths keeps
// the N strongest paths of the sorted report.
func TestGhostThresholdAndTopN(t *testing.T) {
	sys := ghostDoublet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	base := GhostOptions{
		Field: FieldTan(0.5*math.Pi/180, 0), Pupil: Pupil{Z: -5, Radius: 0.25},
		ImageZ: p.BackFocalZ, LambdaNm: Wavelengthd,
	}
	// The Fresnel energies of the doublet, strongest first:
	// (2,0) = 2.671e-3 > (2,1) = 1.514e-4 > (1,0) = 1.007e-4.
	all, err := sys.GhostAnalysis(base)
	if err != nil {
		t.Fatal(err)
	}
	if len(all.Paths) != 3 {
		t.Fatalf("unfiltered report: got %d paths, want 3", len(all.Paths))
	}
	want := []GhostPair{{I: 2, J: 0}, {I: 2, J: 1}, {I: 1, J: 0}}
	for i, w := range want {
		if all.Paths[i].Pair != w {
			t.Errorf("report[%d] = %s, want %s (energy ordering)", i, all.Paths[i].Pair, w)
		}
	}
	for i := 1; i < len(all.Paths); i++ {
		if all.Paths[i-1].Energy < all.Paths[i].Energy {
			t.Errorf("report is not sorted by decreasing energy at %d", i)
		}
	}

	cases := []struct {
		name      string
		minRel    float64
		maxPaths  int
		wantPairs []GhostPair
	}{
		{"no filter", 0, 0, want},
		{"threshold below all", 1e-5, 0, want},
		{"threshold exactly on the weakest", all.Paths[2].Energy, 0, want},
		{"threshold drops the weakest", 1.1e-4, 0, []GhostPair{{I: 2, J: 0}, {I: 2, J: 1}}},
		{"threshold keeps the strongest only", 2e-3, 0, []GhostPair{{I: 2, J: 0}}},
		{"threshold above all", 0.3, 0, nil},
		{"top 2", 0, 2, []GhostPair{{I: 2, J: 0}, {I: 2, J: 1}}},
		{"top 1", 0, 1, []GhostPair{{I: 2, J: 0}}},
		{"top 0 keeps all", 0, 0, want},
		{"threshold and top 1", 1e-5, 1, []GhostPair{{I: 2, J: 0}}},
	}
	for _, c := range cases {
		opt := base
		opt.MinRelative, opt.MaxPaths = c.minRel, c.maxPaths
		rep, err := sys.GhostAnalysis(opt)
		if err != nil {
			t.Fatalf("%s: %v", c.name, err)
		}
		if rep.Enumerated != 3 {
			t.Errorf("%s: Enumerated = %d, want 3 (the enumeration is independent of the filters)", c.name, rep.Enumerated)
		}
		if len(rep.Paths) != len(c.wantPairs) {
			t.Errorf("%s: got %d paths, want %d (%v)", c.name, len(rep.Paths), len(c.wantPairs), rep.Paths)
			continue
		}
		for i, w := range c.wantPairs {
			if rep.Paths[i].Pair != w {
				t.Errorf("%s: path %d = %s, want %s", c.name, i, rep.Paths[i].Pair, w)
			}
		}
	}
}

// TestGhostTieOrderByIndex freezes the tie-break of the ordering: paths of
// equal relative energy are ordered by their surface index pair (I ascending,
// then J ascending). A constant per-surface reflectance makes every path of a
// system have exactly the same energy, which is the degenerate case of the
// ordering.
func TestGhostTieOrderByIndex(t *testing.T) {
	sys := ghostDoublet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	rep, err := sys.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.5*math.Pi/180, 0), Pupil: Pupil{Z: -5, Radius: 0.25},
		ImageZ: p.BackFocalZ, LambdaNm: Wavelengthd,
		Reflectance: []float64{0.5, 0.5, 0.5},
	})
	if err != nil {
		t.Fatal(err)
	}
	want := GhostPairs(3) // (1,0), (2,0), (2,1): index order
	if len(rep.Paths) != len(want) {
		t.Fatalf("got %d paths, want %d", len(rep.Paths), len(want))
	}
	for i := range want {
		if rep.Paths[i].Energy != 0.25 {
			t.Fatalf("path %d: energy %g, want exactly 0.25", i, rep.Paths[i].Energy)
		}
		if rep.Paths[i].Pair != want[i] {
			t.Errorf("ties must be ordered by (I, J): path %d = %s, want %s", i, rep.Paths[i].Pair, want[i])
		}
	}
	// A four-surface system with the same constant reflectance exercises the
	// tie-break over all six pairs of the enumeration.
	four := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 80}, Aperture: CircleAperture(10)},
			{Profile: Sphere{R: -80}, Aperture: CircleAperture(10)},
			{Profile: Sphere{R: 60}, Aperture: CircleAperture(10)},
			{Profile: Sphere{R: -60}, Aperture: CircleAperture(10)},
		},
		[]float64{2, 3, 2},
		[]Medium{Air{}, N_BK7, N_SF5, Air{}, Air{}})
	rep, err = four.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.2*math.Pi/180, 0), Pupil: Pupil{Z: -5, Radius: 0.1},
		ImageZ: 100, LambdaNm: Wavelengthd,
		Reflectance: []float64{0.5, 0.5, 0.5, 0.5},
	})
	if err != nil {
		t.Fatal(err)
	}
	want = GhostPairs(4)
	if len(rep.Paths) != len(want) {
		t.Fatalf("four-surface system: got %d paths, want %d", len(rep.Paths), len(want))
	}
	for i := range want {
		if rep.Paths[i].Pair != want[i] {
			t.Errorf("four-surface ties: path %d = %s, want %s", i, rep.Paths[i].Pair, want[i])
		}
	}
}

// ---------------------------------------------------------------------------
// Acceptance (d): bit-for-bit determinism
// ---------------------------------------------------------------------------

// TestGhostDeterminismBits is acceptance (d): two runs of the same analysis
// produce bit-identical reports (every float compared by its bit pattern, not
// by a tolerance).
func TestGhostDeterminismBits(t *testing.T) {
	sys := ghostDoublet(t)
	p, err := sys.Paraxial(Wavelengthd, math.Inf(-1))
	if err != nil {
		t.Fatal(err)
	}
	run := func() *GhostReport {
		rep, err := sys.GhostAnalysis(GhostOptions{
			Field: FieldTan(0.5*math.Pi/180, 0), Pupil: Pupil{Z: -5, Radius: 0.25},
			ImageZ: p.BackFocalZ, LambdaNm: Wavelengthd, Points: Hexapolar(2),
		})
		if err != nil {
			t.Fatal(err)
		}
		return rep
	}
	a, b := run(), run()
	if diff := ghostBitsEqual(reflect.ValueOf(a).Elem(), reflect.ValueOf(b).Elem(), "report"); len(diff) != 0 {
		t.Errorf("two identical runs differ in %d places, first: %v", len(diff), diff[:min(len(diff), 4)])
	}
	if a.String() != b.String() {
		t.Error("the formatted reports differ")
	}
}

// ghostBitsEqual walks two values of the same type and returns a description of
// every difference, comparing float64 fields by their bit pattern (so ±0 and
// NaN are distinguished) and everything else by equality.
func ghostBitsEqual(a, b reflect.Value, path string) []string {
	if a.Type() != b.Type() {
		return []string{fmt.Sprintf("%s: type %v vs %v", path, a.Type(), b.Type())}
	}
	switch a.Kind() {
	case reflect.Struct:
		var bad []string
		for i := 0; i < a.NumField(); i++ {
			bad = append(bad, ghostBitsEqual(a.Field(i), b.Field(i), path+"."+a.Type().Field(i).Name)...)
		}
		return bad
	case reflect.Slice:
		if a.Len() != b.Len() {
			return []string{fmt.Sprintf("%s: length %d vs %d", path, a.Len(), b.Len())}
		}
		var bad []string
		for i := 0; i < a.Len(); i++ {
			bad = append(bad, ghostBitsEqual(a.Index(i), b.Index(i), fmt.Sprintf("%s[%d]", path, i))...)
		}
		return bad
	case reflect.Float64:
		if math.Float64bits(a.Float()) != math.Float64bits(b.Float()) {
			return []string{fmt.Sprintf("%s: %v (bits %x) vs %v (bits %x)",
				path, a.Float(), math.Float64bits(a.Float()), b.Float(), math.Float64bits(b.Float()))}
		}
		return nil
	case reflect.Ptr:
		if a.IsNil() != b.IsNil() {
			return []string{fmt.Sprintf("%s: nil mismatch", path)}
		}
		if a.IsNil() {
			return nil
		}
		return ghostBitsEqual(a.Elem(), b.Elem(), path)
	}
	if !reflect.DeepEqual(a.Interface(), b.Interface()) {
		return []string{fmt.Sprintf("%s: %v vs %v", path, a.Interface(), b.Interface())}
	}
	return nil
}

// ---------------------------------------------------------------------------
// Unsupported paths, poses and option validation
// ---------------------------------------------------------------------------

// TestGhostUnsupportedPaths checks that the paths the double-reflection model
// cannot describe are reported, ranked last and carry a note instead of
// silently producing numbers: mirrors (their nominal action is a reflection)
// and coordinate-break elements (no profile).
func TestGhostUnsupportedPaths(t *testing.T) {
	mirror := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: -200}, Aperture: CircleAperture(20)},
			{Profile: Sphere{R: -100}, Aperture: CircleAperture(10), Mirror: true},
		},
		[]float64{-100}, []Medium{Air{}, Air{}, Air{}})
	rep, err := mirror.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.001, 0), Pupil: Pupil{Z: -5, Radius: 0.1}, ImageZ: 0, LambdaNm: Wavelengthd,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Paths) != 1 || !rep.Paths[0].Skipped() {
		t.Fatalf("mirror system: got %v, want one skipped path", rep.Paths)
	}
	if rep.Paths[0].Energy != 0 {
		t.Errorf("a skipped path must carry energy 0, got %g", rep.Paths[0].Energy)
	}
	if !math.IsNaN(rep.Paths[0].Paraxial.X) {
		t.Errorf("a skipped path must not report a paraxial landing point, got %s", rep.Paths[0].Paraxial)
	}
	// The threshold that drops everything eliminates the skipped paths too.
	rep, err = mirror.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.001, 0), Pupil: Pupil{Z: -5, Radius: 0.1}, ImageZ: 0, LambdaNm: Wavelengthd,
		MinRelative: 1e-9,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Paths) != 0 {
		t.Errorf("with a positive threshold: got %d paths, want 0", len(rep.Paths))
	}

	brk := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(10)},
			{Break: true, Name: "fold"},
			{Profile: Plane{}, Aperture: CircleAperture(10)},
		},
		[]float64{5, 5}, []Medium{Air{}, Air{}, Air{}, Air{}})
	rep, err = brk.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.001, 0), Pupil: Pupil{Z: -5, Radius: 0.1}, ImageZ: 20, LambdaNm: Wavelengthd,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Paths) != 3 {
		t.Fatalf("break system: got %d paths, want all 3 skipped", len(rep.Paths))
	}
	for _, r := range rep.Paths {
		if !r.Skipped() {
			t.Errorf("%s: a path through a coordinate break must be skipped", r.Pair)
		}
	}
}

// TestGhostPosedSystemRealOnly checks the analysis of a decentred system: the
// real-ray replay follows the frames (so the ghost spot is well defined), while
// the first-order chain is defined for coaxial systems only and reports NaN
// instead of a wrong landing point.
func TestGhostPosedSystemRealOnly(t *testing.T) {
	pose := CoordBreak{Dx: 0.2}
	sys := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 61.7}, Aperture: CircleAperture(12), Pose: &pose},
			{Profile: Sphere{R: -61.7}, Aperture: CircleAperture(12)},
		},
		[]float64{4}, []Medium{Air{}, N_BK7, Air{}})
	rep, err := sys.GhostAnalysis(GhostOptions{
		Field: FieldTan(0.5*math.Pi/180, 0), Pupil: Pupil{Z: -5, Radius: 0.1},
		ImageZ: 60, LambdaNm: Wavelengthd,
	})
	if err != nil {
		t.Fatal(err)
	}
	if len(rep.Paths) != 1 {
		t.Fatalf("got %d paths, want 1", len(rep.Paths))
	}
	r := rep.Paths[0]
	if r.Skipped() {
		t.Fatalf("a posed system is traced with real rays, got: %s", r.Note)
	}
	if !math.IsNaN(r.Paraxial.X) || !math.IsNaN(r.Paraxial.Y) {
		t.Errorf("the paraxial landing of a decentred system must be NaN, got %s", r.Paraxial)
	}
	if r.Spot.N == 0 {
		t.Error("the real-ray replay of a decentred system must land rays")
	}
}

// TestGhostOptionValidation checks the argument validation of the analysis.
func TestGhostOptionValidation(t *testing.T) {
	sys := ghostDoublet(t)
	base := GhostOptions{
		Field: FieldTan(0, 0), Pupil: Pupil{Z: -5, Radius: 0.25},
		ImageZ: 190, LambdaNm: Wavelengthd,
	}
	bad := []struct {
		name string
		mod  func(o *GhostOptions)
	}{
		{"image plane NaN", func(o *GhostOptions) { o.ImageZ = math.NaN() }},
		{"wavelength zero", func(o *GhostOptions) { o.LambdaNm = 0 }},
		{"wavelength infinite", func(o *GhostOptions) { o.LambdaNm = math.Inf(1) }},
		{"pupil radius zero", func(o *GhostOptions) { o.Pupil.Radius = 0 }},
		{"pupil radius negative", func(o *GhostOptions) { o.Pupil.Radius = -1 }},
		{"pupil plane infinite", func(o *GhostOptions) { o.Pupil.Z = math.Inf(-1) }},
		{"pupil offset NaN", func(o *GhostOptions) { o.Pupil.Offset = V(math.NaN(), 0, 0) }},
		{"negative threshold", func(o *GhostOptions) { o.MinRelative = -1e-3 }},
		{"negative MaxPaths", func(o *GhostOptions) { o.MaxPaths = -1 }},
		{"non-finite sampling point", func(o *GhostOptions) { o.Points = []PupilPoint{{PX: math.Inf(1)}} }},
	}
	for _, c := range bad {
		opt := base
		c.mod(&opt)
		if _, err := sys.GhostAnalysis(opt); err == nil {
			t.Errorf("%s: expected an error", c.name)
		}
	}
	// A nil Points list falls back to Hexapolar(1): one reference ray plus its
	// six-point ring, all of them landing.
	rep, err := sys.GhostAnalysis(base)
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range rep.Paths {
		if r.Spot.NTotal != 7 || r.Spot.N != 7 {
			t.Errorf("%s: got %d of %d rays, want 7 of 7 (Hexapolar(1))", r.Pair, r.Spot.N, r.Spot.NTotal)
		}
	}
}
