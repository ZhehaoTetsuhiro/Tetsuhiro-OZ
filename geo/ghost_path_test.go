package geo

import (
	"math"
	"reflect"
	"testing"
)

// ---------------------------------------------------------------------------
// Systems shared by the ghost tests
// ---------------------------------------------------------------------------

// ghostDoublet builds the cemented N-BK7/N-SF5 achromat of the ghost tests: a
// 200 mm thin-lens achromat solution of the classical pair of conditions
//
//	φ₁ + φ₂ = 1/200 mm⁻¹ ,   φ₁/V₁ + φ₂/V₂ = 0 ,   V(BK7) = 64.17, V(SF5) = 32.25
//
// (φ₁ = 0.010052, φ₂ = −0.005052 mm⁻¹) with the radii rounded to
// R₁ = (nd₁−1)/φ₁ = 51.41 mm and R₃ = −(nd₂−1)/φ₂ = 133.15 mm. The crown is
// plane-convex towards the front and the flint plane-concave behind it, so the
// cement is a plane: the two inner surfaces, 1 (cement) and 2 (rear), carry the
// classical ghost pair I = 2, J = 1.
func ghostDoublet(t *testing.T) *System {
	t.Helper()
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 51.41}, Aperture: CircleAperture(12), Name: "crown front"},
			{Profile: Plane{}, Aperture: CircleAperture(12), Name: "cement"},
			{Profile: Sphere{R: 133.15}, Aperture: CircleAperture(12), Name: "flint rear"},
		},
		[]float64{4, 3},
		[]Medium{Air{}, N_BK7, N_SF5, Air{}})
}

// ghostSinglet builds a symmetric biconvex N-BK7 singlet (R = ±61.7 mm,
// t = 4 mm, about 60 mm focal length) for the single-lens ghost: its only pair
// is I = 1, J = 0, the reflection off the rear surface followed by the
// reflection off the front surface.
func ghostSinglet(t *testing.T) *System {
	t.Helper()
	return mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: 61.7}, Aperture: CircleAperture(12), Name: "front"},
			{Profile: Sphere{R: -61.7}, Aperture: CircleAperture(12), Name: "rear"},
		},
		[]float64{4},
		[]Medium{Air{}, N_BK7, Air{}})
}

// ---------------------------------------------------------------------------
// Path enumeration
// ---------------------------------------------------------------------------

// TestGhostPairsDeterministic freezes the enumeration order of the
// double-reflection pairs (I ascending, then J ascending) and checks that a
// repeated call returns the identical sequence.
func TestGhostPairsDeterministic(t *testing.T) {
	cases := []struct {
		n    int
		want []GhostPair
	}{
		{0, nil},
		{1, nil},
		{2, []GhostPair{{I: 1, J: 0}}},
		{3, []GhostPair{{I: 1, J: 0}, {I: 2, J: 0}, {I: 2, J: 1}}},
		{4, []GhostPair{{I: 1, J: 0}, {I: 2, J: 0}, {I: 2, J: 1}, {I: 3, J: 0}, {I: 3, J: 1}, {I: 3, J: 2}}},
	}
	for _, c := range cases {
		got := GhostPairs(c.n)
		if len(got) != len(c.want) {
			t.Fatalf("GhostPairs(%d): got %d pairs, want %d", c.n, len(got), len(c.want))
		}
		for i := range got {
			if got[i] != c.want[i] {
				t.Errorf("GhostPairs(%d)[%d] = %s, want %s", c.n, i, got[i], c.want[i])
			}
		}
		if again := GhostPairs(c.n); !reflect.DeepEqual(got, again) {
			t.Errorf("GhostPairs(%d) is not deterministic: %v vs %v", c.n, got, again)
		}
	}
}

// TestGhostNodesInnerSurfaceDoublet freezes the node sequence of the classical
// two-inner-surface ghost of the cemented doublet: refract through the crown
// front, refract through the cement into the flint, reflect at the rear
// surface, travel back to the cement and reflect there, then leave through the
// rear surface again. The gap indices (Medium) are the side of the surface the
// ray arrives from, which is what the reversed pass swaps.
func TestGhostNodesInnerSurfaceDoublet(t *testing.T) {
	sys := ghostDoublet(t)
	path, err := sys.GhostPathOf(GhostPair{I: 2, J: 1})
	if err != nil {
		t.Fatal(err)
	}
	want := []GhostNode{
		{Surface: 0, Action: GhostRefract, Medium: 0},
		{Surface: 1, Action: GhostRefract, Medium: 1},
		{Surface: 2, Action: GhostReflect, Medium: 2},
		{Surface: 1, Action: GhostReflect, Backward: true, Medium: 2},
		{Surface: 2, Action: GhostRefract, Medium: 2},
	}
	if len(path.Nodes) != len(want) {
		t.Fatalf("ghost %s: got %d nodes, want %d (%v)", path.Pair, len(path.Nodes), len(want), path.Nodes)
	}
	for i := range want {
		if path.Nodes[i] != want[i] {
			t.Errorf("node %d = %s, want %s", i, path.Nodes[i], want[i])
		}
	}
	// The first reflection must be the one at the later surface, the second
	// reflection must travel backwards, and the second pass must cross the
	// first reflection site again.
	if path.Nodes[2].Surface != 2 || path.Nodes[2].Backward {
		t.Errorf("first reflection must be the forward hit of surface 2, got %s", path.Nodes[2])
	}
	if !path.Nodes[3].Backward || path.Nodes[3].Action != GhostReflect {
		t.Errorf("second reflection must be backwards at the cement, got %s", path.Nodes[3])
	}
	if path.Nodes[4].Surface != 2 || path.Nodes[4].Action != GhostRefract {
		t.Errorf("the path must cross the first reflection site again, got %s", path.Nodes[4])
	}
}

// TestGhostNodesMediumBookkeeping walks the nodes of every pair of a
// four-surface system and checks the medium bookkeeping of the reversed pass:
// the ray starts in gap 0, a refraction moves it to the gap on the other side
// of the surface (m−1 for a backwards pass, m+1 for a forwards one), and a
// reflection leaves it where it was. The number of backwards nodes must be
// exactly the surfaces strictly between the two reflection sites.
func TestGhostNodesMediumBookkeeping(t *testing.T) {
	const n = 4
	for _, p := range GhostPairs(n) {
		nodes, err := GhostNodes(n, p)
		if err != nil {
			t.Fatalf("GhostNodes(%d, %s): %v", n, p, err)
		}
		if want := 2*p.I + n - 2*p.J; len(nodes) != want {
			t.Errorf("%s: got %d nodes, want %d", p, len(nodes), want)
		}
		gap := 0
		reflections, backwards := 0, 0
		for i, nd := range nodes {
			if nd.Medium != gap {
				t.Fatalf("%s: node %d (%s) is in gap %d, want %d", p, i, nd, nd.Medium, gap)
			}
			if nd.Backward {
				backwards++
			}
			switch nd.Action {
			case GhostReflect:
				reflections++
			case GhostRefract:
				if nd.Backward {
					gap = nd.Surface
				} else {
					gap = nd.Surface + 1
				}
			default:
				t.Fatalf("%s: node %d has action %v", p, i, nd.Action)
			}
		}
		if reflections != 2 {
			t.Errorf("%s: got %d reflections, want 2", p, reflections)
		}
		if want := p.I - p.J; backwards != want {
			t.Errorf("%s: got %d backwards nodes, want %d (the reverse refractions plus the second reflection at J)",
				p, backwards, want)
		}
		// The last node must leave the ray travelling forwards behind the
		// second reflection site, and the very last surface must be n−1.
		if last := nodes[len(nodes)-1]; last.Backward || last.Surface != n-1 || last.Action != GhostRefract {
			t.Errorf("%s: last node = %s, want a forward refraction at surface %d", p, last, n-1)
		}
	}
}

// TestGhostNodesRangeErrors checks the rejection of impossible pairs and of
// systems whose elements the double-reflection model cannot describe.
func TestGhostNodesRangeErrors(t *testing.T) {
	bad := []GhostPair{{I: 0, J: 0}, {I: 1, J: 1}, {I: 0, J: -1}, {I: 3, J: 1}, {I: -1, J: 0}}
	for _, p := range bad {
		if _, err := GhostNodes(3, p); err == nil {
			t.Errorf("GhostNodes(3, %s): expected an error", p)
		}
	}
	// A mirror surface cannot take part in a Fresnel ghost path.
	mirror := mustSystem(t, Air{},
		[]Surface{
			{Profile: Sphere{R: -200}, Aperture: CircleAperture(20)},
			{Profile: Sphere{R: -100}, Aperture: CircleAperture(10), Mirror: true},
		},
		[]float64{-100}, []Medium{Air{}, Air{}, Air{}})
	if _, err := mirror.GhostPathOf(GhostPair{I: 1, J: 0}); err == nil {
		t.Error("GhostPathOf on a mirror system: expected an error")
	}
	// A coordinate-break element has no profile and no interaction.
	brk := mustSystem(t, Air{},
		[]Surface{
			{Profile: Plane{}, Aperture: CircleAperture(10)},
			{Break: true, Name: "fold"},
			{Profile: Plane{}, Aperture: CircleAperture(10)},
		},
		[]float64{5, 5}, []Medium{Air{}, Air{}, Air{}, Air{}})
	if _, err := brk.GhostPathOf(GhostPair{I: 2, J: 0}); err == nil {
		t.Error("GhostPathOf with a coordinate break: expected an error")
	}
}

// ---------------------------------------------------------------------------
// Relative energy
// ---------------------------------------------------------------------------

// TestGhostReflectanceAnalytic freezes the Fresnel formulae and the
// per-surface override: R = ((n₁−n₂)/(n₁+n₂))² is symmetric in the two media
// (so a path contributes the same R whichever pass hits the surface), the
// amplitude r has the sign of n₁ − n₂, and an override outside [0,1] or a
// reflectance range error is reported instead of being used.
func TestGhostReflectanceAnalytic(t *testing.T) {
	checkClose(t, "r(air, 1.5)", FresnelAmplitude(1, 1.5), (1-1.5)/(1+1.5), 0)
	checkClose(t, "R(air, 1.5)", FresnelReflectance(1, 1.5), 0.04, 1e-15)
	checkClose(t, "R(1.5, air)", FresnelReflectance(1.5, 1), 0.04, 1e-15)
	checkClose(t, "R(1.67271, 1.51680)",
		FresnelReflectance(1.67271, 1.51680),
		math.Pow((1.67271-1.51680)/(1.67271+1.51680), 2), 1e-18)
	if r := FresnelAmplitude(1, 1.5); !(r < 0) {
		t.Errorf("r(air, glass) = %g, want a negative amplitude", r)
	}
	if !math.IsNaN(FresnelReflectance(-1, 1)) {
		t.Error("FresnelReflectance(-1, 1): want NaN for a degenerate denominator")
	}

	sys := ghostDoublet(t)
	lam := Wavelengthd
	ov := []float64{math.NaN(), 0.5, 0.25}
	r, err := sys.GhostReflectance(1, lam, ov)
	if err != nil || r != 0.5 {
		t.Errorf("override R(1) = %g (err %v), want 0.5", r, err)
	}
	r, err = sys.GhostReflectance(0, lam, ov)
	if err != nil {
		t.Fatalf("NaN override must fall back to Fresnel: %v", err)
	}
	checkClose(t, "Fresnel fallback R(0)", r,
		FresnelReflectance(sys.Mediums[0].N(lam), sys.Mediums[1].N(lam)), 0)
	if _, err := sys.GhostReflectance(1, lam, []float64{0, 1.5}); err == nil {
		t.Error("GhostReflectance: an override outside [0,1] must be rejected")
	}
	if _, err := sys.GhostReflectance(7, lam, nil); err == nil {
		t.Error("GhostReflectance: an out-of-range surface index must be rejected")
	}
}
