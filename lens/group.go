package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// ---------------------------------------------------------------------------
// Element-level lens groups: the editable model
// ---------------------------------------------------------------------------
//
// A Prescription describes a system surface by surface (one Card per surface).
// A Group describes it element by element: a lens, a mirror or an aperture
// stop with its own position, thickness, glass, surfaces (radius and equation
// parameters), clear aperture and 3D pose. It is the editable model a GUI or a
// CLI shows, and it converts both ways:
//
//	Group.System()        -> geo.System       (real-ray tracing, always)
//	Group.Prescription()  -> *Prescription    (when the layout is expressible
//	                                           as a card chain)
//	Group.Layout()        -> global-frame vertex position and axes of every
//	                         surface, folds included
//	Prescription.Group()  -> *Group           (the card chain as elements)
//
// 3D. Every element carries a pose (decentre along x and y, tilts about x, y
// and z, in the sense of geo.CoordBreak). A mirror that tilts the beam about x
// or y folds the sequential chain (FoldAuto): the following elements are laid
// out along the reflected beam, which is how a folded system is built. The
// chain is expressed in the geo.System with the geo coordinate-break elements
// (Surface.Break) that package geo documents for chained frames: entering the
// element frame, then either the return break (geo.CoordBreak.Inverse) or the
// fold break that redirects the axis.

// ElementKind classifies an element of a Group. It is a string so that the
// JSON form is readable and a wrong kind is caught by Validate.
type ElementKind string

// Element kinds.
const (
	// ElementLens is a refracting element. With a rear surface it is the
	// classical lens: a front surface, the element's glass, a rear surface and
	// a centre thickness. Without a rear surface it is a single refracting
	// surface (a cemented interface split off a doublet, a field surface, a
	// window face). Both surfaces may be planes, so an ElementLens with two
	// plane surfaces is a plane-parallel plate.
	ElementLens ElementKind = "lens"
	// ElementMirror is a mirror: one reflective surface with a radius and an
	// equation of the same choice as a lens surface. It has no thickness (the
	// substrate is not traversed by the sequential beam); the gap carries the
	// distance to the next element and may be negative (the axial return of a
	// Cassegrain). A mirror tilted about x or y folds the chain, see FoldMode.
	ElementMirror ElementKind = "mirror"
	// ElementStop is an aperture stop: a plane surface with a clear aperture
	// that does not change the beam (the classical "aperture stop" card). The
	// group's Stop index normally points at it, but any element can be the
	// stop.
	ElementStop ElementKind = "stop"
)

// FoldMode controls the redirection of the sequential chain after an element.
type FoldMode int

const (
	// FoldAuto (the zero value) redirects the chain after a mirror that tilts
	// the beam (TiltX or TiltY nonzero): the following elements are laid out
	// along the reflected beam, so a 45 degree mirror really folds the system
	// into 3D. A mirror without a tilt about x or y leaves the chain axis
	// alone; its axial return is written as a negative gap (a Cassegrain).
	FoldAuto FoldMode = iota
	// FoldOff never redirects the chain, whatever the element does. Use it for
	// a tilted mirror that steers the beam out of the sequence (the following
	// elements then stay on the incoming axis) or for an array of tilted
	// elements whose chain is written by hand.
	FoldOff
)

// String labels the fold mode.
func (m FoldMode) String() string {
	if m == FoldOff {
		return "off"
	}
	return "auto"
}

// Face selects one surface of an element.
type Face int

// Faces.
const (
	// Front is the first surface of an element, the one the beam meets first.
	Front Face = iota
	// Rear is the second surface (lenses only).
	Rear
)

// String labels the face.
func (f Face) String() string {
	if f == Rear {
		return "rear"
	}
	return "front"
}

// SurfaceSpec is one surface of an element: the shape vocabulary of Card (the
// vertex radius of curvature, the conic constant, and one of the even asphere
// coefficients, the general polynomial terms, the Forbes Q coefficients, the
// freeform XY polynomial or the Noll-index Zernike list) plus the clear
// semi-aperture. It is a Card without the card's thickness, glass, mirror flag
// and pose, which belong to the element.
//
// RadiusMM == 0 or ±Inf means a plane. SemiAperMM <= 0 means unlimited.
type SurfaceSpec struct {
	RadiusMM   float64        `json:"radius_mm"`
	ConicK     float64        `json:"conic_k,omitempty"`
	Coeffs     []float64      `json:"coeffs,omitempty"`        // A4, A6, A8, ... of geo.EvenAsphere
	PolyTerms  []PolyTermCard `json:"poly_terms,omitempty"`    // geo.PolyAsphere terms
	QAsphere   *QCard         `json:"q_asphere,omitempty"`     // geo.QAsphere (Forbes)
	XYPoly     *XYCard        `json:"xy_polynomial,omitempty"` // geo.XYPolynomial (freeform)
	Zernike    *ZernikeCard   `json:"zernike,omitempty"`       // geo.ZernikeSurface
	SemiAperMM float64        `json:"semi_aperture_mm,omitempty"`
}

// card returns the Card image of the surface, so that the profile builder and
// the asphere validation of the package (which work on cards) apply unchanged.
func (s SurfaceSpec) card() Card {
	return Card{
		RadiusMM: s.RadiusMM, ConicK: s.ConicK, Coeffs: s.Coeffs,
		PolyTerms: s.PolyTerms, QAsphere: s.QAsphere, XYPoly: s.XYPoly, Zernike: s.Zernike,
		SemiAperMM: s.SemiAperMM,
	}
}

// ofCard returns the surface spec of a card: the shape fields only.
func surfaceOfCard(c Card) SurfaceSpec {
	return SurfaceSpec{
		RadiusMM: c.RadiusMM, ConicK: c.ConicK, Coeffs: c.Coeffs,
		PolyTerms: c.PolyTerms, QAsphere: c.QAsphere, XYPoly: c.XYPoly, Zernike: c.Zernike,
		SemiAperMM: c.SemiAperMM,
	}
}

// asphereKinds counts the aspheric descriptions set on the surface (see
// Card.asphereKinds).
func (s SurfaceSpec) asphereKinds() int { return s.card().asphereKinds() }

// profile converts the surface into a geo surface profile.
func (s SurfaceSpec) profile() (geo.Profile, error) { return s.card().profile() }

// isPlane reports whether the surface is flat (radius 0 or ±Inf and no
// polynomial term that would make it curved).
func (s SurfaceSpec) isPlane() bool { return s.RadiusMM == 0 || math.IsInf(s.RadiusMM, 0) }

// validate checks one surface; what names it in the error messages.
func (s SurfaceSpec) validate(what string) error {
	if math.IsNaN(s.RadiusMM) {
		return fmt.Errorf("%s: radius is NaN", what)
	}
	if s.SemiAperMM < 0 {
		return fmt.Errorf("%s: negative semi-aperture %g", what, s.SemiAperMM)
	}
	for j, a := range s.Coeffs {
		if !lensFinite(a) {
			return fmt.Errorf("%s: coefficient A%d is not finite (%v)", what, 2*j+4, a)
		}
	}
	if s.asphereKinds() > 1 {
		return fmt.Errorf("%s: coeffs, poly_terms, q_asphere, xy_polynomial and zernike are mutually exclusive", what)
	}
	if len(s.PolyTerms) > 0 || s.QAsphere != nil || s.XYPoly != nil || s.Zernike != nil {
		if _, err := s.profile(); err != nil {
			return fmt.Errorf("%s: %w", what, err)
		}
	}
	return nil
}

// Element is one editable element of a lens group. Its fields are the ones a
// GUI edits: kind, surfaces (radius and equation), material, thickness,
// position in the chain and 3D pose.
//
// Chain geometry. Elements are laid out in beam order. The first surface of
// element i sits at the distance Element[i-1].GapMM after the last surface of
// element i-1, measured along the chain axis; the chain axis is the global +z
// until a mirror folds it (FoldAuto) and then follows the reflected beam. The
// thickness of a lens is measured along the element's own axis (the element is
// rigid), so a tilted lens keeps its centre thickness.
//
// Media. The medium in front of an element is the GapGlass of the preceding
// element (air for the first). Inside a two-surface element the medium is
// Glass; after the element it is GapGlass again. That is the medium chain of
// the Card model, grouped by element, and it is what makes a cemented doublet
// expressible: lens [R1 R2] with Glass "N-BK7" and GapGlass "F2" at gap 0
// followed by a single-surface lens [R3] is a BK7/F2 doublet.
type Element struct {
	Name string      `json:"name,omitempty"`
	Kind ElementKind `json:"kind"`

	// Front is the first surface; Rear is the second one and is absent for
	// mirrors and stops.
	Front SurfaceSpec  `json:"front"`
	Rear  *SurfaceSpec `json:"rear,omitempty"`

	// Glass is the medium between Front and Rear (a two-surface element); ""
	// or "air" is air, otherwise a name accepted by geo.Glass.
	Glass string `json:"glass,omitempty"`
	// GapGlass is the medium filling the gap between this element and the next
	// one ("" = air).
	GapGlass string `json:"gap_glass,omitempty"`

	// ThicknessMM is the axial distance from Front to Rear along the element's
	// own axis: the centre thickness of a lens. It must be 0 for a
	// single-surface element (a mirror has no thickness in the beam path).
	ThicknessMM float64 `json:"thickness_mm"`
	// GapMM is the axial distance from the element's last surface to the first
	// surface of the next element, measured along the chain axis after the
	// element (which a folding mirror has redirected). It may be negative: an
	// untilted mirror returns the beam along -z and the next element is placed
	// at a negative gap (a Cassegrain). For the last element the gap is the
	// nominal distance to the image plane; it is carried to the Prescription
	// and is not used by the analysis routines, which search for the best
	// focus themselves.
	GapMM float64 `json:"gap_mm"`

	// DecenterXMM/DecenterYMM decentre the element in the chain frame (mm) and
	// TiltXDeg/TiltYDeg/TiltZDeg tilt it about x, y and z (degrees), decentre
	// first, then the tilts: the convention of geo.CoordBreak. The element is
	// rigid (its surfaces move together) and the following elements are not
	// affected.
	DecenterXMM float64 `json:"decenter_x_mm,omitempty"`
	DecenterYMM float64 `json:"decenter_y_mm,omitempty"`
	TiltXDeg    float64 `json:"tilt_x_deg,omitempty"`
	TiltYDeg    float64 `json:"tilt_y_deg,omitempty"`
	TiltZDeg    float64 `json:"tilt_z_deg,omitempty"`

	// Fold selects the chain redirection of a tilted mirror, see FoldMode.
	Fold FoldMode `json:"fold,omitempty"`
}

// LensElement returns a two-surface lens element.
func LensElement(radiusFront, radiusRear, thicknessMM float64, glass string, semiAperMM float64) Element {
	return Element{
		Kind:  ElementLens,
		Front: SurfaceSpec{RadiusMM: radiusFront, SemiAperMM: semiAperMM},
		Rear:  &SurfaceSpec{RadiusMM: radiusRear, SemiAperMM: semiAperMM},
		Glass: glass, ThicknessMM: thicknessMM,
	}
}

// SingleSurface returns a one-surface refracting element (a cemented interface
// or a bare refracting surface).
func SingleSurface(radiusMM float64, semiAperMM float64) Element {
	return Element{Kind: ElementLens, Front: SurfaceSpec{RadiusMM: radiusMM, SemiAperMM: semiAperMM}}
}

// MirrorElement returns a mirror element with the given vertex radius, conic
// constant and clear semi-aperture.
func MirrorElement(radiusMM, conicK, semiAperMM float64) Element {
	return Element{Kind: ElementMirror, Front: SurfaceSpec{RadiusMM: radiusMM, ConicK: conicK, SemiAperMM: semiAperMM}}
}

// StopElement returns an aperture stop of the given clear semi-aperture.
func StopElement(semiAperMM float64) Element {
	return Element{Kind: ElementStop, Front: SurfaceSpec{SemiAperMM: semiAperMM}}
}

// posed reports whether the element is decentred or tilted.
func (e Element) posed() bool {
	return e.DecenterXMM != 0 || e.DecenterYMM != 0 || e.TiltXDeg != 0 || e.TiltYDeg != 0 || e.TiltZDeg != 0
}

// pose returns the coordinate break of the element's pose, nil for an
// undecentred, untilted element.
func (e Element) pose() *geo.CoordBreak {
	if !e.posed() {
		return nil
	}
	const d2r = math.Pi / 180
	return &geo.CoordBreak{
		Dx:    e.DecenterXMM,
		Dy:    e.DecenterYMM,
		TiltX: e.TiltXDeg * d2r,
		TiltY: e.TiltYDeg * d2r,
		TiltZ: e.TiltZDeg * d2r,
	}
}

// folds reports whether the chain is redirected after the element.
func (e Element) folds() bool {
	if e.Fold != FoldAuto {
		return false
	}
	return e.Kind == ElementMirror && (e.TiltXDeg != 0 || e.TiltYDeg != 0)
}

// surfaces returns the number of surfaces of the element.
func (e Element) surfaces() int {
	if e.Rear != nil {
		return 2
	}
	return 1
}

// label names a face of the element for reports and errors.
func (e Element) label(i int, f Face) string {
	if e.Name != "" {
		return fmt.Sprintf("%s %s", e.Name, f)
	}
	return fmt.Sprintf("element %d %s", i, f)
}

// validate checks one element; i is its index in the group.
func (e Element) validate(i int) error {
	where := fmt.Sprintf("lens: element %d (%s)", i, e.Kind)
	if e.Name != "" {
		where = fmt.Sprintf("lens: element %d %q (%s)", i, e.Name, e.Kind)
	}
	switch e.Kind {
	case ElementLens, ElementMirror, ElementStop:
	default:
		return fmt.Errorf("%s: unknown kind %q (known: %q, %q, %q)", where, string(e.Kind), ElementLens, ElementMirror, ElementStop)
	}
	if e.Fold != FoldAuto && e.Fold != FoldOff {
		return fmt.Errorf("%s: unknown fold mode %d", where, int(e.Fold))
	}
	for _, v := range [...]float64{e.DecenterXMM, e.DecenterYMM, e.TiltXDeg, e.TiltYDeg, e.TiltZDeg} {
		if !lensFinite(v) {
			return fmt.Errorf("%s: non-finite decentre or tilt (%v)", where, v)
		}
	}
	if !lensFinite(e.GapMM) {
		return fmt.Errorf("%s: gap is not finite (%v)", where, e.GapMM)
	}
	if !lensFinite(e.ThicknessMM) || e.ThicknessMM < 0 {
		return fmt.Errorf("%s: thickness %v must be finite and >= 0", where, e.ThicknessMM)
	}
	if _, err := resolveGlass(e.GapGlass); err != nil {
		return fmt.Errorf("%s: gap glass: %w", where, err)
	}
	if err := e.Front.validate(where + " front"); err != nil {
		return err
	}
	if e.Rear == nil {
		if e.Kind == ElementLens && !isAirName(e.Glass) {
			return fmt.Errorf("%s: a one-surface lens has no glass of its own (glass %q); set the medium of the gap with the preceding element", where, e.Glass)
		}
		if e.ThicknessMM != 0 {
			return fmt.Errorf("%s: a one-surface element has thickness 0, got %g", where, e.ThicknessMM)
		}
	} else {
		if e.Kind != ElementLens {
			return fmt.Errorf("%s: only a lens has a rear surface", where)
		}
		if err := e.Rear.validate(where + " rear"); err != nil {
			return err
		}
		if _, err := resolveGlass(e.Glass); err != nil {
			return fmt.Errorf("%s: glass: %w", where, err)
		}
		if e.Glass == "" || isAirName(e.Glass) {
			if e.ThicknessMM == 0 {
				return fmt.Errorf("%s: a two-surface element with an air gap must have a positive thickness", where)
			}
		}
	}
	if e.Kind == ElementStop {
		if e.Front.asphereKinds() > 0 {
			return fmt.Errorf("%s: a stop is a plane surface, it cannot carry an aspheric description", where)
		}
		if !e.Front.isPlane() {
			return fmt.Errorf("%s: a stop is a plane surface, got radius %g", where, e.Front.RadiusMM)
		}
		if !isAirName(e.Glass) || !isAirName(e.GapGlass) {
			return fmt.Errorf("%s: a stop is in air (glass %q, gap glass %q)", where, e.Glass, e.GapGlass)
		}
	}
	if e.Kind == ElementMirror && !isAirName(e.Glass) {
		return fmt.Errorf("%s: a mirror has no glass of its own (glass %q)", where, e.Glass)
	}
	return nil
}

// Group is a lens group: an ordered list of editable elements in 3D, with the
// object position, the aperture stop and the global origin. It converts to a
// geo.System for tracing and real-ray analysis, to a Prescription when the
// layout is expressible as a card chain, and to JSON.
type Group struct {
	Name string `json:"name,omitempty"`

	// ObjectZMM is the axial position (mm) of the object point; math.Inf(-1)
	// is an object at infinity. The real-ray measurements take their object
	// from their Bundle; this field is the design object of the group and is
	// carried to the Prescription.
	ObjectZMM float64 `json:"object_z_mm"`

	Elements []Element `json:"elements"`

	// Stop is the index of the aperture-stop element: the stop surface is the
	// first surface of that element. The zero value (0) is the common
	// front-stop default.
	Stop int `json:"stop"`

	// OriginZ is the axial position (mm) of the first element's first surface
	// in the global frame (geo.System.OriginZ). The zero value puts the front
	// vertex at z = 0.
	OriginZ float64 `json:"origin_z_mm,omitempty"`
}

// NewGroup returns an empty group with an object at infinity.
func NewGroup(name string) *Group {
	return &Group{Name: name, ObjectZMM: math.Inf(-1)}
}

// Len returns the number of elements.
func (g *Group) Len() int { return len(g.Elements) }

// Clone returns a deep copy (the elements and every coefficient slice are
// duplicated).
func (g *Group) Clone() *Group {
	out := *g
	out.Elements = make([]Element, len(g.Elements))
	for i := range g.Elements {
		out.Elements[i] = cloneElement(g.Elements[i])
	}
	return &out
}

// cloneElement deep-copies one element.
func cloneElement(e Element) Element {
	out := e
	if e.Front.Coeffs != nil {
		out.Front.Coeffs = append([]float64(nil), e.Front.Coeffs...)
	}
	if e.Front.PolyTerms != nil {
		out.Front.PolyTerms = append([]PolyTermCard(nil), e.Front.PolyTerms...)
	}
	if e.Front.QAsphere != nil {
		q := *e.Front.QAsphere
		q.Coeffs = append([]float64(nil), q.Coeffs...)
		out.Front.QAsphere = &q
	}
	if e.Front.XYPoly != nil {
		x := *e.Front.XYPoly
		x.Terms = append([]XYTermCard(nil), x.Terms...)
		out.Front.XYPoly = &x
	}
	if e.Front.Zernike != nil {
		z := *e.Front.Zernike
		z.Coeffs = append([]float64(nil), z.Coeffs...)
		out.Front.Zernike = &z
	}
	if e.Rear != nil {
		r := cloneElement(Element{Kind: e.Kind, Front: *e.Rear})
		out.Rear = &r.Front
	}
	return out
}

// Validate checks the whole group and returns a descriptive error.
func (g *Group) Validate() error {
	if len(g.Elements) == 0 {
		return fmt.Errorf("lens: group %q has no elements", g.Name)
	}
	if !lensFinite(g.ObjectZMM) && !math.IsInf(g.ObjectZMM, -1) {
		return fmt.Errorf("lens: group %q: object position %v must be finite or -Inf", g.Name, g.ObjectZMM)
	}
	if !lensFinite(g.OriginZ) {
		return fmt.Errorf("lens: group %q: origin z %v is not finite", g.Name, g.OriginZ)
	}
	if g.Stop < 0 || g.Stop >= len(g.Elements) {
		return fmt.Errorf("lens: group %q: stop index %d out of range [0,%d]", g.Name, g.Stop, len(g.Elements)-1)
	}
	for i := range g.Elements {
		if err := g.Elements[i].validate(i); err != nil {
			return err
		}
	}
	return nil
}

// glassBefore returns the name of the medium in front of element i (the gap
// glass of the preceding element, air for the first).
func (g *Group) glassBefore(i int) string {
	if i <= 0 || i > len(g.Elements) {
		return ""
	}
	return g.Elements[i-1].GapGlass
}

// lensFinite reports whether x is a finite number.
func lensFinite(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }

// ---------------------------------------------------------------------------
// Editing: the element list
// ---------------------------------------------------------------------------

// element returns a pointer to element i after a range check.
func (g *Group) element(i int) (*Element, error) {
	if i < 0 || i >= len(g.Elements) {
		return nil, fmt.Errorf("lens: element %d out of range [0,%d]", i, len(g.Elements)-1)
	}
	return &g.Elements[i], nil
}

// ElementAt returns a copy of element i.
func (g *Group) ElementAt(i int) (Element, error) {
	e, err := g.element(i)
	if err != nil {
		return Element{}, err
	}
	return cloneElement(*e), nil
}

// SetElement replaces element i, checking it first. Replacing an element keeps
// the order of the list; the pose, kind and surfaces may all change.
func (g *Group) SetElement(i int, e Element) error {
	if _, err := g.element(i); err != nil {
		return err
	}
	if err := e.validate(i); err != nil {
		return err
	}
	g.Elements[i] = cloneElement(e)
	return nil
}

// AddElement appends an element and returns its index.
func (g *Group) AddElement(e Element) (int, error) {
	if err := e.validate(len(g.Elements)); err != nil {
		return 0, err
	}
	g.Elements = append(g.Elements, cloneElement(e))
	return len(g.Elements) - 1, nil
}

// InsertElement inserts e before element i (i == Len() appends).
func (g *Group) InsertElement(i int, e Element) error {
	if i < 0 || i > len(g.Elements) {
		return fmt.Errorf("lens: insert position %d out of range [0,%d]", i, len(g.Elements))
	}
	if err := e.validate(i); err != nil {
		return err
	}
	g.Elements = append(g.Elements, Element{})
	copy(g.Elements[i+1:], g.Elements[i:])
	g.Elements[i] = cloneElement(e)
	// The stop index follows the insertion.
	if g.Stop >= i {
		g.Stop++
	}
	return nil
}

// RemoveElement removes element i. The stop index follows the removal; a group
// may become empty, which Validate reports.
func (g *Group) RemoveElement(i int) error {
	if _, err := g.element(i); err != nil {
		return err
	}
	g.Elements = append(g.Elements[:i], g.Elements[i+1:]...)
	if g.Stop > i {
		g.Stop--
	}
	if g.Stop >= len(g.Elements) && len(g.Elements) > 0 {
		g.Stop = len(g.Elements) - 1
	}
	if len(g.Elements) == 0 {
		g.Stop = 0
	}
	return nil
}

// MoveElement moves element from to position to, shifting the elements in
// between (to == Len()-1 moves it to the end). The relative order of the other
// elements is unchanged.
func (g *Group) MoveElement(from, to int) error {
	if _, err := g.element(from); err != nil {
		return err
	}
	if to < 0 || to >= len(g.Elements) {
		return fmt.Errorf("lens: move target %d out of range [0,%d]", to, len(g.Elements)-1)
	}
	if from == to {
		return nil
	}
	stop := g.Stop
	e := g.Elements[from]
	g.Elements = append(g.Elements[:from], g.Elements[from+1:]...)
	rest := append([]Element{}, g.Elements...)
	g.Elements = append(g.Elements[:to], append([]Element{e}, rest[to:]...)...)
	// Remap the stop index through the permutation.
	indexOf := func(k int) int {
		switch {
		case k == from:
			return to
		case from < to && k > from && k <= to:
			return k - 1
		case to < from && k >= to && k < from:
			return k + 1
		}
		return k
	}
	g.Stop = indexOf(stop)
	return nil
}

// SetKind changes the kind of element i. Switching to a mirror or a stop drops
// the rear surface and the thickness; switching back to a lens leaves a single
// surface. The element is validated before it is written.
func (g *Group) SetKind(i int, kind ElementKind) error {
	e, err := g.element(i)
	if err != nil {
		return err
	}
	cand := cloneElement(*e)
	cand.Kind = kind
	switch kind {
	case ElementMirror, ElementStop:
		cand.Rear = nil
		cand.ThicknessMM = 0
		cand.Glass = ""
	case ElementLens:
		if cand.Front.isPlane() && cand.Front.asphereKinds() == 0 {
			// A stop becomes a plane plate, which needs a positive thickness.
			cand.ThicknessMM = math.Max(cand.ThicknessMM, 1)
		}
	}
	if err := cand.validate(i); err != nil {
		return err
	}
	*e = cand
	return nil
}

// ---------------------------------------------------------------------------
// Editing: single parameters
// ---------------------------------------------------------------------------

// GroupParamKind selects the group quantity a GroupParam addresses.
type GroupParamKind int

// Parameter kinds. The surface kinds address the front or rear surface of an
// element (GroupParam.Surf); the element kinds use the element index only.
const (
	GPRadius       GroupParamKind = iota // vertex radius R (mm); 0 or ±Inf is a plane
	GPCurvature                          // vertex curvature c = 1/R (1/mm); 0 is a plane
	GPConic                              // conic constant K
	GPCoeff                              // even asphere coefficient Coeffs[Order] (0 = A4)
	GPPolyTerm                           // polynomial term PolyTerms[Order].Coeff
	GPQCoef                              // Forbes Q coefficient QAsphere.Coeffs[Order]
	GPXYTerm                             // freeform term XYPoly.Terms[Order].Coeff
	GPZernike                            // Zernike coefficient Zernike.Coeffs[Order] (Noll j = Order+1)
	GPSemiAperture                       // clear semi-aperture (mm); <= 0 is unlimited
	GPThickness                          // axial distance from the front to the rear surface (mm)
	GPGap                                // axial distance from this element to the next (mm)
	GPPosition                           // axial chain position of the front surface (mm)
	GPDecenterX                          // decentre along x (mm)
	GPDecenterY                          // decentre along y (mm)
	GPTiltX                              // tilt about x (degrees)
	GPTiltY                              // tilt about y (degrees)
	GPTiltZ                              // roll about z (degrees)
)

// SurfRef addresses one surface of a group: the front or the rear face of an
// element. The element-level parameter kinds ignore it.
type SurfRef struct {
	Element int
	Face    Face
}

// String labels the reference, e.g. "e0.front".
func (r SurfRef) String() string { return fmt.Sprintf("e%d.%s", r.Element, r.Face) }

// GroupParam addresses one editable quantity of a group, the editing analogue
// of the optimization VarSpec and the sensitivity Param.
type GroupParam struct {
	Surf  SurfRef
	Kind  GroupParamKind
	Order int // coefficient index for the coefficient kinds
}

// String labels the parameter, e.g. "e0.front.R", "e1.thk", "e2.gap".
func (p GroupParam) String() string {
	switch p.Kind {
	case GPRadius:
		return p.Surf.String() + ".R"
	case GPCurvature:
		return p.Surf.String() + ".c"
	case GPConic:
		return p.Surf.String() + ".K"
	case GPCoeff:
		return fmt.Sprintf("%s.A%d", p.Surf, 2*p.Order+4)
	case GPPolyTerm:
		return fmt.Sprintf("%s.poly%d", p.Surf, p.Order)
	case GPQCoef:
		return fmt.Sprintf("%s.q%d", p.Surf, p.Order)
	case GPXYTerm:
		return fmt.Sprintf("%s.xy%d", p.Surf, p.Order)
	case GPZernike:
		return fmt.Sprintf("%s.z%d", p.Surf, p.Order+1)
	case GPSemiAperture:
		return p.Surf.String() + ".semi"
	case GPThickness:
		return fmt.Sprintf("e%d.thk", p.Surf.Element)
	case GPGap:
		return fmt.Sprintf("e%d.gap", p.Surf.Element)
	case GPPosition:
		return fmt.Sprintf("e%d.z", p.Surf.Element)
	case GPDecenterX:
		return fmt.Sprintf("e%d.dx", p.Surf.Element)
	case GPDecenterY:
		return fmt.Sprintf("e%d.dy", p.Surf.Element)
	case GPTiltX:
		return fmt.Sprintf("e%d.tx", p.Surf.Element)
	case GPTiltY:
		return fmt.Sprintf("e%d.ty", p.Surf.Element)
	case GPTiltZ:
		return fmt.Sprintf("e%d.tz", p.Surf.Element)
	}
	return fmt.Sprintf("e%d.?", p.Surf.Element)
}

// surface returns a pointer to the addressed surface spec.
func (g *Group) surface(r SurfRef) (*SurfaceSpec, error) {
	e, err := g.element(r.Element)
	if err != nil {
		return nil, err
	}
	if r.Face == Rear {
		if e.Rear == nil {
			return nil, fmt.Errorf("lens: element %d (%s) has no rear surface", r.Element, e.Kind)
		}
		return e.Rear, nil
	}
	return &e.Front, nil
}

// GetParam reads one quantity of the group.
func (g *Group) GetParam(p GroupParam) (float64, error) {
	e, err := g.element(p.Surf.Element)
	if err != nil {
		return 0, err
	}
	switch p.Kind {
	case GPThickness:
		return e.ThicknessMM, nil
	case GPGap:
		return e.GapMM, nil
	case GPPosition:
		return g.PositionMM(p.Surf.Element)
	case GPDecenterX:
		return e.DecenterXMM, nil
	case GPDecenterY:
		return e.DecenterYMM, nil
	case GPTiltX:
		return e.TiltXDeg, nil
	case GPTiltY:
		return e.TiltYDeg, nil
	case GPTiltZ:
		return e.TiltZDeg, nil
	}
	s, err := g.surface(p.Surf)
	if err != nil {
		return 0, fmt.Errorf("lens: parameter %s: %w", p, err)
	}
	switch p.Kind {
	case GPRadius:
		return s.RadiusMM, nil
	case GPCurvature:
		if s.isPlane() {
			return 0, nil
		}
		return 1 / s.RadiusMM, nil
	case GPConic:
		return s.ConicK, nil
	case GPSemiAperture:
		return s.SemiAperMM, nil
	case GPCoeff:
		if p.Order < 0 || p.Order >= len(s.Coeffs) {
			return 0, fmt.Errorf("lens: parameter %s: coefficient out of range (surface has %d)", p, len(s.Coeffs))
		}
		return s.Coeffs[p.Order], nil
	case GPPolyTerm:
		if p.Order < 0 || p.Order >= len(s.PolyTerms) {
			return 0, fmt.Errorf("lens: parameter %s: polynomial term out of range (surface has %d)", p, len(s.PolyTerms))
		}
		return s.PolyTerms[p.Order].Coeff, nil
	case GPQCoef:
		if s.QAsphere == nil {
			return 0, fmt.Errorf("lens: parameter %s: the surface has no q_asphere", p)
		}
		if p.Order < 0 || p.Order >= len(s.QAsphere.Coeffs) {
			return 0, fmt.Errorf("lens: parameter %s: q coefficient out of range", p)
		}
		return s.QAsphere.Coeffs[p.Order], nil
	case GPXYTerm:
		if s.XYPoly == nil {
			return 0, fmt.Errorf("lens: parameter %s: the surface has no xy_polynomial", p)
		}
		if p.Order < 0 || p.Order >= len(s.XYPoly.Terms) {
			return 0, fmt.Errorf("lens: parameter %s: freeform term out of range", p)
		}
		return s.XYPoly.Terms[p.Order].Coeff, nil
	case GPZernike:
		if s.Zernike == nil {
			return 0, fmt.Errorf("lens: parameter %s: the surface has no zernike list", p)
		}
		if p.Order < 0 || p.Order >= len(s.Zernike.Coeffs) {
			return 0, fmt.Errorf("lens: parameter %s: Zernike coefficient out of range", p)
		}
		return s.Zernike.Coeffs[p.Order], nil
	}
	return 0, fmt.Errorf("lens: parameter %s: unknown kind %d", p, int(p.Kind))
}

// SetParam writes one quantity of the group after checking it. Setting the
// vertex radius of a surface makes a plane when the value is 0 or ±Inf;
// setting the curvature 0 makes a plane, too (R = ±Inf).
func (g *Group) SetParam(p GroupParam, v float64) error {
	if !lensFinite(v) && !(p.Kind == GPRadius && math.IsInf(v, 0)) {
		return fmt.Errorf("lens: parameter %s: value %v is not finite", p, v)
	}
	e, err := g.element(p.Surf.Element)
	if err != nil {
		return err
	}
	switch p.Kind {
	case GPThickness:
		if v < 0 {
			return fmt.Errorf("lens: parameter %s: thickness %g is negative", p, v)
		}
		if e.Rear == nil && v != 0 {
			return fmt.Errorf("lens: parameter %s: element %d (%s) has one surface, its thickness is 0", p, p.Surf.Element, e.Kind)
		}
		e.ThicknessMM = v
		return nil
	case GPGap:
		e.GapMM = v
		return nil
	case GPPosition:
		return g.SetPositionMM(p.Surf.Element, v)
	case GPDecenterX:
		e.DecenterXMM = v
		return nil
	case GPDecenterY:
		e.DecenterYMM = v
		return nil
	case GPTiltX:
		e.TiltXDeg = v
		return nil
	case GPTiltY:
		e.TiltYDeg = v
		return nil
	case GPTiltZ:
		e.TiltZDeg = v
		return nil
	}
	s, err := g.surface(p.Surf)
	if err != nil {
		return fmt.Errorf("lens: parameter %s: %w", p, err)
	}
	switch p.Kind {
	case GPRadius:
		if v == 0 {
			s.RadiusMM = math.Inf(1)
			return nil
		}
		s.RadiusMM = v
		return nil
	case GPCurvature:
		if v == 0 {
			s.RadiusMM = math.Inf(1)
			return nil
		}
		s.RadiusMM = 1 / v
		return nil
	case GPConic:
		s.ConicK = v
		return nil
	case GPSemiAperture:
		if v < 0 {
			return fmt.Errorf("lens: parameter %s: semi-aperture %g is negative", p, v)
		}
		s.SemiAperMM = v
		return nil
	case GPCoeff:
		if p.Order < 0 || p.Order >= len(s.Coeffs) {
			return fmt.Errorf("lens: parameter %s: coefficient out of range (surface has %d)", p, len(s.Coeffs))
		}
		if s.asphereKinds() > 1 {
			return fmt.Errorf("lens: parameter %s: the surface already carries another aspheric description", p)
		}
		s.Coeffs[p.Order] = v
		return nil
	case GPPolyTerm:
		if p.Order < 0 || p.Order >= len(s.PolyTerms) {
			return fmt.Errorf("lens: parameter %s: polynomial term out of range (surface has %d)", p, len(s.PolyTerms))
		}
		if s.asphereKinds() > 1 {
			return fmt.Errorf("lens: parameter %s: the surface already carries another aspheric description", p)
		}
		s.PolyTerms[p.Order].Coeff = v
		return nil
	case GPQCoef:
		if s.QAsphere == nil {
			return fmt.Errorf("lens: parameter %s: the surface has no q_asphere", p)
		}
		if p.Order < 0 || p.Order >= len(s.QAsphere.Coeffs) {
			return fmt.Errorf("lens: parameter %s: q coefficient out of range", p)
		}
		s.QAsphere.Coeffs[p.Order] = v
		return nil
	case GPXYTerm:
		if s.XYPoly == nil {
			return fmt.Errorf("lens: parameter %s: the surface has no xy_polynomial", p)
		}
		if p.Order < 0 || p.Order >= len(s.XYPoly.Terms) {
			return fmt.Errorf("lens: parameter %s: freeform term out of range", p)
		}
		s.XYPoly.Terms[p.Order].Coeff = v
		return nil
	case GPZernike:
		if s.Zernike == nil {
			return fmt.Errorf("lens: parameter %s: the surface has no zernike list", p)
		}
		if p.Order < 0 || p.Order >= len(s.Zernike.Coeffs) {
			return fmt.Errorf("lens: parameter %s: Zernike coefficient out of range", p)
		}
		s.Zernike.Coeffs[p.Order] = v
		return nil
	}
	return fmt.Errorf("lens: parameter %s: unknown kind %d", p, int(p.Kind))
}

// ---------------------------------------------------------------------------
// Editing: named quantities
// ---------------------------------------------------------------------------

// GapMM returns the axial gap to the next element.
func (g *Group) GapMM(i int) (float64, error) {
	e, err := g.element(i)
	if err != nil {
		return 0, err
	}
	return e.GapMM, nil
}

// SetGapMM sets the axial gap to the next element (mm, may be negative).
func (g *Group) SetGapMM(i int, v float64) error {
	return g.SetParam(GroupParam{Surf: SurfRef{Element: i}, Kind: GPGap}, v)
}

// chainAdvance returns the amount by which the chain advances inside element i
// before the gap of the element starts: the centre thickness of an unposed
// two-surface element. A posed element does not advance the chain (its return
// break restores the frame in which the following gap is measured), so a
// decentred or tilted element does not move the following elements.
func (g *Group) chainAdvance(i int) float64 {
	e := &g.Elements[i]
	if e.Rear == nil || e.posed() {
		return 0
	}
	return e.ThicknessMM
}

// PositionMM returns the axial chain position of the front surface of element
// i (mm): the origin z plus, for every preceding element, its gap and, for an
// unposed two-surface element, its thickness. For an unfolded group this is the
// global z of the element's first vertex (it is what Layout reports); for a
// folded chain it is the distance along the chain, which the fold has turned.
func (g *Group) PositionMM(i int) (float64, error) {
	if _, err := g.element(i); err != nil {
		return 0, err
	}
	z := g.OriginZ
	for k := 0; k < i; k++ {
		z += g.Elements[k].GapMM + g.chainAdvance(k)
	}
	return z, nil
}

// SetPositionMM moves element i so that its front surface sits at the axial
// chain position z (mm): the gap of the preceding element is adjusted, or the
// group origin for the first element. The positions of the other elements are
// unchanged.
func (g *Group) SetPositionMM(i int, z float64) error {
	if _, err := g.element(i); err != nil {
		return err
	}
	if i == 0 {
		g.OriginZ = z
		return nil
	}
	prev := g.OriginZ
	for k := 0; k < i-1; k++ {
		prev += g.Elements[k].GapMM + g.chainAdvance(k)
	}
	g.Elements[i-1].GapMM = z - prev - g.chainAdvance(i-1)
	return nil
}

// SetThicknessMM sets the centre thickness of element i (the axial distance
// from its front to its rear surface).
func (g *Group) SetThicknessMM(i int, v float64) error {
	return g.SetParam(GroupParam{Surf: SurfRef{Element: i}, Kind: GPThickness}, v)
}

// ThicknessMM returns the centre thickness of element i.
func (g *Group) ThicknessMM(i int) (float64, error) {
	return g.GetParam(GroupParam{Surf: SurfRef{Element: i}, Kind: GPThickness})
}

// SetGlass sets the medium between the two surfaces of element i (a glass name
// accepted by geo.Glass, or "" / "air").
func (g *Group) SetGlass(i int, name string) error {
	e, err := g.element(i)
	if err != nil {
		return err
	}
	if _, err := resolveGlass(name); err != nil {
		return err
	}
	if e.Rear == nil && !isAirName(name) {
		return fmt.Errorf("lens: element %d (%s) has one surface and no glass of its own", i, e.Kind)
	}
	if e.Kind == ElementMirror && !isAirName(name) {
		return fmt.Errorf("lens: element %d is a mirror and has no glass", i)
	}
	e.Glass = name
	return nil
}

// SetGapGlass sets the medium filling the gap after element i ("" / "air", or
// a glass name). A cemented interface uses it.
func (g *Group) SetGapGlass(i int, name string) error {
	e, err := g.element(i)
	if err != nil {
		return err
	}
	if _, err := resolveGlass(name); err != nil {
		return err
	}
	if e.Kind == ElementStop && !isAirName(name) {
		return fmt.Errorf("lens: element %d is a stop and must stay in air", i)
	}
	e.GapGlass = name
	return nil
}

// SetFace sets the vertex radius of one face of element i (mm; 0 or ±Inf is a
// plane).
func (g *Group) SetRadiusMM(i int, f Face, v float64) error {
	return g.SetParam(GroupParam{Surf: SurfRef{Element: i, Face: f}, Kind: GPRadius}, v)
}

// SetCurvature sets the vertex curvature c = 1/R of one face (1/mm; 0 is a
// plane, R = ±Inf).
func (g *Group) SetCurvature(i int, f Face, c float64) error {
	return g.SetParam(GroupParam{Surf: SurfRef{Element: i, Face: f}, Kind: GPCurvature}, c)
}

// SetSemiAperture sets the clear semi-aperture of one face of element i (mm;
// <= 0 is unlimited).
func (g *Group) SetSemiAperture(i int, f Face, v float64) error {
	return g.SetParam(GroupParam{Surf: SurfRef{Element: i, Face: f}, Kind: GPSemiAperture}, v)
}

// SetPose sets the whole 3D pose of element i: the decentre (dx, dy) in mm and
// the tilts about x, y and z in degrees.
func (g *Group) SetPose(i int, dx, dy, tiltXdeg, tiltYdeg, tiltZdeg float64) error {
	e, err := g.element(i)
	if err != nil {
		return err
	}
	for _, v := range [...]float64{dx, dy, tiltXdeg, tiltYdeg, tiltZdeg} {
		if !lensFinite(v) {
			return fmt.Errorf("lens: element %d pose: %v is not finite", i, v)
		}
	}
	e.DecenterXMM, e.DecenterYMM = dx, dy
	e.TiltXDeg, e.TiltYDeg, e.TiltZDeg = tiltXdeg, tiltYdeg, tiltZdeg
	return nil
}

// SetFold sets the fold mode of element i.
func (g *Group) SetFold(i int, m FoldMode) error {
	e, err := g.element(i)
	if err != nil {
		return err
	}
	if m != FoldAuto && m != FoldOff {
		return fmt.Errorf("lens: element %d: unknown fold mode %d", i, int(m))
	}
	e.Fold = m
	return nil
}

// String returns a one-line summary per element, useful in CLIs.
func (g *Group) String() string {
	s := fmt.Sprintf("group %q: %d elements, stop %d, origin z=%.4gmm, object %s",
		g.Name, len(g.Elements), g.Stop, g.OriginZ, describeZ(g.ObjectZMM))
	for i, e := range g.Elements {
		r := "plane"
		if !e.Front.isPlane() {
			r = fmt.Sprintf("R=%.4g", e.Front.RadiusMM)
		}
		kind := string(e.Kind)
		rear := ""
		if e.Rear != nil {
			rr := "plane"
			if !e.Rear.isPlane() {
				rr = fmt.Sprintf("R=%.4g", e.Rear.RadiusMM)
			}
			rear = fmt.Sprintf(" + %s", rr)
		}
		pose := ""
		if e.posed() {
			pose = fmt.Sprintf(" [decentre (%g, %g) mm, tilt (%g, %g, %g) deg]", e.DecenterXMM, e.DecenterYMM, e.TiltXDeg, e.TiltYDeg, e.TiltZDeg)
		}
		if e.folds() {
			pose += " [fold]"
		}
		gl := e.Glass
		if isAirName(gl) {
			gl = "air"
		}
		s += fmt.Sprintf("\n  [%d] %s %s%s%s → %.4gmm into %s, gap %.4gmm", i, kind, r, rear, pose, e.ThicknessMM, gl, e.GapMM)
	}
	return s
}

// describeZ formats an object position.
func describeZ(z float64) string {
	if math.IsInf(z, -1) {
		return "at infinity"
	}
	return fmt.Sprintf("at z=%.4gmm", z)
}
