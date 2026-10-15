package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// ---------------------------------------------------------------------------
// The sequential chain: frames, plan, System, Layout, Prescription
// ---------------------------------------------------------------------------
//
// A group is laid out as a chain of rigid frames, exactly the model package geo
// implements with its coordinate-break elements (geo.CoordBreak, see
// coordbreak.go): each element sits in the current frame, and the frame either
// returns to the enclosing axis after the element (the auto-return of a pose)
// or is redirected by the fold break of a tilted mirror. The chain is written
// into the geo.System as Break elements, so the trace, BoundsZ, VertexZ and the
// real-ray analyses of package geo walk the same geometry this file computes.

// gframe is a rigid frame of the chain: the origin o and the orthonormal local
// axes ex, ey, ez in the global frame, so that the local point (x, y, z) sits
// at o + x·ex + y·ey + z·ez.
type gframe struct {
	o          geo.Vec3
	ex, ey, ez geo.Vec3
}

// startFrame returns the frame of the group origin: the first surface sits at
// z = originZ with the global axes.
func startFrame(originZ float64) gframe {
	return gframe{
		o:  geo.Vec3{Z: originZ},
		ex: geo.Vec3{X: 1}, ey: geo.Vec3{Y: 1}, ez: geo.Vec3{Z: 1},
	}
}

// apply returns the frame obtained by the coordinate break cb, with the
// semantics of package geo (P_outer = d + Q·P_inner): the origin moves by the
// displacement measured in this frame (o + R·d) and the axes are rotated by Q
// in this frame (R·Q). The rotation is read off geo.CoordBreak.Apply, so the
// frames agree exactly with System.VertexZ and with the trace.
func (f gframe) apply(cb geo.CoordBreak) gframe {
	d := cb.Apply(geo.Vec3{})
	local := func(v geo.Vec3) geo.Vec3 { return cb.Apply(v).Sub(d) } // Q·v
	// The columns of Q: the images of the local unit vectors.
	colX, colY, colZ := local(geo.V(1, 0, 0)), local(geo.V(0, 1, 0)), local(geo.V(0, 0, 1))
	comb := func(c0, c1, c2 float64) geo.Vec3 {
		return f.ex.Mul(c0).Add(f.ey.Mul(c1)).Add(f.ez.Mul(c2))
	}
	return gframe{
		o:  f.o.Add(comb(cb.Dx, cb.Dy, cb.Dz)),
		ex: comb(colX.X, colX.Y, colX.Z),
		ey: comb(colY.X, colY.Y, colY.Z),
		ez: comb(colZ.X, colZ.Y, colZ.Z),
	}
}

// shift advances the frame along its own axis by d (mm).
func (f gframe) shift(d float64) gframe {
	return gframe{o: f.o.Add(f.ez.Mul(d)), ex: f.ex, ey: f.ey, ez: f.ez}
}

// dirToGlobal maps a local direction to the global frame.
func (f gframe) dirToGlobal(d geo.Vec3) geo.Vec3 {
	return f.ex.Mul(d.X).Add(f.ey.Mul(d.Y)).Add(f.ez.Mul(d.Z))
}

// breakWithAxis returns the coordinate break (zero displacement, zero roll)
// whose rotation maps the local z axis onto the unit direction w:
// Q = Rx(alpha)·Ry(beta) with beta = asin(w.X) and alpha = atan2(-w.Y, w.Z).
func breakWithAxis(w geo.Vec3) geo.CoordBreak {
	beta := math.Asin(math.Max(-1, math.Min(1, w.X)))
	alpha := math.Atan2(-w.Y, w.Z)
	return geo.CoordBreak{TiltX: alpha, TiltY: beta}
}

// foldBreak returns the break that redirects the chain of a mirror: the
// incoming chain axis axisIn (unit, global frame) is reflected at the mirror,
// whose vertex normal is the local +z of the frame f, and the frame after the
// break follows the reflected beam. The break is expressed in f (the frame of
// the mirror surface), so the fold happens at the mirror vertex.
func foldBreak(f gframe, axisIn geo.Vec3) geo.CoordBreak {
	// The incoming axis in the mirror frame.
	local := geo.Vec3{
		X: f.ex.Dot(axisIn),
		Y: f.ey.Dot(axisIn),
		Z: f.ez.Dot(axisIn),
	}
	w := geo.Reflect(local, geo.Vec3{Z: 1})
	return breakWithAxis(w)
}

// ---------------------------------------------------------------------------
// The plan
// ---------------------------------------------------------------------------

// planSurface is one entry of the chain: an optical surface or a pure
// coordinate break.
type planSurface struct {
	profile  geo.Profile
	mirror   bool
	semiAper float64
	name     string
	pose     *geo.CoordBreak // the break of this entry (breaks only)
	brk      bool
	element  int // owning element, -1 for none
	face     Face
}

// plan is the chain of a group: the entries in beam order with the medium
// behind each of them, the gaps, the stop and the global placement of the real
// surfaces.
type plan struct {
	surfaces []planSurface
	mediums  []geo.Medium // mediums[i] is the medium in front of surfaces[i]
	gaps     []float64    // gaps[i] is the distance from surfaces[i] to surfaces[i+1]
	stop     int

	// layout: one entry per real surface, in beam order
	real   []int // index into surfaces
	vertex []geo.Vec3
	axis   []geo.Vec3
	normal []geo.Vec3
	folds  []int // number of chain folds in front of the surface

	exit gframe // beam frame after the last element (its origin and axis)
}

// emit appends one entry with the medium that follows it and the gap to the
// next entry.
func (pl *plan) emit(ps planSurface, medAfter geo.Medium, gap float64) int {
	pl.surfaces = append(pl.surfaces, ps)
	pl.mediums = append(pl.mediums, medAfter)
	pl.gaps = append(pl.gaps, gap)
	return len(pl.surfaces) - 1
}

// record stores the global placement of one real surface.
func (pl *plan) record(k int, f gframe, prof geo.Profile, nfold int) {
	sf := geo.Surface{Profile: prof}
	pl.real = append(pl.real, k)
	pl.vertex = append(pl.vertex, f.o)
	pl.axis = append(pl.axis, f.ez)
	pl.normal = append(pl.normal, f.dirToGlobal(sf.Normal(0, 0)))
	pl.folds = append(pl.folds, nfold)
}

// plan walks the group and returns its chain.
func (g *Group) plan() (*plan, error) {
	if err := g.Validate(); err != nil {
		return nil, err
	}
	pl := &plan{mediums: []geo.Medium{geo.Air{}}}
	c := startFrame(g.OriginZ)
	nfold := 0
	for i := range g.Elements {
		e := &g.Elements[i]
		before, err := resolveGlass(g.glassBefore(i))
		if err != nil {
			return nil, fmt.Errorf("lens: element %d: gap glass: %w", i, err)
		}
		in := c
		ef := in
		var pose *geo.CoordBreak
		if e.posed() {
			pose = e.pose()
			ef = in.apply(*pose)
		}
		// Enter the element frame.
		if pose != nil {
			pl.emit(planSurface{brk: true, pose: pose, element: i, name: fmt.Sprintf("element %d frame", i)}, before, 0)
			c = ef
		}
		// Front surface.
		prof, err := e.Front.profile()
		if err != nil {
			return nil, fmt.Errorf("lens: element %d front: %w", i, err)
		}
		medAfter := before
		switch {
		case e.Rear != nil:
			medAfter, err = resolveGlass(e.Glass)
		default:
			medAfter, err = resolveGlass(e.GapGlass)
		}
		if err != nil {
			return nil, fmt.Errorf("lens: element %d: %w", i, err)
		}
		gapFront := e.GapMM
		switch {
		case e.Rear != nil:
			gapFront = e.ThicknessMM
		case pose != nil:
			// The closing break of the element carries the gap.
			gapFront = 0
		}
		k := pl.emit(planSurface{
			profile: prof, mirror: e.Kind == ElementMirror, semiAper: e.Front.SemiAperMM,
			name: e.label(i, Front), element: i, face: Front,
		}, medAfter, gapFront)
		pl.record(k, ef, prof, nfold)
		frameAfter := ef
		// Rear surface.
		if e.Rear != nil {
			rf := ef.shift(e.ThicknessMM)
			profR, err := e.Rear.profile()
			if err != nil {
				return nil, fmt.Errorf("lens: element %d rear: %w", i, err)
			}
			medRear, err := resolveGlass(e.GapGlass)
			if err != nil {
				return nil, fmt.Errorf("lens: element %d: gap glass: %w", i, err)
			}
			gapRear := e.GapMM
			if pose != nil {
				// The closing break of the element carries the gap.
				gapRear = 0
			}
			k2 := pl.emit(planSurface{
				profile: profR, semiAper: e.Rear.SemiAperMM,
				name: e.label(i, Rear), element: i, face: Rear,
			}, medRear, gapRear)
			pl.record(k2, rf, profR, nfold)
			frameAfter = rf
			medAfter = medRear
		}
		// Close the element: the return break of its pose, or the fold break of
		// a mirror that tilts the beam.
		switch {
		case pose != nil:
			closeCB := pose.Inverse()
			name := fmt.Sprintf("element %d return", i)
			if e.folds() {
				closeCB = foldBreak(frameAfter, in.ez)
				name = fmt.Sprintf("element %d fold", i)
				nfold++
			}
			after := frameAfter.apply(closeCB)
			pl.emit(planSurface{brk: true, pose: &closeCB, element: i, name: name}, medAfter, e.GapMM)
			pl.exit = after
			c = after.shift(e.GapMM)
		default:
			pl.exit = frameAfter
			c = frameAfter.shift(e.GapMM)
		}
	}
	// The stop: the front surface of the stop element.
	pl.stop = -1
	for _, k := range pl.real {
		if pl.surfaces[k].element == g.Stop && pl.surfaces[k].face == Front {
			pl.stop = k
		}
	}
	if pl.stop < 0 {
		return nil, fmt.Errorf("lens: group %q: no stop surface (stop element %d)", g.Name, g.Stop)
	}
	return pl, nil
}

// system builds the geo.System of the plan.
func (pl *plan) system(g *Group) (*geo.System, error) {
	n := len(pl.surfaces)
	surfs := make([]geo.Surface, n)
	for i, ps := range pl.surfaces {
		if ps.brk {
			surfs[i] = geo.Surface{Break: true, Pose: ps.pose, Name: ps.name}
			continue
		}
		sf := geo.Surface{Profile: ps.profile, Mirror: ps.mirror, Name: ps.name, Pose: ps.pose}
		if ps.semiAper > 0 {
			sf.Aperture = geo.CircleAperture(ps.semiAper)
		}
		surfs[i] = sf
	}
	thick := make([]float64, 0, max(n-1, 0))
	if n > 1 {
		thick = append(thick, pl.gaps[:n-1]...)
	}
	mediums := append([]geo.Medium(nil), pl.mediums...)
	sys := &geo.System{
		Surfaces: surfs, Thickness: thick, Mediums: mediums,
		StopIndex: pl.stop, OriginZ: g.OriginZ, Name: g.Name,
	}
	if err := sys.Validate(); err != nil {
		return nil, err
	}
	return sys, nil
}

// System converts the group into a geometry-ready geo.System: one entry per
// surface in beam order, plus the coordinate-break elements that express the
// element poses and the folds (Surface.Break, see the geo documentation of
// Geo.CoordBreak). The surface order is
//
//	element 0 front [rear] [frame/return/fold breaks] element 1 front [rear] ...
//
// and the index of a face is available from SurfaceIndex. Rays are traced with
// System.Trace; the first-order engine rejects a posed or folded system, the
// real-ray routines (ImagePosition, AngularMagnification, Dispersion,
// geo.SpotDiagram, geo.Trace) work for all of them.
func (g *Group) System() (*geo.System, error) {
	pl, err := g.plan()
	if err != nil {
		return nil, err
	}
	return pl.system(g)
}

// SurfaceIndex returns the index, in Group.System(), of one face of element i.
func (g *Group) SurfaceIndex(i int, f Face) (int, error) {
	pl, err := g.plan()
	if err != nil {
		return 0, err
	}
	for _, k := range pl.real {
		if pl.surfaces[k].element == i && pl.surfaces[k].face == f {
			return k, nil
		}
	}
	return 0, fmt.Errorf("lens: element %d has no %s surface", i, f)
}

// LayoutSurface is one optical surface of the group in the global frame, after
// every coordinate break of the chain (folds included).
type LayoutSurface struct {
	Element int      // element index
	Face    Face     // front or rear
	Surface int      // index of the surface in Group.System() (breaks included in the numbering)
	Vertex  geo.Vec3 // global vertex position (mm)
	Axis    geo.Vec3 // unit local +z axis in the global frame: the optical axis of the surface
	Normal  geo.Vec3 // unit surface normal at the vertex in the global frame
	Folds   int      // number of chain folds in front of the surface
}

// Layout returns the global placement of every surface of the group. The chain
// is walked exactly as System() does, so Surface indexes the surface in the
// system of the same group and the vertices agree with System.VertexZ. For a
// system without folds the axes are the global +z; a fold turns them.
func (g *Group) Layout() ([]LayoutSurface, error) {
	pl, err := g.plan()
	if err != nil {
		return nil, err
	}
	out := make([]LayoutSurface, 0, len(pl.real))
	for j, k := range pl.real {
		ps := pl.surfaces[k]
		out = append(out, LayoutSurface{
			Element: ps.element, Face: ps.face, Surface: k,
			Vertex: pl.vertex[j], Axis: pl.axis[j], Normal: pl.normal[j], Folds: pl.folds[j],
		})
	}
	return out, nil
}

// ExitFrame returns the beam frame after the last element: the vertex of the
// last surface and the unit direction of the outgoing beam (the folded axis
// for a folding mirror). It is the frame along which ImagePosition searches for
// the best focus.
func (g *Group) ExitFrame() (vertex, axis geo.Vec3, err error) {
	pl, err := g.plan()
	if err != nil {
		return geo.Vec3{}, geo.Vec3{}, err
	}
	return pl.exit.o, pl.exit.ez, nil
}

// ---------------------------------------------------------------------------
// Prescription: the card form of a group
// ---------------------------------------------------------------------------

// applyPose copies the pose of an element into a card.
func applyPose(c *Card, e *Element) {
	c.DecenterXMM, c.DecenterYMM = e.DecenterXMM, e.DecenterYMM
	c.TiltXDeg, c.TiltYDeg, c.TiltZDeg = e.TiltXDeg, e.TiltYDeg, e.TiltZDeg
}

// Prescription returns the card (surface-by-surface) form of the group. It
// succeeds only when the group is expressible as a card chain:
//
//   - no element may fold the chain except the last one (a card chain runs
//     along one axis; an out-of-plane 3D fold has no card form — use System and
//     Layout instead);
//   - a two-surface element may be decentred and rolled but not tilted about x
//     or y: the card model places each card on the chain axis (auto-return), so
//     the thickness of a tilted element would not follow the element axis.
//
// Everything else round-trips: an axial fold-back is written with negative
// thicknesses (a Cassegrain), a cemented interface with the GapGlass of the
// preceding element, mirrors with the mirror flag and the clear apertures with
// the semi-aperture. ObjectZMM is measured from the group origin, so a non-zero
// OriginZ shifts it. The last card carries the last element's gap (the nominal
// distance to the image plane).
func (g *Group) Prescription() (*Prescription, error) {
	if err := g.Validate(); err != nil {
		return nil, err
	}
	p := &Prescription{Name: g.Name, ObjectZMM: g.ObjectZMM - g.OriginZ, Stop: 0}
	for i := range g.Elements {
		e := &g.Elements[i]
		if e.Rear != nil && (e.TiltXDeg != 0 || e.TiltYDeg != 0) {
			return nil, fmt.Errorf("lens: group %q: element %d is a two-surface element tilted by (%g, %g) degrees: the card model places every card on the chain axis (auto-return), so the element thickness would not follow the element axis",
				g.Name, i, e.TiltXDeg, e.TiltYDeg)
		}
		if e.folds() && i != len(g.Elements)-1 {
			return nil, fmt.Errorf("lens: group %q: element %d is a mirror tilted about x or y by (%g, %g) degrees and folds the beam: the following elements lie on the folded axis, which a card chain cannot express (out-of-plane 3D fold); use Group.System() and Group.Layout()",
				g.Name, i, e.TiltXDeg, e.TiltYDeg)
		}
		if i == g.Stop {
			p.Stop = len(p.Cards)
		}
		front := e.Front.card()
		front.Name = e.Name
		front.Mirror = e.Kind == ElementMirror
		if e.Rear != nil {
			front.Glass = e.Glass
			front.ThicknessMM = e.ThicknessMM
		} else {
			front.Glass = e.GapGlass
			front.ThicknessMM = e.GapMM
		}
		applyPose(&front, e)
		p.Cards = append(p.Cards, front)
		if e.Rear != nil {
			rear := e.Rear.card()
			rear.Glass = e.GapGlass
			rear.ThicknessMM = e.GapMM
			applyPose(&rear, e)
			p.Cards = append(p.Cards, rear)
		}
	}
	if err := p.Validate(); err != nil {
		return nil, err
	}
	return p, nil
}

// Group returns the element-level model of the prescription. Consecutive cards
// that form a glass block are paired into two-surface elements: a cemented
// doublet [R1 R2 R3] with the glasses BK7/F2 becomes a lens [R1 R2] whose Glass
// is BK7 and whose GapGlass is F2 at gap 0, followed by the single-surface
// element R3 — the same surfaces and media as the cards, grouped. An odd card
// at the end of a block stays a single-surface element.
//
// An element carries the name of the first card that forms it (the names of the
// following cards of a pair are not kept: the element has one name).
//
// The conversion is faithful (the two systems have the same surfaces, media and
// thicknesses) except for the cases the element model cannot hold, which are
// reported as errors: a reflective surface inside a glass block (a Mangin
// mirror), and a tilted card pair with a non-zero thickness (the card model
// measures that thickness along the chain axis, the element model along the
// element axis).
func (p *Prescription) Group() (*Group, error) {
	if err := p.Validate(); err != nil {
		return nil, err
	}
	g := &Group{Name: p.Name, ObjectZMM: p.ObjectZMM}
	glassAfter := func(i int) string {
		if i < 0 || i >= len(p.Cards) {
			return ""
		}
		return p.Cards[i].Glass
	}
	// poseOfCard copies the card pose into an element.
	poseOfCard := func(e *Element, c Card) {
		e.DecenterXMM, e.DecenterYMM = c.DecenterXMM, c.DecenterYMM
		e.TiltXDeg, e.TiltYDeg, e.TiltZDeg = c.TiltXDeg, c.TiltYDeg, c.TiltZDeg
	}
	samePose := func(a, b Card) bool {
		return a.DecenterXMM == b.DecenterXMM && a.DecenterYMM == b.DecenterYMM &&
			a.TiltXDeg == b.TiltXDeg && a.TiltYDeg == b.TiltYDeg && a.TiltZDeg == b.TiltZDeg
	}
	i := 0
	for i < len(p.Cards) {
		c := p.Cards[i]
		if isAirName(glassAfter(i)) {
			// A single-surface element.
			e := Element{Name: c.Name, Front: surfaceOfCard(c), GapGlass: glassAfter(i), GapMM: c.ThicknessMM}
			poseOfCard(&e, c)
			switch {
			case c.Mirror:
				e.Kind = ElementMirror
			case i == p.Stop && e.Front.isPlane() && e.Front.asphereKinds() == 0:
				e.Kind = ElementStop
				e.Front = SurfaceSpec{SemiAperMM: c.SemiAperMM}
				e.GapGlass = ""
			default:
				e.Kind = ElementLens
			}
			if i == p.Stop {
				g.Stop = len(g.Elements)
			}
			if _, err := g.AddElement(e); err != nil {
				return nil, err
			}
			i++
			continue
		}
		// A glass block: cards i .. j (j = the first card that ends in air).
		j := i
		for j+1 < len(p.Cards) && !isAirName(glassAfter(j)) {
			j++
		}
		for k := i; k <= j; k += 2 {
			if k+1 <= j {
				a, b := p.Cards[k], p.Cards[k+1]
				if b.Mirror {
					return nil, fmt.Errorf("lens: card %d is reflective inside the glass between cards %d and %d: a reflective surface inside a glass block (a Mangin mirror) is not representable in the element model", k+1, i, j)
				}
				if !samePose(a, b) {
					return nil, fmt.Errorf("lens: cards %d and %d carry different poses: the element model places both surfaces with one pose", k, k+1)
				}
				if (a.TiltXDeg != 0 || a.TiltYDeg != 0) && a.ThicknessMM != 0 {
					return nil, fmt.Errorf("lens: cards %d and %d are tilted by (%g, %g) degrees over a %g mm thickness: the card model measures that thickness along the chain axis, the element model along the element axis; use two poses or a zero thickness",
						k, k+1, a.TiltXDeg, a.TiltYDeg, a.ThicknessMM)
				}
				if k+1 == p.Stop {
					return nil, fmt.Errorf("lens: the stop is card %d, the rear surface of the element formed by cards %d and %d: the element model stops at the front surface of an element", p.Stop, k, k+1)
				}
				rear := surfaceOfCard(b)
				e := Element{
					Name: a.Name,
					Kind: ElementLens, Front: surfaceOfCard(a), Rear: &rear,
					Glass: a.Glass, GapGlass: b.Glass,
					ThicknessMM: a.ThicknessMM, GapMM: b.ThicknessMM,
				}
				poseOfCard(&e, a)
				if k == p.Stop {
					g.Stop = len(g.Elements)
				}
				if _, err := g.AddElement(e); err != nil {
					return nil, err
				}
				continue
			}
			// The odd card that closes the block.
			b := p.Cards[k]
			e := Element{Name: b.Name, Front: surfaceOfCard(b), GapGlass: glassAfter(k), GapMM: b.ThicknessMM}
			poseOfCard(&e, b)
			if b.Mirror {
				e.Kind = ElementMirror
			} else {
				e.Kind = ElementLens
			}
			if k == p.Stop {
				g.Stop = len(g.Elements)
			}
			if _, err := g.AddElement(e); err != nil {
				return nil, err
			}
		}
		i = j + 1
	}
	if err := g.Validate(); err != nil {
		return nil, err
	}
	return g, nil
}
