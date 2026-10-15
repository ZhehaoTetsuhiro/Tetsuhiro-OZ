// Package geo implements the geometrical-optics core of Tetsuhiro-OZ: surface
// shapes, vector ray tracing through sequential systems, paraxial (first
// order) optics, pupil/field sampling, spot diagrams, transverse aberration
// and OPD analysis, third-order (Seidel) aberrations including an aspheric
// contribution (Petzval curvature and distortion included) and the paraxial
// chromatic aberrations (axial and lateral colour).
//
// The package is pure Go standard library, float64 throughout, and is meant to
// be the numerical foundation that the lens/ and coating/ packages (and the
// wave-optics kernel of Tetsuhiro-WOS) are validated against.
//
// # Coordinate system
//
// A right-handed Cartesian frame is used:
//
//	x  points to the right (the "sagittal" / horizontal direction),
//	y  points up (the "meridional" / vertical direction),
//	z  is the optical axis, pointing from the object towards the image
//	   ("left to right" when the system is drawn on paper).
//
// The sign convention for the fields that define the object is
//
//	field angle  (fx, fy) = (tan θx, tan θy): a source above the axis
//	             (θy > 0) sends light travelling in the direction
//	             d = (-fx, -fy, 1)/|...| , i.e. downwards and to the right;
//	object height: the object point lives at z = ObjectZ < 0 (in front of
//	             the first surface, whose vertex is at z = 0).
//
// # Units
//
//	length      millimetre (mm): vertex distances, thicknesses, radii,
//	            apertures, spot radii, OPD values;
//	wavelength  nanometre (nm): every Wavelength field/method is in nm;
//	angles      degree for field angles in the public constructors, radian
//	            for everything returned (half angles, NA, etc.);
//	refractive indices are dimensionless.
//
// # Decentre and tilt (coordinate breaks)
//
// Surfaces may be decentred and tilted with a CoordBreak, the convention of
// Zemax's COORDBRK surface: the origin is first moved by (Dx, Dy, Dz) in the
// enclosing frame, then the axes are tilted about x, about the tilted y, and
// about the tilted z (roll). The rotation is Q = Rx(TiltX)·Ry(TiltY)·Rz(TiltZ)
// with right-handed active rotations, and a point of the new frame is mapped to
// the enclosing frame by
//
//	p_old = d + Q·p_new ,   d = (Dx, Dy, Dz) .
//
// Two ways of using a break are provided:
//
//   - Surface.Pose decentres or tilts one surface. The surface's local frame is
//     the break applied to the chain frame of its position; the next surface is
//     placed on the chain again (auto-return to the optical axis). Thickness[i]
//     is measured along the axis of the chain, not along the tilted surface.
//   - Surface.Break = true makes a pure coordinate-break element: no profile,
//     no interaction, no Step in a trace. It changes the chain for the following
//     elements. A return break is CoordBreak.Inverse() placed directly after it
//     with zero spacing; the pair is the identity (tested to 1e-12 mm).
//
// Analysis under poses: the trace, spot, ray fan and OPD routines are exact in
// the global frame. Their chief reference is the real chief ray (the ray through
// the centre of the aperture stop, System.ChiefRay) whenever the system has a
// break or a pose; the paraxial engine, the Seidel sums and the chromatic
// routines return an error for a decentred or tilted layout, Sub returns an
// error and Reversed returns nil. For such layouts set Pupil.Offset to the
// crossing of the real chief ray (System.ChiefCrossing) to centre the sampling
// on it.
//
// # Surfaces and sign conventions
//
// A surface is described in its own local frame: the vertex is at the local
// origin and the local optical axis is the surface's z axis (the global +z
// axis unless the surface carries a pose). The rotationally symmetric sag is
// z = Sag(r) with r = sqrt(x²+y²) and
//
//   - R > 0: the centre of curvature lies on the +z side of the vertex, i.e.
//     the surface bulges towards the incoming light (convex for a normal
//     left-to-right layout);
//   - R < 0: the centre of curvature lies on the -z side (concave towards the
//     incoming light, "bulging away");
//   - R = ±Inf (or Plane): a flat surface, z ≡ 0.
//
// The unit surface normal returned by Normal is oriented towards +z (the
// normal's z component is always > 0). Ray–surface interaction routines flip
// it as needed, so an incoming ray may hit either face of a surface.
//
// # Freeform surfaces (XY polynomial)
//
// Next to the rotationally symmetric profiles (Plane, Conic, Sphere,
// EvenAsphere, PolyAsphere, QAsphere) the package provides the freeform
// profile XYPolynomial, whose sag in the local frame is
//
//	z(x, y) = conicBase(r) + Σ a_k (x/N)^i_k (y/N)^j_k ,   r = √(x²+y²),
//
// with the optional conic base of revolution (R, K) of Conic, integer
// exponents i_k, j_k ≥ 0 with i_k + j_k ≥ 2 (no piston, no tilt: the vertex is
// at the origin with a horizontal tangent plane) and an optional normalisation
// radius N (N > 0 normalises the coordinates so that the coefficients are sags
// in mm; N = 0 uses the raw coordinates). Terms are summed in a canonical
// order (ascending i+j, then i, then j; see XYPolynomial.Canonical), so the
// evaluation is deterministic bit for bit. The analytic normal is
// n ∝ (−∂z/∂x, −∂z/∂y, 1) and the ray intersection is a Newton iteration with a
// bracketed-bisection fallback that uses the two-dimensional sag — no
// rotational symmetry is assumed anywhere in the solver, so the surface is
// valid off the axis and far from its vertex.
//
// The type implements the FreeformProfile interface; Surface and System
// dispatch to it for the two-dimensional shape (Surface.Sag, Surface.Normal,
// Surface.Intersect, Surface.ZRange), and its Profile methods report the +y
// meridian (x = 0) as the rotationally symmetric surrogate.
//
// Applicability. The paraxial (first-order) engine and the third-order
// (Seidel) sums of the package are defined for rotationally symmetric coaxial
// systems only. A freeform surface whose departure is not a polynomial in r²
// is rejected by them (System.requireCoaxial, like a decentred or tilted
// surface) and has no meaningful first-order description; genuinely freeform
// systems must be analysed with real rays (Trace, SpotDiagram, RayFan, OPDGrid,
// geometric MTF), usually with Surface.Pose / Surface.Break coordinate breaks,
// since freeform elements are almost always decentred or tilted. A
// rotationally symmetric XYPolynomial (XYPolynomial.RotationallySymmetric)
// degenerates to an ordinary rotationally symmetric profile and stays usable
// with the first-order engine.
//
// # Status model
//
// Outcomes of a ray–surface interaction are reported as a Status value
// (StatusOK, StatusHit, StatusVignetted, StatusTIR, StatusEscape,
// StatusDegenerate) rather than by panicking: a vignetted ray, a ray that
// misses a surface and total internal reflection are legitimate physical
// outcomes, not program errors. Trace always returns a TraceResult whose
// Status (and FailedAt index) describes where and why tracing stopped; the
// partial path up to that point remains available in the recorded steps.
//
// # Determinism
//
// All iterative solvers are deterministic: they use fixed iteration counts,
// fixed bracketing grids and fixed tie-breaks (smallest valid intersection
// parameter), so identical inputs always give bit-identical outputs.
package geo
