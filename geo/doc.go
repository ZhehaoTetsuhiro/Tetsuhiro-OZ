// Package geo implements the geometrical-optics core of Tetsuhiro-OZ: surface
// shapes, vector ray tracing through sequential systems, paraxial (first
// order) optics, pupil/field sampling, spot and OPD analysis, and third-order
// (Seidel) aberrations including an aspheric contribution.
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
// # Surfaces and sign conventions
//
// A surface is described in its own local frame: the vertex is at the local
// origin and the local +z axis is the global +z axis (surfaces are currently
// centred and not tilted; see the limitations in the package docs of Surface).
// The rotationally symmetric sag is z = Sag(r) with r = sqrt(x²+y²) and
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
