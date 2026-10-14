package geo

import (
	"fmt"
	"math"
)

// Aperture describes the clear aperture of a surface in the surface's local
// transverse (x, y) plane. Contains decides whether a ray that hit the surface
// at (x, y) passes the aperture; MaxRadius returns the largest distance from
// the aperture centre to a point of the aperture and is used as the scale for
// pupil sampling, vignetting diagnostics and surface search windows.
//
// The aperture of a surface is centred on the vertex; a nil Aperture on a
// Surface means "no aperture" (unlimited), see NoneAperture.
type Aperture interface {
	// Contains reports whether the local point (x, y) lies inside the clear
	// aperture (boundary included).
	Contains(x, y float64) bool
	// MaxRadius returns the maximum radial extent of the aperture (mm).
	MaxRadius() float64
	// String describes the aperture.
	String() string
}

// NoneAperture is an infinite aperture: every point passes. The zero value is
// valid and is what a nil Surface.Aperture behaves like.
type NoneAperture struct{}

// Contains always returns true.
func (NoneAperture) Contains(x, y float64) bool { return true }

// MaxRadius returns +Inf.
func (NoneAperture) MaxRadius() float64 { return math.Inf(1) }

// String implements Aperture.
func (NoneAperture) String() string { return "infinite aperture" }

// CircleAperture is a centred circular clear aperture of radius Radius (mm).
type CircleAperture float64

// Contains reports whether x²+y² <= Radius².
func (a CircleAperture) Contains(x, y float64) bool {
	r := float64(a)
	return r > 0 && x*x+y*y <= r*r
}

// MaxRadius returns the radius (mm).
func (a CircleAperture) MaxRadius() float64 { return float64(a) }

// String implements Aperture.
func (a CircleAperture) String() string { return fmt.Sprintf("circle r=%g mm", float64(a)) }

// RectangleAperture is a centred rectangular clear aperture with half-width
// HalfWidth (along x) and half-height HalfHeight (along y), in mm.
type RectangleAperture struct {
	HalfWidth  float64
	HalfHeight float64
}

// Contains reports whether |x| <= HalfWidth and |y| <= HalfHeight.
func (a RectangleAperture) Contains(x, y float64) bool {
	return a.HalfWidth > 0 && a.HalfHeight > 0 &&
		math.Abs(x) <= a.HalfWidth && math.Abs(y) <= a.HalfHeight
}

// MaxRadius returns the half diagonal sqrt(hw²+hh²), i.e. the farthest corner.
func (a RectangleAperture) MaxRadius() float64 {
	return math.Hypot(a.HalfWidth, a.HalfHeight)
}

// String implements Aperture.
func (a RectangleAperture) String() string {
	return fmt.Sprintf("rectangle %.4g x %.4g mm", 2*a.HalfWidth, 2*a.HalfHeight)
}

// AnnulusAperture is a centred ring-shaped clear aperture with Inner and Outer
// radii in mm (Inner < Outer).
type AnnulusAperture struct {
	Inner float64
	Outer float64
}

// Contains reports whether Inner² <= x²+y² <= Outer².
func (a AnnulusAperture) Contains(x, y float64) bool {
	r2 := x*x + y*y
	return a.Outer > 0 && a.Outer >= a.Inner && r2 <= a.Outer*a.Outer && r2 >= a.Inner*a.Inner
}

// MaxRadius returns the outer radius (mm).
func (a AnnulusAperture) MaxRadius() float64 { return a.Outer }

// String implements Aperture.
func (a AnnulusAperture) String() string {
	return fmt.Sprintf("annulus r=%g..%g mm", a.Inner, a.Outer)
}

// apertureOf returns a usable Aperture for a possibly nil interface value.
func apertureOf(a Aperture) Aperture {
	if a == nil {
		return NoneAperture{}
	}
	return a
}
