package geo

import (
	"fmt"
	"math"
)

// Status reports the outcome of a ray–surface interaction or of a whole trace.
// It is a value, not an error: total internal reflection, vignetting and rays
// that miss a surface are legitimate physical outcomes.
type Status int

// Status values. A trace succeeds while the status is StatusOK; any other
// value terminates the trace at the surface where it occurred.
const (
	// StatusOK: the ray is valid and all interactions so far succeeded.
	StatusOK Status = iota
	// StatusVignetted: the ray hit the surface outside its clear aperture.
	StatusVignetted
	// StatusTIR: total internal reflection at a refracting surface.
	StatusTIR
	// StatusEscape: the ray does not intersect the (rotationally symmetric)
	// surface geometry at all, or it would have to go backwards.
	StatusEscape
	// StatusDegenerate: the ray or the surface is numerically invalid (zero
	// direction, NaN position, undefined sag).
	StatusDegenerate
	// StatusNotConverged: an iterative solver did not reach the required
	// tolerance.
	StatusNotConverged
)

// String returns the status name ("ok", "vignetted", "tir", "escape",
// "degenerate", "not-converged").
func (s Status) String() string {
	switch s {
	case StatusOK:
		return "ok"
	case StatusVignetted:
		return "vignetted"
	case StatusTIR:
		return "tir"
	case StatusEscape:
		return "escape"
	case StatusDegenerate:
		return "degenerate"
	case StatusNotConverged:
		return "not-converged"
	}
	return fmt.Sprintf("status(%d)", int(s))
}

// OK reports whether the status is StatusOK.
func (s Status) OK() bool { return s == StatusOK }

// Ray is a single geometrical ray. Origin is the current point (mm), Dir the
// unit propagation direction, WavelengthNm the vacuum wavelength (nm), OPL the
// optical path length accumulated since the ray was created (mm, sum of
// n·segment) and Path the corresponding geometrical path length (mm).
//
// Status is StatusOK for a freshly created, valid ray and is updated by
// Trace when an interaction fails (see the Status documentation).
type Ray struct {
	Origin       Vec3
	Dir          Vec3
	WavelengthNm float64
	OPL          float64
	Path         float64
	Status       Status
}

// NewRay returns a ray starting at origin with the direction dir, normalized to
// unit length. If dir is zero or not finite the ray is marked StatusDegenerate.
func NewRay(origin, dir Vec3, wavelengthNm float64) Ray {
	r := Ray{Origin: origin, WavelengthNm: wavelengthNm, Status: StatusOK}
	u, ok := dir.Unit()
	if !ok {
		r.Status = StatusDegenerate
		return r
	}
	r.Dir = u
	return r
}

// RayThrough returns a ray from a through b (both in mm).
func RayThrough(a, b Vec3, wavelengthNm float64) Ray {
	return NewRay(a, b.Sub(a), wavelengthNm)
}

// Valid reports whether the ray has a finite origin and a unit direction.
func (r Ray) Valid() bool {
	if !r.Dir.UnitValid() || r.Origin.IsNaN() {
		return false
	}
	return math.Abs(r.Dir.Norm()-1) < 1e-9
}

// At returns the point origin + t*dir of the ray line, t in mm.
func (r Ray) At(t float64) Vec3 { return r.Origin.Add(r.Dir.Mul(t)) }

// UnitValid reports whether v is a finite, normalisable vector (helper shared
// by Valid).
func (v Vec3) UnitValid() bool {
	n := v.Norm()
	return n > 0 && !math.IsInf(n, 0) && !math.IsNaN(n)
}

// PlaneDir returns the direction d restricted to the transverse plane,
// i.e. d with the z component removed.
func (v Vec3) PlaneDir() Vec3 { return Vec3{v.X, v.Y, 0} }

// String formats the ray compactly.
func (r Ray) String() string {
	return fmt.Sprintf("ray[o=%s d=%s λ=%gnm opl=%.6g %s]",
		r.Origin, r.Dir, r.WavelengthNm, r.OPL, r.Status)
}
