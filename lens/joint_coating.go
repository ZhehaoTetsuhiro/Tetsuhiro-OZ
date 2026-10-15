package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/geo"
)

// Coating is a thin-film stack attached to one interface (surface) of a
// prescription. Layers run from the incidence side (the medium in front of the
// surface) to the substrate side (the medium behind it), as in coating.Stack.
// The two bounding media are not stored: they are the geometric media on either
// side of the surface, so the same layers apply to any prescription. An empty
// Layers list is an uncoated interface, which the transfer-matrix evaluation
// reduces to the bare Fresnel interface.
//
// Coatings are kept outside the Prescription (one entry per surface index), so
// the glass-name model of the prescription is unchanged.
type Coating struct {
	Name   string
	Layers []coating.Layer
}

// mediumDispersion presents a geo medium to the coating package as a lossless
// dispersion model with the same refractive index.
type mediumDispersion struct{ m geo.Medium }

// IndexNm implements coating.Dispersion.
func (d mediumDispersion) IndexNm(lambdaNm float64) complex128 {
	return complex(d.m.N(lambdaNm), 0)
}

// Lossless implements coating.Dispersion; geo media carry no extinction.
func (d mediumDispersion) Lossless() bool { return true }

// Model implements coating.Dispersion.
func (d mediumDispersion) Model() string { return "geo medium " + d.m.Name() }

// coatingMaterial wraps a geo medium as a coating material.
func coatingMaterial(m geo.Medium) coating.Material {
	return coating.NewMaterial(m.Name(), mediumDispersion{m: m})
}

// stackAt binds the layers of c to the media before and after a surface.
func (c Coating) stackAt(before, after geo.Medium) coating.Stack {
	return coating.Stack{
		Incident:  coatingMaterial(before),
		Substrate: coatingMaterial(after),
		Layers:    c.Layers,
	}
}

// CoatedInterface returns the power reflectance and transmittance of the
// interface between the media before and after, carrying coating c, at
// wavelength lambdaNm and incidence angle angleDeg (from the surface normal).
// The values are the s and p results of the transfer-matrix evaluation of
// package coating (coating.Stack.Evaluate), so they coincide with a direct
// coating-package computation of the same stack.
func CoatedInterface(c Coating, before, after geo.Medium, lambdaNm, angleDeg float64) (rs, rp, ts, tp float64, err error) {
	st := c.stackAt(before, after)
	s, err := st.Evaluate(lambdaNm, angleDeg, coating.PolS)
	if err != nil {
		return 0, 0, 0, 0, err
	}
	p, err := st.Evaluate(lambdaNm, angleDeg, coating.PolP)
	if err != nil {
		return 0, 0, 0, 0, err
	}
	return s.R, p.R, s.T, p.T, nil
}

// InterfaceHit is the interaction of one traced ray with one interface.
type InterfaceHit struct {
	Surface int     // surface (card) index
	IncDeg  float64 // incidence angle of the ray on the surface (degrees)
	Rs, Rp  float64 // power reflectance of the interface, s and p
	Ts, Tp  float64 // power transmittance of the interface, s and p
}

// RayTransmission is the outcome of TraceCoated for one ray.
type RayTransmission struct {
	Passed bool       // the ray traversed every surface (not vignetted, TIR or escaped)
	Status geo.Status // status of the trace (StatusOK when Passed)
	Hits   []InterfaceHit
	// TS and TP are the single-pass transmittances of the ray for s and p
	// polarization: the product of the interface transmittances along its
	// path. They are 0 when the ray does not pass. Reflections between surfaces
	// (ghost paths) are not included; this is the single-pass convention of the
	// system transmittance.
	TS, TP float64
}

// T returns the transmittance of the ray for polarization pol. For
// coating.PolUnpolarized the s and p powers are averaged, which is exact for
// unpolarized light because the two components do not interfere. For an
// unknown polarization it returns NaN.
func (t RayTransmission) T(pol coating.Polarization) float64 {
	switch pol {
	case coating.PolS:
		return t.TS
	case coating.PolP:
		return t.TP
	case coating.PolUnpolarized:
		return 0.5 * (t.TS + t.TP)
	}
	return math.NaN()
}

// TraceCoated traces ray r through the (uncoated) system sys and evaluates
// every interface it meets with the coating in coatings[surface] (missing or
// empty entries are bare interfaces). For each interaction the incidence angle
// is the one of the traced ray, so R and T follow the actual geometry.
//
// Polarization is carried as two label channels, s and p, each evaluated with
// the incidence angle of the ray; the two products are kept separately. This is
// exact for planar (meridional) ray families, where the plane of incidence is
// the same at every surface. For skew rays the s/p basis rotates from surface
// to surface and the label model is an approximation: the polarization coupling
// of skew rays belongs to the Jones tracing (G5), which is not part of this file.
//
// Mirror surfaces are not part of the transmittance model and are rejected.
func TraceCoated(sys *geo.System, r geo.Ray, coatings []Coating) (RayTransmission, error) {
	res := sys.Trace(r)
	out := RayTransmission{Passed: res.OK(), Status: res.Status, TS: 1, TP: 1}
	for _, st := range res.Steps {
		if st.Status != geo.StatusOK {
			continue
		}
		if st.Mirror {
			return RayTransmission{}, fmt.Errorf("lens: mirror surface %d is outside the transmittance model", st.Surface)
		}
		var c Coating
		if st.Surface < len(coatings) {
			c = coatings[st.Surface]
		}
		before := sys.Mediums[st.Surface]
		after := sys.Mediums[st.Surface+1]
		angle := st.IncAngle * 180 / math.Pi
		rs, rp, ts, tp, err := CoatedInterface(c, before, after, r.WavelengthNm, angle)
		if err != nil {
			return RayTransmission{}, fmt.Errorf("lens: surface %d at %g nm, %.4g°: %w", st.Surface, r.WavelengthNm, angle, err)
		}
		out.Hits = append(out.Hits, InterfaceHit{Surface: st.Surface, IncDeg: angle, Rs: rs, Rp: rp, Ts: ts, Tp: tp})
		out.TS *= ts
		out.TP *= tp
	}
	if !out.Passed {
		out.TS, out.TP = 0, 0
	}
	return out, nil
}
