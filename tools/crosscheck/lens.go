package main

import (
	"encoding/json"
	"fmt"
	"io"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// lensMarginalFrac is the entrance-pupil fraction of the marginal ray: the
// exact rim of the entrance pupil.
const lensMarginalFrac = 1.0

// lensSurface is one surface of the exported prescription. The rayoptics side
// rebuilds its model from exactly these numbers.
type lensSurface struct {
	Name      string  `json:"name"`
	Radius    float64 `json:"radius_mm"`    // vertex radius, 0 for a plane
	Thickness float64 `json:"thickness_mm"` // to the next surface, 0 for the last
	Nd        float64 `json:"nd"`           // index of the medium behind the surface
	Vd        float64 `json:"vd"`           // Abbe number (0 for air)
	Semi      float64 `json:"semi_mm"`      // clear semi-aperture, 0 if unlimited
}

// lensExport is the JSON contract read by lenscheck.py. Lengths are mm, the
// object is at infinity and the wavelength is geo.Wavelengthd.
type lensExport struct {
	Name           string        `json:"name"`
	Wavelength     float64       `json:"wavelength_nm"`
	ObjectZ        string        `json:"object_z_mm"`
	StopIndex      int           `json:"stop_index"`
	StopSemi       float64       `json:"stop_semi_mm"`
	Surfaces       []lensSurface `json:"surfaces"`
	EFL            float64       `json:"efl_mm"`
	BFL            float64       `json:"bfl_mm"`
	FNumber        float64       `json:"fno"`
	EPRadius       float64       `json:"epr_mm"`
	FrontPrincipal float64       `json:"front_principal_z_mm"`
	BackPrincipal  float64       `json:"back_principal_z_mm"`
	VertexLastZ    float64       `json:"vertex_last_z_mm"`
	MarginalFrac   float64       `json:"marginal_pupil_fraction"`
	MarginalHgt    float64       `json:"marginal_height_mm"`
	MarginalStatus string        `json:"marginal_status"`
	MarginalY      float64       `json:"marginal_y_at_focus_mm"`
}

// exportLens computes the first-order data and the marginal-ray intercept of
// an object-at-infinity prescription with the OZ kernels.
//
// The marginal ray is a parallel ray at height frac·EPR, launched in front of
// the system, traced by real ray tracing and intersected with the paraxial
// focal plane z = VertexZ(last) + BFL. That is the definition of the marginal
// ray for an infinite conjugate; a finite-distance launch would add a focus
// error of order f²/d and is therefore not used on either side.
func exportLens(p *lens.Prescription) (lensExport, error) {
	if !math.IsInf(p.ObjectZMM, -1) {
		return lensExport{}, fmt.Errorf("lenscheck: object at z=%g; only an object at infinity is cross-checked", p.ObjectZMM)
	}
	sys, err := p.System()
	if err != nil {
		return lensExport{}, err
	}
	par, err := sys.Paraxial(geo.Wavelengthd, math.Inf(-1))
	if err != nil {
		return lensExport{}, err
	}
	n := len(sys.Surfaces)
	stop := sys.StopIndex
	if !(par.EntrancePupilRadius > 0) || math.IsInf(par.EntrancePupilRadius, 0) {
		return lensExport{}, fmt.Errorf("lenscheck: entrance pupil radius %g is not finite and positive", par.EntrancePupilRadius)
	}
	out := lensExport{
		Name:           p.Name,
		Wavelength:     geo.Wavelengthd,
		ObjectZ:        "inf",
		StopIndex:      stop,
		StopSemi:       p.Cards[stop].SemiAperMM,
		EFL:            par.EFL,
		BFL:            par.BFL,
		FNumber:        par.FNumber,
		EPRadius:       par.EntrancePupilRadius,
		FrontPrincipal: par.FrontPrincipalZ,
		BackPrincipal:  par.BackPrincipalZ,
		VertexLastZ:    par.VertexZ[n-1],
		MarginalFrac:   lensMarginalFrac,
	}
	for i, c := range p.Cards {
		m := sys.Mediums[i+1]
		s := lensSurface{Name: c.Name, Nd: m.N(geo.Wavelengthd), Semi: c.SemiAperMM}
		if !math.IsInf(c.RadiusMM, 0) {
			s.Radius = c.RadiusMM
		}
		if g, ok := m.(geo.Sellmeier); ok {
			s.Vd = g.AbbeNumber()
		}
		if i < n-1 {
			s.Thickness = c.ThicknessMM
		}
		out.Surfaces = append(out.Surfaces, s)
	}

	// The marginal ray is traced without clear apertures, as rayoptics does
	// with apply_vignetting=False. A real ray launched at the entrance-pupil
	// rim does not in general pass the stop rim once aberrations act, so a
	// stop clip there would reflect the definition of the marginal ray rather
	// than the system. The paraxial data above keep the physical apertures.
	free := &geo.System{Surfaces: append([]geo.Surface(nil), sys.Surfaces...),
		Thickness: sys.Thickness, Mediums: sys.Mediums, StopIndex: sys.StopIndex}
	for i := range free.Surfaces {
		free.Surfaces[i].Aperture = geo.NoneAperture{}
	}
	zmin, _ := free.BoundsZ(par.EntrancePupilRadius)
	h := lensMarginalFrac * par.EntrancePupilRadius
	res := free.Trace(geo.NewRay(geo.V(0, h, zmin-1), geo.V(0, 0, 1), geo.Wavelengthd))
	out.MarginalHgt = h
	if !res.OK() {
		out.MarginalStatus = res.Status.String()
		return out, nil
	}
	d := res.Ray.Dir
	t := (par.BackFocalZ - res.Ray.Origin.Z) / d.Z
	out.MarginalStatus = "ok"
	out.MarginalY = res.Ray.Origin.Y + d.Y*t
	return out, nil
}

// runLens loads a prescription, exports its reference data and writes it as
// JSON to w.
func runLens(path string, w io.Writer) error {
	p, err := lens.LoadJSON(path)
	if err != nil {
		return err
	}
	out, err := exportLens(p)
	if err != nil {
		return err
	}
	enc := json.NewEncoder(w)
	enc.SetIndent("", " ")
	return enc.Encode(out)
}
