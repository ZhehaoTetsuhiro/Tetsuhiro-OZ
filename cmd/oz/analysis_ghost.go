package main

// Ghost-path wiring (G10): `oz geo ghost` and POST /api/geo/ghost call
// runGeoGhost below, so the command line and the service share one
// implementation (see analysis.go for the other analyses).
//
// The analysis itself is geo.System.GhostAnalysis: it enumerates the
// double-reflection (ghost) paths of a system, estimates the relative energy
// R_I·R_J of each pair from the normal-incidence Fresnel reflectance, sorts the
// paths by energy and reports the paraxial landing point and the real-ray
// scattered spot of every kept path.

import (
	"errors"
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

// geoGhostRequest is the input of runGeoGhost.
type geoGhostRequest struct {
	Prescription  lens.Prescription `json:"prescription"`
	LambdaNm      float64           `json:"lambda_nm"`       // 0 = d-line
	FieldDeg      float64           `json:"field_deg"`       // meridional field angle
	ImageZMM      float64           `json:"image_z_mm"`      // 0 = paraxial image plane
	PupilRadiusMM float64           `json:"pupil_radius_mm"` // 0 = 1e-3 of the entrance pupil radius
	Rings         int               `json:"rings"`           // sampling rings of Hexapolar (0 = 1)
	MinRelative   float64           `json:"min_relative"`    // energy threshold (0 keeps everything)
	MaxPaths      int               `json:"max_paths"`       // top-N after sorting (0 keeps all)
	// Reflectance overrides the per-surface intensity reflectance (len =
	// surfaces); null entries fall back to the normal-incidence Fresnel value.
	Reflectance []*float64 `json:"reflectance"`
}

// geoGhostRow is one ghost path of the table.
type geoGhostRow struct {
	Pair     string  `json:"pair"` // "reflect(I=2, J=1)"
	I        int     `json:"i"`    // first reflection surface (further along the axis)
	J        int     `json:"j"`    // second reflection surface (J < I)
	Energy   float64 `json:"relative_energy"`
	ReflectI float64 `json:"reflect_i"`
	ReflectJ float64 `json:"reflect_j"`
	// ParaxialXMM/YMM is the paraxial landing point on the image plane; null
	// for a path the first-order engine cannot handle (decentred/tilted).
	ParaxialXMM     *float64 `json:"paraxial_x_mm"`
	ParaxialYMM     *float64 `json:"paraxial_y_mm"`
	ChiefXMM        float64  `json:"chief_x_mm"`
	ChiefYMM        float64  `json:"chief_y_mm"`
	ChiefStatus     string   `json:"chief_status"`
	SpotN           int      `json:"spot_rays"`
	SpotNTotal      int      `json:"spot_rays_total"`
	SpotCentroidXMM float64  `json:"spot_centroid_x_mm"`
	SpotCentroidYMM float64  `json:"spot_centroid_y_mm"`
	SpotRMSMM       float64  `json:"spot_rms_mm"`
	SpotGEOM        float64  `json:"spot_geo_mm"`
	Note            string   `json:"note,omitempty"` // non-empty: path not analysed (mirror, ...)
}

// geoGhostResponse is the answer of `oz geo ghost` and POST /api/geo/ghost.
type geoGhostResponse struct {
	LambdaNm    float64       `json:"lambda_nm"`
	ImageZMM    float64       `json:"image_z_mm"`
	FieldDeg    float64       `json:"field_deg"`
	PupilZMM    float64       `json:"pupil_z_mm"`
	PupilRMM    float64       `json:"pupil_radius_mm"`
	Rings       int           `json:"rings"`
	MinRelative float64       `json:"min_relative"`
	Enumerated  int           `json:"paths_enumerated"`
	Paths       []geoGhostRow `json:"paths"`
}

// runGeoGhost analyses the double-reflection paths of a prescription.
func runGeoGhost(req geoGhostRequest) (*geoGhostResponse, error) {
	lambda := req.LambdaNm
	if lambda == 0 {
		lambda = dLineNm
	}
	if math.IsNaN(lambda) || math.IsInf(lambda, 0) || lambda <= 0 {
		return nil, fmt.Errorf("lambda_nm %v 必须为正且有限", req.LambdaNm)
	}
	if math.IsNaN(req.FieldDeg) || math.IsInf(req.FieldDeg, 0) || math.Abs(req.FieldDeg) >= 90 {
		return nil, fmt.Errorf("field_deg %v 必须是 |θ| < 90 的有限度数", req.FieldDeg)
	}
	if math.IsNaN(req.MinRelative) || req.MinRelative < 0 {
		return nil, fmt.Errorf("min_relative %v 必须是非负有限数", req.MinRelative)
	}
	rings := defI(req.Rings, 1)
	if rings < 1 || rings > maxGhostRings {
		return nil, fmt.Errorf("rings %d 超出 [1, %d]", req.Rings, maxGhostRings)
	}
	p := &req.Prescription
	if err := p.Validate(); err != nil {
		return nil, err
	}
	sys, err := p.System()
	if err != nil {
		return nil, err
	}
	if n := sys.NumSurfaces(); n*(n-1)/2 > maxGhostPaths {
		return nil, fmt.Errorf("面数 %d 的双反射路径有 %d 条，超过上限 %d", n, n*(n-1)/2, maxGhostPaths)
	}
	par, err := p.Paraxial(lambda)
	if err != nil {
		return nil, fmt.Errorf("一阶分析（对偏心/倾斜面不适用）: %w", err)
	}
	imageZ := req.ImageZMM
	if imageZ == 0 {
		if math.IsInf(par.ImageZ, 0) || math.IsNaN(par.ImageZ) {
			return nil, errors.New("近轴像面不存在，请用 image_z_mm 指定像面")
		}
		imageZ = par.ImageZ
	}
	field := geo.FieldTan(0, math.Tan(req.FieldDeg*math.Pi/180))
	if !par.Infinite {
		field = geo.FieldHeight(par.ObjectZ, 0, math.Tan(req.FieldDeg*math.Pi/180)*(par.EntrancePupilZ-par.ObjectZ))
	}
	pup, err := sys.EntrancePupil(field, lambda)
	if err != nil {
		return nil, err
	}
	if req.PupilRadiusMM != 0 {
		if !(req.PupilRadiusMM > 0) || math.IsInf(req.PupilRadiusMM, 0) {
			return nil, fmt.Errorf("pupil_radius_mm %v 必须为正", req.PupilRadiusMM)
		}
		pup.Radius = req.PupilRadiusMM
	} else {
		// Small-aperture sampling about the chief ray: the ghost spot of the
		// reference bundle, not of the full aperture.
		pup.Radius *= 1e-3
	}
	var refl []float64
	if req.Reflectance != nil {
		if len(req.Reflectance) != sys.NumSurfaces() {
			return nil, fmt.Errorf("reflectance 有 %d 项，系统有 %d 个面", len(req.Reflectance), sys.NumSurfaces())
		}
		refl = make([]float64, len(req.Reflectance))
		for i, v := range req.Reflectance {
			if v == nil {
				refl[i] = math.NaN()
				continue
			}
			if math.IsNaN(*v) || *v < 0 || *v > 1 {
				return nil, fmt.Errorf("reflectance[%d] = %v 不在 [0,1]", i, *v)
			}
			refl[i] = *v
		}
	}
	rep, err := sys.GhostAnalysis(geo.GhostOptions{
		Field: field, Pupil: pup, Points: geo.Hexapolar(rings),
		ImageZ: imageZ, LambdaNm: lambda,
		MinRelative: req.MinRelative, MaxPaths: req.MaxPaths, Reflectance: refl,
	})
	if err != nil {
		return nil, err
	}
	out := &geoGhostResponse{
		LambdaNm: lambda, ImageZMM: imageZ, FieldDeg: req.FieldDeg,
		PupilZMM: pup.Z, PupilRMM: pup.Radius, Rings: rings,
		MinRelative: req.MinRelative, Enumerated: rep.Enumerated,
		Paths: make([]geoGhostRow, 0, len(rep.Paths)),
	}
	for _, r := range rep.Paths {
		row := geoGhostRow{
			Pair: r.Pair.String(), I: r.Pair.I, J: r.Pair.J,
			Energy: r.Energy, ReflectI: r.ReflectI, ReflectJ: r.ReflectJ,
			ChiefXMM: r.Chief.X, ChiefYMM: r.Chief.Y, ChiefStatus: r.ChiefStatus.String(),
			SpotN: r.Spot.N, SpotNTotal: r.Spot.NTotal,
			SpotCentroidXMM: r.Spot.Centroid.X, SpotCentroidYMM: r.Spot.Centroid.Y,
			SpotRMSMM: r.Spot.RMS, SpotGEOM: r.Spot.GEO, Note: r.Note,
		}
		if !math.IsNaN(r.Paraxial.X) && !math.IsNaN(r.Paraxial.Y) {
			x, y := r.Paraxial.X, r.Paraxial.Y
			row.ParaxialXMM, row.ParaxialYMM = &x, &y
		}
		out.Paths = append(out.Paths, row)
	}
	return out, nil
}
