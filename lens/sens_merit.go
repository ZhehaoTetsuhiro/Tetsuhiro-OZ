package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// MeritKind selects the image-quality measure of a MeritSpec.
type MeritKind int

const (
	// MeritSpotRMS is the RMS radial spot size (mm) about each field's
	// centroid of hits, pooled over all fields and pupil samples.
	MeritSpotRMS MeritKind = iota
	// MeritWavefrontRMS is the RMS optical path difference (mm) of the exit
	// wavefront relative to the reference sphere centred on the field's
	// centroid of hits (geo.OPDGrid), piston removed, pooled over fields.
	MeritWavefrontRMS
)

// String implements fmt.Stringer.
func (k MeritKind) String() string {
	switch k {
	case MeritSpotRMS:
		return "spot RMS"
	case MeritWavefrontRMS:
		return "wavefront RMS"
	}
	return fmt.Sprintf("MeritKind(%d)", int(k))
}

// MeritSpec defines a scalar merit of a prescription,
//
//	M(p) = sqrt( Σ_k r_k(p)² / N ),
//
// where N is the number of sampled rays over all fields and r_k the per-ray
// residuals of Kind (see Residuals). M is a root-mean-square measure in mm.
//
// The image plane is either fixed at ImageZ (the detector of a built design:
// a perturbation that moves the focus then degrades M) or floats to the
// paraxial image of each evaluated system (FloatingFocus: the design-time
// convention, where focus is always re-found). Fixed returns a copy with
// ImageZ set to the nominal focus of a prescription; every tolerance analysis
// in this package starts from such a fixed merit.
type MeritSpec struct {
	Kind          MeritKind
	LambdaNm      float64  // wavelength (nm); 0 = geo.Wavelengthd
	Fields        []Field  // nil = on-axis and 0.5° (the design-driver default)
	Samples       []Sample // nil = UnitDiskSamples(3, 6)
	ImageZ        float64  // fixed image plane (mm); used unless FloatingFocus
	FloatingFocus bool
}

// withDefaults returns a copy with every zero field replaced by its default.
func (m MeritSpec) withDefaults() MeritSpec {
	if m.LambdaNm <= 0 {
		m.LambdaNm = geo.Wavelengthd
	}
	if m.Fields == nil {
		m.Fields = []Field{{0, 0}, {0, 0.0087}}
	}
	if m.Samples == nil {
		m.Samples = UnitDiskSamples(3, 6)
	}
	return m
}

// Fixed returns a copy of the merit whose image plane is fixed at the nominal
// paraxial image of p, with FloatingFocus cleared.
func (m MeritSpec) Fixed(p *Prescription) (MeritSpec, error) {
	m = m.withDefaults()
	par, err := p.Paraxial(m.LambdaNm)
	if err != nil {
		return m, err
	}
	if math.IsInf(par.ImageZ, 0) || math.IsNaN(par.ImageZ) {
		return m, fmt.Errorf("lens: no finite paraxial image to fix the merit plane (z=%v)", par.ImageZ)
	}
	m.ImageZ = par.ImageZ
	m.FloatingFocus = false
	return m, nil
}

// Len returns the length of the residual vector of the merit.
func (m MeritSpec) Len() int {
	m = m.withDefaults()
	n := len(m.Fields) * len(m.Samples)
	if m.Kind == MeritSpotRMS {
		return 2 * n
	}
	return n
}

// Residuals returns the residual vector of p under the merit (mm): for the
// spot kind two entries per ray (x and y displacement from the field's
// centroid of hits), for the wavefront kind one entry per ray (OPD). A ray
// that does not reach the image plane contributes failPenalty entries, so
// the length is fixed for a given MeritSpec (see Len).
func (m MeritSpec) Residuals(p *Prescription) ([]float64, error) {
	m = m.withDefaults()
	if m.Kind != MeritSpotRMS && m.Kind != MeritWavefrontRMS {
		return nil, fmt.Errorf("lens: unknown merit kind %d", int(m.Kind))
	}
	if !m.FloatingFocus && m.ImageZ == 0 {
		return nil, fmt.Errorf("lens: merit image plane is unset (call Fixed or set FloatingFocus)")
	}
	base, err := p.System()
	if err != nil {
		return nil, err
	}
	par, err := base.Paraxial(m.LambdaNm, p.ObjectZMM)
	if err != nil {
		return nil, err
	}
	imageZ := m.ImageZ
	if m.FloatingFocus {
		imageZ = par.ImageZ
	}
	if math.IsInf(imageZ, 0) || math.IsNaN(imageZ) {
		return nil, fmt.Errorf("lens: merit image plane not finite (z=%v)", imageZ)
	}
	out := make([]float64, 0, m.Len())
	for _, f := range m.Fields {
		switch m.Kind {
		case MeritSpotRMS:
			out, err = appendSpotResiduals(out, base, par, p.ObjectZMM, f, m.Samples, m.LambdaNm, imageZ)
		case MeritWavefrontRMS:
			out = appendOPDResiduals(out, base, p.ObjectZMM, f, m.Samples, m.LambdaNm, imageZ)
		}
		if err != nil {
			return nil, err
		}
	}
	return out, nil
}

// Value returns M(p) in millimetres.
func (m MeritSpec) Value(p *Prescription) (float64, error) {
	r, err := m.Residuals(p)
	if err != nil {
		return 0, err
	}
	var s float64
	for _, v := range r {
		s += v * v
	}
	return math.Sqrt(s / float64(m.withDefaults().rays())), nil
}

// Func adapts the merit to the ScalarFunc signature of Sensitivity.
func (m MeritSpec) Func() ScalarFunc {
	return func(p *Prescription) (float64, error) { return m.Value(p) }
}

// rays returns the number of sampled rays over all fields.
func (m MeritSpec) rays() int { return len(m.Fields) * len(m.Samples) }

// appendSpotResiduals appends the centroid-referenced spot displacements of one
// field to out (two entries per sample, failPenalty for lost rays).
func appendSpotResiduals(out []float64, base *geo.System, par *geo.Paraxial, objZ float64, f Field, samples []Sample, lambdaNm, imageZ float64) ([]float64, error) {
	sys := base.WithImagePlane(imageZ)
	n := len(samples)
	xs := make([]float64, n)
	ys := make([]float64, n)
	ok := make([]bool, n)
	var mx, my float64
	nOK := 0
	for i, s := range samples {
		r, err := launchRay(par, objZ, f, s)
		if err != nil {
			return out, err
		}
		r.WavelengthNm = lambdaNm
		res := sys.Trace(r)
		if !res.OK() {
			continue
		}
		pt := res.Final()
		xs[i], ys[i], ok[i] = pt.X, pt.Y, true
		mx += pt.X
		my += pt.Y
		nOK++
	}
	if nOK > 0 {
		mx /= float64(nOK)
		my /= float64(nOK)
	}
	for i := range samples {
		if !ok[i] {
			out = append(out, failPenalty, failPenalty)
			continue
		}
		out = append(out, xs[i]-mx, ys[i]-my)
	}
	return out, nil
}

// appendOPDResiduals appends the wavefront OPD of one field to out (one entry
// per sample, failPenalty for rays without a valid reference-sphere crossing).
func appendOPDResiduals(out []float64, base *geo.System, objZ float64, f Field, samples []Sample, lambdaNm, imageZ float64) []float64 {
	gf := geoField(objZ, f)
	pts := make([]geo.PupilPoint, len(samples))
	for i, s := range samples {
		pts[i] = geo.PupilPoint{PX: s.Px, PY: s.Py}
	}
	pup, err := base.EntrancePupil(gf, lambdaNm)
	if err == nil {
		var g *geo.OPDGrid
		g, err = base.OPDGrid(gf, pup, pts, imageZ, lambdaNm, geo.OPDOptions{
			Reference:    geo.CentroidReference(),
			RemovePiston: true,
		})
		if err == nil {
			for _, v := range g.OPD {
				if math.IsNaN(v) {
					out = append(out, failPenalty)
				} else {
					out = append(out, v)
				}
			}
			return out
		}
	}
	for range samples {
		out = append(out, failPenalty)
	}
	return out
}

// geoField converts a lens field into the geo convention for the object
// position objZ (−Inf for an object at infinity).
func geoField(objZ float64, f Field) geo.Field {
	if math.IsInf(objZ, 0) {
		return geo.FieldTan(f.TanX, f.TanY)
	}
	return geo.FieldHeight(objZ, f.TanX, f.TanY)
}

// EFLScalar returns the paraxial effective focal length (mm) at lambdaNm as a
// ScalarFunc, the first-order quantity of the analytic checks.
func EFLScalar(lambdaNm float64) ScalarFunc {
	return func(p *Prescription) (float64, error) {
		par, err := p.Paraxial(lambdaNm)
		if err != nil {
			return 0, err
		}
		return par.EFL, nil
	}
}
