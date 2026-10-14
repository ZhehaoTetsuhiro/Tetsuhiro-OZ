package lens

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
)

// Field is a field point. For an object at infinity (TanX, TanY) are the
// tangents of the field angles about x and y; for a finite object they are
// the transverse object coordinates in mm (the convention of package geo).
type Field struct {
	TanX, TanY float64
}

// Sample is a normalised pupil point inside the unit disk (|Px²+Py²| <= 1).
type Sample struct{ Px, Py float64 }

// pupilInset keeps the outermost sampling ring just inside the pupil edge:
// rays that graze the rim of a curved stop surface can drift a few
// micrometres past it (field tilt x sag) or land a floating-point ulp
// outside the exact x²+y² <= r² test, and would otherwise be flagged
// vignetted, injecting artificial failPenalty terms that lock the merit
// function. 0.1% of the pupil radius is far below any physical effect.
const pupilInset = 1e-3

// UnitDiskSamples returns a deterministic quasi-hexapolar ring sampling of
// the unit disk: the centre plus `rings` rings with (ring x perRing) points
// each, with the starting angle staggered between rings. Ray-count control
// for spot residuals: len = 1 + perRing x rings x (rings+1) / 2. The
// outermost ring is inset by pupilInset (see the constant).
func UnitDiskSamples(rings, perRing int) []Sample {
	if rings < 1 {
		rings = 1
	}
	if perRing < 3 {
		perRing = 3
	}
	out := []Sample{{0, 0}}
	for r := 1; r <= rings; r++ {
		radius := float64(r) / float64(rings)
		if r == rings {
			radius *= 1 - pupilInset
		}
		n := r * perRing
		stagger := 0.0
		if r%2 == 0 {
			stagger = math.Pi / float64(n)
		}
		for k := 0; k < n; k++ {
			a := stagger + 2*math.Pi*float64(k)/float64(n)
			out = append(out, Sample{Px: radius * math.Cos(a), Py: radius * math.Sin(a)})
		}
	}
	return out
}

// failPenalty is the residual both components get when a sampled ray does
// not reach the image plane (vignetted, total internal reflection, escaped,
// or no paraxial image exists). It is large compared to typical scaled spot
// residuals, so the optimizer strongly avoids losing rays, but it is finite
// so the merit stays usable for DLS.
const failPenalty = 10.0

// launchRay builds one sampling ray of a field through the entrance pupil
// of the paraxially characterised system. For an object at infinity the
// rays are parallel to the field direction and start on the entrance-pupil
// plane; for a finite object they start at the object point and pass
// through the sampled point of the entrance-pupil plane.
func launchRay(par *geo.Paraxial, objZ float64, f Field, s Sample) (geo.Ray, error) {
	zp := par.EntrancePupilZ
	rp := par.EntrancePupilRadius
	if math.IsInf(zp, 0) || math.IsNaN(zp) {
		return geo.Ray{}, fmt.Errorf("lens: entrance pupil is at infinity (telecentric object space); sampling unsupported")
	}
	if !(rp > 0) || math.IsNaN(rp) || math.IsInf(rp, 0) {
		return geo.Ray{}, fmt.Errorf("lens: entrance pupil radius %v is unusable (set a finite stop aperture)", rp)
	}
	zStop := par.VertexZ[par.StopIndex]
	if par.Infinite {
		dir := geo.V(-f.TanX, -f.TanY, 1)
		u, ok := dir.Unit()
		if !ok {
			return geo.Ray{}, fmt.Errorf("lens: degenerate field direction")
		}
		dz := zp - zStop
		ch := geo.V(u.X/u.Z*dz, u.Y/u.Z*dz, 0)
		origin := geo.V(ch.X+s.Px*rp, ch.Y+s.Py*rp, zp)
		return geo.NewRay(origin, u, 0), nil
	}
	obj := geo.V(f.TanX, f.TanY, objZ)
	pup := geo.V(s.Px*rp, s.Py*rp, zp)
	_ = zStop
	return geo.RayThrough(obj, pup, 0), nil
}

// tracedAtImagePlane returns the system with an image plane appended at the
// paraxial image position, plus that system and the paraxial data.
func tracedAtImagePlane(p *Prescription, lambdaNm float64) (*geo.System, *geo.Paraxial, error) {
	base, err := p.System()
	if err != nil {
		return nil, nil, err
	}
	par, err := base.Paraxial(lambdaNm, p.ObjectZMM)
	if err != nil {
		return nil, nil, err
	}
	img := par.ImageZ
	if math.IsInf(img, 0) || math.IsNaN(img) {
		return nil, nil, fmt.Errorf("lens: no finite paraxial image (z=%v)", img)
	}
	return base.WithImagePlane(img), par, nil
}

// SpotResiduals samples the entrance pupil over the given fields, traces
// every ray to the paraxial image plane and returns the transverse
// displacements from each field's chief ray, divided by targetMM (so every
// residual is "spot displacement in units of the target"). Rays that fail
// to reach the image plane contribute failPenalty in both components.
//
// The residual order is, per field: 2 components per sample, in the order
// of the samples slice. The chief ray of each field is the sample (0,0).
func SpotResiduals(p *Prescription, lambdaNm float64, fields []Field, samples []Sample, targetMM float64) ([]float64, error) {
	if !(targetMM > 0) {
		return nil, fmt.Errorf("lens: spot target must be positive")
	}
	if len(samples) == 0 {
		return nil, fmt.Errorf("lens: no pupil samples")
	}
	sys, par, err := tracedAtImagePlane(p, lambdaNm)
	if err != nil {
		return nil, err
	}
	out := make([]float64, 0, len(fields)*len(samples)*2)
	for _, f := range fields {
		cr, err := launchRay(par, p.ObjectZMM, f, Sample{0, 0})
		if err != nil {
			return nil, err
		}
		cr.WavelengthNm = lambdaNm
		cres := sys.Trace(cr)
		haveChief := cres.OK()
		cx, cy := cres.Final().X, cres.Final().Y
		for _, s := range samples {
			r, err := launchRay(par, p.ObjectZMM, f, s)
			if err != nil {
				return nil, err
			}
			r.WavelengthNm = lambdaNm
			res := sys.Trace(r)
			if !haveChief || !res.OK() {
				out = append(out, failPenalty, failPenalty)
				continue
			}
			pt := res.Final()
			out = append(out, (pt.X-cx)/targetMM, (pt.Y-cy)/targetMM)
		}
	}
	return out, nil
}

// SpotRMS returns the RMS and geometric radius (mm) of one field's spot at
// wavelength lambdaNm together with the number of rays that reached the
// image plane (centroid-subtracted, chief ray included in the sample set).
func SpotRMS(p *Prescription, lambdaNm float64, field Field, samples []Sample) (rms, geoR float64, nOK int, err error) {
	sys, par, err := tracedAtImagePlane(p, lambdaNm)
	if err != nil {
		return 0, 0, 0, err
	}
	xs := make([]float64, 0, len(samples))
	ys := make([]float64, 0, len(samples))
	for _, s := range samples {
		r, err := launchRay(par, p.ObjectZMM, field, s)
		if err != nil {
			return 0, 0, 0, err
		}
		r.WavelengthNm = lambdaNm
		res := sys.Trace(r)
		if !res.OK() {
			continue
		}
		pt := res.Final()
		xs = append(xs, pt.X)
		ys = append(ys, pt.Y)
	}
	n := len(xs)
	if n == 0 {
		return 0, 0, 0, fmt.Errorf("lens: no ray reached the image plane")
	}
	var mx, my float64
	for i := 0; i < n; i++ {
		mx += xs[i]
		my += ys[i]
	}
	mx /= float64(n)
	my /= float64(n)
	var sum, worst float64
	for i := 0; i < n; i++ {
		dx, dy := xs[i]-mx, ys[i]-my
		d2 := dx*dx + dy*dy
		sum += d2
		if d2 > worst {
			worst = d2
		}
	}
	return math.Sqrt(sum / float64(n)), math.Sqrt(worst), n, nil
}

// EFLResidual returns (EFL(lambda) - target) / tolMM.
func EFLResidual(p *Prescription, lambdaNm, targetEFL, tolMM float64) (float64, error) {
	if !(tolMM > 0) {
		return 0, fmt.Errorf("lens: focal tolerance must be positive")
	}
	par, err := p.Paraxial(lambdaNm)
	if err != nil {
		return 0, err
	}
	if math.IsInf(par.EFL, 0) || math.IsNaN(par.EFL) {
		return 0, fmt.Errorf("lens: degenerate EFL %v", par.EFL)
	}
	return (par.EFL - targetEFL) / tolMM, nil
}

// ChromaticResidual returns (EFL(F) - EFL(C)) / tolMM, i.e. the axial
// chromatic focal shift scaled by its tolerance.
func ChromaticResidual(p *Prescription, tolMM float64) (float64, error) {
	if !(tolMM > 0) {
		return 0, fmt.Errorf("lens: chromatic tolerance must be positive")
	}
	parF, err := p.Paraxial(geo.WavelengthF)
	if err != nil {
		return 0, err
	}
	parC, err := p.Paraxial(geo.WavelengthC)
	if err != nil {
		return 0, err
	}
	if math.IsInf(parF.EFL, 0) || math.IsInf(parC.EFL, 0) || math.IsNaN(parF.EFL) || math.IsNaN(parC.EFL) {
		return 0, fmt.Errorf("lens: degenerate EFL at F/C")
	}
	return (parF.EFL - parC.EFL) / tolMM, nil
}
