package lens

import (
	"fmt"

	"tetsuhiro-oz/coating"
	"tetsuhiro-oz/geo"
)

// TransmittanceSpec selects the operating points of Prescription.Transmittance.
// Zero fields take their defaults: Wavelengths {geo.Wavelengthd}; Fields
// on-axis and 0.5° (the design-driver fields); Samples UnitDiskSamples(3, 6);
// Pols s, p and unpolarized. Coatings[i] is the coating of surface i; missing
// or empty entries are bare interfaces.
type TransmittanceSpec struct {
	Wavelengths []float64
	Fields      []Field
	Samples     []Sample
	Pols        []coating.Polarization
	Coatings    []Coating
}

// TransmittanceTable is the system transmittance on the requested grid. The
// index order is [wavelength][field][polarization] (Pols order), with the
// sampled rays of a field traced at each wavelength.
type TransmittanceTable struct {
	Wavelengths []float64
	Fields      []Field
	Pols        []coating.Polarization
	// T is the mean single-pass transmittance (coating and Fresnel losses) of
	// the rays of the field that traverse the system. A field with no
	// passing ray has T = 0.
	T [][][]float64
	// Pass is the fraction of the sampled rays of the field that traverse the
	// system (vignetting and total internal reflection are the losses here).
	Pass [][]float64
	// Throughput = T × Pass: the transmitted fraction of the field's sampled
	// rays including vignetting.
	Throughput [][][]float64
}

// withDefaults returns a copy with every zero field replaced by its default.
func (s TransmittanceSpec) withDefaults() TransmittanceSpec {
	if s.Wavelengths == nil {
		s.Wavelengths = []float64{geo.Wavelengthd}
	}
	if s.Fields == nil {
		s.Fields = []Field{{0, 0}, {0, 0.0087}}
	}
	if s.Samples == nil {
		s.Samples = UnitDiskSamples(3, 6)
	}
	if s.Pols == nil {
		s.Pols = []coating.Polarization{coating.PolS, coating.PolP, coating.PolUnpolarized}
	}
	return s
}

// Transmittance evaluates the system transmittance of p with the coatings of
// spec: at every wavelength the sampled rays of every field are traced through
// the geometry, each interface is evaluated with the thin-film transfer matrix
// at the incidence angle of its ray and the requested polarization, and the
// single-pass interface transmittances are multiplied along each ray (see
// RayTransmission for the convention). The table averages over the rays that
// traverse the system.
func (p *Prescription) Transmittance(spec TransmittanceSpec) (*TransmittanceTable, error) {
	s := spec.withDefaults()
	sys, err := p.System()
	if err != nil {
		return nil, err
	}
	if len(s.Coatings) > len(sys.Surfaces) {
		return nil, fmt.Errorf("lens: %d coatings for a system of %d surfaces", len(s.Coatings), len(sys.Surfaces))
	}
	nW, nF, nP := len(s.Wavelengths), len(s.Fields), len(s.Pols)
	tab := &TransmittanceTable{
		Wavelengths: append([]float64(nil), s.Wavelengths...),
		Fields:      append([]Field(nil), s.Fields...),
		Pols:        append([]coating.Polarization(nil), s.Pols...),
		T:           make([][][]float64, nW),
		Pass:        make([][]float64, nW),
		Throughput:  make([][][]float64, nW),
	}
	for w, lam := range s.Wavelengths {
		par, err := sys.Paraxial(lam, p.ObjectZMM)
		if err != nil {
			return nil, err
		}
		tab.T[w] = make([][]float64, nF)
		tab.Pass[w] = make([]float64, nF)
		tab.Throughput[w] = make([][]float64, nF)
		for f, fld := range s.Fields {
			tab.T[w][f] = make([]float64, nP)
			tab.Throughput[w][f] = make([]float64, nP)
			sum := make([]float64, nP)
			nPass := 0
			for _, smp := range s.Samples {
				r, err := launchRay(par, p.ObjectZMM, fld, smp)
				if err != nil {
					return nil, err
				}
				r.WavelengthNm = lam
				rt, err := TraceCoated(sys, r, s.Coatings)
				if err != nil {
					return nil, err
				}
				if !rt.Passed {
					continue
				}
				nPass++
				for k, pol := range s.Pols {
					sum[k] += rt.T(pol)
				}
			}
			tab.Pass[w][f] = float64(nPass) / float64(len(s.Samples))
			for k := range s.Pols {
				if nPass > 0 {
					tab.T[w][f][k] = sum[k] / float64(nPass)
				}
				tab.Throughput[w][f][k] = tab.T[w][f][k] * tab.Pass[w][f]
			}
		}
	}
	return tab, nil
}
