package lens

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sort"
	"sync"

	"tetsuhiro-oz/geo"
)

// ToleranceSpec describes the manufacturing tolerances of a prescription.
// Every value is the ± half-range of a uniform perturbation: radii,
// thicknesses, conic constants and aspheric coefficients drift per surface,
// while IndexDelta is the glass melt error Δn shared by all glass cards of
// one Monte-Carlo run (one melt per batch — physically the right model and
// a determinism requirement, since the medium hook must render the same
// index for every evaluation of that run).
type ToleranceSpec struct {
	RadiusMM    float64 `json:"radius_mm"`
	ThicknessMM float64 `json:"thickness_mm"`
	ConicK      float64 `json:"conic_k"`
	CoefAbs     float64 `json:"coef_abs"`
	IndexDelta  float64 `json:"index_delta"`
}

// empty reports whether the spec perturbs nothing.
func (s ToleranceSpec) empty() bool {
	return s.RadiusMM == 0 && s.ThicknessMM == 0 && s.ConicK == 0 && s.CoefAbs == 0 && s.IndexDelta == 0
}

// ToleranceCriteria is the per-run pass/fail definition; at least one
// field must be positive. Zero fields are ignored.
type ToleranceCriteria struct {
	TargetEFLMM     float64 `json:"target_efl_mm"`
	MaxFocalErrMM   float64 `json:"max_focal_err_mm"`
	MaxColorShiftMM float64 `json:"max_color_shift_mm"`
	MaxRMSOnAxisMM  float64 `json:"max_rms_on_axis_mm"`
	MaxRMSOffAxisMM float64 `json:"max_rms_off_axis_mm"`
}

func (c ToleranceCriteria) any() bool {
	return c.MaxFocalErrMM > 0 || c.MaxColorShiftMM > 0 || c.MaxRMSOnAxisMM > 0 || c.MaxRMSOffAxisMM > 0
}

// Distribution summarizes one metric over all runs (sorted-percentile
// statistics; P90/P99 are the usual "how bad can it get" numbers).
type Distribution struct {
	Min    float64 `json:"min"`
	Max    float64 `json:"max"`
	Mean   float64 `json:"mean"`
	StdDev float64 `json:"stddev"`
	P50    float64 `json:"p50"`
	P90    float64 `json:"p90"`
	P99    float64 `json:"p99"`
}

// Nominal holds the unperturbed metrics for reference.
type Nominal struct {
	OnAxisRMSMM  float64 `json:"on_axis_rms_mm"`
	OffAxisRMSMM float64 `json:"off_axis_rms_mm"`
	ColorShiftMM float64 `json:"color_shift_mm"`
	FocalErrMM   float64 `json:"focal_err_mm"`
}

// ToleranceReport is the outcome of a Monte-Carlo tolerance analysis.
type ToleranceReport struct {
	Runs         int            `json:"runs"`
	Seed         int64          `json:"seed"`
	Yield        float64        `json:"yield"` // fraction meeting every criterion
	Nominal      Nominal        `json:"nominal"`
	OnAxisRMS    Distribution   `json:"on_axis_rms"`
	OffAxisRMS   Distribution   `json:"off_axis_rms"`
	ColorShift   Distribution   `json:"color_shift"`
	FocalErr     Distribution   `json:"focal_err"`
	FailReasons  map[string]int `json:"fail_reasons"`
	FailedTraces int            `json:"failed_traces"`
}

// ToleranceOptions tunes the run.
type ToleranceOptions struct {
	Runs     int      // default 2000
	Seed     int64    // default 1 (deterministic; change to draw another sample)
	Workers  int      // default min(8, GOMAXPROCS)
	Fields   []Field  // default on-axis + 0.5° (same as the design drivers)
	Samples  []Sample // default UnitDiskSamples(2, 6) (19 rays per field)
	LambdaNm float64  // default geo.Wavelengthd
}

func (o *ToleranceOptions) applyDefaults() {
	if o.Runs <= 0 {
		o.Runs = 2000
	}
	if o.Seed == 0 {
		o.Seed = 1
	}
	if o.Workers <= 0 {
		o.Workers = runtime.GOMAXPROCS(0)
		if o.Workers > 8 {
			o.Workers = 8
		}
	}
	if o.Fields == nil {
		o.Fields = []Field{{0, 0}, {0, 0.0087}}
	}
	if o.Samples == nil {
		o.Samples = UnitDiskSamples(2, 6)
	}
	if o.LambdaNm <= 0 {
		o.LambdaNm = geo.Wavelengthd
	}
}

// meltMedium adds a constant Δn to a base medium (glass melt error).
type meltMedium struct {
	base geo.Medium
	dn   float64
}

// N implements geo.Medium.
func (m meltMedium) N(wavelengthNm float64) float64 { return m.base.N(wavelengthNm) + m.dn }

// Name implements geo.Medium.
func (m meltMedium) Name() string { return fmt.Sprintf("%s%+.5f", m.base.Name(), m.dn) }

// toleranceRunResult is one perturbed evaluation, tagged for aggregation.
type toleranceRunResult struct {
	failed  bool
	focal   float64 // |EFL(d) − target|, NaN when unavailable
	color   float64 // |EFL(F) − EFL(C)|
	onAxis  float64
	offAxis float64
	okFocal bool
	okColor bool
	okOn    bool
	okOff   bool
}

// RunTolerances performs a deterministic Monte-Carlo tolerance analysis of
// the prescription: Runs independently perturbed prescriptions are traced
// and scored against the criteria; the report gives per-metric
// distributions and the yield (fraction meeting every criterion).
//
// Determinism: every run seeds its own PRNG with Seed+run, and results are
// aggregated by run index, so the report is bit-identical regardless of
// worker scheduling.
func RunTolerances(p *Prescription, spec ToleranceSpec, crit ToleranceCriteria, opts *ToleranceOptions) (*ToleranceReport, error) {
	if err := p.Validate(); err != nil {
		return nil, err
	}
	if spec.empty() {
		return nil, fmt.Errorf("lens: tolerance spec perturbs nothing")
	}
	if !crit.any() {
		return nil, fmt.Errorf("lens: at least one tolerance criterion must be positive")
	}
	o := ToleranceOptions{}
	if opts != nil {
		o = *opts
	}
	o.applyDefaults()

	// Nominal reference.
	nom, err := evaluate(p, o.LambdaNm, crit.TargetEFLMM, o.Fields, o.Samples)
	if err != nil {
		return nil, fmt.Errorf("lens: nominal evaluation failed: %w", err)
	}

	results := make([]toleranceRunResult, o.Runs)
	jobs := make(chan int)
	var wg sync.WaitGroup
	for w := 0; w < o.Workers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				rng := rand.New(rand.NewSource(o.Seed + int64(i)))
				q, err := perturbPrescription(p, spec, rng)
				if err != nil {
					results[i] = toleranceRunResult{failed: true}
					continue
				}
				ev, err := evaluate(q, o.LambdaNm, crit.TargetEFLMM, o.Fields, o.Samples)
				if err != nil {
					results[i] = toleranceRunResult{failed: true}
					continue
				}
				results[i] = score(ev, crit)
			}
		}()
	}
	for i := 0; i < o.Runs; i++ {
		jobs <- i
	}
	close(jobs)
	wg.Wait()

	rep := &ToleranceReport{
		Runs:        o.Runs,
		Seed:        o.Seed,
		Nominal:     nom,
		FailReasons: map[string]int{},
	}
	focal := make([]float64, 0, o.Runs)
	color := make([]float64, 0, o.Runs)
	on := make([]float64, 0, o.Runs)
	off := make([]float64, 0, o.Runs)
	pass := 0
	for _, r := range results {
		if r.failed {
			rep.FailedTraces++
			rep.FailReasons["trace"]++
			continue
		}
		focal = append(focal, r.focal)
		color = append(color, r.color)
		on = append(on, r.onAxis)
		off = append(off, r.offAxis)
		if r.okFocal && r.okColor && r.okOn && r.okOff {
			pass++
			continue
		}
		if !r.okFocal {
			rep.FailReasons["focal"]++
		}
		if !r.okColor {
			rep.FailReasons["color"]++
		}
		if !r.okOn {
			rep.FailReasons["spot_on_axis"]++
		}
		if !r.okOff {
			rep.FailReasons["spot_off_axis"]++
		}
	}
	rep.Yield = float64(pass) / float64(o.Runs)
	rep.FocalErr = distribution(focal)
	rep.ColorShift = distribution(color)
	rep.OnAxisRMS = distribution(on)
	rep.OffAxisRMS = distribution(off)
	return rep, nil
}

// evaluate computes the metrics of one prescription.
func evaluate(p *Prescription, lambdaNm, targetEFL float64, fields []Field, samples []Sample) (Nominal, error) {
	var out Nominal
	par, err := p.Paraxial(lambdaNm)
	if err != nil {
		return out, err
	}
	if targetEFL != 0 && !math.IsInf(par.EFL, 0) {
		out.FocalErrMM = math.Abs(par.EFL - targetEFL)
	}
	chrom, err := ChromaticResidual(p, 1.0)
	if err != nil {
		return out, err
	}
	out.ColorShiftMM = math.Abs(chrom)
	on, _, nOn, err := SpotRMS(p, lambdaNm, fields[0], samples)
	if err != nil || nOn == 0 {
		return out, fmt.Errorf("on-axis spot unavailable: %v", err)
	}
	out.OnAxisRMSMM = on
	if len(fields) > 1 {
		off, _, nOff, err := SpotRMS(p, lambdaNm, fields[1], samples)
		if err != nil || nOff == 0 {
			return out, fmt.Errorf("off-axis spot unavailable: %v", err)
		}
		out.OffAxisRMSMM = off
	}
	return out, nil
}

// score applies the criteria to one evaluation.
func score(ev Nominal, c ToleranceCriteria) toleranceRunResult {
	r := toleranceRunResult{
		focal:   ev.FocalErrMM,
		color:   ev.ColorShiftMM,
		onAxis:  ev.OnAxisRMSMM,
		offAxis: ev.OffAxisRMSMM,
	}
	r.okFocal = c.MaxFocalErrMM <= 0 || ev.FocalErrMM <= c.MaxFocalErrMM
	r.okColor = c.MaxColorShiftMM <= 0 || ev.ColorShiftMM <= c.MaxColorShiftMM
	r.okOn = c.MaxRMSOnAxisMM <= 0 || ev.OnAxisRMSMM <= c.MaxRMSOnAxisMM
	r.okOff = c.MaxRMSOffAxisMM <= 0 || ev.OffAxisRMSMM <= c.MaxRMSOffAxisMM
	// A criterion with target 0 and threshold set still needs the metric.
	if c.MaxFocalErrMM > 0 && c.TargetEFLMM == 0 {
		r.okFocal = true // no target: cannot judge, treat as pass
	}
	return r
}

// perturbPrescription clones p and applies one perturbation draw. The glass
// melt Δn is drawn once and installed through the mediumWrap hook so every
// evaluation of this run sees the same index.
func perturbPrescription(p *Prescription, spec ToleranceSpec, rng *rand.Rand) (*Prescription, error) {
	q := p.Clone()
	u := func(scale float64) float64 { return (rng.Float64()*2 - 1) * scale }
	for i := range q.Cards {
		c := &q.Cards[i]
		if spec.RadiusMM > 0 && c.RadiusMM != 0 && !math.IsInf(c.RadiusMM, 0) {
			c.RadiusMM += u(spec.RadiusMM)
			if math.Abs(c.RadiusMM) < 1e-6 {
				return nil, fmt.Errorf("lens: perturbed radius collapsed")
			}
		}
		if spec.ThicknessMM > 0 && i < len(q.Cards)-1 {
			c.ThicknessMM += u(spec.ThicknessMM)
			if c.ThicknessMM < 0.05 {
				c.ThicknessMM = 0.05
			}
		}
		if spec.ConicK > 0 {
			c.ConicK += u(spec.ConicK)
		}
		if spec.CoefAbs > 0 {
			for j := range c.Coeffs {
				c.Coeffs[j] += u(spec.CoefAbs)
			}
		}
	}
	if spec.IndexDelta > 0 {
		dn := u(spec.IndexDelta)
		q.mediumWrap = func(_ int, m geo.Medium) geo.Medium {
			if _, isAir := m.(geo.Air); isAir {
				return m
			}
			return meltMedium{base: m, dn: dn}
		}
	}
	return q, nil
}

// distribution computes the summary statistics (input is copied for
// sorting, so the caller's slice stays intact).
func distribution(xs []float64) Distribution {
	var d Distribution
	if len(xs) == 0 {
		return d
	}
	n := float64(len(xs))
	var sum float64
	for _, x := range xs {
		sum += x
	}
	mean := sum / n
	var sq float64
	for _, x := range xs {
		dv := x - mean
		sq += dv * dv
	}
	sorted := append([]float64(nil), xs...)
	sort.Float64s(sorted)
	return Distribution{
		Min:    sorted[0],
		Max:    sorted[len(sorted)-1],
		Mean:   mean,
		StdDev: math.Sqrt(sq / n),
		P50:    percentile(sorted, 0.50),
		P90:    percentile(sorted, 0.90),
		P99:    percentile(sorted, 0.99),
	}
}

// percentile is the nearest-rank percentile of an ascending slice.
func percentile(sorted []float64, p float64) float64 {
	if len(sorted) == 0 {
		return 0
	}
	rank := int(math.Ceil(p * float64(len(sorted))))
	if rank < 1 {
		rank = 1
	}
	if rank > len(sorted) {
		rank = len(sorted)
	}
	return sorted[rank-1]
}

// String renders the report compactly for CLIs.
func (r *ToleranceReport) String() string {
	return fmt.Sprintf(
		"tolerance: runs=%d yield=%.1f%% (failed traces %d)\n"+
			"  nominal: RMS_on=%.4fmm RMS_off=%.4fmm colour=%.4fmm focal=%.4fmm\n"+
			"  RMS_on  : mean=%.4f p90=%.4f p99=%.4f max=%.4f\n"+
			"  RMS_off : mean=%.4f p90=%.4f p99=%.4f max=%.4f\n"+
			"  colour  : mean=%.4f p90=%.4f p99=%.4f max=%.4f\n"+
			"  focal   : mean=%.4f p90=%.4f p99=%.4f max=%.4f\n"+
			"  fail reasons: %v",
		r.Runs, r.Yield*100, r.FailedTraces,
		r.Nominal.OnAxisRMSMM, r.Nominal.OffAxisRMSMM, r.Nominal.ColorShiftMM, r.Nominal.FocalErrMM,
		r.OnAxisRMS.Mean, r.OnAxisRMS.P90, r.OnAxisRMS.P99, r.OnAxisRMS.Max,
		r.OffAxisRMS.Mean, r.OffAxisRMS.P90, r.OffAxisRMS.P99, r.OffAxisRMS.Max,
		r.ColorShift.Mean, r.ColorShift.P90, r.ColorShift.P99, r.ColorShift.Max,
		r.FocalErr.Mean, r.FocalErr.P90, r.FocalErr.P99, r.FocalErr.Max,
		r.FailReasons)
}
