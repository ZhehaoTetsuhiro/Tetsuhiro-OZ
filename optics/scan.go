package optics

import (
	"encoding/json"
	"fmt"
	"strconv"
	"strings"
)

// ScanSpec describes one parameter sweep over a scene configuration: every
// listed value is substituted into the config, the run is repeated, and one
// scalar metric is collected per value.
type ScanSpec struct {
	// Param is a dotted path into the config: "scene.components[0].params.f",
	// "scene.components[0].roll", "sources[0].params.tilt_x", "wavelength", ...
	Param string `json:"param"`
	// Values are the substituted values (numbers, in the config's own units).
	Values []float64 `json:"values"`
	// Metric selects the reported quantity: power | peak | centroid_x |
	// centroid_y | rms_x | rms_y | strehl | intensity_max.
	Metric string `json:"metric"`
	// Plane is the recorded plane to read (id or label); "" = the first one.
	Plane string `json:"plane,omitempty"`
}

// ScanResult is one point of a sweep.
type ScanResult struct {
	Value  float64 `json:"value"`
	Metric float64 `json:"metric"`
	Plane  string  `json:"plane"`
	// Error is set when this value could not be simulated (invalid config);
	// the sweep continues with the remaining values.
	Error string `json:"error,omitempty"`
}

// ScanReport is the full sweep: one entry per requested value, in order.
type ScanReport struct {
	Param  string       `json:"param"`
	Metric string       `json:"metric"`
	Plane  string       `json:"plane"`
	Points []ScanResult `json:"points"`
	// Warnings collects every run warning, de-duplicated.
	Warnings []string `json:"warnings,omitempty"`
}

// RunScan evaluates cfg once per value of spec.Values, substituting the value
// at spec.Param, and reports spec.Metric of the selected plane. It is a pure
// function of (cfg, spec): identical input gives an identical report.
func RunScan(cfg Config, spec ScanSpec) (*ScanReport, error) {
	if spec.Param == "" {
		return nil, fmt.Errorf("scan: param is required")
	}
	if len(spec.Values) == 0 {
		return nil, fmt.Errorf("scan: values is empty")
	}
	if len(spec.Values) > maxScanPoints {
		return nil, fmt.Errorf("scan: %d values exceeds the limit of %d", len(spec.Values), maxScanPoints)
	}
	if _, err := ParseMethod(cfg.Method); err != nil {
		return nil, err
	}
	rep := &ScanReport{Param: spec.Param, Metric: spec.Metric}
	seen := map[string]bool{}
	for _, v := range spec.Values {
		c := cfg // shallow copy; deep-copied per value below
		one, err := cloneConfig(&c)
		if err != nil {
			return nil, err
		}
		if err := setConfigValue(one, spec.Param, v); err != nil {
			rep.Points = append(rep.Points, ScanResult{Value: v, Error: err.Error()})
			continue
		}
		res, err := Simulate(*one)
		if err != nil {
			rep.Points = append(rep.Points, ScanResult{Value: v, Error: err.Error()})
			continue
		}
		for _, w := range res.Warnings {
			key := w.Code + "|" + w.Message
			if !seen[key] {
				seen[key] = true
				rep.Warnings = append(rep.Warnings, key)
			}
		}
		pl := selectPlane(res, spec.Plane)
		if pl == nil {
			rep.Points = append(rep.Points, ScanResult{Value: v,
				Error: fmt.Sprintf("no plane %q in the run (have %d)", spec.Plane, len(res.Planes))})
			continue
		}
		mv, merr := planeMetric(pl, spec.Metric)
		if merr != nil {
			return nil, merr
		}
		rep.Points = append(rep.Points, ScanResult{Value: v, Metric: mv, Plane: pl.ID})
	}
	return rep, nil
}

// maxScanPoints bounds one sweep: every point is a full simulation, and the
// request has to stay bounded (and the response serialisable).
const maxScanPoints = 512

// selectPlane picks the plane a metric is read from: by id, then by label,
// then the first one.
func selectPlane(res *Result, want string) *Plane {
	if len(res.Planes) == 0 {
		return nil
	}
	if want == "" {
		return res.Planes[0]
	}
	for _, p := range res.Planes {
		if p.ID == want || p.Label == want {
			return p
		}
	}
	return nil
}

// planeMetric extracts one scalar from a plane's statistics. An unknown
// metric name is an error, not a silent zero.
func planeMetric(pl *Plane, metric string) (float64, error) {
	switch metric {
	case "", "power":
		return pl.Stats.Power, nil
	case "peak":
		return pl.Stats.Peak, nil
	case "centroid_x":
		return pl.Stats.CentroidX, nil
	case "centroid_y":
		return pl.Stats.CentroidY, nil
	case "rms_x":
		return pl.Stats.RMSX, nil
	case "rms_y":
		return pl.Stats.RMSY, nil
	case "strehl":
		return pl.Stats.Strehl, nil
	case "intensity_max":
		return pl.Stats.IntensityMax, nil
	case "intensity_min":
		return pl.Stats.IntensityMin, nil
	}
	return 0, fmt.Errorf("scan: unknown metric %q (have: power, peak, centroid_x, centroid_y, rms_x, rms_y, strehl, intensity_max, intensity_min)", metric)
}

// setConfigValue writes v at a dotted path of the configuration. Supported
// paths: top-level scalars (wavelength, method, ...), sources[i].params.k,
// scene.components[i].params.k, scene.components[i].{roll,yaw,pitch},
// scene.components[i].shape.params.k and elements[i].params.k. Anything else
// is refused, so a typo cannot silently do nothing.
func setConfigValue(cfg *Config, path string, v float64) error {
	parts := strings.Split(path, ".")
	if len(parts) == 1 {
		// A top-level scalar field (wavelength).
		switch parts[0] {
		case "wavelength":
			cfg.Wavelength = v
			return nil
		}
		return fmt.Errorf("scan: %q is not a sweepable top-level field (have: wavelength)", path)
	}
	if parts[0] == "grid" {
		// The grid is sampling, not physics: sweeping it would change the
		// window and the pixel size of every run, so the "curve" it produces
		// measures the discretisation, not the system. Refused on purpose.
		return fmt.Errorf("scan: grid.* is not sweepable (the grid is sampling, not physics)")
	}
	// The legacy single-source form ("source.params.k") is addressed by
	// indexRest as if it were "sources[0]".
	if parts[0] == "source" {
		parts = append([]string{"sources[0]"}, parts[1:]...)
	}
	switch {
	case strings.HasPrefix(parts[0], "sources["):
		i, rest, err := indexRest(parts, "sources")
		if err != nil {
			return err
		}
		if i >= len(cfg.Sources) {
			return fmt.Errorf("scan: sources[%d] out of range (%d sources)", i, len(cfg.Sources))
		}
		return setParamValue(&cfg.Sources[i].Params, rest, v)
	case parts[0] == "scene":
		if cfg.Scene == nil {
			return fmt.Errorf("scan: the configuration has no scene")
		}
		// scene.components[i].<...>: indexRest parses the "components[i]"
		// segment, so the path is checked by prefix there.
		i, rest, err := indexRest(parts[1:], "components")
		if err != nil {
			return err
		}
		if i >= len(cfg.Scene.Components) {
			return fmt.Errorf("scan: components[%d] out of range (%d components)", i, len(cfg.Scene.Components))
		}
		c := &cfg.Scene.Components[i]
		if len(rest) == 0 {
			return fmt.Errorf("scan: %q names a component, not a field", path)
		}
		switch rest[0] {
		case "params":
			return setParamValue(&c.Params, rest, v)
		case "shape":
			if len(rest) < 3 || rest[1] != "params" {
				return fmt.Errorf("scan: %q must be ...shape.params.<key>", path)
			}
			if c.Shape == nil {
				c.Shape = &ShapeSpec{Kind: "circle", Params: map[string]any{}}
			}
			if c.Shape.Params == nil {
				c.Shape.Params = map[string]any{}
			}
			c.Shape.Params[rest[2]] = v
			return nil
		case "roll", "yaw", "pitch":
			if len(rest) != 1 {
				return fmt.Errorf("scan: %q is a scalar field", path)
			}
			switch rest[0] {
			case "roll":
				c.Roll = v
			case "yaw":
				c.Yaw = v
			case "pitch":
				c.Pitch = v
			}
			return nil
		case "pos":
			if len(rest) != 2 {
				return fmt.Errorf("scan: %q must be ...pos.{x,y,z}", path)
			}
			switch rest[1] {
			case "x":
				c.Pos.X = v
			case "y":
				c.Pos.Y = v
			case "z":
				c.Pos.Z = v
			default:
				return fmt.Errorf("scan: pos.%s is not a coordinate", rest[1])
			}
			return nil
		}
		return fmt.Errorf("scan: unsupported component path %q", path)
	case strings.HasPrefix(parts[0], "elements["):
		i, rest, err := indexRest(parts, "elements")
		if err != nil {
			return err
		}
		if i >= len(cfg.Elements) {
			return fmt.Errorf("scan: elements[%d] out of range (%d elements)", i, len(cfg.Elements))
		}
		return setParamValue(&cfg.Elements[i].Params, rest, v)
	}
	return fmt.Errorf("scan: unsupported path %q", path)
}

// indexRest parses the segment "name[i]" and returns i plus the remaining
// segments (the ones after the index).
func indexRest(parts []string, name string) (int, []string, error) {
	seg := parts[0]
	if !strings.HasPrefix(seg, name+"[") || !strings.HasSuffix(seg, "]") {
		return 0, nil, fmt.Errorf("scan: expected %s[i], got %q", name, seg)
	}
	i, err := strconv.Atoi(strings.TrimSuffix(strings.TrimPrefix(seg, name+"["), "]"))
	if err != nil || i < 0 {
		return 0, nil, fmt.Errorf("scan: %q is not a valid index", seg)
	}
	return i, parts[1:], nil
}

// setParamValue writes v into a parameter map at the remaining path
// (params.<key> or params.<key>.<sub>).
func setParamValue(m *map[string]any, rest []string, v float64) error {
	if len(rest) < 2 || rest[0] != "params" {
		return fmt.Errorf("scan: expected params.<key>, got %v", rest)
	}
	if *m == nil {
		*m = map[string]any{}
	}
	if len(rest) == 2 {
		(*m)[rest[1]] = v
		return nil
	}
	// Nested object: params.<key>.<sub>.
	sub, ok := (*m)[rest[1]].(map[string]any)
	if !ok {
		sub = map[string]any{}
	}
	cur := sub
	for _, k := range rest[2 : len(rest)-1] {
		next, ok := cur[k].(map[string]any)
		if !ok {
			next = map[string]any{}
		}
		cur[k] = next
		cur = next
	}
	cur[rest[len(rest)-1]] = v
	(*m)[rest[1]] = sub
	return nil
}

// cloneConfig deep-copies a configuration through JSON, so a sweep can mutate
// one value per run without touching the caller's object.
func cloneConfig(cfg *Config) (*Config, error) {
	buf, err := json.Marshal(cfg)
	if err != nil {
		return nil, fmt.Errorf("scan: encode config: %w", err)
	}
	var out Config
	if err := json.Unmarshal(buf, &out); err != nil {
		return nil, fmt.Errorf("scan: decode config: %w", err)
	}
	return &out, nil
}
