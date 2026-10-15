package lens

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
)

// ---------------------------------------------------------------------------
// JSON I/O of a group
// ---------------------------------------------------------------------------
//
// The wire form of a group is the struct itself: every element field carries a
// JSON tag (radius_mm, conic_k, coeffs, poly_terms, q_asphere, xy_polynomial,
// zernike, semi_aperture_mm, thickness_mm, glass, gap_glass, gap_mm, kind,
// decenter_*_mm, tilt_*_deg, fold), the surface shape fields are the ones of a
// Card, and the only special case is the object position, which may be the
// string "inf" (JSON has no infinity), exactly as in Prescription. A round trip
// marshal → unmarshal → marshal is byte-stable.

// groupJSON is the wire form of Group.
type groupJSON struct {
	Name      string          `json:"name,omitempty"`
	ObjectZMM json.RawMessage `json:"object_z_mm,omitempty"`
	Elements  []Element       `json:"elements"`
	Stop      int             `json:"stop"`
	OriginZ   float64         `json:"origin_z_mm,omitempty"`
}

// MarshalJSON implements json.Marshaler. An object at infinity (ObjectZMM =
// math.Inf(-1)) is written as the string "inf".
func (g *Group) MarshalJSON() ([]byte, error) {
	w := groupJSON{Name: g.Name, Elements: g.Elements, Stop: g.Stop, OriginZ: g.OriginZ}
	switch {
	case math.IsInf(g.ObjectZMM, -1):
		w.ObjectZMM = json.RawMessage(`"inf"`)
	case math.IsNaN(g.ObjectZMM) || math.IsInf(g.ObjectZMM, 1):
		return nil, fmt.Errorf("lens: non-finite object position %v", g.ObjectZMM)
	default:
		buf, err := json.Marshal(g.ObjectZMM)
		if err != nil {
			return nil, err
		}
		w.ObjectZMM = buf
	}
	return json.Marshal(w)
}

// UnmarshalJSON implements json.Unmarshaler; a missing object position means an
// object at infinity.
func (g *Group) UnmarshalJSON(data []byte) error {
	var w groupJSON
	if err := json.Unmarshal(data, &w); err != nil {
		return err
	}
	g.Name = w.Name
	g.Elements = w.Elements
	g.Stop = w.Stop
	g.OriginZ = w.OriginZ
	g.ObjectZMM = math.Inf(-1)
	if len(w.ObjectZMM) > 0 && string(w.ObjectZMM) != "null" {
		var f float64
		if err := json.Unmarshal(w.ObjectZMM, &f); err == nil {
			g.ObjectZMM = f
		} else {
			var s string
			if err := json.Unmarshal(w.ObjectZMM, &s); err != nil {
				return fmt.Errorf("lens: object_z_mm must be a number or \"inf\"")
			}
			switch s {
			case "inf", "-inf":
				g.ObjectZMM = math.Inf(-1)
			default:
				return fmt.Errorf("lens: unrecognised object_z_mm %q (use a number or \"inf\")", s)
			}
		}
	}
	return nil
}

// SaveJSON writes the group as pretty JSON.
func (g *Group) SaveJSON(path string) error {
	buf, err := json.MarshalIndent(g, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, append(buf, '\n'), 0o644)
}

// LoadGroupJSON reads a group from a JSON file and validates it.
func LoadGroupJSON(path string) (*Group, error) {
	buf, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var g Group
	if err := json.Unmarshal(buf, &g); err != nil {
		return nil, fmt.Errorf("lens: %s: %w", path, err)
	}
	if err := g.Validate(); err != nil {
		return nil, err
	}
	return &g, nil
}
