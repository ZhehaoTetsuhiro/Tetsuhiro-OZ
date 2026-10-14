package coating

import (
	"encoding/json"
	"fmt"
	"os"
)

// This file adds the JSON persistence of stacks. Materials are referenced by
// their library names (LookupMaterial), so a saved design is readable,
// diffable and stable across versions; custom (non-library) dispersion
// models cannot be represented by name and are rejected on save.

// layerJSON is the wire form of one layer.
type layerJSON struct {
	Mat string  `json:"mat"`
	D   float64 `json:"d_nm"`
}

// stackJSON is the wire form of a Stack.
type stackJSON struct {
	Incident  string      `json:"incident"`
	Substrate string      `json:"substrate"`
	Layers    []layerJSON `json:"layers"`
}

// layerMaterials extracts the layer materials of a stack.
func layerMaterials(s Stack) []Material {
	out := make([]Material, len(s.Layers))
	for i, l := range s.Layers {
		out[i] = l.Material
	}
	return out
}

// MarshalJSON implements json.Marshaler. Every material must be a library
// material; a stack carrying a custom dispersion model fails with an error
// naming the material instead of silently losing information.
func (s Stack) MarshalJSON() ([]byte, error) {
	for _, m := range append([]Material{s.Incident, s.Substrate}, layerMaterials(s)...) {
		if !m.Valid() {
			return nil, fmt.Errorf("coating: cannot serialise stack: invalid material")
		}
		if _, ok := LookupMaterial(m.Name); !ok {
			return nil, fmt.Errorf("coating: cannot serialise stack: material %q is not a library material", m.Name)
		}
	}
	w := stackJSON{Incident: s.Incident.Name, Substrate: s.Substrate.Name}
	for _, l := range s.Layers {
		w.Layers = append(w.Layers, layerJSON{Mat: l.Material.Name, D: l.ThicknessNm})
	}
	return json.Marshal(w)
}

// UnmarshalJSON implements json.Unmarshaler; the reconstructed stack is
// validated (unknown material names, bad thicknesses and missing media are
// reported as errors).
func (s *Stack) UnmarshalJSON(data []byte) error {
	var w stackJSON
	if err := json.Unmarshal(data, &w); err != nil {
		return err
	}
	inc, ok := LookupMaterial(w.Incident)
	if !ok {
		return fmt.Errorf("coating: unknown incident material %q (see Library)", w.Incident)
	}
	sub, ok := LookupMaterial(w.Substrate)
	if !ok {
		return fmt.Errorf("coating: unknown substrate material %q (see Library)", w.Substrate)
	}
	layers := make([]Layer, len(w.Layers))
	for i, l := range w.Layers {
		m, ok := LookupMaterial(l.Mat)
		if !ok {
			return fmt.Errorf("coating: unknown layer %d material %q (see Library)", i, l.Mat)
		}
		layers[i] = Layer{Material: m, ThicknessNm: l.D}
	}
	out := NewStack(inc, sub, layers...)
	if err := out.Validate(); err != nil {
		return err
	}
	*s = out
	return nil
}

// SaveJSON writes the stack as indented JSON.
func (s Stack) SaveJSON(path string) error {
	buf, err := json.MarshalIndent(s, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, append(buf, '\n'), 0o644)
}

// LoadStackJSON reads and validates a stack JSON file.
func LoadStackJSON(path string) (Stack, error) {
	buf, err := os.ReadFile(path)
	if err != nil {
		return Stack{}, err
	}
	var s Stack
	if err := json.Unmarshal(buf, &s); err != nil {
		return Stack{}, fmt.Errorf("coating: %s: %w", path, err)
	}
	return s, nil
}
