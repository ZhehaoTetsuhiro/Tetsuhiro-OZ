package plot

import (
	"os"
	"strings"
	"testing"
)

func TestSVGStructure(t *testing.T) {
	c := NewChart(800, 600)
	c.SetTitle("Test Chart")
	c.SetLabels("X Axis", "Y Axis")
	c.AddLine("Line1", "#ff0000", [][2]float64{{0, 0}, {1, 1}, {2, 2}})
	c.AddScatter("Scatter1", "#00ff00", [][2]float64{{0.5, 0.5}, {1.5, 1.5}})

	svg := c.SVG()

	if !strings.HasPrefix(svg, "<svg") {
		t.Errorf("SVG should start with <svg, got: %s", svg[:10])
	}

	if !strings.HasSuffix(svg, "</svg>\n") && !strings.HasSuffix(svg, "</svg>") {
		t.Errorf("SVG should end with </svg>, got: %s", svg[len(svg)-10:])
	}

	if !strings.Contains(svg, "Line1") {
		t.Error("SVG should contain series name 'Line1'")
	}

	if !strings.Contains(svg, "Scatter1") {
		t.Error("SVG should contain series name 'Scatter1'")
	}

	if !strings.Contains(svg, "#ff0000") {
		t.Error("SVG should contain color '#ff0000'")
	}

	if !strings.Contains(svg, "#00ff00") {
		t.Error("SVG should contain color '#00ff00'")
	}
}

func TestSVGXMLEscape(t *testing.T) {
	c := NewChart(800, 600)
	c.AddLine("<&>Test", "#ff0000", [][2]float64{{0, 0}, {1, 1}})

	svg := c.SVG()

	if strings.Contains(svg, "<&>") {
		t.Error("SVG should not contain unescaped '<&>'")
	}

	if !strings.Contains(svg, "&lt;&amp;&gt;Test") {
		t.Error("SVG should contain escaped '&lt;&amp;&gt;Test'")
	}
}

func TestSVGDeterministic(t *testing.T) {
	c := NewChart(800, 600)
	c.SetTitle("Deterministic Test")
	c.AddLine("Line1", "#ff0000", [][2]float64{{0, 0}, {1, 1}, {2, 2}})
	c.AddScatter("Scatter1", "#00ff00", [][2]float64{{0.5, 0.5}, {1.5, 1.5}})

	svg1 := c.SVG()
	svg2 := c.SVG()

	if svg1 != svg2 {
		t.Error("SVG output should be deterministic")
	}
}

func TestSVGEmptyChart(t *testing.T) {
	c := NewChart(800, 600)

	svg := c.SVG()

	if !strings.HasPrefix(svg, "<svg") {
		t.Errorf("Empty chart SVG should start with <svg, got: %s", svg[:10])
	}

	if !strings.HasSuffix(svg, "</svg>\n") && !strings.HasSuffix(svg, "</svg>") {
		t.Errorf("Empty chart SVG should end with </svg>, got: %s", svg[len(svg)-10:])
	}
}

func TestSVGSinglePointScatter(t *testing.T) {
	c := NewChart(800, 600)
	c.AddScatter("Single", "#ff0000", [][2]float64{{5, 5}})

	svg := c.SVG()

	if !strings.Contains(svg, "<circle") {
		t.Error("Single point scatter should produce at least one <circle>")
	}
}

func TestWriteFile(t *testing.T) {
	c := NewChart(800, 600)
	c.SetTitle("File Test")
	c.AddLine("Line1", "#ff0000", [][2]float64{{0, 0}, {1, 1}})

	expected := c.SVG()

	tmpFile := t.TempDir() + "/test.svg"

	err := c.WriteFile(tmpFile)
	if err != nil {
		t.Fatalf("WriteFile failed: %v", err)
	}

	content, err := os.ReadFile(tmpFile)
	if err != nil {
		t.Fatalf("ReadFile failed: %v", err)
	}

	if string(content) != expected {
		t.Error("File content should match SVG() output")
	}
}
