package plot

import (
	"math"
	"os"
	"strconv"
	"strings"
)

// SeriesKind represents the type of data series.
type SeriesKind int

const (
	// Line represents a line chart series.
	Line SeriesKind = iota
	// Scatter represents a scatter plot series.
	Scatter
)

// series represents a single data series in the chart.
type series struct {
	name  string
	color string
	kind  SeriesKind
	pts   [][2]float64
}

// Chart represents a 2D plot chart.
type Chart struct {
	Width, Height float64
	title         string
	xLabel        string
	yLabel        string
	xmin, xmax    float64
	ymin, ymax    float64
	xAuto, yAuto  bool
	series        []series
}

// NewChart creates a new Chart with the specified width and height.
func NewChart(width, height float64) *Chart {
	return &Chart{
		Width:  width,
		Height: height,
		xAuto:  true,
		yAuto:  true,
	}
}

// SetTitle sets the title of the chart.
func (c *Chart) SetTitle(title string) {
	c.title = title
}

// SetLabels sets the x and y axis labels.
func (c *Chart) SetLabels(xLabel, yLabel string) {
	c.xLabel = xLabel
	c.yLabel = yLabel
}

// SetRange sets the range for the x and y axes.
// Passing 0 for both min and max of an axis enables auto-scaling for that axis.
func (c *Chart) SetRange(xmin, xmax, ymin, ymax float64) {
	if xmin == 0 && xmax == 0 {
		c.xAuto = true
	} else {
		c.xAuto = false
		c.xmin = xmin
		c.xmax = xmax
	}
	if ymin == 0 && ymax == 0 {
		c.yAuto = true
	} else {
		c.yAuto = false
		c.ymin = ymin
		c.ymax = ymax
	}
}

// AddLine adds a line series to the chart.
func (c *Chart) AddLine(name, color string, pts [][2]float64) {
	c.series = append(c.series, series{
		name:  name,
		color: color,
		kind:  Line,
		pts:   pts,
	})
}

// AddScatter adds a scatter series to the chart.
func (c *Chart) AddScatter(name, color string, pts [][2]float64) {
	c.series = append(c.series, series{
		name:  name,
		color: color,
		kind:  Scatter,
		pts:   pts,
	})
}

// escapeXML escapes special XML characters in the input string.
func escapeXML(s string) string {
	var b strings.Builder
	for _, ch := range s {
		switch ch {
		case '<':
			b.WriteString("&lt;")
		case '>':
			b.WriteString("&gt;")
		case '&':
			b.WriteString("&amp;")
		case '"':
			b.WriteString("&quot;")
		case '\'':
			b.WriteString("&apos;")
		default:
			b.WriteRune(ch)
		}
	}
	return b.String()
}

// formatFloat formats a float64 with up to 4 significant digits.
func formatFloat(f float64) string {
	return strconv.FormatFloat(f, 'g', 4, 64)
}

// niceNum calculates a "nice" number for axis ticks.
// (review fix: the parameter must not be named `range` — it is a Go keyword)
func niceNum(span float64, round bool) float64 {
	exp := int(math.Log10(span))
	f := span / math.Pow(10, float64(exp))
	var nf float64
	if round {
		if f < 1.5 {
			nf = 1
		} else if f < 3 {
			nf = 2
		} else if f < 7 {
			nf = 5
		} else {
			nf = 10
		}
	} else {
		if f <= 1 {
			nf = 1
		} else if f <= 2 {
			nf = 2
		} else if f <= 5 {
			nf = 5
		} else {
			nf = 10
		}
	}
	return nf * math.Pow(10, float64(exp))
}

// calcTicks calculates nice tick values for an axis.
func calcTicks(min, max float64) []float64 {
	if min == max {
		return []float64{min}
	}
	rangeVal := max - min
	if rangeVal < 0 {
		rangeVal = -rangeVal
	}
	if rangeVal == 0 {
		return []float64{min}
	}

	// Aim for 5-8 ticks
	tickSpacing := niceNum(rangeVal/6, true)
	if tickSpacing == 0 {
		tickSpacing = 1
	}

	niceMin := math.Floor(min/tickSpacing) * tickSpacing
	niceMax := math.Ceil(max/tickSpacing) * tickSpacing

	var ticks []float64
	for t := niceMin; t <= niceMax+tickSpacing*0.5; t += tickSpacing {
		// Avoid floating point issues
		if t >= min-tickSpacing*0.5 && t <= max+tickSpacing*0.5 {
			ticks = append(ticks, t)
		}
	}

	// Ensure we have at least one tick
	if len(ticks) == 0 {
		ticks = append(ticks, min)
	}

	return ticks
}

// SVG generates the complete SVG document as a string.
func (c *Chart) SVG() string {
	// Determine data ranges
	var dataXMin, dataXMax, dataYMin, dataYMax float64
	hasData := false

	for _, s := range c.series {
		for _, p := range s.pts {
			if !hasData {
				dataXMin = p[0]
				dataXMax = p[0]
				dataYMin = p[1]
				dataYMax = p[1]
				hasData = true
			} else {
				if p[0] < dataXMin {
					dataXMin = p[0]
				}
				if p[0] > dataXMax {
					dataXMax = p[0]
				}
				if p[1] < dataYMin {
					dataYMin = p[1]
				}
				if p[1] > dataYMax {
					dataYMax = p[1]
				}
			}
		}
	}

	// Set axis ranges
	xmin, xmax, ymin, ymax := c.xmin, c.xmax, c.ymin, c.ymax

	if c.xAuto {
		if hasData {
			xmin = dataXMin
			xmax = dataXMax
			if xmin == xmax {
				xmin -= 1
				xmax += 1
			}
		} else {
			xmin = 0
			xmax = 1
		}
	}

	if c.yAuto {
		if hasData {
			ymin = dataYMin
			ymax = dataYMax
			if ymin == ymax {
				ymin -= 1
				ymax += 1
			}
		} else {
			ymin = 0
			ymax = 1
		}
	}

	// Calculate ticks
	xTicks := calcTicks(xmin, xmax)
	yTicks := calcTicks(ymin, ymax)

	// Define plot area margins
	marginLeft := 60.0
	marginRight := 20.0
	marginTop := 40.0
	marginBottom := 50.0

	plotWidth := c.Width - marginLeft - marginRight
	plotHeight := c.Height - marginTop - marginBottom

	if plotWidth <= 0 {
		plotWidth = 1
	}
	if plotHeight <= 0 {
		plotHeight = 1
	}

	// Coordinate transformation functions (clamped to the plot area so a tick
	// that lands a hair outside the data range cannot draw outside the frame)
	clamp := func(v, lo, hi float64) float64 {
		if v < lo {
			return lo
		}
		if v > hi {
			return hi
		}
		return v
	}
	xToPixel := func(x float64) float64 {
		if xmax == xmin {
			return marginLeft + plotWidth/2
		}
		return clamp(marginLeft+(x-xmin)/(xmax-xmin)*plotWidth, marginLeft, marginLeft+plotWidth)
	}

	yToPixel := func(y float64) float64 {
		if ymax == ymin {
			return marginTop + plotHeight/2
		}
		return clamp(marginTop+plotHeight-(y-ymin)/(ymax-ymin)*plotHeight, marginTop, marginTop+plotHeight)
	}

	var b strings.Builder

	// SVG header
	b.WriteString("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 ")
	b.WriteString(formatFloat(c.Width))
	b.WriteString(" ")
	b.WriteString(formatFloat(c.Height))
	b.WriteString("\" style=\"width:100%;height:auto\">\n")

	// Background
	b.WriteString("  <rect width=\"100%\" height=\"100%\" fill=\"#1e1e1e\"/>\n")

	// Title
	if c.title != "" {
		b.WriteString("  <text x=\"")
		b.WriteString(formatFloat(c.Width / 2))
		b.WriteString("\" y=\"20\" text-anchor=\"middle\" fill=\"#ffffff\" font-size=\"14\" font-family=\"sans-serif\">")
		b.WriteString(escapeXML(c.title))
		b.WriteString("</text>\n")
	}

	// Grid lines and ticks
	// X-axis grid and ticks
	for _, tick := range xTicks {
		x := xToPixel(tick)
		// Grid line
		b.WriteString("  <line x1=\"")
		b.WriteString(formatFloat(x))
		b.WriteString("\" y1=\"")
		b.WriteString(formatFloat(marginTop))
		b.WriteString("\" x2=\"")
		b.WriteString(formatFloat(x))
		b.WriteString("\" y2=\"")
		b.WriteString(formatFloat(marginTop + plotHeight))
		b.WriteString("\" stroke=\"#444444\" stroke-width=\"1\"/>\n")
		// Tick label
		b.WriteString("  <text x=\"")
		b.WriteString(formatFloat(x))
		b.WriteString("\" y=\"")
		b.WriteString(formatFloat(marginTop + plotHeight + 15))
		b.WriteString("\" text-anchor=\"middle\" fill=\"#cccccc\" font-size=\"10\" font-family=\"sans-serif\">")
		b.WriteString(escapeXML(formatFloat(tick)))
		b.WriteString("</text>\n")
	}

	// Y-axis grid and ticks
	for _, tick := range yTicks {
		y := yToPixel(tick)
		// Grid line
		b.WriteString("  <line x1=\"")
		b.WriteString(formatFloat(marginLeft))
		b.WriteString("\" y1=\"")
		b.WriteString(formatFloat(y))
		b.WriteString("\" x2=\"")
		b.WriteString(formatFloat(marginLeft + plotWidth))
		b.WriteString("\" y2=\"")
		b.WriteString(formatFloat(y))
		b.WriteString("\" stroke=\"#444444\" stroke-width=\"1\"/>\n")
		// Tick label
		b.WriteString("  <text x=\"")
		b.WriteString(formatFloat(marginLeft - 5))
		b.WriteString("\" y=\"")
		b.WriteString(formatFloat(y + 3))
		b.WriteString("\" text-anchor=\"end\" fill=\"#cccccc\" font-size=\"10\" font-family=\"sans-serif\">")
		b.WriteString(escapeXML(formatFloat(tick)))
		b.WriteString("</text>\n")
	}

	// Axes
	// X-axis
	b.WriteString("  <line x1=\"")
	b.WriteString(formatFloat(marginLeft))
	b.WriteString("\" y1=\"")
	b.WriteString(formatFloat(marginTop + plotHeight))
	b.WriteString("\" x2=\"")
	b.WriteString(formatFloat(marginLeft + plotWidth))
	b.WriteString("\" y2=\"")
	b.WriteString(formatFloat(marginTop + plotHeight))
	b.WriteString("\" stroke=\"#ffffff\" stroke-width=\"1.5\"/>\n")

	// Y-axis
	b.WriteString("  <line x1=\"")
	b.WriteString(formatFloat(marginLeft))
	b.WriteString("\" y1=\"")
	b.WriteString(formatFloat(marginTop))
	b.WriteString("\" x2=\"")
	b.WriteString(formatFloat(marginLeft))
	b.WriteString("\" y2=\"")
	b.WriteString(formatFloat(marginTop + plotHeight))
	b.WriteString("\" stroke=\"#ffffff\" stroke-width=\"1.5\"/>\n")

	// X-axis label
	if c.xLabel != "" {
		b.WriteString("  <text x=\"")
		b.WriteString(formatFloat(marginLeft + plotWidth/2))
		b.WriteString("\" y=\"")
		b.WriteString(formatFloat(c.Height - 10))
		b.WriteString("\" text-anchor=\"middle\" fill=\"#ffffff\" font-size=\"12\" font-family=\"sans-serif\">")
		b.WriteString(escapeXML(c.xLabel))
		b.WriteString("</text>\n")
	}

	// Y-axis label
	if c.yLabel != "" {
		b.WriteString("  <text x=\"15\" y=\"")
		b.WriteString(formatFloat(marginTop + plotHeight/2))
		b.WriteString("\" text-anchor=\"middle\" fill=\"#ffffff\" font-size=\"12\" font-family=\"sans-serif\" transform=\"rotate(-90 15,")
		b.WriteString(formatFloat(marginTop + plotHeight/2))
		b.WriteString(")\">")
		b.WriteString(escapeXML(c.yLabel))
		b.WriteString("</text>\n")
	}

	// Data series
	for _, s := range c.series {
		if s.kind == Line {
			if len(s.pts) > 1 {
				var points []string
				for _, p := range s.pts {
					points = append(points, formatFloat(xToPixel(p[0]))+","+formatFloat(yToPixel(p[1])))
				}
				b.WriteString("  <polyline points=\"")
				b.WriteString(strings.Join(points, " "))
				b.WriteString("\" fill=\"none\" stroke=\"")
				b.WriteString(escapeXML(s.color))
				b.WriteString("\" stroke-width=\"2\"/>\n")
			} else if len(s.pts) == 1 {
				// Single point for line series - draw as circle
				p := s.pts[0]
				b.WriteString("  <circle cx=\"")
				b.WriteString(formatFloat(xToPixel(p[0])))
				b.WriteString("\" cy=\"")
				b.WriteString(formatFloat(yToPixel(p[1])))
				b.WriteString("\" r=\"2.5\" fill=\"")
				b.WriteString(escapeXML(s.color))
				b.WriteString("\"/>\n")
			}
		} else if s.kind == Scatter {
			for _, p := range s.pts {
				b.WriteString("  <circle cx=\"")
				b.WriteString(formatFloat(xToPixel(p[0])))
				b.WriteString("\" cy=\"")
				b.WriteString(formatFloat(yToPixel(p[1])))
				b.WriteString("\" r=\"2.5\" fill=\"")
				b.WriteString(escapeXML(s.color))
				b.WriteString("\"/>\n")
			}
		}
	}

	// Legend (kept inside the plot frame even for very short charts)
	legendX := marginLeft + plotWidth - 120
	if legendX < marginLeft {
		legendX = marginLeft
	}
	legendY := marginTop + 10
	legendRowHeight := 18.0

	for i, s := range c.series {
		y := legendY + float64(i)*legendRowHeight
		// Color box
		b.WriteString("  <rect x=\"")
		b.WriteString(formatFloat(legendX))
		b.WriteString("\" y=\"")
		b.WriteString(formatFloat(y))
		b.WriteString("\" width=\"12\" height=\"12\" fill=\"")
		b.WriteString(escapeXML(s.color))
		b.WriteString("\"/>\n")
		// Name
		b.WriteString("  <text x=\"")
		b.WriteString(formatFloat(legendX + 16))
		b.WriteString("\" y=\"")
		b.WriteString(formatFloat(y + 10))
		b.WriteString("\" fill=\"#ffffff\" font-size=\"10\" font-family=\"sans-serif\">")
		b.WriteString(escapeXML(s.name))
		b.WriteString("</text>\n")
	}

	b.WriteString("</svg>\n")

	return b.String()
}

// WriteFile writes the SVG output to the specified file path.
func (c *Chart) WriteFile(path string) error {
	content := c.SVG()
	return os.WriteFile(path, []byte(content), 0o644)
}
