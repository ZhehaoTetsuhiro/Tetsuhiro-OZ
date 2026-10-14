package coating

import (
	"fmt"
	"testing"
	"time"
)

func TestZZProbe(t *testing.T) {
	if testing.Short() {
		t.Skip("probe")
	}
	report := func(name string, fn func() (Stack, *Result, error), extra func(Stack) string) {
		t0 := time.Now()
		st, res, err := fn()
		dt := time.Since(t0)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		line := fmt.Sprintf("%-16s %5.2fs layers=%2d merit=%.4g starts=%d iters=%d", name, dt.Seconds(), len(st.Layers), res.Merit, res.Starts, res.Iterations)
		if extra != nil {
			line += " | " + extra(st)
		}
		t.Log(line)
	}
	meanR := func(st Stack, from, to float64) string {
		v, err := st.BandAverageR(from, to, 5, 0, PolUnpolarized)
		if err != nil {
			return err.Error()
		}
		return fmt.Sprintf("meanR(%.0f-%.0f)=%.5f", from, to, v)
	}

	report("SingleLayerAR", func() (Stack, *Result, error) {
		return SingleLayerAR(Air, BK7, 550, Options{})
	}, func(st Stack) string {
		s, _ := st.Evaluate(550, 0, PolUnpolarized)
		return fmt.Sprintf("mat=%s R550=%.5f", st.Layers[0].Material.Name, s.R)
	})
	report("VCoat", func() (Stack, *Result, error) { return VCoat(Air, BK7, 550, Options{}) }, func(st Stack) string {
		s, _ := st.Evaluate(550, 0, PolUnpolarized)
		s2, _ := st.Evaluate(470, 0, PolUnpolarized)
		return fmt.Sprintf("R550=%.6f R470=%.4f", s.R, s2.R)
	})
	report("BroadbandAR", func() (Stack, *Result, error) {
		return BroadbandAR(Air, BK7, 420, 680, Options{})
	}, func(st Stack) string {
		v, _ := st.BandAverageR(420, 680, 5, 0, PolUnpolarized)
		mx, _ := st.Reflectance(680, 0, PolUnpolarized)
		return fmt.Sprintf("mean=%0.5f maxR680=%.4f", v, mx)
	})
	report("HighReflector14", func() (Stack, *Result, error) {
		return HighReflector(Air, BK7, 1064, 14, Options{})
	}, func(st Stack) string {
		a, _ := st.Reflectance(1064, 0, PolUnpolarized)
		b, _ := st.Reflectance(532, 0, PolUnpolarized)
		c, _ := st.Reflectance(1000, 0, PolUnpolarized)
		return fmt.Sprintf("R1064=%.6f R532=%.4f R1000=%.6f", a, b, c)
	})
	report("LongWavePass550", func() (Stack, *Result, error) {
		return LongWavePass(Air, BK7, 550, Options{})
	}, func(st Stack) string {
		sb, _ := st.BandAverageT(440, 539, 5, 0, PolUnpolarized)
		pb, _ := st.BandAverageT(567, 688, 5, 0, PolUnpolarized)
		return fmt.Sprintf("stopT=%.5f passT=%.5f", sb, pb)
	})
	report("ShortWavePass550", func() (Stack, *Result, error) {
		return ShortWavePass(Air, BK7, 550, Options{})
	}, func(st Stack) string {
		pb, _ := st.BandAverageT(440, 534, 5, 0, PolUnpolarized)
		sb, _ := st.BandAverageT(561, 688, 5, 0, PolUnpolarized)
		return fmt.Sprintf("passT=%.5f stopT=%.5f", pb, sb)
	})
	report("Splitter50", func() (Stack, *Result, error) {
		return Splitter50(Air, BK7, 550, Options{})
	}, func(st Stack) string {
		a, _ := st.Reflectance(550, 0, PolUnpolarized)
		v, _ := st.BandAverageR(500, 600, 5, 0, PolUnpolarized)
		return fmt.Sprintf("R550=%.5f meanR=%.5f", a, v)
	})
	// Layer growth from a two-layer start.
	report("GrowToBroadband", func() (Stack, *Result, error) {
		targets, err := BandTargets(420, 680, 20, 0, PolUnpolarized, Reflectance, 1e-4, 1, true)
		if err != nil {
			return Stack{}, nil, err
		}
		start := NewStack(Air, BK7, Layer{Material: MgF2, ThicknessNm: 60}, Layer{Material: Ta2O5, ThicknessNm: 80})
		res, err := GrowLayers(start, targets, []Material{SiO2, Ta2O5, TiO2, MgF2}, Options{GrowthRounds: 4, Restarts: 2})
		if err != nil {
			return Stack{}, nil, err
		}
		return res.Stack, res, nil
	}, func(st Stack) string {
		return meanR(st, 420, 680) + fmt.Sprintf(" total=%.0fnm", st.TotalThicknessNm())
	})
}
