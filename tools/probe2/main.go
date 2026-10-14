package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

func idx(g string) float64 { m, _ := geo.Glass(g); return m.N(geo.Wavelengthd) }
func vv(g string) float64  { m, _ := geo.Glass(g); return m.AbbeNumber() }

func seedAchromat(crown, flint string, f float64) []float64 {
	n1, n2 := idx(crown), idx(flint)
	V1, V2 := vv(crown), vv(flint)
	phi := 1 / f
	phi1 := phi * V1 / (V1 + V2)
	phi2 := phi - phi1
	cA := phi1 / (2 * (n1 - 1))
	cB := -cA
	cC := cB - phi2/(n2-1)
	return []float64{cA, cB, cC}
}

func main() {
	// --- rev doublet
	efl, semi := 100.0, 10.0
	base, _, err := lens.DesignAchromat(efl, "N-BK7", "F2", semi, &lens.DesignOptions{SkipPolish: true})
	fmt.Println("design err", err)
	if err == nil {
		for i, c := range base.Cards {
			fmt.Printf("  base %d R=%.4f t=%.3f %s\n", i, c.RadiusMM, c.ThicknessMM, c.Glass)
		}
		par, _ := base.Paraxial(geo.Wavelengthd)
		fmt.Printf("  base EFL=%.6f BFL=%.6f\n", par.EFL, par.BFL)
		n := len(base.Cards)
		cards := []lens.Card{{Name: "stop", RadiusMM: math.Inf(1), ThicknessMM: 3, SemiAperMM: semi}}
		for j := 0; j < n; j++ {
			src := base.Cards[n-1-j]
			t := 40.0
			g := ""
			if j < n-1 {
				t = base.Cards[n-2-j].ThicknessMM
				g = base.Cards[n-2-j].Glass
			}
			cards = append(cards, lens.Card{Name: src.Name, RadiusMM: src.RadiusMM, ThicknessMM: t, Glass: g, SemiAperMM: semi * 1.3})
		}
		p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
		par2, err2 := p.Paraxial(geo.Wavelengthd)
		fmt.Printf("  reversed EFL=%v err=%v (raw)\n", par2.EFL, err2)
		// scan the scale to see the shape of EFL(s)
		orig := make([]float64, 0)
		idx2 := []int{}
		for i, c := range p.Cards {
			if !math.IsInf(c.RadiusMM, 0) && c.RadiusMM != 0 {
				orig = append(orig, c.RadiusMM)
				idx2 = append(idx2, i)
			}
		}
		for _, s := range []float64{0.5, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2} {
			for k, i := range idx2 {
				p.Cards[i].RadiusMM = orig[k] / s
			}
			pp, _ := p.Paraxial(geo.Wavelengthd)
			fmt.Printf("  s=%.2f EFL=%v\n", s, pp.EFL)
		}
	}
	// --- NIR doublet
	{
		const l1, l2 = 1064.0, 1550.0
		efl := 200.0
		seed := seedAchromat("FusedSilica", "N-SF5", efl)
		cards := []lens.Card{
			{Name: "stop", RadiusMM: math.Inf(1), ThicknessMM: 3, SemiAperMM: 20},
			{Name: "S1", RadiusMM: 1 / seed[0], ThicknessMM: 6, Glass: "FusedSilica", SemiAperMM: 26},
			{Name: "S2", RadiusMM: 1 / seed[1], ThicknessMM: 3, Glass: "N-SF5", SemiAperMM: 26},
			{Name: "S3", RadiusMM: 1 / seed[2], ThicknessMM: 60, Glass: "", SemiAperMM: 26},
		}
		p := &lens.Prescription{ObjectZMM: math.Inf(-1), Stop: 0, Cards: cards}
		for _, s := range []float64{0.9, 1.0, 1.1} {
			p.Cards[1].RadiusMM = 1 / (seed[0] * s)
			p.Cards[2].RadiusMM = -1 / (seed[0] * s)
			pa, _ := p.Paraxial(l1)
			pb, _ := p.Paraxial(l2)
			fmt.Printf("  nir s(seed scale)=%.2f EFL1064=%v EFL1550=%v\n", s, pa.EFL, pb.EFL)
		}
	}
}
