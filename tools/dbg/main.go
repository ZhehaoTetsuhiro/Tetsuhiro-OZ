package main

import (
	"fmt"
	"math"

	"tetsuhiro-oz/geo"
	"tetsuhiro-oz/lens"
)

func main() {
	p := &lens.Prescription{
		Name: "galilean 4x", ObjectZMM: math.Inf(-1), Stop: 0,
		Cards: []lens.Card{
			{Name: "stop", RadiusMM: math.Inf(1), ThicknessMM: 5, SemiAperMM: 5},
			{Name: "neg lens", RadiusMM: 2 * (1.5168 - 1) * -25, ThicknessMM: 3, Glass: "N-BK7", SemiAperMM: 8},
			{Name: "neg rear", RadiusMM: -2 * (1.5168 - 1) * -25, ThicknessMM: 75, Glass: "", SemiAperMM: 8},
			{Name: "pos front", RadiusMM: 2 * (1.5168 - 1) * 100, ThicknessMM: 5, Glass: "N-BK7", SemiAperMM: 25},
			{Name: "pos rear", RadiusMM: -2 * (1.5168 - 1) * 100, ThicknessMM: 50, Glass: "", SemiAperMM: 25},
		},
	}
	sys, err := p.System()
	fmt.Println("system err:", err)
	for i, s := range sys.Surfaces {
		fmt.Printf("  surf %d: %-18s vignetted-r=%v z=%v\n", i, s.Label(), s.MaxRadius(), sys.VertexZ(i))
	}
	println("thickness", len(sys.Thickness))
	res := sys.Trace(geo.NewRay(geo.V(0, 0.5, -30), geo.V(0, 0, 1), geo.Wavelengthd))
	fmt.Println("status:", res.Status, "failedAt:", res.FailedAt, "err:", res.Error())
	for _, st := range res.Steps {
		fmt.Printf("  step %d %-12s status=%-10v hit=%s dir=%s\n", st.Surface, st.Label, st.Status, st.Point, st.DirOut)
	}
	par, perr := sys.Paraxial(geo.Wavelengthd, math.Inf(-1))
	fmt.Println("paraxial err:", perr)
	if perr == nil {
		fmt.Printf("  EFL=%v BFL=%v F#=%v ENPz=%v ENPr=%v exitPupil=%v/%v\n", par.EFL, par.BFL, par.FNumber, par.EntrancePupilZ, par.EntrancePupilRadius, par.ExitPupilZ, par.ExitPupilRadius)
	}
}
