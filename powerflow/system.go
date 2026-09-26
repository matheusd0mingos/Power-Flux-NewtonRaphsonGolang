// Package powerflow solves the AC load-flow problem with the Newton-Raphson
// method, optionally warm-started by a few Gauss-Seidel sweeps.
//
// Line data are given in per-unit impedance (R, X) while shunts are given in
// per-unit susceptance.
package powerflow

import "fmt"

// BusType identifies how a bus is modelled in the load flow.
type BusType int

const (
	// Slack is the reference bus (V and Theta specified).
	Slack BusType = 1
	// PV is a voltage-controlled bus (P and V specified).
	PV BusType = 2
	// PQ is a load bus (P and Q specified).
	PQ BusType = 3
)

func (t BusType) String() string {
	switch t {
	case Slack:
		return "VTheta"
	case PV:
		return "PV"
	case PQ:
		return "PQ"
	}
	return fmt.Sprintf("BusType(%d)", int(t))
}

// Bus is a node of the network. Angles are stored in radians.
type Bus struct {
	Number int
	Type   BusType
	P      float64 // specified (or, after solving, calculated) active injection
	Q      float64 // specified (or, after solving, calculated) reactive injection
	V      float64
	Theta  float64
	Shunt  float64 // bus shunt susceptance
}

// Line is a branch between two buses. A tap different from 1 models an
// off-nominal transformer: the origin side sees y, the destiny side y/t².
type Line struct {
	Number  int
	Origin  int // 1-based bus number
	Destiny int // 1-based bus number
	R       float64
	X       float64
	BShunt  float64 // shunt susceptance added at each end
	Tap     float64
}

// System holds the network and its bus admittance matrix.
type System struct {
	Buses []*Bus
	Lines []Line
	Y     [][]complex128
}

// NewSystem builds a system with a flat start (PQ buses at 1∠0, PV buses at
// V∠0) and computes the admittance matrix.
func NewSystem(buses []Bus, lines []Line) (*System, error) {
	s := &System{Lines: lines}
	index := make(map[int]int, len(buses))
	for i, b := range buses {
		if b.Number != i+1 {
			return nil, fmt.Errorf("bus %d found at position %d: buses must be numbered 1..n in order", b.Number, i+1)
		}
		switch b.Type {
		case Slack:
		case PV:
			b.Theta = 0
		case PQ:
			b.V, b.Theta = 1, 0
		default:
			return nil, fmt.Errorf("bus %d: invalid type %d (1=VTheta, 2=PV, 3=PQ)", b.Number, int(b.Type))
		}
		bus := b
		s.Buses = append(s.Buses, &bus)
		index[b.Number] = i
	}
	for _, l := range lines {
		if _, ok := index[l.Origin]; !ok {
			return nil, fmt.Errorf("line %d: unknown origin bus %d", l.Number, l.Origin)
		}
		if _, ok := index[l.Destiny]; !ok {
			return nil, fmt.Errorf("line %d: unknown destiny bus %d", l.Number, l.Destiny)
		}
		if l.R == 0 && l.X == 0 {
			return nil, fmt.Errorf("line %d: zero impedance", l.Number)
		}
	}
	s.Y = AdmittanceMatrix(s.Buses, lines)
	return s, nil
}

// AdmittanceMatrix builds the bus admittance matrix (Y-bus).
func AdmittanceMatrix(buses []*Bus, lines []Line) [][]complex128 {
	n := len(buses)
	y := newComplexMatrix(n)
	for i, b := range buses {
		y[i][i] += complex(0, b.Shunt)
	}
	for _, l := range lines {
		o, d := l.Origin-1, l.Destiny-1
		yl := 1 / complex(l.R, l.X)
		t := complex(tapOf(l), 0)
		sh := complex(0, l.BShunt)

		y[o][o] += yl + sh
		y[d][d] += yl/(t*t) + sh
		y[o][d] -= yl / t
		y[d][o] -= yl / t
	}
	return y
}

func tapOf(l Line) float64 {
	if l.Tap == 0 {
		return 1
	}
	return l.Tap
}

func newComplexMatrix(n int) [][]complex128 {
	m := make([][]complex128, n)
	for i := range m {
		m[i] = make([]complex128, n)
	}
	return m
}
