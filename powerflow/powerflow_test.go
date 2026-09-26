package powerflow

import (
	"math"
	"strings"
	"testing"
)

const busCSV = `bar_number,typebar,P,Q,V,Theta,Cshunt
1,1,0,0,1,0,0
2,2,3,0,1.1,0,0
3,3,-4.5,-0.5,1,0,0
`

const lineCSV = `line_number,origin,destiny,R_line,X_line,Xshunt,TapValue
1,1,2,0,0.1,0,1.00
2,1,3,0,0.1,0,1.05
3,2,3,0,0.2,0,0.98
`

func loadSystem(t *testing.T, buses, lines string) *System {
	t.Helper()
	b, err := ReadBuses(strings.NewReader(buses))
	if err != nil {
		t.Fatal(err)
	}
	l, err := ReadLines(strings.NewReader(lines))
	if err != nil {
		t.Fatal(err)
	}
	s, err := NewSystem(b, l)
	if err != nil {
		t.Fatal(err)
	}
	return s
}

func near(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.6f, want %.6f", name, got, want)
	}
}

// Reference values produced by the original Python implementation
// (sistema_new.py) on the same data.
func TestMatchesPythonReference(t *testing.T) {
	cases := []struct {
		name  string
		gs    bool
		iters int
		want  [3][4]float64 // V, Theta (deg), P, Q
	}{
		{"flat start", false, 3, [3][4]float64{
			{1, 0, 1.49995, 0.13970},
			{1.1, 5.35780, 3, 2.11047},
			{0.97229, -15.83727, -4.5, -0.5},
		}},
		{"gauss-seidel start", true, 2, [3][4]float64{
			{1, 0, 1.49978, 0.13916},
			{1.1, 5.35781, 3, 2.11015},
			{0.97233, -15.83540, -4.5, -0.5},
		}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			s := loadSystem(t, busCSV, lineCSV)
			if c.gs {
				s.GaussSeidel(1e-6, 1)
			}
			res, err := s.NewtonRaphson(DefaultOptions)
			if err != nil {
				t.Fatal(err)
			}
			if !res.Converged || res.Iterations != c.iters {
				t.Fatalf("result = %+v, want convergence in %d iterations", res, c.iters)
			}
			for i, b := range s.Buses {
				near(t, "V", b.V, c.want[i][0], 1e-5)
				near(t, "Theta", b.Theta*180/math.Pi, c.want[i][1], 1e-5)
				near(t, "P", b.P, c.want[i][2], 1e-5)
				near(t, "Q", b.Q, c.want[i][3], 1e-5)
			}
		})
	}
}

func TestAdmittanceMatrix(t *testing.T) {
	s := loadSystem(t, busCSV, lineCSV)
	want := [3][3]float64{
		{-20, 10, 10 / 1.05},
		{10, -15, 5 / 0.98},
		{10 / 1.05, 5 / 0.98, -10/(1.05*1.05) - 5/(0.98*0.98)},
	}
	for i := range want {
		for j := range want[i] {
			near(t, "real(Y)", real(s.Y[i][j]), 0, 1e-12)
			near(t, "imag(Y)", imag(s.Y[i][j]), want[i][j], 1e-9)
		}
	}
}

// Line flows leaving each bus must add up to its injection.
func TestLineFlowsBalance(t *testing.T) {
	s := loadSystem(t, busCSV, lineCSV)
	if _, err := s.NewtonRaphson(Options{MaxIterations: 20, Tolerance: 1e-10}); err != nil {
		t.Fatal(err)
	}
	p, q := s.CalculatePower()
	sumP := make([]float64, len(s.Buses))
	sumQ := make([]float64, len(s.Buses))
	for _, f := range s.LineFlows() {
		sumP[f.From-1] += f.P
		sumQ[f.From-1] += f.Q
	}
	for i := range s.Buses {
		near(t, "ΣP", sumP[i], p[i], 1e-9)
		near(t, "ΣQ", sumQ[i], q[i], 1e-9)
	}
}

// The analytic Jacobian must agree with finite differences of the mismatch.
func TestJacobianFiniteDifference(t *testing.T) {
	s := loadSystem(t, busCSV, lineCSV)
	s.Buses[1].Theta, s.Buses[2].Theta, s.Buses[2].V = 0.1, -0.2, 0.95
	angle, magnitude := s.unknowns()
	j := s.Jacobian()
	const h = 1e-7

	vars := make([]*float64, 0)
	for _, i := range angle {
		vars = append(vars, &s.Buses[i].Theta)
	}
	for _, i := range magnitude {
		vars = append(vars, &s.Buses[i].V)
	}
	for c, x := range vars {
		*x += h
		up := s.Mismatch()
		*x -= 2 * h
		down := s.Mismatch()
		*x += h
		for r := range up {
			// Mismatch = specified - calculated, so ∂calc/∂x = -(Δmismatch)/Δx.
			fd := -(up[r] - down[r]) / (2 * h)
			near(t, "J", j[r][c], fd, 1e-5)
		}
	}
}

// A PQ bus listed before a PV bus must be handled correctly.
func TestBusOrderIndependence(t *testing.T) {
	buses := `bar_number,typebar,P,Q,V,Theta,Cshunt
1,1,0,0,1,0,0
2,3,-4.5,-0.5,1,0,0
3,2,3,0,1.1,0,0
`
	lines := `line_number,origin,destiny,R_line,X_line,Xshunt,TapValue
1,1,3,0,0.1,0,1.00
2,1,2,0,0.1,0,1.05
3,3,2,0,0.2,0,0.98
`
	s := loadSystem(t, buses, lines)
	res, err := s.NewtonRaphson(DefaultOptions)
	if err != nil || !res.Converged {
		t.Fatalf("res=%+v err=%v", res, err)
	}
	near(t, "V PQ", s.Buses[1].V, 0.97229, 1e-5)
	near(t, "Theta PV", s.Buses[2].Theta*180/math.Pi, 5.35780, 1e-5)
}

func TestReadErrors(t *testing.T) {
	if _, err := ReadBuses(strings.NewReader("bar_number,typebar\n1,1\n")); err == nil {
		t.Error("expected missing column error")
	}
	if _, err := ReadLines(strings.NewReader(strings.Replace(lineCSV, "0.1,0,1.00", "x,0,1.00", 1))); err == nil {
		t.Error("expected parse error")
	}
	b, _ := ReadBuses(strings.NewReader(strings.Replace(busCSV, "3,3,", "3,7,", 1)))
	l, _ := ReadLines(strings.NewReader(lineCSV))
	if _, err := NewSystem(b, l); err == nil {
		t.Error("expected invalid bus type error")
	}
}
