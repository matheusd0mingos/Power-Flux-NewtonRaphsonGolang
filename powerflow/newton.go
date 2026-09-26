package powerflow

import (
	"math"
	"math/cmplx"
)

// CalculatePower returns the active and reactive power injected at each bus
// for the current voltages.
func (s *System) CalculatePower() (p, q []float64) {
	n := len(s.Buses)
	p = make([]float64, n)
	q = make([]float64, n)
	for i, bi := range s.Buses {
		for j, bj := range s.Buses {
			g, b := real(s.Y[i][j]), imag(s.Y[i][j])
			if g == 0 && b == 0 {
				continue
			}
			th := bi.Theta - bj.Theta
			sin, cos := math.Sincos(th)
			p[i] += bi.V * bj.V * (g*cos + b*sin)
			q[i] += bi.V * bj.V * (g*sin - b*cos)
		}
	}
	return p, q
}

// unknowns lists the bus indices whose angle (non-slack) and magnitude (PQ)
// are state variables, in bus order.
func (s *System) unknowns() (angle, magnitude []int) {
	for i, b := range s.Buses {
		if b.Type != Slack {
			angle = append(angle, i)
		}
		if b.Type == PQ {
			magnitude = append(magnitude, i)
		}
	}
	return angle, magnitude
}

// Mismatch returns ΔP for every non-slack bus followed by ΔQ for every PQ bus.
func (s *System) Mismatch() []float64 {
	p, q := s.CalculatePower()
	angle, magnitude := s.unknowns()
	m := make([]float64, 0, len(angle)+len(magnitude))
	for _, i := range angle {
		m = append(m, s.Buses[i].P-p[i])
	}
	for _, i := range magnitude {
		m = append(m, s.Buses[i].Q-q[i])
	}
	return m
}

// Jacobian returns the reduced Jacobian
//
//	| H  N |   H = ∂P/∂θ   N = ∂P/∂V
//	| M  L |   M = ∂Q/∂θ   L = ∂Q/∂V
//
// with rows/columns ordered as in Mismatch.
func (s *System) Jacobian() [][]float64 {
	p, q := s.CalculatePower()
	angle, magnitude := s.unknowns()
	na := len(angle)
	size := na + len(magnitude)

	j := make([][]float64, size)
	for r := range j {
		j[r] = make([]float64, size)
	}

	// dP returns ∂P_k/∂θ_m and ∂P_k/∂V_m; dQ the same for Q_k.
	dP := func(k, m int) (dth, dv float64) {
		bk, bm := s.Buses[k], s.Buses[m]
		g, b := real(s.Y[k][m]), imag(s.Y[k][m])
		if k == m {
			return -q[k] - bk.V*bk.V*b, p[k]/bk.V + bk.V*g
		}
		sin, cos := math.Sincos(bk.Theta - bm.Theta)
		return bk.V * bm.V * (g*sin - b*cos), bk.V * (g*cos + b*sin)
	}
	dQ := func(k, m int) (dth, dv float64) {
		bk, bm := s.Buses[k], s.Buses[m]
		g, b := real(s.Y[k][m]), imag(s.Y[k][m])
		if k == m {
			return p[k] - bk.V*bk.V*g, q[k]/bk.V - bk.V*b
		}
		sin, cos := math.Sincos(bk.Theta - bm.Theta)
		return -bk.V * bm.V * (g*cos + b*sin), bk.V * (g*sin - b*cos)
	}

	for r, k := range angle {
		for c, m := range angle {
			j[r][c], _ = dP(k, m)
		}
		for c, m := range magnitude {
			_, j[r][na+c] = dP(k, m)
		}
	}
	for r, k := range magnitude {
		for c, m := range angle {
			j[na+r][c], _ = dQ(k, m)
		}
		for c, m := range magnitude {
			_, j[na+r][na+c] = dQ(k, m)
		}
	}
	return j
}

// Options configures NewtonRaphson.
type Options struct {
	MaxIterations int
	Tolerance     float64
	// OnIteration, if set, is called after every voltage update.
	OnIteration func(iteration int, s *System)
}

// DefaultOptions matches the original Python implementation.
var DefaultOptions = Options{MaxIterations: 10, Tolerance: 3e-3}

// Result summarizes a Newton-Raphson run.
type Result struct {
	Converged  bool
	Iterations int // number of voltage updates performed
	MaxError   float64
}

// NewtonRaphson iterates until every power mismatch is below the tolerance.
// On convergence the calculated slack P/Q and PV Q are written back to the
// buses.
func (s *System) NewtonRaphson(opt Options) (Result, error) {
	angle, magnitude := s.unknowns()
	var res Result
	for it := 0; ; it++ {
		mismatch := s.Mismatch()
		res.Iterations = it
		res.MaxError = maxAbs(mismatch)
		if res.MaxError < opt.Tolerance {
			res.Converged = true
			s.updatePowerValues()
			return res, nil
		}
		if it == opt.MaxIterations {
			return res, nil
		}

		delta, err := solve(s.Jacobian(), mismatch)
		if err != nil {
			return res, err
		}
		for r, i := range angle {
			s.Buses[i].Theta = wrapAngle(s.Buses[i].Theta + delta[r])
		}
		for r, i := range magnitude {
			s.Buses[i].V += delta[len(angle)+r]
		}
		if opt.OnIteration != nil {
			opt.OnIteration(it+1, s)
		}
	}
}

func (s *System) updatePowerValues() {
	p, q := s.CalculatePower()
	for i, b := range s.Buses {
		switch b.Type {
		case Slack:
			b.P, b.Q = p[i], q[i]
		case PV:
			b.Q = q[i]
		}
	}
}

// GaussSeidel performs up to maxIter Gauss-Seidel sweeps to improve the
// initial voltages before Newton-Raphson.
func (s *System) GaussSeidel(tol float64, maxIter int) {
	n := len(s.Buses)
	prev := make([]complex128, n)
	for i, b := range s.Buses {
		prev[i] = cmplx.Rect(b.V, b.Theta)
	}
	next := append([]complex128(nil), prev...)

	for it := 0; it < maxIter; it++ {
		copy(next, prev)
		for i, b := range s.Buses {
			if b.Type == Slack {
				continue
			}
			var sum complex128
			for m := 0; m < n; m++ {
				switch {
				case m < i:
					sum += s.Y[i][m] * next[m]
				case m > i:
					sum += s.Y[i][m] * prev[m]
				}
			}
			q := b.Q
			if b.Type == PV {
				var yv complex128
				for m := 0; m < n; m++ {
					yv += s.Y[i][m] * prev[m]
				}
				q = -imag(cmplx.Conj(prev[i]) * yv)
			}
			v := (complex(b.P, -q)/cmplx.Conj(prev[i]) - sum) / s.Y[i][i]
			if b.Type == PV {
				v = cmplx.Rect(cmplx.Abs(prev[i]), cmplx.Phase(v))
			}
			next[i] = v
		}

		done := true
		for i := range next {
			if cmplx.Abs(next[i]-prev[i]) > tol+tol*cmplx.Abs(prev[i]) {
				done = false
				break
			}
		}
		if done {
			break
		}
		copy(prev, next)
	}

	for i, b := range s.Buses {
		b.V, b.Theta = cmplx.Abs(next[i]), cmplx.Phase(next[i])
	}
}

// Flow is the complex power leaving a bus through a line.
type Flow struct {
	Line     int
	From, To int
	P, Q     float64
}

// LineFlows returns, for each line, the power flow in both directions
// (origin→destiny first), consistent with the tap model of the Y-bus.
func (s *System) LineFlows() []Flow {
	flows := make([]Flow, 0, 2*len(s.Lines))
	for _, l := range s.Lines {
		o, d := s.Buses[l.Origin-1], s.Buses[l.Destiny-1]
		vo, vd := cmplx.Rect(o.V, o.Theta), cmplx.Rect(d.V, d.Theta)
		y := 1 / complex(l.R, l.X)
		t := complex(tapOf(l), 0)
		sh := complex(0, l.BShunt)

		io := (y+sh)*vo - y/t*vd
		id := (y/(t*t)+sh)*vd - y/t*vo
		so, sd := vo*cmplx.Conj(io), vd*cmplx.Conj(id)

		flows = append(flows,
			Flow{Line: l.Number, From: l.Origin, To: l.Destiny, P: real(so), Q: imag(so)},
			Flow{Line: l.Number, From: l.Destiny, To: l.Origin, P: real(sd), Q: imag(sd)},
		)
	}
	return flows
}

func maxAbs(v []float64) float64 {
	m := 0.0
	for _, x := range v {
		m = math.Max(m, math.Abs(x))
	}
	return m
}

// wrapAngle maps an angle to (-π, π].
func wrapAngle(a float64) float64 {
	a = math.Mod(a, 2*math.Pi)
	if a > math.Pi {
		a -= 2 * math.Pi
	} else if a <= -math.Pi {
		a += 2 * math.Pi
	}
	return a
}
