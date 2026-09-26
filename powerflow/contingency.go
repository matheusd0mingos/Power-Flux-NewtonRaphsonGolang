package powerflow

import (
	"context"
	"math"
	"math/cmplx"
)

// Operating summarizes a solved operating point.
type Operating struct {
	MinV, MaxV       float64
	MinVBus, MaxVBus int
	MaxFlow          float64 // largest |S| leaving either end of a line
	MaxFlowLine      int
}

func (s *System) operating() Operating {
	op := Operating{MinV: math.Inf(1), MaxV: math.Inf(-1)}
	for _, b := range s.Buses {
		if b.V < op.MinV {
			op.MinV, op.MinVBus = b.V, b.Number
		}
		if b.V > op.MaxV {
			op.MaxV, op.MaxVBus = b.V, b.Number
		}
	}
	for _, f := range s.LineFlows() {
		if a := cmplx.Abs(complex(f.P, f.Q)); a > op.MaxFlow {
			op.MaxFlow, op.MaxFlowLine = a, f.Line
		}
	}
	return op
}

// ContingencyResult is the outcome of solving the network with one line out.
type ContingencyResult struct {
	Outage    Line
	Islanded  bool // the outage splits the network; not solved
	Converged bool
	Result    Result
	Err       error
	Operating Operating // valid when Converged
}

// Contingencies runs an N-1 analysis: for every line, the load flow is solved
// from a flat start with that line removed. The cases are solved concurrently
// on the given number of worker goroutines (<= 0 means one per CPU).
func Contingencies(ctx context.Context, buses []Bus, lines []Line, opt Options, workers int) ([]ContingencyResult, error) {
	opt.OnIteration = nil // a callback would be called from many goroutines at once
	return parallel(ctx, len(lines), workers, func(_ context.Context, k int) ContingencyResult {
		res := ContingencyResult{Outage: lines[k]}
		remaining := make([]Line, 0, len(lines)-1)
		remaining = append(remaining, lines[:k]...)
		remaining = append(remaining, lines[k+1:]...)

		if !connected(buses, remaining) {
			res.Islanded = true
			return res
		}
		s, err := NewSystem(buses, remaining)
		if err != nil {
			res.Err = err
			return res
		}
		res.Result, res.Err = s.NewtonRaphson(opt)
		res.Converged = res.Err == nil && res.Result.Converged
		if res.Converged {
			res.Operating = s.operating()
		}
		return res
	})
}

// connected reports whether every bus can reach a slack bus.
func connected(buses []Bus, lines []Line) bool {
	adj := make([][]int, len(buses))
	for _, l := range lines {
		o, d := l.Origin-1, l.Destiny-1
		adj[o] = append(adj[o], d)
		adj[d] = append(adj[d], o)
	}
	seen := make([]bool, len(buses))
	var queue []int
	for i, b := range buses {
		if b.Type == Slack {
			seen[i] = true
			queue = append(queue, i)
		}
	}
	for len(queue) > 0 {
		i := queue[0]
		queue = queue[1:]
		for _, j := range adj[i] {
			if !seen[j] {
				seen[j] = true
				queue = append(queue, j)
			}
		}
	}
	for _, ok := range seen {
		if !ok {
			return false
		}
	}
	return true
}

// SweepPoint is the outcome of solving the network at one load level.
type SweepPoint struct {
	Lambda    float64
	Converged bool
	Result    Result
	Err       error
	Voltages  []float64 // bus voltage magnitudes, valid when Converged
	Operating Operating
}

// LoadSweep solves the load flow with the specified P and Q of every non-slack
// bus multiplied by each value in lambdas, concurrently on the given number of
// workers. It traces the PV (nose) curve of the system.
func LoadSweep(ctx context.Context, buses []Bus, lines []Line, lambdas []float64, opt Options, workers int) ([]SweepPoint, error) {
	opt.OnIteration = nil // a callback would be called from many goroutines at once
	return parallel(ctx, len(lambdas), workers, func(_ context.Context, k int) SweepPoint {
		pt := SweepPoint{Lambda: lambdas[k]}
		scaled := make([]Bus, len(buses))
		for i, b := range buses {
			if b.Type != Slack {
				b.P *= pt.Lambda
				b.Q *= pt.Lambda
			}
			scaled[i] = b
		}
		s, err := NewSystem(scaled, lines)
		if err != nil {
			pt.Err = err
			return pt
		}
		pt.Result, pt.Err = s.NewtonRaphson(opt)
		pt.Converged = pt.Err == nil && pt.Result.Converged
		if pt.Converged {
			pt.Voltages = make([]float64, len(s.Buses))
			for i, b := range s.Buses {
				pt.Voltages[i] = b.V
			}
			pt.Operating = s.operating()
		}
		return pt
	})
}

// LoadabilityLimit returns the largest λ reached from the first point through
// consecutive converged points of a sweep sorted by λ, and false if the first
// point did not converge.
func LoadabilityLimit(points []SweepPoint) (float64, bool) {
	limit, ok := 0.0, false
	for _, p := range points {
		if !p.Converged {
			break
		}
		limit, ok = p.Lambda, true
	}
	return limit, ok
}
