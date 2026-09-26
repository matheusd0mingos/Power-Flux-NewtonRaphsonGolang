// Command powerflow runs a Newton-Raphson load flow on the bus and line data
// given as CSV files.
package main

import (
	"context"
	"flag"
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"os"
	"os/signal"
	"runtime"
	"time"

	"github.com/matheusd0mingos/Power-Flux-NewtonRaphsonGolang/powerflow"
)

func main() {
	busFile := flag.String("bus", "data/dadosbarra.csv", "bus data CSV")
	lineFile := flag.String("line", "data/dadoslinha.csv", "line data CSV")
	gs := flag.Bool("gs", false, "run Gauss-Seidel first for better initial values")
	gsIter := flag.Int("gs-iter", 1, "maximum Gauss-Seidel sweeps")
	maxIter := flag.Int("max-iter", powerflow.DefaultOptions.MaxIterations, "maximum Newton-Raphson iterations")
	tol := flag.Float64("tol", powerflow.DefaultOptions.Tolerance, "power mismatch tolerance (p.u.)")
	n1 := flag.Bool("n1", false, "run an N-1 contingency analysis (one line out at a time) in parallel")
	sweepMax := flag.Float64("sweep-max", 0, "if > 0, sweep the load factor λ from -sweep-step up to this value in parallel")
	sweepStep := flag.Float64("sweep-step", 0.05, "load factor step for -sweep-max")
	workers := flag.Int("workers", runtime.NumCPU(), "worker goroutines for -n1 and -sweep-max")
	flag.Parse()

	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt)
	defer stop()

	opt := powerflow.Options{MaxIterations: *maxIter, Tolerance: *tol}
	err := func() error {
		buses, err := powerflow.LoadBusesFile(*busFile)
		if err != nil {
			return err
		}
		lines, err := powerflow.LoadLinesFile(*lineFile)
		if err != nil {
			return err
		}
		switch {
		case *n1:
			return runContingencies(ctx, os.Stdout, buses, lines, opt, *workers)
		case *sweepMax > 0:
			return runSweep(ctx, os.Stdout, buses, lines, *sweepStep, *sweepMax, opt, *workers)
		}
		return run(os.Stdout, buses, lines, *gs, *gsIter, opt)
	}()
	if err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}

func runContingencies(ctx context.Context, w io.Writer, buses []powerflow.Bus, lines []powerflow.Line, opt powerflow.Options, workers int) error {
	start := time.Now()
	results, err := powerflow.Contingencies(ctx, buses, lines, opt, workers)
	if err != nil {
		return err
	}
	fmt.Fprintf(w, "N-1 contingency analysis: %d cases on %d goroutines in %v\n",
		len(results), min(workers, len(results)), time.Since(start).Round(time.Microsecond))
	fmt.Fprintln(w, "-------------------------------------------------------------------------------------")
	fmt.Fprintln(w, "Line out | Buses   | Status       | Iter | Min V (bus)   | Max V (bus)   | Max |S| (line)")
	fmt.Fprintln(w, "-------------------------------------------------------------------------------------")
	for _, r := range results {
		fmt.Fprintf(w, "%-9d| %3d-%-4d| ", r.Outage.Number, r.Outage.Origin, r.Outage.Destiny)
		switch {
		case r.Islanded:
			fmt.Fprintln(w, "islanded")
		case r.Err != nil:
			fmt.Fprintf(w, "failed: %v\n", r.Err)
		case !r.Converged:
			fmt.Fprintf(w, "%-13s| %-5d|\n", "diverged", r.Result.Iterations)
		default:
			op := r.Operating
			fmt.Fprintf(w, "%-13s| %-5d| %.5f (%-3d) | %.5f (%-3d) | %.5f (%d)\n", "converged", r.Result.Iterations,
				op.MinV, op.MinVBus, op.MaxV, op.MaxVBus, op.MaxFlow, op.MaxFlowLine)
		}
	}
	return nil
}

func runSweep(ctx context.Context, w io.Writer, buses []powerflow.Bus, lines []powerflow.Line, step, max float64, opt powerflow.Options, workers int) error {
	if step <= 0 {
		return fmt.Errorf("-sweep-step must be positive")
	}
	var lambdas []float64
	for i := 1; float64(i)*step <= max+1e-9; i++ {
		lambdas = append(lambdas, float64(i)*step)
	}
	start := time.Now()
	points, err := powerflow.LoadSweep(ctx, buses, lines, lambdas, opt, workers)
	if err != nil {
		return err
	}
	fmt.Fprintf(w, "Load sweep (λ scales P and Q of every non-slack bus): %d points on %d goroutines in %v\n",
		len(points), min(workers, len(points)), time.Since(start).Round(time.Microsecond))
	fmt.Fprintln(w, "-------------------------------------------------")
	fmt.Fprintln(w, "λ        | Status    | Iter | Min V (bus)")
	fmt.Fprintln(w, "-------------------------------------------------")
	for _, p := range points {
		if !p.Converged {
			fmt.Fprintf(w, "%-9.3f| %-10s| %-5d|\n", p.Lambda, "diverged", p.Result.Iterations)
			continue
		}
		fmt.Fprintf(w, "%-9.3f| %-10s| %-5d| %.5f (%d)\n", p.Lambda, "converged", p.Result.Iterations,
			p.Operating.MinV, p.Operating.MinVBus)
	}
	if limit, ok := powerflow.LoadabilityLimit(points); ok {
		fmt.Fprintf(w, "\nLoadability limit: λ ≈ %.3f (within one step of %.3f)\n", limit, step)
	} else {
		fmt.Fprintln(w, "\nThe first load level did not converge.")
	}
	return nil
}

func run(w io.Writer, buses []powerflow.Bus, lines []powerflow.Line, gs bool, gsIter int, opt powerflow.Options) error {
	sys, err := powerflow.NewSystem(buses, lines)
	if err != nil {
		return err
	}

	fmt.Fprintln(w, "Y matrix:")
	printMatrix(w, sys.Y)

	if gs {
		sys.GaussSeidel(1e-6, gsIter)
	}

	fmt.Fprintln(w, "\nInitial conditions:")
	printBuses(w, sys, 0)

	opt.OnIteration = func(it int, s *powerflow.System) { printVoltages(w, s, it) }
	res, err := sys.NewtonRaphson(opt)
	if err != nil {
		return err
	}
	if !res.Converged {
		return fmt.Errorf("did not converge after %d iterations (max mismatch %.3g)", res.Iterations, res.MaxError)
	}

	fmt.Fprintf(w, "\nConverged in %d iterations (max mismatch %.3g p.u.)\n\n", res.Iterations, res.MaxError)
	printFlows(w, sys.LineFlows())
	fmt.Fprintln(w)
	printBuses(w, sys, res.Iterations)
	return nil
}

func printMatrix(w io.Writer, y [][]complex128) {
	for _, row := range y {
		for _, v := range row {
			fmt.Fprintf(w, " %9.4f%+9.4fj", real(v), imag(v))
		}
		fmt.Fprintln(w)
	}
}

func printVoltages(w io.Writer, s *powerflow.System, it int) {
	fmt.Fprintf(w, "Iteration %d\n", it)
	fmt.Fprintln(w, "-----------------------------------------------")
	fmt.Fprintln(w, "Bus Number | Voltage Magnitude | Voltage Angle ")
	fmt.Fprintln(w, "-----------------------------------------------")
	for _, b := range s.Buses {
		fmt.Fprintf(w, "%-11d | %-17.5f | %-13.5f\n", b.Number, b.V, b.Theta)
	}
}

func printBuses(w io.Writer, s *powerflow.System, it int) {
	fmt.Fprintf(w, "Number of iterations: %d\n", it)
	fmt.Fprintln(w, "-------------------------------------------------------------------------------------------")
	fmt.Fprintln(w, "Bus Number | Type   | Voltage Magnitude | Voltage Angle (deg) | Active Power | Reactive Power")
	fmt.Fprintln(w, "-------------------------------------------------------------------------------------------")
	for _, b := range s.Buses {
		fmt.Fprintf(w, "%-11d| %-7s| %-18.5f| %-20.5f| %-13.5f| %-14.5f\n",
			b.Number, b.Type, b.V, b.Theta*180/math.Pi, b.P, b.Q)
	}
}

func printFlows(w io.Writer, flows []powerflow.Flow) {
	fmt.Fprintln(w, "Line flows:")
	fmt.Fprintln(w, "----------------------------------------------")
	fmt.Fprintln(w, "Line | From | To   | P (p.u.)    | Q (p.u.)")
	fmt.Fprintln(w, "----------------------------------------------")
	for i, f := range flows {
		fmt.Fprintf(w, "%-5d| %-5d| %-5d| %-12.5f| %-12.5f\n", f.Line, f.From, f.To, f.P, f.Q)
		if i%2 == 1 {
			loss := complex(flows[i-1].P+f.P, flows[i-1].Q+f.Q)
			fmt.Fprintf(w, "     losses: P = %.5f  Q = %.5f  (|S| = %.5f)\n", real(loss), imag(loss), cmplx.Abs(loss))
		}
	}
}
