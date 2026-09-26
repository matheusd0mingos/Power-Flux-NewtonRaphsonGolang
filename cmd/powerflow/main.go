// Command powerflow runs a Newton-Raphson load flow on the bus and line data
// given as CSV files.
package main

import (
	"flag"
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"os"

	"github.com/matheusd0mingos/Power-Flux-NewtonRaphsonGolang/powerflow"
)

func main() {
	busFile := flag.String("bus", "data/dadosbarra.csv", "bus data CSV")
	lineFile := flag.String("line", "data/dadoslinha.csv", "line data CSV")
	gs := flag.Bool("gs", false, "run Gauss-Seidel first for better initial values")
	gsIter := flag.Int("gs-iter", 1, "maximum Gauss-Seidel sweeps")
	maxIter := flag.Int("max-iter", powerflow.DefaultOptions.MaxIterations, "maximum Newton-Raphson iterations")
	tol := flag.Float64("tol", powerflow.DefaultOptions.Tolerance, "power mismatch tolerance (p.u.)")
	flag.Parse()

	if err := run(os.Stdout, *busFile, *lineFile, *gs, *gsIter, powerflow.Options{MaxIterations: *maxIter, Tolerance: *tol}); err != nil {
		fmt.Fprintln(os.Stderr, "error:", err)
		os.Exit(1)
	}
}

func run(w io.Writer, busFile, lineFile string, gs bool, gsIter int, opt powerflow.Options) error {
	buses, err := powerflow.LoadBusesFile(busFile)
	if err != nil {
		return err
	}
	lines, err := powerflow.LoadLinesFile(lineFile)
	if err != nil {
		return err
	}
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
