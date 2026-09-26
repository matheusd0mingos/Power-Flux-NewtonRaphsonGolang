package powerflow

import (
	"context"
	"errors"
	"fmt"
	"math"
	"reflect"
	"strings"
	"testing"
)

func sampleData(t testing.TB) ([]Bus, []Line) {
	t.Helper()
	b, err := ReadBuses(strings.NewReader(busCSV))
	if err != nil {
		t.Fatal(err)
	}
	l, err := ReadLines(strings.NewReader(lineCSV))
	if err != nil {
		t.Fatal(err)
	}
	return b, l
}

// gridSystem builds an n×n meshed grid: bus 1 is the slack, every fifth bus a
// PV generator and the rest PQ loads.
func gridSystem(n int) ([]Bus, []Line) {
	var buses []Bus
	var lines []Line
	for i := range n * n {
		b := Bus{Number: i + 1, Type: PQ, P: -0.1, Q: -0.03, V: 1}
		switch {
		case i == 0:
			b.Type = Slack
		case i%5 == 0:
			b.Type, b.P, b.V = PV, 0.3, 1.02
		}
		buses = append(buses, b)
	}
	add := func(o, d int) {
		lines = append(lines, Line{Number: len(lines) + 1, Origin: o + 1, Destiny: d + 1, R: 0.01, X: 0.05, BShunt: 0.01, Tap: 1})
	}
	for r := range n {
		for c := range n {
			if c+1 < n {
				add(r*n+c, r*n+c+1)
			}
			if r+1 < n {
				add(r*n+c, (r+1)*n+c)
			}
		}
	}
	return buses, lines
}

var tightOptions = Options{MaxIterations: 30, Tolerance: 1e-8}

func TestContingenciesSample(t *testing.T) {
	buses, lines := sampleData(t)
	res, err := Contingencies(context.Background(), buses, lines, DefaultOptions, 0)
	if err != nil {
		t.Fatal(err)
	}
	if len(res) != len(lines) {
		t.Fatalf("got %d results, want %d", len(res), len(lines))
	}
	// Losing line 1-2 is survivable; losing either line feeding bus 3 is not.
	want := []bool{true, false, false}
	for i, r := range res {
		if r.Outage.Number != i+1 || r.Converged != want[i] {
			t.Errorf("outage %d: converged=%v, want %v", r.Outage.Number, r.Converged, want[i])
		}
	}
	near(t, "min V", res[0].Operating.MinV, 0.94498, 1e-5)
}

// Running on one worker or many must give identical results.
func TestParallelMatchesSequential(t *testing.T) {
	buses, lines := gridSystem(5)
	ctx := context.Background()
	seq, err := Contingencies(ctx, buses, lines, tightOptions, 1)
	if err != nil {
		t.Fatal(err)
	}
	par, err := Contingencies(ctx, buses, lines, tightOptions, 8)
	if err != nil {
		t.Fatal(err)
	}
	if !reflect.DeepEqual(seq, par) {
		t.Fatal("parallel contingency results differ from sequential")
	}
	for _, r := range par {
		if !r.Converged {
			t.Errorf("outage of line %d did not converge in a meshed grid", r.Outage.Number)
		}
	}
}

func TestContingencyIslanding(t *testing.T) {
	// Radial feeder 1-2-3: every outage islands a bus.
	buses := []Bus{
		{Number: 1, Type: Slack, V: 1},
		{Number: 2, Type: PQ, P: -0.5},
		{Number: 3, Type: PQ, P: -0.5},
	}
	lines := []Line{
		{Number: 1, Origin: 1, Destiny: 2, X: 0.1, Tap: 1},
		{Number: 2, Origin: 2, Destiny: 3, X: 0.1, Tap: 1},
	}
	res, err := Contingencies(context.Background(), buses, lines, DefaultOptions, 0)
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range res {
		if !r.Islanded || r.Converged {
			t.Errorf("outage %d: %+v, want islanded", r.Outage.Number, r)
		}
	}
}

func TestLoadSweep(t *testing.T) {
	buses, lines := sampleData(t)
	lambdas := []float64{0.5, 1, 1.5, 1.55, 1.6, 2}
	pts, err := LoadSweep(context.Background(), buses, lines, lambdas, Options{MaxIterations: 30, Tolerance: 1e-6}, 0)
	if err != nil {
		t.Fatal(err)
	}
	limit, ok := LoadabilityLimit(pts)
	if !ok || limit != 1.55 {
		t.Errorf("limit = %v (%v), want 1.55", limit, ok)
	}
	// λ = 1 is the base case.
	near(t, "V3 at λ=1", pts[1].Voltages[2], 0.97229, 1e-4)
	// Voltage at the load bus must fall as the load grows.
	for i := 1; i < 4; i++ {
		if pts[i].Voltages[2] >= pts[i-1].Voltages[2] {
			t.Errorf("V3 did not decrease from λ=%v to λ=%v", lambdas[i-1], lambdas[i])
		}
	}
	// The input must not be modified by the sweep.
	if buses[2].P != -4.5 {
		t.Errorf("input bus modified: P = %v", buses[2].P)
	}
}

func TestParallelCancellation(t *testing.T) {
	ctx, cancel := context.WithCancel(context.Background())
	calls := 0
	_, err := parallel(ctx, 1000, 1, func(_ context.Context, i int) int {
		calls++ // single worker, so no race
		if i == 10 {
			cancel()
		}
		return i
	})
	if !errors.Is(err, context.Canceled) {
		t.Fatalf("err = %v, want context.Canceled", err)
	}
	if calls >= 1000 {
		t.Errorf("all %d jobs ran despite cancellation", calls)
	}
}

func TestGridSystemSolves(t *testing.T) {
	buses, lines := gridSystem(10)
	s, err := NewSystem(buses, lines)
	if err != nil {
		t.Fatal(err)
	}
	res, err := s.NewtonRaphson(tightOptions)
	if err != nil || !res.Converged {
		t.Fatalf("res=%+v err=%v", res, err)
	}
	for _, b := range s.Buses {
		if math.IsNaN(b.V) || b.V < 0.8 || b.V > 1.1 {
			t.Fatalf("bus %d: V = %v", b.Number, b.V)
		}
	}
}

// go test -bench Contingencies -run ^$ ./powerflow
func BenchmarkContingencies(b *testing.B) {
	buses, lines := gridSystem(8) // 64 buses, 112 lines → 112 load flows
	for _, w := range []int{1, 2, 4, 8} {
		b.Run(fmt.Sprintf("workers=%d", w), func(b *testing.B) {
			for b.Loop() {
				if _, err := Contingencies(context.Background(), buses, lines, tightOptions, w); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
