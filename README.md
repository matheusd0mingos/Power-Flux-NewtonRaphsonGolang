# Power Flux Newton-Raphson (Go)

A Go rebuild of [Power-Flux-NewtonRaphson](https://github.com/matheusd0mingos/Power-Flux-NewtonRaphson):
an AC load-flow solver using the Newton-Raphson method, with optional
Gauss-Seidel sweeps to improve the initial values. No dependencies beyond the
Go standard library.

## Usage

```sh
go run ./cmd/powerflow                 # uses data/dadosbarra.csv and data/dadoslinha.csv
go run ./cmd/powerflow -gs             # Gauss-Seidel warm start (1 sweep)
go run ./cmd/powerflow -bus my_buses.csv -line my_lines.csv -tol 1e-6 -max-iter 20
go run ./cmd/powerflow -n1                                  # N-1 contingency analysis, in parallel
go run ./cmd/powerflow -sweep-max 2 -tol 1e-6 -max-iter 30  # load sweep / loadability limit, in parallel
```

| Flag        | Default               | Description                              |
|-------------|-----------------------|------------------------------------------|
| `-bus`      | `data/dadosbarra.csv` | Bus data CSV                             |
| `-line`     | `data/dadoslinha.csv` | Line data CSV                            |
| `-gs`       | `false`               | Run Gauss-Seidel before Newton-Raphson   |
| `-gs-iter`  | `1`                   | Maximum Gauss-Seidel sweeps              |
| `-max-iter` | `10`                  | Maximum Newton-Raphson iterations        |
| `-tol`      | `0.003`               | Power mismatch tolerance (p.u.)          |
| `-n1`       | `false`               | Run the N-1 contingency analysis         |
| `-sweep-max`| `0` (off)             | Sweep the load factor λ up to this value |
| `-sweep-step`| `0.05`               | Load factor step                         |
| `-workers`  | number of CPUs        | Goroutines for `-n1` and `-sweep-max`    |

The program prints the admittance matrix, the voltages at each iteration, the
power flow in both directions of every line (with losses) and the final bus
voltages and injections.

## Parallel studies (goroutines)

A single load flow is fast and sequential by nature, but many power system
studies need *hundreds of independent* load flows. Those are spread over a
pool of worker goroutines (`powerflow/parallel.go`). Each job builds its own
`System`, results are written to their own slot (no locks), and Ctrl-C
cancels the run through a `context.Context`.

**N-1 contingency analysis** (`-n1`, `powerflow.Contingencies`): removes each
line in turn and solves from a flat start. It reports whether the case
converges, the min/max bus voltage and the most loaded line. It also detects
outages that island part of the network.

```
N-1 contingency analysis: 3 cases on 3 goroutines in 115µs
Line out | Buses   | Status       | Iter | Min V (bus)   | Max V (bus)   | Max |S| (line)
1        |   1-2   | converged    | 3    | 0.94498 (3  ) | 1.10000 (2  ) | 3.43596 (3)
2        |   1-3   | diverged     | 10   |
3        |   2-3   | diverged     | 10   |
```

**Load sweep** (`-sweep-max`, `powerflow.LoadSweep`): multiplies P and Q of
every non-slack bus by λ, solves every level at the same time and reports
the largest λ that still converges, an estimate of the loadability limit
(the nose of the PV curve). The sample system reaches λ ≈ 1.55.

Speedup on a synthetic 64-bus meshed grid (112 contingencies, 4-core CPU):

```
go test -bench Contingencies -run '^$' ./powerflow
BenchmarkContingencies/workers=1   200 ms/op
BenchmarkContingencies/workers=2   111 ms/op
BenchmarkContingencies/workers=4    78 ms/op
```

## Input data

All values in per unit. Line data are impedances; shunts are susceptances.

`dadosbarra.csv`

| Column       | Meaning                                         |
|--------------|-------------------------------------------------|
| `bar_number` | Bus number (1..n, in order)                     |
| `typebar`    | 1 = VΘ (slack), 2 = PV, 3 = PQ                  |
| `P`, `Q`     | Specified active / reactive injection           |
| `V`          | Voltage magnitude (used for slack and PV buses) |
| `Theta`      | Voltage angle in radians (used for slack bus)   |
| `Cshunt`     | Bus shunt susceptance                           |

`dadoslinha.csv`

| Column        | Meaning                                              |
|---------------|------------------------------------------------------|
| `line_number` | Line identifier                                      |
| `origin`      | Sending bus                                          |
| `destiny`     | Receiving bus                                        |
| `R_line`      | Series resistance                                    |
| `X_line`      | Series reactance                                     |
| `Xshunt`      | Shunt susceptance added at each end                  |
| `TapValue`    | Off-nominal tap (0 or 1 = none), on the destiny side |

## Library

The solver lives in package `powerflow` and can be used directly:

```go
buses, _ := powerflow.LoadBusesFile("data/dadosbarra.csv")
lines, _ := powerflow.LoadLinesFile("data/dadoslinha.csv")
sys, _ := powerflow.NewSystem(buses, lines)
res, err := sys.NewtonRaphson(powerflow.Options{MaxIterations: 20, Tolerance: 1e-8})
flows := sys.LineFlows()
```

## Differences from the Python version

Bus voltages, angles and injections are identical to the Python script on the
same data (checked in `TestMatchesPythonReference`). A few issues were fixed
along the way:

- **Line flows account for the transformer tap**, consistent with the Y
  matrix, so the flows leaving each bus add up to its injection. The Python
  version ignored the tap in the flow calculation.
- **Bus ordering**: angle/voltage corrections are mapped back by bus index,
  so a PQ bus may appear before a PV bus in the file.
- **Bus shunt (`Cshunt`)** is added to the Y matrix once per bus rather than
  once per connected line.
- Linear systems are solved by Gaussian elimination instead of inverting the
  Jacobian, and angles are kept in (-π, π].
- The Gauss-Seidel option is a command-line flag instead of an interactive
  prompt.

## Tests

```sh
go test -race ./...
```

## License

MIT, see [LICENSE](LICENSE).
