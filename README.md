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
```

| Flag        | Default               | Description                              |
|-------------|-----------------------|------------------------------------------|
| `-bus`      | `data/dadosbarra.csv` | Bus data CSV                             |
| `-line`     | `data/dadoslinha.csv` | Line data CSV                            |
| `-gs`       | `false`               | Run Gauss-Seidel before Newton-Raphson   |
| `-gs-iter`  | `1`                   | Maximum Gauss-Seidel sweeps              |
| `-max-iter` | `10`                  | Maximum Newton-Raphson iterations        |
| `-tol`      | `0.003`               | Power mismatch tolerance (p.u.)          |

The program prints the admittance matrix, the voltages at each iteration, the
power flow in both directions of every line (with losses) and the final bus
voltages and injections.

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
go test ./...
```

## License

MIT, see [LICENSE](LICENSE).
