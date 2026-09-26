package powerflow

import (
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// Column headers expected in the bus file (dadosbarra.csv).
var busColumns = []string{"bar_number", "typebar", "P", "Q", "V", "Theta", "Cshunt"}

// Column headers expected in the line file (dadoslinha.csv).
var lineColumns = []string{"line_number", "origin", "destiny", "R_line", "X_line", "Xshunt", "TapValue"}

// LoadBusesFile reads bus data from a CSV file.
func LoadBusesFile(path string) ([]Bus, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	buses, err := ReadBuses(f)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return buses, nil
}

// LoadLinesFile reads line data from a CSV file.
func LoadLinesFile(path string) ([]Line, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	lines, err := ReadLines(f)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", path, err)
	}
	return lines, nil
}

// ReadBuses parses bus data with the columns
// bar_number,typebar,P,Q,V,Theta,Cshunt (Theta in radians).
func ReadBuses(r io.Reader) ([]Bus, error) {
	rows, err := readTable(r, busColumns)
	if err != nil {
		return nil, err
	}
	buses := make([]Bus, 0, len(rows))
	for _, v := range rows {
		buses = append(buses, Bus{
			Number: int(v[0]),
			Type:   BusType(int(v[1])),
			P:      v[2],
			Q:      v[3],
			V:      v[4],
			Theta:  v[5],
			Shunt:  v[6],
		})
	}
	return buses, nil
}

// ReadLines parses line data with the columns
// line_number,origin,destiny,R_line,X_line,Xshunt,TapValue.
func ReadLines(r io.Reader) ([]Line, error) {
	rows, err := readTable(r, lineColumns)
	if err != nil {
		return nil, err
	}
	lines := make([]Line, 0, len(rows))
	for _, v := range rows {
		lines = append(lines, Line{
			Number:  int(v[0]),
			Origin:  int(v[1]),
			Destiny: int(v[2]),
			R:       v[3],
			X:       v[4],
			BShunt:  v[5],
			Tap:     v[6],
		})
	}
	return lines, nil
}

// readTable reads a CSV with a header row and returns the requested columns,
// in the requested order, as floats.
func readTable(r io.Reader, columns []string) ([][]float64, error) {
	cr := csv.NewReader(r)
	cr.TrimLeadingSpace = true
	records, err := cr.ReadAll()
	if err != nil {
		return nil, err
	}
	if len(records) == 0 {
		return nil, fmt.Errorf("empty file")
	}
	pos := make(map[string]int, len(records[0]))
	for i, h := range records[0] {
		pos[strings.TrimSpace(strings.TrimPrefix(h, "\ufeff"))] = i
	}
	idx := make([]int, len(columns))
	for i, c := range columns {
		p, ok := pos[c]
		if !ok {
			return nil, fmt.Errorf("missing column %q", c)
		}
		idx[i] = p
	}

	rows := make([][]float64, 0, len(records)-1)
	for n, rec := range records[1:] {
		row := make([]float64, len(columns))
		for i, p := range idx {
			f, err := strconv.ParseFloat(strings.TrimSpace(rec[p]), 64)
			if err != nil {
				return nil, fmt.Errorf("row %d, column %q: %w", n+2, columns[i], err)
			}
			row[i] = f
		}
		rows = append(rows, row)
	}
	return rows, nil
}
