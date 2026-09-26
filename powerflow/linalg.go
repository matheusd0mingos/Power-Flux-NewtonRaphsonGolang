package powerflow

import (
	"errors"
	"math"
)

// ErrSingular is returned when the Jacobian cannot be factorized.
var ErrSingular = errors.New("singular matrix")

// solve returns x such that a·x = b using Gaussian elimination with partial
// pivoting. a and b are not modified.
func solve(a [][]float64, b []float64) ([]float64, error) {
	n := len(b)
	m := make([][]float64, n)
	for i := range a {
		m[i] = make([]float64, n+1)
		copy(m[i], a[i])
		m[i][n] = b[i]
	}

	for col := 0; col < n; col++ {
		pivot := col
		for r := col + 1; r < n; r++ {
			if math.Abs(m[r][col]) > math.Abs(m[pivot][col]) {
				pivot = r
			}
		}
		if math.Abs(m[pivot][col]) < 1e-14 {
			return nil, ErrSingular
		}
		m[col], m[pivot] = m[pivot], m[col]
		for r := col + 1; r < n; r++ {
			f := m[r][col] / m[col][col]
			if f == 0 {
				continue
			}
			for c := col; c <= n; c++ {
				m[r][c] -= f * m[col][c]
			}
		}
	}

	x := make([]float64, n)
	for r := n - 1; r >= 0; r-- {
		sum := m[r][n]
		for c := r + 1; c < n; c++ {
			sum -= m[r][c] * x[c]
		}
		x[r] = sum / m[r][r]
	}
	return x, nil
}
