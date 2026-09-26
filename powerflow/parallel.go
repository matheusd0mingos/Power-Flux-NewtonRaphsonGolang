package powerflow

import (
	"context"
	"runtime"
	"sync"
)

// parallel runs fn(ctx, i) for i in [0, n) on a pool of worker goroutines and
// returns the results in index order. Each result is written to its own slot,
// so no locking is needed. If ctx is cancelled, pending jobs are skipped and
// ctx.Err() is returned.
func parallel[T any](ctx context.Context, n, workers int, fn func(ctx context.Context, i int) T) ([]T, error) {
	if workers <= 0 {
		workers = runtime.NumCPU()
	}
	workers = min(workers, n)

	results := make([]T, n)
	jobs := make(chan int)
	var wg sync.WaitGroup
	for range workers {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range jobs {
				results[i] = fn(ctx, i)
			}
		}()
	}

feed:
	for i := range n {
		select {
		case jobs <- i:
		case <-ctx.Done():
			break feed
		}
	}
	close(jobs)
	wg.Wait()
	return results, ctx.Err()
}
