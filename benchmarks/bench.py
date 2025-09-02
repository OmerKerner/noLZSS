import time
import sys
import statistics as stats
import argparse
import noLZSS
import random
import tracemalloc

def run_once(data: bytes):
    tracemalloc.start()
    t0 = time.perf_counter()
    factors = noLZSS.factorize(data)
    t1 = time.perf_counter()
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return t1 - t0, factors, peak

def generate_data(size: int) -> bytes:
    # Generate a random nucleotide sequence (A, C, G, T) of the requested size
    seq = ''.join(random.choices('ACGT', k=size))
    return seq.encode('ascii')

def main():
    parser = argparse.ArgumentParser(description="Benchmark noLZSS factorization")
    parser.add_argument("file", nargs="?", help="Input file to benchmark (optional)")
    parser.add_argument("--size", type=int, default=1000000, help="Size of generated data in bytes (default: 1MB)")
    parser.add_argument("--runs", type=int, default=10, help="Number of benchmark runs (default: 10)")
    args = parser.parse_args()

    if args.file:
        with open(args.file, 'rb') as f:
            data = f.read()
    else:
        data = generate_data(args.size)

    print(f"Benchmarking with {len(data):,} bytes of data")
    print(f"Running {args.runs} benchmark iterations...")

    results = [run_once(data) for _ in range(args.runs)]
    times = [t for t, f, m in results]
    counts = [len(f) for t, f, m in results]
    memories = [m for t, f, m in results]

    mean_time = stats.mean(times)
    median_time = stats.median(times)
    stdev = stats.stdev(times) if len(times) > 1 else 0
    min_time = min(times)
    max_time = max(times)

    mean_factors = stats.mean(counts)
    total_factors = sum(counts)

    mean_memory = stats.mean(memories)
    max_memory = max(memories)

    factors_per_sec = mean_factors / mean_time if mean_time > 0 else float('inf')
    time_per_iter_ms = mean_time * 1000
    time_per_factor_ms = (mean_time * 1000) / mean_factors if mean_factors > 0 else float('nan')
    throughput_mb_s = (len(data) / mean_time) / 1e6 if mean_time > 0 else float('inf')

    print("\nResults:")
    print(f"  Mean time per run: {time_per_iter_ms:.2f} ms")
    print(f"  Median time per run: {median_time * 1000:.2f} ms")
    print(f"  Std dev: {stdev * 1000:.2f} ms")
    print(f"  Min time: {min_time * 1000:.2f} ms")
    print(f"  Max time: {max_time * 1000:.2f} ms")
    print(f"  Mean factors per run: {mean_factors:.2f}")
    print(f"  Total factors across runs: {total_factors}")
    print(f"  Factors/second (avg): {factors_per_sec:,.2f}")
    print(f"  Time per factor (avg): {time_per_factor_ms:.6f} ms")
    print(f"  Throughput: {throughput_mb_s:.2f} MB/s")
    print(f"  Mean peak memory usage: {mean_memory / 1024 / 1024:.2f} MB")
    print(f"  Max peak memory usage: {max_memory / 1024 / 1024:.2f} MB")

if __name__ == '__main__':
    main()
