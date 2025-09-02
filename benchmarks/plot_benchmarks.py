#!/usr/bin/env python3
"""
Benchmark plotter for noLZSS factorization.

This script runs the bench.py benchmark with various input sizes
and creates plots of time vs. input size and memory usage vs. input size.
"""

import subprocess
import re
import matplotlib.pyplot as plt
import argparse
import sys
import os
from pathlib import Path

def run_benchmark(size, runs=3):
    """Run bench.py with given size and return parsed results."""
    cmd = [sys.executable, "benchmarks/bench.py", "--size", str(size), "--runs", str(runs)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd())
        if result.returncode != 0:
            print(f"Error running benchmark for size {size}: {result.stderr}")
            return None

        output = result.stdout

        # Parse the output using regex
        patterns = {
            'mean_time': r'Mean time per run: ([\d.]+) ms',
            'mean_memory': r'Mean peak memory usage: ([\d.]+) MB',
            'max_memory': r'Max peak memory usage: ([\d.]+) MB',
            'throughput': r'Throughput: ([\d.]+) MB/s',
            'factors_per_sec': r'Factors/second \(avg\): ([\d,]+)'
        }

        results = {}
        for key, pattern in patterns.items():
            match = re.search(pattern, output)
            if match:
                value = match.group(1).replace(',', '')  # Remove commas from numbers
                results[key] = float(value)
            else:
                print(f"Warning: Could not find {key} in output for size {size}")
                results[key] = None

        return results

    except Exception as e:
        print(f"Error running benchmark for size {size}: {e}")
        return None

def create_plots(sizes, results, output_dir="benchmarks/plots"):
    """Create and save plots of the benchmark results."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Filter out None results
    valid_results = [(size, res) for size, res in zip(sizes, results) if res is not None]
    if not valid_results:
        print("No valid results to plot!")
        return

    sizes_plot = [size for size, _ in valid_results]
    results_plot = [res for _, res in valid_results]

    # Extract metrics
    mean_times = [res['mean_time'] for res in results_plot if res['mean_time'] is not None]
    mean_memories = [res['mean_memory'] for res in results_plot if res['mean_memory'] is not None]
    throughputs = [res['throughput'] for res in results_plot if res['throughput'] is not None]

    # Convert sizes to MB for plotting
    sizes_mb = [size / 1e6 for size in sizes_plot]

    # Plot 1: Time vs Input Size
    plt.figure(figsize=(12, 8))

    plt.subplot(2, 2, 1)
    if mean_times:
        plt.plot(sizes_mb, mean_times, 'bo-', linewidth=2, markersize=8)
        plt.xlabel('Input Size (MB)')
        plt.ylabel('Mean Time per Run (ms)')
        plt.title('Factorization Time vs Input Size')
        plt.grid(True, alpha=0.3)
        plt.xscale('log')
        plt.yscale('log')

    # Plot 2: Memory vs Input Size
    plt.subplot(2, 2, 2)
    if mean_memories:
        plt.plot(sizes_mb, mean_memories, 'ro-', linewidth=2, markersize=8)
        plt.xlabel('Input Size (MB)')
        plt.ylabel('Mean Peak Memory Usage (MB)')
        plt.title('Memory Usage vs Input Size')
        plt.grid(True, alpha=0.3)
        plt.xscale('log')
        plt.yscale('log')

    # Plot 3: Throughput vs Input Size
    plt.subplot(2, 2, 3)
    if throughputs:
        plt.plot(sizes_mb, throughputs, 'go-', linewidth=2, markersize=8)
        plt.xlabel('Input Size (MB)')
        plt.ylabel('Throughput (MB/s)')
        plt.title('Throughput vs Input Size')
        plt.grid(True, alpha=0.3)
        plt.xscale('log')

    # Plot 4: Time per Byte vs Input Size
    plt.subplot(2, 2, 4)
    if mean_times:
        time_per_byte = [t / (s * 1000) for t, s in zip(mean_times, sizes_plot)]  # Convert to seconds per byte
        plt.plot(sizes_mb, time_per_byte, 'mo-', linewidth=2, markersize=8)
        plt.xlabel('Input Size (MB)')
        plt.ylabel('Time per Byte (μs)')
        plt.title('Time per Byte vs Input Size')
        plt.grid(True, alpha=0.3)
        plt.xscale('log')
        plt.yscale('log')

    plt.tight_layout()
    plt.savefig(f"{output_dir}/benchmark_plots.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/benchmark_plots.pdf", bbox_inches='tight')
    plt.show()

    # Print summary table
    print("\nBenchmark Results Summary:")
    print("-" * 80)
    print(f"{'Size (MB)':<10} {'Time (ms)':<12} {'Memory (MB)':<12} {'Throughput (MB/s)':<18}")
    print("-" * 80)
    for size_mb, res in zip(sizes_mb, results_plot):
        time_str = f"{res['mean_time']:.2f}" if res['mean_time'] else "N/A"
        mem_str = f"{res['mean_memory']:.2f}" if res['mean_memory'] else "N/A"
        tp_str = f"{res['throughput']:.2f}" if res['throughput'] else "N/A"
        print(f"{size_mb:<10.1f} {time_str:<12} {mem_str:<12} {tp_str:<18}")

def main():
    parser = argparse.ArgumentParser(description="Create benchmark plots for noLZSS")
    parser.add_argument("--sizes", nargs="+", type=int,
                       default=[10000, 50000, 100000, 500000, 1000000, 5000000],
                       help="Input sizes to benchmark (in bytes)")
    parser.add_argument("--runs", type=int, default=3,
                       help="Number of runs per benchmark")
    parser.add_argument("--output-dir", default="benchmarks/plots",
                       help="Output directory for plots")

    args = parser.parse_args()

    print(f"Running benchmarks for sizes: {args.sizes}")
    print(f"Number of runs per size: {args.runs}")

    results = []
    for size in args.sizes:
        print(f"\nBenchmarking size: {size:,} bytes ({size/1e6:.1f} MB)")
        result = run_benchmark(size, args.runs)
        results.append(result)

    create_plots(args.sizes, results, args.output_dir)
    print(f"\nPlots saved to {args.output_dir}/")

if __name__ == "__main__":
    main()
