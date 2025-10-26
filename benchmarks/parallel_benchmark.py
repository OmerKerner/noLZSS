#!/usr/bin/env python3
"""
Parallel factorization benchmark script for noLZSS.

This script benchmarks the parallel factorization functions with time and memory measurement
across different input sizes and thread counts. Creates plots comparing parallel speedup
and efficiency.

Functions benchmarked:
- parallel_factorize_to_file() - Parallel factorization to binary file
- parallel_factorize_file_to_file() - Parallel factorization from file to binary file
- parallel_factorize_dna_w_rc_to_file() - Parallel DNA factorization with reverse complement
- parallel_factorize_file_dna_w_rc_to_file() - Parallel DNA factorization from file
"""

import time
import tracemalloc
import tempfile
import os
import random
import argparse
import json
import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Any
import statistics as stats
import multiprocessing

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as scipy_stats

try:
    import noLZSS._noLZSS as cpp
except ImportError:
    print("Error: C++ extension not available. Please build the package first.")
    print("Run: pip install -e .")
    exit(1)

# Verify all required functions are available
REQUIRED_FUNCTIONS = [
    'parallel_factorize_to_file',
    'parallel_factorize_file_to_file',
    'parallel_factorize_dna_w_rc_to_file',
    'parallel_factorize_file_dna_w_rc_to_file',
]

missing_functions = [func for func in REQUIRED_FUNCTIONS if not hasattr(cpp, func)]
if missing_functions:
    print(f"Error: C++ extension is missing required functions: {', '.join(missing_functions)}")
    print("This usually means the package needs to be rebuilt.")
    print("Run: pip install -e . --no-build-isolation --force-reinstall")
    exit(1)


def generate_text_content(size: int, dna: bool = False) -> str:
    """
    Generate random text content for benchmarking.
    
    Args:
        size: Size of the text in characters
        dna: Whether to generate DNA sequence (ACGT)
        
    Returns:
        Random text string
    """
    if dna:
        alphabet = 'ACGT'
    else:
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    return ''.join(random.choices(alphabet, k=size))


def create_temp_file(content: str) -> str:
    """Create a temporary text file and return its path."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write(content)
        return f.name


def benchmark_parallel_function(func, *args, runs: int = 3) -> Dict[str, Any]:
    """
    Benchmark a parallel factorization function.
    
    Args:
        func: Function to benchmark
        *args: Arguments to pass to the function
        runs: Number of runs to average
        
    Returns:
        Dictionary with timing and memory statistics
    """
    times = []
    memories = []
    num_factors_list = []
    
    for _ in range(runs):
        tracemalloc.start()
        start_time = time.perf_counter()
        
        try:
            num_factors = func(*args)
            end_time = time.perf_counter()
            current, peak = tracemalloc.get_traced_memory()
            
            times.append(end_time - start_time)
            memories.append(peak)
            num_factors_list.append(num_factors)
            
        except Exception as e:
            print(f"Error in benchmark: {e}")
            return None
        finally:
            tracemalloc.stop()
    
    return {
        'mean_time': stats.mean(times),
        'std_time': stats.stdev(times) if len(times) > 1 else 0,
        'min_time': min(times),
        'max_time': max(times),
        'mean_memory_mb': stats.mean(memories) / (1024 * 1024),
        'max_memory_mb': max(memories) / (1024 * 1024),
        'mean_factors': stats.mean(num_factors_list),
        'all_times': times,
        'all_memories': memories
    }


def run_parallel_benchmark_suite(
    sizes: List[int], 
    thread_counts: List[int], 
    runs: int = 3, 
    output_dir: str = "benchmarks/parallel_results"
) -> Dict[str, Any]:
    """
    Run the complete benchmark suite for parallel factorization functions.
    
    Args:
        sizes: List of input sizes in characters
        thread_counts: List of thread counts to test
        runs: Number of runs per benchmark
        output_dir: Directory to save results
        
    Returns:
        Complete benchmark results
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Test regular text parallel functions
    print("\n=== Benchmarking parallel_factorize_to_file ===")
    results['parallel_factorize_to_file'] = {'sizes': sizes, 'thread_counts': thread_counts, 'results': {}}
    
    for size in sizes:
        print(f"\n  Size: {size:,} characters")
        results['parallel_factorize_to_file']['results'][size] = {}
        
        # Generate text content
        text_content = generate_text_content(size, dna=False)
        
        for num_threads in thread_counts:
            print(f"    Threads: {num_threads}")
            
            # Create temporary output file
            with tempfile.NamedTemporaryFile(delete=False) as f:
                output_path = f.name
            
            try:
                result = benchmark_parallel_function(
                    cpp.parallel_factorize_to_file,
                    text_content,
                    output_path,
                    num_threads,
                    runs=runs
                )
                
                if result:
                    result['input_size'] = size
                    result['num_threads'] = num_threads
                    results['parallel_factorize_to_file']['results'][size][num_threads] = result
                    
                    time_ms = result['mean_time'] * 1000
                    mem_mb = result['mean_memory_mb']
                    print(f"      Time: {time_ms:.2f} ms, Memory: {mem_mb:.2f} MB")
                else:
                    print(f"      Failed")
                    
            finally:
                # Clean up
                try:
                    os.unlink(output_path)
                except:
                    pass
    
    # Test file-to-file parallel functions
    print("\n=== Benchmarking parallel_factorize_file_to_file ===")
    results['parallel_factorize_file_to_file'] = {'sizes': sizes, 'thread_counts': thread_counts, 'results': {}}
    
    for size in sizes:
        print(f"\n  Size: {size:,} characters")
        results['parallel_factorize_file_to_file']['results'][size] = {}
        
        # Generate text content and create file
        text_content = generate_text_content(size, dna=False)
        input_path = create_temp_file(text_content)
        
        try:
            for num_threads in thread_counts:
                print(f"    Threads: {num_threads}")
                
                # Create temporary output file
                with tempfile.NamedTemporaryFile(delete=False) as f:
                    output_path = f.name
                
                try:
                    result = benchmark_parallel_function(
                        cpp.parallel_factorize_file_to_file,
                        input_path,
                        output_path,
                        num_threads,
                        runs=runs
                    )
                    
                    if result:
                        result['input_size'] = size
                        result['num_threads'] = num_threads
                        results['parallel_factorize_file_to_file']['results'][size][num_threads] = result
                        
                        time_ms = result['mean_time'] * 1000
                        mem_mb = result['mean_memory_mb']
                        print(f"      Time: {time_ms:.2f} ms, Memory: {mem_mb:.2f} MB")
                    else:
                        print(f"      Failed")
                        
                finally:
                    # Clean up output file
                    try:
                        os.unlink(output_path)
                    except:
                        pass
        finally:
            # Clean up input file
            try:
                os.unlink(input_path)
            except:
                pass
    
    # Test DNA with reverse complement parallel functions
    print("\n=== Benchmarking parallel_factorize_dna_w_rc_to_file ===")
    results['parallel_factorize_dna_w_rc_to_file'] = {'sizes': sizes, 'thread_counts': thread_counts, 'results': {}}
    
    for size in sizes:
        print(f"\n  Size: {size:,} characters")
        results['parallel_factorize_dna_w_rc_to_file']['results'][size] = {}
        
        # Generate DNA content
        dna_content = generate_text_content(size, dna=True)
        
        for num_threads in thread_counts:
            print(f"    Threads: {num_threads}")
            
            # Create temporary output file
            with tempfile.NamedTemporaryFile(delete=False) as f:
                output_path = f.name
            
            try:
                result = benchmark_parallel_function(
                    cpp.parallel_factorize_dna_w_rc_to_file,
                    dna_content,
                    output_path,
                    num_threads,
                    runs=runs
                )
                
                if result:
                    result['input_size'] = size
                    result['num_threads'] = num_threads
                    results['parallel_factorize_dna_w_rc_to_file']['results'][size][num_threads] = result
                    
                    time_ms = result['mean_time'] * 1000
                    mem_mb = result['mean_memory_mb']
                    print(f"      Time: {time_ms:.2f} ms, Memory: {mem_mb:.2f} MB")
                else:
                    print(f"      Failed")
                    
            finally:
                # Clean up
                try:
                    os.unlink(output_path)
                except:
                    pass
    
    # Test DNA file-to-file with reverse complement parallel functions
    print("\n=== Benchmarking parallel_factorize_file_dna_w_rc_to_file ===")
    results['parallel_factorize_file_dna_w_rc_to_file'] = {'sizes': sizes, 'thread_counts': thread_counts, 'results': {}}
    
    for size in sizes:
        print(f"\n  Size: {size:,} characters")
        results['parallel_factorize_file_dna_w_rc_to_file']['results'][size] = {}
        
        # Generate DNA content and create file
        dna_content = generate_text_content(size, dna=True)
        input_path = create_temp_file(dna_content)
        
        try:
            for num_threads in thread_counts:
                print(f"    Threads: {num_threads}")
                
                # Create temporary output file
                with tempfile.NamedTemporaryFile(delete=False) as f:
                    output_path = f.name
                
                try:
                    result = benchmark_parallel_function(
                        cpp.parallel_factorize_file_dna_w_rc_to_file,
                        input_path,
                        output_path,
                        num_threads,
                        runs=runs
                    )
                    
                    if result:
                        result['input_size'] = size
                        result['num_threads'] = num_threads
                        results['parallel_factorize_file_dna_w_rc_to_file']['results'][size][num_threads] = result
                        
                        time_ms = result['mean_time'] * 1000
                        mem_mb = result['mean_memory_mb']
                        print(f"      Time: {time_ms:.2f} ms, Memory: {mem_mb:.2f} MB")
                    else:
                        print(f"      Failed")
                        
                finally:
                    # Clean up output file
                    try:
                        os.unlink(output_path)
                    except:
                        pass
        finally:
            # Clean up input file
            try:
                os.unlink(input_path)
            except:
                pass
    
    return results


def create_speedup_plots(results: Dict[str, Any], output_dir: str):
    """
    Create speedup and efficiency plots for parallel benchmarks.
    
    Args:
        results: Benchmark results
        output_dir: Directory to save plots
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Parallel Factorization Benchmarks', fontsize=16)
    
    colors = ['blue', 'red', 'green', 'orange']
    markers = ['o', 's', '^', 'D']
    
    for idx, (func_name, data) in enumerate(results.items()):
        if 'results' not in data or not data['results']:
            continue
        
        color = colors[idx % len(colors)]
        marker = markers[idx % len(markers)]
        
        # Select a representative size for detailed analysis
        sizes = list(data['results'].keys())
        if not sizes:
            continue
        representative_size = sizes[len(sizes) // 2]  # Middle size
        
        size_results = data['results'][representative_size]
        thread_counts = sorted(size_results.keys())
        
        if not thread_counts or 1 not in thread_counts:
            continue
        
        # Get baseline (single-threaded) performance
        baseline_time = size_results[1]['mean_time']
        
        # Calculate speedup
        speedups = []
        efficiencies = []
        times = []
        
        for num_threads in thread_counts:
            thread_time = size_results[num_threads]['mean_time']
            speedup = baseline_time / thread_time
            efficiency = speedup / num_threads
            
            speedups.append(speedup)
            efficiencies.append(efficiency * 100)  # Convert to percentage
            times.append(thread_time * 1000)  # Convert to ms
        
        # Plot 1: Speedup vs Thread Count
        ax = axes[0, 0]
        ax.plot(thread_counts, speedups, color=color, marker=marker, 
               label=func_name.replace('_', ' '), linewidth=2)
        # Add ideal speedup line
        if idx == 0:
            ax.plot(thread_counts, thread_counts, 'k--', alpha=0.5, label='Ideal speedup')
        
        # Plot 2: Efficiency vs Thread Count
        ax = axes[0, 1]
        ax.plot(thread_counts, efficiencies, color=color, marker=marker,
               label=func_name.replace('_', ' '), linewidth=2)
        # Add 100% efficiency line
        if idx == 0:
            ax.axhline(y=100, color='k', linestyle='--', alpha=0.5, label='100% efficiency')
        
        # Plot 3: Time vs Thread Count
        ax = axes[1, 0]
        ax.plot(thread_counts, times, color=color, marker=marker,
               label=func_name.replace('_', ' '), linewidth=2)
        
        # Plot 4: Time vs Size for different thread counts
        ax = axes[1, 1]
        # Plot for maximum thread count
        max_threads = max(thread_counts)
        size_list = []
        time_list = []
        for size in sorted(data['results'].keys()):
            if max_threads in data['results'][size]:
                size_list.append(size)
                time_list.append(data['results'][size][max_threads]['mean_time'] * 1000)
        
        if size_list and time_list:
            ax.loglog(size_list, time_list, color=color, marker=marker,
                     label=f"{func_name.replace('_', ' ')} ({max_threads} threads)", linewidth=2)
    
    # Configure Plot 1
    axes[0, 0].set_xlabel('Number of Threads')
    axes[0, 0].set_ylabel('Speedup')
    axes[0, 0].set_title(f'Speedup vs Thread Count (Size: {representative_size:,})')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend()
    
    # Configure Plot 2
    axes[0, 1].set_xlabel('Number of Threads')
    axes[0, 1].set_ylabel('Efficiency (%)')
    axes[0, 1].set_title(f'Efficiency vs Thread Count (Size: {representative_size:,})')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend()
    
    # Configure Plot 3
    axes[1, 0].set_xlabel('Number of Threads')
    axes[1, 0].set_ylabel('Time (ms)')
    axes[1, 0].set_title(f'Execution Time vs Thread Count (Size: {representative_size:,})')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend()
    
    # Configure Plot 4
    axes[1, 1].set_xlabel('Input Size (characters)')
    axes[1, 1].set_ylabel('Time (ms)')
    axes[1, 1].set_title('Execution Time vs Input Size (Max Threads)')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/parallel_benchmark_plots.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/parallel_benchmark_plots.pdf", bbox_inches='tight')
    print(f"\nPlots saved to {output_dir}/parallel_benchmark_plots.png and .pdf")


def save_results(results: Dict[str, Any], output_dir: str):
    """Save benchmark results to files."""
    # Save raw results as JSON
    results_file = Path(output_dir) / "benchmark_results.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save as pickle for easy loading in Python
    results_pickle = Path(output_dir) / "benchmark_results.pkl"
    with open(results_pickle, 'wb') as f:
        pickle.dump(results, f)
    
    print(f"\nResults saved to:")
    print(f"  Raw data (JSON): {results_file}")
    print(f"  Raw data (pickle): {results_pickle}")


def print_summary_table(results: Dict[str, Any]):
    """Print a summary table of the benchmark results."""
    print("\n" + "="*100)
    print("PARALLEL BENCHMARK SUMMARY")
    print("="*100)
    
    for func_name, data in results.items():
        if 'results' not in data or not data['results']:
            continue
            
        print(f"\n{func_name.upper().replace('_', ' ')}")
        print("-" * 80)
        
        # Get a representative size
        sizes = sorted(data['results'].keys())
        for size in sizes:
            print(f"\nSize: {size:,} characters")
            print(f"{'Threads':<10} {'Time (ms)':<12} {'Memory (MB)':<12} {'Speedup':<10} {'Efficiency (%)':<15}")
            print("-" * 80)
            
            size_results = data['results'][size]
            thread_counts = sorted(size_results.keys())
            
            if not thread_counts:
                continue
            
            # Get baseline time (single thread or minimum threads)
            baseline_time = size_results[min(thread_counts)]['mean_time']
            
            for num_threads in thread_counts:
                result = size_results[num_threads]
                time_ms = result['mean_time'] * 1000
                memory_mb = result['mean_memory_mb']
                speedup = baseline_time / result['mean_time']
                efficiency = (speedup / num_threads) * 100
                
                print(f"{num_threads:<10} {time_ms:<12.2f} {memory_mb:<12.2f} {speedup:<10.2f} {efficiency:<15.1f}")


def main():
    """Main function to run the parallel factorization benchmark suite."""
    parser = argparse.ArgumentParser(description="Benchmark parallel factorization functions")
    parser.add_argument("--min-size", type=int, default=10000, 
                       help="Minimum input size in characters (default: 10000)")
    parser.add_argument("--max-size", type=int, default=1000000,
                       help="Maximum input size in characters (default: 1000000)")
    parser.add_argument("--num-sizes", type=int, default=5,
                       help="Number of different sizes to test (default: 5)")
    parser.add_argument("--max-threads", type=int, default=None,
                       help="Maximum number of threads to test (default: CPU count)")
    parser.add_argument("--runs", type=int, default=3,
                       help="Number of runs per benchmark (default: 3)")
    parser.add_argument("--output-dir", default="benchmarks/parallel_results",
                       help="Output directory for results (default: benchmarks/parallel_results)")
    parser.add_argument("--custom-sizes", nargs="+", type=int,
                       help="Custom list of sizes to benchmark")
    parser.add_argument("--custom-threads", nargs="+", type=int,
                       help="Custom list of thread counts to test")
    
    args = parser.parse_args()
    
    # Generate size list
    if args.custom_sizes:
        sizes = args.custom_sizes
    else:
        # Generate logarithmically spaced sizes
        sizes = np.logspace(np.log10(args.min_size), np.log10(args.max_size), args.num_sizes)
        sizes = [int(size) for size in sizes]
    
    # Generate thread count list
    if args.custom_threads:
        thread_counts = args.custom_threads
    else:
        max_threads = args.max_threads or multiprocessing.cpu_count()
        # Test 1, 2, 4, ..., up to max_threads
        thread_counts = [1]
        current = 2
        while current <= max_threads:
            thread_counts.append(current)
            current *= 2
        # Add max_threads if it's not already in the list
        if max_threads not in thread_counts:
            thread_counts.append(max_threads)
        thread_counts.sort()
    
    print("Parallel Factorization Function Benchmark Suite")
    print("="*50)
    print(f"Size range: {min(sizes):,} to {max(sizes):,} characters")
    print(f"Number of sizes: {len(sizes)}")
    print(f"Thread counts: {thread_counts}")
    print(f"Runs per benchmark: {args.runs}")
    print(f"Output directory: {args.output_dir}")
    
    # Run benchmarks
    results = run_parallel_benchmark_suite(sizes, thread_counts, args.runs, args.output_dir)
    
    # Create plots
    create_speedup_plots(results, args.output_dir)
    
    # Save results
    save_results(results, args.output_dir)
    
    # Print summary
    print_summary_table(results)
    
    print(f"\nBenchmark completed! Results saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
