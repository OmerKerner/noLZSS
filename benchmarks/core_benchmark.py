#!/usr/bin/env python3
"""
Core factorization benchmark script for noLZSS.

This script benchmarks the core factorization functions with time and memory measurement
across different input sizes. Creates log-log plots with trend lines and saves
coefficients for use by other scripts.

Functions benchmarked:
- factorize() - Basic string factorization
- factorize_file() - File-based factorization
- count_factors() - Factor counting only
- count_factors_file() - File-based factor counting
- write_factors_binary_file() - Binary file output
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
    'factorize',
    'factorize_file',
    'count_factors',
    'count_factors_file',
    'write_factors_binary_file',
]

missing_functions = [func for func in REQUIRED_FUNCTIONS if not hasattr(cpp, func)]
if missing_functions:
    print(f"Error: C++ extension is missing required functions: {', '.join(missing_functions)}")
    print("This usually means the package needs to be rebuilt.")
    print("Run: pip install -e . --no-build-isolation --force-reinstall")
    exit(1)


def generate_text_content(size: int, alphabet_size: int = 4) -> str:
    """
    Generate random text content for benchmarking.
    
    Args:
        size: Size of the text in characters
        alphabet_size: Size of the alphabet (default: 4 for DNA-like)
        
    Returns:
        Random text string
    """
    if alphabet_size == 4:
        # DNA-like alphabet
        alphabet = 'ACGT'
    elif alphabet_size == 26:
        # English-like alphabet
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    else:
        # Custom alphabet
        alphabet = ''.join(chr(ord('A') + i) for i in range(alphabet_size))
    
    return ''.join(random.choices(alphabet, k=size))


def create_temp_file(content: str) -> str:
    """Create a temporary text file and return its path."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        f.write(content)
        return f.name


def benchmark_function(func, *args, runs: int = 3) -> Dict[str, Any]:
    """
    Benchmark a function that returns factors or a count.
    
    Args:
        func: Function to benchmark
        *args: Arguments to pass to the function
        runs: Number of runs to average
        
    Returns:
        Dictionary with timing and memory statistics
    """
    times = []
    memories = []
    result_values = []
    
    # Convert string arguments to bytes for C++ functions
    converted_args = []
    for arg in args:
        if isinstance(arg, str):
            converted_args.append(arg.encode('utf-8'))
        else:
            converted_args.append(arg)
    
    for _ in range(runs):
        tracemalloc.start()
        start_time = time.perf_counter()
        
        try:
            result = func(*converted_args)
            end_time = time.perf_counter()
            current, peak = tracemalloc.get_traced_memory()
            
            times.append(end_time - start_time)
            memories.append(peak)
            
            # Store result value (either factor count or number of factors)
            if isinstance(result, int):
                result_values.append(result)
            elif isinstance(result, list):
                result_values.append(len(result))
            else:
                result_values.append(0)
            
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
        'mean_result': stats.mean(result_values),
        'all_times': times,
        'all_memories': memories
    }


def benchmark_binary_function(func, input_path: str, runs: int = 3) -> Dict[str, Any]:
    """
    Benchmark a binary file writing function.
    
    Args:
        func: Function to benchmark
        input_path: Path to input text file
        runs: Number of runs to average
        
    Returns:
        Dictionary with timing, memory, and disk space statistics
    """
    times = []
    memories = []
    file_sizes = []
    
    for _ in range(runs):
        # Create temporary binary output file
        with tempfile.NamedTemporaryFile(delete=False) as f:
            binary_path = f.name
        
        try:
            tracemalloc.start()
            start_time = time.perf_counter()
            
            func(input_path, binary_path)
            
            end_time = time.perf_counter()
            current, peak = tracemalloc.get_traced_memory()
            
            # Get file size
            file_size = os.path.getsize(binary_path)
            
            times.append(end_time - start_time)
            memories.append(peak)
            file_sizes.append(file_size)
            
        except Exception as e:
            print(f"Error in binary benchmark: {e}")
            return None
        finally:
            tracemalloc.stop()
            # Clean up temporary file
            try:
                os.unlink(binary_path)
            except:
                pass
    
    return {
        'mean_time': stats.mean(times),
        'std_time': stats.stdev(times) if len(times) > 1 else 0,
        'min_time': min(times),
        'max_time': max(times),
        'mean_memory_mb': stats.mean(memories) / (1024 * 1024),
        'max_memory_mb': max(memories) / (1024 * 1024),
        'mean_file_size_mb': stats.mean(file_sizes) / (1024 * 1024),
        'all_times': times,
        'all_memories': memories,
        'all_file_sizes': file_sizes
    }


def fit_trend_line(x_data: List[float], y_data: List[float], log_scale: bool = True) -> Dict[str, float]:
    """
    Fit a trend line to the data.
    
    Args:
        x_data: X values (e.g., input sizes)
        y_data: Y values (e.g., times or memory)
        log_scale: Whether to fit in log space
        
    Returns:
        Dictionary with slope, intercept, r_squared, and prediction function parameters
    """
    # Filter out zero or negative values for log scale
    if log_scale:
        valid_indices = [i for i, y in enumerate(y_data) if y > 0]
        if len(valid_indices) < 2:
            # Not enough valid data points for log scale, fallback to linear
            log_scale = False
        else:
            x_data = [x_data[i] for i in valid_indices]
            y_data = [y_data[i] for i in valid_indices]
    
    if log_scale:
        # Convert to log space for fitting
        x_log = np.log10(x_data)
        y_log = np.log10(y_data)
        slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(x_log, y_log)
        r_squared = r_value ** 2
        
        return {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'p_value': p_value,
            'std_err': std_err,
            'log_scale': True,
            'equation': f"log10(y) = {slope:.3f} * log10(x) + {intercept:.3f}",
            'power_law': f"y = {10**intercept:.3e} * x^{slope:.3f}"
        }
    else:
        slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(x_data, y_data)
        r_squared = r_value ** 2
        
        return {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_squared,
            'p_value': p_value,
            'std_err': std_err,
            'log_scale': False,
            'equation': f"y = {slope:.3f} * x + {intercept:.3f}"
        }


def predict_from_trend(size: float, trend_params: Dict[str, float]) -> float:
    """Predict value from trend line parameters."""
    if trend_params['log_scale']:
        log_size = np.log10(size)
        log_prediction = trend_params['slope'] * log_size + trend_params['intercept']
        return 10 ** log_prediction
    else:
        return trend_params['slope'] * size + trend_params['intercept']


def run_benchmark_suite(sizes: List[int], runs: int = 3, output_dir: str = "benchmarks/core_results") -> Dict[str, Any]:
    """
    Run the complete benchmark suite for all core factorization functions.
    
    Args:
        sizes: List of input sizes in characters (e.g., [1000, 10000, 100000])
        runs: Number of runs per benchmark
        output_dir: Directory to save results
        
    Returns:
        Complete benchmark results
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Benchmark string-based functions
    string_functions = {
        'factorize': cpp.factorize,
        'count_factors': cpp.count_factors,
    }
    
    for func_name, func in string_functions.items():
        print(f"\n=== Benchmarking {func_name} ===")
        results[func_name] = {'sizes': sizes, 'results': []}
        
        for size in sizes:
            print(f"  Size: {size:,} characters")
            
            # Generate text content
            text_content = generate_text_content(size)
            
            try:
                # Run benchmark
                result = benchmark_function(func, text_content, runs=runs)
                if result:
                    result['input_size'] = size
                    results[func_name]['results'].append(result)
                    
                    # Print quick summary
                    time_ms = result['mean_time'] * 1000
                    mem_mb = result['mean_memory_mb']
                    print(f"    Time: {time_ms:.2f} ms, Memory: {mem_mb:.2f} MB")
                else:
                    print(f"    Failed to benchmark size {size}")
                    
            except Exception as e:
                print(f"    Error: {e}")
    
    # Benchmark file-based functions
    file_functions = {
        'factorize_file': cpp.factorize_file,
        'count_factors_file': cpp.count_factors_file,
    }
    
    for func_name, func in file_functions.items():
        print(f"\n=== Benchmarking {func_name} ===")
        results[func_name] = {'sizes': sizes, 'results': []}
        
        for size in sizes:
            print(f"  Size: {size:,} characters")
            
            # Generate text content and create file
            text_content = generate_text_content(size)
            text_path = create_temp_file(text_content)
            
            try:
                # Run benchmark with appropriate arguments
                if func_name == 'factorize_file':
                    result = benchmark_function(func, text_path, 0, runs=runs)
                else:
                    result = benchmark_function(func, text_path, runs=runs)
                
                if result:
                    result['input_size'] = size
                    results[func_name]['results'].append(result)
                    
                    # Print quick summary
                    time_ms = result['mean_time'] * 1000
                    mem_mb = result['mean_memory_mb']
                    print(f"    Time: {time_ms:.2f} ms, Memory: {mem_mb:.2f} MB")
                else:
                    print(f"    Failed to benchmark size {size}")
                    
            finally:
                # Clean up temporary file
                try:
                    os.unlink(text_path)
                except:
                    pass
    
    # Benchmark binary file writing function
    print(f"\n=== Benchmarking write_factors_binary_file ===")
    results['write_factors_binary_file'] = {'sizes': sizes, 'results': []}
    
    for size in sizes:
        print(f"  Size: {size:,} characters")
        
        # Generate text content and create file
        text_content = generate_text_content(size)
        text_path = create_temp_file(text_content)
        
        try:
            # Run benchmark
            result = benchmark_binary_function(cpp.write_factors_binary_file, text_path, runs=runs)
            
            if result:
                result['input_size'] = size
                results['write_factors_binary_file']['results'].append(result)
                
                # Print quick summary
                time_ms = result['mean_time'] * 1000
                mem_mb = result['mean_memory_mb']
                file_mb = result['mean_file_size_mb']
                print(f"    Time: {time_ms:.2f} ms, Memory: {mem_mb:.2f} MB, File: {file_mb:.2f} MB")
            else:
                print(f"    Failed to benchmark size {size}")
                
        finally:
            # Clean up temporary file
            try:
                os.unlink(text_path)
            except:
                pass
    
    return results


def create_plots_and_trends(results: Dict[str, Any], output_dir: str) -> Dict[str, Dict[str, Any]]:
    """
    Create plots and fit trend lines for all benchmark results.
    
    Args:
        results: Benchmark results from run_benchmark_suite
        output_dir: Directory to save plots and trends
        
    Returns:
        Dictionary with trend line parameters for each function and metric
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    trends = {}
    
    # Set up the plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Core Factorization Function Benchmarks', fontsize=16)
    
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    markers = ['o', 's', '^', 'D', 'v']
    
    # Plot 1: Time vs Size
    ax = axes[0, 0]
    for i, (func_name, data) in enumerate(results.items()):
        if not data['results']:
            continue
            
        sizes = [r['input_size'] for r in data['results']]
        times = [r['mean_time'] * 1000 for r in data['results']]  # Convert to ms
        
        ax.loglog(sizes, times, color=colors[i % len(colors)], 
                 marker=markers[i % len(markers)], label=func_name.replace('_', ' '), linewidth=2)
        
        # Fit trend line
        if len(sizes) >= 2:
            trend = fit_trend_line(sizes, times, log_scale=True)
            trends[f"{func_name}_time"] = trend
            
            # Plot trend line
            x_trend = np.logspace(np.log10(min(sizes)), np.log10(max(sizes)), 100)
            y_trend = [predict_from_trend(x, trend) for x in x_trend]
            ax.loglog(x_trend, y_trend, '--', color=colors[i % len(colors)], alpha=0.7)
    
    ax.set_xlabel('Input Size (characters)')
    ax.set_ylabel('Time (ms)')
    ax.set_title('Execution Time vs Input Size')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Plot 2: Memory vs Size
    ax = axes[0, 1]
    for i, (func_name, data) in enumerate(results.items()):
        if not data['results']:
            continue
            
        sizes = [r['input_size'] for r in data['results']]
        memories = [r['mean_memory_mb'] for r in data['results']]
        
        ax.loglog(sizes, memories, color=colors[i % len(colors)], 
                 marker=markers[i % len(markers)], label=func_name.replace('_', ' '), linewidth=2)
        
        # Fit trend line
        if len(sizes) >= 2:
            trend = fit_trend_line(sizes, memories, log_scale=True)
            trends[f"{func_name}_memory"] = trend
            
            # Plot trend line
            x_trend = np.logspace(np.log10(min(sizes)), np.log10(max(sizes)), 100)
            y_trend = [predict_from_trend(x, trend) for x in x_trend]
            ax.loglog(x_trend, y_trend, '--', color=colors[i % len(colors)], alpha=0.7)
    
    ax.set_xlabel('Input Size (characters)')
    ax.set_ylabel('Peak Memory (MB)')
    ax.set_title('Memory Usage vs Input Size')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Plot 3: Throughput vs Size
    ax = axes[1, 0]
    for i, (func_name, data) in enumerate(results.items()):
        if not data['results']:
            continue
            
        sizes = [r['input_size'] for r in data['results']]
        throughputs = [size / (r['mean_time'] * 1e6) for size, r in zip(sizes, data['results'])]  # MB/s
        
        ax.semilogx(sizes, throughputs, color=colors[i % len(colors)], 
                   marker=markers[i % len(markers)], label=func_name.replace('_', ' '), linewidth=2)
    
    ax.set_xlabel('Input Size (characters)')
    ax.set_ylabel('Throughput (MB/s)')
    ax.set_title('Throughput vs Input Size')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # Plot 4: File Size vs Input Size (for binary function)
    ax = axes[1, 1]
    if 'write_factors_binary_file' in results:
        data = results['write_factors_binary_file']
        if data['results']:
            sizes = [r['input_size'] for r in data['results']]
            file_sizes = [r['mean_file_size_mb'] for r in data['results']]
            
            ax.loglog(sizes, file_sizes, color='blue', marker='o', label='Binary file size', linewidth=2)
            
            # Fit trend line
            if len(sizes) >= 2:
                trend = fit_trend_line(sizes, file_sizes, log_scale=True)
                trends["write_factors_binary_file_disk_space"] = trend
                
                # Plot trend line
                x_trend = np.logspace(np.log10(min(sizes)), np.log10(max(sizes)), 100)
                y_trend = [predict_from_trend(x, trend) for x in x_trend]
                ax.loglog(x_trend, y_trend, '--', color='blue', alpha=0.7)
    
    ax.set_xlabel('Input Size (characters)')
    ax.set_ylabel('File Size (MB)')
    ax.set_title('Disk Space vs Input Size')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/core_benchmark_plots.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/core_benchmark_plots.pdf", bbox_inches='tight')
    print(f"\nPlots saved to {output_dir}/core_benchmark_plots.png and .pdf")
    
    return trends


def save_results_and_trends(results: Dict[str, Any], trends: Dict[str, Any], output_dir: str):
    """Save benchmark results and trend parameters to files."""
    # Save raw results as JSON
    results_file = Path(output_dir) / "benchmark_results.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save trends as JSON
    trends_file = Path(output_dir) / "trend_parameters.json"
    with open(trends_file, 'w') as f:
        json.dump(trends, f, indent=2)
    
    # Save trends as pickle for easy loading in Python
    trends_pickle = Path(output_dir) / "trend_parameters.pkl"
    with open(trends_pickle, 'wb') as f:
        pickle.dump(trends, f)
    
    print(f"\nResults saved to:")
    print(f"  Raw data: {results_file}")
    print(f"  Trends (JSON): {trends_file}")
    print(f"  Trends (pickle): {trends_pickle}")


def print_summary_table(results: Dict[str, Any], trends: Dict[str, Any]):
    """Print a summary table of the benchmark results."""
    print("\n" + "="*100)
    print("BENCHMARK SUMMARY")
    print("="*100)
    
    for func_name, data in results.items():
        if not data['results']:
            continue
            
        print(f"\n{func_name.upper().replace('_', ' ')}")
        print("-" * 80)
        print(f"{'Size (chars)':<15} {'Time (ms)':<12} {'Memory (MB)':<12} {'Throughput (MB/s)':<18}")
        print("-" * 80)
        
        for result in data['results']:
            size = result['input_size']
            time_ms = result['mean_time'] * 1000
            memory_mb = result['mean_memory_mb']
            throughput = result['input_size'] / (result['mean_time'] * 1e6)
            
            print(f"{size:<15,} {time_ms:<12.2f} {memory_mb:<12.2f} {throughput:<18.2f}")
    
    print("\n" + "="*100)
    print("TREND ANALYSIS")
    print("="*100)
    
    for trend_name, trend_params in trends.items():
        print(f"\n{trend_name.upper().replace('_', ' ')}")
        print(f"  Equation: {trend_params['equation']}")
        if 'power_law' in trend_params:
            print(f"  Power law: {trend_params['power_law']}")
        print(f"  R² = {trend_params['r_squared']:.4f}")
        print(f"  P-value = {trend_params['p_value']:.2e}")


def main():
    """Main function to run the core factorization benchmark suite."""
    parser = argparse.ArgumentParser(description="Benchmark core factorization functions with trend analysis")
    parser.add_argument("--min-size", type=int, default=1000, 
                       help="Minimum input size in characters (default: 1000)")
    parser.add_argument("--max-size", type=int, default=1000000,
                       help="Maximum input size in characters (default: 1000000)")
    parser.add_argument("--num-sizes", type=int, default=10,
                       help="Number of different sizes to test (default: 10)")
    parser.add_argument("--runs", type=int, default=3,
                       help="Number of runs per benchmark (default: 3)")
    parser.add_argument("--output-dir", default="benchmarks/core_results",
                       help="Output directory for results (default: benchmarks/core_results)")
    parser.add_argument("--custom-sizes", nargs="+", type=int,
                       help="Custom list of sizes to benchmark (overrides min-size, max-size, num-sizes)")
    
    args = parser.parse_args()
    
    # Generate size list
    if args.custom_sizes:
        sizes = args.custom_sizes
    else:
        # Generate logarithmically spaced sizes
        sizes = np.logspace(np.log10(args.min_size), np.log10(args.max_size), args.num_sizes)
        sizes = [int(size) for size in sizes]
    
    print("Core Factorization Function Benchmark Suite")
    print("="*50)
    print(f"Size range: {min(sizes):,} to {max(sizes):,} characters")
    print(f"Number of sizes: {len(sizes)}")
    print(f"Runs per benchmark: {args.runs}")
    print(f"Output directory: {args.output_dir}")
    print(f"Sizes to test: {[f'{s:,}' for s in sizes]}")
    
    # Run benchmarks
    results = run_benchmark_suite(sizes, args.runs, args.output_dir)
    
    # Create plots and fit trends
    trends = create_plots_and_trends(results, args.output_dir)
    
    # Save results
    save_results_and_trends(results, trends, args.output_dir)
    
    # Print summary
    print_summary_table(results, trends)
    
    print(f"\nBenchmark completed! Results saved to {args.output_dir}/")


if __name__ == "__main__":
    main()
