#!/usr/bin/env python3
"""
Unified benchmark runner for noLZSS.

This script runs all benchmark suites (core, DNA, FASTA, and parallel) with
configurable parameters and generates a comprehensive report.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List, Dict, Any
import json
import time


class BenchmarkRunner:
    """Unified benchmark runner for all noLZSS benchmark suites."""
    
    def __init__(self, output_dir: str = "benchmarks/all_results"):
        """
        Initialize the benchmark runner.
        
        Args:
            output_dir: Base output directory for all benchmark results
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.benchmarks_dir = Path(__file__).parent
        self.results = {}
        
    def run_core_benchmark(self, sizes: List[int], runs: int = 3) -> bool:
        """
        Run the core factorization benchmark.
        
        Args:
            sizes: List of input sizes to test
            runs: Number of runs per benchmark
            
        Returns:
            True if successful, False otherwise
        """
        print("\n" + "="*80)
        print("RUNNING CORE FACTORIZATION BENCHMARK")
        print("="*80)
        
        output_dir = self.output_dir / "core_results"
        
        cmd = [
            sys.executable,
            str(self.benchmarks_dir / "core_benchmark.py"),
            "--custom-sizes"] + [str(s) for s in sizes] + [
            "--runs", str(runs),
            "--output-dir", str(output_dir)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            self.results['core'] = {'status': 'success', 'output_dir': str(output_dir)}
            return True
        except subprocess.CalledProcessError as e:
            print(f"Core benchmark failed: {e}")
            self.results['core'] = {'status': 'failed', 'error': str(e)}
            return False
    
    def run_dna_benchmark(self, sizes: List[int], runs: int = 3) -> bool:
        """
        Run the DNA factorization benchmark.
        
        Args:
            sizes: List of input sizes to test
            runs: Number of runs per benchmark
            
        Returns:
            True if successful, False otherwise
        """
        print("\n" + "="*80)
        print("RUNNING DNA FACTORIZATION BENCHMARK")
        print("="*80)
        
        output_dir = self.output_dir / "dna_results"
        
        cmd = [
            sys.executable,
            str(self.benchmarks_dir / "dna_benchmark.py"),
            "--custom-sizes"] + [str(s) for s in sizes] + [
            "--runs", str(runs),
            "--output-dir", str(output_dir)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            self.results['dna'] = {'status': 'success', 'output_dir': str(output_dir)}
            return True
        except subprocess.CalledProcessError as e:
            print(f"DNA benchmark failed: {e}")
            self.results['dna'] = {'status': 'failed', 'error': str(e)}
            return False
    
    def run_fasta_benchmark(self, sizes: List[int], runs: int = 3) -> bool:
        """
        Run the FASTA factorization benchmark.
        
        Args:
            sizes: List of input sizes to test
            runs: Number of runs per benchmark
            
        Returns:
            True if successful, False otherwise
        """
        print("\n" + "="*80)
        print("RUNNING FASTA FACTORIZATION BENCHMARK")
        print("="*80)
        
        output_dir = self.output_dir / "fasta_results"
        
        cmd = [
            sys.executable,
            str(self.benchmarks_dir / "fasta_benchmark.py"),
            "--custom-sizes"] + [str(s) for s in sizes] + [
            "--runs", str(runs),
            "--output-dir", str(output_dir)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            self.results['fasta'] = {'status': 'success', 'output_dir': str(output_dir)}
            return True
        except subprocess.CalledProcessError as e:
            print(f"FASTA benchmark failed: {e}")
            self.results['fasta'] = {'status': 'failed', 'error': str(e)}
            return False
    
    def run_parallel_benchmark(self, sizes: List[int], threads: List[int], runs: int = 3) -> bool:
        """
        Run the parallel factorization benchmark.
        
        Args:
            sizes: List of input sizes to test
            threads: List of thread counts to test
            runs: Number of runs per benchmark
            
        Returns:
            True if successful, False otherwise
        """
        print("\n" + "="*80)
        print("RUNNING PARALLEL FACTORIZATION BENCHMARK")
        print("="*80)
        
        output_dir = self.output_dir / "parallel_results"
        
        cmd = [
            sys.executable,
            str(self.benchmarks_dir / "parallel_benchmark.py"),
            "--custom-sizes"] + [str(s) for s in sizes] + [
            "--custom-threads"] + [str(t) for t in threads] + [
            "--runs", str(runs),
            "--output-dir", str(output_dir)
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=False)
            self.results['parallel'] = {'status': 'success', 'output_dir': str(output_dir)}
            return True
        except subprocess.CalledProcessError as e:
            print(f"Parallel benchmark failed: {e}")
            self.results['parallel'] = {'status': 'failed', 'error': str(e)}
            return False
    
    def save_summary(self):
        """Save a summary of all benchmark results."""
        summary_file = self.output_dir / "benchmark_summary.json"
        
        with open(summary_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        print(f"\n\nBenchmark summary saved to: {summary_file}")
    
    def print_summary(self):
        """Print a summary of all benchmark results."""
        print("\n" + "="*80)
        print("BENCHMARK SUMMARY")
        print("="*80)
        
        for benchmark_name, result in self.results.items():
            status = result['status']
            if status == 'success':
                print(f"\n{benchmark_name.upper()}: ✓ Success")
                print(f"  Results: {result['output_dir']}")
            else:
                print(f"\n{benchmark_name.upper()}: ✗ Failed")
                if 'error' in result:
                    print(f"  Error: {result['error']}")
        
        success_count = sum(1 for r in self.results.values() if r['status'] == 'success')
        total_count = len(self.results)
        
        print(f"\n{success_count}/{total_count} benchmarks completed successfully")


def main():
    """Main entry point for the unified benchmark runner."""
    parser = argparse.ArgumentParser(
        description="Run all noLZSS benchmark suites",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all benchmarks with default settings
  python benchmarks/run_all_benchmarks.py --all
  
  # Run specific benchmarks with custom sizes
  python benchmarks/run_all_benchmarks.py --core --dna --custom-sizes 1000 10000 100000
  
  # Run parallel benchmark with specific thread counts
  python benchmarks/run_all_benchmarks.py --parallel --custom-threads 1 2 4 8
  
  # Run quick test (smaller sizes, fewer runs)
  python benchmarks/run_all_benchmarks.py --all --quick
        """
    )
    
    # Benchmark selection
    parser.add_argument("--all", action="store_true",
                       help="Run all benchmarks")
    parser.add_argument("--core", action="store_true",
                       help="Run core factorization benchmark")
    parser.add_argument("--dna", action="store_true",
                       help="Run DNA factorization benchmark")
    parser.add_argument("--fasta", action="store_true",
                       help="Run FASTA factorization benchmark")
    parser.add_argument("--parallel", action="store_true",
                       help="Run parallel factorization benchmark")
    
    # Size configuration
    parser.add_argument("--custom-sizes", nargs="+", type=int,
                       help="Custom list of sizes to benchmark")
    parser.add_argument("--quick", action="store_true",
                       help="Quick test with smaller sizes (1k, 10k, 100k)")
    parser.add_argument("--small", action="store_true",
                       help="Small benchmark (1k-100k, 8 sizes)")
    parser.add_argument("--medium", action="store_true",
                       help="Medium benchmark (1k-1M, 10 sizes) - default")
    parser.add_argument("--large", action="store_true",
                       help="Large benchmark (1k-10M, 12 sizes)")
    
    # Thread configuration for parallel benchmarks
    parser.add_argument("--custom-threads", nargs="+", type=int,
                       help="Custom list of thread counts for parallel benchmark")
    
    # Run configuration
    parser.add_argument("--runs", type=int, default=3,
                       help="Number of runs per benchmark (default: 3)")
    parser.add_argument("--output-dir", default="benchmarks/all_results",
                       help="Output directory for all results (default: benchmarks/all_results)")
    
    args = parser.parse_args()
    
    # Determine which benchmarks to run
    run_core = args.all or args.core
    run_dna = args.all or args.dna
    run_fasta = args.all or args.fasta
    run_parallel = args.all or args.parallel
    
    if not (run_core or run_dna or run_fasta or run_parallel):
        parser.error("At least one benchmark must be selected (--all, --core, --dna, --fasta, or --parallel)")
    
    # Determine sizes to test
    if args.custom_sizes:
        sizes = args.custom_sizes
    elif args.quick:
        sizes = [1000, 10000, 100000]
    elif args.small:
        sizes = [int(s) for s in [1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5]]
    elif args.large:
        sizes = [int(s) for s in [1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7]]
    else:  # medium (default)
        sizes = [int(s) for s in [1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6]]
    
    # Determine thread counts for parallel benchmark
    if args.custom_threads:
        threads = args.custom_threads
    else:
        threads = [1, 2, 4]
    
    print("noLZSS Unified Benchmark Runner")
    print("="*80)
    print(f"Benchmarks to run:")
    if run_core:
        print("  - Core factorization")
    if run_dna:
        print("  - DNA factorization")
    if run_fasta:
        print("  - FASTA factorization")
    if run_parallel:
        print("  - Parallel factorization")
    print(f"\nSizes: {[f'{s:,}' for s in sizes]}")
    if run_parallel:
        print(f"Thread counts: {threads}")
    print(f"Runs per benchmark: {args.runs}")
    print(f"Output directory: {args.output_dir}")
    
    # Create runner
    runner = BenchmarkRunner(args.output_dir)
    
    # Track overall timing
    start_time = time.time()
    
    # Run selected benchmarks
    if run_core:
        runner.run_core_benchmark(sizes, args.runs)
    
    if run_dna:
        runner.run_dna_benchmark(sizes, args.runs)
    
    if run_fasta:
        runner.run_fasta_benchmark(sizes, args.runs)
    
    if run_parallel:
        runner.run_parallel_benchmark(sizes, threads, args.runs)
    
    # Calculate total time
    total_time = time.time() - start_time
    minutes = int(total_time // 60)
    seconds = int(total_time % 60)
    
    # Save and print summary
    runner.save_summary()
    runner.print_summary()
    
    print(f"\nTotal benchmark time: {minutes}m {seconds}s")
    
    # Exit with appropriate code
    if all(r['status'] == 'success' for r in runner.results.values()):
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
