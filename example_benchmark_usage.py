#!/usr/bin/env python3
"""
Example script demonstrating how to use the FASTA benchmark system for cluster resource planning.

This script shows how to:
1. Load pre-computed benchmark trends
2. Estimate resources for different input sizes
3. Generate resource tables for cluster job submission
"""

import sys
from pathlib import Path

# Add benchmarks directory to path
sys.path.append(str(Path(__file__).parent / 'benchmarks'))

from fasta_predictor import load_trends, predict_resources, estimate_cluster_resources, generate_resource_table


def main():
    # Path to the benchmark trends (run fasta_benchmark.py first to generate these)
    trends_file = "benchmarks/fasta_results/trend_parameters.pkl"
    
    if not Path(trends_file).exists():
        print(f"Error: Benchmark trends not found at {trends_file}")
        print("Please run: python benchmarks/fasta_benchmark.py")
        return
    
    # Load the trends
    print("Loading benchmark trends...")
    trends = load_trends(trends_file)
    print(f"Loaded trends for {len(trends)} metrics")
    
    # Example 1: Estimate resources for a specific job
    print("\n" + "="*60)
    print("EXAMPLE 1: Single Job Resource Estimation")
    print("="*60)
    
    input_size = 2_500_000  # 2.5 Mbp
    function_name = "factorize_fasta_multiple_dna_w_rc"
    
    estimate = estimate_cluster_resources(trends, input_size, function_name, safety_factor=2.0)
    
    print(f"Input size: {input_size:,} nucleotides ({input_size/1e6:.1f} Mbp)")
    print(f"Function: {function_name}")
    print(f"Estimated time: {estimate['estimated_time_minutes']:.1f} minutes")
    print(f"Safe time (2x safety): {estimate['safe_time_minutes']:.1f} minutes")
    print(f"Estimated memory: {estimate['estimated_memory_gb']:.2f} GB")
    print(f"Cluster memory allocation: {estimate['cluster_memory_gb']} GB")
    
    # Example 2: Resource table for different sizes
    print("\n" + "="*60)
    print("EXAMPLE 2: Resource Table for Multiple Sizes")
    print("="*60)
    
    sizes = [100_000, 500_000, 1_000_000, 5_000_000, 10_000_000]  # 100kbp to 10Mbp
    
    table = generate_resource_table(trends, sizes, function_name)
    print(table)
    
    # Example 3: Compare all functions for a specific size
    print("\n" + "="*60)
    print("EXAMPLE 3: Function Comparison")
    print("="*60)
    
    comparison_size = 1_000_000  # 1 Mbp
    functions = [
        'factorize_fasta_multiple_dna_w_rc',
        'factorize_fasta_multiple_dna_no_rc', 
        'write_factors_binary_file_fasta_multiple_dna_w_rc',
        'write_factors_binary_file_fasta_multiple_dna_no_rc'
    ]
    
    print(f"Resource comparison for {comparison_size:,} nucleotides (1 Mbp):")
    print(f"{'Function':<50} {'Time (min)':<12} {'Memory (GB)':<12} {'Disk (GB)':<12}")
    print("-" * 86)
    
    for func in functions:
        est = estimate_cluster_resources(trends, comparison_size, func)
        time_min = est.get('estimated_time_minutes', 0)
        memory_gb = est.get('estimated_memory_gb', 0)
        disk_gb = est.get('estimated_disk_gb', 0)
        
        func_short = func.replace('factorize_fasta_multiple_dna_', '').replace('write_factors_binary_file_fasta_multiple_dna_', 'binary_')
        print(f"{func_short:<50} {time_min:<12.2f} {memory_gb:<12.3f} {disk_gb:<12.3f}")
    
    # Example 4: Scaling analysis
    print("\n" + "="*60)
    print("EXAMPLE 4: Scaling Analysis")
    print("="*60)
    
    # Show how resources scale with input size
    test_sizes = [10_000, 100_000, 1_000_000, 10_000_000]
    
    print(f"Scaling analysis for {function_name}:")
    print(f"{'Size':<15} {'Time (min)':<12} {'Memory (GB)':<12} {'Time/Mbp':<12}")
    print("-" * 51)
    
    for size in test_sizes:
        est = estimate_cluster_resources(trends, size, function_name)
        time_min = est.get('estimated_time_minutes', 0)
        memory_gb = est.get('estimated_memory_gb', 0)
        time_per_mbp = time_min / (size / 1e6)
        
        size_str = f"{size/1000:.0f}kbp" if size < 1e6 else f"{size/1e6:.0f}Mbp"
        print(f"{size_str:<15} {time_min:<12.2f} {memory_gb:<12.3f} {time_per_mbp:<12.2f}")
    
    print("\n" + "="*60)
    print("CLUSTER JOB SUBMISSION TIPS")
    print("="*60)
    print("1. Use the 'Safe time' estimates for job time limits")
    print("2. Use the 'Cluster memory' values for memory requests")
    print("3. For binary functions, reserve additional disk space")
    print("4. Consider using job arrays for multiple files")
    print("5. Monitor actual usage to refine safety factors")


if __name__ == "__main__":
    main()