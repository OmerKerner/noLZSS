#!/usr/bin/env python3
"""
Utility script for reading and using FASTA benchmark trend parameters.

This script provides functions to load trend parameters saved by fasta_benchmark.py
and make predictions for time, memory, and disk space usage based on input size.

Example usage:
    from benchmarks.fasta_predictor import load_trends, predict_resources
    
    trends = load_trends("benchmarks/fasta_results/trend_parameters.pkl")
    predictions = predict_resources(trends, input_size=100000)
    print(f"Expected time: {predictions['factorize_fasta_multiple_dna_w_rc_time']:.2f} ms")
"""

import json
import pickle
from pathlib import Path
from typing import Dict, Any, Optional
import numpy as np


def load_trends(trend_file: str) -> Dict[str, Any]:
    """
    Load trend parameters from a file.
    
    Args:
        trend_file: Path to trend parameters file (.json or .pkl)
        
    Returns:
        Dictionary with trend parameters
    """
    trend_path = Path(trend_file)
    
    if not trend_path.exists():
        raise FileNotFoundError(f"Trend file not found: {trend_file}")
    
    if trend_path.suffix == '.pkl':
        with open(trend_path, 'rb') as f:
            return pickle.load(f)
    elif trend_path.suffix == '.json':
        with open(trend_path, 'r') as f:
            return json.load(f)
    else:
        raise ValueError(f"Unsupported file format: {trend_path.suffix}. Use .pkl or .json")


def predict_from_trend(size: float, trend_params: Dict[str, float]) -> float:
    """
    Predict value from trend line parameters.
    
    Args:
        size: Input size (nucleotides)
        trend_params: Trend parameters dictionary
        
    Returns:
        Predicted value
    """
    if trend_params.get('log_scale', False):
        log_size = np.log10(size)
        log_prediction = trend_params['slope'] * log_size + trend_params['intercept']
        return 10 ** log_prediction
    else:
        return trend_params['slope'] * size + trend_params['intercept']


def predict_resources(trends: Dict[str, Any], input_size: int, 
                     functions: Optional[list] = None) -> Dict[str, float]:
    """
    Predict resource usage for FASTA functions based on input size.
    
    Args:
        trends: Trend parameters loaded from file
        input_size: Input size in nucleotides
        functions: List of function names to predict for (default: all)
        
    Returns:
        Dictionary with predictions for each function and metric
    """
    if functions is None:
        functions = [
            'factorize_fasta_multiple_dna_w_rc',
            'factorize_fasta_multiple_dna_no_rc',
            'write_factors_binary_file_fasta_multiple_dna_w_rc',
            'write_factors_binary_file_fasta_multiple_dna_no_rc'
        ]
    
    predictions = {}
    
    for func_name in functions:
        # Predict time
        time_key = f"{func_name}_time"
        if time_key in trends:
            time_ms = predict_from_trend(input_size, trends[time_key])
            predictions[time_key] = time_ms
            predictions[f"{func_name}_time_seconds"] = time_ms / 1000
        
        # Predict memory
        memory_key = f"{func_name}_memory"
        if memory_key in trends:
            memory_mb = predict_from_trend(input_size, trends[memory_key])
            predictions[memory_key] = memory_mb
            predictions[f"{func_name}_memory_gb"] = memory_mb / 1024
        
        # Predict disk space (for binary functions only)
        if 'write_factors_binary_file' in func_name:
            disk_key = f"{func_name}_disk_space"
            if disk_key in trends:
                disk_mb = predict_from_trend(input_size, trends[disk_key])
                predictions[disk_key] = disk_mb
                predictions[f"{func_name}_disk_space_gb"] = disk_mb / 1024
    
    return predictions


def estimate_cluster_resources(trends: Dict[str, Any], input_size: int, 
                              function_name: str, safety_factor: float = 1.5) -> Dict[str, Any]:
    """
    Estimate cluster resource requirements with safety factors.
    
    Args:
        trends: Trend parameters loaded from file
        input_size: Input size in nucleotides
        function_name: Function to estimate for
        safety_factor: Safety factor to apply to estimates (default: 1.5)
        
    Returns:
        Dictionary with cluster resource estimates
    """
    predictions = predict_resources(trends, input_size, [function_name])
    
    time_key = f"{function_name}_time_seconds"
    memory_key = f"{function_name}_memory_gb"
    
    estimate = {
        'input_size_nucleotides': input_size,
        'input_size_kbp': input_size / 1000,
        'input_size_mbp': input_size / 1000000,
        'function': function_name,
        'safety_factor': safety_factor
    }
    
    if time_key in predictions:
        base_time = predictions[time_key]
        estimate['estimated_time_seconds'] = base_time
        estimate['estimated_time_minutes'] = base_time / 60
        estimate['estimated_time_hours'] = base_time / 3600
        estimate['safe_time_seconds'] = base_time * safety_factor
        estimate['safe_time_minutes'] = (base_time * safety_factor) / 60
        estimate['safe_time_hours'] = (base_time * safety_factor) / 3600
    
    if memory_key in predictions:
        base_memory = predictions[memory_key]
        estimate['estimated_memory_gb'] = base_memory
        estimate['safe_memory_gb'] = base_memory * safety_factor
        # Round up to common cluster memory allocations
        safe_mem_rounded = round_up_to_power_of_2(base_memory * safety_factor)
        estimate['cluster_memory_gb'] = safe_mem_rounded
    
    # Add disk space for binary functions
    if 'write_factors_binary_file' in function_name:
        disk_key = f"{function_name}_disk_space_gb"
        if disk_key in predictions:
            base_disk = predictions[disk_key]
            estimate['estimated_disk_gb'] = base_disk
            estimate['safe_disk_gb'] = base_disk * safety_factor
    
    return estimate


def round_up_to_power_of_2(value: float) -> int:
    """Round up to the next power of 2 (common for cluster memory allocations)."""
    if value <= 1:
        return 1
    return int(2 ** np.ceil(np.log2(value)))


def generate_resource_table(trends: Dict[str, Any], sizes: list, 
                           function_name: str, output_file: Optional[str] = None) -> str:
    """
    Generate a resource estimation table for different input sizes.
    
    Args:
        trends: Trend parameters loaded from file
        sizes: List of input sizes to estimate for
        function_name: Function to estimate for
        output_file: Optional file to save the table to
        
    Returns:
        Formatted table as string
    """
    header = f"Resource Estimates for {function_name.replace('_', ' ').title()}"
    table = [header, "=" * len(header), ""]
    
    table.append(f"{'Size (kbp)':<12} {'Time (s)':<10} {'Time (min)':<12} {'Memory (GB)':<12} {'Cluster Mem (GB)':<16}")
    if 'write_factors_binary_file' in function_name:
        table[-1] += f" {'Disk (GB)':<12}"
    table.append("-" * (12 + 10 + 12 + 12 + 16 + (12 if 'write_factors_binary_file' in function_name else 0)))
    
    for size in sizes:
        estimate = estimate_cluster_resources(trends, size, function_name)
        
        size_kbp = size / 1000
        time_s = estimate.get('safe_time_seconds', 0)
        time_min = estimate.get('safe_time_minutes', 0)
        memory_gb = estimate.get('estimated_memory_gb', 0)
        cluster_mem = estimate.get('cluster_memory_gb', 0)
        
        row = f"{size_kbp:<12.1f} {time_s:<10.2f} {time_min:<12.2f} {memory_gb:<12.3f} {cluster_mem:<16d}"
        
        if 'write_factors_binary_file' in function_name:
            disk_gb = estimate.get('safe_disk_gb', 0)
            row += f" {disk_gb:<12.3f}"
        
        table.append(row)
    
    table_str = "\n".join(table)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(table_str)
        print(f"Table saved to {output_file}")
    
    return table_str


def main():
    """Example usage of the predictor functions."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Predict resource usage from FASTA benchmark trends")
    parser.add_argument("trend_file", help="Path to trend parameters file (.pkl or .json)")
    parser.add_argument("--size", type=int, help="Input size in nucleotides")
    parser.add_argument("--function", help="Function name to predict for (default: all functions)")
    parser.add_argument("--safety-factor", type=float, default=1.5, help="Safety factor for cluster resources")
    parser.add_argument("--table", nargs="+", type=int, 
                       help="Generate resource table for multiple sizes")
    parser.add_argument("--output", help="Save table to file")
    
    args = parser.parse_args()
    
    # Load trends
    trends = load_trends(args.trend_file)
    print(f"Loaded trends from {args.trend_file}")
    
    if args.table:
        # Generate table for multiple sizes
        if not args.function:
            print("Error: --function is required when using --table")
            return
        
        table = generate_resource_table(trends, args.table, args.function, args.output)
        print("\n" + table)
    else:
        # Single prediction
        if not args.size:
            print("Error: --size is required for single predictions")
            return
            
        if args.function:
            functions = [args.function]
        else:
            functions = None
        
        predictions = predict_resources(trends, args.size, functions)
        
        print(f"\nResource predictions for input size: {args.size:,} nucleotides ({args.size/1000:.1f} kbp)")
        print("-" * 80)
        
        for key, value in predictions.items():
            if 'time' in key and 'seconds' not in key:
                print(f"{key:<50}: {value:.2f} ms")
            elif 'memory' in key and 'gb' not in key:
                print(f"{key:<50}: {value:.3f} MB")
            elif 'disk_space' in key and 'gb' not in key:
                print(f"{key:<50}: {value:.3f} MB")
        
        # Show cluster estimates for each function
        function_names = [
            'factorize_fasta_multiple_dna_w_rc',
            'factorize_fasta_multiple_dna_no_rc',
            'write_factors_binary_file_fasta_multiple_dna_w_rc',
            'write_factors_binary_file_fasta_multiple_dna_no_rc'
        ]
        
        print(f"\nCluster Resource Estimates (safety factor: {args.safety_factor})")
        print("=" * 80)
        
        for func_name in function_names:
            if args.function and func_name != args.function:
                continue
                
            estimate = estimate_cluster_resources(trends, args.size, func_name, args.safety_factor)
            
            print(f"\n{func_name.replace('_', ' ').title()}:")
            print(f"  Estimated time: {estimate.get('estimated_time_minutes', 0):.2f} minutes")
            print(f"  Safe time: {estimate.get('safe_time_minutes', 0):.2f} minutes")
            print(f"  Estimated memory: {estimate.get('estimated_memory_gb', 0):.3f} GB")
            print(f"  Cluster memory: {estimate.get('cluster_memory_gb', 0)} GB")
            
            if 'estimated_disk_gb' in estimate:
                print(f"  Estimated disk: {estimate.get('estimated_disk_gb', 0):.3f} GB")
                print(f"  Safe disk: {estimate.get('safe_disk_gb', 0):.3f} GB")


if __name__ == "__main__":
    main()