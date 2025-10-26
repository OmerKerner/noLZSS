#!/usr/bin/env python3
"""
Example script demonstrating LSF batch factorizer usage.

This example shows how to use the LSF batch factorizer to process
multiple FASTA files on an LSF cluster.
"""

import tempfile
from pathlib import Path

# Note: This is a demonstration script. To actually run it, you need:
# 1. LSF cluster with bsub/bjobs available
# 2. Real FASTA files to process
# 3. Benchmark trends file (optional but recommended)

def show_example_file_list_content():
    """
    Show example file list content for demonstration.
    
    In real usage, create a text file with one FASTA path per line.
    """
    # Example file paths (replace with your actual files)
    example_files = [
        "/path/to/genome1.fasta",
        "/path/to/genome2.fasta",
        "/path/to/genome3.fasta",
    ]
    
    # Example content for a file list
    example_content = "\n".join(example_files) + "\n"
    
    print("Example files.txt content:")
    print(example_content)
    print("\nTo create this file:")
    print("  echo '/path/to/genome1.fasta' > files.txt")
    print("  echo '/path/to/genome2.fasta' >> files.txt")
    print("  echo '/path/to/genome3.fasta' >> files.txt")


def example_basic_usage():
    """Basic usage example."""
    print("="*60)
    print("Example 1: Basic Usage")
    print("="*60)
    
    command = """
python -m noLZSS.genomics.lsf_batch_factorize \\
    --file-list files.txt \\
    --output-dir results \\
    --max-threads 16 \\
    --mode with_reverse_complement
    """
    
    print("Command:")
    print(command)
    print("\nThis will:")
    print("- Process all files listed in files.txt")
    print("- Use with_reverse_complement mode")
    print("- Allow up to 16 threads per job")
    print("- Automatically estimate resources")
    print()


def example_with_benchmark_trends():
    """Example with benchmark trends."""
    print("="*60)
    print("Example 2: With Benchmark Trends")
    print("="*60)
    
    print("\nStep 1: Generate benchmark trends")
    benchmark_command = """
python benchmarks/fasta_benchmark.py \\
    --output-dir benchmarks/fasta_results
    """
    print(benchmark_command)
    
    print("\nStep 2: Use trends for accurate estimation")
    command = """
python -m noLZSS.genomics.lsf_batch_factorize \\
    --file-list files.txt \\
    --output-dir results \\
    --max-threads 16 \\
    --trend-file benchmarks/fasta_results/trend_parameters.pkl \\
    --mode with_reverse_complement
    """
    print(command)
    print()


def example_with_custom_lsf_options():
    """Example with custom LSF options."""
    print("="*60)
    print("Example 3: Custom LSF Options")
    print("="*60)
    
    command = """
python -m noLZSS.genomics.lsf_batch_factorize \\
    --file-list files.txt \\
    --output-dir results \\
    --max-threads 12 \\
    --queue priority \\
    --bsub-args \"-P project_name\" \"-G group_name\" \\
    --safety-factor 2.0 \\
    --check-interval 120 \\
    --log-file batch_process.log
    """
    
    print("Command:")
    print(command)
    print("\nThis will:")
    print("- Submit to 'priority' queue")
    print("- Use project_name and group_name")
    print("- Apply 2x safety factor for resources")
    print("- Check job status every 2 minutes")
    print("- Log to batch_process.log")
    print()


def example_resource_estimation():
    """Example showing resource estimation."""
    print("="*60)
    print("Example 4: Understanding Resource Estimation")
    print("="*60)
    
    print("\nThe script automatically estimates resources based on file size:")
    print()
    print("File Size       | Threads | Memory (estimate)")
    print("----------------|---------|------------------")
    print("< 100 kb        | 1       | 2-4 GB")
    print("100 kb - 1 Mbp  | 2-4     | 4-8 GB")
    print("1 Mbp - 10 Mbp  | up to 8 | 8-16 GB")
    print("> 10 Mbp        | max     | 16-64 GB")
    print()
    print("Actual allocation depends on:")
    print("- Benchmark trends (if available)")
    print("- Safety factor (default 1.5x)")
    print("- Cluster memory increments (2, 4, 8, 16, 32, 64, 128 GB)")
    print()


def example_output_structure():
    """Example showing output structure."""
    print("="*60)
    print("Example 5: Output Structure")
    print("="*60)
    
    print("\nAfter processing, your output directory will contain:")
    print()
    print("results/")
    print("├── with_reverse_complement/")
    print("│   ├── genome1.bin           # Factorization results")
    print("│   ├── genome2.bin")
    print("│   └── genome3.bin")
    print("├── lsf_scripts/")
    print("│   ├── factorize_genome1.sh  # Job scripts")
    print("│   ├── factorize_genome2.sh")
    print("│   └── factorize_genome3.sh")
    print("├── lsf_logs/")
    print("│   ├── factorize_genome1.out # Job stdout")
    print("│   ├── factorize_genome1.err # Job stderr")
    print("│   ├── factorize_genome2.out")
    print("│   ├── factorize_genome2.err")
    print("│   ├── factorize_genome3.out")
    print("│   └── factorize_genome3.err")
    print("└── lsf_job_results.json      # Job status summary")
    print()


def example_check_results():
    """Example of checking results."""
    print("="*60)
    print("Example 6: Checking Results")
    print("="*60)
    
    print("\nView summary:")
    print("  cat results/lsf_job_results.json | jq")
    print()
    print("Check for failures:")
    print("  grep -A5 '\"failed\"' results/lsf_job_results.json")
    print()
    print("View error logs:")
    print("  less results/lsf_logs/factorize_genome1.err")
    print()
    print("Count successful outputs:")
    print("  ls results/with_reverse_complement/*.bin | wc -l")
    print()


def main():
    """Run all examples."""
    print()
    print("*" * 60)
    print("LSF Batch Factorizer - Usage Examples")
    print("*" * 60)
    print()
    
    example_basic_usage()
    example_with_benchmark_trends()
    example_with_custom_lsf_options()
    example_resource_estimation()
    example_output_structure()
    example_check_results()
    
    print("="*60)
    print("For more information, see:")
    print("  docs/LSF_BATCH_FACTORIZER.md")
    print("="*60)
    print()


if __name__ == "__main__":
    main()
