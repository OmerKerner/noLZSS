#!/usr/bin/env python3
"""
Standalone LSF job script for noLZSS factorization.

This script is designed to be self-contained and can be used directly in LSF jobs
without requiring the noLZSS package to be installed on compute nodes.
"""

import sys
import argparse
from pathlib import Path


def main():
    """Execute factorization based on command-line arguments."""
    parser = argparse.ArgumentParser(description="Execute noLZSS factorization")
    parser.add_argument("input_file", type=Path, help="Input FASTA file")
    parser.add_argument("output_file", type=Path, help="Output binary file")
    parser.add_argument("--mode", choices=["with_rc", "no_rc"], 
                       default="with_rc", help="Factorization mode")
    
    args = parser.parse_args()
    
    # Import noLZSS (this requires it to be installed)
    try:
        from noLZSS._noLZSS import (
            write_factors_binary_file_fasta_multiple_dna_w_rc,
            write_factors_binary_file_fasta_multiple_dna_no_rc
        )
    except ImportError as e:
        print(f"Error: Could not import noLZSS: {e}", file=sys.stderr)
        print("Make sure noLZSS is installed in the Python environment", file=sys.stderr)
        sys.exit(1)
    
    # Verify input file exists
    if not args.input_file.exists():
        print(f"Error: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Execute factorization
    print(f"Input: {args.input_file}")
    print(f"Output: {args.output_file}")
    print(f"Mode: {args.mode}")
    
    try:
        if args.mode == "with_rc":
            num_factors = write_factors_binary_file_fasta_multiple_dna_w_rc(
                str(args.input_file), str(args.output_file)
            )
        else:
            num_factors = write_factors_binary_file_fasta_multiple_dna_no_rc(
                str(args.input_file), str(args.output_file)
            )
        
        print(f"Success! Generated {num_factors} factors")
        print(f"Output written to: {args.output_file}")
        
    except Exception as e:
        print(f"Error during factorization: {e}", file=sys.stderr)
        # Clean up partial output
        if args.output_file.exists():
            try:
                args.output_file.unlink()
            except:
                pass
        sys.exit(1)


if __name__ == "__main__":
    main()
