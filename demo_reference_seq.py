#!/usr/bin/env python3
"""
Demonstration of factorize_w_reference_seq functionality.

This script shows how to use the new reference sequence factorization
functions to factorize a target sequence using a reference sequence.
"""

import sys
import os
from pathlib import Path

# Add build directory to path for demonstration
build_path = Path(__file__).parent / "build"
if build_path.exists():
    sys.path.insert(0, str(build_path))
    try:
        import _noLZSS
        print("Using built extension from build directory")
    except ImportError:
        print("Could not import from build directory, trying installed package")
        from noLZSS._noLZSS import factorize_w_reference_seq, factorize_w_reference_seq_file
        _noLZSS = type('Module', (), {
            'factorize_w_reference_seq': factorize_w_reference_seq,
            'factorize_w_reference_seq_file': factorize_w_reference_seq_file
        })()

def demonstrate_basic_usage():
    """Demonstrate basic usage of reference sequence factorization."""
    print("=== Basic Usage Example ===")
    
    # Example: Human chromosome reference and a target sequence
    reference = "ATCGATCGATCGATCGATCG"  # Reference sequence
    target = "GATCGATCGATCGA"          # Target sequence to factorize
    
    print(f"Reference sequence: {reference}")
    print(f"Target sequence:    {target}")
    print()
    
    # Factorize target using reference
    factors = _noLZSS.factorize_w_reference_seq(reference, target)
    
    print(f"Factorization results: {len(factors)} factors")
    print("Factor format: (start_in_target, length, reference_position, is_reverse_complement)")
    print()
    
    for i, (start, length, ref, is_rc) in enumerate(factors):
        print(f"Factor {i+1}: start={start:2d}, length={length}, ref={ref:2d}, RC={is_rc}")
        
        # Show what part of target this factor covers
        target_part = target[start:start+length]
        
        if is_rc:
            print(f"  Target part: '{target_part}' (matches reverse complement in reference)")
        else:
            print(f"  Target part: '{target_part}' (matches forward pattern in reference)")
    
    print()

def demonstrate_file_output():
    """Demonstrate writing factors to a binary file."""
    print("=== File Output Example ===")
    
    reference = "ATCGATCGATCGATCGATCGATCGATCG"
    target = "GATCGATCGATCGATCGA"
    output_file = "/tmp/demo_factors.bin"
    
    print(f"Writing factors to: {output_file}")
    
    num_factors = _noLZSS.factorize_w_reference_seq_file(reference, target, output_file)
    
    if os.path.exists(output_file):
        file_size = os.path.getsize(output_file)
        print(f"✓ Successfully wrote {num_factors} factors ({file_size} bytes)")
        
        # Clean up
        os.unlink(output_file)
    else:
        print("✗ Failed to create output file")
    
    print()

def demonstrate_reverse_complement():
    """Demonstrate reverse complement matching."""
    print("=== Reverse Complement Example ===")
    
    # Create a case where the target matches the reverse complement of reference
    reference = "ATCGATCG"
    target = "CGATCGAT"  # This is the reverse complement of the reference
    
    print(f"Reference:         {reference}")
    print(f"Target:            {target}")
    print(f"Reference RC:      {''.join({'A':'T','T':'A','C':'G','G':'C'}[c] for c in reference[::-1])}")
    print()
    
    factors = _noLZSS.factorize_w_reference_seq(reference, target)
    
    for i, (start, length, ref, is_rc) in enumerate(factors):
        if is_rc:
            print(f"Factor {i+1}: Found reverse complement match!")
            print(f"  Position {start}-{start+length-1} in target matches RC pattern at ref position {ref}")
        else:
            print(f"Factor {i+1}: Forward match at position {start}-{start+length-1}")
    
    print()

def demonstrate_practical_example():
    """Demonstrate a more practical genomics example."""
    print("=== Practical Genomics Example ===")
    
    # Simulate a reference genome segment and a read with sequencing errors
    reference = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    target = "GATCGATCGATCGATCGA"  # A read that maps to part of the reference
    
    print("Scenario: Mapping a sequencing read to a reference genome")
    print(f"Reference genome: {reference}")
    print(f"Sequencing read:  {target}")
    print()
    
    factors = _noLZSS.factorize_w_reference_seq(reference, target)
    
    print("Factorization helps identify:")
    print("1. Which parts of the read match the reference")
    print("2. Whether matches are on forward or reverse strand")
    print("3. Compression ratio relative to reference")
    print()
    
    total_matched_length = sum(length for _, length, _, _ in factors)
    compression_ratio = total_matched_length / len(target)
    
    print(f"Read length: {len(target)} bp")
    print(f"Total matched: {total_matched_length} bp")
    print(f"Match ratio: {compression_ratio:.2%}")
    
    if compression_ratio > 0.8:
        print("✓ Good alignment to reference")
    else:
        print("⚠ Poor alignment - may contain many mutations or sequencing errors")
    
    print()

def main():
    """Run all demonstrations."""
    print("Reference Sequence Factorization Demonstration")
    print("=" * 50)
    print()
    
    try:
        demonstrate_basic_usage()
        demonstrate_file_output()
        demonstrate_reverse_complement()
        demonstrate_practical_example()
        
        print("Demonstration completed successfully!")
        print()
        print("Key benefits of reference sequence factorization:")
        print("• More efficient compression when target shares patterns with reference")
        print("• Identifies both forward and reverse complement matches")
        print("• Useful for genomics applications like read mapping and variant calling")
        print("• Factors start positions are relative to target sequence")
        
    except Exception as e:
        print(f"Error during demonstration: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())