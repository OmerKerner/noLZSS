"""
Example demonstrating per-sequence FASTA factorization.

This example shows how to use the new per-sequence factorization functions
that factorize each sequence in a FASTA file independently without sentinel
concatenation.
"""

import tempfile
import os

# Import the new functions
from noLZSS._noLZSS import (
    factorize_fasta_dna_w_rc_per_sequence,
    factorize_fasta_dna_no_rc_per_sequence,
    count_factors_fasta_dna_w_rc_per_sequence,
    write_factors_binary_file_fasta_dna_w_rc_per_sequence,
    parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence,
)

def main():
    # Create a sample FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        fasta_path = f.name
        f.write(">sequence1\n")
        f.write("ATCGATCGATCG\n")
        f.write(">sequence2\n")
        f.write("GCTAGCTAGCTA\n")
        f.write(">sequence3\n")
        f.write("TTAATTAATTAA\n")
    
    try:
        print("="*70)
        print("Per-Sequence FASTA Factorization Example")
        print("="*70)
        
        # 1. Basic factorization with reverse complement
        print("\n1. Factorize each sequence separately (with RC):")
        per_seq_factors, sequence_ids = factorize_fasta_dna_w_rc_per_sequence(fasta_path)
        
        for i, (seq_id, factors) in enumerate(zip(sequence_ids, per_seq_factors)):
            print(f"   Sequence {i+1} '{seq_id}':")
            print(f"     - {len(factors)} factors")
            if len(factors) > 0:
                # Show first factor
                start, length, ref, is_rc = factors[0]
                rc_str = " (RC)" if is_rc else ""
                print(f"     - First factor: start={start}, len={length}, ref={ref}{rc_str}")
        
        # 2. Compare with no reverse complement
        print("\n2. Factorize without reverse complement:")
        per_seq_factors_no_rc, _ = factorize_fasta_dna_no_rc_per_sequence(fasta_path)
        
        for i, (seq_id, factors) in enumerate(zip(sequence_ids, per_seq_factors_no_rc)):
            print(f"   Sequence {i+1} '{seq_id}': {len(factors)} factors")
        
        print("\n   Note: RC typically produces fewer factors (better compression)")
        
        # 3. Count factors without storing them (memory efficient)
        print("\n3. Count factors (memory efficient):")
        total_count = count_factors_fasta_dna_w_rc_per_sequence(fasta_path)
        print(f"   Total factors across all sequences: {total_count}")
        
        # 4. Write to binary file
        print("\n4. Write factors to binary file:")
        with tempfile.NamedTemporaryFile(suffix='.bin', delete=False) as tmp:
            output_path = tmp.name
        
        count = write_factors_binary_file_fasta_dna_w_rc_per_sequence(fasta_path, output_path)
        file_size = os.path.getsize(output_path)
        print(f"   Wrote {count} factors to {output_path}")
        print(f"   File size: {file_size} bytes")
        os.remove(output_path)
        
        # 5. Parallel processing
        print("\n5. Parallel processing (4 threads):")
        count_parallel = parallel_write_factors_binary_file_fasta_dna_w_rc_per_sequence(
            fasta_path, output_path, num_threads=4
        )
        print(f"   Processed {count_parallel} factors using parallel processing")
        os.remove(output_path)
        
        # 6. Key differences from concatenated version
        print("\n6. Key Differences from Concatenated Version:")
        print("   ✓ No sentinel limitations - can handle unlimited sequences")
        print("   ✓ No cross-sequence matches - cleaner per-sequence results")
        print("   ✓ Better parallelization - distributes sequences across threads")
        print("   ✓ Easier interpretation - factors only reference their own sequence")
        
        print("\n" + "="*70)
        print("Example completed successfully!")
        print("="*70)
        
    finally:
        if os.path.exists(fasta_path):
            os.remove(fasta_path)


if __name__ == "__main__":
    main()
