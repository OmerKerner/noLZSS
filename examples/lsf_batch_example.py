#!/usr/bin/env python3
"""
Example demonstrating the LSF batch factorizer.

This script shows how to use the LSF batch factorizer for processing
multiple FASTA files on an LSF cluster.
"""

import tempfile
from pathlib import Path


def create_example_files(output_dir: Path):
    """Create example FASTA files for demonstration."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Constants for sequence generation
    MEDIUM_SEQUENCE_REPEATS = 30  # ~1 kbp
    LARGE_SEQUENCE_REPEATS = 300  # ~10 kbp
    LARGE_FILE_SCAFFOLD_COUNT = 10
    
    # Create sample FASTA files of different sizes
    files = []
    
    # Small file (~100 bp)
    small_file = output_dir / "small_genome.fasta"
    with open(small_file, 'w') as f:
        f.write(">chr1\n")
        f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
        f.write(">chr2\n")
        f.write("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
    files.append(small_file)
    
    # Medium file (~1 kbp)
    medium_file = output_dir / "medium_genome.fasta"
    with open(medium_file, 'w') as f:
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCG" * MEDIUM_SEQUENCE_REPEATS
        f.write(">chromosome1\n")
        f.write(sequence + "\n")
    files.append(medium_file)
    
    # Large file (~10 kbp)
    large_file = output_dir / "large_genome.fasta"
    with open(large_file, 'w') as f:
        sequence = "ATCGATCGATCGATCGATCGATCGATCGATCG" * LARGE_SEQUENCE_REPEATS
        for i in range(LARGE_FILE_SCAFFOLD_COUNT):
            f.write(f">scaffold_{i}\n")
            f.write(sequence + "\n")
    files.append(large_file)
    
    return files


def create_file_list(files, output_path: Path):
    """Create a file list for batch processing."""
    with open(output_path, 'w') as f:
        for file_path in files:
            f.write(str(file_path) + "\n")
    
    print(f"Created file list: {output_path}")
    print(f"Contains {len(files)} files:")
    for file_path in files:
        print(f"  - {file_path.name}")


def main():
    """Run the example."""
    print("LSF Batch Factorizer Example")
    print("=" * 60)
    
    # Create temporary directory for example
    with tempfile.TemporaryDirectory(prefix="lsf_example_") as temp_dir:
        temp_path = Path(temp_dir)
        
        # Create example FASTA files
        print("\n1. Creating example FASTA files...")
        files = create_example_files(temp_path / "input")
        
        # Create file list
        print("\n2. Creating file list...")
        file_list = temp_path / "file_list.txt"
        create_file_list(files, file_list)
        
        # Output directory
        output_dir = temp_path / "output"
        
        print("\n3. Example command to submit jobs:")
        print("-" * 60)
        print(f"python -m noLZSS.genomics.lsf_batch_factorize \\")
        print(f"    --file-list {file_list} \\")
        print(f"    --output-dir {output_dir} \\")
        print(f"    --mode with_reverse_complement \\")
        print(f"    --queue normal \\")
        print(f"    --dry-run")
        print("-" * 60)
        
        print("\n4. To actually submit (remove --dry-run):")
        print("-" * 60)
        print(f"python -m noLZSS.genomics.lsf_batch_factorize \\")
        print(f"    --file-list {file_list} \\")
        print(f"    --output-dir {output_dir} \\")
        print(f"    --mode with_reverse_complement \\")
        print(f"    --queue normal")
        print("-" * 60)
        
        print("\n5. Alternative: Submit and exit without waiting:")
        print("-" * 60)
        print(f"python -m noLZSS.genomics.lsf_batch_factorize \\")
        print(f"    --file-list {file_list} \\")
        print(f"    --output-dir {output_dir} \\")
        print(f"    --no-wait")
        print("-" * 60)
        
        print("\n6. With custom resource estimation:")
        print("-" * 60)
        print(f"python -m noLZSS.genomics.lsf_batch_factorize \\")
        print(f"    --file-list {file_list} \\")
        print(f"    --output-dir {output_dir} \\")
        print(f"    --trend-file benchmarks/fasta_results/trend_parameters.pkl \\")
        print(f"    --safety-factor 2.0 \\")
        print(f"    --max-threads 4")
        print("-" * 60)
        
        print("\n" + "=" * 60)
        print("Example complete!")
        print("\nFor more information, see: docs/LSF_BATCH_GUIDE.md")


if __name__ == "__main__":
    main()
