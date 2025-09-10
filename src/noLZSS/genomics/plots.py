"""
FASTA file plotting utilities.

This module provides functions for creating plots and visualizations
from FASTA files and their factorizations.
"""

from typing import Union, Optional, Dict, Any
from pathlib import Path

from ..utils import NoLZSSError
from .fasta import _parse_fasta_content
from .sequences import detect_sequence_type


class PlotError(NoLZSSError):
    """Raised when plotting operations fail."""
    pass


def plot_single_seq_accum_factors_from_fasta(
    fasta_filepath: Union[str, Path],
    output_dir: Union[str, Path],
    max_sequences: Optional[int] = None,
    save_factors_text: bool = True,
    save_factors_binary: bool = False
) -> Dict[str, Dict[str, Any]]:
    """
    Process a FASTA file, factorize all sequences, create plots, and save results.

    For each sequence in the FASTA file:
    - Factorizes the sequence
    - Saves factor data (text and/or binary format)
    - Creates and saves a plot of factor lengths

    Args:
        fasta_filepath: Path to input FASTA file
        output_dir: Directory to save all output files
        max_sequences: Maximum number of sequences to process (None for all)
        save_factors_text: Whether to save factors as text files
        save_factors_binary: Whether to save factors as binary files

    Returns:
        Dictionary with processing results for each sequence:
        {
            'sequence_id': {
                'sequence_length': int,
                'num_factors': int,
                'factors_file': str,  # path to saved factors
                'plot_file': str,     # path to saved plot
                'factors': List[Tuple[int, int, int]]  # the factors
            }
        }

    Raises:
        PlotError: If FASTA processing fails
        FileNotFoundError: If input file doesn't exist
    """
    from ..core import factorize, write_factors_binary_file
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend for batch processing
    import re

    fasta_filepath = Path(fasta_filepath)
    output_dir = Path(output_dir)

    if not fasta_filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_filepath}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read FASTA file
    sequences = _parse_fasta_content(fasta_filepath.read_text())

    if not sequences:
        raise PlotError("No sequences found in FASTA file")

    results = {}
    processed_count = 0

    for seq_id, sequence in sequences.items():
        if max_sequences is not None and processed_count >= max_sequences:
            break

        print(f"Processing sequence {seq_id} ({len(sequence)} bp)...")

        # Detect sequence type and validate
        seq_type = detect_sequence_type(sequence)

        if seq_type == 'dna':
            # Validate as nucleotide
            if not re.match(r'^[ACGT]+$', sequence.upper()):
                invalid_chars = set(sequence.upper()) - set('ACGT')
                print(f"  Warning: Skipping {seq_id} - contains invalid nucleotides: {invalid_chars}")
                continue
            sequence = sequence.upper()
            print("  Detected nucleotide sequence")

        elif seq_type == 'protein':
            # Validate as amino acid
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
            if not all(c in valid_aa for c in sequence.upper()):
                invalid_chars = set(sequence.upper()) - valid_aa
                print(f"  Warning: Skipping {seq_id} - contains invalid amino acids: {invalid_chars}")
                continue
            sequence = sequence.upper()
            print("  Detected amino acid sequence")

        else:
            print(f"  Warning: Skipping {seq_id} - unknown sequence type: {seq_type}")
            continue

        # Factorize
        try:
            factors = factorize(sequence.encode('ascii'))
            print(f"  Factorized into {len(factors)} factors")
        except Exception as e:
            print(f"  Warning: Failed to factorize {seq_id}: {e}")
            continue

        # Save factors as text
        factors_text_file = None
        if save_factors_text:
            factors_text_file = output_dir / f"factors_{seq_id}.txt"
            try:
                with open(factors_text_file, 'w') as f:
                    f.write(f"Sequence: {seq_id}\n")
                    f.write(f"Length: {len(sequence)}\n")
                    f.write(f"Number of factors: {len(factors)}\n")
                    f.write("Factors (position, length, reference):\n")
                    for i, (pos, length, ref) in enumerate(factors):
                        f.write(f"{i+1:4d}: ({pos:6d}, {length:4d}, {ref:6d})\n")
                print(f"  Saved factors to {factors_text_file}")
            except Exception as e:
                print(f"  Warning: Failed to save text factors for {seq_id}: {e}")

        # Save factors as binary
        factors_binary_file = None
        if save_factors_binary:
            factors_binary_file = output_dir / f"factors_{seq_id}.bin"
            try:
                # Create a temporary file with just this sequence
                temp_fasta = output_dir / f"temp_{seq_id}.fasta"
                with open(temp_fasta, 'w') as f:
                    f.write(f">{seq_id}\n{sequence}\n")

                write_factors_binary_file(str(temp_fasta), str(factors_binary_file))
                temp_fasta.unlink()  # Clean up temp file
                print(f"  Saved binary factors to {factors_binary_file}")
            except Exception as e:
                print(f"  Warning: Failed to save binary factors for {seq_id}: {e}")

        # Create plot
        plot_file = output_dir / f"plot_{seq_id}.png"
        try:
            from ..utils import plot_factor_lengths
            plot_factor_lengths(factors, save_path=plot_file, show_plot=False)
            print(f"  Saved plot to {plot_file}")
        except Exception as e:
            print(f"  Warning: Failed to create plot for {seq_id}: {e}")
            plot_file = None

        # Store results
        results[seq_id] = {
            'sequence_length': len(sequence),
            'num_factors': len(factors),
            'factors_file': str(factors_text_file) if factors_text_file else None,
            'binary_file': str(factors_binary_file) if factors_binary_file else None,
            'plot_file': str(plot_file) if plot_file else None,
            'factors': factors
        }

        processed_count += 1

    print(f"\nProcessed {len(results)} sequences successfully")
    return results