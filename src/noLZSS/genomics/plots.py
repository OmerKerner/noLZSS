"""
FASTA file plotting utilities.

This module provides functions for creating plots and visualizations
from FASTA files and their factorizations.
"""

from typing import Union, Optional, Dict, Any
from pathlib import Path
import warnings

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


def plot_multiple_seq_self_weizmann_factor_plot_from_fasta(
    fasta_filepath: Union[str, Path],
    name: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True
) -> None:
    """
    Create a Weizmann factor plot for multiple DNA sequences from a FASTA file.

    This function reads a FASTA file containing multiple DNA sequences, factorizes them
    using the multiple DNA with reverse complement algorithm, and creates a specialized
    plot where each factor is represented as a line. The plot shows the relationship
    between factor positions and their reference positions.

    Args:
        fasta_filepath: Path to the FASTA file containing DNA sequences
        name: Optional name for the plot title (defaults to FASTA filename)
        save_path: Optional path to save the plot image
        show_plot: Whether to display the plot

    Raises:
        PlotError: If plotting fails or FASTA file cannot be processed
        FileNotFoundError: If FASTA file doesn't exist
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import numpy as np
    except ImportError:
        warnings.warn("matplotlib is required for plotting. Install with: pip install matplotlib", UserWarning)
        return

    from .._noLZSS import factorize_fasta_multiple_dna_w_rc

    fasta_filepath = Path(fasta_filepath)

    if not fasta_filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_filepath}")

    # Determine plot title
    if name is None:
        name = fasta_filepath.stem

    try:
        # Get factors from FASTA file
        print(f"Reading and factorizing sequences from {fasta_filepath}...")
        factors = factorize_fasta_multiple_dna_w_rc(str(fasta_filepath))

        print(f"Preparing plot for {len(factors)} factors...")
        if not factors:
            raise PlotError("No factors found in FASTA file")

        # Extract factor data
        positions = []
        lengths = []
        refs = []
        is_rcs = []

        for factor in factors:
            if len(factor) == 4:  # (start, length, ref, is_rc) tuple
                start, length, ref, is_rc = factor
            else:  # Assume (start, length, ref) format, default is_rc to False
                start, length, ref = factor
                is_rc = False

            positions.append(start)
            lengths.append(length)
            refs.append(ref)
            is_rcs.append(is_rc)

        # Convert to numpy arrays for easier processing
        positions = np.array(positions)
        lengths = np.array(lengths)
        refs = np.array(refs)
        is_rcs = np.array(is_rcs)

        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 12))

        # Calculate color intensities based on factor lengths
        # Normalize lengths to [0, 1] for color scaling
        if len(lengths) > 1:
            # Use log scale for better visualization of length distribution
            log_lengths = np.log(lengths + 1)  # +1 to avoid log(0)
            norm_lengths = (log_lengths - log_lengths.min()) / (log_lengths.max() - log_lengths.min())
        else:
            norm_lengths = np.array([0.5])  # Default for single factor

        # Additionally consider position to vary color intensity
        if len(positions) > 1:
            norm_positions = (positions - positions.min()) / (positions.max() - positions.min())
        else:
            norm_positions = np.array([0.5])  # Default for single factor

        # Plot each factor as a line
        for i, (pos, length, ref, is_rc, norm_len, norm_pos) in enumerate(zip(positions, lengths, refs, is_rcs, norm_lengths, norm_positions)):
            if is_rc:
                # Reverse complement: red line
                # x_init = pos, x_final = pos + length
                # y_init = ref + length, y_final = ref
                x_coords = [pos, pos + length]
                y_coords = [ref + length, ref]

                # Red color with intensity based on length (more opaque for longer factors)
                alpha = min((norm_len * 0.5 + norm_pos * 0.5), 1.0)  # Combine length and position for alpha
                color = (1.0, 0.0, 0.0, alpha)  # Red with variable alpha
            else:
                # Forward: blue line
                # x_init = pos, x_final = pos + length
                # y_init = ref, y_final = ref + length
                x_coords = [pos, pos + length]
                y_coords = [ref, ref + length]

                # Blue color with intensity based on length (more opaque for longer factors)
                alpha = min((norm_len * 0.5 + norm_pos * 0.5), 1.0)  # Combine length and position for alpha
                color = (0.0, 0.0, 1.0, alpha)  # Blue with variable alpha

            # Plot the line
            ax.plot(x_coords, y_coords, color=color, linewidth=1.5, alpha=alpha)

        # Add diagonal line y=x for reference
        max_val = max(positions.max() + lengths.max(), refs.max() + lengths.max())
        ax.plot([0, max_val], [0, max_val], color='gray', linestyle='--', linewidth=1, alpha=0.5)

        # Set axis labels and title
        ax.set_xlabel(f'Position in concatenated sequence ({name})')
        ax.set_ylabel(f'Reference position ({name})')
        ax.set_title(f'Weizmann Factor Plot - {name}')

        # Make axes equal for better visualization
        ax.set_aspect('equal', adjustable='box')

        # Add grid
        ax.grid(True, alpha=0.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='blue', alpha=0.8, label='Forward factors'),
            Patch(facecolor='red', alpha=0.8, label='Reverse complement factors')
        ]
        ax.legend(handles=legend_elements, loc='upper right')

        # Adjust layout
        plt.tight_layout()

        # Save plot if requested
        if save_path:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")

        # Show plot
        if show_plot:
            plt.show()
        else:
            plt.close(fig)

    except Exception as e:
        raise PlotError(f"Failed to create Weizmann factor plot: {e}")


