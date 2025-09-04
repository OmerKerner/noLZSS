"""
FASTA file parsing and compression utilities.

This module provides functions for reading, parsing, and compressing FASTA files
with proper handling of biological sequences and edge cases.
"""

from typing import Union, List, Tuple, Optional, Dict, Any
import re
from pathlib import Path

from ..utils import validate_input, NoLZSSError
from ..core import factorize
from .sequences import is_dna_sequence, is_protein_sequence, detect_sequence_type


class FASTAError(NoLZSSError):
    """Raised when FASTA file parsing or validation fails."""
    pass


def _parse_fasta_content(content: str) -> Dict[str, str]:
    """
    Parse FASTA format content into a dictionary of sequence IDs to sequences.
    
    Args:
        content: Raw FASTA file content as string
        
    Returns:
        Dictionary mapping sequence IDs to their sequences
        
    Raises:
        FASTAError: If FASTA format is invalid
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    for line_num, line in enumerate(content.splitlines(), 1):
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
            
            # Start new sequence
            header = line[1:].strip()
            if not header:
                raise FASTAError(f"Empty sequence header at line {line_num}")
            current_id = header.split()[0]  # Use first word as ID
            current_seq = []
        else:
            # Sequence line
            if current_id is None:
                raise FASTAError(f"Sequence data before header at line {line_num}")
            # Remove whitespace only (preserve all sequence characters for validation)
            clean_line = re.sub(r'\s', '', line.upper())
            current_seq.append(clean_line)
    
    # Save last sequence
    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)
    
    if not sequences:
        raise FASTAError("No valid sequences found in FASTA file")
    
    return sequences


def read_nucleotide_fasta(filepath: Union[str, Path]) -> List[Tuple[str, List[Tuple[int, int, int]]]]:
    """
    Read and factorize nucleotide sequences from a FASTA file.
    
    Only accepts sequences containing A, C, T, G (case insensitive).
    Sequences are converted to uppercase and factorized.
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        List of (sequence_id, factors) tuples where factors is the LZSS factorization
        
    Raises:
        FASTAError: If file format is invalid or contains invalid nucleotides
        FileNotFoundError: If file doesn't exist
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError as e:
        raise FASTAError(f"File encoding error: {e}")
    
    # Parse FASTA content
    sequences = _parse_fasta_content(content)
    
    results = []
    for seq_id, sequence in sequences.items():
        # Validate as DNA sequence (only A, C, T, G)
        if not re.match(r'^[ACGT]+$', sequence.upper()):
            invalid_chars = set(sequence.upper()) - set('ACGT')
            raise FASTAError(f"Sequence '{seq_id}' contains invalid nucleotides: {invalid_chars}")
        
        # Convert to uppercase
        sequence = sequence.upper()
        
        # Factorize the sequence
        try:
            factors = factorize(sequence.encode('ascii'))
            results.append((seq_id, factors))
        except Exception as e:
            raise FASTAError(f"Failed to factorize sequence '{seq_id}': {e}")
    
    return results


def read_protein_fasta(filepath: Union[str, Path]) -> List[Tuple[str, str]]:
    """
    Read amino acid sequences from a FASTA file.
    
    Only accepts sequences containing canonical amino acids.
    Sequences are converted to uppercase.
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        List of (sequence_id, sequence) tuples
        
    Raises:
        FASTAError: If file format is invalid or contains invalid amino acids
        FileNotFoundError: If file doesn't exist
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError as e:
        raise FASTAError(f"File encoding error: {e}")
    
    # Parse FASTA content
    sequences = _parse_fasta_content(content)
    
    results = []
    # Canonical amino acids (20 standard)
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    
    for seq_id, sequence in sequences.items():
        # Convert to uppercase
        sequence = sequence.upper()
        
        # Validate amino acids
        if not all(c in valid_aa for c in sequence):
            invalid_chars = set(sequence) - valid_aa
            raise FASTAError(f"Sequence '{seq_id}' contains invalid amino acids: {invalid_chars}")
        
        results.append((seq_id, sequence))
    
    return results


def read_fasta_auto(filepath: Union[str, Path]) -> Union[
    List[Tuple[str, List[Tuple[int, int, int]]]],  # For nucleotide
    List[Tuple[str, str]]  # For protein
]:
    """
    Read a FASTA file and automatically detect whether it contains nucleotide or amino acid sequences.
    
    For nucleotide sequences: validates A,C,T,G only and returns factorized results
    For amino acid sequences: validates canonical amino acids and returns sequences
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        For nucleotide FASTA: List of (sequence_id, factors) tuples
        For amino acid FASTA: List of (sequence_id, sequence) tuples
        
    Raises:
        FASTAError: If file format is invalid or sequence type cannot be determined
        FileNotFoundError: If file doesn't exist
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError as e:
        raise FASTAError(f"File encoding error: {e}")
    
    # Parse FASTA content
    sequences = _parse_fasta_content(content)
    
    if not sequences:
        raise FASTAError("No sequences found in FASTA file")
    
    # Sample sequences to detect type
    sample_seq = next(iter(sequences.values()))
    
    # Detect sequence type
    seq_type = detect_sequence_type(sample_seq)
    
    if seq_type == 'dna':
        # Process as nucleotide
        return read_nucleotide_fasta(filepath)
    elif seq_type == 'protein':
        # Process as protein
        return read_protein_fasta(filepath)
    else:
        raise FASTAError(f"Cannot determine sequence type. Detected: {seq_type}. "
                        f"Expected DNA (A,C,T,G) or protein (amino acids) sequences.")


def process_fasta_with_plots(
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
        FASTAError: If FASTA processing fails
        FileNotFoundError: If input file doesn't exist
    """
    from ..core import write_factors_binary_file
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend for batch processing
    
    fasta_filepath = Path(fasta_filepath)
    output_dir = Path(output_dir)
    
    if not fasta_filepath.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_filepath}")
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read FASTA file
    sequences = _parse_fasta_content(fasta_filepath.read_text())
    
    if not sequences:
        raise FASTAError("No sequences found in FASTA file")
    
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
            print(f"  Detected nucleotide sequence")
            
        elif seq_type == 'protein':
            # Validate as amino acid
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
            if not all(c in valid_aa for c in sequence.upper()):
                invalid_chars = set(sequence.upper()) - valid_aa
                print(f"  Warning: Skipping {seq_id} - contains invalid amino acids: {invalid_chars}")
                continue
            sequence = sequence.upper()
            print(f"  Detected amino acid sequence")
            
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

