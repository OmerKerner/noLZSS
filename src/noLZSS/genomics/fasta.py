"""
FASTA file parsing and compression utilities.

This module provides functions for reading, parsing, and compressing FASTA files
with proper handling of biological sequences and edge cases.
"""

from typing import Union, List, Tuple, Dict
import re
from pathlib import Path

from ..utils import NoLZSSError
from ..core import factorize

try:
    from .._noLZSS import write_factors_dna_w_reference_fasta_files_to_binary as _write_factors_dna_w_reference_fasta_files_to_binary
except ImportError:
    _write_factors_dna_w_reference_fasta_files_to_binary = None

from .sequences import detect_sequence_type


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


def write_factors_dna_w_reference_fasta_files_to_binary(reference_fasta_path: Union[str, Path], 
                                                        target_fasta_path: Union[str, Path],
                                                        output_path: Union[str, Path]) -> int:
    """
    Factorize DNA sequences from FASTA files with reference and write factors to binary file.
    
    Reads DNA sequences from reference and target FASTA files, performs noLZSS factorization
    of the target using the reference, and writes the resulting factors to a binary file.
    Specialized for nucleotide sequences (A, C, T, G) with reverse complement matching capability.
    
    Args:
        reference_fasta_path: Path to reference FASTA file containing DNA sequences
        target_fasta_path: Path to target FASTA file containing DNA sequences to factorize
        output_path: Path to output file where binary factors will be written
        
    Returns:
        Number of factors written to the output file
        
    Raises:
        ValueError: If files contain empty sequences or invalid nucleotides
        FileNotFoundError: If FASTA files do not exist
        RuntimeError: If unable to read FASTA files, create output file, or processing errors occur
        FASTAError: If C++ extension is not available
        
    Note:
        - Factor start positions are absolute positions in the combined reference+target string
        - Supports reverse complement matching for DNA sequences (indicated by MSB in reference field)
        - All sequences from both FASTA files are concatenated with sentinel separators
        - Only nucleotides A, C, T, G are allowed (case-insensitive)
        - This function overwrites the output file if it exists
        
    Warning:
        Characters 1-251 are used as sentinel separators and must not appear in sequences.
    """
    # Check if C++ extension is available
    if _write_factors_dna_w_reference_fasta_files_to_binary is None:
        raise FASTAError("C++ extension not available. Cannot process FASTA files with reference sequences.")
    
    # Validate input file paths
    reference_path = Path(reference_fasta_path)
    target_path = Path(target_fasta_path)
    
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference FASTA file not found: {reference_path}")
    if not target_path.exists():
        raise FileNotFoundError(f"Target FASTA file not found: {target_path}")
    
    # Ensure output directory exists
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    return _write_factors_dna_w_reference_fasta_files_to_binary(str(reference_path), str(target_path), str(output_path))

