"""
Test factorization correctness by validating that factors correspond to actual substring matches.

This test uses the viral genome FASTA file to test factorize_fasta_multiple_dna_w_rc
and validates that each factor correctly represents a substring match or reverse complement match.
"""

import pytest
import os
from pathlib import Path
from typing import List

# Import test helpers
try:
    import noLZSS
    from noLZSS.genomics import factorize_fasta_multiple_dna_w_rc, prepare_multiple_dna_sequences_w_rc
    HAS_CPP_EXTENSION = True
except ImportError:
    HAS_CPP_EXTENSION = False


def reverse_complement(dna_sequence: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.
    
    Args:
        dna_sequence: DNA string containing only A, C, T, G
        
    Returns:
        Reverse complement of the input sequence
        
    Raises:
        ValueError: If sequence contains invalid nucleotides
    """
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    try:
        # Build reverse complement by going backwards through sequence
        rc = ''.join(complement_map[nucleotide] for nucleotide in reversed(dna_sequence))
        return rc
    except KeyError as e:
        raise ValueError(f"Invalid nucleotide found: {e.args[0]}")


def parse_fasta_sequences(fasta_path: str) -> List[str]:
    """
    Parse FASTA file and return list of DNA sequences (uppercase, validated).
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        List of DNA sequences from the FASTA file
    """
    sequences = []
    current_sequence = ""
    
    with open(fasta_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence if any
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
            else:
                # Accumulate sequence data, convert to uppercase and validate
                for char in line:
                    if char in 'ACTGactg':
                        current_sequence += char.upper()
                    elif char not in ' \t\n\r':  # Ignore whitespace
                        raise ValueError(f"Invalid nucleotide found: {char}")
        
        # Don't forget the last sequence
        if current_sequence:
            sequences.append(current_sequence)
    
    if not sequences:
        raise ValueError("No sequences found in FASTA file")
    
    return sequences


@pytest.mark.skipif(not HAS_CPP_EXTENSION, reason="C++ extension not available")
@pytest.mark.parametrize("fasta_path", ["tests/resources/test_viral_dna.fna", "tests/resources/test_bacterial_dna.fna"])
class TestFactorizationCorrectness:
    """Test factorization correctness using real genomic data."""
    
    def test_factorization_string_validation(self, fasta_path):
        """
        Test that each factor in the factorization correctly represents a substring match.
        
        This test validates that:
        1. Normal factors match substrings exactly
        2. Reverse complement factors match the reverse complement of reference substrings
        3. No factors are "bad" (i.e., all represent valid matches)
        """
        self.test_fasta_path = fasta_path
        
        # Check if test file exists
        if not os.path.exists(self.test_fasta_path):
            pytest.skip(f"Test FASTA file not found: {self.test_fasta_path}")
        # Parse FASTA file to get individual sequences
        sequences = parse_fasta_sequences(self.test_fasta_path)
        print(f"Parsed {len(sequences)} sequences from FASTA file")
        
        # Prepare sequences using the same method as the C++ function
        concatenated_string, original_length, sentinel_positions = prepare_multiple_dna_sequences_w_rc(sequences)
        print(f"Prepared string: original_length={original_length}, total_length={len(concatenated_string)}")
        
        # Perform factorization
        factors_list, sentinel_indices, sequence_ids = factorize_fasta_multiple_dna_w_rc(self.test_fasta_path)
        print(f"Factorization produced {len(factors_list)} factors with {len(sentinel_indices)} sentinels")
        
        # Validate each factor
        bad_factors = []
        total_factors = len(factors_list)
        
        for i, factor in enumerate(factors_list):
            start, length, ref, is_rc = factor
            
            # Skip validation for factors that might be sentinels or extend beyond original sequence
            if start >= original_length:
                continue
                
            # Ensure factor bounds are valid
            if start + length > original_length:
                print(f"WARNING: Factor {i} extends beyond original sequence: {factor}")
                continue
                
            if ref + length > len(concatenated_string):
                print(f"WARNING: Factor {i} reference extends beyond total sequence: {factor}")
                continue
            
            # Extract the substring at the factor position
            factor_substring = concatenated_string[start:start + length]
            
            # Extract the reference substring
            reference_substring = concatenated_string[ref:ref + length]
            
            # Compute expected substring (reverse complement if needed)
            if is_rc:
                try:
                    expected_substring = reverse_complement(reference_substring)
                except ValueError as e:
                    print(f"ERROR: Factor {i} reverse complement computation failed: {factor}")
                    print(f"  Reference substring: '{reference_substring}'")
                    print(f"  Error: {e}")
                    bad_factors.append((i, factor, f"RC computation failed: {e}"))
                    continue
            else:
                expected_substring = reference_substring
            
            # Compare the substrings
            if factor_substring != expected_substring:
                bad_factors.append((i, factor, f"Mismatch: got '{factor_substring}', expected '{expected_substring}'"))
                print(f"BAD FACTOR {i}: {factor}")
                print(f"  Factor substring:    '{factor_substring}'")
                print(f"  Expected substring:  '{expected_substring}'")
                print(f"  Reference substring: '{reference_substring}'")
                print(f"  Is reverse complement: {is_rc}")
        
        # Report results
        print(f"Validation complete. Bad factors: {len(bad_factors)} out of {total_factors}")
        
        if bad_factors:
            print("\nAll bad factors:")
            for factor_index, factor, reason in bad_factors:
                print(f"  Factor {factor_index}: {factor} - {reason}")
        
        # Assert that no factors are bad
        assert len(bad_factors) == 0, f"Found {len(bad_factors)} bad factors out of {total_factors}"
    
    def test_factorization_coverage(self, fasta_path):
        """
        Test that factors provide complete coverage of the original sequence.
        """
        self.test_fasta_path = fasta_path
        
        # Check if test file exists
        if not os.path.exists(self.test_fasta_path):
            pytest.skip(f"Test FASTA file not found: {self.test_fasta_path}")
        # Parse sequences and prepare
        sequences = parse_fasta_sequences(self.test_fasta_path)
        concatenated_string, original_length, sentinel_positions = prepare_multiple_dna_sequences_w_rc(sequences)
        
        # Perform factorization
        factors_list, sentinel_indices, sequence_ids = factorize_fasta_multiple_dna_w_rc(self.test_fasta_path)
        
        # Filter factors within original sequence and sort by start position
        valid_factors = sorted([f for f in factors_list if f[0] < original_length], key=lambda x: x[0])
        
        if not valid_factors:
            assert False, "No valid factors found within original sequence"
        
        # Check for complete coverage by ensuring no gaps between sorted factors
        current_end = 0
        for start, length, ref, is_rc in valid_factors:
            if start > current_end:
                # Gap found
                uncovered = start - current_end
                assert False, f"Gap of {uncovered} uncovered positions starting at {current_end}"
            current_end = max(current_end, start + length)
        
        # Ensure we reach the end of the original sequence
        assert current_end >= (original_length - 1), f"Sequence not fully covered: reached {current_end} of {original_length}"
        
        print(f"Coverage analysis: Complete coverage verified for {original_length} positions")
    
    def test_factorization_basic_properties(self, fasta_path):
        """
        Test basic properties of the factorization result.
        """
        self.test_fasta_path = fasta_path
        
        # Check if test file exists
        if not os.path.exists(self.test_fasta_path):
            pytest.skip(f"Test FASTA file not found: {self.test_fasta_path}")
        factors_list, sentinel_indices, sequence_ids = factorize_fasta_multiple_dna_w_rc(self.test_fasta_path)
        
        # Basic sanity checks
        assert len(factors_list) > 0, "Should produce at least one factor"
        assert len(sequence_ids) > 0, "Should have at least one sequence ID"
        
        # Check factor properties
        for i, factor in enumerate(factors_list):
            start, length, ref, is_rc = factor
            
            assert length > 0, f"Factor {i} should have positive length: {factor}"
            assert start >= 0, f"Factor {i} should have non-negative start: {factor}"
            assert ref >= 0, f"Factor {i} should have non-negative reference: {factor}"
            assert isinstance(is_rc, bool), f"Factor {i} is_rc should be boolean: {factor}"
        
        # Check sentinel indices are valid
        for sentinel_idx in sentinel_indices:
            assert 0 <= sentinel_idx < len(factors_list), f"Sentinel index {sentinel_idx} out of bounds"
        
        print(f"Basic properties validated for {len(factors_list)} factors")
    
    def test_reverse_complement_functionality(self, fasta_path):
        """
        Test that reverse complement functionality is working as expected.
        """
        self.test_fasta_path = fasta_path
        
        # Check if test file exists
        if not os.path.exists(self.test_fasta_path):
            pytest.skip(f"Test FASTA file not found: {self.test_fasta_path}")
        factors_list, sentinel_indices, sequence_ids = factorize_fasta_multiple_dna_w_rc(self.test_fasta_path)
        
        # Count reverse complement vs normal factors
        rc_factors = sum(1 for factor in factors_list if factor[3])  # is_rc is 4th element
        normal_factors = len(factors_list) - rc_factors
        
        print(f"Factor analysis: {normal_factors} normal factors, {rc_factors} reverse complement factors")
        
        # We expect at least some normal factors
        assert normal_factors > 0, "Should have at least some normal factors"
        
        # Total should match
        assert rc_factors + normal_factors == len(factors_list), "RC + normal should equal total factors"


class TestReverseComplementHelper:
    """Test the reverse complement helper function."""
    
    def test_reverse_complement_basic(self):
        """Test basic reverse complement functionality."""
        assert reverse_complement("A") == "T"
        assert reverse_complement("T") == "A"
        assert reverse_complement("C") == "G"
        assert reverse_complement("G") == "C"
    
    def test_reverse_complement_sequences(self):
        """Test reverse complement on longer sequences."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("ATGC") == "GCAT"
        assert reverse_complement("GCTA") == "TAGC"
    
    def test_reverse_complement_invalid(self):
        """Test reverse complement with invalid nucleotides."""
        with pytest.raises(ValueError, match="Invalid nucleotide"):
            reverse_complement("ATCGN")
        
        with pytest.raises(ValueError, match="Invalid nucleotide"):
            reverse_complement("ATCG123")


class TestFastaParsingHelper:
    """Test the FASTA parsing helper function."""
    
    @pytest.mark.parametrize("fasta_file", ["tests/resources/test_viral_dna.fna", "tests/resources/test_bacterial_dna.fna"])
    def test_parse_fasta_basic(self, fasta_file):
        """Test basic FASTA parsing functionality."""
        # This test only runs if the test file exists
        if not os.path.exists(fasta_file):
            pytest.skip(f"Test FASTA file not found: {fasta_file}")

        sequences = parse_fasta_sequences(fasta_file)

        assert len(sequences) > 0, "Should parse at least one sequence"
        
        for seq in sequences:
            assert len(seq) > 0, "Each sequence should be non-empty"
            assert all(nucleotide in "ACTG" for nucleotide in seq), "All nucleotides should be valid"