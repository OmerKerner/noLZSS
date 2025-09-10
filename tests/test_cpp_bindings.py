import noLZSS
import tempfile
import os
import struct

def check_invariants(text: bytes):
    factors = noLZSS.factorize(text)
    n = len(text)
    covered = 0
    prev_end = 0
    for start, length, ref in factors:
        assert 0 <= start < n
        assert length > 0
        assert start >= prev_end
        prev_end = start + length
        covered += length
    assert covered == n
    return factors

def test_repeated():
    check_invariants(b"aaaaa")

def test_mixed():
    check_invariants(b"abracadabra")

def test_short():
    check_invariants(b"a")

def test_version():
    assert hasattr(noLZSS, "factorize")
    assert hasattr(noLZSS, "__version__")
    assert isinstance(noLZSS.__version__, str)

def test_count_factors():
    """Test count_factors function"""
    text = b"aaaaa"
    factors = noLZSS.factorize(text)
    count = noLZSS.count_factors(text)
    assert count == len(factors)
    
    text2 = b"abracadabra"
    factors2 = noLZSS.factorize(text2)
    count2 = noLZSS.count_factors(text2)
    assert count2 == len(factors2)
    
    # Test with different text
    text3 = b"hello"
    count3 = noLZSS.count_factors(text3)
    factors3 = noLZSS.factorize(text3)
    assert count3 == len(factors3)

def test_count_factors_file():
    """Test count_factors_file function"""
    # Create temporary file
    with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
        f.write(b"aaaaa")
        temp_path = f.name
    
    try:
        count = noLZSS.count_factors_file(temp_path)
        factors = noLZSS.factorize_file(temp_path, 0)
        assert count == len(factors)
        
        # Test with different content
        with open(temp_path, 'wb') as f:
            f.write(b"abracadabra")
        
        count2 = noLZSS.count_factors_file(temp_path)
        factors2 = noLZSS.factorize_file(temp_path, 0)
        assert count2 == len(factors2)
        
    finally:
        os.unlink(temp_path)

def test_factorize_file():
    """Test factorize_file function"""
    # Create temporary file
    with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
        f.write(b"aaaaa")
        temp_path = f.name
    
    try:
        # Test factorize_file
        factors_file = noLZSS.factorize_file(temp_path, 0)
        factors_memory = noLZSS.factorize(b"aaaaa")
        assert factors_file == factors_memory
        
        # Test with reserve_hint
        factors_file_reserved = noLZSS.factorize_file(temp_path, 10)
        assert factors_file_reserved == factors_file
        
        # Test with different content
        with open(temp_path, 'wb') as f:
            f.write(b"abracadabra")
        
        factors_file2 = noLZSS.factorize_file(temp_path, 0)
        factors_memory2 = noLZSS.factorize(b"abracadabra")
        assert factors_file2 == factors_memory2
        
    finally:
        os.unlink(temp_path)

def test_write_factors_binary_file():
    """Test write_factors_binary_file function"""
    # Create input file
    with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
        f.write(b"aaaaa")
        in_path = f.name
    
    # Create output file
    with tempfile.NamedTemporaryFile(delete=False) as f:
        out_path = f.name
    
    try:
        # Write binary factors
        noLZSS.write_factors_binary_file(in_path, out_path)
        
        # Get count from memory factorization for comparison
        factors_memory = noLZSS.factorize(b"aaaaa")
        expected_count = len(factors_memory)
        
        # Read binary file and verify
        with open(out_path, 'rb') as f:
            binary_data = f.read()
        
        # Each factor is 24 bytes (3 uint64_t)
        actual_count = len(binary_data) // 24
        assert actual_count == expected_count
        
        # Parse binary data
        binary_factors = []
        for i in range(actual_count):
            start, length, ref = struct.unpack('<QQQ', binary_data[i*24:(i+1)*24])
            binary_factors.append((start, length, ref))
        
        assert binary_factors == factors_memory
        
        # Test a second call to ensure consistent behavior
        noLZSS.write_factors_binary_file(in_path, out_path)
        with open(out_path, 'rb') as f:
            binary_data2 = f.read()
        actual_count2 = len(binary_data2) // 24
        assert actual_count2 == expected_count
        
    finally:
        os.unlink(in_path)
        os.unlink(out_path)

def test_consistency_across_functions():
    """Test that all functions give consistent results"""
    test_texts = [
        b"a",
        b"aa",
        b"aaa",
        b"aaaa",
        b"aaaaa",
        b"abracadabra",
        b"mississippi"
    ]
    
    for text in test_texts:
        # Get factors from memory
        factors = noLZSS.factorize(text)
        count = noLZSS.count_factors(text)
        
        # Create temporary file
        with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
            f.write(text)
            temp_path = f.name
        
        try:
            # Get factors from file
            factors_file = noLZSS.factorize_file(temp_path, 0)
            count_file = noLZSS.count_factors_file(temp_path)
            
            # All should be consistent
            assert factors == factors_file
            assert count == len(factors)
            assert count_file == len(factors_file)
            assert count == count_file
            
            # Test binary output
            with tempfile.NamedTemporaryFile(delete=False) as f:
                bin_path = f.name
            
            try:
                noLZSS.write_factors_binary_file(temp_path, bin_path)
                # Verify binary file was created and has correct size
                with open(bin_path, 'rb') as f:
                    binary_data = f.read()
                expected_size = len(factors) * 24  # 24 bytes per factor
                assert len(binary_data) == expected_size
            finally:
                os.unlink(bin_path)
                
        finally:
            os.unlink(temp_path)

def test_prepare_multiple_dna_sequences():
    """Test prepare_multiple_dna_sequences function"""
    import noLZSS._noLZSS as cpp
    
    # Test basic functionality with valid DNA sequences
    sequences = ["ATCG", "GCTA"]
    result = cpp.prepare_multiple_dna_sequences(sequences)
    
    assert isinstance(result, tuple)
    assert len(result) == 2
    concatenated, original_length = result
    assert isinstance(concatenated, str)
    assert isinstance(original_length, int)
    
    # Should contain original sequences + sentinels + reverse complements + sentinels
    assert len(concatenated) > len("ATCGGCTA")  # At least original sequences
    assert original_length < len(concatenated)  # Original length should be less than total
    
    # Test with single sequence
    single_result = cpp.prepare_multiple_dna_sequences(["ATCG"])
    assert isinstance(single_result, tuple)
    assert len(single_result) == 2
    
    # Test with mixed case (should be converted to uppercase)
    mixed_result = cpp.prepare_multiple_dna_sequences(["atcg", "GCTA"])
    mixed_concatenated, _ = mixed_result
    # Should not contain lowercase letters
    assert 'a' not in mixed_concatenated
    assert 'c' not in mixed_concatenated
    assert 'g' not in mixed_concatenated
    assert 't' not in mixed_concatenated

def test_prepare_multiple_dna_sequences_validation():
    """Test prepare_multiple_dna_sequences input validation"""
    import noLZSS._noLZSS as cpp
    
    # Test empty sequences list
    result = cpp.prepare_multiple_dna_sequences([])
    assert result == ("", 0)
    
    # Test invalid nucleotides
    try:
        cpp.prepare_multiple_dna_sequences(["ATCGN"])  # N is not valid
        assert False, "Should have raised an exception for invalid nucleotide"
    except RuntimeError as e:
        assert "Invalid nucleotide" in str(e)
    
    # Test sequences with numbers
    try:
        cpp.prepare_multiple_dna_sequences(["ATCG1"])
        assert False, "Should have raised an exception for invalid nucleotide"
    except RuntimeError as e:
        assert "Invalid nucleotide" in str(e)
    
    # Test too many sequences (more than 125)
    try:
        many_sequences = ["ATCG"] * 126
        cpp.prepare_multiple_dna_sequences(many_sequences)
        assert False, "Should have raised an exception for too many sequences"
    except Exception as e:
        assert "Too many sequences" in str(e)

def test_factorize_multiple_dna_w_rc():
    """Test factorize_multiple_dna_w_rc function"""
    import noLZSS._noLZSS as cpp
    
    # Test basic functionality
    sequences = ["ATCG", "CGTA"]
    prepared, _ = cpp.prepare_multiple_dna_sequences(sequences)
    # Convert string to bytes for the factorization function
    prepared_bytes = prepared.encode('latin-1')
    factors = cpp.factorize_multiple_dna_w_rc(prepared_bytes)
    
    assert isinstance(factors, list)
    assert len(factors) > 0
    
    # Each factor should be a tuple with 4 elements (start, length, ref, is_rc)
    for factor in factors:
        assert isinstance(factor, tuple)
        assert len(factor) == 4
        start, length, ref, is_rc = factor
        assert isinstance(start, int)
        assert isinstance(length, int)
        assert isinstance(ref, int)
        assert isinstance(is_rc, bool)
        assert length > 0  # Length should be positive

def test_count_factors_multiple_dna_w_rc():
    """Test count_factors_multiple_dna_w_rc function"""
    import noLZSS._noLZSS as cpp
    
    # Test basic functionality
    sequences = ["ATCG", "CGTA"]
    prepared, _ = cpp.prepare_multiple_dna_sequences(sequences)
    prepared_bytes = prepared.encode('latin-1')
    
    # Get factors and count
    factors = cpp.factorize_multiple_dna_w_rc(prepared_bytes)
    count = cpp.count_factors_multiple_dna_w_rc(prepared_bytes)
    
    # Count should match number of factors
    assert count == len(factors)
    
    # Test with single sequence
    single_seq = ["ATCG"]
    single_prepared, _ = cpp.prepare_multiple_dna_sequences(single_seq)
    single_prepared_bytes = single_prepared.encode('latin-1')
    single_count = cpp.count_factors_multiple_dna_w_rc(single_prepared_bytes)
    assert single_count > 0

def test_factorize_fasta_multiple_dna_w_rc():
    """Test factorize_fasta_multiple_dna_w_rc function"""
    import noLZSS._noLZSS as cpp
    
    # Create a test FASTA file
    fasta_content = """>seq1
ATCGATCG
>seq2
GCTAGCTA
>seq3
TTTTAAAA
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        temp_path = f.name
    
    try:
        # Test factorization
        factors = cpp.factorize_fasta_multiple_dna_w_rc(temp_path)
        
        assert isinstance(factors, list)
        assert len(factors) > 0
        
        # Each factor should be valid
        for factor in factors:
            assert isinstance(factor, tuple)
            assert len(factor) == 4  # (start, length, ref, is_rc)
            start, length, ref, is_rc = factor
            assert isinstance(start, int)
            assert isinstance(length, int)
            assert isinstance(ref, int)
            assert isinstance(is_rc, bool)
            assert length > 0
            
    finally:
        os.unlink(temp_path)

def test_factorize_fasta_multiple_dna_w_rc_validation():
    """Test factorize_fasta_multiple_dna_w_rc input validation"""
    import noLZSS._noLZSS as cpp
    
    # Test with invalid nucleotides in FASTA
    invalid_fasta = """>seq1
ATCGNTCG
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(invalid_fasta)
        temp_path = f.name
    
    try:
        try:
            cpp.factorize_fasta_multiple_dna_w_rc(temp_path)
            assert False, "Should have raised an exception for invalid nucleotide"
        except RuntimeError as e:
            assert "Invalid nucleotide" in str(e)
    finally:
        os.unlink(temp_path)
    
    # Test with non-existent file
    try:
        cpp.factorize_fasta_multiple_dna_w_rc("/non/existent/file.fasta")
        assert False, "Should have raised an exception for non-existent file"
    except RuntimeError as e:
        assert "Cannot open FASTA file" in str(e)
    
    # Test with empty FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write("")
        empty_path = f.name
    
    try:
        try:
            cpp.factorize_fasta_multiple_dna_w_rc(empty_path)
            assert False, "Should have raised an exception for empty FASTA"
        except RuntimeError as e:
            assert "No valid sequences found" in str(e)
    finally:
        os.unlink(empty_path)

def test_dna_w_rc_functions_integration():
    """Test integration between DNA reverse complement functions"""
    import noLZSS._noLZSS as cpp
    
    # Test that single DNA factorization works
    dna_text = "ATCGATCG"
    dna_bytes = dna_text.encode('latin-1')
    single_factors = cpp.factorize_dna_w_rc(dna_bytes)
    single_count = cpp.count_factors_dna_w_rc(dna_bytes)
    
    assert len(single_factors) == single_count
    
    # Test multiple DNA functions with single sequence
    sequences = [dna_text]
    prepared, _ = cpp.prepare_multiple_dna_sequences(sequences)
    prepared_bytes = prepared.encode('latin-1')
    multi_factors = cpp.factorize_multiple_dna_w_rc(prepared_bytes)
    multi_count = cpp.count_factors_multiple_dna_w_rc(prepared_bytes)
    
    assert len(multi_factors) == multi_count
    
    # The multiple sequence version should handle the single sequence case
    assert len(multi_factors) > 0
    assert len(single_factors) > 0

def test_reverse_complement_awareness():
    """Test that reverse complement matching works correctly"""
    import noLZSS._noLZSS as cpp
    
    # Create sequences where reverse complement matches should occur
    # ATCG reverse complement is CGAT
    sequences = ["ATCG", "CGAT", "ATCG"]  # Third sequence should match first
    prepared, _ = cpp.prepare_multiple_dna_sequences(sequences)
    prepared_bytes = prepared.encode('latin-1')
    factors = cpp.factorize_multiple_dna_w_rc(prepared_bytes)
    
    assert len(factors) > 0
    
    # Check that factors have valid references and RC flags
    for factor in factors:
        start, length, ref, is_rc = factor
        # Reference should be valid
        assert ref >= 0
        # RC flag should be boolean
        assert isinstance(is_rc, bool)
        # Length should be positive
        assert length > 0

def test_multiple_dna_binary_output():
    """Test binary output functions for multiple DNA sequences"""
    import noLZSS._noLZSS as cpp
    
    # Create test FASTA file
    fasta_content = """>seq1
ATCGATCG
>seq2
GCTAGCTA
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_path = f.name
    
    with tempfile.NamedTemporaryFile(delete=False) as f:
        binary_path = f.name
    
    try:
        # Test binary output
        num_factors = cpp.write_factors_binary_file_multiple_dna_w_rc(fasta_path, binary_path)
        
        assert num_factors > 0
        
        # Verify binary file was created and has correct size
        with open(binary_path, 'rb') as f:
            binary_data = f.read()
        
        expected_size = num_factors * 24  # 24 bytes per factor (3 * 8 bytes)
        assert len(binary_data) == expected_size
        
        # Compare with regular factorization
        factors = cpp.factorize_fasta_multiple_dna_w_rc(fasta_path)
        assert len(factors) == num_factors
        
    finally:
        os.unlink(fasta_path)
        os.unlink(binary_path)
