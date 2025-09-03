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
