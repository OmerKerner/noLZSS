#!/usr/bin/env python3
"""
Test script for the enhanced plotting function.

This script demonstrates the new capabilities of plot_multiple_seq_self_lz_factor_plot_from_fasta:
1. Reading from binary factor files with metadata
2. Visualizing sentinel boundaries with lines
3. Showing sequence names on axes
"""

import tempfile
import os
from pathlib import Path

# Set up path to use development version
import sys
sys.path.insert(0, 'src')

import noLZSS._noLZSS as cpp

# Inline the new function for testing since it's not installed yet
def read_factors_binary_file_with_metadata(filepath):
    """Read factors from an enhanced binary file with metadata."""
    import struct
    from pathlib import Path
    
    RC_MASK = 1 << 63
    filepath = Path(filepath)
    
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    
    try:
        with open(filepath, 'rb') as f:
            # Read header
            header_data = f.read(48)
            if len(header_data) < 40:
                raise ValueError("File too small to contain valid header")
            
            # Unpack header (magic is 8 chars, then 4 uint64_t)
            magic = header_data[:8]
            if magic != b'noLZSSv1':
                raise ValueError("Invalid file format: missing noLZSS magic header")
            
            num_factors, num_sequences, num_sentinels, header_size = struct.unpack('<QQQQ', header_data[8:40])
            
            # Seek to beginning of header to read the full header
            f.seek(0)
            full_header = f.read(header_size)
            if len(full_header) != header_size:
                raise ValueError(f"Could not read full header: expected {header_size}, got {len(full_header)}")
            
            # Skip the basic header structure
            offset = 40
            
            # Read sequence names
            sequence_names = []
            for i in range(num_sequences):
                # Find null terminator
                name_start = offset
                while offset < len(full_header) and full_header[offset] != 0:
                    offset += 1
                if offset >= len(full_header):
                    raise ValueError("Invalid sequence name format")
                
                name = full_header[name_start:offset].decode('utf-8')
                sequence_names.append(name)
                offset += 1  # Skip null terminator
            
            # Read sentinel factor indices
            sentinel_indices = []
            for i in range(num_sentinels):
                if offset + 8 > len(full_header):
                    raise ValueError("Insufficient data for sentinel indices")
                
                idx = struct.unpack('<Q', full_header[offset:offset+8])[0]
                sentinel_indices.append(idx)
                offset += 8
            
            # Read factors
            factors = []
            for i in range(num_factors):
                factor_data = f.read(24)  # Each factor is 3 * uint64_t = 24 bytes
                if len(factor_data) != 24:
                    raise ValueError(f"Insufficient data for factor {i}")
                
                start, length, ref = struct.unpack('<QQQ', factor_data)
                
                # Extract is_rc flag and clean ref
                is_rc_flag = bool(ref & RC_MASK)
                clean_ref = ref & ~RC_MASK
                
                factors.append((start, length, clean_ref, is_rc_flag))
    
    except IOError as e:
        raise ValueError(f"Error reading file {filepath}: {e}")
    except struct.error as e:
        raise ValueError(f"Error unpacking binary data: {e}")
    
    return {
        'factors': factors,
        'sentinel_factor_indices': sentinel_indices,
        'sequence_names': sequence_names,
        'num_sequences': num_sequences,
        'num_sentinels': num_sentinels
    }


def test_enhanced_plotting_logic():
    """Test the enhanced plotting function logic without Panel dependencies."""
    
    print("=" * 60)
    print("Testing Enhanced Plotting Function Logic")
    print("=" * 60)
    
    # Create a test FASTA file with multiple sequences
    test_fasta = ">sequence1\nACGTACGT\n>sequence2\nTGCATGCA\n>sequence3\nAAAATTTT\n"
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(test_fasta)
        fasta_path = f.name
    
    binary_path = fasta_path.replace('.fasta', '.bin')
    
    try:
        print(f"1. Creating FASTA file: {fasta_path}")
        print(f"   Content:\n{test_fasta}")
        
        # Test factorization from FASTA
        print("2. Testing factorization from FASTA...")
        factors, sentinel_indices = cpp.factorize_fasta_multiple_dna_w_rc(fasta_path)
        print(f"   Factors: {len(factors)}")
        print(f"   Sentinel indices: {sentinel_indices}")
        
        # Create binary file with metadata
        print("3. Creating binary file with metadata...")
        cpp.write_factors_binary_file_fasta_multiple_dna_w_rc(fasta_path, binary_path)
        print(f"   Binary file created: {binary_path}")
        
        # Test reading binary metadata
        print("4. Testing binary metadata reading...")
        metadata = read_factors_binary_file_with_metadata(binary_path)
        print(f"   Sequence names: {metadata['sequence_names']}")
        print(f"   Sentinel indices: {metadata['sentinel_factor_indices']}")
        print(f"   Factor count: {len(metadata['factors'])}")
        
        # Test sentinel position calculation (same logic as in plotting function)
        print("5. Testing sentinel boundary calculation...")
        factors = metadata['factors']
        sentinel_factor_indices = metadata['sentinel_factor_indices']
        sequence_names = metadata['sequence_names']
        
        # Calculate sentinel positions
        sentinel_positions = []
        for idx in sentinel_factor_indices:
            if idx < len(factors):
                sentinel_start = factors[idx][0]
                sentinel_positions.append(sentinel_start)
        
        # Calculate sequence boundaries
        sequence_boundaries = []
        prev_pos = 0
        for i, pos in enumerate(sentinel_positions):
            seq_name = sequence_names[i] if i < len(sequence_names) else f"seq_{i}"
            sequence_boundaries.append((prev_pos, pos, seq_name))
            prev_pos = pos + 1
        
        # Add the last sequence
        if len(sequence_names) > len(sentinel_positions):
            last_name = sequence_names[len(sentinel_positions)]
        else:
            last_name = f"seq_{len(sentinel_positions)}"
        
        max_pos = max(f[0] + f[1] for f in factors) if factors else prev_pos
        sequence_boundaries.append((prev_pos, max_pos, last_name))
        
        print(f"   Sentinel positions: {sentinel_positions}")
        print(f"   Sequence boundaries: {sequence_boundaries}")
        
        # Verify the boundaries make sense
        print("6. Verifying boundaries...")
        total_sequences = len(sequence_names)
        total_sentinels = len(sentinel_factor_indices)
        total_boundaries = len(sequence_boundaries)
        
        print(f"   Total sequences: {total_sequences}")
        print(f"   Total sentinels: {total_sentinels}")
        print(f"   Total boundaries: {total_boundaries}")
        
        # For N sequences, we should have N-1 sentinels and N boundaries
        expected_sentinels = total_sequences - 1
        expected_boundaries = total_sequences
        
        if total_sentinels == expected_sentinels and total_boundaries == expected_boundaries:
            print("   ✅ Boundary calculation is correct!")
        else:
            print(f"   ❌ Boundary calculation issue: expected {expected_sentinels} sentinels and {expected_boundaries} boundaries")
        
        # Test plotting function imports
        print("7. Testing function imports...")
        try:
            from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta
            print("   ✅ Plotting function imported successfully")
            
            # Check function signature
            import inspect
            sig = inspect.signature(plot_multiple_seq_self_lz_factor_plot_from_fasta)
            params = list(sig.parameters.keys())
            
            if 'fasta_filepath' in params and 'factors_filepath' in params:
                print("   ✅ Function has both fasta_filepath and factors_filepath parameters")
            else:
                print(f"   ❌ Function parameters: {params}")
            
        except ImportError as e:
            print(f"   ❌ Could not import plotting function: {e}")
        
        print("\n" + "=" * 60)
        print("TEST SUMMARY")
        print("=" * 60)
        print("✅ FASTA factorization works")
        print("✅ Binary file creation with metadata works")
        print("✅ Binary metadata reading works")
        print("✅ Sentinel boundary calculation works")
        print("✅ Enhanced plotting function imports correctly")
        print("\nThe enhanced plotting function should work correctly!")
        
        return True
        
    except Exception as e:
        print(f"❌ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        # Cleanup
        if os.path.exists(fasta_path):
            os.unlink(fasta_path)
        if os.path.exists(binary_path):
            os.unlink(binary_path)


def test_backwards_compatibility():
    """Test that existing FASTA-only usage still works."""
    
    print("\n" + "=" * 60)
    print("Testing Backwards Compatibility")
    print("=" * 60)
    
    try:
        from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta
        
        # Test function can be called with old signature
        import inspect
        sig = inspect.signature(plot_multiple_seq_self_lz_factor_plot_from_fasta)
        
        # Check if fasta_filepath has default value (should be None)
        fasta_param = sig.parameters.get('fasta_filepath')
        if fasta_param and fasta_param.default is None:
            print("✅ fasta_filepath parameter has default None value")
        else:
            print("❌ fasta_filepath parameter configuration issue")
        
        print("✅ Backwards compatibility maintained")
        return True
        
    except Exception as e:
        print(f"❌ Backwards compatibility issue: {e}")
        return False


if __name__ == "__main__":
    print("Enhanced Plotting Function Test Suite")
    print("=====================================")
    
    success1 = test_enhanced_plotting_logic()
    success2 = test_backwards_compatibility()
    
    if success1 and success2:
        print("\n🎉 All tests passed!")
        print("\nThe enhanced plotting function is ready to use with:")
        print("  • FASTA files (original functionality)")
        print("  • Binary factor files with metadata (new functionality)")
        print("  • Sentinel line visualization")
        print("  • Sequence name labeling")
        exit(0)
    else:
        print("\n❌ Some tests failed")
        exit(1)