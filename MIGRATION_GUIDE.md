# Binary Format Migration Guide: v1 to v2

## Overview

noLZSS v0.2.2 introduces a new binary file format (v2) that moves metadata from the beginning (header) to the end (footer) of binary factor files. This change enables more efficient writing by allowing factors to be written directly without first knowing the total count.

## What Changed?

### Old Format (v1 - Header at Beginning)
```
Bytes:
0-7:    Magic string "noLZSSv1"
8-15:   num_factors (uint64)
16-23:  num_sequences (uint64)
24-31:  num_sentinels (uint64)
32-39:  header_size (uint64)
40+:    Sequence names (null-terminated strings)
X+:     Sentinel indices (uint64 array)
Y+:     Factors (24 bytes each: start, length, ref)
```

### New Format (v2 - Footer at End)
```
Bytes:
0+:     Factors (24 bytes each: start, length, ref)
X+:     Sequence names (null-terminated strings)
Y+:     Sentinel indices (uint64 array)
Z-39:   num_factors (uint64)
Z-31:   num_sequences (uint64)
Z-23:   num_sentinels (uint64)
Z-15:   footer_size (uint64)
Z-7:    Magic string "noLZSSv2"
```

## Why This Change?

### Benefits of Footer Format

1. **Streaming Output**: Factors can be written directly as they're computed, without buffering in memory
2. **Memory Efficiency**: No need to store all factors before writing metadata
3. **Performance**: Reduced memory footprint for large-scale factorizations
4. **Scalability**: Better suited for processing very large genomic datasets

### Backward Compatibility

⚠️ **Important**: The v2 format is NOT backward compatible with v1 readers. Files written with v2 cannot be read by older versions of noLZSS.

## Migration Instructions

### For Users

If you have existing binary files in v1 format, use the provided migration script:

```bash
# Convert individual files
python scripts/migrate_binary_format.py old_file.bin new_file.bin

# In-place conversion (creates backup)
python scripts/migrate_binary_format.py --in-place file.bin

# Batch conversion
for f in *.bin; do
    python scripts/migrate_binary_format.py --in-place "$f"
done
```

### For Developers

If you're directly reading binary files (not using noLZSS utilities):

**Old code (v1):**
```python
with open(filepath, 'rb') as f:
    header = f.read(40)
    magic = header[:8]  # Should be b'noLZSSv1'
    num_factors = struct.unpack('<Q', header[8:16])[0]
    # Read factors after header...
```

**New code (v2):**
```python
with open(filepath, 'rb') as f:
    # Seek to end to read footer
    f.seek(-40, 2)
    footer = f.read(40)
    magic = footer[:8]  # Should be b'noLZSSv2'
    num_factors = struct.unpack('<Q', footer[8:16])[0]
    
    # Seek to beginning to read factors
    f.seek(0)
    # Read factors...
```

## Code Changes Summary

### C++ Changes

**Modified files:**
- `src/cpp/factorizer.hpp`: Renamed `FactorFileHeader` to `FactorFileFooter` with v2 magic
- `src/cpp/factorizer.cpp`: Updated all `write_factors_binary_file*` functions to write footer
- `src/cpp/fasta_processor.cpp`: Updated FASTA binary write functions to write footer

**Key changes:**
- All functions now write factors first, then metadata at the end
- Footer struct includes constructor to ensure proper initialization
- Magic string changed from "noLZSSv1" to "noLZSSv2"

### Python Changes

**Modified files:**
- `src/noLZSS/utils.py`: Updated `read_factors_binary_file` and `read_factors_binary_file_with_metadata`

**Key changes:**
- Both functions now seek to end of file to read footer first
- Footer is read using negative seek offset: `f.seek(-40, 2)`
- Magic string check updated to expect "noLZSSv2"

### Test Changes

**Modified files:**
- `tests/test_utils.py`: Updated mock binary data to use footer format
- `tests/test_cpp_bindings.py`: Updated file size verification logic

**Key changes:**
- Test files now create footer at end instead of header at beginning
- File size calculations updated to account for footer position

## FAQ

### Q: Will old binary files still work?

No, v1 binary files will not be readable by v2 code. You must migrate them using the provided script.

### Q: Can I convert v2 files back to v1?

The migration script only supports v1→v2 conversion. Create a reverse script if needed by swapping the read/write logic.

### Q: What happens if I try to read a v1 file with v2 code?

You'll get an error: "Invalid file format: missing noLZSS magic footer (expected v2 format)"

### Q: Are there any performance differences for reading?

Reading performance is essentially identical. The footer must be read first (one seek to end), then factors are read from the beginning.

### Q: How do I identify file format version?

Check the magic string:
- v1: First 8 bytes contain "noLZSSv1"
- v2: Last 40 bytes (footer) start with "noLZSSv2"

## Version Compatibility Matrix

| noLZSS Version | Creates | Can Read |
|----------------|---------|----------|
| ≤ 0.2.1        | v1      | v1       |
| ≥ 0.2.2        | v2      | v2       |

## Additional Resources

- Migration script: `scripts/migrate_binary_format.py`
- Script documentation: `scripts/README.md`
- Test examples: `tests/test_cpp_bindings.py`, `tests/test_utils.py`

## Support

For issues or questions about migration:
1. Check the migration script's `--help` output
2. Review test files for usage examples
3. Open an issue on GitHub with the "migration" label
