# noLZSS Scripts

This directory contains utility scripts for working with noLZSS.

## migrate_binary_format.py

Migration script to convert old noLZSS binary files (v1 format with header at beginning) to the new v2 format (with footer at end).

### Why Migrate?

The v2 format moves metadata from the beginning to the end of the file. This allows:
- Writing factors directly to the output file without buffering in memory
- More efficient processing of large files
- Streaming factor output during computation

### Usage

**Convert to new file:**
```bash
python scripts/migrate_binary_format.py old_file.bin new_file.bin
```

**In-place conversion (with backup):**
```bash
python scripts/migrate_binary_format.py --in-place file.bin
# Creates file.bin.bak and overwrites file.bin
```

**In-place conversion (no backup):**
```bash
python scripts/migrate_binary_format.py --in-place --no-backup file.bin
```

### Format Differences

**V1 Format (old):**
```
[Header: magic(8) + metadata(32) + sequence_names + sentinel_indices]
[Factors: array of Factor structs (24 bytes each)]
```

**V2 Format (new):**
```
[Factors: array of Factor structs (24 bytes each)]
[Sequence names (null-terminated strings)]
[Sentinel indices (uint64 array)]
[Footer: magic(8) + metadata(32)]
```

The magic string changes from "noLZSSv1" to "noLZSSv2" to distinguish formats.

### Notes

- The script preserves all metadata (factors, sequence names, sentinel indices)
- In-place migration creates a `.bak` backup by default
- Both formats store factors as 24-byte structs (3 × uint64_t: start, length, ref)
- The v2 format is not compatible with v1 readers; migrate all files consistently
