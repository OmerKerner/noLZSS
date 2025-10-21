#!/usr/bin/env python3
"""
Migration script to convert old noLZSS binary files (v1 with header) to new format (v2 with footer).

This script reads binary files created with the old header-at-beginning format and
converts them to the new footer-at-end format.

Usage:
    python scripts/migrate_binary_format.py old_file.bin new_file.bin
    python scripts/migrate_binary_format.py --in-place file.bin
"""

import struct
import sys
import argparse
from pathlib import Path
from typing import Tuple, List, Optional


def read_v1_format(filepath: Path) -> Tuple[bytes, dict]:
    """
    Read an old v1 format binary file (header at beginning).
    
    Returns:
        Tuple of (factors_data, metadata_dict)
    """
    with open(filepath, 'rb') as f:
        # Read basic header (40 bytes)
        header_data = f.read(40)
        if len(header_data) < 40:
            raise ValueError(f"File too small: {len(header_data)} bytes")
        
        magic = header_data[:8]
        if magic != b'noLZSSv1':
            raise ValueError(f"Not a v1 format file (magic: {magic})")
        
        num_factors, num_sequences, num_sentinels, header_size = struct.unpack('<QQQQ', header_data[8:40])
        
        # Seek back and read full header
        f.seek(0)
        full_header = f.read(header_size)
        
        # Parse metadata from header
        offset = 40
        
        # Read sequence names
        sequence_names = []
        for _ in range(num_sequences):
            name_start = offset
            while offset < len(full_header) and full_header[offset] != 0:
                offset += 1
            if offset >= len(full_header):
                break
            name = full_header[name_start:offset].decode('utf-8')
            sequence_names.append(name)
            offset += 1
        
        # Read sentinel indices
        sentinel_indices = []
        for _ in range(num_sentinels):
            if offset + 8 > len(full_header):
                break
            idx = struct.unpack('<Q', full_header[offset:offset+8])[0]
            sentinel_indices.append(idx)
            offset += 8
        
        # Read factors
        factors_data = f.read()
        
        metadata = {
            'num_factors': num_factors,
            'num_sequences': num_sequences,
            'num_sentinels': num_sentinels,
            'sequence_names': sequence_names,
            'sentinel_indices': sentinel_indices
        }
        
        return factors_data, metadata


def write_v2_format(filepath: Path, factors_data: bytes, metadata: dict) -> None:
    """
    Write a new v2 format binary file (footer at end).
    """
    with open(filepath, 'wb') as f:
        # Write factors first
        f.write(factors_data)
        
        # Write sequence names
        for name in metadata['sequence_names']:
            f.write(name.encode('utf-8'))
            f.write(b'\0')
        
        # Write sentinel indices
        for idx in metadata['sentinel_indices']:
            f.write(struct.pack('<Q', idx))
        
        # Calculate footer size
        names_size = sum(len(name) + 1 for name in metadata['sequence_names'])
        sentinels_size = len(metadata['sentinel_indices']) * 8
        footer_size = 48 + names_size + sentinels_size  # 48 bytes for footer struct
        
        # Calculate total_length from factors_data
        total_length = 0
        num_factors = len(factors_data) // 24  # Each factor is 24 bytes
        for i in range(num_factors):
            offset = i * 24
            # Factor is (start, length, ref) each uint64_t
            length = struct.unpack('<Q', factors_data[offset+8:offset+16])[0]
            total_length += length
        
        # Write footer
        magic = b'noLZSSv2'
        footer = struct.pack('<8sQQQQQ',
                           magic,
                           metadata['num_factors'],
                           metadata['num_sequences'],
                           metadata['num_sentinels'],
                           footer_size,
                           total_length)
        f.write(footer)


def migrate_file(input_path: Path, output_path: Path, backup: bool = True) -> None:
    """
    Migrate a single file from v1 to v2 format.
    
    Args:
        input_path: Path to input v1 file
        output_path: Path to output v2 file
        backup: If True and doing in-place migration, create a .bak file
    """
    print(f"Reading v1 format from: {input_path}")
    factors_data, metadata = read_v1_format(input_path)
    
    print(f"  Found {metadata['num_factors']} factors")
    print(f"  Found {metadata['num_sequences']} sequences")
    print(f"  Found {metadata['num_sentinels']} sentinels")
    
    # If in-place and backup requested, create backup
    if input_path == output_path and backup:
        backup_path = input_path.with_suffix(input_path.suffix + '.bak')
        print(f"Creating backup: {backup_path}")
        import shutil
        shutil.copy2(input_path, backup_path)
    
    print(f"Writing v2 format to: {output_path}")
    write_v2_format(output_path, factors_data, metadata)
    print("✓ Migration complete")


def main():
    parser = argparse.ArgumentParser(
        description='Migrate noLZSS binary files from v1 (header) to v2 (footer) format',
        epilog='Example: python scripts/migrate_binary_format.py old.bin new.bin'
    )
    parser.add_argument('input', type=Path, help='Input v1 binary file')
    parser.add_argument('output', type=Path, nargs='?', help='Output v2 binary file (optional for --in-place)')
    parser.add_argument('--in-place', action='store_true', help='Migrate file in-place (overwrite input)')
    parser.add_argument('--no-backup', action='store_true', help='Skip creating .bak file for in-place migration')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.in_place:
        if args.output:
            parser.error("Cannot specify output file with --in-place")
        output_path = args.input
    else:
        if not args.output:
            parser.error("Output file required (or use --in-place)")
        output_path = args.output
    
    # Check input exists
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Migrate
    try:
        migrate_file(args.input, output_path, backup=not args.no_backup)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
