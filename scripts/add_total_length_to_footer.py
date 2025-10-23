#!/usr/bin/env python3
"""
Migration script to add total_length field to existing noLZSS v2 binary files.

This script updates binary files that have the old 40-byte footer to the new 48-byte
footer format that includes the total_length field.

Usage:
    python scripts/add_total_length_to_footer.py old_file.bin new_file.bin
    python scripts/add_total_length_to_footer.py --in-place file.bin
"""

import struct
import sys
import argparse
from pathlib import Path
from typing import Tuple, Dict


def read_old_v2_footer(filepath: Path) -> Tuple[bytes, Dict, int]:
    """
    Read an old v2 format binary file (40-byte footer without total_length).
    
    Returns:
        Tuple of (factors_data, metadata_dict, total_length_calculated)
    """
    with open(filepath, 'rb') as f:
        # Try to read old format footer (40 bytes)
        f.seek(-40, 2)
        footer_data = f.read(40)
        if len(footer_data) < 40:
            raise ValueError(f"File too small: {len(footer_data)} bytes")
        
        magic = footer_data[:8]
        if magic != b'noLZSSv2':
            raise ValueError(f"Not a v2 format file (magic: {magic})")
        
        num_factors, num_sequences, num_sentinels, footer_size = struct.unpack('<QQQQ', footer_data[8:40])
        
        # Check if this is already the new format (48 bytes)
        if footer_size >= 48:
            # Try reading with new format
            f.seek(-48, 2)
            new_footer_data = f.read(48)
            if len(new_footer_data) == 48:
                new_magic = new_footer_data[:8]
                if new_magic == b'noLZSSv2':
                    # This is already new format
                    raise ValueError("File already has total_length field (new format)")
        
        # Read full footer
        f.seek(-footer_size, 2)
        full_footer = f.read(footer_size)
        
        # Parse metadata from footer
        metadata_size = footer_size - 40
        metadata = full_footer[:metadata_size]
        offset = 0
        
        # Read sequence names
        sequence_names = []
        for _ in range(num_sequences):
            name_start = offset
            while offset < len(metadata) and metadata[offset] != 0:
                offset += 1
            if offset >= len(metadata):
                break
            name = metadata[name_start:offset].decode('utf-8')
            sequence_names.append(name)
            offset += 1
        
        # Read sentinel indices
        sentinel_indices = []
        for _ in range(num_sentinels):
            if offset + 8 > len(metadata):
                break
            idx = struct.unpack('<Q', metadata[offset:offset+8])[0]
            sentinel_indices.append(idx)
            offset += 8
        
        # Read factors and calculate total_length
        factors_start = 0
        factors_end = len(full_footer) - footer_size
        f.seek(0)
        factors_data = f.read(num_factors * 24)
        
        total_length = 0
        for i in range(num_factors):
            factor_offset = i * 24
            # Factor is (start, length, ref) each uint64_t
            length = struct.unpack('<Q', factors_data[factor_offset+8:factor_offset+16])[0]
            total_length += length
        
        metadata_dict = {
            'num_factors': num_factors,
            'num_sequences': num_sequences,
            'num_sentinels': num_sentinels,
            'sequence_names': sequence_names,
            'sentinel_indices': sentinel_indices
        }
        
        return factors_data, metadata_dict, total_length


def write_new_v2_format(filepath: Path, factors_data: bytes, metadata: Dict, total_length: int) -> None:
    """
    Write a new v2 format binary file with total_length field (48-byte footer).
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
        
        # Calculate new footer size (48 bytes for struct + metadata)
        names_size = sum(len(name) + 1 for name in metadata['sequence_names'])
        sentinels_size = len(metadata['sentinel_indices']) * 8
        footer_size = 48 + names_size + sentinels_size
        
        # Write new footer with total_length field
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
    Migrate a single file from old v2 format (40-byte footer) to new v2 format (48-byte footer).
    
    Args:
        input_path: Path to input old v2 file
        output_path: Path to output new v2 file
        backup: If True and doing in-place migration, create a .bak file
    """
    print(f"Reading old v2 format from: {input_path}")
    try:
        factors_data, metadata, total_length = read_old_v2_footer(input_path)
    except ValueError as e:
        print(f"Error: {e}")
        print("This script only works on old v2 format files (40-byte footer).")
        sys.exit(1)
    
    print(f"  Found {metadata['num_factors']} factors")
    print(f"  Found {metadata['num_sequences']} sequences")
    print(f"  Found {metadata['num_sentinels']} sentinels")
    print(f"  Calculated total_length: {total_length}")
    
    # If in-place and backup requested, create backup
    if input_path == output_path and backup:
        backup_path = input_path.with_suffix(input_path.suffix + '.bak')
        print(f"Creating backup: {backup_path}")
        import shutil
        shutil.copy2(input_path, backup_path)
    
    print(f"Writing new v2 format with total_length to: {output_path}")
    write_new_v2_format(output_path, factors_data, metadata, total_length)
    print("✓ Migration complete - added total_length field")


def main():
    parser = argparse.ArgumentParser(
        description='Add total_length field to noLZSS v2 binary files',
        epilog='Example: python scripts/add_total_length_to_footer.py old.bin new.bin'
    )
    parser.add_argument('input', type=Path, help='Input v2 binary file (old format)')
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
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
