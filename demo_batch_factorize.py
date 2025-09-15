#!/usr/bin/env python3
"""
Demo script showing batch FASTA factorization capabilities.

This script demonstrates the main features of the batch_factorize module.
"""

import tempfile
from pathlib import Path
import sys
import os

# Add the source directory to path so we can import
sys.path.insert(0, str(Path(__file__).parent / "src"))

try:
    from noLZSS.genomics import batch_factorize
    print("✓ Successfully imported batch_factorize module")
except ImportError as e:
    print(f"✗ Failed to import batch_factorize: {e}")
    sys.exit(1)

def create_demo_files(demo_dir: Path):
    """Create sample FASTA files for demonstration."""
    print(f"\n📁 Creating demo files in {demo_dir}")
    
    # Create sample FASTA files
    files = {
        "sample1.fasta": """>chromosome_1
ATCGATCGATCGATCGATCGATCG
>chromosome_2
GCTAGCTAGCTAGCTAGCTAGCTA
""",
        "sample2.fasta": """>gene_alpha
AAAATTTGGGCCCAAAATTTGGGCCC
>gene_beta
TGCATGCATGCATGCATGCA
""",
        "sample3.fasta": """>sequence_X
ACGTACGTACGTACGT
>sequence_Y
TGCATGCATGCATGCA
>sequence_Z
GGCCGGCCGGCCGGCC
""",
    }
    
    file_paths = []
    for filename, content in files.items():
        file_path = demo_dir / filename
        with open(file_path, 'w') as f:
            f.write(content)
        file_paths.append(file_path)
        print(f"  📄 Created {filename}")
    
    # Create file list
    file_list = demo_dir / "file_list.txt"
    with open(file_list, 'w') as f:
        for file_path in file_paths:
            f.write(str(file_path) + "\n")
    print(f"  📝 Created file_list.txt with {len(file_paths)} files")
    
    return file_paths, file_list

def run_demo():
    """Run the batch factorization demo."""
    print("🧬 noLZSS Batch FASTA Factorization Demo")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory(prefix="nolzss_demo_") as temp_dir:
        demo_dir = Path(temp_dir)
        
        # Create demo files
        file_paths, file_list = create_demo_files(demo_dir)
        
        # Demo 1: Factorize using file list with both modes
        print(f"\n🧪 Demo 1: Processing file list with both modes")
        output_dir1 = demo_dir / "output_both"
        
        try:
            results = batch_factorize.process_file_list(
                file_list=[str(f) for f in file_paths],
                output_dir=output_dir1,
                mode=batch_factorize.FactorizationMode.BOTH,
                skip_existing=False
            )
            
            print(f"  ✓ Processed {len(results)} files")
            for file_path, result in results.items():
                if "error" not in result:
                    modes = [mode for mode, success in result.items() if success]
                    print(f"    📁 {Path(file_path).name}: {', '.join(modes)}")
            
            # Show output structure
            print(f"\n  📂 Output structure:")
            for mode_dir in output_dir1.iterdir():
                if mode_dir.is_dir():
                    files = list(mode_dir.glob("*.bin"))
                    print(f"    {mode_dir.name}/: {len(files)} files")
                    for file in files:
                        size = file.stat().st_size
                        print(f"      📄 {file.name} ({size} bytes)")
                        
        except Exception as e:
            print(f"  ✗ Demo 1 failed: {e}")
        
        # Demo 2: Single file with reverse complement only
        print(f"\n🧪 Demo 2: Single file with reverse complement mode")
        output_dir2 = demo_dir / "output_rc_only"
        
        try:
            results = batch_factorize.process_file_list(
                file_list=[str(file_paths[0])],
                output_dir=output_dir2,
                mode=batch_factorize.FactorizationMode.REVERSE_COMPLEMENT,
                skip_existing=False
            )
            
            print(f"  ✓ Processed {Path(file_paths[0]).name}")
            
            # Show binary file content info
            bin_files = list(output_dir2.rglob("*.bin"))
            for bin_file in bin_files:
                size = bin_file.stat().st_size
                factors_count = size // 24  # Each factor is 24 bytes
                print(f"    📄 {bin_file.name}: {factors_count} factors ({size} bytes)")
                
        except Exception as e:
            print(f"  ✗ Demo 2 failed: {e}")
        
        # Demo 3: Show utility functions
        print(f"\n🧪 Demo 3: Utility functions")
        
        print("  🔗 URL detection:")
        test_urls = [
            "https://example.com/sequence.fasta",
            "/local/path/file.fasta",
            "http://ftp.site.org/data.fa",
            "relative_file.fasta"
        ]
        
        for url in test_urls:
            is_url_result = batch_factorize.is_url(url)
            print(f"    {url}: {'URL' if is_url_result else 'Local'}")
        
        print("  ✅ FASTA validation:")
        for file_path in file_paths[:2]:  # Test first two files
            is_valid = batch_factorize.validate_fasta_file(file_path)
            print(f"    {file_path.name}: {'Valid' if is_valid else 'Invalid'}")
        
        # Demo 4: Show command line usage
        print(f"\n📚 Command Line Usage Examples:")
        print("  # Process file list with both modes:")
        print("  python -m noLZSS.genomics.batch_factorize --file-list files.txt --output-dir results --mode both")
        print()
        print("  # Process individual files:")
        print("  python -m noLZSS.genomics.batch_factorize file1.fasta file2.fasta --output-dir results --mode reverse_complement")
        print()
        print("  # Process remote files:")
        print("  python -m noLZSS.genomics.batch_factorize --file-list remote_urls.txt --output-dir results --download-dir downloads")
        print()
        print("  # Force overwrite existing files:")
        print("  python -m noLZSS.genomics.batch_factorize --file-list files.txt --output-dir results --force")
        
        print(f"\n✨ Demo completed successfully!")
        print(f"   All demo files were created in: {demo_dir}")
        print(f"   (Temporary directory will be cleaned up automatically)")

if __name__ == "__main__":
    run_demo()