#!/usr/bin/env python3
"""
Example usage of the enhanced plot_multiple_seq_self_lz_factor_plot_from_fasta function.

This example demonstrates:
1. Original usage with FASTA files
2. New usage with binary factor files
3. The sentinel visualization features
"""

import tempfile
import os
import sys

# Set up path to use development version
sys.path.insert(0, 'src')

import noLZSS._noLZSS as cpp


def create_example_files():
    """Create example FASTA and binary files for demonstration."""
    
    # Create a multi-sequence FASTA file
    fasta_content = """>Human_Chr21_Fragment
ACGTACGTACGTACGT
>Mouse_Chr19_Fragment  
TGCATGCATGCATGCA
>Zebrafish_Chr3_Fragment
AAAATTTTGGGGCCCC
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_path = f.name
    
    # Create corresponding binary file with metadata
    binary_path = fasta_path.replace('.fasta', '.bin')
    cpp.write_factors_binary_file_fasta_multiple_dna_w_rc(fasta_path, binary_path)
    
    return fasta_path, binary_path


def example_fasta_usage():
    """Example 1: Using the original FASTA input (backwards compatible)."""
    
    print("Example 1: Using FASTA file input (original functionality)")
    print("-" * 60)
    
    fasta_path, binary_path = create_example_files()
    
    try:
        from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta
        
        print(f"FASTA file: {fasta_path}")
        print("Calling plot_multiple_seq_self_lz_factor_plot_from_fasta with FASTA...")
        
        # Original usage - backwards compatible
        # Note: This would show the plot if Panel dependencies were available
        result = plot_multiple_seq_self_lz_factor_plot_from_fasta(
            fasta_filepath=fasta_path,
            name="Example_Multi_Species",
            show_plot=False,  # Don't show since we might not have display
            return_panel=True
        )
        
        print("✅ FASTA plotting function completed successfully!")
        print("   Features enabled:")
        print("   • Multiple sequence factorization")
        print("   • Sentinel boundary detection")
        print("   • Red lines at sequence boundaries")
        print("   • Sequence names on axes")
        
    except ImportError as e:
        print(f"⚠️  Plotting dependencies not available: {e}")
        print("   Install with: pip install 'noLZSS[panel]'")
    except Exception as e:
        print(f"❌ Error in FASTA plotting: {e}")
    
    finally:
        os.unlink(fasta_path)
        os.unlink(binary_path)


def example_binary_usage():
    """Example 2: Using the new binary factor file input."""
    
    print("\nExample 2: Using binary factor file input (new functionality)")
    print("-" * 60)
    
    fasta_path, binary_path = create_example_files()
    
    try:
        from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta
        
        print(f"Binary file: {binary_path}")
        print("Calling plot_multiple_seq_self_lz_factor_plot_from_fasta with binary factors...")
        
        # New usage - binary factor file input
        result = plot_multiple_seq_self_lz_factor_plot_from_fasta(
            factors_filepath=binary_path,  # New parameter!
            name="Example_From_Binary",
            show_plot=False,
            return_panel=True
        )
        
        print("✅ Binary plotting function completed successfully!")
        print("   New features:")
        print("   • Reads pre-computed factors from binary files")
        print("   • Faster than re-factorizing FASTA")
        print("   • Includes sequence names from FASTA headers")
        print("   • Preserves all sentinel information")
        
    except ImportError as e:
        print(f"⚠️  Plotting dependencies not available: {e}")
        print("   Install with: pip install 'noLZSS[panel]'")
    except Exception as e:
        print(f"❌ Error in binary plotting: {e}")
    
    finally:
        os.unlink(fasta_path)
        os.unlink(binary_path)


def show_plotting_enhancements():
    """Show what enhancements have been added to the plotting."""
    
    print("\nPlotting Enhancements Summary")
    print("=" * 60)
    
    print("🆕 NEW FEATURES:")
    print("  1. Binary factor file input support")
    print("     - Faster plotting (no re-factorization needed)")
    print("     - Preserves sequence metadata")
    print("     - Compatible with batch processing workflows")
    
    print("\n  2. Sentinel boundary visualization")
    print("     - Red vertical lines at sequence boundaries")
    print("     - Red horizontal lines at sequence boundaries") 
    print("     - Clear visual separation between sequences")
    
    print("\n  3. Sequence name labeling")
    print("     - Sequence names displayed on X and Y axes")
    print("     - Names extracted from FASTA headers")
    print("     - Enhanced plot titles with sequence information")
    
    print("\n🔄 BACKWARDS COMPATIBILITY:")
    print("  - Original FASTA-only usage still works")
    print("  - Function signature extended, not changed")
    print("  - Existing code continues to work unchanged")
    
    print("\n📋 FUNCTION SIGNATURE:")
    print("  plot_multiple_seq_self_lz_factor_plot_from_fasta(")
    print("      fasta_filepath=None,      # Original parameter")
    print("      factors_filepath=None,    # NEW: Binary factors input")
    print("      name=None,")
    print("      save_path=None,")
    print("      show_plot=True,")
    print("      return_panel=False")
    print("  )")


if __name__ == "__main__":
    print("Enhanced Plotting Function Usage Examples")
    print("=========================================")
    
    example_fasta_usage()
    example_binary_usage()
    show_plotting_enhancements()
    
    print("\n🎯 NEXT STEPS:")
    print("  1. Install Panel dependencies: pip install 'noLZSS[panel]'")
    print("  2. Use the function with your own FASTA or binary files")
    print("  3. Enjoy enhanced visualization with sequence boundaries!")