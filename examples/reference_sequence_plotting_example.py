#!/usr/bin/env python3
"""
Example usage of reference sequence plotting functions.

This script demonstrates how to use the new plotting functions:
- plot_reference_seq_lz_factor_plot_simple()  # matplotlib-based
- plot_reference_seq_lz_factor_plot()         # interactive Panel/Datashader
- Command-line interface via python -m noLZSS.genomics.plots reference-plot

These functions are designed to work with the output of factorize_w_reference_seq().
"""

# Example usage (requires the package to be built with C++ extension)

def example_cli_usage():
    """Example using the command-line interface"""
    
    print("""
    Command-Line Interface Usage:
    ============================
    
    The plotting functions can be accessed via command-line:
    
    # Simple matplotlib plot
    python -m noLZSS.genomics.plots reference-plot \\
        "ATCGATCGATCGATCG" \\
        "ATCGCCCCGATCGAAA" \\
        --reference_name "Reference Genome" \\
        --target_name "Query Sequence" \\
        --save_path "reference_plot.png" \\
        --show_plot
    
    # Interactive Panel plot
    python -m noLZSS.genomics.plots reference-plot \\
        "ATCGATCGATCGATCG" \\
        "ATCGCCCCGATCGAAA" \\
        --reference_name "Reference Genome" \\
        --target_name "Query Sequence" \\
        --save_path "reference_plot_interactive.png" \\
        --interactive \\
        --show_plot
    
    # Using pre-computed factors from binary file
    python -m noLZSS.genomics.plots reference-plot \\
        "ATCGATCGATCGATCG" \\
        "ATCGCCCCGATCGAAA" \\
        --factors_filepath "reference_factors.bin" \\
        --save_path "reference_plot.png"
    
    # View all available plot types
    python -m noLZSS.genomics.plots --help
    
    # View help for reference-plot
    python -m noLZSS.genomics.plots reference-plot --help
    """)

def example_simple_plot():
    """Example using the simple matplotlib-based plotting function"""
    
    from noLZSS.genomics.sequences import factorize_w_reference_seq
    from noLZSS.genomics.plots import plot_reference_seq_lz_factor_plot_simple
    
    # Define reference and target sequences
    reference_seq = "ATCGATCGATCGATCG"    # 16 nucleotides
    target_seq = "ATCGCCCCGATCGAAA"      # 16 nucleotides
    
    print(f"Reference: {reference_seq}")
    print(f"Target: {target_seq}")
    
    # Factorize target sequence using reference
    factors = factorize_w_reference_seq(reference_seq, target_seq)
    print(f"Generated {len(factors)} factors")
    
    # Create and display plot
    plot_reference_seq_lz_factor_plot_simple(
        reference_seq=reference_seq,
        target_seq=target_seq,
        factors=factors,  # Optional: will compute if not provided
        reference_name="Reference Genome",
        target_name="Query Sequence",
        save_path="reference_plot_simple.png",
        show_plot=True
    )

def example_interactive_plot():
    """Example using the interactive Panel/Datashader plotting function"""
    
    from noLZSS.genomics.sequences import factorize_w_reference_seq
    from noLZSS.genomics.plots import plot_reference_seq_lz_factor_plot
    
    # Define reference and target sequences
    reference_seq = "ATCGATCGATCGATCGATCGATCGATCGATCG"  # 32 nucleotides
    target_seq = "ATCGCCCCGATCGAAAATTTTCCCCGGGG"        # 28 nucleotides
    
    print(f"Reference: {reference_seq}")
    print(f"Target: {target_seq}")
    
    # Create interactive plot
    app = plot_reference_seq_lz_factor_plot(
        reference_seq=reference_seq,
        target_seq=target_seq,
        reference_name="Reference Genome",
        target_name="Query Sequence",
        save_path="reference_plot_interactive.png",
        show_plot=True,
        return_panel=True
    )
    
    return app

def example_with_binary_file():
    """Example using a pre-computed binary factors file"""
    
    from noLZSS.genomics.sequences import factorize_w_reference_seq_file
    from noLZSS.genomics.plots import plot_reference_seq_lz_factor_plot_simple
    
    # Create sequences
    reference_seq = "ATCGATCGATCGATCG"
    target_seq = "ATCGCCCCGATCGAAA"
    
    # Save factors to binary file
    binary_file = "reference_factors.bin"
    count = factorize_w_reference_seq_file(reference_seq, target_seq, binary_file)
    print(f"Saved {count} factors to {binary_file}")
    
    # Create plot from binary file
    plot_reference_seq_lz_factor_plot_simple(
        reference_seq=reference_seq,
        target_seq=target_seq,
        factors_filepath=binary_file,  # Use binary file instead of computing
        reference_name="Reference",
        target_name="Target",
        save_path="reference_plot_from_file.png",
        show_plot=True
    )

def example_detailed_explanation():
    """Detailed explanation of the plotting output"""
    
    print("""
    Reference Sequence Factor Plot Explanation:
    ==========================================
    
    The plot shows factorization of a target sequence using a reference sequence.
    The x-axis represents positions in the concatenated sequence (reference + target).
    The y-axis represents reference positions that factors point to.
    
    Plot Elements:
    -------------
    1. Reference Region (left part):
       - Blue background tint
       - Contains the reference sequence (positions 0 to ref_length-1)
    
    2. Target Region (right part):
       - Red/orange background tint  
       - Contains the target sequence (positions ref_length+1 to end)
       - This is where all factors start (factorized region)
    
    3. Sequence Boundary:
       - Green vertical and horizontal lines
       - Separates reference and target regions
    
    4. Factor Lines:
       - Red lines: Forward factors in target region
       - Orange lines: Reverse complement factors in target region
       - Blue lines: Forward factors in reference region (rare)
       - Dark blue lines: Reverse complement factors in reference region (rare)
    
    5. Diagonal Line:
       - Gray dashed line where y=x
       - Reference for self-similarity
    
    Interpretation:
    --------------
    - Each line represents a factor (substring match)
    - Line starts at target position, points to reference position
    - Line length represents factor length
    - Color indicates orientation (forward/reverse complement)
    - Steep lines indicate distant references
    - Lines parallel to diagonal indicate nearby references
    
    This visualization helps understand:
    - How much the target sequence reuses patterns from the reference
    - Which parts of the reference are most reused
    - The orientation of matches (forward vs reverse complement)
    """)

if __name__ == "__main__":
    print("Reference Sequence Plotting Examples")
    print("=" * 40)
    
    try:
        print("\n0. Command-line interface:")
        example_cli_usage()
        
        print("\n1. Simple matplotlib plot:")
        example_simple_plot()
        
        print("\n2. Interactive Panel plot:")
        example_interactive_plot()
        
        print("\n3. Plot from binary file:")
        example_with_binary_file()
        
        print("\n4. Detailed explanation:")
        example_detailed_explanation()
        
    except ImportError as e:
        print(f"Error: {e}")
        print("\nNote: This example requires the noLZSS package to be built with C++ extension.")
        print("Install with: pip install -e .")
        print("For plotting dependencies: pip install 'noLZSS[plotting]' or pip install matplotlib")
        print("For interactive plots: pip install 'noLZSS[panel]'")