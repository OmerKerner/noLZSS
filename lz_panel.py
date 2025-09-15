#!/usr/bin/env python3
"""
Interactive LZ Factor Plot using Datashader/Panel

A high-performance, interactive visualization for non-overlapping Lempel-Ziv-Storer-Szymanski (LZSS) 
factorization results from DNA FASTA files. Uses Datashader for scalable rendering of millions 
of factors and Panel for interactive controls.

Usage:
    python lz_panel.py --fasta path/to/file.fasta
    panel serve lz_panel.py --show

Dependencies:
    numpy pandas holoviews bokeh panel datashader colorcet

Performance Notes:
    - Base (datashaded) layer handles ≥5M segments smoothly
    - Overlay recompute target ≤100ms for ≤2000 bins  
    - All filtering & decimation in vectorized NumPy/Pandas
"""

import argparse
import sys
from pathlib import Path

try:
    import panel as pn
    from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta
except ImportError as e:
    print(f"Error: Missing required dependencies: {e}")
    print("Install with: pip install 'noLZSS[panel]'")
    sys.exit(1)


def create_demo_app():
    """Create a demo app with synthetic data when no FASTA is provided."""
    import numpy as np
    import pandas as pd
    import holoviews as hv
    import datashader as ds
    from holoviews.operation.datashader import datashade, dynspread
    
    hv.extension('bokeh')
    pn.extension()
    
    # Generate synthetic factor data
    np.random.seed(42)
    n_factors = 1000
    
    starts = np.sort(np.random.randint(0, 10000, n_factors))
    lengths = np.random.randint(5, 100, n_factors)
    refs = np.random.randint(0, 10000, n_factors)
    is_rcs = np.random.choice([True, False], n_factors, p=[0.3, 0.7])
    
    # Build coordinates
    x0_vals = starts
    x1_vals = starts + lengths
    y0_vals = np.where(is_rcs, refs + lengths, refs)
    y1_vals = np.where(is_rcs, refs, refs + lengths)
    
    df = pd.DataFrame({
        'x0': x0_vals,
        'y0': y0_vals, 
        'x1': x1_vals,
        'y1': y1_vals,
        'length': lengths,
        'dir': is_rcs.astype(int),
        'start': starts,
        'ref': refs,
        'end': x1_vals,
        'is_rc': is_rcs
    })
    
    # Create segments for forward and reverse complement
    df_fwd = df[df['dir'] == 0]
    df_rc = df[df['dir'] == 1]
    
    segments_fwd = hv.Segments(df_fwd, kdims=['x0','y0','x1','y1'], vdims=['length'])
    segments_rc = hv.Segments(df_rc, kdims=['x0','y0','x1','y1'], vdims=['length'])
    
    # Apply datashader
    shaded_fwd = dynspread(datashade(segments_fwd, aggregator=ds.max('length'), cmap=['white', 'blue']))
    shaded_rc = dynspread(datashade(segments_rc, aggregator=ds.max('length'), cmap=['white', 'red']))
    
    # Add diagonal
    max_val = max(df[['x1', 'y1']].max())
    min_val = min(df[['x0', 'y0']].min())
    diagonal = hv.Curve([(min_val, min_val), (max_val, max_val)]).opts(
        line_dash='dashed', line_color='gray', line_width=1, alpha=0.5
    )
    
    plot = (shaded_fwd * shaded_rc * diagonal).opts(
        width=800, height=800, 
        xlabel='Position in sequence',
        ylabel='Reference position', 
        title='Demo LZ Factor Plot (Synthetic Data)'
    )
    
    info_panel = pn.pane.Markdown(f"""
    ## Demo Mode
    
    Showing **{len(df)} synthetic factors** to demonstrate the interactive visualization.
    
    - **Blue lines**: Forward factors 
    - **Red lines**: Reverse complement factors
    - **Gray dashed**: y=x diagonal reference
    
    To use with real data:
    ```bash
    python lz_panel.py --fasta your_file.fasta
    ```
    """)
    
    return pn.Row(info_panel, plot)


def main():
    parser = argparse.ArgumentParser(
        description="Interactive LZ Factor Plot using Datashader/Panel"
    )
    parser.add_argument(
        "--fasta", 
        type=str,
        help="Path to FASTA file containing DNA sequences"
    )
    parser.add_argument(
        "--name",
        type=str,
        help="Optional name for the plot (defaults to filename)"
    )
    parser.add_argument(
        "--port",
        type=int,
        default=5007,
        help="Port for Panel server (default: 5007)"
    )
    
    args = parser.parse_args()
    
    if args.fasta:
        fasta_path = Path(args.fasta)
        if not fasta_path.exists():
            print(f"Error: FASTA file not found: {fasta_path}")
            sys.exit(1)
        
        print(f"Creating interactive plot for {fasta_path}...")
        
        try:
            app = plot_multiple_seq_self_lz_factor_plot_from_fasta(
                fasta_filepath=fasta_path,
                name=args.name,
                show_plot=False,
                return_panel=True
            )
        except Exception as e:
            print(f"Error creating plot: {e}")
            sys.exit(1)
    else:
        print("No FASTA file specified, creating demo with synthetic data...")
        app = create_demo_app()
    
    print(f"Starting Panel server on port {args.port}...")
    print(f"Open your browser to: http://localhost:{args.port}")
    
    pn.serve(app, show=True, port=args.port)


if __name__ == "__main__":
    main()