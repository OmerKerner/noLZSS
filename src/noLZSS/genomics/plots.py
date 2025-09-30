"""
FASTA file plotting utilities.

This module provides functions for creating plots and visualizations
from FASTA files and their factorizations.
"""

from typing import Union, Optional, Dict, Any, List, Tuple
from pathlib import Path
import warnings
import argparse

from ..utils import NoLZSSError
from .fasta import _parse_fasta_content
from .sequences import detect_sequence_type


class PlotError(NoLZSSError):
    """Raised when plotting operations fail."""
    pass


def plot_single_seq_accum_factors_from_file(
    fasta_filepath: Optional[Union[str, Path]] = None,
    factors_filepath: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    max_sequences: Optional[int] = None,
    save_factors_text: bool = True,
    save_factors_binary: bool = False
) -> Dict[str, Dict[str, Any]]:
    """
    Process a FASTA file or binary factors file, factorize sequences (if needed), create plots, and save results.

    For each sequence:
    - If FASTA file: reads sequences, factorizes them, and saves factor data and plots
    - If binary factors file: reads existing factors and creates plots

    Args:
        fasta_filepath: Path to input FASTA file (mutually exclusive with factors_filepath)
        factors_filepath: Path to binary factors file (mutually exclusive with fasta_filepath)
        output_dir: Directory to save all output files (required for FASTA, optional for binary)
        max_sequences: Maximum number of sequences to process (None for all)
        save_factors_text: Whether to save factors as text files (only for FASTA input)
        save_factors_binary: Whether to save factors as binary files (only for FASTA input)

    Returns:
        Dictionary with processing results for each sequence:
        {
            'sequence_id': {
                'sequence_length': int,
                'num_factors': int,
                'factors_file': str,  # path to saved factors
                'plot_file': str,     # path to saved plot
                'factors': List[Tuple[int, int, int]]  # the factors
            }
        }

    Raises:
        PlotError: If file processing fails
        FileNotFoundError: If input file doesn't exist
        ValueError: If both or neither input files are provided, or if output_dir is missing for FASTA input
    """
    from ..core import factorize, write_factors_binary_file
    from ..utils import read_factors_binary_file_with_metadata
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend for batch processing
    import re

    # Validate input arguments
    if (fasta_filepath is None) == (factors_filepath is None):
        raise ValueError("Exactly one of fasta_filepath or factors_filepath must be provided")

    # Determine input type and file path
    if fasta_filepath is not None:
        input_filepath = Path(fasta_filepath)
        input_type = "fasta"
        if output_dir is None:
            raise ValueError("output_dir is required when processing FASTA files")
        output_dir = Path(output_dir)
    else:
        if factors_filepath is None:
            raise ValueError("Either fasta_filepath or factors_filepath must be provided")
        input_filepath = Path(factors_filepath)
        input_type = "binary"
        if output_dir is None:
            output_dir = input_filepath.parent  # Default to same directory as binary file
        else:
            output_dir = Path(output_dir)

    if not input_filepath.exists():
        raise FileNotFoundError(f"Input file not found: {input_filepath}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    if input_type == "fasta":
        # Process FASTA file (original logic)
        # Read FASTA file
        sequences = _parse_fasta_content(input_filepath.read_text())

        if not sequences:
            raise PlotError("No sequences found in FASTA file")

        processed_count = 0

        for seq_id, sequence in sequences.items():
            if max_sequences is not None and processed_count >= max_sequences:
                break

            print(f"Processing sequence {seq_id} ({len(sequence)} bp)...")

            # Detect sequence type and validate
            seq_type = detect_sequence_type(sequence)

            if seq_type == 'dna':
                # Validate as nucleotide
                if not re.match(r'^[ACGT]+$', sequence.upper()):
                    invalid_chars = set(sequence.upper()) - set('ACGT')
                    print(f"  Warning: Skipping {seq_id} - contains invalid nucleotides: {invalid_chars}")
                    continue
                sequence = sequence.upper()
                print("  Detected nucleotide sequence")

            elif seq_type == 'protein':
                # Validate as amino acid
                valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
                if not all(c in valid_aa for c in sequence.upper()):
                    invalid_chars = set(sequence.upper()) - valid_aa
                    print(f"  Warning: Skipping {seq_id} - contains invalid amino acids: {invalid_chars}")
                    continue
                sequence = sequence.upper()
                print("  Detected amino acid sequence")

            else:
                print(f"  Warning: Skipping {seq_id} - unknown sequence type: {seq_type}")
                continue

            # Factorize
            try:
                factors = factorize(sequence.encode('ascii'))
                print(f"  Factorized into {len(factors)} factors")
            except Exception as e:
                print(f"  Warning: Failed to factorize {seq_id}: {e}")
                continue

            # Save factors as text
            factors_text_file = None
            if save_factors_text:
                factors_text_file = output_dir / f"factors_{seq_id}.txt"
                try:
                    with open(factors_text_file, 'w') as f:
                        f.write(f"Sequence: {seq_id}\n")
                        f.write(f"Length: {len(sequence)}\n")
                        f.write(f"Number of factors: {len(factors)}\n")
                        f.write("Factors (position, length, reference):\n")
                        for i, (pos, length, ref) in enumerate(factors):
                            f.write(f"{i+1:4d}: ({pos:6d}, {length:4d}, {ref:6d})\n")
                    print(f"  Saved factors to {factors_text_file}")
                except Exception as e:
                    print(f"  Warning: Failed to save text factors for {seq_id}: {e}")

            # Save factors as binary
            factors_binary_file = None
            if save_factors_binary:
                factors_binary_file = output_dir / f"factors_{seq_id}.bin"
                try:
                    # Create a temporary file with just this sequence
                    temp_fasta = output_dir / f"temp_{seq_id}.fasta"
                    with open(temp_fasta, 'w') as f:
                        f.write(f">{seq_id}\n{sequence}\n")

                    write_factors_binary_file(str(temp_fasta), str(factors_binary_file))
                    temp_fasta.unlink()  # Clean up temp file
                    print(f"  Saved binary factors to {factors_binary_file}")
                except Exception as e:
                    print(f"  Warning: Failed to save binary factors for {seq_id}: {e}")

            # Create plot
            plot_file = output_dir / f"plot_{seq_id}.png"
            try:
                from ..utils import plot_factor_lengths
                plot_factor_lengths(factors, save_path=plot_file, show_plot=False)
                print(f"  Saved plot to {plot_file}")
            except Exception as e:
                print(f"  Warning: Failed to create plot for {seq_id}: {e}")
                plot_file = None

            # Store results
            results[seq_id] = {
                'sequence_length': len(sequence),
                'num_factors': len(factors),
                'factors_file': str(factors_text_file) if factors_text_file else None,
                'binary_file': str(factors_binary_file) if factors_binary_file else None,
                'plot_file': str(plot_file) if plot_file else None,
                'factors': factors
            }

            processed_count += 1

        print(f"\nProcessed {len(results)} sequences from FASTA successfully")

    else:
        # Process binary factors file
        print(f"Reading factors from binary file {input_filepath}...")
        
        try:
            # Try to read with metadata first (for multi-sequence files)
            metadata = read_factors_binary_file_with_metadata(input_filepath)
            factors = metadata['factors']
            sequence_names = metadata.get('sequence_names', ['sequence'])
            sequence_lengths = metadata.get('sequence_lengths', [])
            sentinel_factor_indices = metadata.get('sentinel_factor_indices', [])
            
            print(f"Loaded {len(factors)} factors with metadata for {len(sequence_names)} sequences")
            
            # For binary files with multiple sequences, we need to split factors by sequence
            if len(sequence_names) > 1 and sentinel_factor_indices:
                # Split factors by sequence using sentinel indices
                factor_groups = []
                start_idx = 0
                
                for sentinel_idx in sentinel_factor_indices:
                    factor_groups.append(factors[start_idx:sentinel_idx])
                    start_idx = sentinel_idx + 1  # Skip the sentinel factor
                
                # Add the last group (after the last sentinel)
                if start_idx < len(factors):
                    factor_groups.append(factors[start_idx:])
                
                # Process each sequence
                for i, (seq_id, seq_factors) in enumerate(zip(sequence_names, factor_groups)):
                    if max_sequences is not None and i >= max_sequences:
                        break
                        
                    print(f"Processing sequence {seq_id} ({len(seq_factors)} factors)...")
                    
                    # Create plot
                    plot_file = output_dir / f"plot_{seq_id}.png"
                    try:
                        from ..utils import plot_factor_lengths
                        plot_factor_lengths(seq_factors, save_path=plot_file, show_plot=False)
                        print(f"  Saved plot to {plot_file}")
                    except Exception as e:
                        print(f"  Warning: Failed to create plot for {seq_id}: {e}")
                        plot_file = None
                    
                    # Store results
                    seq_length = sequence_lengths[i] if i < len(sequence_lengths) else None
                    results[seq_id] = {
                        'sequence_length': seq_length,
                        'num_factors': len(seq_factors),
                        'factors_file': None,  # No text file created from binary input
                        'binary_file': str(input_filepath),  # Original binary file
                        'plot_file': str(plot_file) if plot_file else None,
                        'factors': seq_factors
                    }
            else:
                # Single sequence binary file
                seq_id = sequence_names[0] if sequence_names else input_filepath.stem
                print(f"Processing single sequence {seq_id} ({len(factors)} factors)...")
                
                # Create plot
                plot_file = output_dir / f"plot_{seq_id}.png"
                try:
                    from ..utils import plot_factor_lengths
                    plot_factor_lengths(factors, save_path=plot_file, show_plot=False)
                    print(f"  Saved plot to {plot_file}")
                except Exception as e:
                    print(f"  Warning: Failed to create plot for {seq_id}: {e}")
                    plot_file = None
                
                # Store results
                seq_length = sequence_lengths[0] if sequence_lengths else None
                results[seq_id] = {
                    'sequence_length': seq_length,
                    'num_factors': len(factors),
                    'factors_file': None,  # No text file created from binary input
                    'binary_file': str(input_filepath),  # Original binary file
                    'plot_file': str(plot_file) if plot_file else None,
                    'factors': factors
                }
                
        except Exception as e:
            # Fallback: try to read as simple binary file without metadata
            try:
                from ..utils import read_factors_binary_file
                factors = read_factors_binary_file(input_filepath)
                seq_id = input_filepath.stem
                
                print(f"Loaded {len(factors)} factors from simple binary file")
                print(f"Processing sequence {seq_id} ({len(factors)} factors)...")
                
                # Create plot
                plot_file = output_dir / f"plot_{seq_id}.png"
                try:
                    from ..utils import plot_factor_lengths
                    plot_factor_lengths(factors, save_path=plot_file, show_plot=False)
                    print(f"  Saved plot to {plot_file}")
                except Exception as e:
                    print(f"  Warning: Failed to create plot for {seq_id}: {e}")
                    plot_file = None
                
                # Store results
                results[seq_id] = {
                    'sequence_length': None,  # Unknown from simple binary file
                    'num_factors': len(factors),
                    'factors_file': None,  # No text file created from binary input
                    'binary_file': str(input_filepath),  # Original binary file
                    'plot_file': str(plot_file) if plot_file else None,
                    'factors': factors
                }
                
            except Exception as e2:
                raise PlotError(f"Failed to read binary factors file: {e2}")

        print(f"\nProcessed {len(results)} sequences from binary file successfully")

    return results


def plot_multiple_seq_self_lz_factor_plot_from_file(
    fasta_filepath: Optional[Union[str, Path]] = None,
    factors_filepath: Optional[Union[str, Path]] = None,
    name: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    return_panel: bool = False
) -> Optional["panel.viewable.Viewable"]:
    """
    Create an interactive Datashader/Panel factor plot for multiple DNA sequences from a FASTA file or binary factors file.

    This function reads factors either from a FASTA file (by factorizing multiple DNA sequences)
    or from an enhanced binary factors file with metadata. It creates a high-performance
    interactive plot using Datashader and Panel with level-of-detail rendering, zoom/pan-aware 
    decimation, hover functionality, and sequence boundaries visualization.

    Args:
        fasta_filepath: Path to the FASTA file containing DNA sequences (mutually exclusive with factors_filepath)
        factors_filepath: Path to binary factors file with metadata (mutually exclusive with fasta_filepath)
        name: Optional name for the plot title (defaults to input filename)
        save_path: Optional path to save the plot image (PNG export)
        show_plot: Whether to display/serve the plot
        return_panel: Whether to return the Panel app for embedding

    Returns:
        Panel app if return_panel=True, otherwise None

    Raises:
        PlotError: If plotting fails or input files cannot be processed
        FileNotFoundError: If input file doesn't exist
        ImportError: If required dependencies are missing
        ValueError: If both or neither input files are provided
    """
    # Check for required dependencies
    try:
        import numpy as np
        import pandas as pd
        import holoviews as hv
        import datashader as ds
        import panel as pn
        import colorcet as cc
        from holoviews.operation.datashader import datashade, dynspread
        from holoviews import streams
        import bokeh
    except ImportError as e:
        missing_dep = str(e).split("'")[1] if "'" in str(e) else str(e)
        raise ImportError(
            f"Missing required dependency: {missing_dep}. "
            f"Install with: pip install 'noLZSS[panel]' or "
            f"pip install numpy pandas holoviews bokeh panel datashader colorcet"
        )

    # Initialize extensions
    hv.extension('bokeh')
    pn.extension()

    from .._noLZSS import factorize_fasta_multiple_dna_w_rc
    from ..utils import read_factors_binary_file_with_metadata

    # Validate input arguments
    if (fasta_filepath is None) == (factors_filepath is None):
        raise ValueError("Exactly one of fasta_filepath or factors_filepath must be provided")

    # Determine input type and file path
    if fasta_filepath is not None:
        input_filepath = Path(fasta_filepath)
        input_type = "fasta"
    else:
        if factors_filepath is None:
            raise ValueError("Either fasta_filepath or factors_filepath must be provided")
        input_filepath = Path(factors_filepath)
        input_type = "binary"

    if not input_filepath.exists():
        raise FileNotFoundError(f"Input file not found: {input_filepath}")

    # Determine plot title
    if name is None:
        name = input_filepath.stem

    try:
        # Get factors and metadata based on input type
        if input_type == "fasta":
            print(f"Reading and factorizing sequences from {input_filepath}...")
            factors, sentinel_factor_indices, sequence_names = factorize_fasta_multiple_dna_w_rc(str(input_filepath))
        else:
            print(f"Reading factors from binary file {input_filepath}...")
            metadata = read_factors_binary_file_with_metadata(input_filepath)
            factors = metadata['factors']
            sentinel_factor_indices = metadata['sentinel_factor_indices']
            sequence_names = metadata['sequence_names']

        print(f"Loaded {len(factors)} factors with {len(sentinel_factor_indices)} sentinels")
        print(f"Sequence names: {sequence_names}")
        
        if not factors:
            raise PlotError("No factors found in input file")

        # Build DataFrame with plot coordinates
        print("Building factor DataFrame...")
        x0_vals = []
        y0_vals = []
        x1_vals = []
        y1_vals = []
        lengths = []
        dirs = []
        starts = []
        refs = []
        ends = []
        is_rcs = []

        for factor in factors:
            start, length, ref, is_rc = factor
            
            # Calculate coordinates
            x0 = start
            x1 = start + length
            
            if is_rc:
                # Reverse complement: y0 = ref + length, y1 = ref
                y0 = ref + length
                y1 = ref
                dir_val = 1
            else:
                # Forward: y0 = ref, y1 = ref + length
                y0 = ref
                y1 = ref + length
                dir_val = 0
            
            x0_vals.append(x0)
            y0_vals.append(y0)
            x1_vals.append(x1)
            y1_vals.append(y1)
            lengths.append(length)
            dirs.append(dir_val)
            starts.append(start)
            refs.append(ref)
            ends.append(x1)
            is_rcs.append(is_rc)

        # Create DataFrame
        df = pd.DataFrame({
            'x0': x0_vals,
            'y0': y0_vals,
            'x1': x1_vals,
            'y1': y1_vals,
            'length': lengths,
            'dir': dirs,
            'start': starts,
            'ref': refs,
            'end': ends,
            'is_rc': is_rcs
        })

        print(f"DataFrame created with {len(df)} factors")

        # Calculate sentinel positions for lines and labels
        sentinel_positions = []
        sequence_boundaries = []  # (start_pos, end_pos, sequence_name)
        
        if sentinel_factor_indices:
            # Get positions of sentinel factors
            for idx in sentinel_factor_indices:
                if idx < len(factors):
                    sentinel_start = factors[idx][0]  # start position of sentinel factor
                    sentinel_positions.append(sentinel_start)
            
            # Calculate sequence boundaries
            prev_pos = 0
            for i, pos in enumerate(sentinel_positions):
                seq_name = sequence_names[i] if i < len(sequence_names) else f"seq_{i}"
                sequence_boundaries.append((prev_pos, pos, seq_name))
                prev_pos = pos + 1  # Skip the sentinel itself
            
            # Add the last sequence
            if len(sequence_names) > len(sentinel_positions):
                last_name = sequence_names[len(sentinel_positions)]
            else:
                last_name = f"seq_{len(sentinel_positions)}"
            
            # Find the maximum position for the last sequence
            max_pos = max(max(df['x1']), max(df['y1'])) if len(df) > 0 else prev_pos
            sequence_boundaries.append((prev_pos, max_pos, last_name))
        else:
            # No sentinels - single sequence
            seq_name = sequence_names[0] if sequence_names else "sequence"
            max_pos = max(max(df['x1']), max(df['y1'])) if len(df) > 0 else 1000
            sequence_boundaries.append((0, max_pos, seq_name))

        print(f"Sequence boundaries: {sequence_boundaries}")
        print(f"Sentinel positions: {sentinel_positions}")

        # Define color mapping
        def create_base_layers(df_filtered):
            """Create the base datashaded layers"""
            # Split data by direction
            df_fwd = df_filtered[df_filtered['dir'] == 0]
            df_rc = df_filtered[df_filtered['dir'] == 1]
            
            # Create HoloViews segments
            segments_fwd = hv.Segments(
                df_fwd, 
                kdims=['x0','y0','x1','y1'], 
                vdims=['length','start','ref','end']
            ).opts(color='blue')
            
            segments_rc = hv.Segments(
                df_rc, 
                kdims=['x0','y0','x1','y1'], 
                vdims=['length','start','ref','end']
            ).opts(color='red')
            
            # Apply datashader with max aggregator
            shaded_fwd = dynspread(
                datashade(
                    segments_fwd, 
                    aggregator=ds.max('length'),
                    cmap=['white', 'blue']
                )
            )
            
            shaded_rc = dynspread(
                datashade(
                    segments_rc, 
                    aggregator=ds.max('length'),
                    cmap=['white', 'red']
                )
            )
            
            return shaded_fwd * shaded_rc

        # Create range streams for interactivity
        rangexy = streams.RangeXY()
        
        def create_hover_overlay(x_range, y_range, df_filtered, k_per_bin=1, plot_width=800):
            """Create decimated overlay for hover functionality"""
            if x_range is None or y_range is None:
                return hv.Segments([])
            
            x_min, x_max = x_range
            y_min, y_max = y_range
            
            # Filter to visible range with some padding
            x_pad = (x_max - x_min) * 0.1
            y_pad = (y_max - y_min) * 0.1
            
            visible_mask = (
                (df_filtered['x0'] <= x_max + x_pad) & 
                (df_filtered['x1'] >= x_min - x_pad) &
                (df_filtered['y0'] <= y_max + y_pad) & 
                (df_filtered['y1'] >= y_min - y_pad)
            )
            
            visible_df = df_filtered[visible_mask].copy()
            
            if len(visible_df) == 0:
                return hv.Segments([])
            
            # Screen-space decimation
            nbins = min(plot_width, 2000)
            
            # Calculate midpoints for binning
            visible_df['mid_x'] = (visible_df['x0'] + visible_df['x1']) / 2
            
            # Bin by x-coordinate
            bins = np.linspace(x_min - x_pad, x_max + x_pad, nbins + 1)
            visible_df['bin'] = pd.cut(visible_df['mid_x'], bins, labels=False, include_lowest=True)
            
            # Keep top-k by length per bin
            top_k_df = (visible_df.groupby('bin', group_keys=False)
                        .apply(lambda x: x.nlargest(k_per_bin, 'length'))
                        .reset_index(drop=True))
            
            if len(top_k_df) == 0:
                return hv.Segments([])
            
            # Create hover data with direction labels
            top_k_df['direction'] = top_k_df['is_rc'].map({True: 'reverse-complement', False: 'forward'})
            
            # Create segments with hover info
            segments = hv.Segments(
                top_k_df,
                kdims=['x0','y0','x1','y1'],
                vdims=['start', 'length', 'end', 'ref', 'direction', 'is_rc']
            ).opts(
                tools=['hover'],
                line_width=2,
                alpha=0.9,
                color='is_rc',
                cmap={True: 'red', False: 'blue'},
                hover_tooltips=[
                    ('Start', '@start'),
                    ('Length', '@length'), 
                    ('End', '@end'),
                    ('Reference', '@ref'),
                    ('Direction', '@direction'),
                    ('Is Reverse Complement', '@is_rc')
                ]
            )
            
            return segments

        # Create widgets
        length_range_slider = pn.widgets.IntRangeSlider(
            name="Length Filter",
            start=int(df['length'].min()),
            end=int(df['length'].max()),
            value=(int(df['length'].min()), int(df['length'].max())),
            step=1
        )
        
        show_overlay_checkbox = pn.widgets.Checkbox(
            name="Show hover overlay",
            value=True
        )
        
        k_spinner = pn.widgets.IntInput(
            name="Top-k per pixel bin",
            value=1,
            start=1,
            end=5
        )
        
        colormap_select = pn.widgets.Select(
            name="Colormap",
            value='gray',
            options=['gray', 'viridis', 'plasma', 'inferno']
        )

        # Create dynamic plot function
        def create_plot(length_range, show_overlay, k_per_bin, colormap_name):
            length_min, length_max = length_range
            # Filter by length
            df_filtered = df[
                (df['length'] >= length_min) & 
                (df['length'] <= length_max)
            ].copy()
            
            if len(df_filtered) == 0:
                return hv.Text(0, 0, "No data in range").opts(width=800, height=800)
            
            # Create base layers
            base_plot = create_base_layers(df_filtered)
            
            # Add diagonal y=x line
            max_val = max(df_filtered[['x1', 'y1']].max())
            min_val = min(df_filtered[['x0', 'y0']].min())
            diagonal = hv.Curve([(min_val, min_val), (max_val, max_val)]).opts(
                line_dash='dashed',
                line_color='gray',
                line_width=1,
                alpha=0.5
            )
            
            # Add sentinel lines and sequence labels
            sentinel_elements = []
            
            for pos in sentinel_positions:
                if min_val <= pos <= max_val:
                    # Vertical line at sentinel position
                    v_line = hv.VLine(pos).opts(
                        line_color='red',
                        line_width=2,
                        alpha=0.7,
                        line_dash='solid'
                    )
                    sentinel_elements.append(v_line)
                    
                    # Horizontal line at sentinel position  
                    h_line = hv.HLine(pos).opts(
                        line_color='red',
                        line_width=2,
                        alpha=0.7,
                        line_dash='solid'
                    )
                    sentinel_elements.append(h_line)
            
            # Add sequence name labels
            label_elements = []
            for start_pos, end_pos, seq_name in sequence_boundaries:
                mid_pos = (start_pos + end_pos) / 2
                if min_val <= mid_pos <= max_val:
                    # X-axis label (bottom)
                    x_label = hv.Text(mid_pos, min_val - (max_val - min_val) * 0.05, seq_name).opts(
                        text_color='blue',
                        text_font_size='10pt',
                        text_align='center'
                    )
                    label_elements.append(x_label)
                    
                    # Y-axis label (left side)  
                    y_label = hv.Text(min_val - (max_val - min_val) * 0.05, mid_pos, seq_name).opts(
                        text_color='blue', 
                        text_font_size='10pt',
                        text_align='center',
                        angle=90
                    )
                    label_elements.append(y_label)
            
            # Combine all plot elements
            plot = base_plot * diagonal
            
            # Add sentinel lines
            for element in sentinel_elements:
                plot = plot * element
                
            # Add sequence labels
            for element in label_elements:
                plot = plot * element
            
            # Add hover overlay if requested
            if show_overlay:
                # Use rangexy stream to get current view
                overlay_func = lambda x_range, y_range: create_hover_overlay(
                    x_range, y_range, df_filtered, k_per_bin
                )
                hover_dmap = hv.DynamicMap(overlay_func, streams=[rangexy])
                plot = plot * hover_dmap
            
            # Configure plot options
            plot = plot.opts(
                width=800,
                height=800,
                aspect='equal',
                xlabel=f'Position in concatenated sequence ({name}) - Sequences: {", ".join([b[2] for b in sequence_boundaries])}',
                ylabel=f'Reference position ({name}) - Sequences: {", ".join([b[2] for b in sequence_boundaries])}',
                title=f'LZ Factor Plot - {name} ({len(sequence_boundaries)} sequences)',
                toolbar='above'
            )
            
            return plot

        # Bind widgets to plot function
        interactive_plot = pn.bind(
            create_plot,
            length_range=length_range_slider.param.value,
            show_overlay=show_overlay_checkbox,
            k_per_bin=k_spinner,
            colormap_name=colormap_select
        )

        # Export functionality
        def export_png():
            # This is a placeholder - actual implementation would use bokeh.io.export_png
            print("PNG export not implemented - requires selenium/chromedriver")
            return

        export_button = pn.widgets.Button(name="Export PNG", button_type="primary")
        export_button.on_click(lambda event: export_png())

        # Create Panel app layout
        controls = pn.Column(
            "## Controls",
            length_range_slider,
            show_overlay_checkbox,
            k_spinner,
            colormap_select,
            export_button,
            width=300
        )

        app = pn.Row(
            controls,
            pn.panel(interactive_plot, width=800, height=800)
        )

        # Handle save_path
        if save_path:
            print(f"Note: PNG export to {save_path} requires additional setup (selenium/chromedriver)")

        # Handle display/serving
        if show_plot:
            if return_panel:
                return app
            else:
                # In jupyter notebooks, the app will display automatically
                # For script execution, we need to serve
                try:
                    # Check if we're in a notebook
                    get_ipython()  # noqa: F821
                    return app  # In notebook, just return for display
                except NameError:
                    # Not in notebook, serve the app
                    if __name__ == "__main__":
                        pn.serve(app, show=True, port=5007)
                    else:
                        print("To serve the app, run: panel serve script.py --show")
                        return app
        elif return_panel:
            return app
        else:
            return None

    except Exception as e:
        raise PlotError(f"Failed to create interactive LZ factor plot: {e}")


def plot_reference_seq_lz_factor_plot_simple(
    reference_seq: Union[str, bytes],
    target_seq: Union[str, bytes],
    factors: Optional[List[Tuple[int, int, int, bool]]] = None,
    factors_filepath: Optional[Union[str, Path]] = None,
    reference_name: str = "Reference",
    target_name: str = "Target",
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True
) -> None:
    """
    Create a simple matplotlib factor plot for a DNA sequence factorized with a reference sequence.
    
    This function creates a plot compatible with the output of factorize_w_reference_seq().
    The plot shows the reference sequence at the beginning, concatenated with the target sequence,
    and uses distinct colors for reference vs target regions.
    
    Args:
        reference_seq: Reference DNA sequence (A, C, T, G - case insensitive)
        target_seq: Target DNA sequence (A, C, T, G - case insensitive)
        factors: Optional list of (start, length, ref, is_rc) tuples from factorize_w_reference_seq().
                If None, will compute factors automatically.
        factors_filepath: Optional path to binary factors file (mutually exclusive with factors)
        reference_name: Name for the reference sequence (default: "Reference")
        target_name: Name for the target sequence (default: "Target")
        save_path: Optional path to save the plot image
        show_plot: Whether to display the plot
        
    Raises:
        PlotError: If plotting fails or input sequences are invalid
        ValueError: If both factors and factors_filepath are provided
        ImportError: If matplotlib is not available
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        import numpy as np
    except ImportError as e:
        missing_dep = str(e).split("'")[1] if "'" in str(e) else str(e)
        raise ImportError(
            f"Missing required dependency: {missing_dep}. "
            f"Install with: pip install matplotlib"
        )

    # Validate input arguments
    if factors is not None and factors_filepath is not None:
        raise ValueError("Cannot provide both factors and factors_filepath")
    
    # Convert sequences to strings if they're bytes
    if isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    # Get factors if not provided
    if factors is None:
        if factors_filepath is not None:
            from ..utils import read_factors_binary_file_with_metadata
            metadata = read_factors_binary_file_with_metadata(factors_filepath)
            factors = metadata['factors']
        else:
            from .sequences import factorize_w_reference_seq
            factors = factorize_w_reference_seq(reference_seq, target_seq)
    
    if not factors:
        raise PlotError("No factors available for plotting")
    
    # Calculate sequence boundaries
    ref_length = len(reference_seq)
    target_start = ref_length + 1  # +1 for sentinel
    target_length = len(target_seq)
    total_length = target_start + target_length
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot factors with different colors for forward/reverse and reference/target
    for start, length, ref_pos, is_rc in factors:
        end = start + length
        ref_end = ref_pos + length
        
        # Determine color based on region and orientation
        if start >= target_start:
            # Target region
            color = 'red' if not is_rc else 'darkorange'
            alpha = 0.7
        else:
            # Reference region (shouldn't happen with factorize_w_reference_seq)
            color = 'blue' if not is_rc else 'darkblue'
            alpha = 0.7
        
        # Draw line segment from (start, ref_pos) to (end, ref_end)
        ax.plot([start, end], [ref_pos, ref_end], color=color, alpha=alpha, linewidth=2)
    
    # Add diagonal reference line
    max_pos = max(total_length, max(ref_pos + length for _, length, ref_pos, _ in factors))
    ax.plot([0, max_pos], [0, max_pos], 'gray', linestyle='--', alpha=0.5, linewidth=1)
    
    # Add sequence boundary lines
    ax.axvline(x=target_start - 0.5, color='green', linewidth=3, alpha=0.8, 
               label=f'Boundary: {reference_name}|{target_name}')
    ax.axhline(y=target_start - 0.5, color='green', linewidth=3, alpha=0.8)
    
    # Add sequence region backgrounds
    ax.axvspan(0, ref_length, alpha=0.1, color='blue', label=f'{reference_name} region')
    ax.axvspan(target_start, total_length, alpha=0.1, color='red', label=f'{target_name} region')
    
    # Set labels and title
    ax.set_xlabel(f'Position in concatenated sequence ({reference_name} + {target_name})')
    ax.set_ylabel(f'Reference position')
    ax.set_title(f'Reference Sequence LZ Factor Plot\n{reference_name} vs {target_name}')
    
    # Add legend
    legend_elements = [
        patches.Patch(color='red', alpha=0.7, label=f'{target_name} forward factors'),
        patches.Patch(color='darkorange', alpha=0.7, label=f'{target_name} reverse complement factors'),
        patches.Patch(color='blue', alpha=0.1, label=f'{reference_name} region'),
        patches.Patch(color='red', alpha=0.1, label=f'{target_name} region'),
        patches.Patch(color='green', alpha=0.8, label='Sequence boundary')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    # Set equal aspect ratio and grid
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3)
    
    # Add text annotations for sequence regions
    ax.text(ref_length/2, -max_pos*0.05, reference_name, ha='center', va='top', 
            fontsize=12, color='blue', weight='bold')
    ax.text((target_start + total_length)/2, -max_pos*0.05, target_name, ha='center', va='top',
            fontsize=12, color='red', weight='bold')
    ax.text(-max_pos*0.05, ref_length/2, reference_name, ha='right', va='center',
            fontsize=12, color='blue', weight='bold', rotation=90)
    ax.text(-max_pos*0.05, (target_start + total_length)/2, target_name, ha='right', va='center',
            fontsize=12, color='red', weight='bold', rotation=90)
    
    plt.tight_layout()
    
    # Save plot if requested
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {save_path}")
    
    # Show plot
    if show_plot:
        plt.show()
    else:
        plt.close(fig)


def plot_reference_seq_lz_factor_plot(
    reference_seq: Union[str, bytes],
    target_seq: Union[str, bytes],
    factors: Optional[List[Tuple[int, int, int, bool]]] = None,
    factors_filepath: Optional[Union[str, Path]] = None,
    reference_name: str = "Reference",
    target_name: str = "Target",
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    return_panel: bool = False
) -> Optional["panel.viewable.Viewable"]:
    """
    Create an interactive Datashader/Panel factor plot for a DNA sequence factorized with a reference sequence.
    
    This function creates a plot compatible with the output of factorize_w_reference_seq().
    The plot shows the reference sequence at the beginning, concatenated with the target sequence,
    and uses distinct colors for reference vs target regions.
    
    Args:
        reference_seq: Reference DNA sequence (A, C, T, G - case insensitive)
        target_seq: Target DNA sequence (A, C, T, G - case insensitive)
        factors: Optional list of (start, length, ref, is_rc) tuples from factorize_w_reference_seq().
                If None, will compute factors automatically.
        factors_filepath: Optional path to binary factors file (mutually exclusive with factors)
        reference_name: Name for the reference sequence (default: "Reference")
        target_name: Name for the target sequence (default: "Target")
        save_path: Optional path to save the plot image (PNG export)
        show_plot: Whether to display/serve the plot
        return_panel: Whether to return the Panel app for embedding
        
    Returns:
        Panel app if return_panel=True, otherwise None
        
    Raises:
        PlotError: If plotting fails or input sequences are invalid
        ValueError: If both factors and factors_filepath are provided
        ImportError: If required dependencies are missing
    """
    # Check for required dependencies
    try:
        import numpy as np
        import pandas as pd
        import holoviews as hv
        import datashader as ds
        import panel as pn
        import colorcet as cc
        from holoviews.operation.datashader import datashade, dynspread
        from holoviews import streams
        import bokeh
    except ImportError as e:
        missing_dep = str(e).split("'")[1] if "'" in str(e) else str(e)
        raise ImportError(
            f"Missing required dependency: {missing_dep}. "
            f"Install with: pip install 'noLZSS[panel]' or "
            f"pip install numpy pandas holoviews bokeh panel datashader colorcet"
        )

    # Initialize extensions
    hv.extension('bokeh')
    pn.extension()

    # Validate input arguments
    if factors is not None and factors_filepath is not None:
        raise ValueError("Cannot provide both factors and factors_filepath")
    
    # Convert sequences to strings if they're bytes
    if isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    # Get factors if not provided
    if factors is None:
        if factors_filepath is not None:
            from ..utils import read_factors_binary_file_with_metadata
            metadata = read_factors_binary_file_with_metadata(factors_filepath)
            factors = metadata['factors']
        else:
            from .sequences import factorize_w_reference_seq
            factors = factorize_w_reference_seq(reference_seq, target_seq)
    
    if not factors:
        raise PlotError("No factors available for plotting")
    
    # Calculate sequence boundaries
    ref_length = len(reference_seq)
    target_start = ref_length + 1  # +1 for sentinel
    target_length = len(target_seq)
    total_length = target_start + target_length
    
    # Create sequence boundaries for visualization
    sequence_boundaries = [
        (0, ref_length, reference_name),
        (target_start, total_length, target_name)
    ]
    
    # Convert factors to DataFrame for plotting
    factor_data = []
    for start, length, ref_pos, is_rc in factors:
        end = start + length
        factor_data.append({
            'x0': start,
            'x1': end,
            'y0': ref_pos,
            'y1': ref_pos + length,
            'length': length,
            'is_rc': is_rc,
            'region': 'target' if start >= target_start else 'reference'
        })
    
    df = pd.DataFrame(factor_data)
    
    if len(df) == 0:
        raise PlotError("No valid factors to plot")
    
    # Create plot name
    plot_name = f"{reference_name} vs {target_name}"
    
    # Create interactive plot components
    def create_base_layers(df_filtered):
        """Create base datashader layers for factors"""
        if len(df_filtered) == 0:
            return hv.Text(0, 0, "No data").opts(width=800, height=800)
        
        # Separate reference and target factors for different colors
        df_ref = df_filtered[df_filtered['region'] == 'reference']
        df_target = df_filtered[df_filtered['region'] == 'target']
        
        layers = []
        
        # Reference factors (if any) - use blue tones
        if len(df_ref) > 0:
            ref_forward = df_ref[~df_ref['is_rc']]
            ref_reverse = df_ref[df_ref['is_rc']]
            
            if len(ref_forward) > 0:
                ref_forward_segments = hv.Segments(
                    ref_forward, ['x0', 'y0', 'x1', 'y1']
                ).opts(color='blue', alpha=0.7, line_width=2)
                layers.append(ref_forward_segments)
            
            if len(ref_reverse) > 0:
                ref_reverse_segments = hv.Segments(
                    ref_reverse, ['x0', 'y0', 'x1', 'y1']
                ).opts(color='darkblue', alpha=0.7, line_width=2)
                layers.append(ref_reverse_segments)
        
        # Target factors - use red/orange tones  
        if len(df_target) > 0:
            target_forward = df_target[~df_target['is_rc']]
            target_reverse = df_target[df_target['is_rc']]
            
            if len(target_forward) > 0:
                target_forward_segments = hv.Segments(
                    target_forward, ['x0', 'y0', 'x1', 'y1']
                ).opts(color='red', alpha=0.7, line_width=2)
                layers.append(target_forward_segments)
            
            if len(target_reverse) > 0:
                target_reverse_segments = hv.Segments(
                    target_reverse, ['x0', 'y0', 'x1', 'y1']
                ).opts(color='darkorange', alpha=0.7, line_width=2)
                layers.append(target_reverse_segments)
        
        if not layers:
            return hv.Text(0, 0, "No data").opts(width=800, height=800)
        
        # Combine layers
        plot = layers[0]
        for layer in layers[1:]:
            plot = plot * layer
            
        return plot
    
    def create_hover_overlay(x_range, y_range, df_filtered, k_per_bin):
        """Create hover overlay for detailed factor information"""
        if x_range is None or y_range is None or len(df_filtered) == 0:
            return hv.Text(0, 0, "").opts(alpha=0)
        
        x_min, x_max = x_range
        y_min, y_max = y_range
        
        # Filter data to current view
        view_data = df_filtered[
            (df_filtered['x0'] >= x_min) & (df_filtered['x1'] <= x_max) &
            (df_filtered['y0'] >= y_min) & (df_filtered['y1'] <= y_max)
        ]
        
        if len(view_data) == 0:
            return hv.Text(0, 0, "").opts(alpha=0)
        
        # Sample data if too many points
        if len(view_data) > k_per_bin:
            view_data = view_data.sample(n=k_per_bin)
        
        # Create points for hover
        hover_data = []
        for _, row in view_data.iterrows():
            hover_data.append({
                'x': (row['x0'] + row['x1']) / 2,
                'y': (row['y0'] + row['y1']) / 2,
                'start': int(row['x0']),
                'length': int(row['length']),
                'ref_pos': int(row['y0']),
                'is_rc': row['is_rc'],
                'region': row['region']
            })
        
        if not hover_data:
            return hv.Text(0, 0, "").opts(alpha=0)
        
        hover_df = pd.DataFrame(hover_data)
        
        return hv.Points(
            hover_df, ['x', 'y'], ['start', 'length', 'ref_pos', 'is_rc', 'region']
        ).opts(
            size=8, alpha=0.8, color='yellow',
            tools=['hover'], hover_alpha=1.0
        )
    
    # Set up interactive plot
    rangexy = streams.RangeXY()
    
    # Create dynamic plot function
    def create_plot(length_range, show_overlay, k_per_bin, colormap_name):
        length_min, length_max = length_range
        # Filter by length
        df_filtered = df[
            (df['length'] >= length_min) & 
            (df['length'] <= length_max)
        ].copy()
        
        if len(df_filtered) == 0:
            return hv.Text(0, 0, "No data in range").opts(width=800, height=800)
        
        # Create base layers
        base_plot = create_base_layers(df_filtered)
        
        # Add diagonal y=x line
        max_val = max(df_filtered[['x1', 'y1']].max())
        min_val = min(df_filtered[['x0', 'y0']].min())
        diagonal = hv.Curve([(min_val, min_val), (max_val, max_val)]).opts(
            line_dash='dashed',
            line_color='gray',
            line_width=1,
            alpha=0.5
        )
        
        # Add sequence boundary lines
        boundary_elements = []
        
        # Vertical line separating reference and target
        sep_pos = target_start - 0.5  # Position between ref and target
        if min_val <= sep_pos <= max_val:
            v_line = hv.VLine(sep_pos).opts(
                line_color='green',
                line_width=3,
                alpha=0.8,
                line_dash='solid'
            )
            boundary_elements.append(v_line)
            
            # Horizontal line at same position
            h_line = hv.HLine(sep_pos).opts(
                line_color='green',
                line_width=3,
                alpha=0.8,
                line_dash='solid'
            )
            boundary_elements.append(h_line)
        
        # Add sequence name labels
        label_elements = []
        for start_pos, end_pos, seq_name in sequence_boundaries:
            mid_pos = (start_pos + end_pos) / 2
            if min_val <= mid_pos <= max_val:
                # Determine color based on sequence type
                label_color = 'blue' if seq_name == reference_name else 'red'
                
                # X-axis label (bottom)
                x_label = hv.Text(mid_pos, min_val - (max_val - min_val) * 0.05, seq_name).opts(
                    text_color=label_color,
                    text_font_size='12pt',
                    text_align='center'
                )
                label_elements.append(x_label)
                
                # Y-axis label (left side)  
                y_label = hv.Text(min_val - (max_val - min_val) * 0.05, mid_pos, seq_name).opts(
                    text_color=label_color, 
                    text_font_size='12pt',
                    text_align='center',
                    angle=90
                )
                label_elements.append(y_label)
        
        # Combine all plot elements
        plot = base_plot * diagonal
        
        # Add boundary lines
        for element in boundary_elements:
            plot = plot * element
            
        # Add sequence labels
        for element in label_elements:
            plot = plot * element
        
        # Add hover overlay if requested
        if show_overlay:
            # Use rangexy stream to get current view
            overlay_func = lambda x_range, y_range: create_hover_overlay(
                x_range, y_range, df_filtered, k_per_bin
            )
            hover_dmap = hv.DynamicMap(overlay_func, streams=[rangexy])
            plot = plot * hover_dmap
        
        # Configure plot options
        plot = plot.opts(
            width=800,
            height=800,
            aspect='equal',
            xlabel=f'Position in concatenated sequence ({plot_name})',
            ylabel=f'Reference position ({plot_name})',
            title=f'Reference Sequence LZ Factor Plot - {plot_name}',
            toolbar='above'
        )
        
        return plot
    
    # Create interactive controls
    length_min, length_max = df['length'].min(), df['length'].max()
    
    length_slider = pn.widgets.RangeSlider(
        name='Factor Length Range',
        start=length_min,
        end=length_max,
        value=(length_min, length_max),
        step=1
    )
    
    overlay_toggle = pn.widgets.Toggle(
        name='Show Hover Details',
        value=True
    )
    
    k_spinner = pn.widgets.IntInput(
        name='Max Points for Hover',
        value=min(1000, len(df)),
        start=100,
        end=5000,
        step=100
    )
    
    colormap_select = pn.widgets.Select(
        name='Colormap',
        value='viridis',
        options=['viridis', 'plasma', 'inferno', 'magma', 'cividis']
    )
    
    # Create dynamic plot with controls
    dmap = hv.DynamicMap(
        create_plot,
        kdims=[
            hv.Dimension('length_range', values=[length_slider.param.value]),
            hv.Dimension('show_overlay', values=[overlay_toggle.param.value]),
            hv.Dimension('k_per_bin', values=[k_spinner.param.value]),
            hv.Dimension('colormap_name', values=[colormap_select.param.value])
        ]
    )
    
    # Create layout
    controls = pn.Column(
        pn.pane.Markdown("### Plot Controls"),
        length_slider,
        overlay_toggle,
        k_spinner,
        colormap_select,
        pn.pane.Markdown(f"**Dataset Info:**\n- Total factors: {len(df)}\n- Reference length: {ref_length}\n- Target length: {target_length}"),
        width=300
    )
    
    plot_pane = pn.pane.HoloViews(dmap, width=850, height=850)
    
    app = pn.Row(controls, plot_pane)
    
    # Save plot if requested
    if save_path:
        try:
            # Create a static version for saving
            static_plot = create_plot(
                (length_min, length_max),
                False,  # No hover for static
                1000,
                'viridis'
            )
            
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            
            hv.save(static_plot, save_path, fmt='png', dpi=300)
            print(f"Plot saved to {save_path}")
        except Exception as e:
            warnings.warn(f"Could not save plot: {e}")
    
    # Show or return plot
    if return_panel:
        return app
    elif show_plot:
        try:
            app.show(port=0)  # Auto-select port
        except Exception as e:
            warnings.warn(f"Could not display plot: {e}")
            return app
    
    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run LZSS plots")
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Subparser for cumulative plot
    cumulative_parser = subparsers.add_parser('cumulative', help='Plot cumulative factors')
    cumulative_parser.add_argument('fasta_filepath', help='Path to FASTA file')
    cumulative_parser.add_argument('output_dir', help='Output directory')
    cumulative_parser.add_argument('--max_sequences', type=int, default=None, help='Maximum number of sequences to process')
    cumulative_parser.add_argument('--save_factors_text', action='store_true', help='Save factors as text files')
    cumulative_parser.add_argument('--save_factors_binary', action='store_true', help='Save factors as binary files')

    # Subparser for self-factors-plot
    self_factors_parser = subparsers.add_parser('self-factors-plot', help='Plot self-factors')
    self_factors_parser.add_argument('--fasta_filepath', help='Path to FASTA file')
    self_factors_parser.add_argument('--factors_filepath', help='Path to binary factors file')
    self_factors_parser.add_argument('--name', default=None, help='Name for the plot title')
    self_factors_parser.add_argument('--save_path', default=None, help='Path to save the plot image')
    self_factors_parser.add_argument('--show_plot', action='store_true', default=True, help='Whether to display the plot')
    self_factors_parser.add_argument('--return_panel', action='store_true', help='Whether to return the Panel app')

    # Subparser for reference sequence plotting
    ref_plot_parser = subparsers.add_parser('reference-plot', help='Plot target sequence factorized with reference sequence')
    ref_plot_parser.add_argument('reference_seq', help='Reference DNA sequence')
    ref_plot_parser.add_argument('target_seq', help='Target DNA sequence')
    ref_plot_parser.add_argument('--factors_filepath', help='Path to binary factors file (optional, will compute if not provided)')
    ref_plot_parser.add_argument('--reference_name', default='Reference', help='Name for the reference sequence')
    ref_plot_parser.add_argument('--target_name', default='Target', help='Name for the target sequence')
    ref_plot_parser.add_argument('--save_path', default=None, help='Path to save the plot image')
    ref_plot_parser.add_argument('--show_plot', action='store_true', default=True, help='Whether to display the plot')
    ref_plot_parser.add_argument('--interactive', action='store_true', help='Use interactive Panel/Datashader plot instead of simple matplotlib')
    ref_plot_parser.add_argument('--return_panel', action='store_true', help='Whether to return the Panel app (interactive mode only)')

    args = parser.parse_args()

    if args.command == 'cumulative':
        plot_single_seq_accum_factors_from_file(
            fasta_filepath=args.fasta_filepath,
            output_dir=args.output_dir,
            max_sequences=args.max_sequences,
            save_factors_text=args.save_factors_text,
            save_factors_binary=args.save_factors_binary
        )
    elif args.command == 'self-factors-plot':
        plot_multiple_seq_self_lz_factor_plot_from_file(
            fasta_filepath=args.fasta_filepath,
            factors_filepath=args.factors_filepath,
            name=args.name,
            save_path=args.save_path,
            show_plot=args.show_plot,
            return_panel=args.return_panel
        )
    elif args.command == 'reference-plot':
        # Choose between simple and interactive plotting
        if args.interactive:
            plot_reference_seq_lz_factor_plot(
                reference_seq=args.reference_seq,
                target_seq=args.target_seq,
                factors_filepath=args.factors_filepath,
                reference_name=args.reference_name,
                target_name=args.target_name,
                save_path=args.save_path,
                show_plot=args.show_plot,
                return_panel=args.return_panel
            )
        else:
            plot_reference_seq_lz_factor_plot_simple(
                reference_seq=args.reference_seq,
                target_seq=args.target_seq,
                factors_filepath=args.factors_filepath,
                reference_name=args.reference_name,
                target_name=args.target_name,
                save_path=args.save_path,
                show_plot=args.show_plot
            )