"""
FASTA file plotting utilities.

This module provides functions for creating plots and visualizations
from FASTA files and their factorizations.
"""

from typing import Union, Optional, Dict, Any
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


# Keep old function for backward compatibility

# def plot_multiple_seq_self_weizmann_factor_plot_from_fasta(
#     fasta_filepath: Union[str, Path],
#     name: Optional[str] = None,
#     save_path: Optional[Union[str, Path]] = None,
#     show_plot: bool = True
# ) -> None:
#     """
#     Create a Weizmann factor plot for multiple DNA sequences from a FASTA file.

#     This function reads a FASTA file containing multiple DNA sequences, factorizes them
#     using the multiple DNA with reverse complement algorithm, and creates a specialized
#     plot where each factor is represented as a line. The plot shows the relationship
#     between factor positions and their reference positions.

#     Args:
#         fasta_filepath: Path to the FASTA file containing DNA sequences
#         name: Optional name for the plot title (defaults to FASTA filename)
#         save_path: Optional path to save the plot image
#         show_plot: Whether to display the plot

#     Raises:
#         PlotError: If plotting fails or FASTA file cannot be processed
#         FileNotFoundError: If FASTA file doesn't exist
#     """
#     try:
#         import matplotlib.pyplot as plt
#         import matplotlib.colors as mcolors
#         import numpy as np
#     except ImportError:
#         warnings.warn("matplotlib is required for plotting. Install with: pip install matplotlib", UserWarning)
#         return

#     from .._noLZSS import factorize_fasta_multiple_dna_w_rc

#     fasta_filepath = Path(fasta_filepath)

#     if not fasta_filepath.exists():
#         raise FileNotFoundError(f"FASTA file not found: {fasta_filepath}")

#     # Determine plot title
#     if name is None:
#         name = fasta_filepath.stem

#     try:
#         # Get factors from FASTA file
#         print(f"Reading and factorizing sequences from {fasta_filepath}...")
#         factors = factorize_fasta_multiple_dna_w_rc(str(fasta_filepath))

#         print(f"Preparing plot for {len(factors)} factors...")
#         if not factors:
#             raise PlotError("No factors found in FASTA file")

#         # Extract factor data
#         positions = []
#         lengths = []
#         refs = []
#         is_rcs = []

#         for factor in factors:
#             if len(factor) == 4:  # (start, length, ref, is_rc) tuple
#                 start, length, ref, is_rc = factor
#             else:  # Assume (start, length, ref) format, default is_rc to False
#                 start, length, ref = factor
#                 is_rc = False

#             positions.append(start)
#             lengths.append(length)
#             refs.append(ref)
#             is_rcs.append(is_rc)

#         # Convert to numpy arrays for easier processing
#         positions = np.array(positions)
#         lengths = np.array(lengths)
#         refs = np.array(refs)
#         is_rcs = np.array(is_rcs)

#         # Create the plot
#         fig, ax = plt.subplots(figsize=(12, 12))

#         # Calculate color intensities based on factor lengths
#         # Normalize lengths to [0, 1] for color scaling
#         if len(lengths) > 1:
#             # Use log scale for better visualization of length distribution
#             log_lengths = np.log(lengths + 1)  # +1 to avoid log(0)
#             norm_lengths = (log_lengths - log_lengths.min()) / (log_lengths.max() - log_lengths.min())
#         else:
#             norm_lengths = np.array([0.5])  # Default for single factor

#         # Additionally consider position to vary color intensity
#         if len(positions) > 1:
#             norm_positions = (positions - positions.min()) / (positions.max() - positions.min())
#         else:
#             norm_positions = np.array([0.5])  # Default for single factor

#         # Plot each factor as a line
#         for i, (pos, length, ref, is_rc, norm_len, norm_pos) in enumerate(zip(positions, lengths, refs, is_rcs, norm_lengths, norm_positions)):
#             if is_rc:
#                 # Reverse complement: red line
#                 # x_init = pos, x_final = pos + length
#                 # y_init = ref + length, y_final = ref
#                 x_coords = [pos, pos + length]
#                 y_coords = [ref + length, ref]

#                 # Red color with intensity based on length (more opaque for longer factors)
#                 alpha = min((norm_len * 0.5 + norm_pos * 0.5), 1.0)  # Combine length and position for alpha
#                 color = (1.0, 0.0, 0.0, alpha)  # Red with variable alpha
#             else:
#                 # Forward: blue line
#                 # x_init = pos, x_final = pos + length
#                 # y_init = ref, y_final = ref + length
#                 x_coords = [pos, pos + length]
#                 y_coords = [ref, ref + length]

#                 # Blue color with intensity based on length (more opaque for longer factors)
#                 alpha = min((norm_len * 0.5 + norm_pos * 0.5), 1.0)  # Combine length and position for alpha
#                 color = (0.0, 0.0, 1.0, alpha)  # Blue with variable alpha

#             # Plot the line
#             ax.plot(x_coords, y_coords, color=color, linewidth=1.5, alpha=alpha)

#         # Add diagonal line y=x for reference
#         max_val = max(positions.max() + lengths.max(), refs.max() + lengths.max())
#         ax.plot([0, max_val], [0, max_val], color='gray', linestyle='--', linewidth=1, alpha=0.5)

#         # Set axis labels and title
#         ax.set_xlabel(f'Position in concatenated sequence ({name})')
#         ax.set_ylabel(f'Reference position ({name})')
#         ax.set_title(f'Weizmann Factor Plot - {name}')

#         # Make axes equal for better visualization
#         ax.set_aspect('equal', adjustable='box')

#         # Add grid
#         ax.grid(True, alpha=0.3)

#         # Add legend
#         from matplotlib.patches import Patch
#         legend_elements = [
#             Patch(facecolor='blue', alpha=0.8, label='Forward factors'),
#             Patch(facecolor='red', alpha=0.8, label='Reverse complement factors')
#         ]
#         ax.legend(handles=legend_elements, loc='upper right')

#         # Adjust layout
#         plt.tight_layout()

#         # Save plot if requested
#         if save_path:
#             save_path = Path(save_path)
#             save_path.parent.mkdir(parents=True, exist_ok=True)
#             plt.savefig(save_path, dpi=300, bbox_inches='tight')
#             print(f"Plot saved to {save_path}")

#         # Show plot
#         if show_plot:
#             plt.show()
#         else:
#             plt.close(fig)

#     except Exception as e:
#         raise PlotError(f"Failed to create Weizmann factor plot: {e}")

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