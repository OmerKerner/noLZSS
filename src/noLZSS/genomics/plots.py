"""
FASTA file plotting utilities.

This module provides functions for creating plots and visualizations
from FASTA files and their factorizations.
"""

from typing import Union, Optional, Dict, Any, List, Tuple, Literal, Sequence, cast
from pathlib import Path
import warnings
import argparse

from ..utils import NoLZSSError
from .fasta import _parse_fasta_content
from .sequences import detect_sequence_type



class PlotError(NoLZSSError):
    """Raised when plotting operations fail."""
    pass


def _normalize_reference_factors(
    factors: List[Tuple[int, ...]]
) -> List[Tuple[int, int, int, bool]]:
    """Ensure reference-sequence factors include an explicit reverse-complement flag."""
    normalized: List[Tuple[int, int, int, bool]] = []

    for idx, factor in enumerate(factors):
        if len(factor) == 4:
            start, length, ref, is_rc = factor
            normalized.append((int(start), int(length), int(ref), bool(is_rc)))
        elif len(factor) == 3:
            start, length, ref = factor
            normalized.append((int(start), int(length), int(ref), False))
        else:
            raise PlotError(
                f"Factor at index {idx} has {len(factor)} elements; expected 3 or 4 values"
            )

    return normalized


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
) -> Optional[Any]:
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
        save_path: Optional path to save the plot (supports .html or .png; PNG export
            requires optional selenium-based dependencies)
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
            target_path = Path(save_path)
            try:
                target_path.parent.mkdir(parents=True, exist_ok=True)
                suffix = target_path.suffix.lower()

                if suffix in {".html", ".htm"}:
                    try:
                        from bokeh.resources import INLINE  # type: ignore[import-not-found]
                    except ImportError as exc:  # pragma: no cover - bokeh required earlier
                        raise PlotError(
                            "Saving to HTML requires bokeh to be installed"
                        ) from exc

                    app.save(str(target_path), embed=False, resources=INLINE, title=f"LZ Factor Plot - {name}")
                    print(f"Saved interactive HTML plot to {target_path}")

                elif suffix == ".png":
                    try:
                        from panel.io.save import save as panel_save_png  # type: ignore[import-not-found]
                    except ImportError as exc:
                        raise PlotError(
                            "PNG export requires selenium and a web browser driver. "
                            "Install with: pip install selenium and either "
                            "'conda install -c conda-forge firefox geckodriver' or "
                            "'conda install -c conda-forge chromium chromedriver'"
                        ) from exc

                    try:
                        panel_save_png(app, filename=str(target_path), as_png=True)
                    except Exception as exc:
                        raise PlotError(
                            "Failed to export PNG. Ensure selenium and a compatible web driver are available."
                        ) from exc

                    print(f"Saved PNG snapshot to {target_path}")

                else:
                    raise PlotError(
                        f"Unsupported save_path extension '{suffix}'. Use .html or .png"
                    )

            except Exception as exc:
                raise PlotError(f"Failed to save plot to {target_path}: {exc}") from exc

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


def plot_multiple_seq_self_lz_factor_plot_simple(
    fasta_filepath: Optional[Union[str, Path]] = None,
    factors_filepath: Optional[Union[str, Path]] = None,
    name: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True
) -> None:
    """
    Create a simple matplotlib factor plot for multiple DNA sequences from a FASTA file or binary factors file.

    This function reads factors either from a FASTA file (by factorizing multiple DNA sequences)
    or from an enhanced binary factors file with metadata. It creates a static plot using matplotlib
    with sequence boundaries visualization - a simplified alternative to the interactive Panel/Datashader version.

    Args:
        fasta_filepath: Path to the FASTA file containing DNA sequences (mutually exclusive with factors_filepath)
        factors_filepath: Path to binary factors file with metadata (mutually exclusive with fasta_filepath)
        name: Optional name for the plot title (defaults to input filename)
        save_path: Optional path to save the plot image (PNG, PDF, SVG, etc.)
        show_plot: Whether to display the plot

    Raises:
        PlotError: If plotting fails or input files cannot be processed
        FileNotFoundError: If input file doesn't exist
        ImportError: If matplotlib is not available
        ValueError: If both or neither input files are provided
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
    except ImportError as e:
        missing_dep = str(e).split("'")[1] if "'" in str(e) else str(e)
        raise ImportError(
            f"Missing required dependency: {missing_dep}. "
            f"Install with: pip install matplotlib"
        )

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
            max_x = max(factor[0] + factor[1] for factor in factors)  # max end position (start + length)
            max_y = max(factor[2] + (factor[1] if not factor[3] else 0) for factor in factors)  # max ref position
            max_pos = max(max_x, max_y)
            sequence_boundaries.append((prev_pos, max_pos, last_name))
        else:
            # No sentinels - single sequence
            seq_name = sequence_names[0] if sequence_names else "sequence"
            max_x = max(factor[0] + factor[1] for factor in factors)
            max_y = max(factor[2] + (factor[1] if not factor[3] else 0) for factor in factors)
            max_pos = max(max_x, max_y)
            sequence_boundaries.append((0, max_pos, seq_name))

        print(f"Sequence boundaries: {sequence_boundaries}")
        print(f"Sentinel positions: {sentinel_positions}")

        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # Plot factors with different colors for forward/reverse
        for start, length, ref, is_rc in factors:
            x0 = start
            x1 = start + length
            
            if is_rc:
                # Reverse complement: y0 = ref + length, y1 = ref
                y0 = ref + length
                y1 = ref
                color = 'red'
            else:
                # Forward: y0 = ref, y1 = ref + length
                y0 = ref
                y1 = ref + length
                color = 'blue'
            
            # Draw line segment
            ax.plot([x0, x1], [y0, y1], color=color, alpha=0.6, linewidth=1.5)
        
        # Add diagonal reference line
        max_pos = max(sequence_boundaries[-1][1] if sequence_boundaries else 1000,
                     max(factor[0] + factor[1] for factor in factors))
        ax.plot([0, max_pos], [0, max_pos], 'gray', linestyle='--', alpha=0.5, linewidth=1)
        
        # Add sentinel boundary lines
        for pos in sentinel_positions:
            ax.axvline(x=pos, color='green', linewidth=2, alpha=0.7, linestyle='solid')
            ax.axhline(y=pos, color='green', linewidth=2, alpha=0.7, linestyle='solid')
        
        # Add sequence region backgrounds (alternating light colors)
        colors_cycle = ['lightblue', 'lightcoral', 'lightgreen', 'lightyellow', 'lightpink']
        for idx, (start_pos, end_pos, seq_name) in enumerate(sequence_boundaries):
            bg_color = colors_cycle[idx % len(colors_cycle)]
            ax.axvspan(start_pos, end_pos, alpha=0.1, color=bg_color)
            ax.axhspan(start_pos, end_pos, alpha=0.1, color=bg_color)
        
        # Add sequence name labels
        label_offset = max_pos * 0.02
        for start_pos, end_pos, seq_name in sequence_boundaries:
            mid_pos = (start_pos + end_pos) / 2
            # X-axis label (bottom)
            ax.text(mid_pos, -label_offset, seq_name, ha='center', va='top',
                   fontsize=10, color='darkblue', weight='bold')
            # Y-axis label (left side)
            ax.text(-label_offset, mid_pos, seq_name, ha='right', va='center',
                   fontsize=10, color='darkblue', weight='bold', rotation=90)
        
        # Set labels and title
        sequence_list = ", ".join([b[2] for b in sequence_boundaries])
        ax.set_xlabel(f'Position in concatenated sequence ({name})', fontsize=12)
        ax.set_ylabel(f'Reference position ({name})', fontsize=12)
        ax.set_title(f'LZ Factor Plot - {name}\nSequences: {sequence_list}', fontsize=14, weight='bold')
        
        # Add legend
        legend_elements = [
            patches.Patch(color='blue', alpha=0.6, label='Forward factors'),
            patches.Patch(color='red', alpha=0.6, label='Reverse complement factors'),
            patches.Patch(color='green', alpha=0.7, label='Sequence boundaries')
        ]
        ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
        
        # Set equal aspect ratio, grid, and ensure axes start at zero
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(-label_offset * 2, max_pos)
        ax.set_ylim(-label_offset * 2, max_pos)
        
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

    except Exception as e:
        raise PlotError(f"Failed to create matplotlib LZ factor plot: {e}")


def plot_reference_seq_lz_factor_plot_simple(
    reference_seq: Optional[Union[str, bytes]] = None,
    target_seq: Optional[Union[str, bytes]] = None,
    factors: Optional[List[Tuple[int, int, int, bool]]] = None,
    factors_filepath: Optional[Union[str, Path]] = None,
    reference_name: str = "Reference",
    target_name: str = "Target",
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    factorization_mode: Literal["dna", "general"] = "dna"
) -> None:
    """
    Create a simple matplotlib factor plot for a sequence factorized with a reference sequence.
    
    This function creates a plot compatible with the outputs of factorize_dna_w_reference_seq()
    or the general factorize_w_reference() wrapper.
    The plot shows the reference sequence at the beginning, concatenated with the target sequence,
    and uses distinct colors for reference vs target regions.
    
    Args:
        reference_seq: Reference DNA sequence (A, C, T, G - case insensitive) or
            general ASCII text when ``factorization_mode`` is "general". Optional if
            ``factors_filepath`` is provided (parameters will be inferred from factors).
        target_seq: Target DNA sequence (A, C, T, G - case insensitive) or
            general ASCII text when ``factorization_mode`` is "general". Optional if
            ``factors_filepath`` is provided (parameters will be inferred from factors).
        factors: Optional list of (start, length, ref, is_rc) tuples from
            factorize_dna_w_reference_seq() or factorize_w_reference(). If None,
            the function will compute factors automatically based on
            ``factorization_mode``.
        factors_filepath: Optional path to binary factors file (mutually exclusive
            with ``factors``). When provided and sequences are not, parameters are
            inferred from the factors: first factor start = target_start, last factor
            end = total_length.
        reference_name: Name for the reference sequence (default: "Reference")
        target_name: Name for the target sequence (default: "Target")
        save_path: Optional path to save the plot image
        show_plot: Whether to display the plot
        factorization_mode: Choose "dna" for reverse-complement-aware
            factorization or "general" for ASCII/general sequences without
            reverse complements
        
    Raises:
        PlotError: If plotting fails or input sequences are invalid
        ValueError: If both factors and factors_filepath are provided, or if
            sequences are not provided and factors_filepath is not provided
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
    
    # If sequences are not provided, factors_filepath must be provided
    sequences_provided = reference_seq is not None and target_seq is not None
    if not sequences_provided and factors_filepath is None and factors is None:
        raise ValueError("Either provide reference_seq and target_seq, or provide factors_filepath or factors")

    factorization_mode_normalized = factorization_mode.lower()
    if factorization_mode_normalized not in {"dna", "general"}:
        raise ValueError("factorization_mode must be 'dna' or 'general'")
    
    # Convert sequences to strings if they're bytes
    if reference_seq is not None and isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if target_seq is not None and isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    raw_factors: Optional[List[Tuple[int, ...]]] = None
    metadata: Optional[Dict[str, Any]] = None
    sequence_names: List[str] = []
    sentinel_factor_indices: List[int] = []

    # Get factors if not provided
    if factors is None:
        if factors_filepath is not None:
            from ..utils import read_factors_binary_file_with_metadata

            try:
                metadata = read_factors_binary_file_with_metadata(factors_filepath)
                raw_factors = metadata.get('factors', [])
                sequence_names = metadata.get('sequence_names', [])
                sentinel_factor_indices = metadata.get('sentinel_factor_indices', [])
            except NoLZSSError as exc:
                raise PlotError(
                    f"Failed to read factors metadata from {factors_filepath}: {exc}"
                ) from exc
        else:
            if factorization_mode_normalized == "dna":
                from .sequences import factorize_dna_w_reference_seq

                assert reference_seq is not None and target_seq is not None
                raw_factors = factorize_dna_w_reference_seq(
                    cast(str, reference_seq),
                    cast(str, target_seq)
                )
            else:
                from ..core import factorize_w_reference

                assert reference_seq is not None and target_seq is not None
                raw_factors = factorize_w_reference(
                    cast(str, reference_seq),
                    cast(str, target_seq)
                )
    else:
        raw_factors = list(factors)

    if not raw_factors:
        raise PlotError("No factors available for plotting")

    factors = _normalize_reference_factors(raw_factors)
    
    # Calculate sequence boundaries
    if sequences_provided:
        # Use provided sequences
        assert reference_seq is not None and target_seq is not None  # Type hint for mypy
        ref_length = len(reference_seq)
        target_start = ref_length + 1  # +1 for sentinel
        target_length = len(target_seq)
        total_length = target_start + target_length
    else:
        # Infer boundaries from factors
        if not factors:
            raise PlotError("Cannot infer sequence boundaries from empty factors")
        
        # First factor starts at target_start, last factor ends at total_length
        target_start = min(start for start, length, ref, is_rc in factors)
        total_length = max(start + length for start, length, ref, is_rc in factors)
        ref_length = target_start - 1  # -1 because target_start = ref_length + 1
        target_length = total_length - target_start
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot factors with different colors for forward/reverse and reference/target
    for start, length, ref_pos, is_rc in factors:
        end = start + length
        
        # Calculate coordinates based on direction
        if is_rc:
            # Reverse complement: y0 = ref + length, y1 = ref
            y_start = ref_pos + length
            y_end = ref_pos
        else:
            # Forward: y0 = ref, y1 = ref + length
            y_start = ref_pos
            y_end = ref_pos + length
        
        # Determine color based on region and orientation
        if start >= target_start:
            # Target region
            color = 'red' if not is_rc else 'darkorange'
            alpha = 0.7
        else:
            # Reference region (shouldn't happen with factorize_dna_w_reference_seq)
            color = 'blue' if not is_rc else 'darkblue'
            alpha = 0.7
        
        # Draw line segment
        ax.plot([start, end], [y_start, y_end], color=color, alpha=alpha, linewidth=2)
    
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
        patches.Patch(color='blue', alpha=0.7, label=f'{target_name} forward factors'),
        patches.Patch(color='darkred', alpha=0.7, label=f'{target_name} reverse complement factors'),
        patches.Patch(color='blue', alpha=0.1, label=f'{reference_name} region'),
        patches.Patch(color='red', alpha=0.1, label=f'{target_name} region'),
        patches.Patch(color='green', alpha=0.8, label='Sequence boundary')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    # Set equal aspect ratio, grid, and ensure axes start at zero
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max_pos)
    ax.set_ylim(0, max_pos)

    # Add text annotations for sequence regions within the visible range
    label_offset = max(max_pos * 0.03, 1.0)
    ax.text(ref_length / 2, label_offset, reference_name, ha='center', va='bottom',
        fontsize=12, color='blue', weight='bold')
    ax.text((target_start + total_length) / 2, label_offset, target_name, ha='center', va='bottom',
        fontsize=12, color='red', weight='bold')
    ax.text(label_offset, ref_length / 2, reference_name, ha='left', va='center',
        fontsize=12, color='blue', weight='bold', rotation=90)
    ax.text(label_offset, (target_start + total_length) / 2, target_name, ha='left', va='center',
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
    reference_seq: Optional[Union[str, bytes]] = None,
    target_seq: Optional[Union[str, bytes]] = None,
    factors: Optional[List[Tuple[int, int, int, bool]]] = None,
    factors_filepath: Optional[Union[str, Path]] = None,
    reference_name: str = "Reference",
    target_name: str = "Target",
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    return_panel: bool = False,
    factorization_mode: Literal["dna", "general"] = "dna"
) -> Optional[Any]:
    """
    Create an interactive Datashader/Panel factor plot for a sequence factorized with a reference sequence.
    
    This function creates a plot compatible with the outputs of factorize_dna_w_reference_seq()
    or the general factorize_w_reference() wrapper.
    The plot shows the reference sequence at the beginning, concatenated with the target sequence,
    and uses distinct colors for reference vs target regions.
    
    Args:
        reference_seq: Reference DNA sequence (A, C, T, G - case insensitive) or
            general ASCII text when ``factorization_mode`` is "general". Optional if
            ``factors_filepath`` is provided (parameters will be inferred from factors).
        target_seq: Target DNA sequence (A, C, T, G - case insensitive) or
            general ASCII text when ``factorization_mode`` is "general". Optional if
            ``factors_filepath`` is provided (parameters will be inferred from factors).
        factors: Optional list of (start, length, ref, is_rc) tuples from
            factorize_dna_w_reference_seq() or factorize_w_reference(). If None,
            the function will compute factors automatically based on
            ``factorization_mode``.
        factors_filepath: Optional path to binary factors file (mutually exclusive
            with ``factors``). When provided and sequences are not, parameters are
            inferred from the factors: first factor start = target_start, last factor
            end = total_length.
        reference_name: Name for the reference sequence (default: "Reference")
        target_name: Name for the target sequence (default: "Target")
        save_path: Optional path to save the plot image (PNG export)
        show_plot: Whether to display/serve the plot
        return_panel: Whether to return the Panel app for embedding
        factorization_mode: Choose "dna" for reverse-complement-aware
            factorization or "general" for ASCII/general sequences without
            reverse complements
        
    Returns:
        Panel app if return_panel=True, otherwise None
        
    Raises:
        PlotError: If plotting fails or input sequences are invalid
        ValueError: If both factors and factors_filepath are provided, or if
            sequences are not provided and factors_filepath is not provided
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
    
    # If sequences are not provided, factors_filepath must be provided
    sequences_provided = reference_seq is not None and target_seq is not None
    if not sequences_provided and factors_filepath is None and factors is None:
        raise ValueError("Either provide reference_seq and target_seq, or provide factors_filepath or factors")

    factorization_mode_normalized = factorization_mode.lower()
    if factorization_mode_normalized not in {"dna", "general"}:
        raise ValueError("factorization_mode must be 'dna' or 'general'")

    # Convert sequences to strings if they're bytes
    if reference_seq is not None and isinstance(reference_seq, bytes):
        reference_seq = reference_seq.decode('ascii')
    if target_seq is not None and isinstance(target_seq, bytes):
        target_seq = target_seq.decode('ascii')
    
    raw_factors: Optional[List[Tuple[int, ...]]] = None
    metadata: Optional[Dict[str, Any]] = None
    sequence_names: List[str] = []
    sentinel_factor_indices: List[int] = []

    # Get factors if not provided
    if factors is None:
        if factors_filepath is not None:
            from ..utils import read_factors_binary_file_with_metadata

            try:
                metadata = read_factors_binary_file_with_metadata(factors_filepath)
                raw_factors = metadata.get('factors', [])
                sequence_names = metadata.get('sequence_names', [])
                sentinel_factor_indices = metadata.get('sentinel_factor_indices', [])
            except NoLZSSError as exc:
                raise PlotError(
                    f"Failed to read factors metadata from {factors_filepath}: {exc}"
                ) from exc
        else:
            if factorization_mode_normalized == "dna":
                from .sequences import factorize_dna_w_reference_seq

                assert reference_seq is not None and target_seq is not None
                raw_factors = factorize_dna_w_reference_seq(
                    cast(str, reference_seq),
                    cast(str, target_seq)
                )
            else:
                from ..core import factorize_w_reference

                assert reference_seq is not None and target_seq is not None
                raw_factors = factorize_w_reference(
                    cast(str, reference_seq),
                    cast(str, target_seq)
                )
    else:
        raw_factors = list(factors)

    if not raw_factors:
        raise PlotError("No factors available for plotting")

    factors = _normalize_reference_factors(raw_factors)

    valid_sentinels = [idx for idx in sentinel_factor_indices if 0 <= idx < len(factors)]
    sentinel_positions = sorted(factors[idx][0] for idx in valid_sentinels)

    factors_for_plot = [factor for idx, factor in enumerate(factors) if idx not in valid_sentinels]
    if not factors_for_plot:
        factors_for_plot = factors

    reference_label = reference_name
    target_label = target_name

    if not sequences_provided and sequence_names:
        if len(sequence_names) >= 1 and reference_name == "Reference":
            reference_label = sequence_names[0]
        if len(sequence_names) >= 2 and target_name == "Target":
            target_label = sequence_names[1]

    max_end = max((start + length) for start, length, _, _ in factors_for_plot) if factors_for_plot else 0

    if sequences_provided:
        assert reference_seq is not None and target_seq is not None  # mypy hint
        ref_length = len(reference_seq)
        target_length = len(target_seq)
        target_start = ref_length + 1
        total_length = target_start + target_length
        sequence_boundaries = [
            (0, ref_length, reference_label),
            (target_start, total_length, target_label)
        ]
    elif sentinel_positions:
        sequence_boundaries = []
        prev_pos = 0
        for idx, pos in enumerate(sentinel_positions):
            seq_name = sequence_names[idx] if idx < len(sequence_names) else (
                reference_label if idx == 0 else f"sequence_{idx}"
            )
            sequence_boundaries.append((prev_pos, pos, seq_name))
            prev_pos = pos + 1

        max_pos = max(max_end, prev_pos)
        last_index = len(sequence_boundaries)
        if sequence_names and last_index < len(sequence_names):
            last_name = sequence_names[last_index]
        elif last_index == 0:
            last_name = reference_label
        elif last_index == 1:
            last_name = target_label
        else:
            last_name = f"sequence_{last_index}"

        sequence_boundaries.append((prev_pos, max_pos, last_name))

        total_length = int(sequence_boundaries[-1][1])
        ref_length = int(sequence_boundaries[0][1]) if sequence_boundaries else 0
        target_start = int(sequence_boundaries[1][0]) if len(sequence_boundaries) > 1 else ref_length + 1
        target_length = max(total_length - target_start, 0)
    else:
        if not factors_for_plot:
            raise PlotError("No factors available for plotting")

        target_start = min(start for start, length, _, _ in factors_for_plot)
        total_length = max(start + length for start, length, _, _ in factors_for_plot)
        ref_length = max(target_start - 1, 0)
        target_length = max(total_length - target_start, 0)
        sequence_boundaries = [
            (0, ref_length, reference_label),
            (target_start, total_length, target_label)
        ]

    if sequence_boundaries:
        reference_label = sequence_boundaries[0][2]
        if len(sequence_boundaries) > 1:
            target_label = sequence_boundaries[1][2]

    if not factors_for_plot:
        raise PlotError("No valid factors to plot")

    factor_data = []
    for start, length, ref_pos, is_rc in factors_for_plot:
        end = start + length

        # Calculate coordinates based on direction
        if is_rc:
            # Reverse complement: y0 = ref + length, y1 = ref
            y0 = ref_pos + length
            y1 = ref_pos
        else:
            # Forward: y0 = ref, y1 = ref + length
            y0 = ref_pos
            y1 = ref_pos + length

        if start >= target_start:
            region = 'target'
            color = 'darkred' if is_rc else 'blue'
        else:
            region = 'reference'
            color = 'red' if is_rc else 'blue'

        factor_data.append({
            'x0': start,
            'x1': end,
            'y0': y0,
            'y1': y1,
            'length': length,
            'is_rc': is_rc,
            'region': region,
            'color': color,
            'start': start,
            'end': end,
            'ref_pos': ref_pos
        })

    df = pd.DataFrame(factor_data)

    if df.empty:
        raise PlotError("No valid factors to plot")

    global_x_max = max(float(df[['x0', 'x1']].max().max()), float(total_length))
    raw_y_max = float(df[['y0', 'y1']].max().max())
    global_y_max = max(raw_y_max, float(target_start))
    global_max = max(global_x_max, global_y_max)

    plot_name = (
        f"{reference_label} vs {target_label}"
        if len(sequence_boundaries) >= 2
        else reference_label
    )
    
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
                    ref_forward,
                    kdims=['x0', 'y0', 'x1', 'y1'],
                    vdims=['length', 'is_rc', 'region', 'start', 'end', 'ref_pos'],
                    label=f"{reference_label} forward"
                ).opts(color='blue', alpha=0.7, line_width=2)
                layers.append(ref_forward_segments)
            
            if len(ref_reverse) > 0:
                ref_reverse_segments = hv.Segments(
                    ref_reverse,
                    kdims=['x0', 'y0', 'x1', 'y1'],
                    vdims=['length', 'is_rc', 'region', 'start', 'end', 'ref_pos'],
                    label=f"{reference_label} reverse"
                ).opts(color='red', alpha=0.7, line_width=2)
                layers.append(ref_reverse_segments)
        
        # Target factors - use red/orange tones  
        if len(df_target) > 0:
            target_forward = df_target[~df_target['is_rc']]
            target_reverse = df_target[df_target['is_rc']]
            
            if len(target_forward) > 0:
                target_forward_segments = hv.Segments(
                    target_forward,
                    kdims=['x0', 'y0', 'x1', 'y1'],
                    vdims=['length', 'is_rc', 'region', 'start', 'end', 'ref_pos'],
                    label=f"{target_label} forward"
                ).opts(color='blue', alpha=0.7, line_width=2)
                layers.append(target_forward_segments)
            
            if len(target_reverse) > 0:
                target_reverse_segments = hv.Segments(
                    target_reverse,
                    kdims=['x0', 'y0', 'x1', 'y1'],
                    vdims=['length', 'is_rc', 'region', 'start', 'end', 'ref_pos'],
                    label=f"{target_label} reverse"
                ).opts(color='darkred', alpha=0.7, line_width=2)
                layers.append(target_reverse_segments)
        
        if not layers:
            return hv.Text(0, 0, "No data").opts(width=800, height=800)
        
        # Combine layers
        overlay = hv.Overlay(layers).opts(show_legend=True)
        return overlay
    
    def create_hover_overlay(x_range, y_range, df_filtered, k_per_bin):
        """Create hover overlay for detailed factor information"""
        if x_range is None or y_range is None or len(df_filtered) == 0:
            return hv.Segments([])
        
        x_min, x_max = x_range
        y_min, y_max = y_range
        
        # Filter data to current view (allow partial overlaps)
        view_mask = (
            (df_filtered['x1'] >= x_min) &
            (df_filtered['x0'] <= x_max) &
            (df_filtered['y1'] >= y_min) &
            (df_filtered['y0'] <= y_max)
        )
        view_data = df_filtered[view_mask].copy()

        if len(view_data) == 0:
            return hv.Segments([])

        if len(view_data) > k_per_bin:
            view_data = view_data.nlargest(int(k_per_bin), 'length')

        view_data['ref_pos'] = view_data['y0']

        hover_segments = hv.Segments(
            view_data,
            kdims=['x0', 'y0', 'x1', 'y1'],
            vdims=['start', 'end', 'length', 'ref_pos', 'is_rc', 'region', 'color']
        ).opts(
            color='color',
            line_width=4,
            alpha=0.9,
            tools=['hover'],
            show_legend=False
        )

        return hover_segments
    
    # Set up interactive plot
    initial_padding = max(global_max * 0.05, 10.0)
    padded_min = -initial_padding
    rangexy = streams.RangeXY(x_range=(padded_min, global_x_max), y_range=(padded_min, global_max))
    
    # Create dynamic plot function
    def create_plot(length_range, show_overlay, k_per_bin, colormap_name):
        length_min: float
        length_max: float

        if isinstance(length_range, Sequence) and len(length_range) == 2:
            length_min, length_max = length_range[0], length_range[1]
        else:
            start = getattr(length_range, "start", None)
            end = getattr(length_range, "end", None)
            if start is None or end is None:
                raise TypeError(
                    "length_range must be a (min, max) sequence or expose 'start' and 'end' attributes"
                )
            length_min, length_max = start, end

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
        filtered_max = float(df_filtered[['x1', 'y1']].max().max())
        max_val = max(global_max, filtered_max)
        min_val = 0.0
        range_span = max(max_val - min_val, 1.0)
        label_offset = max(range_span * 0.03, 1.0)

        diagonal = hv.Curve([(min_val, min_val), (max_val, max_val)]).opts(
            line_dash='dashed',
            line_color='gray',
            line_width=1,
            alpha=0.5
        )
        
        # Add sequence boundary lines
        boundary_elements = []
        if sentinel_positions:
            boundary_positions = sentinel_positions
        else:
            boundary_positions = [start_pos - 0.5 for start_pos, _, _ in sequence_boundaries[1:]]

        for pos in boundary_positions:
            if min_val <= pos <= max_val:
                boundary_elements.append(
                    hv.VLine(pos).opts(
                        line_color='green',
                        line_width=3,
                        alpha=0.8,
                        line_dash='solid'
                    )
                )
                boundary_elements.append(
                    hv.HLine(pos).opts(
                        line_color='green',
                        line_width=3,
                        alpha=0.8,
                        line_dash='solid'
                    )
                )

        # Add sequence name labels
        label_elements = []
        sequence_colors = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#9467bd']
        for idx, (start_pos, end_pos, seq_name) in enumerate(sequence_boundaries):
            mid_pos = (start_pos + end_pos) / 2
            if min_val <= mid_pos <= max_val:
                label_color = sequence_colors[idx % len(sequence_colors)]

                x_label = hv.Text(mid_pos, min_val - label_offset, seq_name).opts(
                    text_color=label_color,
                    text_font_size='12pt',
                    text_align='center',
                    text_baseline='top'
                )
                label_elements.append(x_label)

                y_label = hv.Text(min_val - label_offset, mid_pos, seq_name).opts(
                    text_color=label_color,
                    text_font_size='12pt',
                    text_align='center',
                    text_baseline='middle',
                    angle=90
                )
                label_elements.append(y_label)
        
        # Combine all plot elements
        plot = cast(Any, base_plot) * diagonal
        
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
        lower_bound = min(-label_offset * 1.5, 0.0)

        plot = plot.opts(
            width=800,
            height=800,
            aspect='equal',
            xlabel=f'Position in concatenated sequence ({plot_name})',
            ylabel=f'Reference position ({plot_name})',
            title=f'Reference Sequence LZ Factor Plot - {plot_name}',
            toolbar='above',
            xlim=(lower_bound, global_x_max),
            ylim=(lower_bound, global_max),
            legend_position='top_left'
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
    
    # Bind widgets to plotting function
    interactive_plot = pn.bind(
        create_plot,
        length_range=length_slider,
        show_overlay=overlay_toggle,
        k_per_bin=k_spinner,
        colormap_name=colormap_select
    )
    
    # Create layout
    sequence_info_lines = []
    for start_pos, end_pos, seq_name in sequence_boundaries:
        seq_length = int(end_pos - start_pos)
        sequence_info_lines.append(f"- {seq_name}: length {seq_length}")

    dataset_info = "\n".join([
        "**Dataset Info:**",
        f"- Total factors: {len(df)}",
        *sequence_info_lines
    ])

    controls = pn.Column(
        pn.pane.Markdown("### Plot Controls"),
        length_slider,
        overlay_toggle,
        k_spinner,
        colormap_select,
        pn.pane.Markdown(dataset_info),
        width=300
    )
    
    plot_pane = pn.pane.HoloViews(interactive_plot, width=850, height=850)
    
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


def plot_factor_length_ccdf(
    factors_filepath: Union[str, Path],
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    separate: bool = True
) -> None:
    """
    Create an empirical CCDF plot of factor lengths on log-log axes from a binary factors file.

    This function reads factors from a binary file and plots the complementary cumulative
    distribution function (CCDF) of factor lengths. Forward and reverse complement factors
    can be plotted separately or together on the same axes with different colors.

    Args:
        factors_filepath: Path to binary factors file with metadata
        save_path: Optional path to save the plot image (PNG, PDF, SVG, etc.)
        show_plot: Whether to display the plot
        separate: Whether to plot forward and reverse complement factors separately
            (default: True). If False, both are plotted on the same axes with different colors.

    Raises:
        PlotError: If file reading or plotting fails
        FileNotFoundError: If factors file doesn't exist
        ImportError: If matplotlib is not available
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        missing_dep = str(e).split("'")[1] if "'" in str(e) else str(e)
        raise ImportError(
            f"Missing required dependency: {missing_dep}. "
            f"Install with: pip install matplotlib"
        )

    from ..utils import read_factors_binary_file_with_metadata

    # Validate input
    factors_filepath = Path(factors_filepath)
    if not factors_filepath.exists():
        raise FileNotFoundError(f"Factors file not found: {factors_filepath}")

    try:
        # Read factors from binary file
        metadata = read_factors_binary_file_with_metadata(factors_filepath)
        factors = metadata['factors']

        if not factors:
            raise PlotError("No factors found in the binary file")

        # Extract lengths and separate by direction
        forward_lengths = []
        reverse_lengths = []

        for factor in factors:
            if len(factor) >= 3:
                start, length, ref = factor[:3]
                is_rc = factor[3] if len(factor) >= 4 else False

                if is_rc:
                    reverse_lengths.append(length)
                else:
                    forward_lengths.append(length)

        if not forward_lengths and not reverse_lengths:
            raise PlotError("No valid factors with lengths found")

        # Function to compute empirical CCDF
        def compute_ccdf(lengths):
            if not lengths:
                return np.array([]), np.array([])
            sorted_lengths = np.sort(lengths)
            unique_lengths = np.unique(sorted_lengths)
            ccdf_values = []

            total_count = len(lengths)
            for length in unique_lengths:
                # Count how many factors have length >= this value
                count_ge = np.sum(sorted_lengths >= length)
                ccdf_values.append(count_ge / total_count)

            return unique_lengths, np.array(ccdf_values)

        # Create the plot
        if separate:
            # Create subplots
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

            # Forward factors
            if forward_lengths:
                lengths_fwd, ccdf_fwd = compute_ccdf(forward_lengths)
                ax1.loglog(lengths_fwd, ccdf_fwd, 'o-', color='blue', alpha=0.7,
                          label=f'Forward (n={len(forward_lengths)})')
                ax1.set_xlabel('Factor Length')
                ax1.set_ylabel('CCDF (P(Length ≥ x))')
                ax1.set_title('Forward Factors')
                ax1.grid(True, alpha=0.3)
                ax1.legend()

            # Reverse complement factors
            if reverse_lengths:
                lengths_rev, ccdf_rev = compute_ccdf(reverse_lengths)
                ax2.loglog(lengths_rev, ccdf_rev, 'o-', color='red', alpha=0.7,
                          label=f'Reverse Complement (n={len(reverse_lengths)})')
                ax2.set_xlabel('Factor Length')
                ax2.set_ylabel('CCDF (P(Length ≥ x))')
                ax2.set_title('Reverse Complement Factors')
                ax2.grid(True, alpha=0.3)
                ax2.legend()

            plt.suptitle(f'Factor Length CCDF - {factors_filepath.stem}', fontsize=14, weight='bold')
            plt.tight_layout()

        else:
            # Plot on same axes
            fig, ax = plt.subplots(figsize=(10, 8))

            # Forward factors
            if forward_lengths:
                lengths_fwd, ccdf_fwd = compute_ccdf(forward_lengths)
                ax.loglog(lengths_fwd, ccdf_fwd, 'o-', color='blue', alpha=0.7,
                         label=f'Forward (n={len(forward_lengths)})')

            # Reverse complement factors
            if reverse_lengths:
                lengths_rev, ccdf_rev = compute_ccdf(reverse_lengths)
                ax.loglog(lengths_rev, ccdf_rev, 'o-', color='red', alpha=0.7,
                         label=f'Reverse Complement (n={len(reverse_lengths)})')

            ax.set_xlabel('Factor Length')
            ax.set_ylabel('CCDF (P(Length ≥ x))')
            ax.set_title(f'Factor Length CCDF - {factors_filepath.stem}', fontsize=14, weight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend()

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

    except Exception as e:
        raise PlotError(f"Failed to create factor length CCDF plot: {e}")


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
    self_factors_parser = subparsers.add_parser('self-factors-plot', help='Plot self-factors (simple matplotlib by default, use --interactive for Panel/Datashader)')
    self_factors_parser.add_argument('--fasta_filepath', help='Path to FASTA file')
    self_factors_parser.add_argument('--factors_filepath', help='Path to binary factors file')
    self_factors_parser.add_argument('--name', default=None, help='Name for the plot title')
    self_factors_parser.add_argument('--save_path', default=None, help='Path to save the plot image')
    self_factors_parser.add_argument('--no-show', action='store_true', help='Do not display the plot')
    self_factors_parser.add_argument('--interactive', action='store_true', help='Use interactive Panel/Datashader plot instead of simple matplotlib')
    self_factors_parser.add_argument('--return_panel', action='store_true', help='Whether to return the Panel app (interactive mode only)')

    # Subparser for factor length CCDF
    ccdf_parser = subparsers.add_parser('factor-length-ccdf', help='Plot empirical CCDF of factor lengths on log-log axes')
    ccdf_parser.add_argument('factors_filepath', help='Path to binary factors file')
    ccdf_parser.add_argument('--save_path', default=None, help='Path to save the plot image')
    ccdf_parser.add_argument('--no-show', action='store_true', help='Do not display the plot')
    ccdf_parser.add_argument('--no-separate', action='store_true', help='Plot forward and reverse complement factors on the same axes')

    # Subparser for reference sequence plotting
    ref_plot_parser = subparsers.add_parser('reference-plot', help='Plot target sequence factorized with reference sequence')
    ref_plot_parser.add_argument('reference_seq', nargs='?', help='Reference sequence (DNA by default, optional if factors_filepath provided)')
    ref_plot_parser.add_argument('target_seq', nargs='?', help='Target sequence (DNA by default, optional if factors_filepath provided)')
    ref_plot_parser.add_argument('--factors_filepath', help='Path to binary factors file (optional, will compute if not provided)')
    ref_plot_parser.add_argument('--reference_name', default='Reference', help='Name for the reference sequence')
    ref_plot_parser.add_argument('--target_name', default='Target', help='Name for the target sequence')
    ref_plot_parser.add_argument('--save_path', default=None, help='Path to save the plot image')
    ref_plot_parser.add_argument('--show_plot', action='store_true', default=True, help='Whether to display the plot')
    ref_plot_parser.add_argument('--interactive', action='store_true', help='Use interactive Panel/Datashader plot instead of simple matplotlib')
    ref_plot_parser.add_argument('--return_panel', action='store_true', help='Whether to return the Panel app (interactive mode only)')
    ref_plot_parser.add_argument(
        '--factorization-mode',
        choices=['dna', 'general'],
        default='dna',
        help='Choose "dna" for reverse-complement-aware factorization or "general" for arbitrary strings'
    )

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
        # Choose between simple and interactive plotting
        if args.interactive:
            plot_multiple_seq_self_lz_factor_plot_from_file(
                fasta_filepath=args.fasta_filepath,
                factors_filepath=args.factors_filepath,
                name=args.name,
                save_path=args.save_path,
                show_plot=not args.no_show,
                return_panel=args.return_panel
            )
        else:
            plot_multiple_seq_self_lz_factor_plot_simple(
                fasta_filepath=args.fasta_filepath,
                factors_filepath=args.factors_filepath,
                name=args.name,
                save_path=args.save_path,
                show_plot=not args.no_show
            )
    elif args.command == 'factor-length-ccdf':
        plot_factor_length_ccdf(
            factors_filepath=args.factors_filepath,
            save_path=args.save_path,
            show_plot=not args.no_show,
            separate=not args.no_separate
        )
    elif args.command == 'reference-plot':
        # Validate that sequences are provided if factors_filepath is not
        if args.factors_filepath is None and (args.reference_seq is None or args.target_seq is None):
            parser.error("reference_seq and target_seq are required when factors_filepath is not provided")
        
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
                return_panel=args.return_panel,
                factorization_mode=args.factorization_mode
            )

        else:
            plot_reference_seq_lz_factor_plot_simple(
                reference_seq=args.reference_seq,
                target_seq=args.target_seq,
                factors_filepath=args.factors_filepath,
                reference_name=args.reference_name,
                target_name=args.target_name,
                save_path=args.save_path,
                show_plot=args.show_plot,
                factorization_mode=args.factorization_mode
            )


