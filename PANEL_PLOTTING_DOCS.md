# Interactive LZ Factor Plotting with Datashader/Panel

This module provides a high-performance, interactive visualization for non-overlapping Lempel-Ziv-Storer-Szymanski (LZSS) factorization results from DNA FASTA files.

## Installation

Install the base package:
```bash
pip install noLZSS
```

For interactive plotting functionality, install with panel dependencies:
```bash
pip install 'noLZSS[panel]'
```

Or install dependencies manually:
```bash
pip install numpy pandas holoviews bokeh panel datashader colorcet
```

## Usage

### Python API

```python
from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta

# Basic usage - creates interactive plot
app = plot_multiple_seq_self_lz_factor_plot_from_fasta(
    fasta_filepath="sequences.fasta",
    show_plot=True,
    return_panel=True
)

# For Jupyter notebooks, display directly
app

# For scripts, serve the app
import panel as pn
pn.serve(app, show=True)
```

### Command Line Interface

```bash
# Interactive plot from FASTA file
python lz_panel.py --fasta sequences.fasta

# Demo mode with synthetic data
python lz_panel.py

# Panel server mode
panel serve lz_panel.py --show

# Custom port
python lz_panel.py --fasta sequences.fasta --port 8080
```

## Features

### Visual Semantics

- **Blue lines**: Forward factors (y0=ref, y1=ref+length)
- **Red lines**: Reverse complement factors (y0=ref+length, y1=ref)  
- **Gray dashed line**: y=x diagonal reference
- **Equal aspect ratio**: Synchronized x/y ranges for proper visualization

### Interactive Controls

- **Length Filter**: IntRangeSlider to filter factors by length
- **Show Hover Overlay**: Toggle decimated overlay with tooltips
- **Top-k per Pixel**: Select 1-5 factors per screen pixel for hover
- **Colormap**: Choose intensity colormap for datashaded rendering
- **Export PNG**: Export current view (requires selenium/chromedriver)

### Performance Features

- **Datashader Rendering**: Handles millions of factors smoothly
- **Level-of-Detail**: Zoom/pan-aware decimation for overlay
- **Screen-Space Binning**: Max 2000 bins for optimal performance
- **Vectorized Operations**: NumPy/Pandas for fast filtering

### Hover Information

When overlay is enabled, hover reveals:
- **Start**: Factor start position
- **Length**: Factor length  
- **End**: Factor end position
- **Reference**: Reference position
- **Direction**: "forward" or "reverse-complement"

## Function Signatures

### Main Function

```python
def plot_multiple_seq_self_lz_factor_plot_from_fasta(
    fasta_filepath: Union[str, Path],
    name: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    show_plot: bool = True,
    return_panel: bool = False
) -> Optional[panel.viewable.Viewable]:
```

**Parameters:**
- `fasta_filepath`: Path to FASTA file containing DNA sequences
- `name`: Optional plot title (defaults to filename)
- `save_path`: Optional PNG export path
- `show_plot`: Whether to display/serve the plot
- `return_panel`: Whether to return Panel app for embedding

**Returns:**
- Panel app if `return_panel=True`, otherwise `None`

### Backward Compatibility

The original Matplotlib function remains available:

```python
from noLZSS.genomics.plots import plot_multiple_seq_self_weizmann_factor_plot_from_fasta

plot_multiple_seq_self_weizmann_factor_plot_from_fasta(
    fasta_filepath="sequences.fasta",
    show_plot=True
)
```

## Error Handling

The function provides comprehensive error handling:

- **ImportError**: Missing dependencies with install instructions
- **FileNotFoundError**: FASTA file not found
- **PlotError**: No factors found or processing failure

Example:
```python
try:
    app = plot_multiple_seq_self_lz_factor_plot_from_fasta("file.fasta")
except ImportError as e:
    print(f"Install dependencies: {e}")
except FileNotFoundError as e:
    print(f"File issue: {e}")
except PlotError as e:
    print(f"Processing error: {e}")
```

## Performance Guidelines

### Recommended Usage

- **Large datasets**: Use datashader base layers for overview
- **Interactive exploration**: Enable hover overlay for detailed inspection
- **Memory efficiency**: Use file-based processing for huge FASTA files
- **Export**: PNG export requires additional selenium setup

### Performance Targets

- **Base layer**: ≥5M segments rendered smoothly via datashader
- **Overlay**: ≤100ms recompute for ≤2000 bins on zoom/pan
- **Memory**: Vectorized operations minimize memory overhead

## Examples

### Jupyter Notebook

```python
import panel as pn
from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta

# Enable Panel in Jupyter
pn.extension()

# Create and display interactive plot
app = plot_multiple_seq_self_lz_factor_plot_from_fasta(
    fasta_filepath="example.fasta",
    name="Example Genome",
    return_panel=True
)

app  # Display in notebook cell
```

### Standalone Script

```python
import panel as pn
from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta

def main():
    app = plot_multiple_seq_self_lz_factor_plot_from_fasta(
        fasta_filepath="large_genome.fasta",
        show_plot=False,
        return_panel=True
    )
    
    print("Starting server at http://localhost:5007")
    pn.serve(app, show=True, port=5007)

if __name__ == "__main__":
    main()
```

### Batch Processing

```python
import os
from pathlib import Path
from noLZSS.genomics.plots import plot_multiple_seq_self_lz_factor_plot_from_fasta

# Process multiple FASTA files
fasta_dir = Path("genomes/")
output_dir = Path("plots/")

for fasta_file in fasta_dir.glob("*.fasta"):
    try:
        # Save PNG export for each file
        plot_multiple_seq_self_lz_factor_plot_from_fasta(
            fasta_filepath=fasta_file,
            save_path=output_dir / f"{fasta_file.stem}_plot.png",
            show_plot=False
        )
        print(f"Processed: {fasta_file.name}")
    except Exception as e:
        print(f"Failed to process {fasta_file.name}: {e}")
```

## Technical Details

### Data Flow

1. **Input**: FASTA file with multiple DNA sequences
2. **Factorization**: `factorize_fasta_multiple_dna_w_rc()` C++ function
3. **DataFrame**: Convert factors to pandas DataFrame with plot coordinates
4. **Datashader**: Render base layers with max('length') aggregator
5. **Overlay**: Decimated overlay with hover tooltips via HoloViews streams
6. **Controls**: Panel widgets for interactive filtering and configuration

### Coordinate System

Factor coordinates are calculated as:
- **x0**: start position
- **x1**: start + length
- **Forward factors**: y0=ref, y1=ref+length (blue upward lines)
- **Reverse complement**: y0=ref+length, y1=ref (red downward lines)

### Dependencies

- **Core**: numpy, pandas
- **Visualization**: holoviews, bokeh, datashader
- **Interface**: panel, colorcet
- **Optional**: selenium (for PNG export)

## Troubleshooting

### Common Issues

1. **ImportError**: Install `pip install 'noLZSS[panel]'`
2. **Empty plot**: Check FASTA file format and content
3. **Performance issues**: Reduce overlay k-value or disable overlay for huge datasets
4. **Export fails**: Install selenium and chromedriver for PNG export

### Debug Mode

Enable verbose output:
```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Function will print progress messages
app = plot_multiple_seq_self_lz_factor_plot_from_fasta("file.fasta")
```