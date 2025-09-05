#!/usr/bin/env python3
"""
Phylogenetic Tree Visualization Script using baltic_bokeh

This script creates interactive phylogenetic tree visualizations from Newick files
and metadata, using the baltic_bokeh package for web-based interactive plots.

Usage:
    python bokeh-phy.py --tree TREE_FILE --metadata METADATA_FILE [OPTIONS]

Example:
    python bokeh-phy.py --tree results/vgsc_focal.fasta.treefile --metadata fastas/vgsc_focal.metadata.tsv --output tree_plot.html
"""

import argparse
import sys
import os
from pathlib import Path

import baltic as bt
import baltic_bokeh as bt_bokeh
import pandas as pd
import numpy as np
import plotly.express as px


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Create interactive phylogenetic tree visualizations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --tree tree.newick --metadata samples.tsv
  %(prog)s --tree results/vgsc_focal.fasta.treefile --metadata fastas/vgsc_focal.metadata.tsv --output my_tree.html
  %(prog)s --tree tree.newick --metadata samples.tsv --color-column species --width 1200 --height 800
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--tree', '-t',
        required=True,
        help='Path to Newick tree file (e.g., .treefile from IQ-TREE)'
    )
    
    parser.add_argument(
        '--metadata', '-m',
        required=True,
        help='Path to metadata TSV file with sample information'
    )
    
    # Optional arguments
    parser.add_argument(
        '--output', '-o',
        default='phylo_tree.html',
        help='Output HTML file name (default: phylo_tree.html)'
    )
    
    parser.add_argument(
        '--color-column', '-c',
        default='taxon',
        help='Metadata column to use for tip coloring (default: taxon)'
    )
    
    parser.add_argument(
        '--hover-columns',
        nargs='+',
        default=['taxon', 'country'],
        help='Metadata columns to show in hover tooltips (default: taxon country)'
    )
    
    parser.add_argument(
        '--width',
        type=int,
        default=1000,
        help='Plot width in pixels (default: 1000)'
    )
    
    parser.add_argument(
        '--height',
        type=int,
        default=1000,
        help='Plot height in pixels (default: 1000)'
    )
    
    parser.add_argument(
        '--tree-type',
        choices=['c', 'r', 'u'],
        default='c',
        help='Tree layout type: c=circular, r=rectangular, u=unrooted (default: c)'
    )
    
    parser.add_argument(
        '--marker',
        choices=['circle', 'hex', 'square', 'triangle'],
        default='hex',
        help='Tip marker shape (default: hex)'
    )
    
    parser.add_argument(
        '--marker-size',
        type=int,
        default=15,
        help='Tip marker size (default: 15)'
    )
    
    parser.add_argument(
        '--max-branch-length',
        type=float,
        default=0.005,
        help='Maximum branch length for display (default: 0.005)'
    )
    
    parser.add_argument(
        '--min-branch-length',
        type=float,
        default=5e-5,
        help='Minimum branch length for display (default: 5e-5)'
    )
    
    parser.add_argument(
        '--backend',
        choices=['webgl', 'canvas', 'svg'],
        default='webgl',
        help='Rendering backend (default: webgl)'
    )
    
    parser.add_argument(
        '--custom-colors',
        help='Path to JSON file with custom color mapping (optional)'
    )
    
    parser.add_argument(
        '--show-plot',
        action='store_true',
        help='Show plot in browser after creating (default: False)'
    )
    
    return parser.parse_args()


def get_default_colors():
    """Get default color palette for common taxa."""
    TAXON_PALETTE = px.colors.qualitative.Vivid
    return {
        "gambiae": TAXON_PALETTE[1],
        "coluzzii": TAXON_PALETTE[0], 
        "arabiensis": TAXON_PALETTE[2],
        "merus": TAXON_PALETTE[3],
        "melas": TAXON_PALETTE[4],
        "quadriannulatus": TAXON_PALETTE[5],
        "fontenillei": TAXON_PALETTE[6],
        "gcx1": TAXON_PALETTE[7],
        "gcx2": TAXON_PALETTE[8],
        "gcx3": TAXON_PALETTE[9],
        "gcx4": TAXON_PALETTE[10],
        "bissau": TAXON_PALETTE[7],
        "pwani": TAXON_PALETTE[9],
        "unassigned": "black",
    }


def load_custom_colors(color_file):
    """Load custom color mapping from JSON file."""
    import json
    try:
        with open(color_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Warning: Could not load custom colors from {color_file}: {e}")
        print("Using default colors instead.")
        return get_default_colors()


def validate_inputs(tree_file, metadata_file):
    """Validate that input files exist and are readable."""
    if not os.path.exists(tree_file):
        raise FileNotFoundError(f"Tree file not found: {tree_file}")
    
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")
    
    # Check if tree file is readable
    try:
        with open(tree_file, 'r') as f:
            content = f.read().strip()
            if not content:
                raise ValueError(f"Tree file is empty: {tree_file}")
    except Exception as e:
        raise ValueError(f"Cannot read tree file {tree_file}: {e}")


def load_data(tree_file, metadata_file, color_column, hover_columns):
    """Load and validate tree and metadata."""
    print(f"Loading tree from: {tree_file}")
    try:
        tree = bt.loadNewick(tree_file)
    except Exception as e:
        raise ValueError(f"Error loading tree: {e}")
    
    print(f"Loading metadata from: {metadata_file}")
    try:
        # Try different separators
        for sep in ['\t', ',', ' ']:
            try:
                metadata = pd.read_csv(metadata_file, sep=sep)
                if len(metadata.columns) > 1:  # Successfully parsed multiple columns
                    break
            except:
                continue
        else:
            raise ValueError("Could not parse metadata file with any common separator")
            
    except Exception as e:
        raise ValueError(f"Error loading metadata: {e}")
    
    print(f"Loaded tree with {len(tree.getExternal())} tips")
    print(f"Loaded metadata for {len(metadata)} samples")
    
    # Validate required columns
    if color_column not in metadata.columns:
        available_cols = list(metadata.columns)
        raise ValueError(f"Color column '{color_column}' not found in metadata. Available columns: {available_cols}")
    
    # Check hover columns
    missing_hover = [col for col in hover_columns if col not in metadata.columns]
    if missing_hover:
        print(f"Warning: Hover columns not found in metadata: {missing_hover}")
        hover_columns = [col for col in hover_columns if col in metadata.columns]
    
    return tree, metadata, hover_columns


def create_plot(tree, metadata, args, color_map):
    """Create the phylogenetic tree plot."""
    print("Creating phylogenetic tree plot...")
    
    try:
        p = bt_bokeh.plotTree(
            tree,
            type=args.tree_type,
            df_metadata=metadata,
            color_column=args.color_column,
            color_discrete_map=color_map,
            size=args.marker_size,
            plot_width=args.width,
            plot_height=args.height,
            output_backend=args.backend,
            marker=args.marker,
            marker_line_color='white',
            marker_line_width=0,
            hover_data=args.hover_columns,
            max_branch_length=args.max_branch_length,
            min_branch_length=args.min_branch_length
        )
        return p
    except Exception as e:
        raise RuntimeError(f"Error creating plot: {e}")


def save_plot(plot, output_file, show_plot=False):
    """Save plot to HTML file and optionally show in browser."""
    from bokeh.plotting import output_file, save, show
    
    print(f"Saving plot to: {output_file}")
    
    # Ensure output directory exists
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    output_file(output_file)
    
    if show_plot:
        show(plot)
    else:
        save(plot)
    
    print(f"✅ Plot saved successfully to {output_file}")
    if not show_plot:
        print(f"   Open in browser: file://{os.path.abspath(output_file)}")


def main():
    """Main function."""
    args = parse_arguments()
    
    try:
        # Validate inputs
        validate_inputs(args.tree, args.metadata)
        
        # Load data
        tree, metadata, hover_columns = load_data(
            args.tree, args.metadata, args.color_column, hover_columns=args.hover_columns
        )
        
        # Update hover columns based on what's actually available
        args.hover_columns = hover_columns
        
        # Get color mapping
        if args.custom_colors:
            color_map = load_custom_colors(args.custom_colors)
        else:
            color_map = get_default_colors()
        
        # Create plot
        plot = create_plot(tree, metadata, args, color_map)
        
        # Save plot
        save_plot(plot, args.output, args.show_plot)
        
    except KeyboardInterrupt:
        print("\n❌ Interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()