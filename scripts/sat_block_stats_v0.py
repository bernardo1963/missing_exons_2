#!/usr/bin/env python3
# sat_block_stats_v0.py  created with help from  DeepSeek  16ago2025   v.  16set2025
# Used to perform the analyses of the ms. "Triplex DNA and inverted repeats cause long-read sequencing bias against satellite DNA", by AB Carvalho, B Kim, F Uno

# sat_block_stats_v0.py SRR29479668_mms50_mtc5_blocks.txt  --target size --histo --mon_bp 5 
# sat_block_stats_v0.py SRR29479668_mms50_mtc5_blocks.txt  --target size --histo --mon_bp 5 --orientation F   # only gets the matches in F (ie, will get AAC, but not GTT)

# Default PNG output
# sat_block_stats_v0.py input.txt --target copies --histo

# SVG vector output
# sat_block_stats_v0.py input.txt --target size --histo --graph_format svg

# PDF output with monomer size filter
# sat_block_stats_v0.py input.txt --target size --histo --mon_bp 12 --graph_format pdf

import pandas as pd
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from matplotlib import rcParams

def calculate_statistics(df, stat_column):
    """Calculate statistics for the specified column grouped by monomer."""
    grouped = df.groupby('monomer')[stat_column]
    
    stats = grouped.agg([
        ('N', 'count'),
        ('sum', 'sum'),
        ('mean', 'mean'),
        ('median', 'median'),
        ('95percentile', lambda x: np.percentile(x, 95)),
        ('max', 'max')
    ]).reset_index()
    
    mon_bp = df.groupby('monomer')['mon_bp'].first().reset_index()
    stats = pd.merge(stats, mon_bp, on='monomer')
    stats = stats[['monomer', 'mon_bp', 'N', 'sum', 'mean', 'median', '95percentile', 'max']]
    return stats

def plot_histograms(df, stat_column, mon_bp_filter=None, bins='auto', density=False, xlog=False, ylog=False):
    """Plot histograms for each monomer's target variable distribution."""
    if mon_bp_filter is not None:
        df = df[df['mon_bp'] == mon_bp_filter]
        if df.empty:
            print(f"No data found for monomers with size {mon_bp_filter}", file=sys.stderr)
            return None
    
    monomers = df['monomer'].unique()
    if not len(monomers):
        print("No data available for plotting", file=sys.stderr)
        return None
    
    rcParams.update({
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 9
    })
    
    cols = min(3, len(monomers))
    rows = (len(monomers) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(15, 5*rows))
    
    if len(monomers) == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for i, (monomer, ax) in enumerate(zip(monomers, axes)):
        data = df[df['monomer'] == monomer][stat_column]
        if len(data) < 2:
            ax.text(0.5, 0.5, f"Not enough data\nfor {monomer}", 
                    ha='center', va='center')
            ax.set_title(monomer)
            continue
        
        # Filter out zeros if using log scale
        plot_data = data
        if xlog:
            plot_data = plot_data[plot_data > 0]
            if len(plot_data) == 0:
                ax.text(0.5, 0.5, f"No positive values\nfor {monomer}", 
                        ha='center', va='center')
                ax.set_title(monomer)
                continue
        
        counts, bins, patches = ax.hist(
            plot_data, 
            bins=bins,
            histtype='step', 
            linewidth=1.5, 
            edgecolor='black',
            fill=False,
            density=density,
            log=ylog
        )
        
        if xlog:
            ax.set_xscale('log')
        
        median = np.median(plot_data)
        mean = np.mean(plot_data)
        percentile95 = np.percentile(plot_data, 95)
        
        ax.axvline(median, color='red', linestyle='--', label=f'Median: {median:.1f}')
        ax.axvline(mean, color='green', linestyle='-.', label=f'Mean: {mean:.1f}')
        ax.axvline(percentile95, color='purple', linestyle=':', label=f'95%: {percentile95:.1f}')
        
        ax.set_title(f"{monomer} (bp={df[df['monomer']==monomer]['mon_bp'].iloc[0]})")
        ax.set_xlabel(stat_column)
        ax.set_ylabel('Density' if density else 'Count')
        ax.legend()
    
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()
    return fig

def main():
    parser = argparse.ArgumentParser(
        description='Calculate statistics for satellite block distributions.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('input_file', help='Path to the input tabular data file')
    parser.add_argument('--target', required=True, choices=['copies', 'size'], 
                      help='Column to calculate statistics for')
    parser.add_argument('--output', help='Optional output file to save results')
    parser.add_argument('--histo', action='store_true', 
                      help='Generate histograms of the distribution')
    parser.add_argument('--mon_bp', type=int, 
                      help='Filter results to monomers of specific size (bp)')
    parser.add_argument('--graph_format', default='png', 
                      choices=['png', 'svg', 'pdf', 'jpg', 'tiff', 'eps'],
                      help='Output format for the histograms')
    parser.add_argument('--bins', type=str, default='auto',
                      help='Number of bins for histograms ("auto" by default)')
    parser.add_argument('--density', action='store_true',
                      help='Show density instead of counts in histograms')
    parser.add_argument('--Xlog', action='store_true',
                      help='Use logarithmic scale for X-axis')
    parser.add_argument('--Ylog', action='store_true',
                      help='Use logarithmic scale for Y-axis')
    parser.add_argument('--orientation', default='all', choices=['F', 'R', 'FR', 'all'], 
                      help='Column to calculate statistics for')    
    args = parser.parse_args()
    
    try:
        # df = pd.read_csv(args.input_file, sep='\s+')
        df = pd.read_csv(args.input_file, sep='\s+')
        if args.orientation != 'all':
            df = df[df['orient'] == args.orientation]
  
        if args.target not in df.columns:
            raise ValueError(f"Column '{args.target}' not found in input file")
        
        if args.mon_bp is not None:
            df = df[df['mon_bp'] == args.mon_bp]
            if df.empty:
                raise ValueError(f"No data found for monomers with size {args.mon_bp}")
        
        stats = calculate_statistics(df, args.target)

        # First, let's define the base column names
        base_columns = ['monomer', 'mon_bp', 'N', 'sum', 'mean', 'median', '95percentile', 'max']
        # Conditionally modify the column names based on orientation
        if args.orientation == 'all':
            output_columns = base_columns
        else:
            # Add suffix to all columns except the first two
            output_columns = base_columns[:2] + [f"{col}_{args.orientation}" for col in base_columns[2:]]


        if args.output:
            stats.to_csv(args.output, sep='\t', index=False, 
                        float_format='%.1f',
                        # columns=['monomer', 'mon_bp', 'N', 'sum', 'mean', 'median', '95percentile', 'max'])
                         columns=base_columns,  # Use the original column names to select data
                         header=output_columns)  # Use the modified names for the header                        
            print(f"Results saved to {args.output}", file=sys.stdout)
        else:
            print(
                stats.to_string(
                    index=False, 
                    float_format='%.1f',
                    columns=base_columns,
                    header=output_columns
                ),
                file=sys.stdout
            )

 
        if args.histo:
            bins = args.bins if not args.bins.isdigit() else int(args.bins)
            fig = plot_histograms(
                df, 
                args.target, 
                args.mon_bp,
                bins=bins,
                density=args.density,
                xlog=args.Xlog,
                ylog=args.Ylog
            )
            if fig:
                plot_filename = f"{args.input_file.rsplit('.', 1)[0]}_{args.target}_histograms.{args.graph_format}"
                fig.savefig(plot_filename, format=args.graph_format, dpi=300)
                print(f"Histograms saved to {plot_filename}", file=sys.stderr)
                
                if not args.output and 'DISPLAY' in os.environ:
                    try:
                        import matplotlib.pyplot as plt_backend
                        plt_backend.switch_backend('TkAgg')
                        plt_backend.show()
                    except Exception as e:
                        print(f"Could not display plot: {e}", file=sys.stderr)
                
                plt.close(fig)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()