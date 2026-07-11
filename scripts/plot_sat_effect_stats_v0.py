#!/usr/bin/env python3.9
# plot_sat_effect_stats_v0.py  written by Bernardo w/ DeepSeek help.   5oct2025    v. 28jun2026
# Used to perform the analyses of the ms. "Triplex DNA and inverted repeats cause long-read sequencing bias against satellite DNA", by AB Carvalho, B Kim, F Uno
# to analyse the effect of the satellites in LILAP assembly and HiFi read: how much do they enter inside sat blocks?] 
# usage:plot_sat_effect_stats_v0.py --input your_data.txt --variable HiFi_reads --stats median --graph_min -1000 --graph_max 5000 --monomer "AATATAT AAGAG"

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import math
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']  # Specify Arial as preferred sans-serif font
plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.grid'] = False # suppress background grid


def round_to_next_thousand(x, direction='up'):
    """
    Round to the next 1000, with direction control.
    For positive numbers: round up to next 1000
    For negative numbers: round down to next -1000
    """
    if x == 0:
        return 0
    
    if x > 0:
        if direction == 'up':
            return math.ceil(x / 1000) * 1000
        else:
            return math.floor(x / 1000) * 1000
    else:
        if direction == 'up':
            return math.ceil(x / 1000) * 1000
        else:
            return math.floor(x / 1000) * 1000

def format_table(stats_df):
    """
    Format the statistics table with aligned columns for human readability
    """
    # Calculate column widths
    monomer_width = max(stats_df['monomer'].str.len().max(), len('monomer'))
    n_width = max(stats_df['n'].astype(str).str.len().max(), len('n'))
    
    # For float columns, find the maximum width needed
    mean_vals = stats_df['mean'].apply(lambda x: f"{x:.2f}")
    median_vals = stats_df['median'].apply(lambda x: f"{x:.2f}")
    # For min and max, use integer formatting
    min_vals = stats_df['min'].apply(lambda x: f"{int(x)}")
    max_vals = stats_df['max'].apply(lambda x: f"{int(x)}")
    
    mean_width = max(mean_vals.str.len().max(), len('mean'))
    median_width = max(median_vals.str.len().max(), len('median'))
    min_width = max(min_vals.str.len().max(), len('min'))
    max_width = max(max_vals.str.len().max(), len('max'))
    
    # Create header
    header = (f"{'monomer':<{monomer_width}} "
              f"{'n':>{n_width}} "
              f"{'mean':>{mean_width}} "
              f"{'median':>{median_width}} "
              f"{'min':>{min_width}} "
              f"{'max':>{max_width}}")
    
    # Create separator line
    separator = (f"{'-' * monomer_width} "
                 f"{'-' * n_width} "
                 f"{'-' * mean_width} "
                 f"{'-' * median_width} "
                 f"{'-' * min_width} "
                 f"{'-' * max_width}")
    
    # Create rows
    rows = []
    for _, row in stats_df.iterrows():
        row_str = (f"{row['monomer']:<{monomer_width}} "
                   f"{row['n']:>{n_width}} "
                   f"{row['mean']:>{mean_width}.2f} "
                   f"{row['median']:>{median_width}.2f} "
                   f"{int(row['min']):>{min_width}} "
                   f"{int(row['max']):>{max_width}}")
        rows.append(row_str)
    
    return header, separator, rows

def main():
    parser = argparse.ArgumentParser(description='Plot statistics for monomers')
    parser.add_argument('--input', required=True, help='Input dataset file')
    parser.add_argument('--variable', required=True, 
                       help='Column name to calculate statistics on')
    parser.add_argument('--stats', required=True, 
                       choices=['mean', 'median', 'min', 'max'],
                       help='Statistic to plot')
    parser.add_argument('--graph_min', type=float,
                       help='Starting point for horizontal lines. If not provided, will be computed from data.')
    parser.add_argument('--graph_max', type=float,
                       help='Maximum x-axis value. If not provided, will be computed from data.')
    parser.add_argument('--monomer', type=str,
                       help='Space-separated list of monomers to plot')
    parser.add_argument('--order', type=str,
                       help='Space-separated list of monomers in desired order')
    parser.add_argument('--graph', choices=['none', 'svg', 'png', 'pdf'], 
                       default='png', help='Output format for graph (default: png)')
    parser.add_argument('--graph_suffix', type=str, default='',
                       help='Suffix to add to graph filename to avoid overwriting')
    parser.add_argument('--line_weight', type=float, default=2,
                       help='Line weight for horizontal lines (default: 2)')
    parser.add_argument('--line_color', type=str, default='blue',
                       help='Line color for horizontal lines (default: blue)')
    parser.add_argument('--Xlabel', type=str,
                       help='Custom x-axis label (default: based on variable and stats)')
    parser.add_argument('--font_size', type=float, default=10,
                       help='Font size for graph text (default: 10)')
    parser.add_argument("--figsize_W", type=float, default=10, help="fig size WIDTH  (in inches)")
    parser.add_argument("--figsize_H", type=float, default=6 , help="fig size HEIGHT (in inches)")
                       
    args = parser.parse_args()
    
    if args.monomer == None:       # new in 4may2026
        args.monomer = args.order
    
    
    
    # Read the data with more robust error handling
    try:
        # First try with whitespace separator (handles multiple spaces)
        df = pd.read_csv(args.input, sep='\s+', engine='python', na_values=['.', 'NA', 'N/A', ''])
        print(f"Successfully read file with whitespace separator")
    except Exception as e:
        print(f"Error reading file with whitespace separator: {e}")
        print("Trying with tab separator...")
        try:
            df = pd.read_csv(args.input, sep='\t', na_values=['.', 'NA', 'N/A', ''])
            print(f"Successfully read file with tab separator")
        except Exception as e2:
            print(f"Error reading file with tab separator: {e2}")
            sys.exit(1)
    
    # Check if required columns exist
    required_columns = ['monomer', args.variable]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: Missing columns in dataset: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)
    
    # Filter data by monomer if specified
    if args.monomer:
        monomers_to_plot = args.monomer.split()
        df = df[df['monomer'].isin(monomers_to_plot)]
    
    if df.empty:
        print("Error: No data found for the specified monomers")
        sys.exit(1)
    
    # Convert variable column to numeric, handling non-numeric values and missing values
    # First make a copy of the original for reporting
    original_count = len(df)
    df[args.variable] = pd.to_numeric(df[args.variable], errors='coerce')
    
    # Count how many values were converted to NaN
    missing_count = df[args.variable].isna().sum()
    if missing_count > 0:
        print(f"Warning: {missing_count} out of {original_count} values in '{args.variable}' column are missing or non-numeric and will be excluded")
    
    df = df.dropna(subset=[args.variable])
    
    if df.empty:
        print("Error: No valid numeric data found after removing missing values")
        sys.exit(1)
    
    # Calculate statistics for each monomer
    stats_df = df.groupby('monomer')[args.variable].agg([
        ('n', 'count'),
        ('mean', 'mean'),
        ('median', 'median'),
        ('min', 'min'),
        ('max', 'max')
    ]).reset_index()
    
    # Apply monomer ordering if specified
    if args.order:
        order_list = args.order.split()
        # Filter to only include monomers that exist in the data
        order_list = [m for m in order_list if m in stats_df['monomer'].values]
        # order_list.reverse() 
        # Add any missing monomers that are in data but not in order list
        missing_monomers = [m for m in stats_df['monomer'].values if m not in order_list]
        order_list.extend(missing_monomers)
        
        # Create a categorical type with the specified order
        stats_df['monomer'] = pd.Categorical(stats_df['monomer'], categories=order_list, ordered=True)
        stats_df = stats_df.sort_values('monomer')
    
    # Print the formatted statistics table to stdout
    header, separator, rows = format_table(stats_df)
    print(header)
    print(separator)
    for row in rows:
        print(row)
    
    # Get the statistic values for plotting
    stat_values = stats_df[args.stats].values
    monomers = stats_df['monomer'].values
    counts = stats_df['n'].values
    
    # Calculate graph_min and graph_max if not provided
    if args.graph_max is None:
        max_val = max(stat_values)
        args.graph_max = round_to_next_thousand(max_val, 'up')
        print(f"\nAuto-calculated graph_max: {args.graph_max}")
    
    if args.graph_min is None:
        min_val = min(stat_values)
        args.graph_min = round_to_next_thousand(min_val, 'down')
        print(f"Auto-calculated graph_min: {args.graph_min}")
    
    # Create the plot
    plt.figure(figsize=(args.figsize_W, args.figsize_H))
    
    # Set font size for all text elements
    plt.rc('font', size=args.font_size)
    plt.rc('axes', titlesize=args.font_size)
    plt.rc('axes', labelsize=args.font_size)
    plt.rc('xtick', labelsize=args.font_size)
    plt.rc('ytick', labelsize=args.font_size)
    plt.rc('legend', fontsize=args.font_size)
    
    # Create horizontal lines
    y_positions = range(len(monomers))
    
    for i, (monomer, stat_val) in enumerate(zip(monomers, stat_values)):
        plt.hlines(y=i, xmin=args.graph_min, xmax=stat_val, 
                  linewidth=args.line_weight, color=args.line_color, alpha=0.7)
    
    # Add monomer labels with count on the left
    monomer_labels = [f"{monomer} ({count})" for monomer, count in zip(monomers, counts)]
    plt.yticks(y_positions, monomer_labels)
    
    # Add vertical dashed line at x=0
    plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='x=0')
    
    # Set plot limits and labels
    plt.xlim(args.graph_min, args.graph_max)
    plt.ylim(-0.5, len(monomers) - 0.5)
    
    # Set x-axis label
    if args.Xlabel:
        plt.xlabel(args.Xlabel, fontsize=args.font_size)
    else:
        plt.xlabel(f'{args.stats.capitalize()} of {args.variable}', fontsize=args.font_size)
    
    plt.ylabel('Monomer', fontsize=args.font_size)
    # plt.grid(True, alpha=0.3)
    plt.grid(False)  # suppress backfround grid  
    # Adjust layout to ensure labels are visible
    plt.tight_layout()
    
    # Save or show the graph
    if args.graph != 'none':
        # Generate output filename based on input filename
        # base_name = os.path.splitext(args.input)[0]
        base_name = os.path.splitext(os.path.basename(args.input))[0] # new in 4may2026
        suffix = f"_{args.graph_suffix}" if args.graph_suffix else ""
        output_file = f"{base_name}{suffix}.{args.graph}"
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nGraph saved as {output_file}")
    
    # Only try to show if we're likely in an interactive environment
    if args.graph == 'none':
        # Check if we're likely in an interactive environment
        if 'DISPLAY' in os.environ:
            try:
                plt.show()
            except Exception as e:
                print(f"Note: Could not display interactive plot: {e}")
                print("Plot was generated but cannot be shown.")
        else:
            print("No display available. Plot was generated but cannot be shown.")
            print("Use --graph to save the plot in a file format (svg, png, pdf).")

if __name__ == '__main__':
    main()