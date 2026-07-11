#!/usr/bin/env python3.9

#############################################################################################################################################
# visualize_repeats_censor5.py  :  v.  7jul2026
# streamlined repeat classification. Sat DNA is now identificed by a regex; TEs are now the uncclassified stuff, after  aseries of elif
#
# written by Fabiana Uno
#
# The script processes a Censor output file (.map), classifies sequence regions (e.g., simple repeats, TEs, and Ycds), generates a
# graphical visualization, and outputs sequence composition statistics.
#############################################################################################################################################


import argparse
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import numpy as np
from scipy import ndimage
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']  # Specify Arial as preferred sans-serif font
plt.rcParams['pdf.fonttype'] = 42 
plt.rcParams['ps.fonttype'] = 42

def calculate_at_content(sequence):
    """Calculate AT% for a simple repeat sequence."""
    match = re.match(r'^\((.*?)\)n$', sequence)
    if not match:
        return None  # Not a simple repeat
    repeat_bases = match.group(1)
    at_count = repeat_bases.count('A') + repeat_bases.count('T')
    total_count = len(repeat_bases)
    return (at_count / total_count) * 100 if total_count > 0 else None


def parse_depth_file(depth_file):
    """Parse the depth/coverage file and return position and depth arrays."""
    positions = []
    depths = []
    
    if not depth_file or not os.path.exists(depth_file):
        return positions, depths
    
    print(f"Parsing depth file: {depth_file}")
    try:
        with open(depth_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        position = int(parts[0])
                        depth = int(parts[1])
                        # Filter positions within the specified range
                        if position >= args.Xmin and position <= args.Xmax:
                            positions.append(position)
                            depths.append(depth)
                    except ValueError as e:
                        print(f"Warning: Invalid data in depth file line {line_num}: {line} - {e}")
                else:
                    print(f"Warning: Insufficient columns in depth file line {line_num}: {line}")
    
    except Exception as e:
        print(f"Error parsing depth file {depth_file}: {e}")
        return [], []
    
    print(f"Successfully parsed {len(positions)} depth entries from {depth_file}")
    return positions, depths


def parse_cds_file(cds_file):
    """Parse the CDS/exon file and return a list of exon information."""
    exons = []
    if not cds_file or not os.path.exists(cds_file):
        return exons

    print(f"Parsing CDS/exon file: {cds_file}")
    try:
        with open(cds_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue

                parts = line.split()
                if len(parts) >= 3:
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        exon_id = parts[2]
                        if start >= args.Xmin and end <= args.Xmax:  # new in 6jun2025
                            exons.append({'start': start, 'end': end, 'id': exon_id})
                    except ValueError as e:
                        print(f"Warning: Invalid data in CDS file line {line_num}: {line} - {e}")

                elif len(parts) == 2:   # new in 10jun2025  allow unlabeled exons
                    try:
                        start = int(parts[0])
                        end = int(parts[1])
                        exon_id = " "
                        if start >= args.Xmin and end <= args.Xmax:  # new in 6jun2025
                            exons.append({'start': start, 'end': end, 'id': exon_id})
                    except ValueError as e:
                        print(f"Warning: Invalid data in CDS file line {line_num}: {line} - {e}")

                else:
                    print(f"Warning: Insufficient columns in CDS file line {line_num}: {line}")

    except Exception as e:
        print(f"Error parsing CDS file {cds_file}: {e}")
        return []

    print(f"Successfully parsed {len(exons)} exons from {cds_file}")
    return exons




def process_censor_output(censor_file, total_length):
    """Parses the Censor output file and classifies sequence regions."""
    classifications = []
    satellite_counts = {}  # Tracks satellite counts for summaries
    debug_log = []

    with open(censor_file, 'r') as file:
        for line in file:
            if line.strip():
                try:
                    columns = line.split()
                    start = int(columns[1])
                    end = int(columns[2])
                    repeat = columns[3]  # 4th column
                    # escaped_repeat = re.escape(repeat)  # new in 6set2025  the parenthisis in (AAGAG)n were interfering 
                    bases = end - start + 1                   

                    if re.match(fr"({args.exclude_regex})", repeat):   # skips undesired entries in censor, eg, anaother gene
                        continue
                    if  start < args.Xmin or end > args.Xmax:  # new in 6jun2025
                        continue

                    # if re.match(r'^Ycds_', repeat):  # Handle Ycds sequences
                    if re.match(fr"({args.gene_regex})", repeat):  # Handle gene names  w/ dynamic regex
                        classifications.append({
                            'type': 'gene',  # Rename for clarity
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"gene: {repeat} Start: {start} End: {end} Bases: {bases}")

                    elif bases < args.min_repeat_len :   # skips small sat or TE
                        continue

                    elif re.match(fr"({args.repeat_regex})", repeat):  # Handle a specific repeat name  w/ dynamic regex    e.g.,  1.688_GK_42F7/2-X
                    # elif re.match(fr"({args.repeat_regex})", escaped_repeat):  # Handle a specific repeat name  w/ dynamic regex    e.g.,  1.688_GK_42F7/2-X
                        classifications.append({
                            'type': 'repeat_regex',
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"repeat_regex: {repeat} Start: {start} End: {end} Bases: {bases}")
                    
                    elif args.repeat_exact == repeat:  # Handle a specific repeat name   e.g., (AAGAG)n
                        classifications.append({
                            'type': 'repeat_regex',
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"repeat_regex: {repeat} Start: {start} End: {end} Bases: {bases}")



                    elif repeat == "NNNN_block":  # Handle NNNN_block
                        classifications.append({
                            'type': 'NNNN_block',
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"NNNN_block: {repeat} Start: {start} End: {end} Bases: {bases}")

                    elif re.match("\([ATGC]+\)n", repeat):    #new in v3. simple repeat
                        at_content = calculate_at_content(repeat)
                        classifications.append({
                            'type': 'simple_repeat',
                            'start': start,
                            'end': end,
                            'at_content': at_content
                        })
                        if repeat in satellite_counts:
                            satellite_counts[repeat] += bases
                        else:
                            satellite_counts[repeat] = bases
                        debug_log.append(f"Simple Repeat: {repeat} Start: {start} End: {end} Bases: {bases}")

                    else:  # Transposable elements
                        # if re.match(r'^[A-Za-z0-9_-]+$', repeat):
                        classifications.append({
                            'type': 'TE',
                            'start': start,
                            'end': end
                        })
                        debug_log.append(f"TE: {repeat} Start: {start} End: {end} Bases: {bases}")

                    # else:  # Unclassified  will never reach here, will be cathed before , as a TE
                        # debug_log.append(f"Unclassified: {repeat} Start: {start} End: {end} Bases: {bases}")

                except (ValueError, IndexError) as e:
                    debug_log.append(f"Skipping malformed line: {line.strip()} Error: {e}")

    # Output debug log
    with open("debug_log.txt", "w") as debug_file:
        debug_file.write("\n".join(debug_log))

    # Prepare satellite summary
    satellite_summary = [
        (repeat, bases, (bases / total_length) * 100)
        for repeat, bases in satellite_counts.items()
    ]
    satellite_summary = sorted(satellite_summary, key=lambda x: x[2], reverse=True)
    print("\ninput file:", args.censor)
    # Print Top 10 satellite
    print("Top 10 Satellite DNA Summary:")
    print(f"{'Repeat':<15} {'Bases':<10} {'Percentage':<10}")
    for repeat, bases, percentage in satellite_summary[:10]:
        print(f"{repeat:<15} {bases:<10} {percentage:.1f}%")

    return classifications, satellite_summary


def summarize_categories(classifications, total_length, satellite_summary):
    """Summarize the sequence composition by category (single-copy, TE, satellite)."""
    global total_repeat_regex_bases
    total_sat_bases = 0
    total_te_bases = 0
    total_ycds_bases = 0
    total_repeat_regex_bases = 0
    total_NNNN_bases = 0
    covered_bases = 0  # Tracks bases covered by satellites, TEs, and gene

    for entry in classifications:
        start = entry['start']
        end = entry['end']
        bases = end - start + 1

        if entry['type'] == 'simple_repeat':
            total_sat_bases += bases
        elif entry['type'] == 'gene':
            total_ycds_bases += bases
        elif entry['type'] == 'repeat_regex':
            total_repeat_regex_bases += bases
        elif entry['type'] == 'TE':
            total_te_bases += bases
        elif entry['type'] == 'NNNN_block':
            total_NNNN_bases += bases
        covered_bases += bases

    # Calculate single-copy bases as the remainder
    single_copy_bases = total_length - covered_bases

    # Get the most common satellite DNA
    if satellite_summary:
        most_common_satellite, _, most_common_percentage = satellite_summary[0]
        details = f"{most_common_satellite} {most_common_percentage:.1f}%"
    else:
        details = "None"

    # Print the summary
    print("Sequence Composition Summary:")
    try:
        print(f"{'Category':<15} {'Bases':<10} {'Percentage':<15} {'Details':<20}")
        print(f"{'single-copy':<15} {single_copy_bases + total_ycds_bases:<10} {100 * (single_copy_bases + total_ycds_bases) / total_length:.1f}%")
        print(f"{'TE':<15} {total_te_bases:<10} {100 * total_te_bases / total_length:.1f}%")
        print(f"{'sat':<15} {total_sat_bases:<10} {100 * total_sat_bases / total_length:.1f}%           {details}")
        if args.repeat_regex != 'not_used':
            print(f"{'repeat_regex':<15} {total_repeat_regex_bases:<10} {100 * total_repeat_regex_bases / total_length:.1f}%")
        if args.repeat_exact != 'not_used':
            print(f"{'repeat_exact':<15} {total_repeat_regex_bases:<10} {100 * total_repeat_regex_bases / total_length:.1f}%")
    except ZeroDivisionError:
        print("Error: attempted division by zero in the Sequence Composition Summary")


def create_visualization(classifications, total_length, output_file, exons, cds_yaxis, depth_positions, depth_values):
    """Creates the graphical output using Matplotlib with depth plotting."""
    # Determine if we have depth data to adjust figure height
    has_depth = len(depth_positions) > 0
    fig_height = args.figsize_H_depth if has_depth else args.figsize_H_nodepth
    
    if has_depth:
        # Create subplot with separate axes for depth and annotations
        fig, (ax_depth, ax_anno) = plt.subplots(2, 1, figsize=(args.figsize_W, fig_height), height_ratios=[3, 1], sharex=True)
        if args.repeat_exact != 'not_used':
            if total_repeat_regex_bases < 1000:
                cleaned_sat =  args.repeat_exact.strip('()n') + "  " + str(total_repeat_regex_bases) + " bp"   # new in 31may2026
            else:
                cleaned_sat =  args.repeat_exact.strip('()n')  + "  " +  f"{total_repeat_regex_bases/1000:.1f}" + " kb"            
         
        if args.title != '' and args.title != 'sat_info' and args.title != '.' :
            ax_depth.set_title(args.title, fontsize=args.title_fontsize, color=args.title_color) 
        elif args.title == 'sat_info' and args.repeat_exact != 'not_used': 
            ax_depth.set_title(cleaned_sat, fontsize=args.title_fontsize,color=args.repeat_color)           
        if args.title != 'sat_info' and args.repeat_exact != 'not_used' and args.print_sat == 'yes' : # new in 21sep2025 / 27jun2026  args.print_sat
            ax_depth.text(0.97, 1.025, cleaned_sat, transform=ax_depth.transAxes,fontsize=9, ha='left', va='baseline', color=args.repeat_color)            
            # ax_depth.annotate(cleaned_sat,xy=(0.97, 1.025),xycoords='axes fraction',fontsize=9,ha='left',va='baseline',color=args.repeat_color,annotation_clip=False)



            
        # Set up the plot limits
        ax_depth.set_xlim(args.Xmin, total_length)
        ax_anno.set_xlim(args.Xmin, total_length)
        
        # Plot depth in top subplot (simple line, with optional smoothing)
        if args.Ymax==0:
            max_depth = max(depth_values) if depth_values else 0
            ax_depth.set_ylim(0, max_depth * 1.05)  # Small padding above depth
        else:
            max_depth = args.Ymax if depth_values else 0
            ax_depth.set_ylim(0, max_depth)
        if depth_positions and depth_values:
            # Apply smoothing if requested
            if args.smooth_depth > 0:
                # Convert to numpy arrays for smoothing
                depth_array = np.array(depth_values)
                smoothed_values = ndimage.gaussian_filter1d(depth_array, sigma=args.smooth_depth)
                ax_depth.plot(depth_positions, smoothed_values, color='blue', linewidth=args.linewidth, zorder=4)
                print(f"Applied Gaussian smoothing with sigma={args.smooth_depth}")
            else:
                ax_depth.plot(depth_positions, depth_values, color='blue', linewidth=args.linewidth, zorder=4)
        
        # Configure depth subplot (simple, clean)
        ax_depth.set_ylabel(args.Ylabel, fontsize=args.label_fontsize)
        ax_depth.tick_params(axis="y", labelsize=args.scale_fontsize)
        # Remove all ticks from depth subplot
        ax_depth.set_xticks([])
        ax_depth.tick_params(axis="x", which="both", length=0)  # Hide any remaining x-ticks
        
        # Set up annotation subplot (bottom) - increased separation
        ax_anno.set_ylim(0, 1)
        repeat_y_base = 0.6
        repeat_height = 0.2
        baseline_y = 0.7
        exon_y_base = 0.1
        exon_height = 0.12
        
        # Use annotation subplot for plotting elements
        ax = ax_anno
        
    else:
        # Original single subplot layout without depth
        fig, ax = plt.subplots(figsize=(args.figsize_W, fig_height))
        ax.set_xlim(args.Xmin, total_length)
        ax.set_ylim(0, 1)
        repeat_y_base = 0.375
        repeat_height = 0.25
        baseline_y = 0.5
        exon_y_base = cds_yaxis
        exon_height = 0.1
        if args.no_background:               
            fig.patch.set_visible(False)  # Removes figure background
            ax.patch.set_visible(False)   # Removes axis background
            
    # Plot classified regions (middle layer)
    for entry in classifications:
        start = entry['start']
        end = entry['end']
        width = end - start

        if entry['type'] == 'simple_repeat':
            at_content = entry['at_content']
            color = (
                0 + (0.7 - 0) * (1 - at_content / 100),  # Red component
                0 + (0.9 - 0) * (1 - at_content / 100),  # Green component
                0.5 + (1 - 0.5) * (1 - at_content / 100)  # Blue component
            )
            ax.add_patch(
                patches.Rectangle(
                    (start, repeat_y_base), width, repeat_height, facecolor=color, zorder=2
                )
            )
        elif entry['type'] == 'TE':
            ax.add_patch(
                patches.Rectangle(
                    (start, repeat_y_base), width, repeat_height,
                    facecolor=args.TE_color, edgecolor=None, linewidth=0, zorder=2
                )
            )
        elif entry['type'] == 'repeat_regex':
            ax.add_patch(
                patches.Rectangle(
                    (start, repeat_y_base), width, repeat_height,
                    facecolor=args.repeat_color, edgecolor=None, linewidth=0, zorder=3
                )
            )
        elif entry['type'] == 'NNNN_block':
            nnnn_y = repeat_y_base + repeat_height * 0.3
            nnnn_height = repeat_height * 0.4
            ax.add_patch(
                patches.Rectangle(
                    (start, nnnn_y), width, nnnn_height,
                    facecolor="lightgreen", edgecolor=None, linewidth=0, zorder=2
                )
            )
        elif entry['type'] == 'gene' and len(exons) == 0:  # no external exon location file
            ax.add_patch(
                patches.Rectangle(
                    (start, exon_y_base + exon_height + 0.15), width, exon_height,
                    facecolor=args.gene_color, edgecolor=None, linewidth=0, zorder=3
                )
            )

    # Plot exons from external CDS file (bottom layer)
    if exons:
        print(f"Adding {len(exons)} exons to visualization")
        for exon in exons:
            start = exon['start']
            end = exon['end']
            width = end - start
            exon_id = exon['id']

            # Draw exon rectangle
            ax.add_patch(
                patches.Rectangle(
                    (start, exon_y_base), width, exon_height,
                    facecolor=args.gene_color, edgecolor=None, linewidth=0, zorder=3
                )
            )

            # Add exon label above the rectangle
            mid_point = (start + end) / 2
            if has_depth:
                label_y = exon_y_base + exon_height + 0.05  # Small gap above exon rectangle
            else:
                # label_y = exon_y_base + exon_height + 0.15
                label_y = exon_y_base + exon_height + 0.05

            ax.text(mid_point, label_y, f"{exon_id}",
                   ha='center', va='bottom', fontsize=args.exon_fontsize, fontweight=args.exon_fontweight)

    # Draw baseline for unmatched sequences (single-copy)
    ax.plot([args.Xmin, total_length], [baseline_y, baseline_y], color="black", lw=0.5, zorder=1)

    # Configure plot appearance
    if has_depth:
        # Configure annotation subplot - no left spine
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)  # Remove left spine from annotation subplot
        ax.spines['bottom'].set_visible(True)  # Ensure bottom spine is visible
        ax.set_yticks([])  # No y-axis labels or spine
        
        # Configure depth subplot - keep left spine only here
        ax_depth.spines['top'].set_visible(False)
        ax_depth.spines['right'].set_visible(False)
        ax_depth.spines['bottom'].set_visible(False)
        ax_depth.spines['left'].set_visible(True)  # Keep left spine only for depth plot
        
        # Set x-axis formatting only on bottom subplot
        tick_interval_minor = args.tick_minor
        tick_interval_major = args.tick_major
        minor_ticks = np.arange(args.Xmin, total_length, tick_interval_minor)
        major_ticks = np.arange(args.Xmin, total_length, tick_interval_major)
        # major_labels = [f"{int(tick/1000)} kb" for tick in major_ticks]
        major_labels = [f"{int(tick/1000)}" for tick in major_ticks] #put "kb" only in the  last label
        major_labels[-1] = major_labels[-1].strip() + " kb" 

        
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_xticks(major_ticks)
        ax.set_xticklabels(major_labels, fontsize=args.scale_fontsize)
        
        # Customize tick appearance
        ax.tick_params(axis="x", which="minor", length=4, width=0.5)
        ax.tick_params(axis="x", which="major", length=6, width=1.0)
        
    else:
        # Original single plot configuration
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Set x-axis formatting
        tick_interval_minor = args.tick_minor
        tick_interval_major = args.tick_major
        minor_ticks = np.arange(args.Xmin, total_length, tick_interval_minor)
        major_ticks = np.arange(args.Xmin, total_length, tick_interval_major)
        # major_labels = [f"{int(tick/1000)} kb" for tick in major_ticks]
        major_labels = [f"{int(tick/1000)}" for tick in major_ticks] #put "kb" only in the  last label
        major_labels[-1] = major_labels[-1].strip() + " kb" 

        ax.set_xticks(minor_ticks, minor=True)
        ax.set_xticks(major_ticks)
        ax.set_xticklabels(major_labels, fontsize=args.scale_fontsize)
        
        # Customize tick appearance
        ax.tick_params(axis="x", which="minor", length=4, width=0.5)
        ax.tick_params(axis="x", which="major", length=6, width=1.0)
        ax.set_yticks([])
        
        if args.no_scale:
            ax.xaxis.set_visible(False)
            ax.spines['bottom'].set_visible(False)
            # plt.subplots_adjust(bottom=0.4)  # Reduce bottom margin
            # plt.subplots_adjust(top=0.5)  # Reduce bottom margin
            # plt.subplots_adjust(left=0.5)  # Reduce bottom margin

    # Adjust subplot spacing - minimal gap
    if has_depth:
        plt.subplots_adjust(hspace=0.02)  # Very minimal space between subplots
    
    if args.no_scale:
        plt.tight_layout(pad=2.1, rect=[0, 0, 1, 1])  # rect: [left, bottom, right, top]
    else:
        plt.tight_layout()
    
    if args.no_background:     
        plt.savefig(output_file, dpi=600,  facecolor='none', edgecolor='none')
    else:
        plt.savefig(output_file, dpi=600)           
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Censor output and classify repeats.")
    parser.add_argument("--censor", required=True, help="Path to the Censor output file.")
    parser.add_argument("--len", type=int, required=True, help="Total length of the FASTA sequence.")
    parser.add_argument("--output", required=True, help="graph output file for visualization. File extension dtermines formnat PNG/PDF/SVG")
    parser.add_argument("--len_graph", type=int, default=0,  help="Total length of the FASTA sequence for the graph.")
    parser.add_argument("--tick_minor", type=int, default=5000,  help="tick_interval_minor (unlabelled marks) for the graph.")
    parser.add_argument("--tick_major", type=int, default=10000,  help="tick_interval_major (labelled marks) for the graph.")
    parser.add_argument("--gene_regex", type=str, default='Ycds',  help="gene names regex to detect gene hits.")
    parser.add_argument("--gene_color", type=str, default='red',  help="color that marks gene sequences")
    parser.add_argument("--repeat_regex", type=str, default='not_used',  help="repeat name  regex to detect specific  repeats. eg, 1.688_GK.")
    parser.add_argument("--repeat_exact", type=str, default='not_used',  help="repeat name        to detect specific  repeats. eg, (AAGAG)n.")
    parser.add_argument("--repeat_color", type=str, default='orange',  help="color that marks specific  repeats. eg, 1.688_GK.")
    parser.add_argument("--TE_color", type=str, default='lightcoral',  help="color that marks trnasposable elements.")  
    parser.add_argument("--exclude_regex", type=str, default='not_used',  help="  regex to ignore specific  repeats")
    parser.add_argument("--min_repeat_len", type=int, default=0,  help="repeats smaller than this value will be ignored")
    # NEW PARAMETERS FOR EXTERNAL CDS FILE
    parser.add_argument('--cds_file', type=str, default='', help='optional file containing exon positions (start end exon_id) to display on the plot')
    parser.add_argument('--cds_Yaxis', type=float, default=0.65, help='Y-axis position where exon rectangles will be displayed (default: 0.65)')
    parser.add_argument('--Xmin', type=int, default=0,  help="start of the graph (bp)")
    parser.add_argument('--Xmax', type=int, default=0,  help="end of the graph (bp)")

    # NEW PARAMETER FOR DEPTH FILE
    parser.add_argument('--depth_file', type=str, default='', help='optional file containing depth/coverage data (position depth) to display as a line plot')
    parser.add_argument('--smooth_depth', type=float, default=0, help='smoothing factor for depth plot (0 = no smoothing, higher values = more smoothing, e.g., 1.0-5.0)')
    parser.add_argument('--Ymax', type=float, default=0, help='Y max  for depth plot (if not set, the pgm will use the max depth of the data')    
    parser.add_argument('--title', type=str, default='', help='graph title')
    # parser.add_argument('--subtitle', type=str, default='', help='graphg subtitle')
    parser.add_argument('--no_scale', action='store_true', help='Disable printing the scale')
    parser.add_argument('--no_background', action='store_true', help='Disable white background in the output file')

    parser.add_argument("--exon_fontsize", type=float, default=10, help="font size of the exon labels")
    parser.add_argument("--exon_fontweight", type=str, default='normal', help="font weight of the exon labels   normal/bold/heavy/light/ultrabold/ultralight")
    parser.add_argument("--scale_fontsize", type=float, default=10, help="font size of the X-axis scale")
    parser.add_argument('--Ylabel', type=str, default='Coverage', help='label of the Y-axis')
    parser.add_argument("--label_fontsize", type=float, default=10, help="font size of the axes labels")
    parser.add_argument("--linewidth", type=float, default=1.2, help="linewidth of the coverage data")
    parser.add_argument("--title_fontsize", type=float, default=9, help="font size of the title")   # 6jul2026
    parser.add_argument("--title_color", type=str, default="black", help="font color of the title")   # 7jul2026


    parser.add_argument("--figsize_W", type=float, default=10, help="fig size WIDTH (in inches)")
    parser.add_argument("--figsize_H_depth", type=float, default=5, help="fig size HEIGHT when depth data is plotted (in inches)")
    parser.add_argument("--figsize_H_nodepth", type=float, default=2, help="fig size HEIGHT  when NO depth data is plotted(in inches)")

    parser.add_argument("--print_sat", type=str, default="yes", help="print satellite annd block size")  #27jun2026



    args = parser.parse_args()
    if args.len_graph == 0:
        args.len_graph = args.len
    if args.Xmax == 0:
        args.Xmax = args.len_graph

    # Parse external CDS file if provided
    exons = []
    if args.cds_file:
        exons = parse_cds_file(args.cds_file)

    # Parse depth file if provided
    depth_positions, depth_values = [], []
    if args.depth_file:
        depth_positions, depth_values = parse_depth_file(args.depth_file)

    classifications, satellite_summary = process_censor_output(args.censor, args.len)
    summarize_categories(classifications, args.len, satellite_summary)
    create_visualization(classifications, args.Xmax, args.output, exons, args.cds_Yaxis, depth_positions, depth_values)
