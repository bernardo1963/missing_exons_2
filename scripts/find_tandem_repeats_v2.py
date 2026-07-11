#!/usr/bin/env python3
# find_tandem_repeats_v2.py  created with help of  DeepSeek   26jul2025   v. 1oct2025
# Used to perform the analyses of the ms. "Triplex DNA and inverted repeats cause long-read sequencing bias against satellite DNA", by AB Carvalho, B Kim, F Uno
# based on the regex r'([ATGC]{2}?)\1+'  from McGinty, Lyskov and Mirkin 2025  NAR  https://doi.org/10.1093/nar/gkaf619
# usage: find_tandem_repeats_v2.py  Mitf_generegion_HiFi.fasta --min_tandem_copies_block 4  --print_mode summary
#        find_tandem_repeats_v2.py  LILAP_hifiasm_15sep2025.fasta --min_match_size 30 --min_tandem_copies_block 5  --print_mode blocks 
#        cat Mitf_generegion_HiFi.fasta | find_tandem_repeats_v2.py --min_tandem_copies_block 4  --print_mode summary

# output:
# sequenceID                  monomer  mon_bp  copies   start     end    size  orig_match
# Mitf_generegion_HiFi        AAATTAT       7      60    2870    3289     420     AATTATA
# Mitf_generegion_HiFi        AATTATC       7      15    3290    3394     105     AATTATC
# Mitf_generegion_HiFi        AAATTAT       7       4    3395    3422      28     AATTATA


import re
import argparse
from Bio import SeqIO
from collections import defaultdict
import gzip
import pandas as pd
import sys
import tempfile
from io import StringIO

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence"""
    complement = str.maketrans('ATCGN', 'TAGCN')
    return sequence.translate(complement)[::-1]
    
def get_lex_smallest_monomer(monomer):
    """Returns the lexicographically smallest rotation of a repeat monomer, including the reverse-complement"""
    # For dinucleotides, just compare forward and reverse
    # if len(monomer) == 2:
    # return min(monomer, monomer[::-1])
    # For longer monomers, find all rotations and return the smallest
    rotations_F = [monomer[i:] + monomer[:i] for i in range(len(monomer))]
    # adding the reverse complement
    monomer_rc = reverse_complement(monomer)
    rotations_R = [monomer_rc[i:] + monomer_rc[:i] for i in range(len(monomer_rc))]
    rotations = rotations_F + rotations_R
    canonical_monomer = min(rotations)
    orientation=""
    if canonical_monomer in rotations_F:
        orientation += "F"
    if canonical_monomer in rotations_R:
        orientation += "R"        
    return (canonical_monomer,orientation)

def find_tandem_repeats(seq_file, min_tandem_copies_block=4):    
    """
    Finds tandem repeats in DNA sequences with configurable minimum repeats. monomer len: 2 bp - 30 bp
    Reports the lexicographically smallest version of each repeat monomer.
    """
    if args.print_mode == "blocks" and args.merge_blocks != -1:
        data = []  # will be used to store the output lines 
    pattern = re.compile(r'(([ATGCN]{2,30}?)\2{%d,})' % (min_tandem_copies_block-1))
    overall_cumulative_cp_dict = defaultdict(int) # int is the factory function that returns 0 for new keys
    overall_max_cp_dict = defaultdict(int)     
    if args.print_mode == 'blocks' and args.merge_blocks == -1:
        print("{:<25} {:>20} {:>8} {:>8} {:>8} {:>8} {:>8} {:>6} {:>6}".format("sequenceID","monomer","mon_bp","copies","start","end","size","orient","seq_bp"))
    if args.print_mode == 'summary':
        print("\n\nSummary:")
    if args.print_mode == 'summary' or args.print_mode == 'summary_short' :
        print("{:<25} {:>20} {:>8}{:>8} {:>8} {:>8} {:>8} {:>15} {:>12}".format("sequenceID","monomer","mon_bp","tot_cp","tot_bp","max_cp","max_bp","max_start","max_end"))    
    seq_counter = 0
    
    # MODIFIED: Handle both file input and stdin piping with proper format detection
    if seq_file == '-':
        # Read from stdin - we need to handle the non-seekable stream
        # Read all data from stdin into memory
        stdin_data = sys.stdin.read()
        
        # Use StringIO to create a file-like object that we can seek on
        file_handle = StringIO(stdin_data)
        
        # Detect format by first character
        first_char = file_handle.read(1)
        file_handle.seek(0)  # Rewind to start
        file_format = 'fasta' if first_char == '>' else 'fastq'
        records = SeqIO.parse(file_handle, file_format)
        
    else:
        # Read from file (original logic)
        file_handle = gzip.open(seq_file, 'rt') if seq_file.endswith('.gz') else open(seq_file, 'r')
        # Detect format by first character
        first_char = file_handle.read(1)
        file_handle.seek(0)  # Rewind to start
        file_format = 'fasta' if first_char == '>' else 'fastq'
        records = SeqIO.parse(file_handle, file_format)
    
    seq_counter = 0
    for record in records:
        seq_counter += 1
        cumulative_cp_dict = defaultdict(int) # int is the factory function that returns 0 for new keys
        max_cp_dict = defaultdict(int) 
        max_cp_start_dict = defaultdict(int)
        max_cp_end_dict = defaultdict(int)
        sequence = str(record.seq).upper()
        seq_bp = len(record) 
        for match in pattern.finditer(sequence):
            full_match = match.group(1)
            original_monomer = match.group(2)
            canonical_monomer,orientation = get_lex_smallest_monomer(original_monomer)
            if len(canonical_monomer) == 2 and re.match(r'(AA|GG|CC|TT|NN)', canonical_monomer):
                canonical_monomer = canonical_monomer[0]
            copy_number = len(full_match) // len(canonical_monomer)
            start_pos = match.start(1) + 1  # 1-based position
            end_pos = match.end(1)          # already behaves as 1-based
            if len(full_match) >= args.min_match_size:
                if args.print_mode == 'blocks':
                    if  args.merge_blocks == -1:
                        print("{:<25} {:>20} {:>8} {:>8} {:>8} {:>8} {:>8} {:>6} {:>6}".format(record.id,canonical_monomer,len(canonical_monomer),copy_number,start_pos,end_pos,len(full_match),orientation,seq_bp))
                    else:
                        data.append((record.id, canonical_monomer, len(canonical_monomer),copy_number, start_pos, end_pos, len(full_match),orientation, seq_bp))
                elif args.print_mode == 'censor':
                    if canonical_monomer == 'N':
                        sat_censor_string = 'NNNN_block'
                    else:
                        sat_censor_string = "(" + canonical_monomer + ")n"
                    output_line = f"{record.id}\t{start_pos}\t{end_pos}\t{sat_censor_string}\t{len(full_match)}\t0\tc\t1.0000\t0\t0\t0\t0"  
                    print(output_line) 
                elif args.print_mode == 'summary' or args.print_mode == 'summary_short': 
                    cumulative_cp_dict[canonical_monomer] += copy_number
                    if copy_number > max_cp_dict[canonical_monomer]:
                        max_cp_dict[canonical_monomer] = copy_number
                        max_cp_start_dict[canonical_monomer] = start_pos
                        max_cp_end_dict[canonical_monomer] = end_pos               
                    # max_cp_dict[canonical_monomer] = max(max_cp_dict[canonical_monomer], copy_number)
                    overall_cumulative_cp_dict[canonical_monomer] += copy_number
                    overall_max_cp_dict[canonical_monomer] = max(overall_max_cp_dict[canonical_monomer], copy_number)          
                    
        if args.print_mode == 'summary' or args.print_mode == 'summary_short' :
            flag_sat_block_found = False
            sorted_dict = sorted(cumulative_cp_dict.items(), key=lambda item: item[1], reverse=True)  # to list in decreasing value 
            for sat , total_copy in sorted_dict:
                if max_cp_dict[sat] >= args.min_tandem_copies_summary:            
                    flag_sat_block_found = True
            if flag_sat_block_found :   # avoid printing headers with no filtered data                
                for sat , total_copy in sorted_dict:
                    if max_cp_dict[sat] >= args.min_tandem_copies_summary:
                        total_bp = total_copy*len(sat)
                        max_bp = max_cp_dict[sat]*len(sat)                    
                        print("{:<25} {:>20} {:>8} {:>8} {:>8} {:>8} {:>8} {:>12} {:>12}".format(record.id,sat,len(sat),total_copy,total_bp,max_cp_dict[sat],max_bp,max_cp_start_dict[sat],max_cp_end_dict[sat]))
                if args.print_mode == 'summary':
                    print('')    
    if seq_counter > 1 and  args.print_mode == 'summary':
        # MODIFIED: Handle filename for stdin case
        if seq_file == '-':
            fasta_name_short = "stdin"
        else:
            fasta_name_short = seq_file.replace('.fasta','')
        overall_sorted_dict = sorted(overall_cumulative_cp_dict.items(), key=lambda item: item[1], reverse=True)  # to list in decreasing value 
        print("\n\nOverall summary across",seq_counter,"sequences:")
        print("{:<25} {:>20} {:>8} {:>8} {:>8} {:>8} {:>8}".format("sequenceID","monomer","mon_bp","tot_cp","tot_bp","max_cp","max_bp"))
        for sat , total_copy in overall_sorted_dict:
            if overall_max_cp_dict[sat] >= args.min_tandem_copies_summary:
                total_bp = total_copy*len(sat)
                max_bp = overall_max_cp_dict[sat]*len(sat)                    
                print("{:<25} {:>20} {:>8} {:>8} {:>8} {:>8} {:>8}".format(fasta_name_short,sat,len(sat),total_copy,total_bp,overall_max_cp_dict[sat],max_bp))
    
    # MODIFIED: Always close the file handle
    file_handle.close()

    if args.print_mode == "blocks" and args.merge_blocks != -1:
        original_df = pd.DataFrame(data, columns=['sequenceID', 'monomer', 'mon_bp', 'copies', 'start', 'end', 'size', 'orient', 'seq_bp'])
        # Sort the dataframe to ensure proper ordering
        original_df = original_df.sort_values(['sequenceID', 'monomer', 'start', 'orient']).reset_index(drop=True)
        merged_df = merge_blocks(original_df, args.merge_blocks)
        # Reorder columns to have mrged blocks (m_blks) at the end
        columns = [col for col in merged_df.columns if col != 'm_blks'] + ['m_blks']
        merged_df = merged_df[columns]
        # Output results
        format_output(merged_df)
    
  
        
def format_output(df):
    """Format the output with nice spacing between columns"""
    # Print header
    header_format = "{:<25} {:>20} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}"
    print(header_format.format(
        "sequenceID", "monomer", "mon_bp", "copies", "start", "end", "size", "orient", "seq_bp", "m_blks"
    ))
    
    # Print each row
    for _, row in df.iterrows():
        row_format = "{:<25} {:>20} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}"
        print(row_format.format(
            str(row['sequenceID']),
            str(row['monomer']),
            int(row['mon_bp']),
            int(row['copies']),
            int(row['start']),
            int(row['end']),
            int(row['size']),
            str(row['orient']),
            int(row['seq_bp']),
            int(row['m_blks'])
        ))       
        
        


def merge_blocks(df, merge_blocks):
    """merge blocks that are close and have the same sequenceID, monomer, and orientation"""
    
    # Initialize the combined blocks column
    df['m_blks'] = 1
    
    # List to store merged blocks
    merged_blocks = []
    
    # Group by sequenceID, monomer, and orient since we only check adjacent blocks within these groups
    current_group = None
    current_blocks = []
    
    for idx, row in df.iterrows():
        # Check if we're starting a new group
        group_key = (row['sequenceID'], row['monomer'], row['orient'])
        
        if current_group is None:
            current_group = group_key
            current_blocks = [row]
        elif group_key == current_group:
            # Same group, check if we can merge with the last block
            last_block = current_blocks[-1]
            
            # Check distance condition
            distance = row['start'] - last_block['end']
            
            if distance <= merge_blocks:
                # merge the blocks - update the last block in current_blocks
                merged_block = last_block.copy()
                merged_block['end'] = row['end']
                merged_block['size'] = row['end'] - merged_block['start'] + 1
                merged_block['copies'] += row['copies']
                merged_block['m_blks'] += 1
                
                # Replace the last block with the merged version
                current_blocks[-1] = merged_block
            else:
                # Cannot merge, add as new block
                current_blocks.append(row)
        else:
            # New group, save current blocks and start new group
            merged_blocks.extend(current_blocks)
            current_group = group_key
            current_blocks = [row]
    
    # Don't forget the last group
    if current_blocks:
        merged_blocks.extend(current_blocks)
    
    # Create new DataFrame
    result_df = pd.DataFrame(merged_blocks)
    # Reset index
    result_df.reset_index(drop=True, inplace=True)
    return result_df




def main():
    parser = argparse.ArgumentParser(description='Find tandem repeats in DNA sequences',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # MODIFIED: Make seq_file optional and accept '-' for stdin
    parser.add_argument('seq_file', nargs='?', help='Input file containing DNA sequences (fasta or fastq). Use - or omit for stdin')
    parser.add_argument('-m', '--min_tandem_copies_block',type=int,default=4, help='Minimum number of tandem repeats to report')
    parser.add_argument('--min_match_size',type=int,default=20, help='Minimum size of the match to report in individua sequences(avoids reporting (AT)4')
    parser.add_argument('-o', '--output',help='Output file (default: print to stdout)')
    parser.add_argument('--print_mode',type=str, default='blocks' , help='summary/blocks/censor/summary_short   summary prints the sum per satellite type. summary_short is more suitable to additional  processing (skips emply lines, etc)  ')
    parser.add_argument('--min_tandem_copies_summary',type=int,default=50, help='Minimum tandem_copies to be reported in the summary (avoids small blocks')
    parser.add_argument('--merge_blocks',type=int,default=-1, help='maximun distance between two satellite blocks that will be merged. Only works with --mode blocks')
    global args
    
    args = parser.parse_args()
    
    # MODIFIED: Handle default case (stdin)
    if args.seq_file is None:
        args.seq_file = '-'
    
    if args.output:
        with open(args.output, 'w') as f:
            import sys
            sys.stdout = f
            find_tandem_repeats(args.seq_file, args.min_tandem_copies_block)
    else:
        find_tandem_repeats(args.seq_file, args.min_tandem_copies_block)

if __name__ == "__main__":
    main()