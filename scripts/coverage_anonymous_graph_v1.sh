#!/bin/bash
# coverage_anonymous_graph_v1.sh  Bernardo  20ago2025 v 31may2026
#  Used to perform the analyses of the ms. "Triplex DNA and inverted repeats cause long-read sequencing bias against satellite DNA", by AB Carvalho, B Kim, F Uno
# WARNING: paths of bam and fasta files are hard-coded, and should be edited to conform the user's server 
# WARNING:coverage_anonymous_graph_v1.sh is a wrapper script that invokes other programs like censor, visualize_repeats_censor5.py . seqtk, etrc. These programs should be available in the path of the user. 
# usage: coverage_anonymous_graph_v1.sh   ptg000003l  0 30000     



set -e

function PrintUsage() {
   echo 'Usage: coverage_anonymous_graph_v1.sh -r HiFi -a hifiasm ptg000045l 0 221980 '
   echo 'Usage: coverage_anonymous_graph_v1.sh -r HiFi -a ONT_hifiasm  -F AAGAG_1.blocks.txt -V "--Ymax 100"  '   
   exit 1
}

function process_data() {
    local seqID="$1"
    local st="$2"
    local en="$3"
	seqname="${seqID}_${st}_${en}"
	# echo $seqname
	if [[ $name_suffix = "" ]]; then
		seqname="${seqID}_${st}_${en}"
	else
		seqname="${seqID}_${st}_${en}_$name_suffix"
	fi	
    [[ "$batch_file" != "" ]] && seqname="${seqname}_$monomer"     # avoids over-writing when different satellites occur in the same contig (and coordinates are the same; eg, full contig)
	region="${seqID}:${st}-${en}"
	echo "point 0   region:" $region
	echo -e "$seqID\t$st\t$en" | seqtk subseq $genomic_fasta -    > $seqname.fasta
	find_tandem_repeats_v2.py  $seqname.fasta --min_match_size 30 --min_tandem_copies_summary 5  --print_mode blocks > $seqname.repeats_mms30mts5.txt
	[[ $debug == 3 ]] && return    # only wants  find_tandem_repeats_v2.py output, with the same cut-offs of the one used in the graphs  7jul226
	# size=`fasta_size3.awk $seqname.fasta`  
	size=`seqtk comp $seqname.fasta | awk '{print $2}'`  
	samtools depth $genomic_bam_file -a  -r $region | awk -v offset=$st '{print $2-offset,$3}' > $seqname.depth
	[[ $debug == -1 ]] && return   # I want only the depth data
	# /home6/tools/miniconda3/envs/censor/bin/censor  $seqname.fasta -lib dro -lib inv   -mode norm   -bprm 'cpus=100'  -bprm '-filter=none'  -nomasked   -nofound -minsim 0.90 &> /dev/null
	# /home6/tools/miniconda3/envs/censor/bin/censor  $seqname.fasta  -lib diptera2   -mode norm   -bprm 'cpus=100'  -bprm '-filter=none'  -nomasked   -nofound -minsim 0.90 &> /dev/null
	if [[ $assembly = "CHM13" ]]; then	
		/home6/tools/miniconda3/envs/censor/bin/censor  $seqname.fasta  -lib hum     -mode norm  -bprm 'cpus=100'   -nomasked   -nofound -minsim 0.90 &> /dev/null
	elif [[ $assembly = "simulated" ]]; then	
		/home6/tools/miniconda3/envs/censor/bin/censor  $seqname.fasta  -lib hum -lib diptera2    -mode norm  -bprm 'cpus=100'   -nomasked   -nofound -minsim 0.90 &> /dev/null	
	else
		/home6/tools/miniconda3/envs/censor/bin/censor  $seqname.fasta  -lib diptera2 -mode norm -bprm 'cpus=100'   -nomasked   -nofound -minsim 0.90 &> /dev/null
	fi	
	find_tandem_repeats_v2.py  $seqname.fasta --min_match_size 30 --min_tandem_copies_summary 4  --print_mode summary_short > $seqname.repeats_summary_short.txt  # name changed in 7jul2026  from $seqname.repeats.txt  to $seqname.repeats_summary_short.txt   summary_short mode is a bit weird
	echo "point 1"
	if [[ -z "$satellite_force" ]]; then
		satellite=$(awk '(NR==2){print $2}' "$seqname.repeats_summary_short.txt")
		# Check if the command above failed?
		if [[ -z "$satellite" ]]; then
			echo "WARNING: Could not extract satellite from $seqname.repeats_summary_short.txt" >&2
		fi
	else
		satellite="$satellite_force"
	fi	
	echo "point 2   satellite: $satellite"
	set +e
	# sat_info=`grep -w $satellite $satellite_data_file`
	sat_info=$(grep -w "$satellite" "$satellite_data_file" 2>/dev/null || echo "$satellite no non-B data")
	set -e
	echo "point 3"
	sat_composition=`awk '(NR>1){sat_array[$2]=$5; tot_sat += $5}; END{PROCINFO["sorted_in"]="@val_num_desc"; for (sat in sat_array){i++; sat_perc=int(100*sat_array[sat]/tot_sat); if ((sat_perc>10)||(i<=3)){printf("%s %s%  ",sat,sat_perc)}}}' $seqname.repeats_summary_short.txt`
	echo "point 4" 
	find_tandem_repeats_v2.py  $seqname.fasta $find_repeats_parameters  --print_mode censor >> $seqname.fasta.map
	
	echo "point 5"
	if [[ "$mRNA_bam_file" != "no_mRNA" ]]; then
		samtools view -h  $mRNA_bam_file   $region | bedtools bamtobed -bed12 -i - > $seqname.bed12
		echo "point 6"
		# awk -v offset=$st '{internal_offset= $2 -offset
			# n = split($11, blockSizes, ","); split($12, blockStarts, ",");
			# for (i=1; i<=n; i++){exon_start = internal_offset + blockStarts[i]; exon_end = exon_start + blockSizes[i]; 
			# if ((exon_start>=0)&&(array_gene[$4]=="")){ print exon_start, exon_end, gensub(/-R[A-Z]\>/,"",1,$4);array_gene[$4] ++}
			# else if (exon_start>=0) print exon_start, exon_end}}'  $seqname.bed12  > $seqname.exon_locations	
		 # 8jul2026:
		awk -v offset=$st '{internal_offset= $2 -offset
			n = split($11, blockSizes, ","); split($12, blockStarts, ",");
			for (i=1; i<=n; i++){exon_start = internal_offset + blockStarts[i]; exon_end = exon_start + blockSizes[i]; 			
			if (exon_start>=0){print exon_start, exon_end, gensub(/-R[A-Z]\>/,"",1,$4)}}}'  $seqname.bed12  > $seqname.exon_locations.raw			
        awk -v size=$size '{gene_array[$3]++; if ($1>size){next}};(gene_array[$3]==1){print $0}; (gene_array[$3]>1){print $1,$2}'  $seqname.exon_locations.raw >  $seqname.exon_locations					
		if [[ "$flip_sequence" == "yes" ]]; then 
			cds_file_cmd="--cds_file "$seqname.exon_locations.rc
		else
			cds_file_cmd="--cds_file "$seqname.exon_locations
		fi
	else
			cds_file_cmd=""
	fi
	

	if [[ $graph_name2 = "" ]]; then
		if [[ $graph_format == "" ]]; then
			graph_name="$seqname.png"
		else
			graph_name="$seqname.$graph_format"
		fi          # <-- closes the inner if
	else
		graph_name="$graph_name2"
	fi
	visualize_FU_parameters2=$visualize_FU_parameters
	[[ $flag_get_repeat_exact = "yes" ]] && visualize_FU_parameters2="$visualize_FU_parameters2 --repeat_exact (${monomer})n"
	if [[ "$flip_sequence" == "yes" ]]; then
		awk -v size=$size '{print $1,(size-$3),(size-$2),$4,$5,$6,$7,$8,$9,$10,$11,$12}'   $seqname.fasta.map  >  $seqname.fasta.map.rc
		tac  $seqname.depth | awk -v size=$size '{print (size-$1),$2}' > $seqname.depth.rc
		# awk -v size=$size '{print (size-$2),(size-$1),$3}'   $seqname.exon_locations   >  $seqname.exon_locations.rc
		
		awk -v size=$size '($1>size){next}; {print (size-$2),(size-$1),$3}'   $seqname.exon_locations.raw | tac |  awk -v size=$size '{gene_array[$3]++};(gene_array[$3]==1){print $0}; (gene_array[$3]>1){print $1,$2}'  >  $seqname.exon_locations.rc  # 8jul2028

		visualize_repeats_censor5.py  --censor $seqname.fasta.map.rc --min_repeat_len 30  --len $size --output "rc."$graph_name  --depth_file $seqname.depth.rc  $cds_file_cmd  --gene_color black --title "$region       $sat_info  $sat_composition"  $visualize_FU_parameters2  # 2>/dev/null
	else
		visualize_repeats_censor5.py  --censor $seqname.fasta.map    --min_repeat_len 30  --len $size --output $graph_name       --depth_file $seqname.depth     $cds_file_cmd  --gene_color black --title "$region       $sat_info  $sat_composition"  $visualize_FU_parameters2  # 2>/dev/null
	fi

	# echo $seqname  $debug
	[[ $debug == -1 ]] && set -x &&  rm -f *.repeats_mms30mts5.txt # only depth
	[[ $debug == -2 ]] && set -x &&  rm -f censor.*.log $seqname.fasta.aln $seqname.fasta.idx debug_log.txt  $seqname.bed12 $seqname.depth* $seqname.exon_locations*  $seqname.fasta.map*   $seqname.fasta $seqname.repeats_mms30mts5.txt $seqname.repeats_summary_short.txt # only the graph
	[[ $debug == 0 ]] && set -x  &&  rm -f *.log *.aln *.idx debug*  *.bed12 *.depth*
	[[ $debug == 1 ]] && set -x  &&  rm -f *.log *.aln *.idx debug*  *.bed12 
	[[ $debug == 3 ]] && set -x  &&  rm -f $seqname.fasta  # only wants find_tandem_repeats_v2.py output, with the same cut-offs of the one used in the graphs  7jul226 
	if ! [ "$debug" -eq "$debug" ] 2>/dev/null; then # non-numeric debug, for fine control of what will be deleted.  -D "rm *.log"
        eval "$debug"
    fi
	echo "point 7"	
	}


function process_from_file_partial_contig() {  
    local seqID="no_seqID"
    local st="no_st"
    local en="no_en"   
	# Check if file exists
    if [[ ! -f "$batch_file" ]]; then
        echo "Error: File $batch_file not found!" >&2
        exit 1
    fi
    # Read each line from file and process
    while IFS= read -r line; do
        # Skip empty lines and lines starting with #  (comments)
        [[ -z "$line" || "$line" =~ ^# || "$line" =~ sequenceID ]] && continue
        # Extract columns 1 (seqID) 5 (sat_start)  6 sat_end)  9 (contig_size)
        seqID=$(echo "$line" | awk '{print $1}')
		sat_start=$(echo "$line" | awk '{print $5}')
		sat_end=$(echo "$line" | awk '{print $6}')
		contig_size=$(echo "$line" | awk '{print $9}')
		monomer=$(echo "$line" | awk '{print $2}')
		# modified 4oct2025 to: implement $flanking_bp  ; preserve coordintes when -g 0
        if (( sat_start - flanking_bp > 0 )); then st=$((sat_start - flanking_bp)); else st=0; fi
        if (( $contig_size < $((sat_end + flanking_bp)) )); then en=$contig_size; else en=$((sat_end + flanking_bp)); fi
		if (( contig_size <= 100000 && flanking_bp > 0 )); then process_data "$seqID" 0 "$contig_size"; else process_data "$seqID" "$st" "$en"; fi	
        # if (( $sat_start - 50000 > 0 )); then st=$((sat_start - 50000)); else st=0; fi
        # if (( $contig_size < $((sat_end + 50000)) )); then en=$contig_size; else en=$((sat_end + 50000)); fi
		# if (( $contig_size <= 100000 )); then process_data "$seqID" 0 "$contig_size"; else process_data "$seqID" "$st" "$en"; fi  # modified in 28set2025
    done < "$batch_file"
	}





function process_from_file_full_contig() {  
    local seqID="no_seqID"
    local st=0
    local en="no_en"   
	# Check if file exists
    if [[ ! -f "$batch_file" ]]; then
        echo "Error: File $batch_file not found!" >&2
        exit 1
    fi
    
    # Read each line from file and process
    while IFS= read -r line; do
        # Skip empty lines and lines starting with #  (comments)
        [[ -z "$line" || "$line" =~ ^#  || "$line" =~ sequenceID ]] && continue
        # Extract columns 1 (seqID)   9 (contig_size)
        seqID=$(echo "$line" | awk '{print $1}')
        st=0
		en=$(echo "$line" | awk '{print $9}')
		monomer=$(echo "$line" | awk '{print $2}')
        # Process the data
		echo  $seqID   $st $en
		process_data "$seqID" "$st" "$en"		
    done < "$batch_file"
	}


function process_from_args() {
    if [[ $positional_args_count -ne 3 ]]; then
        echo "Error: Expected 3 positional arguments when -F is not used" >&2
        echo "Usage: $0 [-F filename] contig start end" >&2
        exit 1
    fi
    process_data "$1" "$2" "$3"
	}


# main code
if [ "$#" == "0" ]; then
	PrintUsage
	exit 1
fi

while getopts   "hr:d:b:m:s:S:a:g:F:V:n:R:N:D:f:v:t:" OPTION   
do
   case $OPTION in
        h) PrintUsage
         ;;
        r) reads=$OPTARG 
         ;;
        a) assembly=$OPTARG   
         ;;
        b) genomic_bam_file2=$OPTARG   # needs to adjust genomic_fasta_file 
         ;;
        m) mRNA_bam_file2=$OPTARG   
         ;;		 
        # s) satellite_force=$OPTARG   # force a satelliete, instead of relying ion automatic detection. Useful when the target sat is minoritary 
         # ;;   		 
        S) satellite_data_file=$OPTARG
         ;;     
		D) debug=$OPTARG  # debug mode  options: -2 erase everything except  the graph file    -1 only produce the depth file  0: default: removes intermediate  files       1:  do not erase temporary files  2:  do not erase temporary files + set -x    ,
		;;
		n) graph_name2=$OPTARG  # output graph name
		;;
		N) name_suffix=$OPTARG  # suffix that will be added to all file names (eh, LILAP_HiFi)
		;;		
		F) batch_file=$OPTARG  # file in the format of ONThifiasm_mms50_mtc5_blocks.txt that will provide the contig and the size (instead of providing the coordinates as positional parameters)
		;;       
		g) flanking_bp=$OPTARG  # Must be used with -F. modify its behavior: we get  say, 50kb before start, and 50kb after end, replacing thsese values if we get illegal valus (limits below 0 or above the ctg size
		;;
		V) visualize_FU_parameters=$OPTARG  # additional  parameters that will be passed to FU program visualize_repeats_censor5.py  "--Ymax 100 " 
		;;	
		R) find_repeats_parameters=$OPTARG  #   parameters that will be passed to  program find_tandem_repeats_v2.py   " --min_match_size 30 --min_tandem_copies_block 5 " 
		;;						
		f) graph_format=$OPTARG  #  svg/png/pdf   default: png
		;;
		v) flip_sequence=$OPTARG  #   default: not flip 
		;;		
		# t) graph_title=$OPTARG  #   default:   "$region       $sat_info  $sat_composition"
		# ;;				
	esac
done

# Shift off the processed options  Needed to use the positional argunmtes correctly. eg, in coverage_anonymous_graph_v0.sh -r HiFi -a hifiasm  ptg000003l  0 30000    , ptg000003l should be b$1 again
shift $((OPTIND - 1))
positional_args_count=$#

# setting default values if not specified by the user:
[[ $assembly == "" ]] && assembly="hifiasmL2"
[[ $reads == "" ]] && reads="HiFi"
[[ $find_repeats_parameters == "" ]] && find_repeats_parameters=" --min_match_size 30 --min_tandem_copies_block 5 "
[[ $debug == "" ]] && debug=0
[[ $debug == 2 ]] && set -x
if [[ "$flip_sequence" == "" ]]; then
	flip_sequence="no"
else
    flip_sequence="yes"	
fi	
# [[ $graph_title == "" ]] && graph_title="\"$region       $sat_info  $sat_composition\""

# setting the genomic_bam_file:
if [[ "$genomic_bam_file2" == "" ]]; then
	if [[ $assembly = "hifiasm" ]] && [[ $reads = "HiFi" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/HiFi_hifiasmL0p_primary.sorted.bam"
		genomic_fasta="/draft10/bernardo/YGS_2_ms/mel_HiFi/genome/HiFi_hifiasmL0/HiFi_hifiasmL0_p.fasta" 		
	fi
	if [[ $assembly = "hifiasmL2" ]] && [[ $reads = "HiFi" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/HiFi_hifiasmL2primary.sorted.bam"
		genomic_fasta="/draft10/bernardo/YGS_2_ms/mel_HiFi/genome/HiFi_hifiasmL2primary/HiFi_hifiasmL2primary.fasta" 		
	fi
    if [[ $assembly = "R6" ]] && [[ $reads = "HiFi" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/HiFi_R6_primary.sorted.bam"
		genomic_fasta="/home6/bernardo/projects/missing_exons_2/generegion_cov/dmel-all-chromosome-r6.64.fasta.gz" 
	fi
    if [[ $assembly = "ONT" ]] && [[ $reads = "HiFi" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/HiFi_AXFlye6YFlyecanu_primary.sorted.bam"
		genomic_fasta="/home6/bernardo/projects/missing_exons_2/generegion_cov/AXFlye6_YFlyecanu_14jun24_clean2.fasta2.gz"  
	fi
    if [[ $assembly = "ONT_hifiasm" ]] && [[ $reads = "HiFi" ]]; then
		genomic_bam_file="/home6/bernardo//projects/missing_exons_2/sat_assembly_ONT_hifiasm/anonymous/HiFi_ONThifiasm.sorted.bam"
		genomic_fasta="/draft10/bernardo/YGS_2_ms/mel/genome/ONT45kY20_hifiasm_2sep2025.fasta"  
	fi
    if [[ $assembly = "ONT_hifiasm" ]] && [[ $reads = "ONT" ]]; then
		genomic_bam_file="/home6/bernardo//projects/missing_exons_2/sat_assembly_ONT_hifiasm/anonymous/ONT_ONThifiasm.sorted.bam"
		genomic_fasta="/draft10/bernardo/YGS_2_ms/mel/genome/ONT45kY20_hifiasm_2sep2025.fasta"  
	fi	
    if [[ $assembly = "LILAP_hifiasm" ]] && [[ $reads = "HiFi" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/sat_assembly_LILAP_hifiasm/anonymous/HiFi_LILAPhifiasm.sorted.bam"
		genomic_fasta="/draft10/bernardo/YGS_2_ms/mel_LILAP/assemblies/LILAP_hifiasm_15sep2025.fasta"  
	fi		
    if [[ $assembly = "LILAP_hifiasm" ]] && [[ $reads = "ONT" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/sat_assembly_ONT_hifiasm/anonymous/ONT_LILAPhifiasm.sorted.bam"
		genomic_fasta="/draft10/bernardo/YGS_2_ms/mel_LILAP/assemblies/LILAP_hifiasm_15sep2025.fasta"  
	fi	
    if [[ $assembly = "ONT_hifiasm" ]] && [[ $reads = "LILAP" ]]; then
		genomic_bam_file="/home6/bernardo//projects/missing_exons_2/sat_assembly_ONT_hifiasm/anonymous/LILAP_ONThifiasm.sorted.bam"
		genomic_fasta="/draft10/bernardo/YGS_2_ms/mel/genome/ONT45kY20_hifiasm_2sep2025.fasta" 
	fi
    if [[ $assembly = "CS_ONT_hifiasm" ]] && [[ $reads = "CS_HiFi" ]]; then
		genomic_bam_file="/draft6/bernardo/drosophila/CantonS_T2T/bam_files/CS_HF_ONThifiasm.sorted.bam"
		genomic_fasta="/draft6/bernardo/drosophila/CantonS_T2T/genome_completeness/CS_ONT98k_hifiasmL2primary.fasta" 
	fi
    if [[ $assembly = "CS_nT2T" ]] && [[ $reads = "CS_HiFi" ]]; then
		genomic_bam_file="/draft6/bernardo/drosophila/CantonS_T2T/bam_files/CS_HF_nT2T.sorted.bam"
		genomic_fasta="/draft6/bernardo/drosophila/CantonS_T2T/genome/GCA_048772135.1_ASM4877213v1_genomic.fna" 
	fi	
    if [[ $assembly = "CHM13" ]] && [[ $reads = "HiFi20k" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/human/HiFi20k_CHM13.sorted.bam"
		genomic_fasta="/home6/bernardo/projects/missing_exons_2/human/CHM13_XA.fasta.gz" 
	fi	
    if [[ $assembly = "CHM13" ]] && [[ $reads = "HiFi10k" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/human/HiFi10k_CHM13.sorted.bam"
		genomic_fasta="/home6/bernardo/projects/missing_exons_2/human/CHM13_XA.fasta.gz" 
	fi	
    if [[ $assembly = "CHM13" ]] && [[ $reads = "HiFi10k20k" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/human/HiFi10k20k_CHM13.sorted.bam"
		genomic_fasta="/home6/bernardo/projects/missing_exons_2/human/CHM13_XA.fasta.gz" 
	fi
    if [[ $assembly = "simulated" ]] && [[ $reads = "hypothesis1" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/human/simulated_reads/simulated_hyp1.sorted.bam"
		genomic_fasta="/home6/bernardo/projects/missing_exons_2/human/simulated_reads/simulated_genome.fasta" 
	fi	
    if [[ $assembly = "simulated" ]] && [[ $reads = "hypothesis2" ]]; then
		genomic_bam_file="/home6/bernardo/projects/missing_exons_2/human/simulated_reads/simulated_hyp2.sorted.bam"
		genomic_fasta="/home6/bernardo/projects/missing_exons_2/human/simulated_reads/simulated_genome.fasta" 
	fi	
	[[ $genomic_fasta == "" ]] && echo "Wrong combination of assembly and reads:" $assembly  $reads   && exit 1
fi




# setting the mRNA_bam_file  :
if [[ "$mRNA_bam_file2" == "" ]]; then
  [[ $assembly = "hifiasm" ]]     && mRNA_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/mRNA_hifiasmL0p.sorted.bam"  
  [[ $assembly = "hifiasmL2" ]]   && mRNA_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/mRNA_hifiasmL2primary.sorted.bam"
  [[ $assembly = "R6" ]]          && mRNA_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/mRNA_R6.sorted.bam" 
  [[ $assembly = "ONT" ]]         && mRNA_bam_file="/home6/bernardo/projects/missing_exons_2/generegion_cov/mRNA_HiFi_AXFlye6YFlyecanu.sorted.bam"
  [[ $assembly = "ONT_hifiasm" ]] && mRNA_bam_file="/home6/bernardo//projects/missing_exons_2/sat_assembly_ONT_hifiasm/anonymous/mRNA_ONThifiasm.sorted.bam" 	
  [[ $assembly = "LILAP_hifiasm" ]] && mRNA_bam_file="/home6/bernardo/projects/missing_exons_2/sat_assembly_LILAP_hifiasm/anonymous/mRNA_LILAPhifiasm.sorted.bam" 	 
  [[ $assembly = "CS_ONT_hifiasm" ]] && mRNA_bam_file="/draft6/bernardo/drosophila/CantonS_T2T/bam_files/CS_mRNA_ONThifiasm.sorted.bam" 	
  [[ $assembly = "CS_nT2T" ]] && mRNA_bam_file="/draft6/bernardo/drosophila/CantonS_T2T/bam_files/CS_mRNA_nT2T.sorted.bam" 	
  [[ $assembly = "CHM13" ]] && mRNA_bam_file="no_mRNA"
  [[ $assembly = "simulated" ]] && mRNA_bam_file="no_mRNA"
  [[ $assembly == "" ]] && echo " mRNA alignment not identified for the chosen assembly:" $assembly " Check assembly name"     && exit 1
fi



# setting the satellite_data_file  :
[[ "$satellite_data_file" == "" ]]  && satellite_data_file="/home6/bernardo/projects/missing_exons_2/non-B/sat_19ago25_nonBdata_raw3.txt"


# getting the info for repeat_exact parameter for visualize_FU 
if [[ "$batch_file" != "" && ! "$visualize_FU_parameters" =~ "repeat_exact" ]]; then
    flag_get_repeat_exact="yes"
	echo "flag_get_repeat_exact=yes" 
fi

# Execute appropriate processing based on -F flag
if   [[ "$batch_file" == "" && "$flanking_bp" == "" ]]; then
	 process_from_args "$@"
elif [[ "$batch_file" != "" && "$flanking_bp" == "" ]]; then
	 process_from_file_full_contig "$batch_file"
elif [[ "$batch_file" != "" && "$flanking_bp" != "" ]]; then
	 process_from_file_partial_contig "$batch_file"
else
     echo "wrong combination of -F and -g.  -g can only be used if -F is also set"	
fi


