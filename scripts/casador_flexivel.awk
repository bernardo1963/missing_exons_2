#! /bin/gawk -f
# casador_flexivel.awk    inspired on kmer_mask4.awk    , Bernardo   30oct2017  v.14mar2025
#  Used to perform the analyses of the ms. "Triplex DNA and inverted repeats cause long-read sequencing bias against satellite DNA", by AB Carvalho, B Kim, F Uno
# combines two tabular files using an index chosen by thye  add_file_match_field  and  main_file_match_field parameters
#usage:  casador_flexivel.awk    add_file=prot.accession2taxid  main_file_match_field=2   add_file_match_field=2   Git3_refseqX.minusK1.m8 > Git3_refseqX.minusK1.m8.taxified 


function load_add_file(dummy){
    if (system("test -f " add_file_name)) {print "ERROR: kmer file       " add_file_name    " not found. Check name of the file" > "/dev/stderr" ;  exit }
	# if (add_file_output_fields != ""){add_number_output_fields=split(add_file_output_fields,add_file_output_array)}
	if (add_file_output_fields != ""){add_number_output_fields=split(add_file_output_fields,add_file_output_array,/[, ]/)}
	print "\n" strftime() "  loading add_file ", add_file_name   > "/dev/stderr"  
	print "add file lines found (in millions):" > "/dev/stderr"
	while ( (getline line < add_file_name) > 0 ) {	#loads the add_file file
		total_add ++; add_value=""
		if (total_add==1){
				add_num_fields=split(line,b) #gets the number of fields of the add_file (for printing missing values)					
				}  
		split(line,a) 
		add_index=a[add_file_match_field]
		if (add_file_output_fields==""){add_value=line }
		if (add_file_output_fields ~/[0-9]/){
			for (j=1; j<=add_number_output_fields ; j++){add_value=add_value " " a[add_file_output_array[j]]}
			}
		if (add_file_output_fields=="no_index"){   #the index will not be included in add_value, to avoid the duplication in the output 
			for (j=1; j<=add_num_fields ; j++){
				if (j != add_file_match_field) {add_value=add_value " " a[j]}
				}
			}
		add_array[add_index]=add_value
		if ((total_add % 1000000)==0 ){millions_add_info=total_add / 1000000 ;   printf "%i " , millions_add_info  > "/dev/stderr" ; fflush() }
		}
	print "\n" strftime() "  add_file " add_file_name " loaded." total_add " add info found \n"   > "/dev/stderr"   #prints to stderr so it will not go to a file if output is redirected
	close(add_file_name)
	# for (i in add_array){print i, add_array[i]}; exit
return}




	
function print_usage(dummy){
	print "\n usage:  casador_flexivel.awk    add_file=prot.accession2taxid  main_file_match_field=2   add_file_match_field=2   Git3_refseqX.minusK1.m8 > Git3_refseqX.minusK1.m8.taxified"   > "/dev/stderr"
	print "\n usage:  casador_flexivel.awk    add_file=prot.accession2taxid  main_file_match_field=2   add_file_match_field=2 add_file_output_fields="3"  Git3_refseqX.minusK1.m8 > Git3_refseqX.minusK1.m8.taxified"   > "/dev/stderr"
	print "\n usage:  casador_flexivel.awk    add_file_output_fields=no_index  add_file=temp_acc.txt2 temp_puncs.txt  "
	print "\n usage:  casador_flexivel.awk    add_file=prot.accession2taxid  main_file_match_field=2   add_file_match_field=2 add_file_output_fields=\"3 5\"  Git3_refseqX.minusK1.m8 > Git3_refseqX.minusK1.m8.taxified"   > "/dev/stderr"
	
	return}



BEGIN{
	  OFS="\t"
      start_time =  strftime() 
	  main_file_match_field=1 ;   add_file_match_field=1 
      for (i = 0; i < ARGC; i++){
	    if (ARGV[i] ~ /add_file=/)               {add_file_name=gensub(/.+=/,"","g",ARGV[i])}  #deletes everything until the =   "fasta1_file=melR6_frag_A.PB.fasta" became "melR6_frag_A.PB.fasta"
	    if (ARGV[i] ~ /main_file_match_field=/)  {main_file_match_field=gensub(/.+=/,"","g",ARGV[i])}
	    if (ARGV[i] ~ /add_file_match_field=/)   {add_file_match_field =gensub(/.+=/,"","g",ARGV[i])}
	    if (ARGV[i] ~ /add_file_output_fields=/) {add_file_output_fields =gensub(/.+=/,"","g",ARGV[i])}  # desired fields, either "2 3" or 2,3
	    if (ARGV[i] ~ /line2=/) {line2 =gensub(/.+=/,"","g",ARGV[i])}   # 'snnnnnnnnn'    add a dummy line to force SYSTAT to consider variable string or numeric 	 
	    if (ARGV[i] ~ /add_missing_char=/) {add_missing_char =gensub(/.+=/,"","g",ARGV[i])}   # defalt is "."
		}
      if (add_file_name==""){print "ERROR: add_file not specified." ; print_usage(dummy); exit }
	  if (add_missing_char==""){add_missing_char="."}
	  if (ARGC==1){print_usage(dummy); exit }  
      load_add_file(dummy)
	  main_missing_value=""; add_missing_value=""
      }


( NR==1) {
		main_num_fields=split($0,b)     #gets the number of fileds of the main_file (for printing missing values)
		for (i=1; i<=main_num_fields ; i++){main_missing_value= main_missing_value ". "}   #probably will not be used , at leas in blobtools taxify
		# if (add_file_output_fields ~/[0-9]/){for (i=1; i<=add_number_output_fields ; i++){add_missing_value= add_missing_value ". "}}
		# else{for (i=1; i<=add_num_fields ; i++){add_missing_value= add_missing_value ". "}}
		if (add_file_output_fields ~/[0-9]/){for (i=1; i<=add_number_output_fields ; i++){add_missing_value= add_missing_value  add_missing_char " "}}
		else{for (i=1; i<=add_num_fields ; i++){add_missing_value= add_missing_value  add_missing_char " "}}
		}  
( (NR==2) && (line2 != "") ){
	for (i=1 ; i <= length(line2) ; i++){
		if (substr(line2,i,1)== "s"){printf "string\t"}
		else                        {printf "-1\t"}		#default is numeric (more robust to mispecifictions.
		}
	print ""
	}

(1==1){
	   split($0,main_file)
	   main_file_index=main_file[main_file_match_field]
       if(main_file_index in add_array){print $0, add_array[main_file_index]}
	   else{print $0, add_missing_value}
      }


    
#     dynamic regex:    
#     BEGIN { digits_regexp = "[[:digit:]]+" }
#      $0 ~ digits_regexp    { print }

# 
#     if ($0 ~ /^[ATGCNatgcn][ATGCNatgcn][ATGCNatgcn]+/) {seq_original = seq_original $0 ; print "sequence data",NR}				#sequence data
#     if ($0 == "+") {print $0 ; print "fastq qual defline"}							                # fastq qual defline
#     if (($0 ~ /^@[A-Za-z0-9]+/) && (trace_number >  0)) {trace_number ++ ; process_trace(seq_original); seq_original = "";  curr_defline = $0 ; print "new sequence", NR}   #fastq sequence defline    
#     if (($0 ~ /^@[A-Za-z0-9]+/) && (trace_number == 0)) {curr_defline = $0 ; trace_number ++ ; print "first fastq sequence defline",NR}   #first fastq sequence defline

# LC_ALL=C explanation: needed for this regex to work: /^[!-\]]/     awk4 does not need LC_ALL=C
#  http://www.catonmat.net/blog/ten-awk-tips-tricks-and-pitfalls/
# http://unix.stackexchange.com/questions/87745/what-does-lc-all-c-do

# I will assume the following format, in order to handle fastq and fasta:
# 
# ^>              : fasta defline
# ^@[A-Za-z0-9]+  : fastq sequence defline
# ^\+$            : fastq qual defline
# ^[ATGCN]{10}    : sequence data (fasta and fastq)
# ^[!-\]]{10}     : fastq qual data  (ascii 33 - ascii 93)  CANNOT BE USED, because range inlcudes ATGC!

