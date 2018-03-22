############################################
# File Name : downsample_unique
#
# Purpose : Downsample a bam file keeping only a specified number of unique reads.
#
# Creation Date : 27-02-2016
#
# Last Modified : Mon Feb 29 12:29:00 2016
#
# Created By : ahaber
#
# Usage : 
############################################

set -e # quit on error

paired=false
output_bam=''
input_bam=''
n_reads=0
rlp="reads_list"

while getopts 'po:i:n:r:' flag; do
  case "${flag}" in
    p) paired=true ;;
    o) output_bam="${OPTARG}" ;;
	i) input_bam="${OPTARG}" ;;
	n) n_reads="${OPTARG}" ;;
	r) rlp="${OPTARG}" ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done

echo "Downsample settings"
echo "-------------------------------------"
echo "Input bam: "$input_bam
echo "N reads: "$n_reads
echo "Paired-end: "$paired
echo "Output .bam:" $output_bam
echo "Reads list prefix: " $rlp
echo "-------------------------------------"

[ $# -eq 0 ] && { echo "Usage: $0 -i reads.bam -n n_reads -p -o output.bam"; exit 1; }

echo "Keeping only QC-passed reads that are primary alignments"

# Explaining the filtering. 

# Star manual:
# For multi-mappers, all alignments except one are marked with 0x100 (secondary alignment) in
# the FLAG (column 2 of the SAM). The unmarked alignment is either the best one (i.e. highest
# scoring), or is randomly selected from the alignments of equal quality.

# -F 512 flag keeps anything that does not have the 'fail vendor QC' sam flag bit set. 
# -F 256 flag keeps anything that does not have the 'not primary alignment' sam flag bit set. 
# for both, use the sum, 768.

all_suf="_all.txt"
final_suf="_processed.txt"
if [ "$paired" = true ]; then
    echo "----------------"
    echo "PAIRED is TRUE"
    echo "----------------"
    echo "Selecting $n_reads unique, QC-passed reads from $input_bam [Linux core-utils]"
	samtools view -h -F 768 $input_bam | cut -f1 | sort | uniq  > $rlp$all_suf


    all_paired_suf="_paired_all.txt"
    filt_paired_suf="_paired_filtered.txt"
    p1_suf="_p1.txt"
    p2_suf="_p2.txt"

	echo "Collecting paired reads [Linux core-utils]"
	awk '{ printf("%s",$0); n++; if(n%2==0) { printf("\n");} else { printf("\t\t");} }' \
	    $rlp$all_suf  > $rlp$all_paired_suf

	echo "Picking pairs"
	N_PAIRS=0
	echo $N_PAIRS
	n_reads_int=${n_reads%.*}
	echo "hi"
	echo $n_reads_int
	N_PAIRS=$(($n_reads_int/2))
	echo "Selecting $N_PAIRS pairs"

	shuf -n $N_PAIRS $rlp$all_paired_suf > $rlp$filt_paired_suf

	awk '{print $1}' $rlp$filt_paired_suf > $rlp$p1_suf
	awk '{print $2}' $rlp$filt_paired_suf > $rlp$p2_suf
	cat $rlp$p1_suf $rlp$p2_suf > $rlp$final_suf
else
   	echo "PAIRED is FALSE"
   	echo "Selecting $n_reads unique, QC-passed reads from $input_bam [Linux core-utils]"
	samtools view -h -F 768 $input_bam | cut -f1 | sort | uniq | shuf -n $n_reads > $rlp$final_suf
fi


echo "Writing out these reads to --> " $output_bam [Picard tools]
#java -jar -Xmx8g /seq/software/picard-public/current/picard.jar FilterSamReads I=$input_bam O=$output_bam \
java -jar -Xmx8g /seq/software/picard-public/2.16.0/picard.jar FilterSamReads I=$input_bam O=$output_bam \
    READ_LIST_FILE=$rlp$final_suf  FILTER=includeReadList WRITE_READS_FILES="false"



