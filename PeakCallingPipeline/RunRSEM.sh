#!/bin/bash


##This command takes in a bam file (with some extra information) to run RSEM on it.
##The result is used to find 'covered' genes (genes with 5' reads on them), which is used in calculating the false positive
##Inputs:
##1) Name of downsampled bam file
##2) Name used for temporary fastq file
##3) Name used for the output from RSEM
##4) The location of the index for RSEM
##5) The location of the jar file for PICARD
##
##Assumes have java installed, have the jar file for PICARD on your computer, and have RSEM installed/ on your path 
##
##

BAM=$1
FASTQ=$2
RSEM_OUT=$3
RSEM_INDEX=$4
PICARD=$5

rsem=rsem-calculate-expression

echo "Convert bam into fastq file!"
java -Xmx8g -jare $PICARD SamToFastq INPUT=$BAM OUTPUT=$FASTQ 

echo "Run RSEM!"
$rsem --no-bam-output --quiet --bowtie-chunkmbs 512 $RSEM_INDEX $RSEM_OUT
