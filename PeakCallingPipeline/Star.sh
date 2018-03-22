#! /bin/bash

##Require Star to be installed and on your path

##Pass for arguments:
##1) Name of output
##2) Fastq file
##3) Name of STAR reference 
##4) Number of bases to clip as 5' end
genome=$3
NumClip=$4
STAR --runThreadN 4 --genomeDir $genome --outFileNamePrefix $1 --readFilesIn $2 --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM Unsorted --clip5pNbases $NumClip  0 
