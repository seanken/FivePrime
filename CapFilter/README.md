
Capfilter is a method devloped by the Megraw lab at Oregon state. In order to run our pipeline on multiple methods (not just NanoCAGE and CAGE) we slightly modified the code, though the vast majority of the work on the code was performed by the Megraw lab (see LICENSE).

The modified version of Capfilter we include here can be run as:

perl CapFilter.pl peak --bed --cutoff $cut --out $OutputFile --soft $Soft --barcodelen $Barlen $BedName $BamName $GenomeName

Here, cutoff is the percent cutoff used by the Capfilter (we suggested setting cut equal to 20). OutputFile is the name of the output file used. BedName is the name of the bed file output by the peak calling algorithm, BamName is the name of the bam file used to produce the bed file, and GenomeName is the location of the fasta genome file (here located in ../data/hg19.fa). Finally there are two other parameters (added for our purposes): Soft, which is the number of bases soft mapped with Star, and Barlen, the number of bases at the 5' end of the reads not mapped to the genome (used for UMI, barcodes, etc). 


