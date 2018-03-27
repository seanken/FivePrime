Need to have installed:
1) samtools (tested with 1.3.1)
2) python (tested with 2.7.9)
3) bedtools


For python, you need to following packages:
1) numpy (tested with version ...)
2) pandas (tested with version ...)
3) scipy (tested with version ...)

How to use:

The pipeline takes in mapped reads from any fiveprime method and produces peaks that are putative eRNA (as well as peaks that will be used for other analysis later on). This takes 5 arguments:
1) A ctss file (output by the standard peak calling pipeline)
2) The name used to save results (will be used as a prefix for all files output)
3) The name of the bam file used/ from which the ctss have been extracted
4) The length used for soft mapping in the bam
5) The number of bp at the beginning of each reads that come from non-genomic sources (umi, barcodes, etc)


For example, let CTSS.txt be the CTSS file for our method (see the standard Peak Calling pipeline for how to calculate this), reads.bam be the read bam file, 12 and 9 be the parameters described in the capfilter description, and output/Nano be the name we want to use for saving the results. Then you simply need to run:

./getEnhancers.sh CTSS.txt output/Nano reads.bam 12 9

And the results will be output into the output directory with the prefix Nano. The important outputs are:
1) Nano.enhance.bed: A Bed file of putative enhancers
2) Nano.TC.bed: A set of 'tage clusters' (an alternative form of peak calling).
3) Nano.results.csv: A file with information about how many of the putative enhancers overlap various sources of data.
