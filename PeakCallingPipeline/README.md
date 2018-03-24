#Peak Calling Pipeline

Code for the pipeline that takes in bamfiles from a given five prime method (mapped to the genome with STAR) and outputs peaks.

#Pre Pipeline

Before using the pipeline to call peaks, users must map the 5' RNA-Seq reads to the genome of interest. We performed this using star, and example script is included in this directory (Star.sh). The script takes 4 arguements: name of output, input file (a fastq file), a path to the genome being mapped to (in STAR index format), and the number of bases to soft clip at the 5' end.

For the soft clipping, this number should be at least as big as the number of bases used for barcoding, umi, etc at the 3' end of a given method. For template switching methods, an extra base should be added to account for the extra G often found in these methods (this is particularly important if one wants to use our modified version of Capfilter to filter peaks).

#Running the Pipeline

The actual pipeline can be run with the SmallPipeline.sh command. This takes in the bam output by star and runs Paraclu on it. It takes in numerous arguments:
1) The input bam file
2) The name of the output directory (will create this directory if it does not exist).
3) The number of reads to downsample to (most of the analysis in our paper used 20000000)
4) - 7) Numerous filtering parameters (see our manuscript for details). More spcifically:
4) Min Value
5) Min density rise
6) Min pos with data
7) Min Sum

The results are output into the directory specified by argument 2. The output directory contains the folloding files:
1) downsampled.bam: The downsampled Bam 
2) paraclu.bed: The final filtered bed  
3) paraclu.ctss: The file of all ctss (the start sites of all reads in downsampled.bed, used for peak calling) 
4) - 6) paraclu.raw.bed  paraclu.temp  rlp_processed.txt: Intermediate files



#After Running the Pipeline

The files output by this method can be used for various things. The bam and bed can be used to run Capfilter (located in ../Capfilter) to clean up peak calls for strand invasion based methods. The bam can be run through RSEM to get a list of covered genes (used for finding false negaitves). This output list of covered genes can then be combined with the methods in ../AnalysisingMethods to calculate false postives, false negatives, etc. The BAM file can be used to explore 5' to 3' bias using the methods in ../BAM_Metrics. Finally, the output ctss and bam file can be used by the methods in ../eRNA_Pipeline to try to find putative eRNA. 


