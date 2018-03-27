Code for analysising peaks in the methods comparision part.

Require:
1) R (tested with 3.3)
2) bedtools

The code takes in the results of peak calling and other data, spits out the results in terms of TP, FP, and FN. The main method, getScore.sh, takes in 6 arguments:
1) Peak file (a bed file)
2) The output directory
3) pointer to the output of RSEM (the genes.results file)
4) The bed file of genes
5) The bed file of Dnase peaks
6) The TSS peak files

For example, for a peak file (peak.bed) and an RSEM file (rsem.genes.results), can run (using the files in the sources folder)

./getScore.sh peak.bed output rsem.genes.results sources/Genes.bed sources/dnase.bed sources/TSS.bed

The output is put into the output directory, in the form of numerous bed files
1) ParaClu_false_neg.bed: The false negatives
2) ParaClu_false_pos.bed: The false positives
3) ParaClu_true_pos.bed: The true positives
4) Other intermediate bed files


