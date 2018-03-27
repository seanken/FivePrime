#! /bin/bash

##
##Need installed:
##1) R (tested with version 3.3)
##2) bedtools
##

PEAKS=$1 #Bed file of peaks
EVAL=$2 #directory to save analysis in
RSEM=$3 ##location of genes.results files from RSEM
GENES=$4 ##bed file of genes
DNASE=$5 ## bed file of dnase
TSS=$6 ##Bed file of annotated TSS


echo "Make evaluation directory"

mkdir $EVAL

echo "get true positives"
bedtools window -w 100 -u -a $PEAKS -b $TSS > ${EVAL}/ParaClu_true_pos.bed
echo "Number TP:"
wc -l ${EVAL}/ParaClu_true_pos.bed 
echo ""

echo "get false positives"
bedtools window -w 100 -v -a $PEAKS -b $TSS > ${EVAL}/notTP.bed
bedtools window -w 100 -u -a ${EVAL}/notTP.bed -b $GENES > ${EVAL}/ParaClu_false_pos.bed
echo "Number FP:"
wc -l ${EVAL}/ParaClu_false_pos.bed 
echo ""

echo "get false negatives"
#Rscript getCovered.R $GENES $RSEM $EVAL
python temp.py $GENES $RSEM $EVAL
#cp $RSEM ${EVAL}/covered.bed
echo Get covered
#head $TSS
#head ${EVAL}/covered.bed
wc -l ${EVAL}/rsem.unsorted.bed
echo "sort"
bedtools sort -i ${EVAL}/rsem.unsorted.bed > ${EVAL}/covered.sorted.bed
echo "merge"
bedtools merge -i ${EVAL}/covered.sorted.bed -c 4,5,6,11,12,13 -o first,first,first,max,max,max > ${EVAL}/temp.bed
mv ${EVAL}/temp.bed ${EVAL}/covered.bed

awk '{if($12 > 1){print $0}}' $EVAL/covered.bed  > ${EVEL}/temp.bed
mv ${EVAL}/temp.bed $EVAL/covered.bed

bedtools intersect -a $TSS -b ${EVAL}/covered.bed -u > ${EVAL}/temp.bed
cp ${EVAL}/temp.bed ${EVAL}/must_call_set.bed
echo Get DNASE
bedtools window -w 100 -u -a ${EVAL}/temp.bed -b $DNASE > ${EVAL}/must_call_set.bed
echo Get FN!
bedtools window -w 100 -v -a ${EVAL}/must_call_set.bed -b $PEAKS  > ${EVAL}/ParaClu_false_neg.bed
echo "Number FN:"
wc -l ${EVAL}/ParaClu_false_neg.bed


