#! /bin/bash

##
##Need installed:
##1) R (tested with version 3.3)
##2) bedtools
##

PEAKS=$1 #Bed file of peaks
EVAL=$2 #directory to save analysis in
COVEREDGENES=$RSEM ##location of genes.results files from RSEM
GENES=$4 ##bed file of genes
DNASE=$5 ## bed file of dnase
TSS=$6 ##Bed file of annotated TSS


echo "Make evaluation directory"

mkdir $EVAL

echo "get true positives"
bedtools window -w 100 -u -a $PEAKS -b $TSS > ${EVAL}/ParaClu_true_pos.bed
echo "Number TP:"
wc -l ${EVAL}/ParaClu_true_pos.bed | echo

echo "get false positives"
bedtools window -w 100 -v -a $PEAKS -b $TSS > ${EVAL}/notTP.bed
bedtools window -w 100 -u -a ${EVAL}/notTP.bed -b $GENES > ${EVAL}/ParaClu_false_pos.bed
echo "Number FP:"
wc -l ${EVAL}/ParaClu_false_pos.bed | echo


echo "get false negatives"
Rscript getCovered.R $GENE $RSEM
bedtools intersect -u -a $TSS -b ${DIR}/coverage/covered.bed > ${EVAL}/temp.bed
bedtools intersect -u -a ${EVAL}/temp.bed -b $DNASE > ${EVAL}/must_call_set.bed
bedtools window -v -w 100 -a ${EVAL}/must_call_set.bed -b $PEAKS > ${EVAL}/ParaClu_false_neg.bed
echo "Number FN:"
wc -l ${EVAL}/ParaClu_false_neg.bed | echo


