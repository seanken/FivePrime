#!/bin/bash


INPUT=$1
OUTPUT_DIR=$2
NUM=$3
MinVal=$4
MIN_WIDTH=$5
MAX_WIDTH=$6
min_density_rise=$7
min_pos_with_data=$8
min_sum=$9



DS=${OUTPUT_DIR}/downsampled.bam
CTSS=${OUTPUT_DIR}/paraclu.ctss

echo "Make directory"
mkdir $OUTPUT_DIR
echo $OUTPUT_DIR

echo "Dowsample to $NUM"
./downsample_unique.sh -o ${OUTPUT_DIR}/downsampled.bam -i $INPUT -n $NUM -r ${OUTPUT_DIR}/rlp

echo "Make CTSS"
./make_ctss3.sh $DS $CTSS

echo "Run Paraclu"
PARACLU_OUT=${OUTPUT_DIR}/paraclu
./RunParaclu.sh $CTSS $MinVal $PARACLU_OUT $MIN_WIDTH $MAX_WIDTH $min_density_rise $min_pos_with_data $min_sum



