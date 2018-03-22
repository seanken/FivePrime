#!/bin/bash

##Runs the ctss generated earlier through paraclu, then filters the results
##Takes in CTSS location, MinVal, and OUTPUT location for paraclu (outputs files $OUTPUT.temp, $OUTPUT.raw.bed, and $OUTPUT.bed
INPUT=$1
MinVal=$2
OUTPUT=$3
MIN_WIDTH=$4
MAX_WIDTH=$5
min_density_rise=$6
min_pos_with_data=$7
min_sum=$8

echo "Run Paraclu!"
paraclu/paraclu $MinVal $INPUT > ${OUTPUT}.temp

echo "Convert output to bed format"
paraclu-to-bed.sh ${OUTPUT}.temp > ${OUTPUT}.raw.bed

echo "Filter!!"
python Filters.py $OUTPUT $MIN_WIDTH $MAX_WIDTH $min_density_rise $min_pos_with_data $min_sum
echo "Paraclu done!"

