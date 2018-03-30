#!/bin/bash

##Runs the ctss generated earlier through paraclu, then filters the results
##Takes in CTSS location, MinVal, and OUTPUT location for paraclu (outputs files $OUTPUT.temp, $OUTPUT.raw.bed, and $OUTPUT.bed
INPUT=$1
MinVal=$2
OUTPUT=$3
MIN_WIDTH=3
MAX_WIDTH=300
min_density_rise=$4
min_pos_with_data=$5
min_sum=$6

echo "Run Paraclu!"
paraclu/paraclu $MinVal $INPUT > ${OUTPUT}.temp

echo "Convert output to bed format"
paraclu-to-bed.sh ${OUTPUT}.temp > ${OUTPUT}.raw.bed
head ${OUTPUT}.raw.bed
wc -l ${OUTPUT}.raw.bed
echo "Filter!!"
echo "Filter by length"
cat ${OUTPUT}.raw.bed | awk '{if($3-$2 < 300){print $0}}' | awk '{if($3-$2>3){print $0}}' > ${OUTPUT}.temp.bed
mv ${OUTPUT}.temp.bed ${OUTPUT}.raw.bed
head ${OUTPUT}.raw.bed
wc -l ${OUTPUT}.raw.bed
echo "Sort!"
bedtools sort -i ${OUTPUT}.raw.bed > ${OUTPUT}.temp.bed
mv ${OUTPUT}.temp.bed ${OUTPUT}.raw.bed
head ${OUTPUT}.raw.bed
wc -l ${OUTPUT}.raw.bed
echo "Merge!"
bedtools merge -s -c 4,5,6,7,8  -o max,max,distinct,min,max  -i ${OUTPUT}.raw.bed > ${OUTPUT}.temp.bed
mv ${OUTPUT}.temp.bed ${OUTPUT}.raw.bed
head ${OUTPUT}.raw.bed
wc -l ${OUTPUT}.raw.bed
echo "Filter!"
echo $OUTPUT $MIN_WIDTH $MAX_WIDTH $min_density_rise $min_pos_with_data $min_sum
python Filters.py $OUTPUT $MIN_WIDTH $MAX_WIDTH $min_density_rise $min_pos_with_data $min_sum


echo "Paraclu done!"

