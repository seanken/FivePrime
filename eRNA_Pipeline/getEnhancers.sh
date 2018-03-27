###
##Need to have installed:
##1) samtools (tested with 1.3.1)
##2) python (tested with 2.7.9)
##3) bedtools
##############################

CTSSFile=$1
saveName=$2
BamName=$3
Soft=$4
Barlen=$5

GENES="../data/Gene.bed"
echo "Get Peaks"
python easyCall.py $CTSSFile ${saveName}.bed

echo "Reorder"
bedtools sort -i ${saveName}.bed > ${saveName}.temp.bed
cp ${saveName}.temp.bed ${saveName}.bed


echo "Run Capfilter"
GenomeName="../data/hg19.fa"
cut=20
OutputFile=${saveName}.temp.bed
perl ../CapFilter/CapFilter.pl peak --bed --cutoff $cut --out $OutputFile --soft $Soft --barcodelen $Barlen ${saveName}.bed $BamName $GenomeName
#cp ${saveName}.temp.bed ${saveName}.mess.bed
python FixCap.py ${saveName}.temp.bed $saveName.bed $saveName.bed

cp ${saveName}.bed ${saveName}.TC.bed

echo "Get plus and minus"

bedtools window -v -w 300 -a ${saveName}.bed -b $GENES | grep -F + > ${saveName}.plus.bed
bedtools window -v -w 300 -a ${saveName}.bed -b $GENES | grep -v -F + > ${saveName}.minus.bed
bedtools intersect -v -a ${saveName}.plus.bed -b ${saveName}.minus.bed > $saveName.bed
bedtools intersect -v -b ${saveName}.plus.bed -a ${saveName}.minus.bed >> $saveName.bed
bedtools sort -i ${saveName}.bed > ${saveName}.temp.bed
mv ${saveName}.temp.bed $saveName.bed

echo "Get possible enhancer pairs"
bedtools window -w 400 -a ${saveName}.bed -b ${saveName}.bed > ${saveName}.paired.bed 
python Cleanup.py ${saveName}

sed -i 's/\.0//g' ${saveName}.enhance.bed
sed -i 's/\.5//g' ${saveName}.enhance.bed

echo "Get QC"
Tot_enhancers=$(wc -l ${saveName}.enhance.bed| awk '{print $1}')
intersect_Dnase=$(bedtools intersect -u -a ${saveName}.enhance.bed -b sources/dnase.narrow.bed | wc -l)
intersect_H3k27ac=$(bedtools intersect -u -a ${saveName}.enhance.bed -b sources/H3k27ac_narrow.bed | wc -l)
intersect_broad_H3k27ac=$(bedtools intersect -u -a ${saveName}.enhance.bed -b sources/H3k27ac.broad.bed | wc -l)
intersect_both=$(bedtools intersect -u -a ${saveName}.enhance.bed -b sources/intersect_broad.bed | wc -l)
intersect_eRNA=$(bedtools intersect -u -a ${saveName}.enhance.bed -b sources/eRNA_encode_sort.bed | wc -l)
echo Total Reads,Overlap dnase,Overlap H3k27ac, Overlap H3k27ac broad, Overlap both narrow,Overlap eRNA > ${saveName}.results.csv
echo $Tot_enhancers,$intersect_Dnase,$intersect_H3k27ac,$intersect_broad_H3k27ac,$intersect_both,$intersect_eRNA >> ${saveName}.results.csv


rm debug.txt


