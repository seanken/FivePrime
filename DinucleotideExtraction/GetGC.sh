
CTSSFile=$1 ##The location of the CTSS file, see that standard peak calling pipeline for how to calculate
saveName=$2 ##Self explanatory, will be used as a prefix
TCFile=$3 ##TC peaks, see 	
offBy=$4 ##The soft clipping number minus the barcod length (see description for Capfilter for more details)


echo $CTSSFile
echo $saveName

echo "Clean Data!"


TCGene=${saveName}.TC.genes.bed
GENES=/home/unix/ahaber/ref/bed/genes_UCSC_hg19.allIds.bed
GENES="../data/Gene.bed"

echo $GENES

bedtools window -w 100 -a $TCFile -b $GENES -u > $TCGene

echo "Remove low coverage"

awk '($5) >= 10' $TCGene > ${saveName}.temp.bed
mv ${saveName}.temp.bed $TCGene

echo "Make CTSS bed"

CTSSBed=${saveName}.ctss.bed
Rscript makeCTSSBed.R $CTSSFile $CTSSBed
bedtools intersect -b $CTSSBed -a $TCGene -wa -wb -s > ${saveName}.temp.bed
#bedtools intersect -b $CTSSBed -a $TCfile -wa -wb -s > output/${saveName}.temp.bed
#mv output/${saveName}.temp.bed $CTSSBed
Rscript GetMode.R $CTSSBed ${saveName}.temp.bed

Rscript CleanData_Para.R $CTSSBed ${saveName}.bed $offBy

echo "Get Sequence!"
bedtools getfasta -fi ~/FivePrime/hg19/hg19.fa -bed ${saveName}.bed -s > ${saveName}.fa

echo "Make table!"
Rscript makeTable.R ${saveName}.fa ${saveName}.table.txt $saveName
