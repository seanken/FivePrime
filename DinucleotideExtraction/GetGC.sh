
CTSSFile=$1
saveName=$2
TCFile=$3
##Amount CTSS off by due to soft mapping
offBy=$4


echo $CTSSFile
echo $saveName

echo "Clean Data!"


TCGene=${saveName}.TC.genes.bed
GENES=/home/unix/ahaber/ref/bed/genes_UCSC_hg19.allIds.bed
GENES="../data/Gene.bed"


bedtools window -w 100 -a $TCfile -b $GENES -u > $TCGene

awk '($5) >= 10' $TCGene > ${saveName}.temp.bed
mv ${saveName}.temp.bed $TCGene

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
