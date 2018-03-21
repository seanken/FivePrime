DIR=$1
NAME=$2_exon
CTSS=$DIR/ParaClu/ctss_tmp
COVERED=$DIR/coverage/all_genes.bed.isoforms.results
GENES=withExons.bed


mkdir $NAME

echo $NAME

head $CTSS

echo Make Bed!
Rscript makeCTSSBed.R $NAME $CTSS

head $GENES

echo Intersect!
bedtools intersect -s -a $NAME/${NAME}.ctss.bed -b $GENES -wa -wb > $NAME/${NAME}.comb.bed

echo Make Table!
Rscript makeTable_exon.R $NAME $COVERED


