DIR=$1
NAME=$2_fin
CTSS=$DIR/ParaClu/ctss_tmp
GENES=$DIR/coverage/covered_genes.bed

mkdir $NAME

echo $NAME

head $CTSS

echo Make Bed!
Rscript makeCTSSBed.R $NAME $CTSS

head $GENES

echo Intersect!
bedtools intersect -s -a $NAME/${NAME}.ctss.bed -b $GENES -wa -wb > $NAME/${NAME}.comb.bed

echo Make Table!
Rscript makeTable.R $NAME


