#!/bin/sh
#
# Copyright 2012 K.K.DNAFORM
# This file is part of idr_paraclu program.
# Idr_paraclu is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, any later version.
#
# Idr_paraclu is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar. If not, see <http://www.gnu.org/licenses/>.
#
# Make CTSS file
#
# Usage :
#		make_ctss3.sh <input>
# Parameter :  
#		input		input file
# Date :
#		Creating	2011/08/15
#
##########################################################
function usage()
{
 cat <<EOF
usage: $0 -i BAMFILE
EOF
 exit 1;
}



var=$1;
out=$2;
foo=$var
file=${foo##*/}
base=${file%%.*}
#echo $base

TMPFILE="/tmp/$(basename $0).$RANDOM.txt"
   samtools view  -F 4 -u -q 10 -b $var |  bamToBed -i stdin > $TMPFILE
   cat  ${TMPFILE} \
| awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$5}}' \
| sort -k1,1 -k2,2n \
| groupBy -i stdin -g 1,2 -c 3 -o count \
| awk -v x="$base" 'BEGIN{OFS="\t"}{print $1,"+",$2,$3}' > $out

cat  ${TMPFILE} \
| awk 'BEGIN{OFS="\t"}{if($6=="-"){print $1,$3,$5}}' \
| sort -k1,1 -k2,2n \
| groupBy -i stdin -g 1,2 -c 3 -o count \
| awk -v x="$base" 'BEGIN{OFS="\t"}{print $1,"-",$2-1,$3}' >> $out

rm $TMPFILE
   


