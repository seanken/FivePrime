#! /bin/sh

# This script simplifies the output of paraclu.  It omits clusters
# that are too long, are singletons, or have too low "strength".
# Then, it removes clusters that are contained in larger clusters.  It
# writes the remaining clusters in BED format.

# Column 4 in the output is: (number of sites):(number of tags).

# Column 5 in the output is the cluster strength.  Using the notation
# of PMID:18032727, this is: log2[(max d) / (min d)].  Infinite
# strength is written as 1000.

if [ $# -lt 1 ]
then
    echo "Usage: $(basename $0) minClusterStrength maxClusterLength file" 1>&2
    exit 2
fi

# Convert to BED format.
cat $1 |
awk '
BEGIN { OFS="\t"; OFMT = "%.3g" }
!/^#/ {
  print $1, $3, $4+1, $5, $6, $2, $7, $8
}
'  

