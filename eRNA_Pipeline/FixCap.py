import numpy as np
import scipy as sp
import pandas as pd
import sys

capfilename=sys.argv[1]
filename=sys.argv[2]
saveFile=sys.argv[3]

print(saveFile)

cap=pd.read_csv(capfilename,"\t",names=["Chr","start","end","Cluster","Count","strand","Number_CTSS"],low_memory=False)
no_cap=pd.read_csv(filename,"\t",names=["Chr","start","end","Cluster","Count","strand","Number_CTSS"],low_memory=False)

no_cap.index=no_cap["Cluster"]
cap.index=cap["Cluster"]

print(cap.head())
print(no_cap.head())

comb=no_cap.join(cap,"Cluster",how="right",rsuffix="none")
names=comb.columns.values[[0,1,2,3,4,5,6]]
print(comb.head())
comb=comb[names]
print(comb.head())
comb.to_csv(saveFile,index=False,header=False,sep="\t")

