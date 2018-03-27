import numpy as np
import scipy as sp
import pandas as pd
import sys

##For calculating peak width
def GetWidth(x):
	width=max(x)-min(x)
	return(width)

print("Ready CTSS")

##load CTSS

CTSSFile=sys.argv[1]
savefile=sys.argv[2]

print(CTSSFile)
print(savefile)

ctss=pd.read_csv(CTSSFile,"\t",names=["Chr","strand","Location","Count"],low_memory=False)

##remove those with low counts
#ctss=ctss[ctss["Count"] >0]

print("Cluster the CTSS!")

loc_prev=0
clust_cur=0
clust=["clust0" for i in ctss["Count"]]
pos=0;

for loc in ctss["Location"]:
	#if abs(loc-loc_prev)> 40:
	if abs(loc-loc_prev)> 20:
		clust_cur=clust_cur+1
	clust[pos]="clust_"+str(clust_cur)
	pos=pos+1
	loc_prev=loc


ctss["Clust"]=clust	

ctss=ctss.groupby(["Clust","Chr","strand"],as_index=False).agg({'Location':['max','min'],'Count':['sum','count']})

ctss.columns=["Cluster","Chr","strand","Count","Number_CTSS","end","start"]

ctss=ctss[["Cluster","Chr","strand","Count","Number_CTSS","start","end"]]

ctss=ctss[["Chr","start","end","Cluster","Count","strand","Number_CTSS"]]
ctss["Number_CTSS"]=ctss["end"]-ctss["start"]

ctss=ctss[ctss["Count"]>2]
ctss=ctss[ctss["Number_CTSS"]<300]


ctss=ctss[ctss["Chr"]!="MT"]
ctss["end"]=ctss["end"]+1

print("Save the file!")

ctss.to_csv(savefile,index=False,header=False,sep="\t")





