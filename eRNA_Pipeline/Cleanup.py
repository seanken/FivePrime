import math
import sys
import pandas as pd
import numpy as np

##Collapse eRNA with overlapping TC
def Collapse(pairs):
	startTC={}
	toRemove=[]
	comb=pd.DataFrame(columns=pairs.columns)
	print(comb.head())
	print(pairs.head())
	pairs.index=range(pairs.shape[0])
	for i in range(0,pairs.shape[0]):
		##If neither cluster occured before
		if (not pairs["Cluster1"][i] in startTC.keys()) & (not pairs["Cluster2"][i] in startTC.keys()):
			pos=comb.shape[0]
			comb.loc[comb.shape[0]]=list(pairs.loc[i])
			startTC[pairs["Cluster1"][i]]=pos
			startTC[pairs["Cluster2"][i]]=pos
			continue;
		
		##If both clusters occured before
		if (pairs["Cluster1"][i] in startTC.keys()) & (pairs["Cluster2"][i] in startTC.keys()):
			pos1=startTC[pairs["Cluster1"][i]]
			pos2=startTC[pairs["Cluster2"][i]]
			start1=min(comb["start1"][pos1],comb["start1"][pos2])
			start2=min(comb["start2"][pos1],comb["start2"][pos2])
			end1=max(comb["end1"][pos1],comb["end1"][pos2])
			end2=max(comb["end2"][pos1],comb["end2"][pos2])
			if end1 > start2:
				continue
			comb["end1"][pos1]=end1
			comb["end2"][pos1]=end2
			comb["start1"][pos1]=start1
			comb["start2"][pos1]=start2
			comb["Count1"][pos1]=comb["Count1"][pos1]+comb["Count1"][pos2]
			comb["Count2"][pos1]=comb["Count2"][pos1]+comb["Count2"][pos2]
			for clust in startTC.keys():
				if startTC[clust]==pos2:
					startTC[clust]==pos1
			comb["Count1"][pos2]=0
			comb["Count2"][pos2]=0
			continue;

		##If one cluster occured before, but not the other
		pos=0
		if pairs["Cluster1"][i] in startTC.keys():
			pos=startTC[pairs["Cluster1"][i]]

		if pairs["Cluster2"][i] in startTC.keys():
			pos=startTC[pairs["Cluster2"][i]]

		start1=min(pairs["start1"][i],comb["start1"][pos]);
		start2=min(pairs["start2"][i],comb["start2"][pos]);
		
		end1=max(pairs["end1"][i],comb["end1"][pos]);
		end2=max(pairs["end2"][i],comb["end2"][pos]);
		
		if start2 < end1:
			continue;

		comb["end1"][pos]=end1
		comb["end2"][pos]=end2
		comb["start1"][pos]=start1
		comb["start2"][pos]=start2
		startTC[pairs["Cluster1"][i]]=pos
		startTC[pairs["Cluster2"][i]]=pos
		comb["Count1"][pos]=comb["Count1"][pos]+pairs["Count1"][i]
		comb["Count2"][pos]=comb["Count2"][pos]+pairs["Count2"][i]
	return(comb)
print("Load pairs file")
loadfile=sys.argv[1]+".paired.bed"
enhancerFile=sys.argv[1]+".enhance.bed"

cols=["Chr1","start1","end1","Cluster1","Count1","strand1","Number_CTSS1","Chr2","start2","end2","Cluster2","Count2","strand2","Number_CTSS2"]

pairs=pd.read_csv(loadfile,"\t",names=cols,low_memory=False)

print("Clean Up")
pairs=pairs[pairs["strand1"]=="-"]
print(pairs.shape)
print(pairs.head())
pairs=pairs[pairs["strand2"]=="+"]
print(pairs.shape)
pairs=pairs[pairs["start1"] < pairs["start2"]]
print(pairs.shape)
print("Collapse!")
pairs=Collapse(pairs)

print("Remove Promoters!")
pairs["totalExpr"]=pairs["Count1"]+pairs["Count2"]

pairs=pairs[pairs["Count1"] > 1]
pairs=pairs[pairs["Count2"] > 1]

##quantity used in paper to remove promoter like pairs
pairs["D"]=abs((pairs["Count1"]-pairs["Count2"])/pairs["totalExpr"])
pairs=pairs[pairs["D"]<.8]

pairs=pairs.sort_values(by='totalExpr', ascending=0)


print(pairs.shape)
print(pairs.head())
pairs.to_csv(loadfile,index=False,header=False,sep="\t")

pairs["Chr"]=pairs["Chr1"]
pairs["mid"]=(pairs["end1"]+pairs["start2"])/2.0#np.floor((pairs["start1"]+pairs["end1"])/2.0)
pairs["start"]=pairs["mid"]-200
pairs["end"]=pairs["mid"]+200
pairs=pairs[["Chr","start","end"]]
print(pairs.head())
pairs.to_csv(enhancerFile,index=False,header=False,sep="\t")
