library(useful)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

name=args[1]
cover=args[2]
dat=read.table(paste(name,"/",name,".comb.bed",sep=""),sep="\t")


print(head(dat))

dat["Length"]=dat[,"V14"]

dat["Dist"]=dat[,2]-dat[,8]
dat[dat[,12]=="-","Dist"]=dat[dat[,12]=="-",9]-dat[dat[,12]=="-",3]

dat=dat[dat$Length>500,]

#dat=dat[order(dat$Dist),]

#dat=dat[!duplicated(dat[,2]),]
print(head(dat))

dat["Percent"]=(dat[,"Dist"]+dat[,"V13"])/dat[,"Length"]

dat["Num"]=as.integer(dat[,"Percent"]*100)


print(head(dat))

dat=dat[,c("V5","V10","Num","Length")]

colnames(dat)=c("Count","Gene","Num","Length")

print(head(dat))

tab <- dat %>% group_by(Gene,Num) %>% summarise(Count=sum(Count))

tab=data.frame(tab)

head(tab)

tab <- tab %>% spread(Num,Count,fill=0)

tab=data.frame(tab)

corner(tab)
#save(tab,file="temp.Robj")
rownames(tab)=tab[,1]
cov=read.table(cover,sep="\t",header=T)
cov=cov[cov$TPM>1,]

rownames(cov)=cov[,1]
inter=intersect(rownames(tab),rownames(cov))
tab=tab[inter,]
tab["Gene"]=cov[rownames(tab),"gene_id"]
tab["TPM"]=cov[rownames(tab),"TPM"]

#tab=tab[,2:101]

tab["total"]=apply(tab[,2:101],1,sum)
tab=tab[order(tab[,"TPM"],decreasing=T),]
print(head(tab))
tab=tab[!duplicated(tab[,"Gene"]),]
print(head(tab))
tab=tab[,2:103]
#tab=tab[1:1000,]
#tab=tab[tab[,"total"]>100,]
for(col in colnames(tab)){tab[col]=tab[,col]/tab[,"total"]}

corner(tab)
mn=apply(tab,2,mean)
mn=mn[1:100]

print(sum(mn))
print(head(mn))

lst=c(name,mn)

write(lst,paste(name,"/",name,".results.txt",sep=""),sep="\n")

