library(useful)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE)

name=args[1]

dat=read.table(paste(name,"/",name,".comb.bed",sep=""),sep="\t")


print(head(dat))

dat["Length"]=dat[,9]-dat[,8]

dat["Dist"]=dat[,2]-dat[,8]
dat[dat[,12]=="-","Dist"]=dat[dat[,12]=="-",9]-dat[dat[,12]=="-",3]

dat=dat[dat$Length>500,]

#dat=dat[order(dat$Dist),]

#dat=dat[!duplicated(dat[,2]),]
print(head(dat))

dat["Percent"]=dat[,"Dist"]/dat[,"Length"]

dat["Num"]=as.integer(dat[,"Percent"]*100)


print(head(dat))

dat=dat[,c("V5","V10","Num")]

colnames(dat)=c("Count","Gene","Num")

print(head(dat))

tab <- dat %>% group_by(Gene,Num) %>% summarise(Count=sum(Count))

tab=data.frame(tab)

head(tab)

tab <- tab %>% spread(Num,Count,fill=0)

tab=data.frame(tab)

corner(tab)

tab=tab[,2:101]

tab["total"]=apply(tab,1,sum)
tab=tab[order(tab[,"total"],decreasing=T),]
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

