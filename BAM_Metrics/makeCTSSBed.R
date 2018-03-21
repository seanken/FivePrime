args = commandArgs(trailingOnly=TRUE)

name=args[1]
ctss=args[2]

print(ctss)
print(name)

dat=read.table(ctss,sep="\t")


colnames(dat)[1:4]=c("Chr","Strad","Start","Count")

dat["end"]=dat[,"Start"]+1


dat["Name"]="None"

dat=dat[,c(1,3,5,6,4,2)]

print("Save!")
write.table(dat,paste(name,"/",name,".ctss.bed",sep=""),sep="\t",row.names=F,col.names=F,quote=F)

