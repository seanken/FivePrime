args = commandArgs(trailingOnly=TRUE)
ctssfile=args[1]
ctssbed=args[2]

ctss=read.table(ctssfile)
colnames(ctss)=c("Chr","Sign","Start","Num")
ctss["End"]=ctss[,"Start"]+1
ctss["Name"]="NoName"
ctss=ctss[,c("Chr","Start","End","Name","Num","Sign")]

ctss=ctss[ctss[,1] %in% c(1:22,"X","Y"),]

write.table(ctss,ctssbed,sep="\t",row.names=F,col.names=F,quote=F)

