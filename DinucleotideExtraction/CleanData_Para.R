##Make CTSS file into nice bed
args = commandArgs(trailingOnly=TRUE)

ctssfile=args[1]
outfile=args[2]
offBy=as.integer(args[3])
print(ctssfile)
ctss=read.table(ctssfile)
print(dim(ctss))
print(head(ctss))
ctss[ctss[,6]=="+",2]=ctss[ctss[,6]=="+",2]-offBy
ctss[ctss[,6]=="-",2]=ctss[ctss[,6]=="-",2]+offBy

ctss["start"]=ctss[,2]-50
ctss["end"]=ctss[,2]+51
ctss=ctss[ctss[,"start"]>0,]
ctss["name"]="NoName"

#ctss=ctss[ctss[,5]>10,]
#ctss=ctss[,c(1,5,6,7,4,2)]
ctss=ctss[,c(1,7,8,4,5,6)]
ctss=ctss[ctss[,1] %in% c(1:22,"X","Y"),]

write.table(ctss,outfile,sep="\t",quote=F,col.names=F,row.names=F)

