library(dplyr)
library(tidyr)
args = commandArgs(trailingOnly=TRUE)

saveFile=args[1]
ctssFile=args[2]


dat=read.table(ctssFile)

tab <- dat %>% group_by(V4,V1,V6) %>% summarise(start=V9[which.max(V12)],end=V10[which.max(V12)],val=max(V12))
tab=data.frame(tab)
tab=tab[,c(2,4,5,1,6,3)]

write.table(tab,saveFile,sep="\t",row.names=F,col.names=F,quote=F)
