args = commandArgs(trailingOnly=TRUE)
fasta=args[1]
outfile=args[2]
method=args[3]

print("load data and reformat!")

lst=scan(fasta,what="")

lst_plus=lst[(grep("+",lst,fixed=T)+1)]
lst_minus=lst[(grep("(-)",lst,fixed=T)+1)]

#temp=lst_plus
#lst_plus=c()
#for(l in temp){lst_plus<-c(lst_plus,substr(l,1,100))}


#temp=lst_minus
#lst_minus=c()
#for(l in temp){lst_minus<-c(lst_minus,substr(l,2,101))}

lst=c(lst_plus,lst_minus)
#lst=lst_minus
lst=toupper(lst)
print(length(lst_plus))
print(length(lst_minus))

tab=matrix(0,ncol=4,nrow=101)
tab=data.frame(tab)
colnames(tab)=c("A","T","G","C")

print("Make table!")

#lst=sample(lst,1000)

for(l in lst){for(i in 1:101){cur=substr(l,i,i);if(cur %in% colnames(tab)){tab[i,cur]=tab[i,cur]+1} }}
tab["method"]=method
tab["row"]=1:dim(tab)[1]

print("save")
write.table(tab,outfile,sep=",",col.names=F,row.names=F,quote=F)

