
args = commandArgs(trailingOnly=TRUE)
geneFile=args[1]
direct=args[2]

print(geneFile)
print(direct)

print("Get covered genes from RSEM output")
RSEM_Output=paste(direct,"coverage","all_genes.bed.genes.results",sep="/")
print(RSEM_Output)
rsem=read.table(RSEM_Output,header=T)
rownames(rsem)=rsem[,"gene_id"]
rsem=rsem[rsem$TPM > 1,]

covered_genes=as.character(rownames(rsem))


print("Extract covered genes location")
gene_location=read.table(geneFile)
rownames(gene_location)=gene_location[,4]
gene_location=gene_location[covered_genes,]

print("Save covered genes")
savefile=paste(direct,"coverage","covered.bed",sep="/")
print(savefile)
write.table(gene_location,savefile,sep="\t",col.names=F,row.names=F,quote=F)
