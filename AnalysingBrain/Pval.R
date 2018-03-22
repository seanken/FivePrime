library(reshape)



cleanUp<-function(dat,nam,normed=F)
{
tab=cast(dat[dat[,"Number_Peaks_in_Gene"]>1,],formula=Gene~Position_in_Gene,value=nam,fun.aggregate=function(x){l=sum(as.numeric(x));if(is.na(l)){return(0)}else{return(l)}})

rownames(tab)<-tab[,"Gene"]
tab=tab[2:dim(tab)[2]]

return(tab);


}


meanInd<-function(lst)
{
n=length(lst)
tot=0;
sm=0;
for(i in 1:n){sm=sm+i*lst[i];tot=tot+lst[i]}
if(tot==0){return(NA)}
return(sm/tot)
}


getPval<-function(dat,samp1_name,samp2_name,exactNum=-1,getGeneMeans=F,scaled=-1,normed=F,getMeans=F,scor_test="mean")
{
if(exactNum>0)
{
dat=dat[dat[,"Number_Peaks_in_Gene"]==exactNum,]
}

samp1<-cleanUp(dat,samp1_name,normed)
samp2<-cleanUp(dat,samp2_name,normed)

score=0


library(useful)
print(corner(samp1))
n=dim(samp1)[1]
diffMean=c()
print(meanInd(samp1[1,]))
diffMean=apply(samp1,1,meanInd)-apply(samp2,1,meanInd)
scal=rep(TRUE,dim(samp1)[1])

if(scaled>0)
{
sm=apply(samp1,1,sum)
sm2=apply(samp2,1,sum)
scal=(sm>scaled &sm2>scaled)
}

if(scor_test=="mean")
{
x=apply(samp1,1,meanInd)
y=apply(samp2,1,meanInd)
}
else{
x=apply(samp1,1,function(x){return(sum(x>0))})
y=apply(samp2,1,function(x){return(sum(x>0))})
}

print(head(x))
I=scal&!is.na(x)&!is.na(y)
x=x[I]
y=y[I]
geneNames=rownames(samp1)[I]
geneMeans=cbind(x,y)

print(dim(geneMeans))
geneMeans=data.frame(geneMeans)
rownames(geneMeans)=geneNames
if(normed)
{
print(dim(dat))
val=cast(dat,formula=~Gene,value="Number_Peaks_in_Gene",fun.aggregate=mean)
nams=names(val)
val=data.frame(val)
colnames(val)=nams
#print(corner(val))
print(length(geneNames))
print(setdiff(geneNames,names(val)))
val=val[1,geneNames]
val=as.numeric(val[1,])
print(head(geneMeans))
geneMeans["temp"]=val
geneMeans["x"]=geneMeans["x"]/geneMeans["temp"]
geneMeans["y"]=geneMeans["y"]/geneMeans["temp"]
geneMeans=geneMeans[c(1,2)]
print(head(geneMeans))
x=geneMeans[,1]
y=geneMeans[,2]
}

mn=(mean((x-y)))
SE=sd((x-y))/sqrt(length(x))
SD=sd((x-y))
print(x[1:10])
print(y[1:10])
if(getGeneMeans)
{
return(geneMeans)
}
res=wilcox.test(x,y,paired=T,alternative="less")
pv=res$p.value
print(res$estimate)
if(getMeans){pv=c(pv,mn,SD,SE)}
return(pv)


}


testPval<-function(dat,samp1_name,samp2_name,numRep,scaled=-1)
{
#dat=read.table("combined.txt",sep="\t",header=T)
pv=getPval(dat,samp1_name,samp2_name,numRep=numRep,scaled=scaled)
return(pv)


}

if(!interactive())
{
dat=read.table("combined.txt",sep="\t",header=T)
lst=colnames(dat)[7:14]
for(i in lst)
{
for(j in lst)
{
numRep=10000
input=commandArgs(TRUE)
samp1_name=i#nput[1]#strsplit(input,split=" ")[[1]][1]
samp2_name=j#input[2]#strsplit(input,split=" ")[[1]][2]
print("got args!")
print(input)
print(samp1_name)
print(samp2_name)
if(i!=j)
{
pv=testPval(dat,samp1_name,samp2_name,numRep=numRep)
save_string=paste(samp1_name,samp2_name,toString(pv),sep=" ")
print(save_string)
write(save_string,"pvals.txt",append=T)
}
}}
}
