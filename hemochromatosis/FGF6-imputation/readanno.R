#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
file=paste(args[1],".hg19_multianno.csv",sep="")
mPred=as.numeric(args[2]);
data<-read.table(file,head=T,sep=",",check.names = F)
RiskAAC<-function(data,mPred){
  risk1<-which(unlist(apply(data[,grep("pred|gwas",colnames(data))],1,function(x) length(grep("D|P|H|M|Name",as.character(unlist(x))))))>mPred)
  rlt=unique(c(risk1))
  return(rlt)
}
Riskloci<-RiskAAC(data,mPred)
rlt<-data[Riskloci,]
ID=paste(rlt$Chr,":",rlt$Start,sep="")
rlt<-data.frame(rlt[,1:2],rlt[,4:5],ID,rlt[,6:7],rlt[,grep("pred|gwasCatalog",colnames(rlt))])
write.table(rlt,file=paste(args[1],"mPred",mPred,"anno.txt",sep="."),col.names=F,row.names=F,quote=F,sep="\t")
write.table(ID,file=paste(args[1],"mPred",mPred,"anno.snp",sep="."),col.names=F,row.names=F,quote=F,sep="\t")


which(rlt$G)
risk1<-which(unlist(apply(xdata[,grep("pred",colnames(xdata))],1,function(x) length(grep("D|P|H|M",as.character(unlist(x))))))>mPred)
