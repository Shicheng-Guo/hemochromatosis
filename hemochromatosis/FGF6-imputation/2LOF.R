#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ssum<-function(vector){
options(warn=-1)
rlt<-c()
  for(i in seq(1,length(vector),by=2)){
    tmp<-(vector[i]>=1)+(vector[i+1]>=1)
    rlt<-c(rlt,tmp)
  }
rlt
}

cssum<-function(vector){
options(warn=-1)
rlt<-c()
  for(i in seq(1,length(vector),by=2)){
    tmp<-(vector[i])*(vector[i+1])
    rlt<-c(rlt,tmp)
  }
rlt
}

restru<-function(data){
  rlt<-c()
  for(i in 1:ncol(data)){
    yy<-lapply(strsplit(as.character(data[,i]),"[|]"),function(x) t(x))
    zz<-matrix(unlist(yy),ncol=2,byrow=T)
    rlt<-cbind(rlt,zz)
  }
  rlt
}

vcfile=as.character(args[1])  
samfile=as.character(args[2]) 
#vcfile="chr12.update.vcf"
#samfile="/gpfs/home/guosa/hpc/hemochromatosis/FGF6-imputation/PheTyp7_Iron_C1.phen"

#samfile="/gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_ENA_rev2_SampleIDs.Michigen.txt"
#vcfile="chr22.update.vcf"
#samfile="/gpfs/home/guosa/hpc/project/pmrp/phen/IBDCH_Phetyp6_SSc_rev2_SampleIDs.Michigen.txt"

print(vcfile);
print(samfile);
chr=unlist(strsplit(vcfile,"[.]"))[1]
fileout=unlist(strsplit(samfile,"/"))
samfilehead=fileout[length(fileout)]
output<-paste(samfilehead,vcfile,"pvalue.txt",sep=".")

head<-readLines(vcfile)
data<-read.table(vcfile,head=F,sep="\t",check.names=F)
head<-head[-(grep("#CHROM",head)+1):-(length(head))]
vcfnames<-unlist(strsplit(head[length(head)],"\t"))
names(data)<-vcfnames

phen<-read.table(samfile,head=T)
samRank1<-as.character(colnames(data))
samRank2<-as.character(paste(phen[,2],phen[,2],sep="_"))
newphen<-phen[match(samRank1,samRank2),]
case<-which(newphen[,6]=="2")
con<-which(newphen[,6]=="1")
rlt<-c()

GeneNameList<-lapply(unique(data$INFO),function(x) unlist(strsplit(as.character(x),split="TAG3=|;"))[2])
GeneList<-unique(unlist(GeneNameList))

add2LOF<-function(sxx,syy){
  A<-sum(ssum(sxx)>1)
  B<-sum(ssum(sxx)<=1)
  C<-sum(ssum(syy)>1)
  D<-sum(ssum(syy)<=1)
  z<-matrix(c(A,B,C,D),2,2)
  try(test<-fisher.test(z, alternative = "greater"))
  P=test$p.value
  OR=((A*D)+0.5)/((B*C)+0.5)
  return(c(A,B,C,D,OR,P))
}

mult2LOF<-function(sxx,syy){
  Av<-cssum(sxx)
  Dv<-cssum(syy)
  try(test<-t.test(Av,Dv))
  P=test$p.value
  beta=test$statistic
  return(c(beta,P))
}


for(i in GeneList){
  sel<-which(unlist(lapply(data$INFO,function(x) sum(unlist(strsplit(as.character(x),split="TAG3=|;"))%in% i)>0)))
  ca<-data[sel,case]
  co<-data[sel,con]
  head(ca)
  head(co)
  xx<-data.matrix(restru(ca))
  yy<-data.matrix(restru(co))
  class(xx)="numeric"
  class(yy)="numeric"
  rm<-which(rowMeans(yy)>=0.3)
  if(length(rm)>0){
  xx<-xx[-rm,]
  yy<-yy[-rm,]  
  print(i)
  }
  xx<-matrix(xx,ncol=2*ncol(ca))
  yy<-matrix(yy,ncol=2*ncol(co))
  if(nrow(xx)>0){
  sxx<-colSums(xx)
  syy<-colSums(yy)

  rlt1<-add2LOF(sxx,syy)
  rlt2<-mult2LOF(sxx,syy)
  temp<-data.frame(i,chr,nrow(yy),A=rlt1[1],B=rlt1[2],C=rlt1[3],D=rlt1[4],OR=rlt1[5],P1=rlt1[6],Beta=rlt2[1],P2=rlt2[2])
  print(temp)
  rlt<-rbind(rlt,temp)
  }
  if(PRINT=T){
	try=data[sel,10:ncol(data)]
	xx<-data.matrix(restru(try))
	class(xx)="numeric"
	sxx<-colSums(xx)
	vxx<-c()
	for(i in seq(1,length(sxx),2)){
	vxx<-c(vxx,sxx[i]+sxx[i+1])
	}
	names(vxx)<-colnames(data)[10:ncol(data)]
	phen2<-data.frame(phen[match(names(vxx),samRank2),],Gene=vxx)
	table(phen2$phen,phen2$Gene)
	newdata=subset(phen2,phen>0)
	table(newdata$phen,newdata$Gene)
	write.table(newdata,file=paste(i,samfile,"raw.txt",sep="."),sep="\t",quote=F,row.names=F,col.names=T)
	table(newdata2$phen,newdata2$Gene)
	glmrlt<-glm(phen~Gene+sex,newdata2,family = gaussian)
	}	
}
colnames(rlt)<-c("GeneSymbol","chr","NSNP","Case+","Case-","Normal+","Normal-","OR","P","Beta","P")
write.table(rlt,file=output,sep="\t",quote=F,col.names=T,row.names = F)




