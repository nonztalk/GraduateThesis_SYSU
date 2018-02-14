library(limma)

args<-commandArgs(trailingOnly = TRUE)
sample<-args[1]

fileInfo<-read.table("fileInfo.txt",header=F,stringsAsFactors=F)
tryCatch({
	RG<<-read.maimages(fileInfo[,1], source="agilent")
	}, error = function(err) {
		print("R, G, RB, GB should be defined")
		RG<<-read.maimages(fileInfo[,1], source="agilent", columns = list(R = "rMeanSignal", G = "gMeanSignal", Rb = "rBGUsed", Gb = "gBGUsed"))
})
RGnorm<-backgroundCorrect(RG,method="normexp",normexp.method = "rma")
RGnorm<-normalizeBetweenArrays(RGnorm,method="quantile")
data<-data.frame(RGnorm$genes$ProbeName,RGnorm$M,stringsAsFactors=F)
colnames(data)<-c("probeId",fileInfo[,2])
annot<-read.table(file=paste("/data2/ArrayReAnnotation/FastaFile/Agilent/",strsplit(sample,"_")[[1]][1],".annot.uniq.strand.txt",sep=""),header=F,stringsAsFactors=F)
exp1<-lapply(unique(annot[,2]), function (x){colMeans(data[data[,1]%in%annot[which(annot[,2]==x),1],][,-1])})
exp<-do.call(rbind,exp1)
rownames(exp)<-unique(annot[,2])
exp<-exp[exp[,1]!="NaN",]
save(exp,file=paste(sample,".exp.R",sep=""))
exp<-data.frame(GeneSymbol=rownames(exp),exp,stringsAsFactors=F)
write.table(exp,paste(sample,"_express.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")