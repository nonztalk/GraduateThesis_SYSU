args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

data0<-readLines(con=paste(sample,"_series_matrix.txt",sep=""))
skipPos<-grep("series_matrix_table_begin", data0)
data<-read.table(file=paste(sample,"_series_matrix.txt",sep=""),header=T,stringsAsFactors=F,sep="\t",skip=skipPos,fill=T)
data <- data[1:nrow(data)-1, ]
annot<-read.table(file=paste("/data2/ArrayReAnnotation/FastaFile/Affy/",strsplit(sample,"-")[[1]][2],".annot.uniq.strand.txt",sep=""),header=F,stringsAsFactors=F)
exp1<-lapply(unique(annot[,2]), function (x){colMeans(data[data[,1]%in%annot[which(annot[,2]==x),1],][,-1])})
exp<-do.call(rbind,exp1)
rownames(exp)<-unique(annot[,2])
exp<-exp[exp[,1]!="NaN",]
save(exp,file=paste(strsplit(sample,"-")[[1]][2],"_",strsplit(sample, "-")[[1]][1],".exp.R",sep=""))
exp<-data.frame(GeneSymbol=rownames(exp),exp,stringsAsFactors=F)
write.table(exp,paste(strsplit(sample,"-")[[1]][2],"_",strsplit(sample, "-")[[1]][1],"_express.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")