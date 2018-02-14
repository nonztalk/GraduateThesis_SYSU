args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
AnnotationFile <- args[2]

data <- read.table(paste(sample, "_norm.txt", sep = ""), header = T)
head(data)
fileInfo <- read.table("fileInfo.txt", header = F, stringsAsFactors = F)
annot<-read.table(file = paste("/data2/ArrayReAnnotation/FastaFile/Affy/", AnnotationFile,".probe_fasta.annot.uniq.strand.txt",sep = ""),header=F,stringsAsFactors=F)
head(annot)
colnames(data)<-c("probeset_id",fileInfo[,2])
exp1 <- lapply(unique(annot[,2]), function (x){colMeans(data[data[,1] %in% annot[which(annot[,2]==x),1],][,-1])})
exp <- do.call(rbind,exp1)
rownames(exp) <- unique(annot[,2])
head(exp)
exp <- exp[exp[,1]!="NaN",]
exp <- log2(exp)
save(exp,file=paste(sample,".exp.R",sep=""))
write.table(exp,paste(sample,"_express.txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
