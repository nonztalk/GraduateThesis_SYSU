library(affy)

args <- commandArgs(trailingOnly = T)
sample <- args[1]

Data <- ReadAffy()
eset <- rma(Data)
write.exprs(eset, file = paste(sample, "_express.txt", sep = ""))
data <- read.table(paste(sample, "_express.txt", sep = ""), header = T, sep = "\t", row.names = 1)
fileInfo <- read.table("fileInfo.txt", header = F)
colnames(data) <- fileInfo[, 2]

filename <- paste(strsplit(sample, "_")[[1]][1], ".annot.uniq.strand.txt", sep = "")
annot <- read.table(file = paste("/data2/ArrayReAnnotation/FastaFile/Affyprocess/", filename, sep = ""), header = F, stringsAsFactors = F)
exp1 <- lapply(unique(annot[, 2]), function(x){
	colMeans(data[annot[which(annot[, 2] == x), 1], ])
	})
exp <- do.call(rbind,exp1)
rownames(exp) <- unique(annot[,2])
save(exp, file = paste(sample, ".exp.R", sep = ""))
exp <- data.frame(GeneSymbol = rownames(exp), exp, stringsAsFactors = F)
write.table(exp, paste(sample, "_express.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
