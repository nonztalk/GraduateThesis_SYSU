library(limma)

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
sampleGPL <- args[2]

tryCatch({
	rawData <<- read.ilmn(files = sample, other.columns = "Detection")
	}, error = function(err) {
		print("Wrong probeId")
		rawData <<- read.ilmn(files = sample, probeid = "ID", other.columns = "Detection")
})

normData <- backgroundCorrect(rawData, method = "normexp", normexp.method = "rma")
normData <- normalizeBetweenArrays(normData, method = "quantile")

expressed <- rowSums(normData$other$Detection < 0.05) >= 2
normData <- normData[expressed, ]
data <- data.frame(probeId = rownames(normData), normData$E, row.names = 1:length(rownames(normData)))

annot <- read.table(file = paste("/data2/ArrayReAnnotation/FastaFile/Illumina", sampleGPL, ".annot.uniq.strand.txt",sep = ""), header = F, stringsAsFactors = F)
exp1 <- lapply(unique(annot[,2]), function(x) {colMeans(data[data[,1] %in% annot[which(annot[,2] == x), 1], ][, -1])})
exp <- do.call(rbind,exp1)
rownames(exp) <- unique(annot[,2])
exp <- exp[exp[,1] != "NaN", ]

filename <- paste(sampleGPL, "_", strsplit(sample, "_")[[1]][1], sep = "")
save(exp,file = paste(filename, ".exp.R", sep = ""))
exp <- data.frame(GeneSymbol = rownames(exp), exp, stringsAsFactors = F)
write.table(exp, paste(filename, "_express.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t") 
