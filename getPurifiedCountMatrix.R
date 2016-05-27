# original count matrix
countMatrix  <- as.matrix(read.delim("/mnt/projects/fikret/results/anduril/execute/htseqExprMatrix/countArray/all.csv", row.names = 1, stringsAsFactors = F, check.names = F))

# purified DTC "counts"
load("/mnt/projects/fikret/results/isopure/ISOpureS2model.RData")
tpm.isopure <- ISOpureS2model$cc_cancerprofiles
counts.isopure <- round(tpm.isopure)

# merge
countMatrix <- countMatrix[,!colnames(countMatrix) %in% colnames(counts.isopure)]
countMatrix <- merge(countMatrix, counts.isopure, by="row.names")
colnames(countMatrix)[1] <- "ID"

# write
write.table(countMatrix, "/mnt/projects/fikret/results/counts.ISOpure.csv", row.names = F, col.names = T, quote = F, sep = "\t")
