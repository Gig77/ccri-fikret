options(warn=1)
library("DESeq2")

samples <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples[nrow(samples)+1,] <- c("01-15-BM100", "01-15-BM100.count", "BM100")
samples[nrow(samples)+1,] <- c("02-15-BM30", "02-15-BM30.count", "BM30")
samples[nrow(samples)+1,] <- c("03-15-BM0", "03-15-BM0.count", "BM0")
samples[nrow(samples)+1,] <- c("04-6-BM100", "04-6-BM100.count", "BM100")
samples[nrow(samples)+1,] <- c("05-6-BM30", "05-6-BM30.count", "BM30")
samples[nrow(samples)+1,] <- c("06-6-BM0", "06-6-BM0.count", "BM0")
samples[nrow(samples)+1,] <- c("T-93-0644-06-Dx", "T-93-0644-06-Dx.count", "T")
samples[nrow(samples)+1,] <- c("M-93-0644-06-Dx", "M-93-0644-06-Dx.count", "M")
samples[nrow(samples)+1,] <- c("D-93-0644-06-Dx", "D-93-0644-06-Dx.count", "D")
samples[nrow(samples)+1,] <- c("T-94-0919-05-Dx", "T-94-0919-05-Dx.count", "T")
samples[nrow(samples)+1,] <- c("M-94-0919-05-Dx", "M-94-0919-05-Dx.count", "M")
samples[nrow(samples)+1,] <- c("D-94-0919-05-Dx", "D-94-0919-05-Dx.count", "D")

#=====================================================================================
# scatter plots
#=====================================================================================

#--------
# DESeq2 - counts
#--------
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples, directory="~/fikret/results/htseq", design=~condition)
counts.raw <- as.data.frame(counts(cds))

cds <- estimateSizeFactors(cds)
sizeFactors(cds)
counts <- as.data.frame(counts(cds, normalized=T))

#-------
# cell line spike ins
#-------

png(file="~/fikret/results/scatter_STA-NB-15_STA-NB-6.png", width=4000, height=3000, res=300)
par(mfrow=c(2,3))

#plot(ifelse(counts[,"01-15-BM100"] >= 1, log2(counts[,"01-15-BM100"]), 0), ifelse(counts[,"02-15-BM30"] >= 1, log2(counts[,"02-15-BM30"]), 0), cex=0.1, col=rgb(0,0,0,0.2), xlab="log2 rc 01-15-BM100", ylab="log2 rc 02-15-BM30")
#smoothScatter(log2(counts[,"01-15-BM100"]), log2(counts[,"02-15-BM30"]), nbin=400, bandwidth=0.2)

plot(log2(counts[,"01-15-BM100"]), log2(counts[,"02-15-BM30"]), cex=0.1, col=rgb(0,0,0,0.2), xlab="log2 rc 01-15-BM100", ylab="log2 rc 02-15-BM30")
plot(log2(counts[,"01-15-BM100"]), log2(counts[,"03-15-BM0"]), cex=0.1, col=rgb(0,0,0,0.2), xlab="log2 rc 01-15-BM100", ylab="log2 rc 03-15-BM0")
plot(log2(counts[,"02-15-BM30"]), log2(counts[,"03-15-BM0"]), cex=0.1, col=rgb(0,0,0,0.2), xlab="log2 rc 02-15-BM30", ylab="log2 rc 03-15-BM0")

plot(log2(counts[,"04-6-BM100"]), log2(counts[,"05-6-BM30"]), cex=0.1, col=rgb(0,0,0,0.2), xlab="log2 rc 04-6-BM100", ylab="log2 rc 05-6-BM30")
plot(log2(counts[,"04-6-BM100"]), log2(counts[,"06-6-BM0"]), cex=0.1, col=rgb(0,0,0,0.2), xlab="log2 rc 04-6-BM100", ylab="log2 rc 06-6-BM0")
plot(log2(counts[,"05-6-BM30"]), log2(counts[,"06-6-BM0"]), cex=0.1, col=rgb(0,0,0,0.2), xlab="log2 rc 05-6-BM30", ylab="log2 rc 06-6-BM0")

dev.off()

png(file="~/fikret/results/scatter_644_919.png", width=4000, height=3000, res=300)

#-------
# primary samples
#-------

# ---- 644

par(mfrow=c(2,3))

d <- data.frame(x = log2(counts[,"T-93-0644-06-Dx"]), y = log2(counts[,"D-93-0644-06-Dx"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) T-93-0644-06-Dx", ylab="log2(nrc) D-93-0644-06-Dx", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

d <- data.frame(x = log2(counts[,"T-93-0644-06-Dx"]), y = log2(counts[,"M-93-0644-06-Dx"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) T-93-0644-06-Dx", ylab="log2(nrc) M-93-0644-06-Dx", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

tf <- 0.43
d <- data.frame(x = log2(counts[,"T-93-0644-06-Dx"]), y = log2((counts[,"D-93-0644-06-Dx"] - (1-tf) * counts[,"M-93-0644-06-Dx"]) / tf))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y) & d$y > 0,]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) T-93-0644-06-Dx", ylab="log2(nrc) D-93-0644-06-Dx adjusted", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

# ---- 919

d <- data.frame(x = log2(counts[,"T-94-0919-05-Dx"]), y = log2(counts[,"D-94-0919-05-Dx"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) T-94-0919-05-Dx", ylab="log2(nrc) D-94-0919-05-Dx", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

d <- data.frame(x = log2(counts[,"T-94-0919-05-Dx"]), y = log2(counts[,"M-94-0919-05-Dx"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) T-94-0919-05-Dx", ylab="log2(nrc) M-94-0919-05-Dx", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

tf <- 0.9
d <- data.frame(x = log2(counts[,"T-94-0919-05-Dx"]), y = log2((counts[,"D-94-0919-05-Dx"] - (1-tf) * counts[,"M-94-0919-05-Dx"]) / tf))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y) & d$y > 0,]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) T-94-0919-05-Dx", ylab="log2(nrc) D-94-0919-05-Dx adjusted", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

dev.off()

#=====================================================================================
# differentially expressed genes
#=====================================================================================

library("biomaRt")
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # GRCh37, v75
set.seed(343)

#------
# DEGs primary samples
#------

# -- tumor vs. bone marrow
samples.TvsM <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.TvsM[nrow(samples.TvsM)+1,] <- c("T-93-0644-06-Dx", "T-93-0644-06-Dx.count", "T")
samples.TvsM[nrow(samples.TvsM)+1,] <- c("T-94-0919-05-Dx", "T-94-0919-05-Dx.count", "T")
samples.TvsM[nrow(samples.TvsM)+1,] <- c("M-93-0644-06-Dx", "M-93-0644-06-Dx.count", "M")
samples.TvsM[nrow(samples.TvsM)+1,] <- c("M-94-0919-05-Dx", "M-94-0919-05-Dx.count", "M")
samples.DvsM$condition <- relevel(as.factor(samples.DvsM$condition), "M")

cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.TvsM, directory="~/fikret/results/htseq", design=~condition)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
counts.norm <- as.data.frame(counts(cds, normalized=T))
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
res.annotated <- merge(res.annotated, counts.norm, by.x="id", by.y="row.names", all.x=T)
write.table(res.annotated, file="~/fikret/results/patients-919-644.T-vs-M.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# -- tumor vs. disseminated tumor cells in bone marrow
samples.TvsD <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.TvsD[nrow(samples.TvsD)+1,] <- c("T-93-0644-06-Dx", "T-93-0644-06-Dx.count", "T")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("T-94-0919-05-Dx", "T-94-0919-05-Dx.count", "T")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("D-94-0919-05-Dx", "D-94-0919-05-Dx.count", "D")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("D-93-0644-06-Dx", "D-93-0644-06-Dx.count", "D")
samples.DvsM$condition <- relevel(as.factor(samples.DvsM$condition), "D")

rm(cds, dds, counts.norm, res, res.df, genes, res.annotated)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.TvsD, directory="~/fikret/results/htseq", design=~condition)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
counts.norm <- as.data.frame(counts(cds, normalized=T))
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
res.annotated <- merge(res.annotated, counts.norm, by.x="id", by.y="row.names", all.x=T)
write.table(res.annotated, file="~/fikret/results/patients-919-644.T-vs-D.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# -- disseminated tumor cells in bone marrow vs. bone marrow
samples.DvsM <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.DvsM[nrow(samples.DvsM)+1,] <- c("D-94-0919-05-Dx", "D-94-0919-05-Dx.count", "D")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("D-93-0644-06-Dx", "D-93-0644-06-Dx.count", "D")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("M-93-0644-06-Dx", "M-93-0644-06-Dx.count", "M")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("M-94-0919-05-Dx", "M-94-0919-05-Dx.count", "M")
samples.DvsM$condition <- relevel(as.factor(samples.DvsM$condition), "M")

rm(cds, dds, counts.norm, res, res.df, genes, res.annotated)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.DvsM, directory="~/fikret/results/htseq", design=~condition)
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
counts.norm <- as.data.frame(counts(cds, normalized=T))
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
res.annotated <- merge(res.annotated, counts.norm, by.x="id", by.y="row.names", all.x=T)
write.table(res.annotated, file="~/fikret/results/patients-919-644.D-vs-M.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#------
# DEGs cell line spike ins
#------

# -- BM100 vs. BM30
samples.BM100vsBM30 <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("01-15-BM100", "01-15-BM100.count", "BM100")
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("04-6-BM100", "04-6-BM100.count", "BM100")
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("02-15-BM30", "02-15-BM30.count", "BM30")
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("05-6-BM30", "05-6-BM30.count", "BM30")

rm(cds, dds, res, res.df, genes, res.annotated)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.BM100vsBM30, directory="~/fikret/results/htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/fikret/results/BM100-vs-BM30.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# -- BM100 vs. BM0
samples.BM100vsBM0 <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("01-15-BM100", "01-15-BM100.count", "BM100")
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("04-6-BM100", "04-6-BM100.count", "BM100")
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("03-15-BM0", "03-15-BM0.count", "BM0")
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("06-6-BM0", "06-6-BM0.count", "BM0")

rm(cds, dds, res, res.df, genes, res.annotated)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.BM100vsBM0, directory="~/fikret/results/htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/fikret/results/BM100-vs-BM0.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# -- BM30 vs. BM0
samples.BM30vsBM0 <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("02-15-BM30", "02-15-BM30.count", "BM30")
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("05-6-BM30", "05-6-BM30.count", "BM30")
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("03-15-BM0", "03-15-BM0.count", "BM0")
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("06-6-BM0", "06-6-BM0.count", "BM0")

rm(cds, dds, res, res.df, genes, res.annotated)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.BM30vsBM0, directory="~/fikret/results/htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/fikret/results/BM30-vs-BM0.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

