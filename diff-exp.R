options(warn=1)
library("DESeq2")

samples <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples[nrow(samples)+1,] <- c("01-15-BM100", "C3JP8ACXX_01-15-BM100_14s001408-1-1_Rifatbegovic_lane114s001408_sequence.count", "BM100")
samples[nrow(samples)+1,] <- c("02-15-BM30", "C3JP8ACXX_02-15-BM30_14s001409-1-1_Rifatbegovic_lane114s001409_sequence.count", "BM30")
samples[nrow(samples)+1,] <- c("03-15-BM0", "C3JP8ACXX_03-15-BM0_14s001410-1-1_Rifatbegovic_lane114s001410_sequence.count", "BM0")
samples[nrow(samples)+1,] <- c("04-6-BM100", "C3JP8ACXX_04-6-BM100_14s001411-1-1_Rifatbegovic_lane114s001411_sequence.count", "BM100")
samples[nrow(samples)+1,] <- c("05-6-BM30", "C3JP8ACXX_05-6-BM30_14s001412-1-1_Rifatbegovic_lane114s001412_sequence.count", "BM30")
samples[nrow(samples)+1,] <- c("06-6-BM0", "C3JP8ACXX_06-6-BM0_14s001413-1-1_Rifatbegovic_lane114s001413_sequence.count", "BM0")
samples[nrow(samples)+1,] <- c("07-644-T", "C3JP8ACXX_07-644-T_14s001414-1-1_Rifatbegovic_lane214s001414_sequence.count", "T")
samples[nrow(samples)+1,] <- c("08-644-M", "C3JP8ACXX_08-644-M_14s001415-1-1_Rifatbegovic_lane214s001415_sequence.count", "M")
samples[nrow(samples)+1,] <- c("09-644-D", "C3JP8ACXX_09-644-D_14s001416-1-1_Rifatbegovic_lane214s001416_sequence.count", "D")
samples[nrow(samples)+1,] <- c("10-919-T", "C3JP8ACXX_10-919-T_14s001417-1-1_Rifatbegovic_lane214s001417_sequence.count", "T")
samples[nrow(samples)+1,] <- c("11-919-M", "C3JP8ACXX_11-919-M_14s001418-1-1_Rifatbegovic_lane214s001418_sequence.count", "M")
samples[nrow(samples)+1,] <- c("12-919-D", "C3JP8ACXX_12-919-D_14s001419-1-1_Rifatbegovic_lane214s001419_sequence.count", "D")

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

d <- data.frame(x = log2(counts[,"07-644-T"]), y = log2(counts[,"09-644-D"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) 07-644-T", ylab="log2(nrc) 09-644-D", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

d <- data.frame(x = log2(counts[,"07-644-T"]), y = log2(counts[,"08-644-M"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) 07-644-T", ylab="log2(nrc) 08-644-M", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

tf <- 0.43
d <- data.frame(x = log2(counts[,"07-644-T"]), y = log2((counts[,"09-644-D"] - (1-tf) * counts[,"08-644-M"]) / tf))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y) & d$y > 0,]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) 07-644-T", ylab="log2(nrc) 09-644-D adjusted", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

# ---- 919

d <- data.frame(x = log2(counts[,"10-919-T"]), y = log2(counts[,"12-919-D"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) 10-919-T", ylab="log2(nrc) 12-919-D", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

d <- data.frame(x = log2(counts[,"10-919-T"]), y = log2(counts[,"11-919-M"]))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y),]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) 10-919-T", ylab="log2(nrc) 11-919-M", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

tf <- 0.9
d <- data.frame(x = log2(counts[,"10-919-T"]), y = log2((counts[,"12-919-D"] - (1-tf) * counts[,"11-919-M"]) / tf))
d <- d[complete.cases(d) & !is.infinite(d$x) & !is.infinite(d$y) & d$y > 0,]
plot(d$x, d$y, cex=0.1, col=rgb(0,0,0,0.2), xlab="log2(nrc) 10-919-T", ylab="log2(nrc) 12-919-D adjusted", main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
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
samples.TvsM[nrow(samples.TvsM)+1,] <- c("07-644-T", "C3JP8ACXX_07-644-T_14s001414-1-1_Rifatbegovic_lane214s001414_sequence.count", "T")
samples.TvsM[nrow(samples.TvsM)+1,] <- c("10-919-T", "C3JP8ACXX_10-919-T_14s001417-1-1_Rifatbegovic_lane214s001417_sequence.count", "T")
samples.TvsM[nrow(samples.TvsM)+1,] <- c("08-644-M", "C3JP8ACXX_08-644-M_14s001415-1-1_Rifatbegovic_lane214s001415_sequence.count", "M")
samples.TvsM[nrow(samples.TvsM)+1,] <- c("11-919-M", "C3JP8ACXX_11-919-M_14s001418-1-1_Rifatbegovic_lane214s001418_sequence.count", "M")

cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.TvsM, directory="~/fikret/results/htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/fikret/results/patients-919-644.T-vs-M.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# -- tumor vs. disseminated tumor cells in bone marrow
samples.TvsD <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.TvsD[nrow(samples.TvsD)+1,] <- c("07-644-T", "C3JP8ACXX_07-644-T_14s001414-1-1_Rifatbegovic_lane214s001414_sequence.count", "T")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("10-919-T", "C3JP8ACXX_10-919-T_14s001417-1-1_Rifatbegovic_lane214s001417_sequence.count", "T")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("12-919-D", "C3JP8ACXX_12-919-D_14s001419-1-1_Rifatbegovic_lane214s001419_sequence.count", "D")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("09-644-D", "C3JP8ACXX_09-644-D_14s001416-1-1_Rifatbegovic_lane214s001416_sequence.count", "D")

rm(cds, dds, res, res.df, genes, res.annotated)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.TvsD, directory="~/fikret/results/htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/fikret/results/patients-919-644.T-vs-D.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

# -- disseminated tumor cells in bone marrow vs. bone marrow
samples.DvsM <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.DvsM[nrow(samples.DvsM)+1,] <- c("12-919-D", "C3JP8ACXX_12-919-D_14s001419-1-1_Rifatbegovic_lane214s001419_sequence.count", "D")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("09-644-D", "C3JP8ACXX_09-644-D_14s001416-1-1_Rifatbegovic_lane214s001416_sequence.count", "D")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("08-644-M", "C3JP8ACXX_08-644-M_14s001415-1-1_Rifatbegovic_lane214s001415_sequence.count", "M")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("11-919-M", "C3JP8ACXX_11-919-M_14s001418-1-1_Rifatbegovic_lane214s001418_sequence.count", "M")

rm(cds, dds, res, res.df, genes, res.annotated)
cds <- DESeqDataSetFromHTSeqCount(sampleTable = samples.DvsM, directory="~/fikret/results/htseq", design=~condition)
dds <- DESeq(cds)
res <- results(dds)
res <- res[order(res$padj),]
res.df <- as.data.frame(res)
res.df$id <- rownames(res.df)

genes <- getGene(res.df$id, "ensembl_gene_id", mart)
res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T)
write.table(res.annotated, file="~/fikret/results/patients-919-644.D-vs-M.deseq2.tsv", col.names=T, row.names=F, sep="\t", quote=F)

#------
# DEGs cell line spike ins
#------

# -- BM100 vs. BM30
samples.BM100vsBM30 <- data.frame(name=character(), file=character(), condition=character(), stringsAsFactors=F)
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("01-15-BM100", "C3JP8ACXX_01-15-BM100_14s001408-1-1_Rifatbegovic_lane114s001408_sequence.count", "BM100")
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("04-6-BM100", "C3JP8ACXX_04-6-BM100_14s001411-1-1_Rifatbegovic_lane114s001411_sequence.count", "BM100")
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("02-15-BM30", "C3JP8ACXX_02-15-BM30_14s001409-1-1_Rifatbegovic_lane114s001409_sequence.count", "BM30")
samples.BM100vsBM30[nrow(samples.BM100vsBM30)+1,] <- c("05-6-BM30", "C3JP8ACXX_05-6-BM30_14s001412-1-1_Rifatbegovic_lane114s001412_sequence.count", "BM30")

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
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("01-15-BM100", "C3JP8ACXX_01-15-BM100_14s001408-1-1_Rifatbegovic_lane114s001408_sequence.count", "BM100")
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("04-6-BM100", "C3JP8ACXX_04-6-BM100_14s001411-1-1_Rifatbegovic_lane114s001411_sequence.count", "BM100")
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("03-15-BM0", "C3JP8ACXX_03-15-BM0_14s001410-1-1_Rifatbegovic_lane114s001410_sequence.count", "BM0")
samples.BM100vsBM0[nrow(samples.BM100vsBM0)+1,] <- c("06-6-BM0", "C3JP8ACXX_06-6-BM0_14s001413-1-1_Rifatbegovic_lane114s001413_sequence.count", "BM0")

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
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("02-15-BM30", "C3JP8ACXX_02-15-BM30_14s001409-1-1_Rifatbegovic_lane114s001409_sequence.count", "BM30")
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("05-6-BM30", "C3JP8ACXX_05-6-BM30_14s001412-1-1_Rifatbegovic_lane114s001412_sequence.count", "BM30")
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("03-15-BM0", "C3JP8ACXX_03-15-BM0_14s001410-1-1_Rifatbegovic_lane114s001410_sequence.count", "BM0")
samples.BM30vsBM0[nrow(samples.BM30vsBM0)+1,] <- c("06-6-BM0", "C3JP8ACXX_06-6-BM0_14s001413-1-1_Rifatbegovic_lane114s001413_sequence.count", "BM0")

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

