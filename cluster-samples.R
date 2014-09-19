options(warn=1)
library("DESeq2")
library("RColorBrewer")
library("gplots")

files <- list.files(path="~/fikret/results/htseq/", pattern=".count$")

# remove samples from initial proof-of-principle project
files <- files[!files %in% c("01-15-BM100.count", "02-15-BM30.count", "03-15-BM0.count", "04-6-BM100.count", "05-6-BM30.count", 
				             "06-6-BM0.count", "07-644-T.count", "08-644-M.count", "09-644-D.count", "10-919-T.count", "11-919-M.count", "12-919-D.count")]

# build sample data frame
names <- unlist(strsplit(files, ".count"))
site <- sapply(strsplit(names, "-"), "[[", 1)
site[site=="D2"] <- "D"
site <- factor(site, levels=c("M", "D", "T"))
samples <- data.frame(name=names, file=files, site=site, stringsAsFactors=F)

# transform counts into normalized values
cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples, directory="~/fikret/results/htseq", design=~site)
dds <- DESeq(cds)
rld <- rlog(dds)

# heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pdf("~/fikret/results/sample-dist.heatmap.pdf", width=12, height=12)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

# PCA
pdf("~/fikret/results/sample-dist.pca.pdf")
p <- plotPCA(rld, intgroup=c("site"))
p <- update(p, panel = function(x, y, ...) {
			lattice::panel.xyplot(x, y, ...);
			lattice::ltext(x=x, y=y, labels=rownames(colData(rld)), pos=1, offset=1, cex=0.4)
		})
print(p)
dev.off()

rld <- rlow(dds)
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
samples.TvsD[nrow(samples.TvsD)+1,] <- c("07-644-T", "07-644-T.count", "T")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("10-919-T", "10-919-T.count", "T")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("12-919-D", "12-919-D.count", "D")
samples.TvsD[nrow(samples.TvsD)+1,] <- c("09-644-D", "09-644-D.count", "D")
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
samples.DvsM[nrow(samples.DvsM)+1,] <- c("12-919-D", "12-919-D.count", "D")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("09-644-D", "09-644-D.count", "D")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("08-644-M", "08-644-M.count", "M")
samples.DvsM[nrow(samples.DvsM)+1,] <- c("11-919-M", "11-919-M.count", "M")
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

