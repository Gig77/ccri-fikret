options(warn=1)
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("edgeR")


samples <- read.delim("~/fikret/results/qlucore/sample-annotations.txt")
samples$NameAnn <- paste(samples$Name, samples$TCC) # annotated name
samples$NameAnn <- paste(samples$NameAnn, samples$MYCN) # annotated name
samples$NameAnn <- paste(samples$NameAnn, samples$Coverage) # annotated name
samples <- samples[,c("NameAnn", "Filename", colnames(samples)[!colnames(samples) %in% c("NameAnn", "Filename")])] # reorder columns compatible with DeSeqDataSetFromHTSeqCount()
samples$Filename <- gsub(".gsnap.filtered.bam", ".count", samples$Filename)

cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples, directory="~/fikret/results/htseq", design=~1)

# keep only expressed genes
#cds <- cds[rowSums(cpm(DGEList(counts=counts(cds))) > 10 ) >= 2,] 

#---------------------
# DESeq2 rlog
#---------------------
rld <- rlog(cds)

# heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pdf("~/fikret/results/sample-dist.heatmap.deseq2.rlog.pdf", width=12, height=12)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

pdf("~/fikret/results/sample-dist.heatmap.deseq2.rlog.customclustering.pdf", width=12, height=12)
heatmap.2(as.matrix(distsRL), hclustfun=function(x) hclust(x, method="average"), distfun=function(x) as.dist(x), trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

# PCA
pdf("~/fikret/results/sample-dist.pca.deseq2.rlog.pdf")
p <- plotPCA(rld, intgroup=c("site"), col=c("red", "blue", "black"))
p$panel.args.common$cex <- 0.7
p <- update(p, panel = function(x, y, ...) {
			lattice::panel.xyplot(x, y, ...);
			lattice::ltext(x=x, y=y, labels=rownames(colData(rld)), pos=1, offset=1, cex=0.3)
		})
print(p)
dev.off()

#---------------------
# DeSeq2 VTS
#---------------------

#maxs <- apply(counts(cds), 1, max) ; cds <- cds[maxs > 20,]
rld <- varianceStabilizingTransformation(cds)
#rld <- rlog(cds.trimmed, blind=T, fast=F)

# heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pdf("~/fikret/results/sample-dist.heatmap.deseq2.VST.pdf", width=12, height=12)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

#---------------------
# Voom
#---------------------

dge <- DGEList(counts=counts(cds))
dge.norm <- calcNormFactors(dge, method="TMM")
y <- voom(dge.norm)

# heatmap
distsRL <- dist(t(y$E))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pdf("~/fikret/results/sample-dist.heatmap.voom.pdf", width=12, height=12)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()
