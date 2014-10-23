options(warn=1)
library("DESeq2")
library("RColorBrewer")
library("gplots")

files <- list.files(path="~/fikret/results/htseq/", pattern=".count$")

# remove samples from initial proof-of-principle project
files <- files[!files %in% c("01-15-BM100.count", "02-15-BM30.count", "03-15-BM0.count", "04-6-BM100.count", "05-6-BM30.count", "06-6-BM0.count")]
names <- unlist(strsplit(files, ".count"))
site <- sapply(strsplit(names, "-"), "[[", 1)
site[site=="D2"] <- "D"
site[site=="T2"] <- "T"
site <- factor(site, levels=c("M", "D", "T"))
samples <- data.frame(name=names, file=files, site=site, stringsAsFactors=F)

#-----------------------
# WITH EXCLUDED SAMPLES
#-----------------------
# transform counts into normalized values
cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples, directory="~/fikret/results/htseq", design=~site)
#dds <- DESeq(cds)
rld <- rlog(cds)

# heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pdf("~/fikret/results/sample-dist.heatmap.with-outliers.pdf.part", width=12, height=12)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

# PCA
pdf("~/fikret/results/sample-dist.pca.with-outliers.pdf.part")
p <- plotPCA(rld, intgroup=c("site"), col=c("red", "blue", "black"))
p$panel.args.common$cex <- 0.7
p <- update(p, panel = function(x, y, ...) {
			lattice::panel.xyplot(x, y, ...);
			lattice::ltext(x=x, y=y, labels=rownames(colData(rld)), pos=1, offset=1, cex=0.3)
		})
print(p)
dev.off()

#-----------------------
# WITHOUT EXCLUDED SAMPLES
#-----------------------
samples.ex <- samples[grep("-x", samples$name, invert=T),] # remove excluded samples

# build sample data frame

# transform counts into normalized values
cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples.ex, directory="~/fikret/results/htseq", design=~site)
#dds <- DESeq(cds)
rld <- rlog(cds)

# heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
pdf("~/fikret/results/sample-dist.heatmap.without-outliers.pdf.part", width=12, height=12)
heatmap.2(mat, trace="none", col=rev(hmcol), margin=c(13, 13))
dev.off()

# PCA
pdf("~/fikret/results/sample-dist.pca.without-outliers.pdf.part")
p <- plotPCA(rld, intgroup=c("site"), col=c("red", "blue", "black"))
p$panel.args.common$cex <- 0.7
p <- update(p, panel = function(x, y, ...) {
			lattice::panel.xyplot(x, y, ...);
			lattice::ltext(x=x, y=y, labels=rownames(colData(rld)), pos=1, offset=1, cex=0.3)
		})
print(p)
dev.off()
