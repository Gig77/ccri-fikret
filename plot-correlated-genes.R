options(warn=1)
library("DESeq2")

files <- list.files(path="~/fikret/results/htseq/", pattern=".count$")
	
# remove samples from initial proof-of-principle project and samples that failed QC
files <- files[!files %in% c("01-15-BM100.count", "02-15-BM30.count", "03-15-BM0.count", "04-6-BM100.count", "05-6-BM30.count", "06-6-BM0.count")]
names <- unlist(strsplit(files, ".count"))
samples <- data.frame(name=names, file=files, stringsAsFactors=F)
samples.ex <- samples[grep("-x", samples$name, invert=T),] # remove excluded samples
	
# get counts of remaining samples
cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples.ex, directory="~/fikret/results/htseq", design=~1)

# regularized log transformation (takes 2-3 hours to compute)
countfile <- "~/fikret/results/deseq/counts.deseq2.rlogMat.RData"
if(file.exists(countfile)) {
	load(countfile)
} else {
	rld <- rlog(cds)
	rlogMat <- assay(rld)
	
	save(rlogMat, file=countfile)
}

# annotate genes HUGO gene names
biomartfile <- "~/generic/data/ensembl/genes.GRCh37v75.biomart.RData"
if(file.exists(biomartfile)) {
	load(biomartfile)
} else {
	library("biomaRt")
	mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # GRCh37, v75
	genes <- getGene(res.df$id, "ensembl_gene_id", mart)
	save(genes, file=biomartfile)
}

geneids <- data.frame(ensemblid = rownames(rlogMat), stringsAsFactors=F)
geneids <- merge(geneids, genes[,1:2], by.x="ensemblid", by.y="ensembl_gene_id", all.x=T) # add gene annotation
geneids <- geneids[!duplicated(geneids$ensemblid),] # remove duplicate entries due to ensembl ids mapping to multiple symbols
rownames(rlogMat) <- ifelse(is.na(geneids$hgnc_symbol) | geneids$hgnc_symbol == "", geneids$ensemblid, geneids$hgnc_symbol)

# only DTCs
rlogMat.dtc <- rlogMat[,grep("D-", colnames(rlogMat))]

# remove unexpressed genes to speed up correlation analysis
maxs <- apply(counts(cds), 1, max)
rlogMat.dtc.filtered <- rlogMat.dtc[maxs > 10,]

# compute all-vs-all correlation matrix
cormat <- cor(t(rlogMat.dtc.filtered))

# remove genes with low variance
#vars <- apply(rlogMat.dtc, 1, var)
#rlogMat.dtc.filtered <- rlogMat.dtc[vars > quantile(vars, 0.9),]


#------------------------
# MAGEC2 SCATTER PLOTS
#------------------------

# extract top correlated genes (positively and negatively)
magec2 <- cormat["MAGEC2",]
pos <- magec2[magec2 >= 0.6]; pos <- pos[order(pos, decreasing=T)]
neg <- magec2[magec2 <= -0.6]; neg <- neg[order(neg, decreasing=F)]

pdf("~/fikret/results/deseq/magec2.scatter.pos.pdf", height=10, width=14)
layout(cbind(matrix(1:16, nrow=4, byrow=T), c(17,17,17,17)))
par(mar=c(2,4,1,1))
samples <- as.factor(colnames(rlogMat.dtc))
for(i in c(2:length(pos))) {
	x <- rlogMat.dtc["MAGEC2",]
	y <- rlogMat.dtc[names(pos)[i],]
	fit <- lm(y~x)
	test <- cor.test(x,y)
	plot(x, y, xlab="", ylab=names(pos)[i], col=samples, pch=as.integer(samples) %% 25, ylim=c(min(0, y), max(y)), cex.lab=1.3)
	abline(fit, col="red")
	text(4.5, min(0, y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0)
	if ((i-1) %% 16 == 0) {
		plot.new()
		legend("topleft", legend=samples, col=samples, pch=as.integer(samples) %% 25)
	}
}
dev.off()

pdf("~/fikret/results/deseq/magec2.scatter.neg.pdf", height=10, width=14)
layout(cbind(matrix(1:16, nrow=4, byrow=T), c(17,17,17,17)))
par(mar=c(2,4,1,1))
samples <- as.factor(colnames(rlogMat.dtc))
for(i in c(2:length(neg))) {
	x <- rlogMat.dtc["MAGEC2",]
	y <- rlogMat.dtc[names(neg)[i],]
	fit <- lm(y~x)
	test <- cor.test(x,y)
	plot(x, y, xlab="", ylab=names(neg)[i], col=samples, pch=as.integer(samples) %% 25, ylim=c(min(0, y), max(y)), cex.lab=1.3)
	abline(fit, col="red")
	text(4.5, min(0, y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0)
	if ((i-1) %% 16 == 0) {
		plot.new()
		legend("topleft", legend=samples, col=samples, pch=as.integer(samples) %% 25)
	}
}
dev.off()

#------------------------
# MYCN SCATTER PLOTS
#------------------------

# extract top correlated genes (positively and negatively)
mycn <- cormat["MYCN",]
pos <- mycn[mycn >= 0.6]; pos <- pos[order(pos, decreasing=T)]
neg <- mycn[mycn <= -0.6]; neg <- neg[order(neg, decreasing=F)]

write.table(pos, file="~/fikret/results/deseq/mycn-corr-pos.tsv", col.names=T, row.names=T, sep="\t", quote=F)
write.table(neg, file="~/fikret/results/deseq/mycn-corr-neg.tsv", col.names=T, row.names=T, sep="\t", quote=F)

pdf("~/fikret/results/deseq/mycn.scatter.pos.pdf", height=10, width=14)
layout(cbind(matrix(1:16, nrow=4, byrow=T), c(17,17,17,17)))
par(mar=c(2,4,1,1))
samples <- as.factor(colnames(rlogMat.dtc))
for(i in c(2:length(pos))) {
	x <- rlogMat.dtc["MYCN",]
	y <- rlogMat.dtc[names(pos)[i],]
	fit <- lm(y~x)
	test <- cor.test(x,y)
	plot(x, y, xlab="", ylab=names(pos)[i], col=samples, pch=as.integer(samples) %% 25, ylim=c(min(0, y), max(y)), cex.lab=1.3)
	abline(fit, col="red")
	text(10, min(0, y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0)
	if ((i-1) %% 16 == 0) {
		plot.new()
		legend("topleft", legend=samples, col=samples, pch=as.integer(samples) %% 25)
	}
}
dev.off()

pdf("~/fikret/results/deseq/mycn.scatter.neg.pdf", height=10, width=14)
layout(cbind(matrix(1:16, nrow=4, byrow=T), c(17,17,17,17)))
par(mar=c(2,4,1,1))
samples <- as.factor(colnames(rlogMat.dtc))
for(i in c(2:length(neg))) {
	x <- rlogMat.dtc["MYCN",]
	y <- rlogMat.dtc[names(neg)[i],]
	fit <- lm(y~x)
	test <- cor.test(x,y)
	plot(x, y, xlab="", ylab=names(neg)[i], col=samples, pch=as.integer(samples) %% 25, ylim=c(min(0, y), max(y)), cex.lab=1.3)
	abline(fit, col="red")
	text(10, min(0, y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0)
	if ((i-1) %% 16 == 0) {
		plot.new()
		legend("topleft", legend=samples, col=samples, pch=as.integer(samples) %% 25)
	}
}
dev.off()

#------------------------
# MAGEC2 AND MYCN LINE PLOT
#------------------------

pdf("~/fikret/results/deseq/magec2.vs.mycn.pdf", height=7, width=14)
samples <- c(grep("-Dx", colnames(rlogMat.dtc)), grep("-Dx", colnames(rlogMat.dtc), invert=T)) 
par(mar=c(8,4,1,2))
matplot(t(rlogMat.dtc[c("MAGEC2", "MYCN"), samples]), type="l", xaxt='n', ylab="MAGEC2                                 MYCN")
par(new=T)
matplot(t(rlogMat.dtc[c("MAGEC2", "MYCN"), samples]), pch=19, xaxt='n', yaxt='n', xlab="", ylab="")
axis(1, at=c(1:length(colnames(rlogMat.dtc))), las=2, cex.axis=0.8, labels=colnames(rlogMat.dtc)[samples])
dev.off()

rlogMat.dtc.corr <- rlogMat.dtc[c(names(pos), names(neg)),]
matplot(rlogMat.dtc.corr, type="l")

magec2 <- rlogMat.dtc["ENSG00000046774",]

magec2.cor (cor(rlogMat.dtc))
# annotate genes HUGO gene names
biomartfile <- "~/generic/data/ensembl/genes.GRCh37v75.biomart.RData"
if(file.exists(biomartfile)) {
	load(biomartfile)
} else {
	library("biomaRt")
	mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl") # GRCh37, v75
	genes <- getGene(res.df$id, "ensembl_gene_id", mart)
	save(genes, file=biomartfile)
}

# compute row-wise (gene-gene) correlation matrix
# WARNING: takes > 40G memory!

# pull out correlations with gene of interest
# MAGEC2: ENSG00000046774
# only DTCs!
# correlate only gene of interest!


res.annotated <- merge(res.df, genes[,1:3], by.x="id", by.y="ensembl_gene_id", all.x=T) # add gene annotation

