samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv")
mito <- read.delim("/mnt/projects/fikret/data/mito_genes.tsv")
rrna <- read.delim("/mnt/projects/fikret/data/rrna_genes.tsv")

fpm <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv", row.names=1)
dim(fpm)

# what is going on with mitochondrial genes??
#sums <- sort(apply(fpm[rownames(fpm) %in% mito$Ensembl,], 2, sum))
#pdf("/mnt/projects/fikret/results/expr-mito-genes-sum-fpms.pdf")
#plot(sums, xaxt = "n", main="Mitochondrial gene expression", ylab="sum(fpm)")
#axis(1, at=1:length(sums), labels=names(sums), las=2, cex.axis=0.4)
#dev.off()

names.dtc <- grep("^D", colnames(fpm), value=T)
#names.dtc <- grep("NA$", invert=T, names.dtc, value=T)
names.dtc
names.mnc <- grep("^M", colnames(fpm), value=T)
names.mnc <- names.mnc[!names.mnc %in% c("M17r", "M05d", "M16r", "M08d", "M26r", "M25r", "M14r")] # MNCs contaminated with DTCs
names.mnc <- names.mnc[!names.mnc %in% c("M06r", "M33r", "M22r", "M03d", "M08r")]  # outliers in PCA of gene-independent expression metrics
names.mnc <- names.mnc[!names.mnc %in% c("M23r")]  # outliers in PCA of gene expression values
names.mnc
names.tum <- grep("^T", colnames(fpm), value=T)
names.tum <- names.tum[!names.tum %in% c("T69d", "T01d", "T03d1", "T03d2", "T04d", "T04d2", "TWNB6NCT", "TWNB6F", "T06rm", "T46d", "T58d", "T67d", "T71c", "T03d3", "T73NA", "T70c", "T3027", "T3752", "T0588", "T2376", "T4860", "T3231", "T2887", "T0364", "T2383", "T1578", "T0454", "T2421", "T1804", "T4436", "T1095", "T2617", "T0317", "T0206")]
names.tum

# only genes expressed in DTC, MNC, or TUM
min.fpm <- 1
max.fpm <- Inf
expr.filt <- fpm[(apply(fpm[,names.dtc], 1, mean) >= min.fpm & apply(fpm[,names.dtc], 1, mean) <= max.fpm) | 
                 (apply(fpm[,names.mnc], 1, mean) >= min.fpm & apply(fpm[,names.mnc], 1, mean) <= max.fpm) |
                 (apply(fpm[,names.tum], 1, mean) >= min.fpm & apply(fpm[,names.tum], 1, mean) <= max.fpm),]
dim(expr.filt)

# exclude mitochondrial and rRNA genes
#expr.filt <- expr.filt[!rownames(expr.filt) %in% mito$Ensembl,] ; dim(expr.filt)
#expr.filt <- expr.filt[!rownames(expr.filt) %in% rrna$Ensembl,] ; dim(expr.filt)

# remove genes with extreme outlier counts in at least one sample
expr.filt <- expr.filt[!apply(expr.filt, 1, max) >= 20000,] ; dim(expr.filt)

# eliminate values < 1
expr.filt <- expr.filt + 1 

# some final quality checks
sums <- sort(apply(expr.filt[,c(names.dtc, names.mnc)], 2, sum))
plot(sums, xaxt = "n", ylim=c(0,max(sums)), ylab="sum(fpm)")
axis(1, at=1:length(sums), labels=names(sums), las=2, cex.axis=0.5)

# get going
expr.filt.dtc <- expr.filt[,names.dtc] ; sum(expr.filt.dtc)
expr.filt.mnc <- expr.filt[,names.mnc] ; sum(expr.filt.mnc)

library(ISOpureR)
set.seed(123)
ISOpureS1model <- ISOpure.step1.CPE(expr.filt.dtc, expr.filt.mnc)
save(ISOpureS1model, file = "/mnt/projects/fikret/results/isopure/ISOpureS1model.RData")

# correlate with experimentally measured purities
samples <- samples[grepl("^D", samples$Alias),c("Alias", "Infiltration")]
samples <- samples[!is.na(samples$Infiltration),]

purities <- data.frame(sample=names.dtc, purity.computational=ISOpureS1model$alphapurities*100, purity.experimental=samples$Infiltration[match(names.dtc, samples$Alias)])
purities <- purities[order(purities$purity.computational, decreasing=T),]
write.table(purities, "/mnt/projects/fikret/results/isopure/dtc.purity-estimates.isopure.txt", col.names = T, row.names = F, sep = "\t", quote = F)

pdf("/mnt/projects/fikret/results/isopure/dtc.purity-estimates.woContaminatedMNC.fpm.expressed-in-mnc-or-dtc.isopure.pdf")
plot(purities$purity.computational, purities$purity.experimental, xlim=c(0,100), ylim=c(0,100), main=sprintf("Computational vs. experimental\ninfiltration estimates of DTCs (%d genes)", nrow(expr.filt.dtc)), xlab="Computational infiltration rate", ylab="Experimental infiltration rate", cex=2, pch=20)
abline(0, 1, lty=2)
plot(purities$purity.computational, purities$purity.experimental, xlim=c(0,100), ylim=c(0,100), main=sprintf("Computational vs. experimental\ninfiltration estimates of DTCs (%d genes)", nrow(expr.filt.dtc)), xlab="Computational infiltration rate", ylab="Experimental infiltration rate", cex=0)
abline(0, 1, lty=2)
text(purities$purity.computational, purities$purity.experimental, purities$sample, cex=0.5)
fit <- lm(purity.experimental~purity.computational, purities)
#abline(fit, col="red")
dev.off()

# 2nd step: patient profile estimation (PPE)
set.seed(456)
ISOpureS2model <- ISOpure.step2.PPE(expr.filt.dtc, expr.filt.mnc, ISOpureS1model)
ISOpureS2model
rownames(ISOpureS2model$cc_cancerprofiles) <- rownames(expr.filt.dtc)
colnames(ISOpureS2model$cc_cancerprofiles) <- colnames(expr.filt.dtc)
names(ISOpureS2model$alphapurities) <- colnames(expr.filt.dtc)
save(ISOpureS2model, file = "/mnt/projects/fikret/results/isopure/ISOpureS2model.RData")

# scatter plot
expr.filt.dtc.isopure <- ISOpureS2model$cc_cancerprofiles
rownames(expr.filt.dtc.isopure) <- rownames(expr.filt.dtc)
colnames(expr.filt.dtc.isopure) <- colnames(expr.filt.dtc)

s <- "D03d"
d <- data.frame(raw=expr.filt.dtc[,"D35d"], isopure=expr.filt.dtc.isopure[,s])
purity <- ISOpureS2model$alphapurities[colnames(expr.filt.dtc) == s]
r.squared <- summary(lm(log(d[,"isopure"])~log(d[,"raw"])))$r.squared
print(ggplot(data=d, aes(x=raw, y=isopure)) +
  geom_point(alpha=0.3, size=0.7) + 
  scale_x_log10(limits=c(1,1000), breaks=c(1,10,100,1000)) + 
  scale_y_log10(limits=c(1,1000), breaks=c(1,10,100,1000)) +
  stat_density2d(aes(fill=..level.., alpha=..level..), geom = "polygon", colour = "black", size=0.1) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "red", "yellow"))(100)) +
  geom_abline(intercept=0, slope=1, size=0.3) +
  guides(alpha="none") +
  labs(color="Density", fill="Density", x=sprintf("%s raw", s), y=sprintf("%s ISOpure", s)) +
  theme_bw() +
  coord_fixed() +
  geom_text(aes(x=1000, y=1, label=sprintf("Inf=%.2f R=%.2f", purity, r.squared), size=4, hjust=1, vjust=0), show_guide = FALSE))

# try VST transformed expression matrix

library(DESeq2)
counts <- read.delim("/mnt/projects/fikret/results/anduril/execute/htseqExprMatrix/countArray/all.csv", row.names = 1)
samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv", stringsAsFactors = F)
rownames(samples) <- samples$Alias
samples <- samples[match(names(counts), samples$Alias),]
deseq <- DESeqDataSetFromMatrix(countData = counts, colData=samples, design = ~1)	
vst <- assay(varianceStabilizingTransformation(deseq)) ; dim(vst)

names.dtc <- grep("^D", colnames(counts), value=T)
names.dtc
names.mnc <- grep("^M", colnames(counts), value=T)
names.mnc <- names.mnc[!names.mnc %in% c("M17r", "M05d", "M16r", "M08d", "M26r", "M25r", "M14r")]  # MNCs contaminated with DTCs
names.mnc
names.tum <- grep("^T", colnames(counts), value=T)
names.tum

#vst.filt <- vst[rowMeans(vst[,names.tum]) > quantile(rowMeans(vst[,names.tum]), 0.99) |
#                rowMeans(vst[,names.mnc]) > quantile(rowMeans(vst[,names.mnc]), 0.99),]
vst.filt <- vst[apply(vst[,names.dtc], 1, mean) >= 2 |
                apply(vst[,names.mnc], 1, mean) >= 2,]
dim(vst.filt)
vst.linear <- 2^vst.filt
vst.linear[vst.linear<1] <- 1

library(ISOpureR)
set.seed(123)
ISOpureS1model <- ISOpure.step1.CPE(vst.linear[,names.dtc], vst.linear[,names.mnc])
ISOpureS1model

purities <- data.frame(sample=names.dtc, purity.computational=ISOpureS1model$alphapurities*100, purity.experimental=samples$Infiltration[match(names.dtc, samples$Alias)])
pdf("/mnt/projects/fikret/results/isopure/dtc.purity-estimates.woContaminatedMNC.vst.expressed-in-mnc-or-dtc.isopure.pdf")
plot(purities$purity.computational, purities$purity.experimental, xlim=c(0,100), ylim=c(0,100), main=sprintf("Computational vs. experimental\ninfiltration estimates of DTCs (%d genes)", nrow(vst.filt)), xlab="Computational infiltration rate", ylab="Experimental infiltration rate", cex=0)
abline(0, 1, lty=2)
text(purities$purity.computational, purities$purity.experimental, purities$sample, cex=0.5)
dev.off()
