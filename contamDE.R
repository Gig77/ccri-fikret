#install.packages("nloptr")
#install.packages("/data_synology/software/contamDE-1.0.tar.gz", repos = NULL, type="source")
library(contamDE)

# demo data set
#install.packages("ARTIVA")
library(ARTIVA)
data(prostate)
d = contamDE(prostate[,-1], R=2, match=TRUE)

# fikret's data
counts <- read.delim("/mnt/projects/fikret/results/anduril/execute/htseqExprMatrix/countArray/all.csv", row.names = 1)
#counts.matched <- counts[,c("M01d", "M02d", "M04d", "M06d", "M07d", "M07r", "M08r", "M10r", "M33r", 
#                            "D01d", "D02d", "D04d", "D06d", "D07d", "D07r", "D08r", "D10r", "D33r")]
counts.matched <- counts[,c("M01d", "M02d", "M04d", "M11r", 
                            "D01d", "D02d", "D04d", "D11r")]
#counts.matched <- counts.matched[rowSums(counts.matched) >= 10,]
counts.matched.filt <- counts.matched[rowSums(counts.matched) >= 10,]

d <- contamDE(counts.matched.filt, R=2, match=TRUE)
d$W

# try VST
library(DESeq2)
samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv", stringsAsFactors = F)
rownames(samples) <- samples$Alias
samples <- samples[match(names(counts.matched), samples$Alias),]
deseq <- DESeqDataSetFromMatrix(countData = counts.matched, colData=samples, design = ~1)	
vst <- assay(varianceStabilizingTransformation(deseq)) ; dim(vst)
#vst <- vst[apply(vst, 1, sd) != 0 & rowMeans(vst) > quantile(rowMeans(vst), 0.8),] ; dim(vst)
vst <- vst[apply(vst, 1, sd) > quantile(apply(vst, 1, sd), 0.9),] ; dim(vst)
d <- contamDE(2^vst, R=2, match=TRUE)

# try voomanized expression matrix
library(edgeR)
m <- DGEList(counts=counts.matched.filt)
m <- calcNormFactors(m, method="TMM")
y <- voom(m)
d <- contamDE(2^y$E, R=2, match=TRUE)
d$W

