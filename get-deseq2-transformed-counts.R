options(warn=1)
library("DESeq2")

files <- list.files(path="~/fikret/results/htseq/", pattern=".count$")

# remove samples from initial proof-of-principle project
files <- files[!files %in% c("01-15-BM100.count", "02-15-BM30.count", "03-15-BM0.count", "04-6-BM100.count", "05-6-BM30.count", "06-6-BM0.count")]
names <- unlist(strsplit(files, ".count"))
samples <- data.frame(name=names, file=files, stringsAsFactors=F)

# remove excluded samples
samples.ex <- samples[grep("-x", samples$name, invert=T),] # remove excluded samples

# add annotations
ann <- read.delim(file="~/fikret/scripts/cluster-samples.R")
# transform counts into normalized values
cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples.ex, directory="~/fikret/results/htseq", design=~1)

# regularized log transformation
rld <- rlog(cds)
rlogMat <- assay(rld)

# variance stabilizing transformation
#vsd <- varianceStabilizingTransformation(cds)
#vstMat <- assay(vsd)

# annotate gene names
#rlogMat <- read.delim("~/fikret/results/deseq/normalized-counts.deseq2.rlog.tsv", check.names=F)
load("~/generic/data/ensembl/genes.GRCh37v75.biomart.RData")
gMat.hgnc <- merge(rlogMat, genes[,c("ensembl_gene_id", "hgnc_symbol")], by.x="gene", by.y="ensembl_gene_id", all.x=T)
n <- names(rlogMat.hgnc)

# write table
write.table(rlogMat.hgnc[n[c(1,length(n),2:(length(n)-1))]], file="~/fikret/results/deseq/normalized-counts.deseq2.rlog.hgnc.tsv", col.names=T, row.names=F, sep="\t", quote=F)

