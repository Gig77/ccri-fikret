options(warn=1)
library("DESeq2")

files <- list.files(path="~/fikret/results/htseq/", pattern=".count$")

# remove samples from initial proof-of-principle project
files <- files[!files %in% c("01-15-BM100.count", "02-15-BM30.count", "03-15-BM0.count", "04-6-BM100.count", "05-6-BM30.count", "06-6-BM0.count")]
names <- unlist(strsplit(files, ".count"))
samples <- data.frame(name=names, file=files, stringsAsFactors=F)

# remove excluded samples
samples.ex <- samples[grep("-x", samples$name, invert=T),] # remove excluded samples

# transform counts into normalized values
cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples.ex, directory="~/fikret/results/htseq", design=~1)

# regularized log transformation
rld <- rlog(cds)
rlogMat <- assay(rld)

# variance stabilizing transformation
#vsd <- varianceStabilizingTransformation(cds)
#vstMat <- assay(vsd)

# write table
