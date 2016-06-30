library(edgeR)
counts <- read.delim("/mnt/projects/fikret/results/anduril/execute/htseqExprMatrixDGEA/countArray/all.csv", row.names = 1)
groups <- read.delim("/mnt/projects/fikret/data/sample_groups.csv")
samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv")
samples <- samples[samples$Alias %in% colnames(counts),]

samples$group <- rep(NA, nrow(samples))
samples$group[samples$Alias %in% unlist(strsplit(groups$Members[groups$ID=="DTC"], ","))] <- "DTC"
samples$group[samples$Alias %in% unlist(strsplit(groups$Members[groups$ID=="TUM"], ","))] <- "TUM"
samples$group[is.na(samples$group)] <- "Other"
samples$group <- factor(samples$group, levels=c("DTC", "TUM", "Other"))

y <- DGEList(counts=counts)
y <- calcNormFactors(y)

design <- model.matrix(~0+TumorFraction+group, data=samples)
rownames(design) <- colnames(y)

#y <- estimateGLMCommonDisp(y, design)
#y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
DTCvTUM <- glmLRT(fit, contrast=makeContrasts(groupDTC+TumorFraction-groupTUM, levels=design))
topTags(DTCvTUM)
