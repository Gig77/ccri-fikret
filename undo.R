#source("https://bioconductor.org/biocLite.R")
#biocLite("UNDO")
library(UNDO)

# vignette example
data(NumericalMixMCF7HS27)
X <- NumericalMixMCF7HS27
data(NumericalMixingMatrix)
A <- NumericalMixingMatrix
data(PureMCF7HS27)
S <- exprs(PureMCF7HS27)
two_source_deconv(X,lowper=0.4,highper=0.1,epsilon1=0.01,epsilon2=0.01,A,S[,1],S[,2],return=0)

# real data

counts <- read.delim("/mnt/projects/fikret/results/anduril/execute/htseqExprMatrix/countArray/all.csv", row.names = 1)
counts <- counts[,c("D01d", "D02d", "D04d", "D06d", "D07d", "D07r", "D08r", "D10r", "D33r")]

library(edgeR)
m <- DGEList(counts=counts)
m <- calcNormFactors(m, method="TMM")
y <- voom(m)
logCPM <- y$E
logCPM <- logCPM + abs(min(logCPM)) # avoid negative values
logCPM <- logCPM[rowSums(logCPM)>0,]
X <- new("ExpressionSet", exprs=logCPM)

two_source_deconv(X, lowper=0.4, highper=0.1, epsilon1=0.01, epsilon2=0.01, return = 0)

# garbage results.... darn