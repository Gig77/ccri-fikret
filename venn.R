options(warn=1)
library(VennDiagram)

t.vs.m <- read.delim("/mnt/projects/fikret/results/patients-919-644.T-vs-M.deseq2.tsv", check.names=F, stringsAsFactors=F)
t.vs.d <- read.delim("/mnt/projects/fikret/results/patients-919-644.T-vs-D.deseq2.tsv", check.names=F, stringsAsFactors=F)
d.vs.m <- read.delim("/mnt/projects/fikret/results/patients-919-644.D-vs-M.deseq2.tsv", check.names=F, stringsAsFactors=F)
set1 <- t.vs.m[!is.na(t.vs.m$padj) & t.vs.m$padj<0.05,"id"]
set2 <- t.vs.d[!is.na(t.vs.d$padj) & t.vs.d$padj<0.05, "id"]
set3 <- d.vs.m[!is.na(d.vs.m$padj) & d.vs.m$padj<0.05, "id"]

pdf("/mnt/projects/fikret/results/venn.pdf")
grid.draw(draw.triple.venn(length(set1), length(set2), length(set3), length(intersect(set1, set2)), length(intersect(set2, set3)), length(intersect(set1, set3)), length(intersect(intersect(set1, set2), set3)), category=c(sprintf("T-vs-M\n(%d)", length(set1)), sprintf("T-vs-D\n(%d)", length(set2)), sprintf("D-vs-M (%d)", length(set3)))))
dev.off()
