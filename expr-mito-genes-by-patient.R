library(reshape)
library(ggplot2)

e <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv") ; names(e)[1] <- "ensembl"
mito <- read.delim("/mnt/projects/fikret/data/mito_genes.tsv")
e <- e[e$ensembl %in% mito$Ensembl,]
e.long <- melt(e, "ensembl") ; names(e.long) <- c("ensembl", "sample", "expression")
e.long$tp <- factor(gsub("([MDT])(\\d+)([dr])", "\\1\\3", e.long$sample), levels=c("Td", "Md", "Mr", "Dr", "Dd"))
e.long$patient <- gsub("([MDT])(\\d+)([dr])", "\\2", e.long$sample)

pdf("/mnt/projects/fikret/results/expr-mito-genes-by-patient.pdf", width=14)
p <- ggplot(e.long[e.long$patient %in% c("02", "05", "06", "07"),], aes(x=tp, y=expression)) + 
  theme_bw() +
  facet_wrap(~ patient, nrow=1) +
  geom_boxplot() +
  geom_point(alpha=0.3, color="blue") +
  geom_line(aes(group=ensembl), alpha=0.3, color="blue")
print(p)
dev.off()