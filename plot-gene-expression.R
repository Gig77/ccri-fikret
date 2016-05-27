library(reshape)
library(ggplot2)

# read data

tpm <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv")
names(tpm)[1] <- "ensembl"
tpm <- melt(tpm, "ensembl")
names(tpm) <- c("ensembl", "sample", "expression")
tpm$tp <- factor(gsub("([MDT])(\\d+).*", "\\1", tpm$sample), levels=c("T", "M", "D"))
tpm$patient <- gsub("([MDT])(\\d+)([dr])", "\\2", tpm$sample)

ggplot(tpm[tpm$sample %in% c("M33d", "D33d", "M11r", "D11r") & 
           tpm$ensembl %in% c("ENSG00000158578", "ENSG00000163554", "ENSG00000112077", "ENSG00000133742",
                              "ENSG00000157005", "ENSG00000104725", "ENSG00000124882", "ENSG00000169429",
                              "ENSG00000109132", "ENSG00000162706"),], 
       aes(x=patient, y=expression)) + 
  theme_bw() +
  facet_wrap(~ tp+ensembl, nrow=3) +
  geom_point() +
  geom_line(aes(group=ensembl))
