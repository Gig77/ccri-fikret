m <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqAnnotated_DTCdx_vs_DTCrel/table.csv")
m <- m[!duplicated(m$ids),]
rownames(m) <- m$ids
s <- read.delim("/mnt/projects/fikret/data/sample_key.csv", check.names = F, stringsAsFactors = F)

gene.names <- c("RPL18", "RPL14", "RPS29", "RPS19", "RPL8", "FAU", "RPL37A", "RPS11")
genes <- m$ids[m$Gene %in% c(gene.names)]
dia <- grep("^D.*d\\d?$", colnames(m), value=T)
rel <- grep("^D.*r\\d?$", colnames(m), value=T)

library(reshape)
d <- melt(m[genes,], id.vars="Gene", measure.vars=c(dia, rel))
d$group <- ifelse(d$variable %in% dia, "dia", "rel")
d$Gene <- factor(d$Gene, levels=c(gene.names))

ggplot(d, aes(x=group, y=value)) + 
  facet_wrap(~Gene) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.1)) + 
  theme_bw()

