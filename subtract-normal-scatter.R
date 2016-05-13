samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv")
counts <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv", row.names = 1)

# ---- D06r (53% infiltration)
tf <- samples$Infiltration[samples$Alias=="D06r"] / 100 # tumor fraction

# to improve this analysis, we should remove genes NOT differentially expressed b/w tumor and BM
# and correlate only the remaining genes!

d <- data.frame(T06d = counts[,"T06d"]+1, D06r = counts[,"D06r"]+1)
d.adj <- data.frame(T06d = counts[,"T06d"]+1, D06r.adj = ((counts[,"D06r"] - (1-tf) * counts[,"M06r"]) / tf) + 1)
max.value <- max(d, d.adj)

par(mfrow=c(1,2))

plot(d[,1], d[,2], log="xy", cex=0.2, col=rgb(0,0,0,0.2), xlab=names(d)[1], ylab=names(d)[2], xlim=c(1, max.value), ylim=c(1, max.value), main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

plot(d.adj[,1], d.adj[,2], log="xy", cex=0.2, col=rgb(0,0,0,0.2), xlab=names(d.adj)[1], ylab=names(d.adj)[2], xlim=c(1, max.value), ylim=c(1, max.value), main=sprintf("n=%d, Spearman=%.2f", nrow(d.adj), cor(d.adj, method="spearman")[2]))
abline(0, 1, col="black")

# fancy scatter with contour map

tf <- samples$Infiltration[samples$Alias=="D06r"] / 100
d <- data.frame(patient="p06 53%", source="raw", gene=rownames(counts), tumor = counts[,"T06d"]+1, dtc = counts[,"D06r"]+1)
d <- rbind(d, data.frame(patient="p06 53%", source="adjusted", gene=rownames(counts), tumor = counts[,"T06d"]+1, dtc = ((counts[,"D06r"] - (1-tf) * counts[,"M06r"]) / tf) + 1))
tf <- samples$Infiltration[samples$Alias=="D02d"] / 100
d <- rbind(d, data.frame(patient="p02 59%", source="raw", gene=rownames(counts), tumor = counts[,"T02d"]+1, dtc = counts[,"D02d"]+1))
d <- rbind(d, data.frame(patient="p02 59%", source="adjusted", gene=rownames(counts), tumor = counts[,"T02d"]+1, dtc = ((counts[,"D02d"] - (1-tf) * counts[,"M02d"]) / tf) + 1))
d$source = factor(d$source, levels=c("raw", "adjusted"))
d$dtc[d$dtc<1] <- 1
d <- d[d$tumor >= 2 | d$dtc >= 2,]

library(ggplot2)
png("/mnt/projects/fikret/results/expression-tumor-vs-dtc-raw-adjusted.png", width=2400, height=2000, res=300)
p <- ggplot(data=d, aes(tumor, dtc)) +
  facet_wrap(~patient+source) +
  geom_point(alpha=0.3, size=0.7) + 
  scale_x_log10(limits=c(1,1000)) + 
  scale_y_log10(limits=c(1,1000)) +
  stat_density2d(aes(fill=..level.., alpha=..level..), geom = "polygon", colour = "black", size=0.1) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "red", "yellow"))(100)) +
#  geom_smooth(method=lm, linetype=2, colour="red", se=F) +
  geom_abline(intercept=0, slope=1) +
  guides(alpha="none") +
  labs(color="Density", fill="Density", x="Normalized expression tumor (tpm)", y="Normalized expression DTC (tpm)") + 
  theme_bw() +
  coord_fixed() +
  ggtitle("Tumor vs. DTC")
print(p)
dev.off()
#  theme(legend.position=c(0,1), legend.justification=c(0,1)) +
