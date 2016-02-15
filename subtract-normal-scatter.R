counts <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv")

# ---- D06r (53% infiltration)

# to improve this analysis, we should remove genes NOT differentially expressed b/w tumor and BM
# and correlate only the remaining genes!

par(mfrow=c(1,2))

d <- data.frame(T06d = counts[,"T06d"], D06r = counts[,"D06r"])
plot(d[,1], d[,2], cex=0.1, col=rgb(0,0,0,0.2), xlab=names(d)[1], ylab=names(d)[2], main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

tf <- 0.53  # tumor fraction
d <- data.frame(T06d = counts[,"T06d"], D06r.adj = log((exp(counts[,"D06r"]) - (1-tf) * exp(counts[,"M06r"])) / tf))
d <- d[complete.cases(d) & !is.infinite(d[,1]) & !is.infinite(d[,2]),]
plot(d[,1], d[,2], cex=0.1, col=rgb(0,0,0,0.2), xlab=names(d)[1], ylab=names(d)[2], main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")
