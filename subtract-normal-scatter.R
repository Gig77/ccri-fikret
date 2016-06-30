# fancy scatter with contour map

samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv")
fpm <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv", row.names = 1)
fpm <- fpm[rowSums(fpm) >= 100,]

# filter for genes diff. expressed b/w TUM and MNC
TvM <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseq_TUMvsMNC/results.csv", stringsAsFactors = F)
TvM <- TvM[!is.na(TvM$qTUMvsMNC) & TvM$qTUMvsMNC <= 1e-5 & abs(TvM$fcTUMvsMNC) >= 1,]
fpm <- fpm[rownames(fpm) %in% TvM$ids,]

# get ISOpure results (= purified DTC profiles)
load("/mnt/projects/fikret/results/isopure/ISOpureS2model.RData")
fpm.isopure <- ISOpureS2model$cc_cancerprofiles

# only genes present in both matrices
shared <- intersect(rownames(fpm), rownames(fpm.isopure))
fpm <- fpm[shared,]
fpm.isopure <- as.data.frame(fpm.isopure[shared,])

tfs.isopure <- ISOpureS2model$alphapurities
name.ref <- "D35d"

# DTCs without matching MNCs (excluded or missing): 
name.dtc.nomnc <- c("05d", "08d", "14r", "16r", "17r", "18r", "25r", "26r", "35d", 
                    "36NA", "37d", "39d", "40d", "42d", "43d", "61d", "69r", "70r", "74NA")
name.dtc.withmnc <- c("01d", "02d", "03d", "04d", "06d", "06r", "07d", "07r", "08r", "10d", "10r", 
                      "11r", "12r", "15r", "20r", "21r", "22r", "23r", "24r", "27d", "30d", "31d",
                      "33d", "33r", "28d")
library(ggplot2)
for (s in c(name.dtc.withmnc, name.dtc.nomnc)) {
  print(s)

  d <- data.frame(src=character(0), gene=character(0), dtc_ref=numeric(0), dtc_sample=numeric(0))
  d <- rbind(d, data.frame(src="Raw", gene=rownames(fpm), dtc_ref = fpm[,name.ref]+1, dtc_sample = fpm[,paste0("D", s)]+1))
  
  tf.exp <- samples$Infiltration[samples$Alias==paste0("D", s)] / 100
#  if (s %in% name.dtc.withmnc) {
#    src.exp <- paste0("Purified (exp. ", tf.exp*100, "%)")
#    d <- rbind(d, data.frame(src=src.exp, gene=rownames(fpm), dtc_ref = fpm[,name.ref]+1, dtc_sample = ((fpm[,paste0("D", s)] - (1-tf.exp) * fpm[,paste0("M", s)]) / tf.exp) + 1))
#  }

  tf.isopure <- tfs.isopure[paste0("D", s)]
  if (!is.na(tf.isopure)) {
    src.isopure <- sprintf("Purified (assuming %.0f%% IR)", tf.isopure*100)
    d <- rbind(d, data.frame(src=src.isopure, gene=rownames(fpm.isopure), dtc_ref = fpm[,name.ref]+1, dtc_sample = fpm.isopure[,paste0("D", s)]))
  }
  d <- d[d$dtc_ref > 1 & d$dtc_sample > 1,]
  d$src <- factor(d$src, levels=unique(d$src))
  
  # compute fit
  fit <- data.frame(src=character(0), p=numeric(0), R=numeric(0))
  fit.mix <- lm(log(d$dtc_sample[grepl("Raw", d$src)])~log(d$dtc_ref[grepl("Raw", d$src)]))
  fit <- rbind(fit, data.frame(src="Raw", p=anova(fit.mix)$'Pr(>F)'[1], R=summary(fit.mix)$r.squared))
#  if (s %in% name.dtc.withmnc) {
#    fit.exp <- lm(log(d$dtc_sample[grepl("exp.", d$src)])~log(d$dtc_ref[grepl("exp.", d$src)]))
#    fit <- rbind(fit, data.frame(src=src.exp, p=anova(fit.exp)$'Pr(>F)'[1], R=summary(fit.exp)$r.squared))
#  }
  if (!is.na(tf.isopure)) {
    fit.isopure <- lm(log(d$dtc_sample[grepl("Purified", d$src)])~log(d$dtc_ref[grepl("Purified", d$src)]))
    fit <- rbind(fit, data.frame(src=src.isopure, p=anova(fit.isopure)$'Pr(>F)'[1], R=summary(fit.isopure)$r.squared))
  }
  fit$src <- factor(fit$src, levels=unique(fit$src))
  
  png(paste0("/mnt/projects/fikret/results/scatterplots/expression-", "D", s, "-vs-", name.ref, "-raw-and-adjusted.png"), width=2400, height=2000, res=300)
  print(ggplot(data=d, aes(dtc_ref, dtc_sample)) +
    facet_wrap(~src) +
    geom_point(alpha=0.3, size=0.7) + 
    scale_x_log10(limits=c(1,1000), breaks=c(1,10,100,1000)) + 
    scale_y_log10(limits=c(1,1000), breaks=c(1,10,100,1000)) +
    stat_density2d(aes(fill=..level.., alpha=..level..), geom = "polygon", colour = "black", size=0.1) +
    scale_fill_gradientn(colours = colorRampPalette(c("blue", "red", "yellow"))(100)) +
    #  geom_smooth(method=lm, linetype=2, colour="red", se=F) +
    geom_abline(intercept=0, slope=1, size=0.3) +
    guides(alpha="none") +
    labs(color="Density", fill="Density", 
         x=paste0("Normalized expression ", name.ref, " (fpm)"), 
         y=paste0("Normalized expression ", "D", s, " (fpm)")) + 
    theme_bw() +
    coord_fixed() +
    ggtitle(paste0("D", s, " (", tf.exp*100, "% IR) vs. ", name.ref, " (96% IR)")) +
    geom_text(aes(x=1000, y=1, label=sprintf("R=%.2f", R), size=4, hjust=1, vjust=0), data=fit, show_guide = FALSE))
  dev.off()
}

# ===============================================================================================

# D07r vs D07r2

fpm <- read.delim("/mnt/projects/fikret/results/exprMatrix.deseq2.fpm.csv", row.names = 1)
TvM <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseq_TUM_vs_MNC/results.csv", stringsAsFactors = F)
TvM <- TvM[!is.na(TvM$qTUM_vs_MNC) & TvM$qTUM_vs_MNC <= 1e-5 & abs(TvM$fcTUM_vs_MNC) >= 1,]
#names.mnc <- grep("^M", colnames(fpm), value=T)
#names.mnc <- names.mnc[!names.mnc %in% c("M17r", "M05d", "M16r", "M08d", "M26r", "M25r", "M14r")] # MNCs contaminated with DTCs
#names.mnc <- names.mnc[!names.mnc %in% c("M06r", "M33r", "M22r", "M03d", "M08r")]  # outliers in PCA of gene-independent expression metrics
#names.mnc <- names.mnc[!names.mnc %in% c("M23r")]  # outliers in PCA of gene expression values
#names.tum <- grep("^T", colnames(fpm), value=T)
#names.tum <- names.tum[!names.tum %in% c("T69d", "T01d", "T03d1", "T03d2", "T04d", "T04d2", "TWNB6NCT", "TWNB6F", "T06rm", "T46d", "T58d", "T67d", "T71c", "T03d3", "T73NA", "T70c", "T3027", "T3752", "T0588", "T2376", "T4860", "T3231", "T2887", "T0364", "T2383", "T1578", "T0454", "T2421", "T1804", "T4436", "T1095", "T2617", "T0317", "T0206")]
fpm <- fpm[rownames(fpm) %in% TvM$ids,]
fpm <- fpm[rowSums(fpm) >= 100,]
fpm <- fpm + 1

# get ISOpure results (= purified DTC profiles)
load("/mnt/projects/fikret/results/isopure/ISOpureS2model.RData")
fpm.isopure <- ISOpureS2model$cc_cancerprofiles
tfs.isopure <- ISOpureS2model$alphapurities

# only genes present in both matrices
shared <- intersect(rownames(fpm), rownames(fpm.isopure))
fpm <- fpm[shared,]
fpm.isopure <- as.data.frame(fpm.isopure[shared,])
fpm.isopure <- fpm.isopure + 1

# build dataframe(s)
d <- data.frame(src=character(0), gene=character(0), dtc_ref=numeric(0), dtc_sample=numeric(0))
fit <- data.frame(src=character(0), p=numeric(0), R=numeric(0))

x <- fpm[,"D07r"] ; y <- fpm[,"D07r2"]
src <- "Raw BE vs. raw AE"
d <- rbind(d, data.frame(src=src, gene=rownames(fpm), dtc_ref = x, dtc_sample = y))
fit <- rbind(fit, data.frame(src=src, p=anova(lm(log(y)~log(x)))$'Pr(>F)'[1], R=summary(lm(log(y)~log(x)))$r.squared))

x <- fpm[,"D07r"] ; y <- fpm.isopure[,"D07r2"]
src <- sprintf("Purified BE vs. raw AE", tfs.isopure[paste0("D07r2")]*100)
d <- rbind(d, data.frame(src=src, gene=rownames(fpm), dtc_ref = x, dtc_sample = y))
fit <- rbind(fit, data.frame(src=src, p=anova(lm(log(y)~log(x)))$'Pr(>F)'[1], R=summary(lm(log(y)~log(x)))$r.squared))

# purified vs purified: VERY high correlation (R=0.99)
#x <- fpm.isopure[,"D07r"] ; y <- fpm.isopure[,"D07r2"]
#src <- sprintf("Purified (%.0f%% IR) vs. purified (%.0f%% IR)", tfs.isopure[paste0("D07r")]*100, tfs.isopure[paste0("D07r2")]*100)
#d <- rbind(d, data.frame(src=src, gene=rownames(fpm), dtc_ref = x, dtc_sample = y))
#fit <- rbind(fit, data.frame(src=src, p=anova(lm(log(y)~log(x)))$'Pr(>F)'[1], R=summary(lm(log(y)~log(x)))$r.squared))

# purified vs purified (control, different patients): also VERY high correlation (R=0.96); hm.....
#x <- fpm.isopure[,"D70r"] ; y <- fpm.isopure[,"D07r2"]
#src <- sprintf("D70r (%.0f%% IR) vs. D07r2 (%.0f%% IR)", tfs.isopure[paste0("D70r")]*100, tfs.isopure[paste0("D07r2")]*100)
#d <- rbind(d, data.frame(src=src, gene=rownames(fpm), dtc_ref = x, dtc_sample = y))
#fit <- rbind(fit, data.frame(src=src, p=anova(lm(log(y)~log(x)))$'Pr(>F)'[1], R=summary(lm(log(y)~log(x)))$r.squared))

d <- d[d$dtc_ref > 1 & d$dtc_sample > 1,]
d$src <- factor(d$src, levels=unique(d$src))
fit$src <- factor(fit$src, levels=unique(fit$src))

png(paste0("/mnt/projects/fikret/results/scatterplots/expression-", "D07r-vs-D07r2-before-after-enrichment.png"), width=2400, height=2000, res=300)
#pdf(paste0("/mnt/projects/fikret/results/scatterplots/expression-", "D07r-vs-D07r2-before-after-enrichment.pdf"))
print(ggplot(data=d, aes(x=dtc_ref, y=dtc_sample)) +
  facet_wrap(~src, nrow = 1) +
  geom_point(alpha=0.3, size=0.6) + 
  scale_x_log10(limits=c(1,1000), breaks=c(1,10,100,1000)) + 
  scale_y_log10(limits=c(1,1000), breaks=c(1,10,100,1000)) +
  stat_density2d(aes(fill=..level.., alpha=..level..), geom = "polygon", colour = "black", size=0.1, bins=9) +
  scale_alpha_continuous(range=c(0.1,1)) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "red", "yellow"))(100)) +
  geom_abline(intercept=0, slope=1, size=0.3) +
  guides(alpha="none") +
  labs(color="Density", fill="Density", x="D07r after enrichment (91% IR)", y="D07r before enrichment (36% IR)") +
  theme_bw() +
  coord_fixed() +
  geom_text(aes(x=1000, y=1, label=sprintf("R=%.2f", R), size=4, hjust=1, vjust=0), data=fit, show_guide = FALSE))
dev.off()


#===================================================================================================
# old code

samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv")
fpm <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv", row.names = 1)

# ---- D06r (53% infiltration)
tf <- samples$Infiltration[samples$Alias=="D06r"] / 100 # tumor fraction

# to improve this analysis, we should remove genes NOT differentially expressed b/w tumor and BM
# and correlate only the remaining genes!

d <- data.frame(T06d = fpm[,"T06d"]+1, D06r = fpm[,"D06r"]+1)
d.adj <- data.frame(T06d = fpm[,"T06d"]+1, D06r.adj = ((fpm[,"D06r"] - (1-tf) * fpm[,"M06r"]) / tf) + 1)
max.value <- max(d, d.adj)

par(mfrow=c(1,2))

plot(d[,1], d[,2], log="xy", cex=0.2, col=rgb(0,0,0,0.2), xlab=names(d)[1], ylab=names(d)[2], xlim=c(1, max.value), ylim=c(1, max.value), main=sprintf("n=%d, Spearman=%.2f", nrow(d), cor(d, method="spearman")[2]))
abline(0, 1, col="black")

plot(d.adj[,1], d.adj[,2], log="xy", cex=0.2, col=rgb(0,0,0,0.2), xlab=names(d.adj)[1], ylab=names(d.adj)[2], xlim=c(1, max.value), ylim=c(1, max.value), main=sprintf("n=%d, Spearman=%.2f", nrow(d.adj), cor(d.adj, method="spearman")[2]))
abline(0, 1, col="black")
