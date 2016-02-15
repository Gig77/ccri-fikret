tpm <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseqExprMatrix/expr.csv", row.names=1)
tpm <- 2^tpm # transform back to non-log scale

names.dtc <- grep("^D", colnames(tpm), value=T) ; names.dtc
names.mnc <- grep("^M", colnames(tpm), value=T) ; names.mnc
names.tum <- grep("^T", colnames(tpm), value=T) ; names.tum

# select top-N most differentially expressed genes between tumor and normal 
TvsM <- read.delim("/mnt/projects/fikret/results/anduril/execute/deseq_DTC_highiniltration_vs_DTC_lowiniltration/results.csv")
TvsM <- TvsM[!is.na(TvsM$qDTC_highiniltration_vs_DTC_lowiniltration) & !is.na(TvsM$fcDTC_highiniltration_vs_DTC_lowiniltration),]
TvsM <- TvsM[order(TvsM$qDTC_highiniltration_vs_DTC_lowiniltration),]
topUp <- as.character(TvsM$ids[TvsM$fcDTC_highiniltration_vs_DTC_lowiniltration >= 2][1:200])
topUp <- topUp[!is.na(topUp)]
topDn <- as.character(TvsM$ids[TvsM$fcDTC_highiniltration_vs_DTC_lowiniltration <= -2][1:200])
topDn <- topDn[!is.na(topDn)]

# test
tpm.filt <- tpm[topDn,]
df <- data.frame(group="dtc", value=rowSums(tpm.filt[,names.dtc]))
df <- rbind(df, data.frame(group="mnc", value=rowSums(tpm.filt[,names.mnc])))
boxplot(value~group, df, log="y")

# re-order tpm to have genes most different between tumor and normal at the top
#tpm.filt <- tpm[c(topUp, topDn),]
tpm.filt <- tpm[rowSums(tpm[,names.tum])>=10,]
dim(tpm.filt)
tpm.filt[tpm.filt<1] <- 1
apply(tpm.filt, 2, sum)

tpm.filt.dtc <- tpm.filt[,names.dtc]
tpm.filt.mnc <- tpm.filt[,names.mnc]
sum(tpm.filt.dtc)
sum(tpm.filt.mnc)

library(ISOpureR)
set.seed(123)
ISOpureS1model <- ISOpure.step1.CPE(tpm.filt.dtc, tpm.filt.mnc)
ISOpureS1model

# correlate with experimentally measured purities
samples <- read.delim("/mnt/projects/fikret/data/sample_key.csv")
samples <- samples[grepl("^D", samples$Alias),c("Alias", "Infiltration")]
samples <- samples[!is.na(samples$Infiltration),]

purities <- data.frame(sample=colnames(tpm.filt.dtc), purity.computational=ISOpureS1model$alphapurities*100, purity.experimental=samples$Infiltration[match(colnames(tpm.filt.dtc), samples$Alias)])
purities <- purities[order(purities$purity.computational, decreasing=T),]

pdf("/mnt/projects/fikret/results/dtc.purity-estimates.isopure.pdf")
plot(purities$purity.computational, purities$purity.experimental, xlim=c(0,100), ylim=c(0,100), main=sprintf("Computational vs. experimental\npurity estimates of DTCs (%d genes)", nrow(tpm.filt.dtc)), xlab="Computational (% purity)", ylab="Experimental (% purity)", cex=0)
text(purities$purity.computational, purities$purity.experimental, purities$sample, cex=0.5)
fit <- lm(purity.experimental~purity.computational, purities)
abline(fit, col="red")
dev.off()

#write.table(purities, "/mnt/projects/fikret/results/dtc.purity-estimates.isopure.txt", col.names = T, row.names = F, sep = "\t", quote = F)

