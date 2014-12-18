options(warn=1)

TUvsDTC <- read.delim("~/fikret/results/deseq/TU-vs-DTC.diagnosis.full.tsv")
TUvsBM <- read.delim("~/fikret/results/deseq/TU-vs-BM.diagnosis.full.tsv")
DTCvsBM <- read.delim("~/fikret/results/deseq/DTC-vs-BM.diagnosis.tsv")

merged <- merge(TUvsDTC, TUvsBM, by="id", suffixes=c(".TUvsDTC", ".TUvsBM"), all.x=T)
merged <- merge(merged, DTCvsBM, by="id", suffixes=c(".TUvsDTC", ".DTCvsBM"), all.x=T)

names(merged)[68:80] <- paste0(names(merged)[68:80], ".DTCvsBM")
merged <- merged[,!names(merged) %in% c("hgnc_symbol.TUvsBM", "description.TUvsBM", "chromosome_name.TUvsBM", "start_position.TUvsBM", "end_position.TUvsBM",
						                "hgnc_symbol.DTCvsBM", "description.DTCvsBM", "chromosome_name.DTCvsBM", "start_position.DTCvsBM", "end_position.DTCvsBM")]

# boxplots to estimate tumor cell content
deg <- rowMeans(merged[,15:22]) < 20 & rowMeans(merged[,51:62]) > 200
data <- merged[deg,23:34]
medians <- apply(data, 2, median)
data <- data[order(medians)]
colnames(data) <- gsub(".TUvsDTC", "", colnames(data))
pdf("~/fikret/results/tumor-cell-content.pdf")
par(mar=c(8,4,1,1))
boxplot(data, ylim=c(0,300), las=2, cex.axis=0.7, ylab=sprintf("No. of reads for %d BM marker genes", nrow(data)))
dev.off()

write.table(merged[(!is.na(merged$padj.TUvsDTC) & merged$padj.TUvsDTC < 0.1) | (!is.na(merged$padj.TUvsBM) & merged$padj.TUvsBM < 0.1) | (!is.na(merged$padj.DTCvsBM) & merged$padj.DTCvsBM < 0.1),], file="~/fikret/results/deseq/TU-vs-DTC.TU-vs-BM.DTC-vs-BM.diagnosis.tsv", row.names=F, col.names=T, quote=F, sep="\t")
