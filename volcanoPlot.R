volcanoPlot <- function
(
    expr, sigthresh=0.001, lfcthresh=2, fccol=2, pcol=3, qcol=4, geneNameCol=2, labelTopN=30, 
    legendpos="bottomright", excludeGenes = "", sampleSize=1000, cexLabel=0.4, cexPoint=1, 
    minQ=1e-200, labelNonameGenesWithId=TRUE, onlyGenesWithNames=FALSE, title=""
) 
{
  # get input data and parameters
  expr <- read.delim(expr, stringsAsFactors = F, check.names = F)
  
  # prepare data
  res <- data.frame(id=expr[,1], log2FC=expr[,fccol], pvalue=expr[,pcol], qValue=pmax(expr[,qcol], minQ), stringsAsFactors = F)
  res <- res[!is.na(res$qValue) & !is.na(res$log2FC),]
  
  if (get.input(cf, "geneNames") != "") {
    geneNames <- CSV.read(get.input(cf, "geneNames"))
    res <- merge(res, geneNames[,c(1,geneNameCol)], by.x=1, by.y=1, all.x=T)
    res <- res[!duplicated(res$id),]
    names(res)[5] <- "Name"
    if (onlyGenesWithNames) {
      res <- res[!is.na(res$Name) & res$Name != "",]
    }
    if (labelNonameGenesWithId) {
      res$Name[is.na(res$Name) | res$Name == ""] <- res$id[is.na(res$Name) | res$Name == ""]
    }
  } else {
    res$Name <- res[,1]
  }
  
  # subsample uninteresting fraction of genes to reduce plot size
  res.interesting   <- subset(res, qValue <= sigthresh | abs(log2FC) >= lfcthresh)
  res.uninteresting <- subset(res, qValue > sigthresh & abs(log2FC) < lfcthresh)
  res.uninteresting.sample <- res.uninteresting[sample(1:nrow(res.uninteresting), min(sampleSize, nrow(res.uninteresting)), replace = F),]
  res.sample <- rbind(res.interesting, res.uninteresting.sample)
  
  if (excludeGenes != "") {
    exclude <- unlist(strsplit(excludeGenes,','))
    res.sample <- res.sample[!res.sample$Name %in% exclude,]
  }
  
  with(subset(res.sample, qValue >  sigthresh & abs(log2FC) <  lfcthresh), plot(log2FC, -log10(qValue), pch=20, col="lightgray", xlim=c(min(res.sample$log2FC)-1,max(res.sample$log2FC)+1), ylim=range(-log10(res.sample$qValue)), cex=cexPoint, main=title))
  with(subset(res.sample, qValue <= sigthresh & abs(log2FC) <  lfcthresh), points(log2FC, -log10(qValue), pch=20, col="black", cex=cexPoint))
  with(subset(res.sample, qValue >  sigthresh & abs(log2FC) >= lfcthresh), points(log2FC, -log10(qValue), pch=20, col="blue", cex=cexPoint))
  with(subset(res.sample, qValue <= sigthresh & abs(log2FC) >= lfcthresh), points(log2FC, -log10(qValue), pch=20, col="red", cex=cexPoint))
  
  print(sprintf("Significantly up-regulated: %d genes", sum(res$qValue <= sigthresh & res$log2FC >= lfcthresh)))
  print(sprintf("Significantly down-regulated: %d genes", sum(res$qValue <= sigthresh & res$log2FC <= lfcthresh)))
  
  # label points
  if (labelTopN > 0) {
    library(calibrate)
    
    res.sample.sorted <- res.sample[order(res.sample$qValue),]
    res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$qValue) & res.sample.sorted$log2FC > 0,]
    label.up.sig <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]
    
    res.sample.sorted <- res.sample[order(-res.sample$log2FC),]
    res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$qValue) & res.sample.sorted$log2FC > 0,]
    label.up.fc <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]
    
    with(res.sample.sorted[res.sample.sorted$id %in% c(label.up.sig, label.up.fc),], textxy(log2FC, -log10(qValue), labs=Name, cex=cexLabel, offset=0.7))
    
    res.sample.sorted <- res.sample[order(res.sample$qValue),]
    res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$qValue) & res.sample.sorted$log2FC < 0,]
    label.dn.sig <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]
    
    res.sample.sorted <- res.sample[order(res.sample$log2FC),]
    res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$qValue) & res.sample.sorted$log2FC < 0,]
    label.dn.fc <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]
    
    with(res.sample.sorted[res.sample.sorted$id %in% c(label.dn.sig, label.dn.fc),], textxy(log2FC, -log10(qValue), labs=Name, cex=cexLabel, offset=0.7))
  }

  legend(legendpos, xjust=1, yjust=1, legend=c(paste("qValue <= ",sigthresh,sep=""), paste("|Log2FC| >= ",lfcthresh,sep=""), "both"), pch=20, col=c("black","blue","red"), bg="white")
}

pdf("/mnt/projects/fikret/results/volcanoPlots.DTCvsTUM.DTCdxvsDTCrel.pdf")

volcanoPlot(expr="/mnt/projects/fikret/results/anduril/execute/deseq_DTC_vs_TUM/results.csv", 
            title="DTC vs. TUM", 
            sigthresh=0.01, 
            lfcthresh=1,
            labelTopN=20)

volcanoPlot(expr="/mnt/projects/fikret/results/anduril/execute/deseq_DTCdx_vs_DTCrel/results.csv", 
            title="DTCdx vs. DTCrel", 
            sigthresh=0.01, 
            lfcthresh=1,
            labelTopN=20)

dev.off()
