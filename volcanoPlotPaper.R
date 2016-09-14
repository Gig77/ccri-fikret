volcanoPlot <- function
(
    expr, sigthresh=0.001, lfcthresh=2, fccol="fc", pcol="p", qcol="q", geneNameCol=1, labelTopN=30, 
    legendpos="bottomright", excludeGenes = "", nonSigSampleSize=1000, cexLabel=0.4, cexPoint=1, 
    minQ=1e-200, labelNonameGenesWithId=TRUE, onlyGenesWithNames=FALSE, title="", labelGenes = ""
) 
{
  # get input data and parameters
  expr <- read.delim(expr, stringsAsFactors = F, check.names = F)
  
  # prepare data
  res <- data.frame(id=expr[,1], log2FC=expr[,fccol], pvalue=expr[,pcol], qValue=pmax(expr[,qcol], minQ), Name=expr[,geneNameCol], stringsAsFactors = F)
  res <- res[!is.na(res$qValue) & !is.na(res$log2FC),]
  if (onlyGenesWithNames) {
    res <- res[!is.na(res$Name) & res$Name != "",]
  }
  if (labelNonameGenesWithId) {
    res$Name[is.na(res$Name) | res$Name == ""] <- res$id[is.na(res$Name) | res$Name == ""]
  }

  # subsample uninteresting fraction of genes to reduce plot size
  res.interesting   <- subset(res, qValue <= sigthresh | abs(log2FC) >= lfcthresh)
  res.uninteresting <- subset(res, qValue > sigthresh & abs(log2FC) < lfcthresh)
  res.uninteresting.sample <- res.uninteresting[sample(1:nrow(res.uninteresting), min(nonSigSampleSize, nrow(res.uninteresting)), replace = F),]
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
  if (labelGenes != "") {
    labelGenes <- unlist(strsplit(labelGenes, ","))
    with(res[res$Name %in% labelGenes,], textxy(log2FC, -log10(qValue), labs=Name, cex=cexLabel, offset=0.7))
  } else if (labelTopN > 0) {
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

#-----------

pdf("/mnt/projects/fikret/results/volcanoPlot.DTCdx_vs_DTCrel.forPaper.pdf")

volcanoPlot(expr="/mnt/projects/fikret/results/anduril/execute/deseqAnnotated_DTCdx_vs_DTCrel/table.csv", 
            title="DTCdx vs. DTCrel", 
            sigthresh=0.2, 
            lfcthresh=0.5,
            geneNameCol="Gene",
            labelGenes="LPPR3,BBC3,POLR2E,ZNF428,SIRT6,ADAMTS10,TMEM259,GPX4,CACNG7,CSNK1G2,STK11,ANKRD24,TMEM160,C19orf43,ALKBH7,ICAM5,OAZ1,ZNF414,SLC25A41,ARRDC2,PVRL2,RPL18,CADM4,POLRMT,TDRD12,GLTSCR2,MAP2K2,RPS19,TRIP10,ZNF444,GRAMD1A,RNF126,SLC39A3,MAGEA4,MAGEB2,MAGEA1,MAGEC2,MAGEA11,MAGEB1,MAGEA8,MAGEA9")

dev.off()

#-----------

pdf("/mnt/projects/fikret/results/volcanoPlot.DTC_vs_TUM.forPaper.pdf")

volcanoPlot(expr="/mnt/projects/fikret/results/anduril/execute/deseqAnnotated_DTC_vs_TUM/table.csv", 
            title="DTC vs. TUM", 
            sigthresh=0.001, 
            lfcthresh=2,
            geneNameCol="Gene",
            labelGenes="LPL,MSX1,LRPAP1,JAG2,APP,CCND2,JAG1,PDGFA,THBD,SPP1,STC1,KCNJ8,TIMP1,COL5A2,FSTL1,SLCO2A1,LUM,COL3A1,POSTN")

dev.off()

#-----------

pdf("/mnt/projects/fikret/results/volcanoPlot.DTC_vs_TUM.forPaper_2ndplot.pdf")

volcanoPlot(expr="/mnt/projects/fikret/results/anduril/execute/deseqAnnotated_DTC_vs_TUM/table.csv", 
            title="DTC vs. TUM", 
            sigthresh=0.001, 
            lfcthresh=2,
            geneNameCol="Gene",
            labelGenes="MGP,COL6A3,ACTA2,PRRX1,TAGLN,ECM1")

dev.off()

