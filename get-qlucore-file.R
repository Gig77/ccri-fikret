# read input data
m <- read.delim("/mnt/projects/fikret/results/anduril/execute/qcReport-samplesClusterHeatmap/vst.csv", check.names = F, stringsAsFactors = F)
ann.var <- read.delim("/mnt/projects/fikret/data/ensembl.homo_sapiens_75_37.geneAnnotations.tsv", check.names = F, stringsAsFactors = F)
ann.sample <- read.delim("/mnt/projects/fikret/data/sample_key.csv", check.names = F, stringsAsFactors = F)
groups <- read.delim("/mnt/projects/fikret/results/anduril/execute/sampleGroupsWithoutExcludedSamples/table.csv", check.names = F, stringsAsFactors = F)

# add group memberships as columns to sample annotation
for(g in groups$ID) {
  ann.sample[,g] <- NA
  members <- unlist(strsplit(groups$Members[groups$ID==g], ","))
  ann.sample[,g] <- ifelse(ann.sample$Alias %in% members, "yes", "no")
}

# add collapsed GO terms to variable annotation
ensembl2go.file <- "/mnt/projects/schwann/results/ensembl2go.RData"
if(file.exists(ensembl2go.file)) {
	load(ensembl2go.file)
} else {
	library(biomaRt)
	mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
	ensembl2go <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006", "definition_1006", "go_linkage_type"), mart = mart)
	save(ensembl2go, file=ensembl2go.file)
}
ensembl2go <- ensembl2go[ensembl2go$go_id != "",]
ensembl2go$id_term <- paste(ensembl2go$go_id, ensembl2go$name_1006)
ensembl2go.agg <- ensembl2go[,c("ensembl_gene_id", "id_term")]
ensembl2go.agg <- aggregate(id_term~ensembl_gene_id, paste, collapse=", ", data=ensembl2go.agg)
names(ensembl2go.agg) <- c("Ensembl Gene ID", "GO")
ann.var <- merge(ann.var, ensembl2go.agg, all.x = T)
ann.var.names <- c("Ensembl", "Description", "Band", "HGNC", "Chr", "Start", "End", "GO")
names(ann.var) <- ann.var.names

# annotate expression matrix
sample.names <- colnames(m)[-1]
m.ann <- merge(m, ann.var[,ann.var.names], all.x = T)
m.ann <- m.ann[,c(ann.var.names, sample.names)]

# combine sample annotation with expression matrix into single data frame
ann.sample <- ann.sample[match(sample.names, ann.sample$Alias),]
ann.sample.transposed <- as.data.frame(t(ann.sample), stringsAsFactors = F)
colnames(ann.sample.transposed) <- ann.sample$Alias
ann.sample.transposed$Attribute <- rownames(ann.sample.transposed)

qlucore <- merge(m.ann, ann.sample.transposed, by.x=c("Ensembl", sample.names), by.y=c("Attribute", sample.names), all=T)
qlucore <- qlucore[order(qlucore$Description, na.last=F),c(ann.var.names, sample.names)]
qlucore <- rbind(qlucore[qlucore$Ensembl %in% ann.sample.transposed$Attribute,], c(ann.var.names, rep(NA, length(sample.names))), qlucore[!qlucore$Ensembl %in% ann.sample.transposed$Attribute,])
qlucore <- cbind(qlucore[,ann.var.names], data.frame(Annotation=NA), qlucore[,sample.names])
qlucore$Annotation <- as.character(ifelse(qlucore$Ensembl %in% ann.sample.transposed$Attribute, qlucore$Ensembl, NA))
qlucore$Ensembl <- ifelse(qlucore$Ensembl %in% ann.sample.transposed$Attribute, NA, qlucore$Ensembl)

# put annotation with sample id and alias as first row
qlucore <- qlucore[order(qlucore$Annotation %in% c("ID", "Alias"), decreasing = T),]

write.table(qlucore, "/mnt/projects/fikret/results/qlucore/NB-DTC.qlucore.expression-matrix.tsv", quote=F, sep="\t", col.names=F, row.names=F)
