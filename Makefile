export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=/mnt/projects/fikret
TRIM_BEFORE_BASE=1
COUNT_STAT_YMAX=60

SAMPLES=03-15-BM0 02-15-BM30 01-15-BM100 \
		06-6-BM0 05-6-BM30 04-6-BM100 \
		T-93-0644-06-Dx D-93-0644-06-Dx M-93-0644-06-Dx \
		T-94-0919-05-Dx D-94-0919-05-Dx M-94-0919-05-Dx \
		T-98-1261-02-Dx D-98-1261-02-Dx M-98-1261-02-Dx \
		T-95-0029-01-Dx D-95-0029-01-Dx M-95-0029-01-Dx \
		T-03-4304-03-Dx T2-03-4304-03-Dx D-03-4304-03-Dx M-03-4304-03-Dx \
		T-99-1650-04-Dx T2-99-1650-04-Dx D-99-1650-04-Dx M-99-1650-04-Dx \
		D-92-0551-07-Dx M-92-0551-07-Dx \
		D-93-0849-10-R M-93-0849-10-R \
		D-93-1141-07-R D2-93-1141-07-R M-93-1141-07-R \
		D-13-0131-14-R M-13-0131-14-R \
		D-14-1234-12-R M-14-1234-12-R \
		D-92-0841-17-R M-92-0841-17-R-x \
		D-06-3040-29-Dx M-06-3040-29-Dx \
		D-10-2922-23-R M-10-2922-23-R \
		D-07-0342-19-R M-07-0342-19-R \
		D-11-2208-18-R M-11-2208-18-R \
		D-98-1897-13-R M-98-1897-13-R \
		D-98-0062-22-R M-98-0062-22-R \
		D-10-0681-28-Dx M-10-0681-28-Dx \
		D-02-0103-21-R M-02-0103-21-R \
		D-95-2120-06-R M-95-2120-06-R D-95-2120-06-R-2 M-95-2120-06-R-2 \
		D-99-3072-24-R M-99-3072-24-R \
		D-10-2326-20-R M-10-2326-20-R \
		D-92-0825-25-R M-92-0825-25-R \
		D-98-1987-30-Dx M-98-1987-30-Dx \
		D-93-0612-08-Dx M-93-0612-08-Dx \
		D-95-0614-08-R M-95-0614-08-R \
		D-12-2162-16-R M-12-2162-16-R \
		D-95-0581-27-Dx M-95-0581-27-Dx \
		D-13-2083-15-R M-13-2083-15-R \
		D-09-0606-11-R M-09-0606-11-R \
		\
		D-14-0613-26-R M-14-0613-26-R-x \
		D-94-1911-33-R-x M-94-1911-33-R \
		D-92-1283-10-Dx-x M-92-1283-10-Dx \
		D-92-1182-33-Dx-x M-92-1182-33-Dx \
		D-95-1415-31-Dx-x M-95-1415-31-Dx \
		D-03-5050-32-Dx-x M-03-5050-32-Dx-x \
		D-11-4162-34-R-x M-11-4162-34-R \
		D-95-2493-35-Dx-x


all: all-default deseq sample-dist.pca.without-outliers.pdf

include /mnt/projects/generic/scripts/rna-seq/default.mk

clean: clean-default clean-deseq
	rm -f sample-dist.*.pdf

#=========================
# sample clustering
#=========================
 
sample-dist.pca.without-outliers.pdf: $(foreach S, $(SAMPLES), htseq/$S.count) /mnt/projects/fikret/scripts/cluster-samples.R
	Rscript /mnt/projects/fikret/scripts/cluster-samples.R
	mv sample-dist.pca.without-outliers.pdf.part sample-dist.pca.without-outliers.pdf
	mv sample-dist.heatmap.without-outliers.pdf.part sample-dist.heatmap.without-outliers.pdf
	mv sample-dist.pca.with-outliers.pdf.part sample-dist.pca.with-outliers.pdf
	mv sample-dist.heatmap.with-outliers.pdf.part sample-dist.heatmap.with-outliers.pdf

#=========================
# DEG analysis with DeSeq2
#=========================
.PHONY: deseq clean-deseq

deseq: deseq/TU-vs-DTC.diagnosis.paired.tsv deseq/TU-vs-DTC.diagnosis.full.tsv deseq/TU-vs-BM.diagnosis.paired.tsv deseq/TU-vs-BM.diagnosis.full.tsv deseq/DTC-vs-BM.diagnosis.tsv deseq/DTC-dia-vs-DTC-rel.tsv deseq/BM-dia-vs-BM-rel.tsv deseq/DTC-vs-BM.relapse.tsv deseq/DTC-rel-vs-BM-dia.tsv

clean-deseq:
	rm deseq/*
	rmdir deseq

deseq/TU-vs-DTC.diagnosis.paired.tsv: /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control    D-99-1650-04-Dx.count,D-94-0919-05-Dx.count,D-93-0644-06-Dx.count,D-95-0029-01-Dx.count,D-03-4304-03-Dx.count \
		--experiment T-99-1650-04-Dx.count,T2-99-1650-04-Dx.count,T-94-0919-05-Dx.count,T-93-0644-06-Dx.count,T-95-0029-01-Dx.count,T-03-4304-03-Dx.count,T2-03-4304-03-Dx.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/TU-vs-DTC.diagnosis.full.tsv: /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control    D-99-1650-04-Dx.count,D-94-0919-05-Dx.count,D-93-0644-06-Dx.count,D-95-0029-01-Dx.count,D-03-4304-03-Dx.count,D-92-0551-07-Dx.count,D-93-0612-08-Dx.count,D-10-0681-28-Dx.count,D-06-3040-29-Dx.count,D-98-1987-30-Dx.count,D-98-1261-02-Dx.count,D-95-0581-27-Dx.count \
		--experiment T-99-1650-04-Dx.count,T2-99-1650-04-Dx.count,T-94-0919-05-Dx.count,T-93-0644-06-Dx.count,T-95-0029-01-Dx.count,T-03-4304-03-Dx.count,T2-03-4304-03-Dx.count,T-98-1261-02-Dx.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/TU-vs-BM.diagnosis.paired.tsv: /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control    M-99-1650-04-Dx.count,M-94-0919-05-Dx.count,M-93-0644-06-Dx.count,M-95-0029-01-Dx.count,M-03-4304-03-Dx.count \
		--experiment T-99-1650-04-Dx.count,T2-99-1650-04-Dx.count,T-94-0919-05-Dx.count,T-93-0644-06-Dx.count,T-95-0029-01-Dx.count,T-03-4304-03-Dx.count,T2-03-4304-03-Dx.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/TU-vs-BM.diagnosis.full.tsv: /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control    M-99-1650-04-Dx.count,M-94-0919-05-Dx.count,M-93-0644-06-Dx.count,M-95-0029-01-Dx.count,M-03-4304-03-Dx.count,M-92-0551-07-Dx.count,M-93-0612-08-Dx.count,M-10-0681-28-Dx.count,M-06-3040-29-Dx.count,M-98-1987-30-Dx.count,M-98-1261-02-Dx.count,M-95-0581-27-Dx.count \
		--experiment T-99-1650-04-Dx.count,T2-99-1650-04-Dx.count,T-94-0919-05-Dx.count,T-93-0644-06-Dx.count,T-95-0029-01-Dx.count,T-03-4304-03-Dx.count,T2-03-4304-03-Dx.count,T-98-1261-02-Dx.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/DTC-vs-BM.diagnosis.tsv: /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control    M-99-1650-04-Dx.count,M-94-0919-05-Dx.count,M-93-0644-06-Dx.count,M-95-0029-01-Dx.count,M-03-4304-03-Dx.count,M-92-0551-07-Dx.count,M-93-0612-08-Dx.count,M-10-0681-28-Dx.count,M-06-3040-29-Dx.count,M-98-1987-30-Dx.count,M-98-1261-02-Dx.count,M-95-0581-27-Dx.count \
		--experiment D-99-1650-04-Dx.count,D-94-0919-05-Dx.count,D-93-0644-06-Dx.count,D-95-0029-01-Dx.count,D-03-4304-03-Dx.count,D-92-0551-07-Dx.count,D-93-0612-08-Dx.count,D-10-0681-28-Dx.count,D-06-3040-29-Dx.count,D-98-1987-30-Dx.count,D-98-1261-02-Dx.count,D-95-0581-27-Dx.count \
		--output-tsv $@.part
	mv $@.part $@
	
deseq/DTC-dia-vs-DTC-rel.tsv:  /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control D-99-1650-04-Dx.count,D-95-0029-01-Dx.count,D-03-4304-03-Dx.count,D-94-0919-05-Dx.count,D-93-0644-06-Dx.count,D-92-0551-07-Dx.count,D-93-0612-08-Dx.count,D-10-0681-28-Dx.count,D-06-3040-29-Dx.count,D-98-1987-30-Dx.count,D-98-1261-02-Dx.count,D-95-0581-27-Dx.count \
		--experiment D-95-2120-06-R.count,D-95-2120-06-R-2.count,D-93-1141-07-R.count,D2-93-1141-07-R.count,D-95-0614-08-R.count,D-09-0606-11-R.count,D-14-1234-12-R.count,D-98-1897-13-R.count,D-13-0131-14-R.count,D-13-2083-15-R.count,D-12-2162-16-R.count,D-92-0841-17-R.count,D-11-2208-18-R.count,D-07-0342-19-R.count,D-10-2326-20-R.count,D-02-0103-21-R.count,D-98-0062-22-R.count,D-10-2922-23-R.count,D-99-3072-24-R.count,D-92-0825-25-R.count,D-93-0849-10-R.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/BM-dia-vs-BM-rel.tsv:  /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control M-99-1650-04-Dx.count,M-94-0919-05-Dx.count,M-93-0644-06-Dx.count,M-95-0029-01-Dx.count,M-03-4304-03-Dx.count,M-92-0551-07-Dx.count,M-93-0612-08-Dx.count,M-10-0681-28-Dx.count,M-06-3040-29-Dx.count,M-98-1987-30-Dx.count,M-98-1261-02-Dx.count,M-95-0581-27-Dx.count \
		--experiment M-95-2120-06-R.count,M-95-2120-06-R-2.count,M-93-1141-07-R.count,M-95-0614-08-R.count,M-09-0606-11-R.count,M-14-1234-12-R.count,M-98-1897-13-R.count,M-13-0131-14-R.count,M-13-2083-15-R.count,M-12-2162-16-R.count,M-11-2208-18-R.count,M-07-0342-19-R.count,M-10-2326-20-R.count,M-02-0103-21-R.count,M-98-0062-22-R.count,M-10-2922-23-R.count,M-99-3072-24-R.count,M-92-0825-25-R.count,M-93-0849-10-R.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/DTC-vs-BM.relapse.tsv: /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control    M-95-2120-06-R.count,M-95-2120-06-R-2.count,M-93-1141-07-R.count,M-95-0614-08-R.count,M-09-0606-11-R.count,M-14-1234-12-R.count,M-98-1897-13-R.count,M-13-0131-14-R.count,M-13-2083-15-R.count,M-12-2162-16-R.count,M-11-2208-18-R.count,M-07-0342-19-R.count,M-10-2326-20-R.count,M-02-0103-21-R.count,M-98-0062-22-R.count,M-10-2922-23-R.count,M-99-3072-24-R.count,M-92-0825-25-R.count,M-93-0849-10-R.count \
		--experiment D-95-2120-06-R.count,D-95-2120-06-R-2.count,D-93-1141-07-R.count,D2-93-1141-07-R.count,D-95-0614-08-R.count,D-09-0606-11-R.count,D-14-1234-12-R.count,D-98-1897-13-R.count,D-13-0131-14-R.count,D-13-2083-15-R.count,D-12-2162-16-R.count,D-92-0841-17-R.count,D-11-2208-18-R.count,D-07-0342-19-R.count,D-10-2326-20-R.count,D-02-0103-21-R.count,D-98-0062-22-R.count,D-10-2922-23-R.count,D-99-3072-24-R.count,D-92-0825-25-R.count,D-93-0849-10-R.count \
		--output-tsv $@.part
	mv $@.part $@

deseq/DTC-rel-vs-BM-dia.tsv:  /mnt/projects/generic/scripts/rna-seq/diff-exp.R $(foreach S, $(SAMPLES), htseq/$S.count)
	mkdir -p deseq
	Rscript /mnt/projects/generic/scripts/rna-seq/diff-exp.R \
		--control M-99-1650-04-Dx.count,M-94-0919-05-Dx.count,M-93-0644-06-Dx.count,M-95-0029-01-Dx.count,M-03-4304-03-Dx.count,M-92-0551-07-Dx.count,M-93-0612-08-Dx.count,M-10-0681-28-Dx.count,M-06-3040-29-Dx.count,M-98-1987-30-Dx.count,M-98-1261-02-Dx.count,M-95-0581-27-Dx.count \
		--experiment D-95-2120-06-R.count,D-95-2120-06-R-2.count,D-93-1141-07-R.count,D2-93-1141-07-R.count,D-95-0614-08-R.count,D-09-0606-11-R.count,D-14-1234-12-R.count,D-98-1897-13-R.count,D-13-0131-14-R.count,D-13-2083-15-R.count,D-12-2162-16-R.count,D-92-0841-17-R.count,D-11-2208-18-R.count,D-07-0342-19-R.count,D-10-2326-20-R.count,D-02-0103-21-R.count,D-98-0062-22-R.count,D-10-2922-23-R.count,D-99-3072-24-R.count,D-92-0825-25-R.count,D-93-0849-10-R.count \
		--output-tsv $@.part
	mv $@.part $@
	