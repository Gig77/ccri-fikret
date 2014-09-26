export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=~/fikret
TRIM_BEFORE_BASE=1
COUNT_STAT_YMAX=60

SAMPLES=01-15-BM100	02-15-BM30 03-15-BM0 \
		04-6-BM100 05-6-BM30 06-6-BM0 \
		07-644-T 08-644-M 09-644-D \
		10-919-T 11-919-M 12-919-D \
		D-1992-0551 D-1992-1182 D-1992-1283 D-1993-0849 D-1993-1141 D-1994-1911 D-1995-0029 D-1995-1415 D-2013-0131 D-2014-1234 D2-1993-1141 M-1195-0029 \
		M-1955-1415 M-1992-0551 M-1992-1182 M-1992-1283 M-1993-0849 M-1993-1141 M-1994-1911 M-2013-0131 M-2014-1234 \
		T-1995-0029 \
		D-1992-0841 M-2006-3040 M-2010-2922 D-1998-1261 M-2007-0342 D-1995-2493 M-1992-0841 M-2011-2208 D-1998-1897 M-1998-1261 D-1998-0062 D-2010-0681 T-1998-126 \
		D-2002-0103 D-2003-5050 D-2006-3040 D-2010-2922 M-1998-1897 D-2007-0342 M-1998-0062 M-2003-5050 D-2011-2208 M-2002-0103 M-2010-0681 \
		D-1995-2120 M-1995-2120 D-1999-3072 M-1999-3072 D-2014-0613 M-2014-0613 D-2011-4162 M-2011-4162 D-2010-2326 M-2010-2326 D-1992-0825 M-1992-0825 \
		T-2003-4304 T2-2003-4304 D-2003-4304 M-2003-4304 D-1998-1987 M-1998-1987 D-1993-0612 M-1993-0612 D-1995-0614 M-1995-0614 D-2012-2162 M-2012-2162 \
		D-1995-0581 M-1995-0581 D-1995-2120 M-1995-2120 D-2013-2083 M-2013-2083
		 
all: gsnap htseq qc blast fastqc

sample-dist.pca.pdf sample-dist.heatmap.pdf: $(foreach S, $(SAMPLES), htseq/$S.count) ~/fikret/scripts/cluster-samples.R
	Rscript ~/fikret/scripts/cluster-samples.R
	mv sample-dist.pca.pdf.part sample-dist.pca.pdf
	mv sample-dist.heatmap.pdf.part sample-dist.heatmap.pdf

include ~/generic/scripts/rna-seq/gsnap.mk
include ~/generic/scripts/rna-seq/htseq.mk
include ~/generic/scripts/rna-seq/qc.mk
include ~/generic/scripts/rna-seq/blast.mk
include ~/generic/scripts/rna-seq/fastqc.mk

