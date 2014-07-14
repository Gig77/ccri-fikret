export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

# build GMAP database
# cd ~/chrisi/data/gsnap
# ~/tools/gmap-2014-05-15/util/gmap_build -d g1k_v37 ~/generic/data/broad/human_g1k_v37.fasta
# cat ~/generic/data/broad/human_g1k_v37.fasta ~/chrisi/results/etv6-runx1.breakpoint.fa > ~/chrisi/data/gsnap/human_g1k_v37_etv6runx1.fasta 
# ~/tools/gmap-2014-05-15/bin/gmap_build --dir=/data/christian/chrisi/data/current/gsnap --db=g1k_v37_etv6runx1 ~/chrisi/data/gsnap/human_g1k_v37_etv6runx1.fasta
# rm ~/chrisi/data/gsnap/human_g1k_v37_etv6runx1.fasta
#
# gunzip -c snp138.txt.gz | ~/tools/gmap-2014-05-15/bin/dbsnp_iit -w 1 > g1k_v37.snp138.txt
# cd g1k_v37_etv6runx1
# cat g1k_v37.snp138.txt | ~/tools/gmap-2014-05-15/bin/iit_store -o g1k_v37.snp138
# ~/tools/gmap-2014-05-15/bin/snpindex -d g1k_v37_etv6runx1 -D g1k_v37_etv6runx1 -v snpfile138
# mv g1k_v37_etv6runx1/g1k_v37.snp138.iit g1k_v37_etv6runx1/g1k_v37_etv6runx1.maps

# for splicesites
# ?? gunzip -c refGene.txt.gz | psl_splicesites -s 1 > mm10splicesites #ok
# cat g1k_v37.splicesites | ~/tools/gmap-2014-05-15/bin/iit_store -o g1k_v37_etv6runx1.maps/g1k_v37.splicesites

SAMPLES=C3JP8ACXX_01-15-BM100_14s001408-1-1_Rifatbegovic_lane114s001408_sequence \
		C3JP8ACXX_02-15-BM30_14s001409-1-1_Rifatbegovic_lane114s001409_sequence \
		C3JP8ACXX_03-15-BM0_14s001410-1-1_Rifatbegovic_lane114s001410_sequence \
		C3JP8ACXX_04-6-BM100_14s001411-1-1_Rifatbegovic_lane114s001411_sequence \
		C3JP8ACXX_05-6-BM30_14s001412-1-1_Rifatbegovic_lane114s001412_sequence \
		C3JP8ACXX_06-6-BM0_14s001413-1-1_Rifatbegovic_lane114s001413_sequence \
		C3JP8ACXX_07-644-T_14s001414-1-1_Rifatbegovic_lane214s001414_sequence \
		C3JP8ACXX_08-644-M_14s001415-1-1_Rifatbegovic_lane214s001415_sequence \
		C3JP8ACXX_09-644-D_14s001416-1-1_Rifatbegovic_lane214s001416_sequence \
		C3JP8ACXX_10-919-T_14s001417-1-1_Rifatbegovic_lane214s001417_sequence \
		C3JP8ACXX_11-919-M_14s001418-1-1_Rifatbegovic_lane214s001418_sequence \
		C3JP8ACXX_12-919-D_14s001419-1-1_Rifatbegovic_lane214s001419_sequence

all: $(foreach S, $(SAMPLES), htseq/$S.count flagstat/$S.samtools.flagstat)

#----------------------------------------------
# BAM to FASTQ, align with GSNAP, sort & index
#----------------------------------------------
gsnap/%.gsnap.bam: ~/fikret/data/fastq/%.txt.gz
	~/tools/gmap-2014-05-15/src/gsnap \
			--db=g1k_v37_etv6runx1 \
			--dir=/data/christian/chrisi/data/current/gsnap/g1k_v37_etv6runx1 \
			--format=sam \
			--npaths=1 \
			--quiet-if-excessive \
			--nofails \
			--batch=4  \
			--quality-protocol=sanger \
			--print-snps \
			--nthreads=24 \
			--input-buffer-size=5000 \
			--use-splicing=g1k_v37.splicesites \
			--use-snps=g1k_v37.snp138 \
			--genome-unk-mismatch=0 \
			--gunzip \
			$< \
		| ~/tools/samtools-0.1.19/samtools view -Shb - \
		| ~/tools/samtools-0.1.19/samtools sort -m 3000000000 - $@ \
		2>&1 | $(LOG)
	mv $@.bam $@
	~/tools/samtools-0.1.19/samtools index $@ 2>&1 | $(LOG)

flagstat/%.samtools.flagstat: gsnap/%.gsnap.bam
	samtools flagstat $< 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

htseq/%.count: gsnap/%.gsnap.bam ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz
	~/tools/HTSeq-0.6.1/scripts/htseq-count -f bam -t exon -s no $< ~/chrisi/data/ensembl/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz 2>&1 1>$@.part | $(LOG)
	mv $@.part $@
	
#----
# RSeQC
#----

quality: rseqc/allpatients.rRNA.count rseqc/allpatients.read-distribution.txt rseqc/allpatients.splice_junction.pdf rseqc/allpatients.splicing_events.pdf rseqc/allpatients.geneBodyCoverage.pdf rseqc/allpatients.DupRate_plot.pdf

rseqc/allpatients.geneBodyCoverage.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.geneBodyCoverage.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/%.rseqc.geneBodyCoverage.pdf: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/geneBody_coverage.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed -o rseqc/$*.rseqc

rseqc/allpatients.DupRate_plot.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.DupRate_plot.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/%.rseqc.DupRate_plot.pdf: gsnap/%.gsnap.bam
	~/tools/RSeQC-2.3.9/bin/read_duplication.py -i $< -o rseqc/$*.rseqc

rseqc/allpatients.read-distribution.txt: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.read-distribution.txt)
	rm -f $@.part
	for S in $^ ; do echo $$S >> $@.part; cat $$S >> $@.part ; done
	mv $@.part $@

rseqc/%.rseqc.read-distribution.txt: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/read_distribution.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed > $@.part
	mv $@.part $@

rseqc/allpatients.splice_junction.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.splice_junction.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/allpatients.splicing_events.pdf: $(foreach S, $(SAMPLES), rseqc/$S.rseqc.splice_events.pdf)
	gs -dBATCH -dNOPAUSE -q -dAutoRotatePages=/None -sDEVICE=pdfwrite -sOutputFile=$@.part $^
	mv $@.part $@

rseqc/%.rseqc.splice_events.pdf rseqc/%.rseqc.splice_junction.pdf: gsnap/%.gsnap.bam ~/generic/data/rseqc/hg19_Ensembl.bed
	~/tools/RSeQC-2.3.9/bin/junction_annotation.py -i $< -r ~/generic/data/rseqc/hg19_Ensembl.bed -o rseqc/$*.rseqc

rseqc/allpatients.rRNA.count: $(foreach S, $(SAMPLES), rseqc/$S.rRNA.count)
	rm -f $@.part
	for S in $^ ; do echo $$S >> $@.part; cat $$S >> $@.part ; done
	mv $@.part $@

rseqc/%.rRNA.count: gsnap/%.gsnap.bam flagstat/%.samtools.flagstat
	echo "Total number of reads:" `head -1 flagstat/$*.samtools.flagstat | cut -f 1 -d ' '` > $@.part
	echo "Number of reads mapping to rRNA genes: " `samtools view $< -L ~/generic/data/ensembl/rRNA.ensembl.biomart.GRCh37.p13.bed | wc -l` >> $@.part
	mv $@.part $@	
