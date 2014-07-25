export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=~/fikret
TRIM_BEFORE_BASE=1

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

all: gsnap htseq qc blast

include ~/generic/scripts/rna-seq/gsnap.mk
include ~/generic/scripts/rna-seq/htseq.mk
include ~/generic/scripts/rna-seq/qc.mk
include ~/generic/scripts/rna-seq/blast.mk

