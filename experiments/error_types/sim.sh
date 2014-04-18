#!/bin/sh

#
# Run Mason on the human genome to obtain 91% of the input data.
# Then run mason on the lambda, E. coli and mouse genomes to
# obtain 3% each of the input data.
#

FASTA_DIR=$HOME/fasta

HG19_FA=$FASTA_DIR/hg19.fa
CONTAM1_FA=$FASTA_DIR/lambda_virus.fa
CONTAM2_FA=$FASTA_DIR/e_coli.fa
CONTAM3_FA=$FASTA_DIR/mm10.fa

# 3% contamination from each of the three types
NUM_READS=910000
NUM_CONTAM_READS=30000

if ! which mason ; then
	echo "Add mason to path first"
fi

mason illumina \
	--num-reads $NUM_READS \
	--num-haplotypes 2 \
	--include-read-information \
	--seed 345 \
	--simulate-qualities \
	--read-length 100 \
	--output-file hg19.sim.fq \
	$HG19_FA

mason illumina \
	--num-reads $NUM_CONTAM_READS \
	--num-haplotypes 1 \
	--include-read-information \
	--seed 254 \
	--simulate-qualities \
	--read-length 100 \
	--output-file contam1.sim.fq \
	$CONTAM1_FA

mason illumina \
	--num-reads $NUM_CONTAM_READS \
	--num-haplotypes 1 \
	--include-read-information \
	--seed 298 \
	--simulate-qualities \
	--read-length 100 \
	--output-file contam2.sim.fq \
	$CONTAM2_FA

mason illumina \
	--num-reads $NUM_CONTAM_READS \
	--num-haplotypes 2 \
	--include-read-information \
	--seed 046 \
	--simulate-qualities \
	--read-length 100 \
	--output-file contam3.sim.fq \
	$CONTAM3_FA
