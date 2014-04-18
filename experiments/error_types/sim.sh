#!/bin/sh

#
# Run Mason on the human genome to obtain 91% of the input data.
# Then run mason on the lambda, E. coli and mouse genomes to
# obtain 3% each of the input data.
#

if [ -z "$TS_HOME" ] ; then
    echo "Set TS_HOME first; should contain software/mason/mason"
    exit 1
fi

if [ -z "$TS_REFS" ] ; then
    echo "Set TS_REFS first; should contain hg19.fa, lambda_virus.fa, e_coli.fa and mm10.fa"
    exit 1
fi

HG19_FA=$TS_REFS/hg19.fa
CONTAM1_FA=$TS_REFS/lambda_virus.fa
CONTAM2_FA=$TS_REFS/e_coli.fa
CONTAM3_FA=$TS_REFS/mm10.fa

# 3% contamination from each of the three types
NUM_READS=910000
NUM_CONTAM_READS=30000

$TS_HOME/software/mason/mason illumina \
	--num-reads $NUM_READS \
	--num-haplotypes 2 \
	--include-read-information \
	--seed 345 \
	--simulate-qualities \
	--read-length 100 \
	--output-file hg19.sim.fq \
	$HG19_FA

$TS_HOME/software/mason/mason illumina \
	--num-reads $NUM_CONTAM_READS \
	--num-haplotypes 1 \
	--include-read-information \
	--seed 254 \
	--simulate-qualities \
	--read-length 100 \
	--output-file .contam1.sim.fq \
	$CONTAM1_FA

sed 's/contig=/contig=lambda-/' < .contam1.sim.fq > contam1.sim.fq
rm -f .contam1.sim.fq

$TS_HOME/software/mason/mason illumina \
	--num-reads $NUM_CONTAM_READS \
	--num-haplotypes 1 \
	--include-read-information \
	--seed 298 \
	--simulate-qualities \
	--read-length 100 \
	--output-file .contam2.sim.fq \
	$CONTAM2_FA

sed 's/contig=/contig=e_coli-/' < .contam2.sim.fq > contam2.sim.fq
rm -f .contam2.sim.fq

$TS_HOME/software/mason/mason illumina \
	--num-reads $NUM_CONTAM_READS \
	--num-haplotypes 2 \
	--include-read-information \
	--seed 046 \
	--simulate-qualities \
	--read-length 100 \
	--output-file .contam3.sim.fq \
	$CONTAM3_FA

sed 's/contig=/contig=mouse-/' < .contam3.sim.fq > contam3.sim.fq
rm -f .contam3.sim.fq
