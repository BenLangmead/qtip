#!/bin/sh

#
# Run Bowtie 2 in --sensitive and --very-sensitive modes to align reads
# from the simulations.
#

if [ -z "$TS_HOME" ] ; then
    echo "Set TS_HOME first; should contain software/bowtie2/bowtie2"
    exit 1
fi
if [ -z "$TS_INDEXES" ] ; then
    echo "Set TS_INDEXES first; should contain hg19.fa.*.bt2"
    exit 1
fi

$TS_HOME/software/bowtie2/bowtie2 \
	--sam-no-qname-trunc \
	-x $TS_INDEXES/hg19.fa \
	-U hg19.sim.fq,contam1.sim.fq,contam2.sim.fq,contam3.sim.fq \
	--sensitive \
	-S sim.sensitive.sam \
	-p 12

$TS_HOME/software/bowtie2/bowtie2 \
	--sam-no-qname-trunc \
	-x $TS_INDEXES/hg19.fa \
	-U hg19.sim.fq,contam1.sim.fq,contam2.sim.fq,contam3.sim.fq \
	--very-sensitive \
	-S sim.very-sensitive.sam \
	-p 12
