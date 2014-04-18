#!/bin/sh

BT2_IDX_DIR=$HOME/bowtie2_indexes
BT2_IDX=$BT2_IDX_DIR/hg19.fa

if ! which bowtie2 ; then
    echo "Add bowtie2 to PATH first"
fi

$BT2_EXE \
	--sam-no-qname-trunc \
	-x $BT2_IDX \
	-U hg19.sim.fq,contam1.sim.fq,contam2.sim.fq,contam3.sim.fq \
	-S sim.sam \
	-p 12
