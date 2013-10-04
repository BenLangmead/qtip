#!/bin/sh

BT2_HOME=$HOME/git/bowtie2/bowtie2

REF=$HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa
IDX=$HOME/bowtie2_indexes/Homo_sapiens.GRCh37.60.dna.fa
UNP_IN=$HOME/reads/SRA_HISEQ2000_FC1.shuffle.2M.1.fastq

SAM_OUT=test_ts_human.sam
SAM_TSV_OUT=test_ts_human.sam.tsv
TRAIN_TSV_OUT=test_ts_human

python ../bin/ts.py \
	--ref $REF \
	--pickle-ref Homo_sapiens.GRCh37.60.dna.fa.pickle \
	--U $UNP_IN \
	--S $SAM_OUT \
	--fastq \
	--num-reads 50000 \
	--bt2-exe $BT2_HOME/bowtie2 \
	--bt2-args "-x $IDX" \
	--upto 50000 \
	--save-sam-tsv $SAM_TSV_OUT \
	--save-training-tsv $TRAIN_TSV_OUT
