#!/bin/sh

GAP=30
NUM_INPUT_READS=500000
NUM_TRAINING_READS=100000
BT2_HOME=$(HOME)/git/bowtie2/bowtie2
BT2_EXE=$(BT2_HOME)/bowtie2
REF=$(HOME)/bowtie2_indexes/hg19.fa
FA=$(HOME)/bowtie2_indexes/hg19.fa
BT2_ARGS=-p 8 -x $(REF) --sensitive
TS_PARAMS=--ref $(REF) --pickle-ref $(REF).pickle --num-reads $(NUM_TRAINING_READS)
BIN=../../bin
TS=python $(BIN)/ts.py

.PHONY: reads
reads: r_0/r_0.fq r_12/r_1.fq r_12/r_2.fq

.PHONY: clean
clean:
	rm -rf r_0 r_12

.PHONY: clean_training
clean_training:
	rm -f r_0/input.sam r_0/training.*
	rm -f r_12/input.sam r_12/training.*

r_0/r_0.fq: $(FA)
	mkdir -p r_0
	mason illumina -hn 2 -i -s 322 -sq -n 100 -N $(NUM_INPUT_READS) -o $@ $(FA)

r_12/r_1.fq r12/r_2.fq: $(FA)
	mkdir -p r_12
	mason illumina -hn 2 -i -s 535 -sq -mp -rn 2 -ll 375 -le 100 -n 100 -N $(NUM_INPUT_READS) -o r_12/r.fq $(FA)

r_0/training.pickle: r_0/r_0.fq
	$(TS) $(TS_PARAMS) --input-from-mason --bt2-exe $(BT2_EXE) --write-all --U $< --out r_0 -- $(BT2_ARGS)

r_12/training.pickle: r_12/r_1.fq r_12/r_2.fq
	$(TS) $(TS_PARAMS) --input-from-mason --bt2-exe $(BT2_EXE) --write-all --m1 r_12/r_1.fq --m2 r_12/r_2.fq --out r_12 -- $(BT2_ARGS)
