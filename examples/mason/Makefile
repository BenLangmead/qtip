#!/bin/sh

TS_PARAMS=--ref $(REF) --pickle-ref $(REF).pickle --num-reads 200000
GAP=30
BT2_HOME=~/git/bowtie2/bowtie2
BT2_EXE=$(BT2_HOME)/bowtie2
REF=~/bowtie2_indexes/hg19.fa
FA=~/bowtie2_indexes/hg19.fa
BT2_ARGS=-p 8
REP=../../data/rep/hg19.fa.out.gz

all: sats

.PHONY: clean
clean:
	rm -f *.sac *.sat *.sam *.roc *.err *.training.tsv *.sam.tsv *.training *.pdf

%.sac: %.sam
	perl ../../bin/correct.pl --no-chi --gap=$(GAP) -- $< /dev/null $@

%.sat: %.sac
	perl ../../bin/tabulate_sam.pl $< | \
		python $(HOME)/Documents/workspace/toolbox/python/bin/rep.py --rep=$(REP) --sat-in=/dev/stdin --sat-out=$@

define ill_reads

# Tabulated alignments
sats: \
    r0_$1.bt2_vs.sat  r0_$1.bt2_s.sat  r0_$1.bt2_f.sat  r0_$1.bt2_vf.sat \
    r12_$1.bt2_vs.sat r12_$1.bt2_s.sat r12_$1.bt2_f.sat r12_$1.bt2_vf.sat
#    r0_$1.bt2_lvs.sat  r0_$1.bt2_ls.sat  r0_$1.bt2_lf.sat  r0_$1.bt2_lvf.sat \
#    r12_$1.bt2_lvs.sat r12_$1.bt2_ls.sat r12_$1.bt2_lf.sat r12_$1.bt2_lvf.sat

# Reads

r0_$1.fq: $$(FA)
	mason illumina -hn 2 -i -s 322 -sq -n $2 -N $3 -o .$$@ $$(FA)
	perl ../../bin/mason2art.pl .$$@ > $$@

# Alignments at various levels of sensitivity

r0_$1.bt2_vs.sam: r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_vs.sam.tsv \
      --save-training-tsv r0_$1.bt2_vs \
      --U $$< \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) \
      --very-sensitive"

r0_$1.bt2_lvs.sam: r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_vs.sam.tsv \
      --save-training-tsv r0_$1.bt2_vs \
      --U $$< \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --very-sensitive-local"

r0_$1.bt2_s.sam:  r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_s.sam.tsv \
      --save-training-tsv r0_$1.bt2_s \
      --U $$< \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --sensitive"

r0_$1.bt2_ls.sam:  r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_s.sam.tsv \
      --save-training-tsv r0_$1.bt2_s \
      --U $$< \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --sensitive-local"

r0_$1.bt2_f.sam:  r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_f.sam.tsv \
      --save-training-tsv r0_$1.bt2_f \
	  --U $$< \
	  --S $$@ \
      --no-mapq \
	  --bt2-exe $$(BT2_EXE) \
	  --bt2-args "$$(BT2_ARGS) -x $$(FA) --fast"

r0_$1.bt2_lf.sam:  r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_f.sam.tsv \
      --save-training-tsv r0_$1.bt2_f \
	  --U $$< \
	  --S $$@ \
      --no-mapq \
	  --bt2-exe $$(BT2_EXE) \
	  --bt2-args "$$(BT2_ARGS) -x $$(FA) --fast-local"

r0_$1.bt2_vf.sam: r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_vf.sam.tsv \
      --save-training-tsv r0_$1.bt2_vf \
      --U $$< \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --very-fast"

r0_$1.bt2_lvf.sam: r0_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r0_$1.bt2_vf.sam.tsv \
      --save-training-tsv r0_$1.bt2_vf \
      --U $$< \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --very-fast-local"

r1_$1.fq: $(FA)
	mason illumina -hn 2 -i -s 535 -sq -mp -rn 2 -ll 375 -le 100 -n $2 -N $3 -o .tmppair.fq $$(FA)
	perl ../../bin/mason2art.pl .tmppair_1.fq > $$@
	perl ../../bin/mason2art.pl .tmppair_2.fq > $$(@:r1_%=r2_%)

r12_$1.bt2_vs.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_vs.sam.tsv \
      --save-training-tsv r12_$1.bt2_vs \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --very-sensitive"

r12_$1.bt2_lvs.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_vs.sam.tsv \
      --save-training-tsv r12_$1.bt2_vs \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --very-sensitive-local"

r12_$1.bt2_s.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_s.sam.tsv \
      --save-training-tsv r12_$1.bt2_s \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --sensitive"

r12_$1.bt2_ls.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_s.sam.tsv \
      --save-training-tsv r12_$1.bt2_s \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) --sensitive-local"

r12_$1.bt2_f.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_f.sam.tsv \
      --save-training-tsv r12_$1.bt2_f \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) \
      --fast"

r12_$1.bt2_lf.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_f.sam.tsv \
      --save-training-tsv r12_$1.bt2_f \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) \
      --fast-local"

r12_$1.bt2_vf.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_vf.sam.tsv \
      --save-training-tsv r12_$1.bt2_vf \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) \
      --very-fast"

r12_$1.bt2_lvf.sam: r1_$1.fq $$(FA).1.bt2
	python ../../bin/ts.py $$(TS_PARAMS) \
      --save-sam-tsv r12_$1.bt2_vf.sam.tsv \
      --save-training-tsv r12_$1.bt2_vf \
      --m1 $$< \
      --m2 $$(<:r1_%=r2_%) \
      --S $$@ \
      --no-mapq \
      --bt2-exe $$(BT2_EXE) \
      --bt2-args "$$(BT2_ARGS) -x $$(FA) \
      --very-fast-local"

endef

$(eval $(call ill_reads,ill_mason_50_500k,50,500000))
$(eval $(call ill_reads,ill_mason_75_500k,75,500000))
$(eval $(call ill_reads,ill_mason_100_500k,100,500000))
$(eval $(call ill_reads,ill_mason_125_500k,125,500000))
$(eval $(call ill_reads,ill_mason_150_500k,150,500000))
$(eval $(call ill_reads,ill_mason_175_500k,175,500000))
$(eval $(call ill_reads,ill_mason_200_500k,200,500000))
$(eval $(call ill_reads,ill_mason_225_500k,225,500000))
$(eval $(call ill_reads,ill_mason_250_500k,250,500000))
$(eval $(call ill_reads,ill_mason_275_500k,275,500000))
$(eval $(call ill_reads,ill_mason_300_500k,300,500000))
