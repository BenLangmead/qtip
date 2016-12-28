MASON_DIR=../../software/mason
MASON=$(MASON_DIR)/mason
BT2_DIR=../../software/bowtie2
BT2=$(BT2_DIR)/bowtie2
MASON_CONVERT=python ../../bin/mason_convert.py
TOOL=qtip
NTHREADS=1

all: full_e2e/final.sam \
     extra_e2e/trial0/final.sam \
     final_full_e2e.sam \
     test_random_forest \
     test_gradient_boost \
     test_reweight_mapq \
     test_reweight_ratio \
     test_max_leaf_nodes

final_full_e2e.sam: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq \
		-- -p $(NTHREADS) --reorder

full_e2e/final.sam: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--output-directory full_e2e \
		--keep-intermediates \
		--write-orig-mapq \
		--write-precise-mapq \
		--assess-accuracy \
		--threads $(NTHREADS) \
		-- -p $(NTHREADS) --reorder

extra_e2e/trial0/final.sam: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--output-directory extra_e2e \
		--keep-intermediates \
		--predict-for-training \
		--try-include-mapq \
		--subsampling-series 0.1,0.5,1.0 \
		--trials 2 \
		--write-orig-mapq \
		--write-precise-mapq \
		--assess-accuracy \
		--collapse \
		-- -p $(NTHREADS) --reorder

#
# Try a few different machine learning models
#
test_random_forest: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--output-directory $@ \
		--keep-intermediates \
		--predict-for-training \
		--try-include-mapq \
		--subsampling-series 0.1,0.5,1.0 \
		--trials 2 \
		--write-orig-mapq \
		--write-precise-mapq \
		--assess-accuracy \
		--model-family RandomForest \
		--threads $(NTHREADS) \
		-- -p $(NTHREADS) --reorder

test_max_leaf_nodes: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--output-directory $@ \
		--model-family RandomForest \
		--max-leaf-nodes 20 \
		--threads $(NTHREADS) \
		-- -p $(NTHREADS) --reorder

test_reweight_mapq: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--output-directory $@ \
		--reweight-mapq-offset 20.0 \
		--reweight-mapq \
		--assess-accuracy \
		--model-family RandomForest \
		--threads $(NTHREADS) \
		-- -p $(NTHREADS) --reorder

test_reweight_ratio: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--output-directory $@ \
		--reweight-ratio 10.0 \
		--assess-accuracy \
		--model-family RandomForest \
		--threads $(NTHREADS) \
		-- -p $(NTHREADS) --reorder

test_gradient_boost: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../$(TOOL) \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--seed 34 \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--output-directory $@ \
		--keep-intermediates \
		--predict-for-training \
		--try-include-mapq \
		--subsampling-series 0.5,1.0 \
		--trials 3 \
		--write-orig-mapq \
		--write-precise-mapq \
		--assess-accuracy \
		--model-family GradientBoosting \
		--verbose \
		--threads $(NTHREADS) \
		-- -p $(NTHREADS) --reorder

.r_1.fastq: .rtmp_1.fastq
	$(MASON_CONVERT) --in1 $< --in2 $(<:_1.fastq=_2.fastq) --out1 $@ --out2 $(@:_1.fastq=_2.fastq)

.rtmp_1.fastq: lambda_virus.fa $(MASON)
	$(MASON) \
		illumina \
		-i \
		-s 1 \
		-sq \
		-mp \
		-rn 2 \
		-ll 300 \
		-le 40 \
		-n 8 \
		-N 250000 \
		-o .rtmp.fastq \
		lambda_virus.fa

$(BT2):
	$(MAKE) -C $(BT2_DIR)

$(MASON):
	$(MAKE) -C ../../software/mason

lambda_virus.fa: $(BT2)
	cp `find $(BT2_DIR) -name $@` .

lambda_virus.1.bt2:
	cp `find $(BT2_DIR) -name lambda_virus.*.bt2` .

.PHONY: clean
clean:
	rm -rf full_e2e extra_e2e test_random_forest test_gradient_boost test_reweight_mapq test_reweight_ratio test_max_leaf_nodes
	rm -f .r*.fastq
	rm -f lambda_virus.*
