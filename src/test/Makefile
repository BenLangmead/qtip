MASON_DIR=../../software/mason
MASON=$(MASON_DIR)/mason
BT2_DIR=../../software/bowtie2
BT2=$(BT2_DIR)/bowtie2
MASON_CONVERT=python ../../bin/mason_convert.py

all: full_e2e/final.sam extra_e2e/final.sam final_full_e2e.sam test_random_forest test_gradient_boost

final_full_e2e.sam: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../qsim \
		--ref lambda_virus.fa \
		--m1 .r_1.fastq \
		--m2 .r_2.fastq \
		--index lambda_virus \
		--aligner bowtie2 \
		--bt2-exe $(BT2) \
		--vanilla-output $@ \
		--write-orig-mapq \
		--write-precise-mapq

full_e2e/final.sam: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../qsim \
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
		--assess-accuracy

extra_e2e/final.sam: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../qsim \
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
		--try-exclude-mapq \
		--subsampling-series 0.1,0.5,1.0 \
		--trials 2 \
		--write-orig-mapq \
		--write-precise-mapq \
		--assess-accuracy \
		--no-collapse

#
# Try a few different machine learning models
#
test_random_forest: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../qsim \
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
		--try-exclude-mapq \
		--subsampling-series 0.1,0.5,1.0 \
		--trials 2 \
		--write-orig-mapq \
		--write-precise-mapq \
		--assess-accuracy \
		--model-family RandomForest \
		--model-params 50:8:0.5

test_gradient_boost: .r_1.fastq .r_2.fastq lambda_virus.fa lambda_virus.1.bt2 $(BT2)
	python ../qsim \
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
		--try-exclude-mapq \
		--subsampling-series 0.5,1.0 \
		--trials 3 \
		--write-orig-mapq \
		--write-precise-mapq \
		--assess-accuracy \
		--model-family GradientBoosting \
		--verbose

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
		-N 25000 \
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
	rm -rf full_e2e extra_e2e test_random_forest test_gradient_boost
	rm -f .r*.fastq
	rm -f lambda_virus.*
