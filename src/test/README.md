To generate unpaired lambda_virus reads:

$TS_HOME/software/mason/mason \
    illumina \
    -i \
    -s 1 \
    -sq \
    -n 10 \
    -N 50000 \
    -o $HOME/.lambda_mason_unp.fastq \
    $BT2_HOME/example/reference/lambda_virus.fa

python $TS_HOME/bin/mason_convert.py \
    --in1 $HOME/.lambda_mason_unp.fastq \
    --out1 $HOME/lambda_mason_unp.fastq

To generate paired-end lambda_virus reads:

$TS_HOME/software/mason/mason \
    illumina \
    -i \
    -s 1 \
    -sq \
    -mp \
    -rn 2 \
    -ll 300 \
    -le 40 \
    -n 10 \
    -N 25000 \
    -o $HOME/.lambda_mason_paired.fastq \
    $BT2_HOME/example/reference/lambda_virus.fa

python $TS_HOME/bin/mason_convert.py \
    --in1 $HOME/.lambda_mason_paired_1.fastq \
    --in2 $HOME/.lambda_mason_paired_2.fastq \
    --out1 $HOME/lambda_mason_paired_1.fastq \
    --out2 $HOME/lambda_mason_paired_2.fastq
