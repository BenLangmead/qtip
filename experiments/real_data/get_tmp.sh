#!/bin/sh

for d in https://s3.amazonaws.com/ss-reads/ERR194147_1M/ERR194147_1.first100M.fastq.000000 ; do
    bs=`basename $d`
    wget -O `echo $bs | sed 's/\.[0-9]*$//'` $d
done

# I need a nice clean command to go from FASTQ to SAM and recalibrated SAM