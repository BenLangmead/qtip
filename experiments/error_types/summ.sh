#!/bin/sh

#
# Summarize the errors in the output files from align.sh
#

python summ_errs.py --in sim.sensitive.sam > results.sensitive.txt
python summ_errs.py --in sim.very-sensitive.sam > results.very-sensitive.txt
