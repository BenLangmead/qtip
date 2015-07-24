#!/bin/sh

OUTPUT_DIR=$HOME/git/ts-paper/preprint/figures
mkdir -p $OUTPUT_DIR

#
# Various lengths
#

python paper_plot.py ill_various_length_r0_bt2s_mason_ill_50_100k.pred_ \
                     ill_various_length_r0_bt2s_mason_ill_100_100k.pred_ \
                     ill_various_length_r0_bt2s_mason_ill_150_100k.pred_ \
                     ill_various_length_r0_bt2s_mason_ill_250_100k.pred_ \
                     ill_various_length_r12_bt2s50_mason_ill_50_100k.pred_ \
                     ill_various_length_r12_bt2s100_mason_ill_100_100k.pred_ \
                     ill_various_length_r12_bt2s150_mason_ill_150_100k.pred_ \
                     ill_various_length_r12_bt2s250_mason_ill_250_100k.pred_

cp tmp.eps $OUTPUT_DIR/various_length.eps

#
# Various lengths (local)
#

python paper_plot.py ill_various_length_r0_bt2sl_mason_ill_50_100k.pred_ \
                     ill_various_length_r0_bt2sl_mason_ill_100_100k.pred_ \
                     ill_various_length_r0_bt2sl_mason_ill_150_100k.pred_ \
                     ill_various_length_r0_bt2sl_mason_ill_250_100k.pred_ \
                     ill_various_length_r12_bt2sl50_mason_ill_50_100k.pred_ \
                     ill_various_length_r12_bt2sl100_mason_ill_100_100k.pred_ \
                     ill_various_length_r12_bt2sl150_mason_ill_150_100k.pred_ \
                     ill_various_length_r12_bt2sl250_mason_ill_250_100k.pred_

cp tmp.eps $OUTPUT_DIR/various_length_local.eps

#
# Various sensitivities
#

for SIZE in 100 250 ; do
for LOCAL in '' 'l' ; do

python paper_plot.py various_sensitivities_r0_bt2vf${LOCAL}_mason_ill_${SIZE}_100k.pred_ \
                     various_sensitivities_r0_bt2f${LOCAL}_mason_ill_${SIZE}_100k.pred_ \
                     various_sensitivities_r0_bt2s${LOCAL}_mason_ill_${SIZE}_100k.pred_ \
                     various_sensitivities_r0_bt2vs${LOCAL}_mason_ill_${SIZE}_100k.pred_ \
                     various_sensitivities_r12_bt2vf${LOCAL}${SIZE}_mason_ill_${SIZE}_100k.pred_ \
                     various_sensitivities_r12_bt2f${LOCAL}${SIZE}_mason_ill_${SIZE}_100k.pred_ \
                     various_sensitivities_r12_bt2s${LOCAL}${SIZE}_mason_ill_${SIZE}_100k.pred_ \
                     various_sensitivities_r12_bt2vs${LOCAL}${SIZE}_mason_ill_${SIZE}_100k.pred_

LOCAL_STR="e2e"
if [ "$LOCAL" = "l" ] ; then
    LOCAL_STR="local"
fi
cp tmp.eps $OUTPUT_DIR/various_sensitivities_${SIZE}_${LOCAL_STR}.eps

done
done

#
# Various genomes
#

for SIZE in 100 250 ; do
python paper_plot.py various_genomes_r0_bt2s_hg_mason_ill_hg_${SIZE}_100k.pred_ \
                     various_genomes_r0_bt2s_mm_mason_ill_mm_${SIZE}_100k.pred_ \
                     various_genomes_r0_bt2s_zm_mason_ill_zm_${SIZE}_100k.pred_ \
                     various_genomes_r12_bt2s${SIZE}_hg_mason_ill_hg_${SIZE}_100k.pred_ \
                     various_genomes_r12_bt2s${SIZE}_mm_mason_ill_mm_${SIZE}_100k.pred_ \
                     various_genomes_r12_bt2s${SIZE}_zm_mason_ill_zm_${SIZE}_100k.pred_

cp tmp.eps $OUTPUT_DIR/various_genomes_${SIZE}.eps
done

#
# Various simulators
#

for SIZE in 100 250 ; do
python paper_plot.py various_simulators_r0_bt2s_art_ill_hg_${SIZE}_100k.pred_ \
                     various_simulators_r0_bt2s_mason_ill_hg_${SIZE}_100k.pred_ \
                     various_simulators_r0_bt2s_wgsim_ill_hg_${SIZE}_100k.pred_ \
                     various_simulators_r12_bt2s${SIZE}_art_ill_hg_${SIZE}_100k.pred_ \
                     various_simulators_r12_bt2s${SIZE}_mason_ill_hg_${SIZE}_100k.pred_ \
                     various_simulators_r12_bt2s${SIZE}_wgsim_ill_hg_${SIZE}_100k.pred_

cp tmp.eps $OUTPUT_DIR/various_simulators_${SIZE}.eps
done

#
# Various aligners
#

python paper_plot.py various_aligners_r0_bt2sl_mason_ill_100_100k.pred_ \
                     various_aligners_r0_bt2sl_mason_ill_250_100k.pred_ \
                     various_aligners_r0_bwamem_mason_ill_100_100k.pred_ \
                     various_aligners_r0_bwamem_mason_ill_250_100k.pred_ \
                     various_aligners_r12_bt2sl100_mason_ill_100_100k.pred_ \
                     various_aligners_r12_bt2sl250_mason_ill_250_100k.pred_ \
                     various_aligners_r12_bwamem100_mason_ill_100_100k.pred_ \
                     various_aligners_r12_bwamem250_mason_ill_250_100k.pred_

cp tmp.eps $OUTPUT_DIR/various_aligners.eps
