import os
import sys
import pandas as pd
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import logging

qsim_dir = os.environ['HOME'] + '/git/mapq'
simdata_dir = os.path.join(qsim_dir, 'data', 'simulated_reads')
roc_table_base = os.path.join(simdata_dir, 'summary', 'roc_table', 'test')
sub_table_base = os.path.join(simdata_dir, 'summary', 'subsampling_tables', 'test')


def read_roc(fn):
    tab = pd.read_csv(fn, sep='\t')
    tab['cum'] = tab['cum_cor'] + tab['cum_incor']
    tab['pcor'] = 1.0 - 10.0 ** (-0.1 * tab['mapq'])
    tab['sse'] = tab['incor'] * tab['pcor'] * tab['pcor'] + \
                 tab['cor'] * (1.0 - tab['pcor']) * (1.0 - tab['pcor'])
    tab['cum_sse'] = tab['sse'].cumsum()
    tab['n'] = tab['cor'] + tab['incor']
    return tab


def cum_diff(c1, ci1, m1, c2, ci2, m2):
    i, j = 0, 0
    cum_incor1, cum_incor2, cum, mapq1, mapq2 = [], [], [], [], []
    c1, c2 = [0] + c1, [0] + c2
    ci1, ci2 = [0] + ci1, [0] + ci2
    m1, m2 = [m1[0]] + m1, [m2[0]] + m2
    while i < len(c1) and j < len(c2):
        if c1[i] == c2[j]:
            cum_incor1.append(ci1[i])
            mapq1.append(m1[i])
            cum_incor2.append(ci2[j])
            mapq2.append(m2[j])
            cum.append(c1[i])
            i += 1
            j += 1
        elif c1[i] < c2[j]:
            assert i > 0 and j > 0
            cum_incor1.append(ci1[i])
            mapq1.append(m1[i])
            prev, nxt = ci2[j-1], ci2[j]
            along = c1[i] - c2[j-1]
            tot = c2[j] - c2[j-1]
            frac = float(nxt - prev) * along / tot
            cum_incor2.append(prev + frac)
            mapq2.append(m2[j])
            cum.append(c1[i])
            i += 1
        else:
            assert i > 0 and j > 0
            cum_incor2.append(ci2[j])
            mapq2.append(m2[j])
            prev, nxt = ci1[i-1], ci1[i]
            along = c2[j] - c1[i-1]
            tot = c1[i] - c1[i-1]
            frac = float(nxt - prev) * along / tot
            cum_incor1.append(prev + frac)
            mapq1.append(m1[i])
            cum.append(c2[j])
            j += 1
    assert i == len(c1) and j == len(c2)
    ret = pd.DataFrame.from_dict({'cum': cum,
                                  'cum_incor1': cum_incor1,
                                  'cum_incor2': cum_incor2,
                                  'mapq1': mapq1,
                                  'mapq2': mapq2})
    ret['diff'] = ret['cum_incor1'] - ret['cum_incor2']
    return ret  # negative means 1 is better


def cum_plots(prefix):
    roc_table_fn = os.path.join(roc_table_base, prefix + 'roc_table.tsv')
    roc_table_orig_fn = os.path.join(roc_table_base, prefix + 'roc_table_orig.tsv')
    roc = read_roc(roc_table_fn)
    roc_orig = read_roc(roc_table_orig_fn)

    # Get cum_diff for # incorrect
    sse_df = cum_diff(roc.cum.tolist(),
                      roc.cum_sse.tolist(),
                      roc.mapq.tolist(),
                      roc_orig.cum.tolist(),
                      roc_orig.cum_sse.tolist(),
                      roc_orig.mapq.tolist())
    incor_df = cum_diff(roc.cum.tolist(),
                        roc.cum_incor.tolist(),
                        roc.mapq.tolist(),
                        roc_orig.cum.tolist(),
                        roc_orig.cum_incor.tolist(),
                        roc_orig.mapq.tolist())

    sse_cum = sse_df['cum'].tolist()
    sse_incor1, sse_incor2 = sse_df['cum_incor1'].tolist(), sse_df['cum_incor2'].tolist()
    incor_cum = incor_df['cum'].tolist()
    incor_incor1, incor_incor2 = incor_df['cum_incor1'].tolist(), incor_df['cum_incor2'].tolist()
    mapq = sse_df['mapq1'].tolist()
    mapq_orig = sse_df['mapq2'].tolist()

    return sse_cum, sse_incor1, sse_incor2, incor_cum, incor_incor1, incor_incor2, mapq, mapq_orig


def go():
    prefixes = sys.argv[1:]

    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG)

    robjects.r('''source('paper_plot.R')''')
    plot_cums = robjects.globalenv['plot_cums']
    grd = importr('grDevices')

    columns = 2
    if True:
        grd.setEPS()
        grd.postscript(file='tmp.eps', width=columns * 5, height=3 * len(prefixes)/columns)
    else:
        grd.pdf(file='tmp.pdf', width=columns * 5, height=3 * len(prefixes)/columns)
    robjects.r('''layout(matrix(1:%d, ncol = %d), widths = 1, heights = rep(c(10,1,1,10), %d), respect = F)''' %
               (4 * len(prefixes), columns, len(prefixes)))
    for prefix in prefixes:
        logging.info('Prefix: %s' % prefix)
        sse_cum, sse_incor1, sse_incor2, incor_cum, incor_incor1, incor_incor2, mapq, mapq_orig = cum_plots(prefix)
        sse_cum = robjects.FloatVector(sse_cum)
        sse_incor1 = robjects.FloatVector(sse_incor1)
        sse_incor2 = robjects.FloatVector(sse_incor2)
        incor_cum = robjects.FloatVector(incor_cum)
        incor_incor1 = robjects.FloatVector(incor_incor1)
        incor_incor2 = robjects.FloatVector(incor_incor2)
        mapq = robjects.FloatVector(mapq)
        mapq_orig = robjects.FloatVector(mapq_orig)
        plot_cums(columns, sse_cum, sse_incor1, sse_incor2, incor_cum, incor_incor1, incor_incor2, mapq, mapq_orig)
    grd.dev_off()

if __name__ == '__main__':
    go()
