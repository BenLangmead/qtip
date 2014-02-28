"""
predict.py
"""

import os
import sys
import random
import logging
import numpy as np
import matplotlib.pyplot as plt
import math
from collections import defaultdict
from pandas import DataFrame
from pandas.io.parsers import read_csv
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from scipy.interpolate import splrep, splev
from matplotlib.backends.backend_pdf import PdfPages


def read_dataset(prefix):
    dfs = {}
    for name, suf, paired in [('u', '_unp.csv', False), ('m', '_mates.csv', False), ('c', '_conc.csv', True)]:
        fn = prefix + suf
        if os.path.exists(fn):
            dfs[name] = read_csv(fn, quoting=2)
        elif os.path.exists(fn + '.gz'):
            dfs[name] = read_csv(fn + '.gz', quoting=2, compression='gzip')
        elif os.path.exists(fn + '.bz2'):
            dfs[name] = read_csv(fn + '.bz2', quoting=2, compression='bz2')
        else:
            raise RuntimeError('No such file: "%s"' % fn)

    for df in dfs.itervalues():
        if df['correct'].count() == len(df['correct']):
            df['correct'] = df['correct'].map(lambda x: 1 if x == 'T' else 0)

    for df in [dfs['u'], dfs['m']]:
        diffv = df.maxv - df.minv
        df['bestnorm'] = (1.0 * df.best - df.minv) / diffv
        secbestnorm = (1.0 * df.secbest - df.minv) / diffv
        df['secbestnorm'] = secbestnorm.fillna(secbestnorm.min())
        df['diffnorm'] = df.bestnorm - df.secbestnorm

    conc = dfs['c']

    diffv_1 = conc.maxv1 - conc.minv1
    conc['bestnorm1'] = (1.0 * conc.best1 - conc.minv1) / diffv_1
    secbestnorm_1 = (1.0 * conc.secbest1 - conc.minv1) / diffv_1
    conc['secbestnorm1'] = secbestnorm_1.fillna(secbestnorm_1.min())
    conc['diffnorm1'] = conc['bestnorm1'] - conc['secbestnorm1']

    diffv_2 = conc.maxv2 - conc.minv2
    conc['bestnorm2'] = (1.0 * conc.best2 - conc.minv2) / diffv_2
    secbestnorm_2 = (1.0 * conc.secbest2 - conc.minv2) / diffv_2
    conc['secbestnorm2'] = secbestnorm_2.fillna(secbestnorm_2.min())
    conc['diffnorm2'] = conc['bestnorm2'] - conc['secbestnorm2']

    diff_conc = diffv_1 + diffv_2
    conc['bestnormconc'] = (1.0 * conc.best1 + conc.best2 - conc.minv1 - conc.minv2) / diff_conc
    secbest_conc = (1.0 * conc.secbest1 + conc.secbest2 - conc.minv1 - conc.minv2) / diff_conc
    conc['secbestnormconc'] = secbest_conc.fillna(secbest_conc.min())
    conc['diffnormconc'] = conc['bestnormconc'] - conc['secbestnormconc']

    return dfs


def pcor_to_mapq(p):
    return map(lambda x: -10.0 * math.log10(1.0 - x) if x < 1.0 else float('inf'), p)


def mapq_to_pcor(p):
    return map(lambda x: 1.0 - (10.0 ** (-0.1 * x)) if x < float('inf') else 1.0, p)


def round_pcor(pcor):
    assert not any(map(lambda x: x < 0.0 or x > 1.0, pcor))
    return mapq_to_pcor(map(round, (pcor_to_mapq(pcor))))


def _tally_cor_per(level, cor):
    tally = defaultdict(lambda: [0, 0])
    for p, c in zip(level, cor):
        c = 0 if c else 1
        tally[p][c] += 1
    return tally


def ranking_error(pcor, cor, rounded=False):
    assert len(pcor) == len(cor)
    if rounded:
        pcor = round_pcor(pcor)
    tally = _tally_cor_per(pcor, cor)
    err, sofar = 0, 0
    # from low-confidence to high confidence
    for p in sorted(tally.iterkeys()):
        ncor, nincor = tally[p]
        ntot = ncor + nincor
        assert ntot > 0
        if nincor > 0:
            # spread error over this grouping of tied pcors
            frac = float(nincor) / ntot
            assert frac <= 1.0
            err += frac * sum(xrange(sofar, sofar + ntot))
        sofar += ntot
    return err


def ranking_error_by_quantile(pcor, cor, nquantiles=10, rounded=False):
    assert len(pcor) == len(cor)
    if rounded:
        pcor = round_pcor(pcor)
    tally = _tally_cor_per(pcor, cor)
    errs, sofar = [0] * nquantiles, 0
    # from low-confidence to high confidence
    cur_quantile = 0
    for p in sorted(tally.iterkeys()):
        ncor, nincor = tally[p]
        ntot = ncor + nincor
        assert ntot > 0
        if nincor > 0:
            # spread error over this grouping of tied pcors
            frac = float(nincor) / ntot
            assert frac <= 1.0
            errs[cur_quantile] += frac * sum(xrange(sofar, sofar + ntot))
        sofar += ntot
        cur_quantile = (sofar * nquantiles) // len(pcor)
    return errs


def smoothfit_error(pcor, cor, rounded=False, s=125.0, plot=False, log2ize=False):
    """ Calculate the smooth-fit error.  Possibly also plot the smooth
        fit. """
    assert len(pcor) == len(cor)
    if rounded:
        pcor = mapq_to_pcor(map(round, (pcor_to_mapq(pcor))))
    tally = _tally_cor_per(pcor, cor)
    # create parallel pcor and cor lists
    sorted_pcor, sorted_cor = [], []
    for p in sorted(tally.iterkeys(), reverse=True):
        ncor, nincor = tally[p]
        ntot = ncor + nincor
        sorted_pcor.extend([p] * ntot)
        # use fractional correctness so results are deterministic
        sorted_cor.extend([float(ncor)/ntot] * ntot)
    m = len(sorted_cor)
    assert m == len(sorted_pcor) == len(pcor)
    y = splev(xrange(m), splrep(xrange(m), sorted_cor, s=s))
    if plot:
        ln = len(pcor)
        x = np.linspace(1.0/ln, (ln - 1.0) / ln, ln - 1.0)
        if log2ize:
            x = -np.log2(x[::-1])
        maxx = math.ceil(np.log2(ln))
        fig = plt.figure(figsize=(12, 4))
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        axes.plot(x, sorted_pcor[:-1], color='r', label='Pcor')
        axes.plot(x, y[:-1], color='k', label='Smoothed cor')
        axes.set_xlabel('Index')
        if log2ize:
            axes.set_xticklabels(map(lambda z: 2 ** z, np.linspace(-1, -maxx, maxx)), rotation=90)
        axes.set_ylabel('Probability')
        axes.set_title('Smooth-fit plot')
        axes.legend(loc=3)
        axes.grid(True)
    y_null = [np.mean(pcor)] * m
    # smaller is better
    return mean_squared_error(sorted_pcor, y) / mean_squared_error(sorted_pcor, y_null), y


def bucket_error(mapq, cor):
    """ Rounded given mapqs off to nearest integer.  For each bin, if
        there is at least one incorrect alignment in the bin, compare
        predicted MAPQ with "actual". """
    tally = _tally_cor_per(mapq, cor)
    err, n = 0, 0
    for m, t in tally.iteritems():
        ncor, nincor = t
        ntot = ncor + nincor
        assert ntot > 0
        if nincor > 0 and ncor > 0:
            pincor = float(nincor) / ntot
            actual = -10.0 * math.log10(pincor)
            err += ((actual - m) ** 2)
            n += 1
    return float(err)/n  # mean squared error


def _drop_rate_cum_sum(pcor, cor):
    """ Create cumulative sum of # of incorrect alignments starting at
        high pcor and moving to low pcor.  For blocks of equal pcors,
        "charge" a fractional amount for the incorrect alignments in
        that stratum. """
    tally = _tally_cor_per(pcor, cor)
    cumsum, cumsums = 0, []
    for p, tup in sorted(tally.iteritems(), reverse=True):
        ncor, nincor = tup
        for _ in xrange(ncor + nincor):
            cumsum += (float(nincor) / (ncor + nincor))
            cumsums.append(cumsum)
    return cumsums


def _plot_drop_rate(pcor, cor, pcor2=None, log2ize=False):
    """ Create a plot showing drop rate on x axis versus cumulative
        number of incorrect alignments on y axis.  Intended as an
        easier-to-read-all-at-once alternative to a ROC-like plot. """
    cumsum = _drop_rate_cum_sum(pcor, cor)
    cumsum2 = None if pcor2 is None else _drop_rate_cum_sum(pcor2, cor)
    ln = len(pcor)
    x_log = np.linspace(1.0/ln, (ln - 1.0) / ln, ln - 1.0)
    if log2ize:
        x_log = -np.log2(x_log[::-1])
    maxx = math.ceil(np.log2(ln))
    fig = plt.figure(figsize=(10, 4))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.plot(x_log, cumsum[:-1], color='r', label='Model')
    if cumsum2 is not None:
        axes.plot(x_log, cumsum2[:-1], color='k', label='Original')
    axes.set_xlabel('MAPQ drop rate')
    if log2ize:
        axes.set_xticklabels(map(lambda x: 2 ** x, np.linspace(-1, -maxx, maxx)), rotation=90)
    axes.set_ylabel('Cumulative # incorrect alignments')
    axes.set_title('Drop rate comparison on -log2 scale')
    axes.legend(loc=2)
    axes.grid(True)
    return fig


def summarize_performance(pcor_test, mapq_test, pcor_orig_test, mapq_orig_test, y_test):

    # ranking error
    rank_err = ranking_error(pcor_test, y_test)
    rank_err_rounded = ranking_error(pcor_test, y_test, rounded=True)
    rank_err_orig = ranking_error(mapq_orig_test, y_test)
    rank_err_improve = (100.0 * (rank_err_orig - rank_err_rounded)) / max(rank_err_orig, 1)

    # smooth-fit error
    sfit_err = smoothfit_error(pcor_test, y_test)
    sfit_err_rounded = smoothfit_error(pcor_test, y_test, rounded=True)
    sfit_err_orig = smoothfit_error(pcor_orig_test, y_test)
    sfit_err_improve = ((100.0 * (sfit_err_orig - sfit_err_rounded))/max(sfit_err_orig, 1))

    # bucket error
    bucket_err = bucket_error(mapq_test, y_test)
    bucket_err_orig = bucket_error(mapq_orig_test, y_test)
    bucket_err_improve = ((100.0 * (bucket_err_orig - bucket_err))/max(bucket_err_orig, 1))

    logging.info('  Original errors: %.2e/%.3f/%.2e' % (rank_err_orig, sfit_err_orig, bucket_err_orig))
    logging.info('  Errors with model: %.2e/%.3f/%.2e' % (rank_err, sfit_err, bucket_err))
    logging.info('  Errors with model (rounded): %.2e/%.3f' % (rank_err_rounded, sfit_err_rounded))
    logging.info('  Improvement: %0.2f%%/%0.2f%%/%0.2f%%' % (rank_err_improve, sfit_err_improve, bucket_err_improve))

    return rank_err_improve, sfit_err_improve, bucket_err_improve


def go(args):
    # Load training and test
    logging.info('Loading training data')
    train_dfs = read_dataset(args.training_prefix)
    logging.info('Loading test data')
    test_dfs = read_dataset(args.test_prefix)

    logging.info('Fitting models')
    datasets = [('Unpaired', 'u', train_dfs['u'], test_dfs['u'], False),
                ('Mate', 'm', train_dfs['m'], test_dfs['m'], False),
                ('Concordant', 'c', train_dfs['c'], test_dfs['c'], True)]
    models = [('RFR', lambda: RandomForestRegressor(n_estimators=args.random_forest_num_trees,
                                                    max_depth=args.random_forest_max_depth,
                                                    random_state=random.randint(0, sys.maxint),
                                                    max_features='auto'))]
    for dataset_name, dataset_shortname, train, test, paired in datasets:
        if test.shape[0] == 0:
            logging.info('  No test data; skipping...')
            continue
        if train.shape[0] == 0:
            logging.error('  Test data exists, but there is no training data...')
            raise RuntimeError('Test data exists, but there is no training data')
        if paired:
            x_train = train[['bestnorm1', 'diffnorm1', 'diffnormconc']].values
            x_test = test[['bestnorm1', 'diffnorm1', 'diffnormconc']].values
            mapq_orig_train = train['mapq1']
            mapq_orig_test = test['mapq1']
        else:
            x_train = train[['bestnorm', 'diffnorm']].values
            x_test = test[['bestnorm', 'diffnorm']].values
            mapq_orig_train = train['mapq']
            mapq_orig_test = test['mapq']
        pcor_orig_test = mapq_to_pcor(mapq_orig_test)
        y_train = map(lambda x: x == 1, train.correct)
        y_test = map(lambda x: x == 1, test.correct)
        test_data_has_correct = test['correct'].count() == len(test['correct'])
        for model_name, model_gen in models:
            nsamples = x_train.shape[0]
            logging.info('  Fitting "%s" on test dataset "%s" (%d samples)' % (model_name, dataset_name, nsamples))
            model = model_gen()
            model.fit(x_train, y_train)
            pcor_test = model.predict(x_test)
            mapq_test = pcor_to_mapq(pcor_test)
            rank_err_improves, sfit_err_improves, bucket_err_improves = [], [], []
            rank_err_quantile, rank_err_quantile_orig = None, None
            if test_data_has_correct:
                rank_err_improve, sfit_err_improve, bucket_err_improve = \
                    summarize_performance(pcor_test, mapq_test, pcor_orig_test, mapq_orig_test, y_test)
                rank_err_improves.append(rank_err_improve)
                sfit_err_improves.append(sfit_err_improve)
                bucket_err_improves.append(bucket_err_improve)
                rank_err_quantile = ranking_error_by_quantile(pcor_test, y_test, rounded=True)
                rank_err_quantile_orig = ranking_error_by_quantile(pcor_orig_test, y_test)
            logging.info('  Writing results')
            result_test_df = DataFrame.from_items([('pcor', pcor_test), ('mapq', mapq_test),
                                                   ('orig', mapq_orig_test), ('correct', y_test)])
            result_test_df.to_csv(args.test_prefix + '_' + dataset_shortname + '_' + model_name + '.csv', index=False)
            if args.downsample:
                logging.info('  Trying various downsampling fractions')
                fractions = [0.5, 0.2, 0.15, 0.1, 0.08, 0.05, 0.035, 0.025, 0.02, 0.01]
                for fraction in fractions:
                    nsamples = int(fraction * x_train.shape[0])
                    logging.info('    Trying %0.03f%% (%d samples)' % (fraction * 100, nsamples))
                    samp = random.sample(xrange(x_train.shape[0]), nsamples)
                    x_train_samp = x_train[samp, ]
                    y_train_samp = [y_train[i] for i in samp]
                    model = model_gen()
                    model.fit(x_train_samp, y_train_samp)
                    pcor_test = model.predict(x_test)
                    mapq_test = pcor_to_mapq(pcor_test)
                    if test_data_has_correct:
                        rank_err_improve, sfit_err_improve, bucket_err_improve = \
                            summarize_performance(pcor_test, mapq_test, pcor_orig_test, mapq_orig_test, y_test)
                        rank_err_improves.append(rank_err_improve)
                        sfit_err_improves.append(sfit_err_improve)
                        bucket_err_improves.append(bucket_err_improve)
                    logging.info('      Writing downsampled results')
                    result_test_df = DataFrame.from_items([('pcor', pcor_test), ('mapq', mapq_test),
                                                           ('orig', mapq_orig_test), ('correct', y_test)])
                    result_test_df.to_csv('%s_%s_%s_down%02d.csv' % (args.test_prefix, dataset_shortname, model_name,
                                                                     int(fraction*100)), index=False)
                # write entire
                if test_data_has_correct:
                    downsample_series_df = DataFrame.from_items([
                        ('fraction', [1.0] + fractions),
                        ('rank_err_improvement', rank_err_improves),
                        ('smoothfit_err_improvement', sfit_err_improves),
                        ('bucket_err_improvement', bucket_err_improves)])
                    downsample_series_df.to_csv('%s_%s_%s_down_series.csv' % (args.test_prefix, dataset_shortname,
                                                                              model_name), index=False)

            # Difference in ranking error between model and original, broken
            # down by decile.  Helpful for distinguishing whether differences
            # are in the lower or higher quantiles.
            if test_data_has_correct:
                quantiles = np.linspace(0.0, 0.9, 10)
                ranking_error_quantile_df = DataFrame.from_items([
                    ('quantiles', quantiles),
                    ('model', rank_err_quantile),
                    ('original', rank_err_quantile_orig)])
                ranking_error_quantile_df['model_minus_original'] =\
                    ranking_error_quantile_df['model'] - ranking_error_quantile_df['original']
                ranking_error_quantile_df.to_csv('%s_%s_%s_rank_err_q.csv' % (args.test_prefix, dataset_shortname,
                                                                              model_name), index=False)

            if args.training_results:
                logging.info('  Fitting "%s" on training dataset "%s"' % (model_name, dataset_name))
                pcor_train = model.predict(x_train)
                mapq_train = pcor_to_mapq(pcor_train)
                logging.info('  Writing results')
                result_train_df = DataFrame.from_items([('pcor', pcor_train), ('mapq', mapq_train),
                                                        ('orig', mapq_orig_train), ('correct', y_train)])
                result_train_df.to_csv(args.training_prefix + '_' + dataset_shortname + '_' + model_name + '.csv',
                                       index=False)

    logging.info('Done')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Fit models, etc.')

    parser.add_argument('--training-prefix', metavar='path', type=str, default='training',
                        help='Prefix for files with training data')
    parser.add_argument('--test-prefix', metavar='path', type=str, default='test',
                        help='Prefix for files with test data')
    parser.add_argument('--training-results', action='store_const', const=True, default=False,
                        help='Use model to fit on training data and emit results')
    parser.add_argument('--downsample', action='store_const', const=True, default=False,
                        help='Do a series of downsamples to see how well model does with less data')
    parser.add_argument('--random-forest-num-trees', metavar='int', type=int, default=35, required=False,
                        help='Number of trees to use in random forest regression')
    parser.add_argument('--random-forest-max-depth', metavar='int', type=int, default=6, required=False,
                        help='Maximum depth of trees in random forest')
    parser.add_argument('--seed', metavar='int', type=int, default=99099, required=False,
                        help='Integer to initialize pseudo-random generator')
    parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG if args.verbose else logging.INFO)

    random.seed(args.seed)
    if args.profile:
        import cProfile
        cProfile.run('go(args)')
    else:
        go(args)
