"""
Given a directory with output from ts.py, predict new MAPQs.
"""

__author__ = 'langmead'

import pandas
import os
import sys
import math
import random
import matplotlib.pyplot as plt
import numpy as np
import logging
import gc
import cPickle
from sklearn import cross_validation
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor, GradientBoostingRegressor
from collections import defaultdict, Counter


VERSION = '0.1.0'


def pcor_to_mapq(p):
    """ Convert probability correct (pcor) to mapping quality (MAPQ) """
    return abs(-10.0 * math.log10(1.0 - p)) if p < 1.0 else float('inf')


def pcor_to_mapq_iter(p):
    """ Apply pcor_to_mapq to every element of iterable """
    return map(pcor_to_mapq, p)


def mapq_to_pcor(p):
    """ Convert mapping quality (MAPQ) to probability correct (pcor) """
    return (1.0 - (10.0 ** (-0.1 * p))) if p < float('inf') else 1.0


def mapq_to_pcor_iter(p):
    """ Apply mapq_to_pcor to every element of iterable """
    return map(mapq_to_pcor, p)


def round_pcor_iter(pcor):
    """ Modify pcor so that pcors correspond to nearest integer-rounded MAPQ """
    assert not any(map(lambda x: x < 0.0 or x > 1.0, pcor))
    return mapq_to_pcor_iter(map(round, (pcor_to_mapq_iter(pcor))))


def read_dataset(prefix):
    """ Read in training and test data having with given path prefix. """
    dfs = {}
    for name, suf, paired in [('u', '_unp.csv', False), ('d', '_disc.csv', True),
                              ('c', '_conc.csv', True), ('b', '_bad_end.csv', False)]:
        fn = prefix + suf
        if os.path.exists(fn):
            dfs[name] = pandas.io.parsers.read_csv(fn, quoting=2)
        elif os.path.exists(fn + '.gz'):
            dfs[name] = pandas.io.parsers.read_csv(fn + '.gz', quoting=2, compression='gzip')
        elif os.path.exists(fn + '.bz2'):
            dfs[name] = pandas.io.parsers.read_csv(fn + '.bz2', quoting=2, compression='bz2')
        else:
            raise RuntimeError('No such file: "%s"' % fn)

    for df in dfs.itervalues():
        if df['correct'].count() == len(df['correct']):
            df['correct'] = df['correct'].map(lambda x: 1 if x == 'T' else 0)

    def _standardize(df, nm, best_nm, secbest_nm, mn, diff):
        df[nm] = (df[best_nm] - df[secbest_nm]) / diff
        df[nm] = df[nm].fillna(np.nanmax(df[nm])).fillna(0)
        df[best_nm] = (df[best_nm].astype(float) - mn) / diff
        df[best_nm] = df[best_nm].fillna(np.nanmax(df[best_nm])).fillna(0)
        assert not any([math.isnan(x) for x in df[nm]])
        assert not any([math.isnan(x) for x in df[best_nm]])

    for df in [dfs['u'], dfs['b']]:
        # TODO: set minv/maxv properly if not defined
        if df.shape[0] == 0:
            continue
        df['minv'] = df['minv'].fillna(df['best'].min())
        df['maxv'] = df['maxv'].fillna(df['best'].max())
        _standardize(df, 'diff', 'best', 'secbest', df.minv, df.maxv - df.minv)

    for df in [dfs['d']]:
        if df.shape[0] == 0:
            continue
        df['minv1'] = df['minv1'].fillna(df['best1'].min())
        df['maxv1'] = df['maxv1'].fillna(df['best1'].max())
        _standardize(df, 'diff1', 'best1', 'secbest1', df.minv1, df.maxv1 - df.minv1)

    for df in [dfs['c']]:
        # TODO: set minv1/maxv1/minv2/maxv2 properly if not defined
        if df.shape[0] == 0:
            continue
        df['minv1'] = df['minv1'].fillna(df['best1'].min())
        df['maxv1'] = df['maxv1'].fillna(df['best1'].max())
        df['minv2'] = df['minv2'].fillna(df['best2'].min())
        df['maxv2'] = df['maxv2'].fillna(df['best2'].max())
        _standardize(df, 'diff1', 'best1', 'secbest1', df.minv1, df.maxv1 - df.minv1)
        _standardize(df, 'diff2', 'best2', 'secbest2', df.minv2, df.maxv2 - df.minv2)
        minconc, maxconc = df.minv1 + df.minv2, df.maxv1 + df.maxv2
        _standardize(df, 'diff_conc', 'bestconc', 'secbestconc', minconc, maxconc - minconc)
        df['best_min12'] = df[['best1', 'best2']].min(axis=1)
        df['diff_min12'] = df[['diff1', 'diff2']].min(axis=1)
        df['diff_min12conc'] = df[['diff1', 'diff2', 'diff_conc']].min(axis=1)
        df['fraglen_z'] = ((df['fraglen'] - df['fraglen'].mean()) / df['fraglen'].std()).abs()

    return dfs


def tally_cor_per(level, cor):
    """ Some predictions from our model are the same; this helper function
        gathers all the correct/incorrect information for each group of equal
        predictions. """
    tally = defaultdict(lambda: [0, 0])
    for p, c in zip(level, cor):
        c = 0 if c else 1
        tally[p][c] += 1
    return tally


def ranking_error(pcor, cor, rounded=False):
    """ Return the ranking error given a list of pcors and a parallel list of
        correct/incorrect booleans.  Round off to nearest MAPQ first if
        rounded=True.  """
    assert len(pcor) == len(cor)
    if rounded:
        pcor = round_pcor_iter(pcor)
    tally = tally_cor_per(pcor, cor)
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


def mse_error(pcor, cor, rounded=False):
    """ Return the mean squared error between the pcors (in [0, 1]) and and
        cors (in {0, 1}).  This is a measure of how close our predictions are
        to true correct/incorrect. """
    if rounded:
        pcor = round_pcor_iter(pcor)
    return np.mean((pcor - cor) ** 2) / np.mean((np.mean(cor) - cor) ** 2)


def cum_squared_error(pcor, cor, rounded=False):
    """ Return the cumulative squared error between pcors and cors, sorted
        from highest to lowest pcor. """
    if rounded:
        pcor = round_pcor_iter(pcor)
    pcor_order = sorted(range(len(pcor)), key=lambda x: pcor[x], reverse=True)
    pcor = np.array([pcor[x] for x in pcor_order])
    cor = np.array([cor[x] for x in pcor_order])
    assert all(pcor[i] >= pcor[i+1] for i in xrange(len(pcor)-1))
    return np.cumsum((pcor - cor) ** 2)


def plot_drop_rate_v_squared_error(pcor, cor, pcor2=None, log2ize=False, rasterize=False):
    cumsum = cum_squared_error(pcor, cor)
    cumsum2 = None if pcor2 is None else cum_squared_error(pcor2, cor)
    assert len(cumsum) == len(cumsum2)
    ln = len(pcor)
    x_log = np.linspace(1.0/ln, (ln - 1.0) / ln, ln)
    assert len(x_log) == len(cumsum)
    if log2ize:
        x_log = -np.log2(x_log[::-1])
    maxx = math.ceil(np.log2(ln))
    fig = plt.figure(figsize=(10, 4))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.plot(x_log, cumsum, color='r', label='Model', rasterized=rasterize)
    if cumsum2 is not None:
        axes.plot(x_log, cumsum2, color='k', label='Original', rasterized=rasterize)
    axes.set_xlabel('MAPQ drop rate')
    if log2ize:
        axes.set_xticklabels(map(lambda x: 2 ** x, np.linspace(-1, -maxx, maxx)), rotation=90)
    axes.set_ylabel('Cumulative squared error')
    axes.set_title('Cumulative MSE comparison %s' % ('on -log2 scale' if log2ize else ''))
    axes.legend(loc=2)
    axes.grid(True)
    return fig


def plot_drop_rate_v_squared_error_difference(pcor, pcor2, cor, log2ize=False):
    cumsum1 = cum_squared_error(pcor, cor)
    cumsum2 = cum_squared_error(pcor2, cor)
    ln = len(pcor)
    x_log = np.linspace(1.0/ln, (ln - 1.0) / ln, ln - 1.0)
    if log2ize:
        x_log = -np.log2(x_log[::-1])
    maxx = math.ceil(np.log2(ln))
    fig = plt.figure(figsize=(10, 4))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    cumsum_diff = np.subtract(cumsum1, cumsum2)
    xs, ys = x_log, cumsum_diff[:-1]
    axes.plot(xs, ys, color='k', label='Model')
    axes.set_xlabel('MAPQ drop rate')
    if log2ize:
        axes.set_xticklabels(map(lambda x: 2 ** x, np.linspace(-1, -maxx, maxx)), rotation=90)
    axes.set_ylabel('Difference in cumulative MSE')
    axes.set_title('Drop rate comparison on -log2 scale')
    axes.fill_between(xs, ys, where=ys>=0, interpolate=True, color='red')
    axes.fill_between(xs, ys, where=ys<=0, interpolate=True, color='green')
    axes.grid(True)
    return fig


def drop_rate_cum_sum(pcor, cor):
    """ Generate a vector giving (something like) the cumulative sum of
        incorrect alignments up to each p-value, from high to low p-value.
        When many p-values are tied, we distribute the incorrect fraction over
        all the elements in that range so that the answer doesn't depend on
        how equal p-values are ordered. """
    tally = tally_cor_per(pcor, cor)
    cumsum, cumsums = 0, []
    for p, tup in sorted(tally.iteritems(), reverse=True):
        ncor, nincor = tup
        for i in xrange(ncor + nincor):
            cumsum += (float(nincor) / (ncor + nincor))
            cumsums.append(cumsum)
    return cumsums


def plot_drop_rate(pcor, cor, pcor2=None, log2ize=False, rasterize=False):
    """  """
    cumsum = drop_rate_cum_sum(pcor, cor)
    cumsum2 = None if pcor2 is None else drop_rate_cum_sum(pcor2, cor)
    ln = len(pcor)
    x_log = np.linspace(1.0/ln, (ln - 1.0) / ln, ln - 1.0)
    if log2ize:
        x_log = -np.log2(x_log[::-1])
    maxx = math.ceil(np.log2(ln))
    fig = plt.figure(figsize=(10, 4))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.plot(x_log, cumsum[:-1], color='r', label='Model', rasterized=rasterize)
    if cumsum2 is not None:
        axes.plot(x_log, cumsum2[:-1], color='k', label='Original', rasterized=rasterize)
    axes.set_xlabel('MAPQ drop rate')
    if log2ize:
        axes.set_xticklabels(map(lambda x: 2 ** x, np.linspace(-1, -maxx, maxx)), rotation=90)
    axes.set_ylabel('Cumulative # incorrect alignments')
    axes.set_title('Drop rate comparison on -log2 scale')
    axes.legend(loc=2)
    axes.grid(True)
    return fig


def plot_drop_rate_difference(pcor, pcor2, cor, log2ize=False, rasterize=False):
    cumsum1, cumsum2 = drop_rate_cum_sum(pcor, cor), drop_rate_cum_sum(pcor2, cor)
    ln = len(pcor)
    x_log = np.linspace(1.0/ln, (ln - 1.0) / ln, ln - 1.0)
    if log2ize:
        x_log = -np.log2(x_log[::-1])
    maxx = math.ceil(np.log2(ln))
    fig = plt.figure(figsize=(10, 4))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    cumsum_diff = np.subtract(cumsum1, cumsum2)
    xs, ys = x_log, cumsum_diff[:-1]
    axes.plot(xs, ys, color='k', label='Model', rasterized=rasterize)
    axes.set_xlabel('MAPQ drop rate')
    if log2ize:
        axes.set_xticklabels(map(lambda x: 2 ** x, np.linspace(-1, -maxx, maxx)), rotation=90)
    axes.set_ylabel('Difference in cumulative # incorrect alignments')
    axes.set_title('Drop rate comparison on -log2 scale')
    axes.fill_between(xs, ys, where=ys>=0, interpolate=True, color='red')
    axes.fill_between(xs, ys, where=ys<=0, interpolate=True, color='green')
    axes.grid(True)
    return fig


def plot_drop_rate_v_squared_error(pcor, cor, pcor2=None, log2ize=False, rasterize=False):
    """  """
    cumsum = cum_squared_error(pcor, cor)
    cumsum2 = None if pcor2 is None else cum_squared_error(pcor2, cor)
    ln = len(pcor)
    x_log = np.linspace(1.0/ln, (ln - 1.0) / ln, ln - 1.0)
    if log2ize:
        x_log = -np.log2(x_log[::-1])
    maxx = math.ceil(np.log2(ln))
    fig = plt.figure(figsize=(10, 4))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.plot(x_log, cumsum[:-1], color='r', label='Model', rasterized=rasterize)
    if cumsum2 is not None:
        axes.plot(x_log, cumsum2[:-1], color='k', label='Original', rasterized=rasterize)
    axes.set_xlabel('MAPQ drop rate')
    if log2ize:
        axes.set_xticklabels(map(lambda x: 2 ** x, np.linspace(-1, -maxx, maxx)), rotation=90)
    axes.set_ylabel('Cumulative squared error')
    axes.set_title('Drop rate comparison on -log2 scale')
    axes.legend(loc=2)
    axes.grid(True)
    return fig


def bucket_error_plot(mapq_lists, labs, colors, cor, title=None):
    """ Plot predicted MAPQ versus actual MAPQ, both rounded to
        nearest integer.  Omit MAPQs where the number of incorrect
        alignments is 0 (making MAPQ infinity) """
    observed_mapq_lists = []
    actual_mapq_lists = []
    mx = 0
    for mapq_list in mapq_lists:
        # round MAPQs off to integers
        mapq = map(round, mapq_list)
        # stratify by mapq & count correct/incorrect
        tally = tally_cor_per(mapq, cor)
        observed_mapq_lists.append([])
        actual_mapq_lists.append([])
        for p, tup in sorted(tally.iteritems()):
            ncor, nincor = tup
            if nincor > 0:  # ignore buckets where all are correct
                ntot = ncor + nincor
                observed_mapq_lists[-1].append(p)
                actual_mapq = pcor_to_mapq(float(ncor) / ntot)
                actual_mapq_lists[-1].append(actual_mapq)
                mx = max(max(mx, p), actual_mapq)
    fig = plt.figure(figsize=(8, 8))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes.set_xlabel('Reported MAPQ')
    axes.set_ylabel('Actual MAPQ')
    if title is not None:
        axes.set_title(title)
    for o, a, l, c in zip(observed_mapq_lists, actual_mapq_lists, labs, colors):
        axes.scatter(o, a, label=l, color=c)
    axes.plot([0, mx], [0, mx], color='r')
    axes.set_xlim((0, mx))
    axes.set_ylim((0, mx))
    axes.legend(loc=2)
    axes.grid(True)
    return fig


def log2ize_p_minus(p, mx=10.0):
    """ Given a probability p, rescale it so that ps are mapped
        to values in [0, mx] and the right-hand end is magnified. """
    return abs(math.log(1 - (1 - (2 ** -mx)) * p, 2))


def unlog2ize_p_minus(p, mx=10.0):
    """ Inverse of log2ize_p_minus """
    return abs((2 ** -p) - 1)/(1 - (2 ** -mx))


def quantile_error_plot(mapq_lists, labs, colors, cor, title=None, quantiles=10,
                        already_sorted=False, exclude_mapq_gt=None, prob_scale=False,
                        log2ize=False):
    # One problem is how to determine the quantiles when there are multiple methods being compared
    estimated_lists = []
    actual_lists = []
    n_filtered = []
    mse_errors = []
    mx = 0
    scale = (lambda x: x)
    if log2ize:
        scale = log2ize_p_minus
    elif not prob_scale:
        scale = pcor_to_mapq
    for mapq_list in mapq_lists:
        # remove MAPQs greater than ceiling
        n_before = len(mapq_list)
        if exclude_mapq_gt is not None:
            mapq_list = filter(lambda x: x <= exclude_mapq_gt, mapq_list)
        n_filtered.append(n_before - len(mapq_list))
        # stratify by mapq & count correct/incorrect
        n = len(mapq_list)
        assert n >= quantiles
        if not already_sorted:
            sorted_order = sorted(range(n), key=lambda x: mapq_list[x])
            mapq_sorted = [mapq_list[x] for x in sorted_order]
            pcor_sorted = np.array(mapq_to_pcor_iter(mapq_sorted))
            cor_sorted = np.array([cor[x] for x in sorted_order])
        else:
            mapq_sorted, pcor_sorted, cor_sorted = mapq_list, mapq_to_pcor_iter(mapq_list), cor
        mse_errors.append(mse_error(pcor_sorted, cor_sorted))
        partition_starts = [int(round(i*n/quantiles)) for i in xrange(quantiles+1)]
        estimated_lists.append([])
        actual_lists.append([])
        for i in xrange(quantiles):
            st, en = partition_starts[i], partition_starts[i+1]
            ncor = sum(cor_sorted[st:en])
            ncor_est = sum(pcor_sorted[st:en])
            assert ncor_est < (en-st)
            estimated_lists[-1].append(scale(float(ncor_est)/(en-st)))
            if ncor == (en-st) and not prob_scale:
                actual_lists[-1].append(None)
            else:
                actual_lists[-1].append(scale(float(ncor)/(en-st)))
        mx = max(mx, max(max(estimated_lists[-1]), max(actual_lists[-1])))
    assert not math.isinf(mx)
    fig = plt.figure(figsize=(6, 6))
    axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    if prob_scale:
        axes.set_xlabel('Predicted probability correct')
        axes.set_ylabel('Actual fraction correct')
    else:
        axes.set_xlabel('Reported MAPQ')
        axes.set_ylabel('Actual MAPQ')
    if title is not None:
        axes.set_title(title)
    for o, a, l, c, mse in zip(estimated_lists, actual_lists, labs, colors, mse_errors):
        assert len(o) == len(a)
        a = [mx * 1.04 if m is None else m for m in a]
        axes.scatter(o, a, label='%s (MSE=%0.3f)' % (l, mse), color=c, alpha=0.5, s=60)
    axes.plot([-mx, 2*mx], [-mx, 2*mx], color='r')  # red line at y=x
    axes.set_xlim((-.02*mx, 1.02*mx))
    axes.set_ylim((-.02*mx, 1.02*mx))
    fracs = [0.0, 0.5, 0.9, 0.95, 0.99, 0.999, 1.0]
    axes.xaxis.set_ticks(map(log2ize_p_minus, fracs))
    axes.yaxis.set_ticks(map(log2ize_p_minus, fracs))
    axes.xaxis.set_ticklabels(map(str, fracs))
    axes.yaxis.set_ticklabels(map(str, fracs))
    axes.legend(loc=4)  # lower-right corner


def plot_subsampling_series(seriess, labs=None, colors=None):
    fig = plt.figure(figsize=(14, 12))

    if labs is None:
        labs = [str(i+1) for i in xrange(len(seriess))]
    if colors is None:
        colors = ['r', 'b', 'g', 'c', 'm'][:len(seriess)]
        assert len(colors) == len(seriess)

    ax1 = fig.add_subplot(3, 1, 1)
    min_rank_err, max_rank_err = float('inf'), 0.0
    for series, lab, color in zip(seriess, labs, colors):
        ax1.plot(series.fraction, series.rank_err_round, 'o-', color=color, label=lab)
        max_rank_err = max(max_rank_err, max(series.rank_err_round))
        min_rank_err = min(min_rank_err, min(series.rank_err_round))
    max_rank_err = min(2, max_rank_err)
    ax1.set_xlabel('Fraction')
    ax1.set_ylabel('Round rank err')
    ax1.set_ylim([min(0.8, min_rank_err * 0.98), max(1.2, max_rank_err * 1.02)])
    ax1.legend(loc=1)
    ax1.grid(True)

    ax2 = fig.add_subplot(3, 1, 2)
    max_mse = 0.0
    for series, lab, color in zip(seriess, labs, colors):
        ax2.plot(series.fraction, series.mse_err_round, 'o-', color=color, label=lab)
        max_mse = max(max_mse, max(series.mse_err_round))
    max_mse = min(2, max_mse)
    ax2.set_xlabel('Fraction')
    ax2.set_ylabel('MSE')
    ax2.set_ylim([0, max(1.0, max_mse * 1.02)])
    ax2.legend(loc=1)
    ax2.grid(True)

    ax3 = fig.add_subplot(3, 1, 3)
    max_mapq_avg = 0.0
    for series, lab, color in zip(seriess, labs, colors):
        ax3.plot(series.fraction, series.mapq_avg, 'o-', color=color, label=lab)
        max_mapq_avg = max(max_mapq_avg, max(series.mapq_avg))
    ax3.set_xlabel('Fraction')
    ax3.set_ylabel('Average MAPQ')
    ax3.set_ylim([0., max(80., max_mapq_avg * 1.02)])
    ax3.legend(loc=1)
    ax3.grid(True)

    plt.close()


def plot_fit(model, x_lim=(0.0, 1.0), y_lim=(0.0, 1.0), dx=0.01, dy=0.01, zmin=0.0, zmax=60.0):
    grid = np.mgrid[slice(y_lim[0], y_lim[1] + dy, dy),
                    slice(x_lim[0], x_lim[1] + dx, dx)]
    dim_x = 1 + ((x_lim[1] - x_lim[0]) / dx)
    dim_y = 1 + ((y_lim[1] - y_lim[0]) / dy)
    assert grid.shape == (2, dim_x, dim_y), "%s, (2, %d, %d)" % (str(grid.shape), dim_x, dim_y)
    xy = np.column_stack([grid[0].flatten(), grid[1].flatten()])
    z = np.array(pcor_to_mapq_iter(MapqFit.postprocess_predictions(model.predict(xy), '')), dtype=np.float32).reshape((dim_x, dim_y))

    # x and y are bounds, so z should be the value *inside* those bounds.
    # Therefore, remove the last value from the z array.
    z = z[:-1, :-1].transpose()
    z = z[::-1, :]
    if zmin is None:
        zmin = z.min()
    if zmax is None:
        zmax = z.max()

    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    cmap = plt.get_cmap('Blues')

    fitplt = plt.figure(figsize=(10, 10))
    axes = fitplt.add_subplot(111)

    extents = [x_lim[0], x_lim[1], y_lim[0], y_lim[1]]

    im = axes.imshow(z, cmap=cmap, extent=extents, vmax=zmax, vmin=zmin, origin='upper')

    cbar = fitplt.colorbar(im)
    cbar.set_label('MAPQ')

    # set the limits of the plot to the limits of the data
    axes.set_xlabel('Best score, 0=min, 1=max')
    axes.set_ylabel('Difference between best and second-best score, 0=min, 1=max')
    axes.set_title('Predicted MAPQ for X, Y in [0, 1]')
    return fitplt, axes


class MapqPredictions:

    def __init__(self):
        # all these lists are parallel
        self.pcor = np.array([])  # predicted pcors
        self.mapq_orig = []  # original mapping qualities
        self.category = []  # categories of alignments
        self.names = None  # names of reads
        self.data = None  # data that gave rise to predictions
        self.correct = None  # whether or not alignment is correct
        self.pcor_hist = None
        self.mapq = None
        self.correct_end, self.correct_run = 0, 0
        self.pcor_orig = None
        self.mapq_avg, self.mapq_orig_avg = 0., 0.
        self.mapq_std, self.mapq_orig_std = 0., 0.
        self.rank_err_orig = None
        self.rank_err = None
        self.rank_err_round = None
        self.mse_err_orig = None
        self.mse_err_raw = None
        self.mse_err = None
        self.mse_err_round = None

    def add_pcors(self, pcor, mapq_orig, category, names=None, data=None, correct=None):
        """ Add a new batch of predictions """
        self.pcor = np.append(self.pcor, pcor)
        self.mapq_orig.extend(mapq_orig)
        self.category.extend([category] * len(pcor))
        if data is not None:
            if self.data is None:
                self.data = []
            self.data.extend(data)
        if correct is not None:
            if self.correct is None:
                self.correct = correct
            else:
                self.correct = np.append(self.correct, correct)
        if names is not None:
            if self.names is None:
                self.names = []
            self.names.extend(names)

    def incorrect_indexes(self):
        """ Return indexes of in correct alignments in order
            from highest to lowest predicted pcor """
        assert self.correct is not None
        return [x for x in xrange(len(self.correct)-1, -1, -1) if not self.correct[x]]

    def top_incorrect(self, n=50):
        assert self.data is not None
        return [self.data[x] for x in self.incorrect_indexes()[:n]]

    def summarize_incorrect(self, n=50):
        assert self.correct is not None
        incor_idx = self.incorrect_indexes()[:n]
        summ_dict = dict()
        summ_dict['category'] = [self.category[x] for x in incor_idx]
        summ_dict['mapq'] = [self.mapq[x] for x in incor_idx]
        summ_dict['mapq_orig'] = [self.mapq_orig[x] for x in incor_idx]
        if self.names is not None:
            summ_dict['names'] = [self.names[x] for x in incor_idx]
        if self.data is not None:
            summ_dict['data'] = map(lambda x: ','.join(map(lambda y: '%0.3f' % y, x)),
                                    [self.data[x] for x in incor_idx])
        summ_dict['correct'] = [self.correct[x] for x in incor_idx]
        return pandas.DataFrame.from_dict(summ_dict)

    def finalize(self, verbose=False):

        # put pcor, mapq, mapq_orig, pcor_orig in order by pcor
        if verbose:
            logging.info('  Sorting data')
        pcor_order = sorted(range(len(self.pcor)), key=lambda x: self.pcor[x])
        self.pcor = self.pcor[pcor_order]
        self.mapq_orig = [self.mapq_orig[x] for x in pcor_order]
        self.category = [self.category[x] for x in pcor_order]
        if self.data is not None:
            self.data = [self.data[x] for x in pcor_order]
        if self.names is not None:
            self.names = [self.names[x] for x in pcor_order]

        # make pcor histogram
        self.pcor_hist = Counter(self.pcor)

        self.mapq = mapq = pcor_to_mapq_iter(self.pcor)
        self.pcor_orig = pcor_orig = mapq_to_pcor_iter(self.mapq_orig)

        pcor, mapq_orig = self.pcor, self.mapq_orig

        # calculate error measures and other measures
        if self.correct is not None:

            # calculate # of highest pcors and max # pcors in a row that
            # correspond to correct alignments
            if verbose:
                logging.info('  Finding correct runs')
            self.correct = correct = self.correct[pcor_order]
            end, run = True, 0
            for i in xrange(len(correct)-1, -1, -1):
                if correct[i] and end:
                    self.correct_end += 1
                elif end:
                    end = False
                run += 1 if correct[i] else -run
                self.correct_run = max(self.correct_run, run)

            # ranking error; +1 is to avoid division-by-zero when a dataset
            # is perfectly ranked
            if verbose:
                logging.info('  Calculating rank error')
            self.rank_err_orig = ranking_error(pcor_orig, correct) + 1
            self.rank_err = (ranking_error(pcor, correct) + 1) / self.rank_err_orig
            self.rank_err_round = (ranking_error(pcor, correct, rounded=True) + 1) / self.rank_err_orig
            if verbose:
                logging.info('    Done: %0.4f, %0.4fr' % (self.rank_err, self.rank_err_round))

            if verbose:
                logging.info('  Calculating MSE')
            self.mse_err_orig = mse_error(pcor_orig, correct)
            self.mse_err_raw = mse_error(pcor, correct)
            self.mse_err = self.mse_err_raw / self.mse_err_orig
            self.mse_err_round = mse_error(pcor, correct, rounded=True) / self.mse_err_orig
            if verbose:
                logging.info('    Done: %0.4f, %0.4fr (raw:%0.4f, orig raw:%0.4f)' % \
                             (self.mse_err, self.mse_err_round,
                              self.mse_err_raw, self.mse_err_orig))

        # summary statistics over pcors and mapqs
        if verbose:
            logging.info('  Calculating MAPQ summaries')
        self.mapq_avg, self.mapq_orig_avg = float(np.mean(mapq)), float(np.mean(mapq_orig))
        self.mapq_std, self.mapq_orig_std = float(np.std(mapq)), float(np.std(mapq_orig))

    def summary(self, include_distinct=False):
        """ Return a concise summary of prediction performance """
        if len(self.pcor) == 0:
            return ''
        ret = []
        if self.correct is not None:
            ret.append("rerr=%0.2f,%0.2fr" % (self.rank_err, self.rank_err_round))
            ret.append("mse=%0.2f,%0.2fr" % (self.mse_err, self.mse_err_round))
        ret.append("q=%0.2f,%0.2f,%0.2f" % (self.mapq_avg, self.mapq_std, max(self.mapq)))
        ret.append("end=%d,%0.2f%%" % (self.correct_end,
                                       100.0 * self.correct_end / len(self.pcor)))
        ret.append("cor=%d,%0.2f%%" % (self.correct_run,
                                       100.0 * self.correct_run / len(self.pcor)))
        if include_distinct:
            ret.append("distinct=%d,%0.2f%%" % (len(self.pcor_hist),
                                                100.0 * len(self.pcor_hist) / len(self.pcor)))
        return ' '.join(ret)


class MapqFit:

    @staticmethod
    def _df_to_mat(data, shortname):
        """ Convert a data frame read with read_dataset into a matrix
            suitable for use with scikit-learn, and parallel vectors
            giving the original MAPQ predictions and whether or not the
            alignments are correct. """
        if shortname == 'c':
            # extract relevant paired-end features from training data
            labs = ['best1',  # score of best alignment for mate
                    'best2',  # score of best alignment for opposite mate
                    'diff1',  # difference for mate
                    'diff2',  # difference for opposite mate
                    'diff_conc',  # concordant difference
                    'fraglen_z']  # # stddevs diff for fraglen
            mapq_header = 'mapq1'
        elif shortname == 'd':
            # extract relevant discordant paired-end features
            labs = ['best1', 'diff1']
            mapq_header = 'mapq1'
        else:
            # extract relevant unpaired features
            labs = ['best', 'diff']
            mapq_header = 'mapq'
        for lab in labs:
            assert not np.isnan(data[lab]).any()
        data_mat = data[labs].values
        correct = np.array(map(lambda x: x == 1, data['correct']))
        assert not np.isinf(data_mat).any() and not np.isnan(data_mat).any()
        return data_mat, data[mapq_header], correct

    @staticmethod
    def _downsample(x_train, mapq_orig_train, y_train, sample_fraction):
        """ Return a random subset of the data, MAPQs and labels.  Size of
            subset given by sample_fraction. """
        n_training_samples = x_train.shape[0]
        if sample_fraction < 1.0:
            sample_indexes = random.sample(xrange(n_training_samples), int(round(n_training_samples * sample_fraction)))
            x_train = x_train[sample_indexes, ]
            mapq_orig_train = mapq_orig_train[sample_indexes]
            y_train = np.array([y_train[i] for i in sample_indexes])
        return x_train, mapq_orig_train, y_train

    def summary(self):
        ret = []
        if self.predict_training:
            ret.append('train: ' + self.pred_overall_train.summary())
        ret.append(' test: ' + self.pred_overall.summary())
        return '\n'.join(ret)

    def pcor(self):
        """ Return list of predicted pcors (probability correct) """
        return self.pred_overall.pcor

    def pcor_orig(self):
        """ Return list of pcors (probability correct) reported by
            the aligner (converted from integer-rounded MAPQs) """
        return self.pred_overall.pcor_orig

    def mapq(self):
        """ Return list of predicted MAPQs """
        return self.pred_overall.mapq

    def mapq_orig(self):
        """ Return list of the original MAPQs reported by the aligner """
        return self.pred_overall.mapq_orig

    def correct(self):
        """ Return list indicating which alignments are correct; None if that
            information isn't available. """
        return self.pred_overall.correct

    def data(self):
        """ Return list of input data given to predictor. """
        return self.pred_overall.data

    def names(self):
        """ Return list of names. """
        return self.pred_overall.names

    def category(self):
        """ Return list of categories. """
        return self.pred_overall.category

    @staticmethod
    def postprocess_predictions(pcor_test, dataset_name, max_pcor=0.999999):
        """ Deal with pcors equal to 1.0, which cause infinite MAPQs """
        mn, mx = min(pcor_test), max(pcor_test)
        if mx >= 1.0:
            if mn == mx:
                logging.warning('All data points for %s are predicted correct; results unreliable' % dataset_name)
                pcor_test = [max_pcor] * len(pcor_test)
            max_noninf_pcor_test = max(filter(lambda x: x < 1.0, pcor_test))
            pcor_test = [max_noninf_pcor_test + 1e-6 if p >= 1.0 else p for p in pcor_test]
        return np.maximum(np.minimum(pcor_test, max_pcor), 0.)

    @staticmethod
    def _fit(mf_gen, x_train, y_train, dataset_shortname):
        """ Use cross validation to pick the best model from a
            collection of possible models (model_family) """
        mf = mf_gen()
        score_avgs, score_stds = [], []
        while True:
            params, pred = mf.next_predictor()
            if pred is None:
                break
            scores = cross_validation.cross_val_score(pred, x_train, y_train)
            score_avg, score_std = float(np.mean(scores)), float(np.std(scores))
            logging.debug("%s, avg=%0.3f, std=%0.3f, %s" % (dataset_shortname, score_avg, score_std, str(params)))
            score_avgs.append(score_avg)
            score_stds.append(score_std)
            mf.set_score(score_avg)
        best_params, best_pred = mf.best_predictor()
        logging.debug("BEST: %s, avg=%0.3f, %s" % (dataset_shortname, max(score_avgs), str(best_params)))
        assert best_pred is not None
        return best_pred, best_params, score_avgs, score_stds

    def __init__(self, dr, model_gen, sample_fraction=1.0,
                 sample_fractions=None,
                 random_seed=628599, predict_training=False,
                 verbose=False, keep_names=False, keep_data=False):

        def log_msg(s):
            if verbose:
                logging.info(s)

        self.sample_fraction = sample_fraction
        self.random_seed = random_seed
        self.predict_training = predict_training
        if sample_fractions is None:
            sample_fractions = {}

        # read datasets
        log_msg('Reading datasets')
        self.train_dfs = train_dfs = read_dataset(os.path.join(dr, 'training'))
        self.test_dfs = test_dfs = read_dataset(os.path.join(dr, 'test'))

        self.trained_models = {}
        self.crossval_avg = {}
        self.crossval_std = {}
        self.datasets = zip('dbcu',
                            ['Discordant', 'Bad-end', 'Concordant', 'Unpaired'],
                            [True, False, True, False])
        self.pred_overall = MapqPredictions()
        self.pred = {}
        self.pred_overall_train = MapqPredictions()
        self.pred_train = {}

        for dataset_shortname, dataset_name, paired in self.datasets:

            train = train_dfs[dataset_shortname]
            if train.shape[0] == 0:
                continue

            log_msg('Fitting %s training data' % dataset_name)

            # seed pseudo-random generators
            random.seed(random_seed)
            np.random.seed(random_seed)

            # extract features, convert to matrix
            x_train, mapq_orig_train, y_train = self._df_to_mat(train, dataset_shortname)

            # optionally downsample
            frac = sample_fractions[dataset_shortname] if dataset_shortname in sample_fractions else sample_fraction
            log_msg('Sampling %0.2f%% of %d rows of %s training data' % (100.0 * frac, train.shape[0], dataset_name))
            x_train, mapq_orig_train, y_train = \
                self._downsample(x_train, mapq_orig_train, y_train, frac)
            log_msg('  Now has %d rows' % x_train.shape[0])

            # use cross-validation to pick a model
            self.trained_models[dataset_shortname], params, avgs, stds = \
                self._fit(model_gen, x_train, y_train, dataset_shortname)

            self.trained_params = params
            self.crossval_avg[dataset_shortname] = max(avgs)

            # fit all the training data with the model
            self.trained_models[dataset_shortname].fit(x_train, y_train)

            names = None
            names_colname = 'name1' if paired else 'name'

            for training in [False, True]:

                if training and not predict_training:
                    continue

                log_msg('Making %s predictions for %s data' % (dataset_name, 'training' if training else 'test'))

                # set up some variables based on training/test
                test = train_dfs[dataset_shortname] if training else test_dfs[dataset_shortname]
                pred_overall = self.pred_overall_train if training else self.pred_overall
                pred = self.pred_train if training else self.pred
                pred[dataset_shortname] = MapqPredictions()

                # get test matrix
                x_test, mapq_orig_test, y_test = self._df_to_mat(test, dataset_shortname)
                y_test = np.array(map(lambda x: x == 1, test.correct))

                # predict
                pcor = self.trained_models[dataset_shortname].predict(x_test)  # make predictions
                pcor = np.array(self.postprocess_predictions(pcor, dataset_name))
                if keep_names and names_colname in test:
                    names = test[names_colname]
                data = None
                if keep_data:
                    data = x_test.tolist()
                for pr in [pred_overall, pred[dataset_shortname]]:
                    pr.add_pcors(pcor, mapq_orig_test, dataset_shortname, names=names, data=data, correct=y_test)

        # finalize the predictions
        log_msg('Finalizing results for overall test data (%d alignments)' % len(self.pred_overall.pcor))
        self.pred_overall.finalize(verbose=verbose)
        for shortname, pred in self.pred.iteritems():
            log_msg('Finalizing results for "%s" test data (%d alignments)' % (shortname, len(pred.pcor)))
            pred.finalize(verbose=verbose)
        if predict_training:
            log_msg('Finalizing results for overall training data (%d alignments)' % len(self.pred_overall_train))
            self.pred_overall_train.finalize(verbose=verbose)
            for shortname, pred in self.pred_train.iteritems():
                log_msg('Finalizing results for "%s" training data (%d alignments)' % (shortname, len(pred.pcor)))
                pred.finalize(verbose=verbose)
        log_msg('Done')


class ModelFamily(object):
    """ Encapsulates a model family and a simple interface for naively
        searching the space of hyperparameters. """

    def __init__(self, new_predictor, params):
        """
        new_predictor: function that takes set of params and returns
                       new predictor with those params
        params: list of lists
        """
        self.new_predictor = new_predictor
        self.params = params
        self.center = tuple([int(round(len(x)/2)) for x in params])
        self.scores = {}
        self.last_params = None
        self.neighborhood_scores = None
        self._setup_neighborhood(self.center)

    def _setup_neighborhood(self, center):
        self.neighborhood_scores = {}
        for i in xrange(len(self.params)):
            if center[i] > 0:
                neighbor = list(center[:])
                neighbor[i] -= 1
                neighbor = tuple(neighbor)
                self.neighborhood_scores[neighbor] = self.scores.get(neighbor, None)
            if center[i] < len(self.params[i])-1:
                neighbor = list(center[:])
                neighbor[i] += 1
                neighbor = tuple(neighbor)
                self.neighborhood_scores[neighbor] = self.scores.get(neighbor, None)

    def _best_neighbor(self):
        assert not self._has_unexplored_neighbors()
        assert len(self.neighborhood_scores) > 0
        prioritized = sorted([(v, k) for k, v in self.neighborhood_scores.iteritems()])
        return prioritized[-1][0], prioritized[-1][1]  # highest cross-val score

    def _has_unexplored_neighbors(self):
        return any(map(lambda x: x is None, self.neighborhood_scores.itervalues()))

    def _get_unexplored_neighbor(self):
        for k, v in self.neighborhood_scores.iteritems():
            if v is None:
                return k
        raise RuntimeError('_get_unexplored_neighbor called with no unexplored neighbors')

    def _idxs_to_params(self, idxs):
        return [self.params[i][j] for i, j in enumerate(idxs)]

    def next_predictor(self):
        if self.center not in self.scores:
            self.last_params = self.center
            return self._idxs_to_params(self.center), self.new_predictor(self.center)
        while True:
            if self._has_unexplored_neighbors():
                # more neighbors to be explored
                params = self._get_unexplored_neighbor()
                if params in self.scores:
                    self.neighborhood_scores[params] = self.scores[params]
                    continue
                self.last_params = params
                return self._idxs_to_params(params), self.new_predictor(params)
            else:
                # just finished exploring neighborhood
                assert self.center in self.scores
                best_score, best_neighbor = self._best_neighbor()
                if best_score > self.scores[self.center]:
                    # improved on center
                    self.center = best_neighbor
                    self._setup_neighborhood(self.center)
                    continue
                else:
                    # didn't improve on best!
                    self.last_params = None
                    return None, None

    def set_score(self, score):
        assert self.last_params is not None
        assert self.last_params not in self.scores
        self.scores[self.last_params] = score

    def best_predictor(self):
        assert len(self.scores) > 0
        prioritized = sorted([(v, k) for k, v in self.scores.iteritems()])
        return self._idxs_to_params(prioritized[-1][1]), self.new_predictor(prioritized[-1][1])


def random_forest_models(random_seed=33):
    # These perform OK but not as well as the extremely random trees
    def _gen(params):
        return RandomForestRegressor(n_estimators=params[0], max_depth=params[1],
                                     random_state=random_seed, max_features='auto')
    return lambda: ModelFamily(_gen, [range(5, 100, 10), range(3, 6)])


def extra_trees_models(random_seed=33):
    # These perform quite well
    def _gen(params):
        return ExtraTreesRegressor(n_estimators=params[0], max_depth=params[1],
                                   random_state=random_seed, max_features='auto')
    return lambda: ModelFamily(_gen, [range(5, 105, 5), range(3, 20)])


def gradient_boosting_models(random_seed=33):
    # These perform OK but not as well as the extremely random trees
    def _gen(params):
        return GradientBoostingRegressor(n_estimators=params[0], max_depth=params[1],
                                         learning_rate=params[2], random_state=random_seed,
                                         max_features='auto', loss='ls')
    return lambda: ModelFamily(_gen, [range(5, 40, 10), range(2, 6), [0.1, 0.2, 0.3, 0.5, 0.75, 1.0]])


def adaboost_models(random_seed=33):
    # These seem to perform quite poorly
    def _gen(params):
        return GradientBoostingRegressor(n_estimators=params[0], learning_rate=params[1],
                                         loss='linear', random_state=random_seed)
    return lambda: ModelFamily(_gen, [[10, 15, 20, 25, 30, 35, 40], [1.0, 1.5, 2.0]])


def model_family(args):
    if args.model_family == 'RandomForest':
        return random_forest_models()
    elif args.model_family == 'ExtraTrees':
        return extra_trees_models()
    elif args.model_family == 'GradientBoosting':
        return gradient_boosting_models()
    elif args.model_family == 'AdaBoost':
        return adaboost_models()
    else:
        raise RuntimeError('Bad value for --model-family: "%s"' % args.model_family)


def mkdir_quiet(dr):
    # Create output directory if needed
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def make_plots(pred, odir, args, prefix=''):
    mkdir_quiet(odir)
    fmat = args.plot_format
    if args.plot_cum_incorrect or args.plot_all:
        logging.info(prefix + 'Making cumulative-incorrect plots')
        assert pred.correct is not None
        plot_drop_rate(pred.pcor, pred.correct, pcor2=pred.pcor_orig, log2ize=False, rasterize=args.rasterize).savefig(
            os.path.join(odir, 'drop_rate.' + fmat))
        plt.close()
        assert len(plt.get_fignums()) == 0
        plot_drop_rate(pred.pcor, pred.correct, pcor2=pred.pcor_orig, log2ize=True, rasterize=args.rasterize).savefig(
            os.path.join(odir, 'drop_rate_log2.' + fmat))
        plt.close()
        assert len(plt.get_fignums()) == 0
        plot_drop_rate_difference(pred.pcor, pred.pcor_orig, pred.correct, log2ize=False, rasterize=args.rasterize).savefig(
            os.path.join(odir, 'drop_rate_diff.' + fmat))
        plt.close()
        assert len(plt.get_fignums()) == 0
        plot_drop_rate_difference(pred.pcor, pred.pcor_orig, pred.correct, log2ize=True, rasterize=args.rasterize).savefig(
            os.path.join(odir, 'drop_rate_diff_log2.' + fmat))
        plt.close()
        assert len(plt.get_fignums()) == 0
    if args.plot_mapq_buckets or args.plot_all:
        logging.info(prefix + 'Making MAPQ bucket plots')
        bucket_error_plot([pred.mapq, pred.mapq_orig], ['Predicted', 'Original'], ['b', 'g'], pred.correct).savefig(
            os.path.join(odir, 'mapq_buckets.' + fmat))
        plt.close()
        assert len(plt.get_fignums()) == 0


def go(args):
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG)

    logging.info('Instantiating model family')
    fam = model_family(args)

    odir = args.output_directory
    mkdir_quiet(odir)

    if args.subsampling_series is not None:
        logging.info('Doing subsampling series')
        ss_odir = os.path.join(odir, 'subsampled')
        fractions = []
        preds = [list() for _ in xrange(args.subsampling_replicates)]
        fits = [list() for _ in xrange(args.subsampling_replicates)]
        for fraction in map(float, args.subsampling_series.split(',')):
            logging.info('  Fraction=%0.3f' % fraction)
            for repl in xrange(1, args.subsampling_replicates+1):
                my_seed = hash(str(args.seed + repl))
                gc.collect()
                logging.info('    Replicate=%d' % repl)
                my_odir = os.path.join(ss_odir, '%0.3f' % fraction, str(repl))
                mkdir_quiet(my_odir)
                my_fit_fn = os.path.join(my_odir, 'fit.pkl')
                if os.path.exists(my_fit_fn) and not args.overwrite_fit:
                    logging.info('      Loading predictions from file')
                    with open(my_fit_fn, 'rb') as fh:
                        ss_fit = cPickle.load(fh)
                else:
                    logging.info('      Fitting and making predictions')
                    ss_fit = MapqFit(args.input_directory, fam, fraction, verbose=args.verbose, random_seed=my_seed)
                    if args.serialize_fit:
                        logging.info('      Serializing fit object')
                        with open(my_fit_fn, 'wb') as ofh:
                            cPickle.dump(ss_fit, ofh, 2)
                fits[repl-1].append(ss_fit)
                preds[repl-1].append(ss_fit.pred_overall)
                logging.info('      Making plots')
                make_plots(ss_fit.pred_overall, my_odir, args, prefix='        ')
            gc.collect()
            fractions.append(fraction)
        perf_dicts = [{'fraction': fractions,
                       'rank_err': [x.rank_err for x in pred],
                       'rank_err_round': [x.rank_err_round for x in pred],
                       'mse_err': [x.mse_err for x in pred],
                       'mse_err_round': [x.mse_err_round for x in pred],
                       'mapq_avg': [x.mapq_avg for x in pred],
                       'mapq_std': [x.mapq_std for x in pred],
                       'params': [str(x.trained_params) for x in fit]} for pred, fit in zip(preds, fits)]
        dfs = [pandas.DataFrame.from_dict(perf_dict) for perf_dict in perf_dicts]
        for i, df in enumerate(dfs):
            df.to_csv(os.path.join(odir, 'subsampling_series_%d.tsv' % (i+1)), sep='\t', index=False)
        plot_subsampling_series(dfs).savefig(os.path.join(odir, 'subsampling_series.' + args.plot_format))
        plt.close()
        assert len(plt.get_fignums()) == 0

    fit_fn = os.path.join(odir, 'fit.pkl')
    if os.path.exists(fit_fn) and not args.overwrite_fit:
        logging.info('Loading fit from file')
        with open(fit_fn, 'rb') as fh:
            fit = cPickle.load(fh)
    else:
        logging.info('Fitting and making predictions')
        fit = MapqFit(args.input_directory, fam, args.subsampling_fraction, verbose=args.verbose, random_seed=args.seed)
        print fit.summary()
        if args.serialize_fit:
            logging.info('Serializing fit object')
            with open(fit_fn, 'wb') as ofh:
                cPickle.dump(fit, ofh, 2)

    logging.info('Making plots')
    make_plots(fit.pred_overall, odir, args, prefix='  ')
    logging.info('Done')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Fit model, make predictions.')

    # Output-related options
    parser.add_argument('--input-directory', metavar='path', type=str, required=True,
                        help='Directory with output from tandem simulator')
    parser.add_argument('--output-directory', metavar='path', type=str, required=True,
                        help='Directory to write summaries and plots to')

    # Model-related options
    parser.add_argument('--model-family', metavar='family', type=str, required=False,
                        default='ExtraTrees', help='{RandomForest | ExtraTrees | GradientBoosting | AdaBoost}')
    parser.add_argument('--overwrite-fit', action='store_const', const=True, default=False,
                        help='Re-fit the model even if a fit is already present in --output-directory')
    parser.add_argument('--serialize-fit', action='store_const', const=True, default=False,
                        help='Write fit model to a pickle file')
    parser.add_argument('--compression-effort', metavar='int', type=int, default=1,
                        help='How hard to try to compress the model when writing to .pkl file')

    # What to generate
    parser.add_argument('--subsampling-fraction', metavar='float', type=float, default=1.0,
                        help='Subsample the training down to this fraction before fitting model')
    parser.add_argument('--subsampling-series', metavar='floats', type=str,
                        help='Comma separated list of subsampling fractions to try')
    parser.add_argument('--subsampling-replicates', metavar='int', type=int, default=1,
                        help='Number of times to repeat fiting/prediction for each subsampling fraction')
    parser.add_argument('--plot-cum-incorrect', action='store_const', const=True, default=False,
                        help='Make cumulative-incorrect plots, including on -log and normal scale, and for predicted '
                             'and difference')
    parser.add_argument('--plot-mapq-buckets', action='store_const', const=True, default=False,
                        help='Plot expected vs actual MAPQ')
    parser.add_argument('--plot-all', action='store_const', const=True, default=False,
                        help='Like specifying all option beginning with --plot')
    parser.add_argument('--plot-format', metavar='format', type=str, default='png',
                        help='Extension (and image format) for plot: {pdf, png, eps, jpg, ...}')

    # Plotting options
    # This doesn't seem to help anything!
    parser.add_argument('--rasterize', action='store_const', const=True, default=False,
                        help='Rasterize complicated lines and curves in plots')

    # Other options
    parser.add_argument('--seed', metavar='int', type=int, default=6277, help='Pseudo-random seed')
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Run unit tests')
    parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
    parser.add_argument('--profile', action='store_const', const=True, default=False, help='Output profiling data')

    if '--version' in sys.argv:
        print 'Tandem predictor, version ' + VERSION
        sys.exit(0)

    if '--test' in sys.argv:
        import unittest

        class Test(unittest.TestCase):
            def test_model_family(self):
                def model_gen(params):
                    assert len(params) == 3
                    return lambda: -sum(map(lambda x: x**2, params))
                mf = ModelFamily(model_gen, [[-2, -1, 0, 1, 2, 3, 4], [-5, -4, -3, -2, -1, 0], [-1, 0, 10]])
                while True:
                    params, fn = mf.next_predictor()
                    if fn is None:
                        break
                    mf.set_score(fn())
                best_params, best_pred = mf.best_predictor()
                self.assertEqual(0, best_pred())

        unittest.main(argv=[sys.argv[0]])
        sys.exit(0)

    myargs = parser.parse_args()

    pr = None
    if myargs.profile:
        import cProfile
        import pstats
        import StringIO
        pr = cProfile.Profile()
        pr.enable()
    go(myargs)
    if myargs.profile:
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(30)
        print s.getvalue()