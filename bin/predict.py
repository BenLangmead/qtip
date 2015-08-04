"""
Given a directory with output from ts.py, predict new MAPQs.

TODO: normalizers that deal naturally with various-length reads
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
import copy
from itertools import imap, izip
from sklearn import cross_validation
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from collections import defaultdict, Counter


VERSION = '0.1.0'


def pcor_to_mapq_np(pcor):
    old = np.seterr(divide='ignore')
    ret = np.abs(-10.0 * np.log10(1.0 - pcor))
    np.seterr(**old)
    return ret


def mapq_to_pcor_np(mapq):
    return 1.0 - 10.0 ** (-0.1 * mapq)


def round_pcor_np(pcor):
    return mapq_to_pcor_np(np.round(pcor_to_mapq_np(pcor)))


def pcor_to_mapq(p):
    """ Convert probability correct (pcor) to mapping quality (MAPQ) """
    return int(round(abs(-10.0 * math.log10(1.0 - p)) if p < 1.0 else float('inf')))


def mapq_to_pcor(p):
    """ Convert mapping quality (MAPQ) to probability correct (pcor) """
    return (1.0 - (10.0 ** (-0.1 * p))) if p < float('inf') else 1.0


class Normalizers(object):
    """ Values that help us to normalize a new dataset prior to making MAPQ
        predictions.  "Normalize" here means that the alignment scores and
        fragment lengths are converted to standardized units.  There is a
        weakness here: when a new dataset contains reads of various lengths,
        we don't necessarily want simple, scalar min/max alignment scores
        here, since these will normalize reads of various lengths as though
        they all live in the same range of relevant alignment scores. """

    def __init__(self):
        self.minv = None
        self.maxv = None
        self.minv_1 = None
        self.maxv_1 = None
        self.minv_2 = None
        self.maxv_2 = None


class AlignmentTableReader(object):

    """ Reads a table of information describing alignments.  These are tables
        output by ts.py.  Tables might describe training or test alignments.
    """

    #            short name   suffix
    datasets = [('u',         '_unp.csv'),
                ('d',         '_disc.csv'),
                ('c',         '_conc.csv'),
                ('b',         '_bad_end.csv')]

    def __init__(self, prefix, use_normalizers=True, learn_normalizers=False, normalizers=None):
        self.prefix = prefix
        self.normalizers = {}
        self.dfs = {}
        self.readers = {}
        self.use_normalizers = use_normalizers
        assert learn_normalizers or normalizers is not None or not use_normalizers
        if learn_normalizers:
            assert use_normalizers
            for sn, suf in self.datasets:
                fn = self.prefix + suf
                if any(map(os.path.exists, [fn, fn + '.gz', fn + '.bz2'])):
                    self.dfs[sn] = self._fn_to_iterator(fn, chunksize=None)

                    def _new_iter(_sn):
                        def _inner():
                            return iter([self.dfs[_sn].copy()])
                        return _inner

                    self.readers[sn] = _new_iter(sn)
            self._extract_normalizers(self.dfs)
        else:
            self.normalizers = normalizers
            for sn, suf in self.datasets:
                fn = self.prefix + suf
                if any(map(os.path.exists, [fn, fn + '.gz', fn + '.bz2'])):

                    def _new_iter(_fn):
                        def _inner():
                            return self._fn_to_iterator(_fn)
                        return _inner

                    self.readers[sn] = _new_iter(fn)

    @staticmethod
    def _fn_to_iterator(fn, chunksize=50000):
        if os.path.exists(fn):
            return pandas.io.parsers.read_csv(fn, quoting=2, chunksize=chunksize)
        elif os.path.exists(fn + '.gz'):
            return pandas.io.parsers.read_csv(fn + '.gz', quoting=2, chunksize=chunksize, compression='gzip')
        elif os.path.exists(fn + '.bz2'):
            return pandas.io.parsers.read_csv(fn + '.bz2', quoting=2, chunksize=chunksize, compression='bz2')
        else:
            raise RuntimeError('No such file: "%s"' % fn)

    def _extract_normalizers(self, dfs):
        """ Calculates normalization factors based on the alignment table.
            Stores results in self.normalizers dictionary. """
        for sn, df in dfs.iteritems():
            if df is None:
                continue
            if sn == 'c':
                norm_conc = self.normalizers['c'] = Normalizers()
                norm_conc.minv_1 = df['best1_1'].min()
                norm_conc.maxv_1 = df['best1_1'].max()
                norm_conc.minv_2 = df['best1_2'].min()
                norm_conc.maxv_2 = df['best1_2'].max()
            elif sn == 'd':
                norm_disc = self.normalizers['d'] = Normalizers()
                norm_disc.minv_1 = df['best1_1'].min()
                norm_disc.maxv_1 = df['best1_1'].max()
            elif sn in 'bu':
                norm = self.normalizers[sn] = Normalizers()
                norm.minv = df['best1'].min()
                norm.maxv = df['best1'].max()
            else:
                raise RuntimeError('Bad shortname: "%s"' % sn)

    def _postprocess_data_frame(self, df, sn):
        """ Applies normalization to a data frame containing a chunk of rows
            from the alignment table.  Requires that we have normalizers,
            which may or may not have been learned from this alignment table.
        """
        def _normalize_and_standardize(_df, nm, best_nm, secbest_nm, mn, diff):
            _df[nm] = (_df[best_nm].astype(float) - _df[secbest_nm]) / diff
            _df[nm] = _df[nm].fillna(np.nanmax(_df[nm])).fillna(0)
            _df[best_nm] = (_df[best_nm].astype(float) - mn) / diff
            _df[best_nm] = _df[best_nm].fillna(np.nanmax(_df[best_nm])).fillna(0)
            assert not any([math.isnan(x) for x in _df[nm]])
            assert not any([math.isnan(x) for x in _df[best_nm]])

        def _standardize(_df, nm, best_nm, secbest_nm):
            _df[nm] = _df[best_nm].astype(float) - _df[secbest_nm]
            _df[nm] = _df[nm].fillna(np.nanmax(_df[nm])).fillna(0)
            assert not any([math.isnan(x) for x in _df[nm]])

        if df.shape[0] == 0:
            return

        if self.use_normalizers:
            assert sn in self.normalizers
            norm = self.normalizers[sn]

        # Turn the correct column into 0/1
        if df['correct'].count() == len(df['correct']):
            df['correct'] = df['correct'].map(lambda x: 1 if x == 'T' else 0)

        if sn in 'c':
            # normalize alignment scores
            if self.use_normalizers:
                minv_1 = df['minv_1'].fillna(norm.minv_1)
                maxv_1 = df['maxv_1'].fillna(norm.maxv_1)
                minv_2 = df['minv_2'].fillna(norm.minv_2)
                maxv_2 = df['maxv_2'].fillna(norm.maxv_2)
                _normalize_and_standardize(df, 'diff_1', 'best1_1', 'best2_1', minv_1, maxv_1 - minv_1)
                _normalize_and_standardize(df, 'diff_2', 'best1_2', 'best2_2', minv_2, maxv_2 - minv_2)
            else:
                _standardize(df, 'diff_1', 'best1_1', 'best2_1')
                _standardize(df, 'diff_2', 'best1_2', 'best2_2')

            # normalize concordant alignment scores
            if math.isnan(df['best1conc'].sum()) or math.isnan(df['best2conc'].sum()):
                # assume the worst
                logging.warning('Difference of concordants not available, so using minimum of mates')
                df['diff_conc'] = df[['diff_1', 'diff_2']].min(axis=1)
            elif self.use_normalizers:
                minconc, maxconc = norm.minv_1 + norm.minv_2, norm.maxv_1 + norm.maxv_2
                _normalize_and_standardize(df, 'diff_conc', 'best1conc', 'best2conc', minconc, maxconc - minconc)
            else:
                _standardize(df, 'diff_conc', 'best1conc', 'best2conc')

        elif sn == 'd':
            if self.use_normalizers:
                minv_1 = df['minv_1'].fillna(norm.minv_1)
                maxv_1 = df['maxv_1'].fillna(norm.maxv_1)
                _normalize_and_standardize(df, 'diff_1', 'best1_1', 'best2_1', minv_1, maxv_1 - minv_1)
            else:
                _standardize(df, 'diff_1', 'best1_1', 'best2_1')

        else:
            assert sn in 'ub'
            if self.use_normalizers:
                minv = df['minv'].fillna(norm.minv)
                maxv = df['maxv'].fillna(norm.maxv)
                _normalize_and_standardize(df, 'diff', 'best1', 'best2', minv, maxv - minv)
            else:
                _standardize(df, 'diff', 'best1', 'best2')

        return df

    def dataset_iter(self, sn):
        assert sn in self.readers
        return imap(lambda x: self._postprocess_data_frame(x, sn), self.readers[sn]())

    def __contains__(self, o):
        return o in self.readers


def tuples_to_unpaired_matrix(tups, normalizers):
    """ Convert a list of UnpairedTuples for unpaired alignments into a
        matrix suitable for use with a predictor in predict.py. """
    pass


def tuples_to_bad_end_matrix(tups, normalizers):
    """ Convert a list of UnpairedTuples for bad-end alignments into a
        matrix suitable for use with a predictor in predict.py. """
    pass


def tuples_to_concordant_matrix(ptups, normalizers):
    """ Convert a list of PairedTuples for concordant alignments into a matrix
        suitable for use with a predictor in predict.py. """
    pass


def tuples_to_discordant_matrix(ptups, normalizers):
    """ Convert a list of PairedTuples for discordant alignments into a matrix
        suitable for use with a predictor in predict.py. """
    pass


def tally_cor_per(level, cor):
    """ Some predictions from our model are the same; this helper function
        gathers all the correct/incorrect information for each group of equal
        predictions. """
    tally = defaultdict(lambda: [0, 0])
    for p, c in izip(level, cor):
        c = 0 if c else 1
        tally[p][c] += 1
    return tally


def auc(pcor, cor, rounded=False):
    """ Calculate area under curve given predictions and correct/incorrect
        information. """
    if rounded:
        pcor = round_pcor_np(pcor)
    area, tot_cor, tot_incor = 0, 0, 0
    last_tot_cor, last_tot_incor = 0, 0
    for pcor, ci in sorted(tally_cor_per(pcor, cor).iteritems(), reverse=True):
        tot_cor += ci[0]
        tot_incor += ci[1]
        cor_diff = tot_cor - last_tot_cor
        incor_diff = tot_incor - last_tot_incor
        if incor_diff > 0:
            area += (0.5 * cor_diff * incor_diff)
            area += last_tot_cor * incor_diff
        last_tot_cor, last_tot_incor = tot_cor, tot_incor
    return area


def roc_table(pcor, cor, rounded=False, mapqize=False):
    """ Return the ranking error given a list of pcors and a parallel list of
        correct/incorrect booleans.  Round off to nearest MAPQ first if
        rounded=True.  """
    assert len(pcor) == len(cor)
    if rounded:
        pcor = round_pcor_np(pcor)
    tally = tally_cor_per(pcor, cor)
    cum_cor, cum_incor = 0, 0
    mapqs, cors, incors, cum_cors, cum_incors = [], [], [], [], []
    for p in sorted(tally.iterkeys(), reverse=True):
        ncor, nincor = tally[p]
        cum_cor += ncor
        cum_incor += nincor
        mapqs.append(pcor_to_mapq(p) if mapqize else p)
        cors.append(ncor)
        incors.append(nincor)
        cum_cors.append(cum_cor)
        cum_incors.append(cum_incor)
    return pandas.DataFrame.from_dict({'mapq': mapqs, 'cor': cors, 'incor': incors,
                                       'cum_cor': cum_cors, 'cum_incor': cum_incors})


def ranking_error(pcor, cor, rounded=False):
    """ Return the ranking error given a list of pcors and a parallel list of
        correct/incorrect booleans.  Round off to nearest MAPQ first if
        rounded=True.  """
    assert len(pcor) == len(cor)
    if rounded:
        pcor = round_pcor_np(pcor)
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


def mseor(pcor, cor, rounded=False):
    """ Return the mean squared error between the pcors (in [0, 1]) and and
        cors (in {0, 1}).  This is a measure of how close our predictions are
        to true correct/incorrect. """
    if rounded:
        pcor = round_pcor_np(pcor)
    return np.mean((pcor - cor) ** 2) / np.mean((np.mean(cor) - cor) ** 2)


def cum_squared_error(pcor, cor, rounded=False):
    """ Return the cumulative squared error between pcors and cors, sorted
        from highest to lowest pcor. """
    if rounded:
        pcor = round_pcor_np(pcor)
    pcor_order = sorted(range(len(pcor)), key=lambda x: pcor[x], reverse=True)
    pcor = np.array([pcor[x] for x in pcor_order])
    cor = np.array([cor[x] for x in pcor_order])
    assert all(pcor[i] >= pcor[i+1] for i in xrange(len(pcor)-1))
    return np.cumsum((pcor - cor) ** 2)


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
    mse = []
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
            pcor_sorted = np.array(mapq_to_pcor_np(mapq_sorted))
            cor_sorted = np.array([cor[x] for x in sorted_order])
        else:
            mapq_sorted, pcor_sorted, cor_sorted = mapq_list, mapq_to_pcor_np(mapq_list), cor
        mse.append(mseor(pcor_sorted, cor_sorted))
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
    for o, a, l, c, mse in zip(estimated_lists, actual_lists, labs, colors, mse):
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

    ax1 = fig.add_subplot(4, 1, 1)
    min_rank_err, max_rank_err = float('inf'), float('-inf')
    for series, lab, color in zip(seriess, labs, colors):
        ax1.plot(series.fraction, series.rank_err_diff_round_pct, 'o-', color=color, label=lab)
        max_rank_err = max(max_rank_err, max(series.rank_err_diff_round_pct))
        min_rank_err = min(min_rank_err, min(series.rank_err_diff_round_pct))
    ax1.axhline(y=0, linewidth=2, color='k')
    ax1.set_xlabel('Subsampling fraction')
    ax1.set_ylabel('% diff, ranking error')
    ax1.set_ylim([min(-20.0, min_rank_err), max(20.0, max_rank_err)])
    ax1.legend(loc=1)
    ax1.grid(True)

    ax2 = fig.add_subplot(4, 1, 2)
    min_mse, max_mse = float('inf'), float('-inf')
    for series, lab, color in zip(seriess, labs, colors):
        ax2.plot(series.fraction, series.mse_diff_round_pct, 'o-', color=color, label=lab)
        max_mse = max(max_mse, max(series.mse_diff_round_pct))
        min_mse = min(min_mse, min(series.mse_diff_round_pct))
    ax2.axhline(y=0, linewidth=2, color='k')
    ax2.set_xlabel('Subsampling fraction')
    ax2.set_ylabel('% diff, MSE')
    ax2.set_ylim([min(-50.0, min_mse), max(50.0, max_mse)])
    ax2.legend(loc=1)
    ax2.grid(True)

    ax3 = fig.add_subplot(4, 1, 3)
    min_auc, max_auc = float('inf'), float('-inf')
    for series, lab, color in zip(seriess, labs, colors):
        ax3.plot(series.fraction, series.auc_diff_round_pct, 'o-', color=color, label=lab)
        max_auc = max(max_auc, max(series.auc_diff_round_pct))
        min_auc = min(min_auc, min(series.auc_diff_round_pct))
    ax3.axhline(y=0, linewidth=2, color='k')
    ax3.set_xlabel('Subsampling fraction')
    ax3.set_ylabel('% diff, AUC')
    ax3.set_ylim([min(-1.1, min_auc), max(1.1, max_auc)])
    ax3.legend(loc=1)
    ax3.grid(True)

    ax4 = fig.add_subplot(4, 1, 4)
    max_mapq_avg = 0.0
    for series, lab, color in zip(seriess, labs, colors):
        ax4.plot(series.fraction, series.mapq_avg, 'o-', color=color, label=lab)
        max_mapq_avg = max(max_mapq_avg, max(series.mapq_avg))
    ax4.axhline(y=0, linewidth=2, color='k')
    ax4.set_xlabel('Subsampling fraction')
    ax4.set_ylabel('Average MAPQ')
    ax4.set_ylim([0., max(80., max_mapq_avg * 1.02)])
    ax4.legend(loc=1)
    ax4.grid(True)

    return fig


def plot_fit(model, x_lim=(0.0, 1.0), y_lim=(0.0, 1.0), dx=0.01, dy=0.01, zmin=0.0, zmax=60.0):
    grid = np.mgrid[slice(y_lim[0], y_lim[1] + dy, dy),
                    slice(x_lim[0], x_lim[1] + dx, dx)]
    dim_x = 1 + ((x_lim[1] - x_lim[0]) / dx)
    dim_y = 1 + ((y_lim[1] - y_lim[0]) / dy)
    assert grid.shape == (2, dim_x, dim_y), "%s, (2, %d, %d)" % (str(grid.shape), dim_x, dim_y)
    xy = np.column_stack([grid[0].flatten(), grid[1].flatten()])
    z = np.array(pcor_to_mapq_np(MapqFit.postprocess_predictions(model.predict(xy), '')), dtype=np.float32).reshape((dim_x, dim_y))

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
    """ Encpsulates mapq predictions for a dataset for evaluation purposes. """

    def __init__(self):
        # all these lists are parallel
        self.pcor = np.array([])  # predicted pcors
        self.mapq_orig = np.array([])  # original mapping qualities
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
        self.rank_err_raw = None
        self.rank_err_raw_round = None
        self.rank_err_diff = None
        self.rank_err_diff_pct = None
        self.rank_err_diff_round = None
        self.rank_err_diff_round_pct = None
        self.auc_orig = None
        self.auc_raw = None
        self.auc_raw_round = None
        self.auc_diff = None
        self.auc_diff_pct = None
        self.auc_diff_round = None
        self.auc_diff_round_pct = None
        self.mse_orig = None
        self.mse_raw = None
        self.mse_raw_round = None
        self.mse_diff = None
        self.mse_diff_pct = None
        self.mse_diff_round = None
        self.mse_diff_round_pct = None
        self.mse = None
        self.mse_round = None

    def add_pcors(self, pcor, mapq_orig, category, names=None, data=None, correct=None):
        """ Add a new batch of predictions """
        self.pcor = np.append(self.pcor, pcor)
        self.mapq_orig = np.append(self.mapq_orig, mapq_orig)
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
        if verbose:
            logging.info('  Histogramming pcors')
        self.pcor_hist = Counter(self.pcor)

        if verbose:
            logging.info('  Sorting data')
        pcor_order = np.argsort(self.pcor)

        if verbose:
            logging.info('  Reordering data')
        self.pcor = self.pcor[pcor_order]
        self.mapq_orig = self.mapq_orig[pcor_order]
        self.category = [self.category[x] for x in pcor_order]
        if self.data is not None:
            self.data = [self.data[x] for x in pcor_order]
        if self.names is not None:
            self.names = [self.names[x] for x in pcor_order]

        if verbose:
            logging.info('  Converting between pcor and mapq')
        self.mapq = mapq = np.abs(-10.0 * np.log10(1.0 - self.pcor))
        self.pcor_orig = pcor_orig = 1.0 - 10.0 ** (-0.1 * self.mapq_orig)

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
            self.rank_err_raw = ranking_error(pcor, correct) + 1
            self.rank_err_raw_round = ranking_error(pcor, correct, rounded=True) + 1
            self.rank_err_diff = self.rank_err_raw - self.rank_err_orig
            self.rank_err_diff_pct = 100.0 * self.rank_err_diff / self.rank_err_orig
            self.rank_err_diff_round = self.rank_err_raw_round - self.rank_err_orig
            self.rank_err_diff_round_pct = 100.0 * self.rank_err_diff_round / self.rank_err_orig
            self.rank_err = self.rank_err_raw / self.rank_err_orig
            self.rank_err_round = self.rank_err_raw_round / self.rank_err_orig
            if verbose:
                logging.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.rank_err_diff_pct,
                                                                       self.rank_err_diff_round_pct))

            if verbose:
                logging.info('  Calculating AUC')
            self.auc_orig = auc(pcor_orig, correct)
            self.auc_raw = auc(pcor, correct)
            self.auc_raw_round = auc(pcor, correct, rounded=True)
            self.auc_diff = self.auc_raw - self.auc_orig
            self.auc_diff_round = self.auc_raw_round - self.auc_orig
            if self.auc_orig == 0.:
                if self.auc_diff > 0.:
                    self.auc_diff_pct = float('inf')
                else:
                    self.auc_diff_pct = 0.0
                if self.auc_diff_round > 0.:
                    self.auc_diff_round_pct = float('inf')
                else:
                    self.auc_diff_round_pct = 0
            else:
                self.auc_diff_pct = 100.0 * self.auc_diff / self.auc_orig
                self.auc_diff_round_pct = 100.0 * self.auc_diff_round / self.auc_orig
            if verbose:
                logging.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.auc_diff_pct, self.auc_diff_round_pct))

            if verbose:
                logging.info('  Calculating MSE')
            self.mse_orig = mseor(pcor_orig, correct)
            self.mse_raw = mseor(pcor, correct)
            self.mse_raw_round = mseor(pcor, correct, rounded=True)
            self.mse_diff = self.mse_raw - self.mse_orig
            self.mse_diff_pct = 100.0 * self.mse_diff / self.mse_orig
            self.mse_diff_round = self.mse_raw_round - self.mse_orig
            self.mse_diff_round_pct = 100.0 * self.mse_diff_round / self.mse_orig
            self.mse = self.mse_raw / self.mse_orig
            self.mse_round = mseor(pcor, correct, rounded=True) / self.mse_orig
            if verbose:
                logging.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.mse_diff_pct, self.mse_diff_round_pct))

        # summary statistics over pcors and mapqs
        if verbose:
            logging.info('  Calculating MAPQ summaries')
        self.mapq_avg, self.mapq_orig_avg = float(np.mean(mapq)), float(np.mean(mapq_orig))
        self.mapq_std, self.mapq_orig_std = float(np.std(mapq)), float(np.std(mapq_orig))


def parse_sam(fh, ofh, alignment_class, normalizers, fit, ival=10000):

    def _rewrite(_ln, _mapq):
        toks = _ln.split()
        toks[4] = str(round(_mapq))
        # TODO: could optionally include decimal MAPQ in extra field
        ofh.write('\t'.join(toks) + '\n')

    last_al, last_ln = None, None
    nal, nignored, npair, nunp = 0, 0, 0, 0
    for ln in fh:
        if ln[0] == '@':
            ofh.write(ln)
            continue  # skip headers

        # Parse SAM record
        al = alignment_class()
        al.parse(ln)
        nal += 1
        if al.flags >= 2048:
            nignored += 1
            ofh.write(ln)
            continue
        elif al.paired:
            npair += 1
        else:
            nunp += 1

        if (nal % ival) == 0:
            logging.info('      # alignments parsed: %d (%d paired, %d unpaired, %d ignored)' %
                         (nal, npair, nunp, nignored))

        if al.paired and last_al is not None:
            assert al.concordant == last_al.concordant
            assert al.discordant == last_al.discordant
            assert al.mate1 != last_al.mate1
            mate1, mate2 = (al, last_al) if al.mate1 else (last_al, al)
            ln1, ln2 = (ln, last_ln) if al.mate1 else (last_ln, ln)
            # if only one end aligned, ensure 'al' is the one that aligned
            if not al.is_aligned() and last_al.is_aligned():
                al, last_al = last_al, al
                ln, last_ln = last_ln, ln
        else:
            mate1, mate2 = None, None

        if al.is_aligned():
            if mate1 is not None:
                if al.concordant:
                    mapq1, mapq2 = fit.predict_concordant(mate1, mate2)
                    ofh.write(_rewrite(ln1, mapq1))
                    ofh.write(_rewrite(ln2, mapq2))
                elif al.discordant:
                    mapq1, mapq2 = fit.predict_discordant(mate1, mate2)
                    ofh.write(_rewrite(ln1, mapq1))
                    ofh.write(_rewrite(ln2, mapq2))
                else:
                    mapq = fit.predict_bad_end(al, last_al)
                    ofh.write(_rewrite(ln, mapq))
            elif not al.paired:
                mapq = fit.predict_unpaired(al)
                ofh.write(_rewrite(ln, mapq))
        else:
            ofh.write(ln)

        if mate1 is None and al.paired:
            last_al, last_ln = al, ln
        else:
            last_al, last_ln = None, None


class MapqFit:

    @staticmethod
    def _df_to_mat(data, shortname, remove_labels=None, include_ztzs=True):
        """ Convert a data frame read with read_dataset into a matrix
            suitable for use with scikit-learn, and parallel vectors
            giving the original MAPQ predictions and whether or not the
            alignments are correct. """
        if remove_labels is None:
            remove_labels = set()
        if shortname == 'c':
            # extract relevant paired-end features from training data
            labs = ['best1_1',  # score of best alignment for mate
                    'best1_2',  # score of best alignment for opposite mate
                    'diff_1',  # difference for mate
                    'diff_2']  # difference for opposite mate
            if 'diff_conc' in data:
                labs.append('diff_conc')  # concordant difference
            #labs.append('fraglen_z')  # # stddevs diff for fraglen
            labs.append('fraglen')
            if data['rdlen_1'].nunique() > 1:
                labs.append('rdlen_1')
            mapq_header = 'mapq_1'
        elif shortname == 'd':
            # extract relevant discordant paired-end features
            labs = ['best1_1', 'diff_1']
            if data['rdlen_1'].nunique() > 1:
                labs.append('rdlen_1')
            mapq_header = 'mapq_1'
        else:
            # extract relevant unpaired features
            labs = ['best1', 'diff']
            if data['rdlen'].nunique() > 1:
                labs.append('rdlen')
            mapq_header = 'mapq'
        if include_ztzs:
            for colname in sorted(data.columns.values):
                if colname.startswith('ztz'):
                    labs.append(colname)
        for lab in labs:
            assert not np.isnan(data[lab]).any()
        for to_remove in remove_labels:
            if to_remove in labs:
                labs.remove(to_remove)
        data_mat = data[labs].values
        correct = np.array(map(lambda x: x == 1, data['correct']))
        assert not np.isinf(data_mat).any() and not np.isnan(data_mat).any()
        return data_mat, data[mapq_header], correct, labs

    @staticmethod
    def _downsample(x_train, mapq_orig_train, y_train, sample_fraction):
        """ Return a random subset of the data, MAPQs and labels.  Size of
            subset given by sample_fraction. """
        n_training_samples = x_train.shape[0]
        if sample_fraction < 1.0:
            sample_indexes = random.sample(xrange(n_training_samples), int(round(n_training_samples * sample_fraction)))
            x_train = x_train[sample_indexes, ]
            mapq_orig_train = mapq_orig_train.iloc[sample_indexes]
            y_train = np.array([y_train[i] for i in sample_indexes])
        return x_train, mapq_orig_train, y_train

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
    def _crossval_fit(mf_gen, x_train, y_train, dataset_shortname, use_oob=True):
        """ Use cross validation to pick the best model from a
            collection of possible models (model_family) """
        mf = mf_gen()
        scores = []

        def _oob_score(pred_):
            pred_.fit(x_train, y_train)
            return pred_.oob_score_

        def _crossval_score(pred_):
            scores_cv = cross_validation.cross_val_score(pred_, x_train, y_train)
            return float(np.mean(scores_cv))

        while True:
            params, pred = mf.next_predictor()
            if pred is None:
                break
            score = _oob_score(pred) if use_oob else _crossval_score(pred)
            scores.append(score)
            best = mf.set_score(score)
            logging.debug("%s, %s=%0.3f, %s%s" % (dataset_shortname, 'oob' if use_oob else 'score', score, str(params), ' *' if best else ''))
        best_params, best_pred = mf.best_predictor()
        logging.debug("BEST: %s, avg=%0.3f, %s" % (dataset_shortname, max(scores), str(best_params)))
        assert best_pred is not None
        return best_pred, best_params, scores

    datasets = zip('dbcu', ['Discordant', 'Bad-end', 'Concordant', 'Unpaired'], [True, False, True, False])

    def _fit(self, dfs, logger=logging.info, sample_fraction=1.0, sample_fractions=None, remove_labels=None,
             include_ztzs=True):

        if sample_fractions is None:
            sample_fractions = {}

        for ds, ds_long, paired in self.datasets:
            if ds not in dfs:
                continue  # empty
            train = pandas.concat([x for x in dfs.dataset_iter(ds)])
            if train.shape[0] == 0:
                continue  # empty

            logger('Fitting %s training data with random seed %d' % (ds_long, self.random_seed))
            # seed pseudo-random generators
            random.seed(self.random_seed)
            np.random.seed(self.random_seed)
            # extract features, convert to matrix
            x_train, mapq_orig_train, y_train, col_names = self._df_to_mat(train, ds, remove_labels=remove_labels,
                                                                           include_ztzs=include_ztzs)
            self.col_names[ds] = col_names
            # optionally downsample
            frac = sample_fractions[ds] if ds in sample_fractions else sample_fraction
            if frac < 1.0:
                logger('Sampling %0.2f%% of %d rows of %s training data' % (100.0 * frac, train.shape[0], ds_long))
                x_train, mapq_orig_train, y_train = \
                    self._downsample(x_train, mapq_orig_train, y_train, frac)
                logger('  Now has %d rows' % x_train.shape[0])

            # use cross-validation to pick a model
            self.trained_models[ds], params, scores = \
                self._crossval_fit(self.model_gen, x_train, y_train, ds)

            self.trained_params = params
            self.crossval_avg[ds] = max(scores)

            # fit all the training data with the model
            self.trained_models[ds].fit(x_train, y_train)

    def predict(self, dfs,
                keep_names=False, keep_data=False, keep_per_category=False, verbose=False,
                logger=logging.info, remove_labels=None):

        pred_overall = MapqPredictions()
        pred_per_category = {}

        for ds, ds_long, paired in self.datasets:
            if ds not in dfs:
                continue
            names = None
            names_colname = 'name_1' if paired else 'name'
            if keep_per_category:
                pred_per_category[ds] = MapqPredictions()
            nchunk = 0
            for test_chunk in dfs.dataset_iter(ds):
                nchunk += 1
                logger('  Making predictions for %s chunk %d, %d rows' % (ds_long, nchunk, test_chunk.shape[0]))
                x_test, mapq_orig_test, y_test, col_names = self._df_to_mat(test_chunk, ds, remove_labels=remove_labels)
                y_test = np.array(map(lambda c: c == 1, test_chunk.correct))
                pcor = self.trained_models[ds].predict(x_test)  # make predictions
                pcor = np.array(self.postprocess_predictions(pcor, ds_long))
                if keep_names and names_colname in test_chunk:
                    names = test_chunk[names_colname]
                data = x_test.tolist() if keep_data else None
                for prd in [pred_overall, pred_per_category[ds]] if keep_per_category else [pred_overall]:
                    prd.add_pcors(pcor, mapq_orig_test, ds, names=names, data=data, correct=y_test)

        logger('Finalizing results for overall test data (%d alignments)' % len(pred_overall.pcor))
        pred_overall.finalize(verbose=verbose)
        if len(pred_per_category) > 1:
            for ds, pred in pred_per_category.iteritems():
                logger('Finalizing results for "%s" test data (%d alignments)' % (ds, len(pred.pcor)))
                pred.finalize(verbose=verbose)
        logger('Done')

        if keep_per_category:
            return pred_overall, pred_per_category
        else:
            return pred_overall

    def predict_unpaired(self, al):
        """ Given unpaided alignment, use fit model to predict a MAPQ for it """
        # TODO: need to normalize best, secbest the same way we normally do
        # Many aspects:
        # 1. Need to parse and extract raw fields from SAM
        # 2. Need to get normalizers from the training data (should we copy
        #    the normalizers to the MapqFit object so we don't have to keep
        #    anything besides the picked fit?)
        # 3.
        assert 'u' in self.trained_models
        secbest = al.secondBestScore
        if hasattr(al, 'thirdBestScore'):
            secbest = max(secbest, al.thirdBestScore)
        x_test = [al.bestScore, secbest]
        pcor = self.trained_models['u'].predict(x_test)  # make predictions

    def rewrite_sam(self, fn_in, fn_out):
        """ Read SAM input line by line.  Output new lines
        """
        with open(fn_in) as ifh:
            with open(fn_out) as ofh:
                for ln in ifh:
                    if ln[0] != '@':
                        toks = ln.split()

                    ofh.write(ln)

    def __init__(self, dfs, model_gen, random_seed=628599,
                 logger=logging.info, sample_fraction=1.0, sample_fractions=None,
                 remove_labels=None, include_ztzs=True):

        assert random_seed >= 0
        self.model_gen = model_gen
        self.random_seed = random_seed
        self.trained_models = {}
        self.crossval_avg = {}
        self.crossval_std = {}
        self.col_names = {}
        self.trained_params = None
        self._fit(dfs, logger=logger, sample_fraction=sample_fraction, sample_fractions=sample_fractions,
                  remove_labels=remove_labels, include_ztzs=include_ztzs)


class ModelFamily(object):
    """ Encapsulates a model family and a simple interface for naively
        searching the space of hyperparameters. """

    def __init__(self, new_predictor, params, round_to, min_separation, start_in_middle=True):
        """
        new_predictor: function that takes set of params and returns
                       new predictor with those params
        params: list of lists
        """
        self.new_predictor = new_predictor  # returns new predictor given parameters
        self.params = params  # space of possible paramter choices
        self.last_params = None  # remember last set of params used for predictor
        self.round_to = round_to
        self.min_separation = min_separation  # have to improve by at least this much to declare new best
        self.best = float('-inf')
        self.best_rounded = float('-inf')
        self.best_translated_params = None
        if start_in_middle:
            center = tuple([int(round(len(x) / 2)) for x in params])
        else:
            center = tuple([0] * len(params))
        self.workset = {center}
        self.added_to_workset = copy.copy(self.workset)
        self._add_neighbors_to_workset(center)

    def _add_neighbors_to_workset(self, center):
        for i in range(len(self.params)):
            if center[i] > 0:
                neighbor = list(center[:])
                neighbor[i] -= 1  # add next-lowest neighbor
                neighbor = tuple(neighbor)
                if neighbor not in self.added_to_workset:
                    self.added_to_workset.add(neighbor)
                    self.workset.add(neighbor)
            if center[i] < len(self.params[i])-1:
                neighbor = list(center[:])
                neighbor[i] += 1  # add next-highest neighbor
                neighbor = tuple(neighbor)
                if neighbor not in self.added_to_workset:
                    self.added_to_workset.add(neighbor)
                    self.workset.add(neighbor)

    def _idxs_to_params(self, idxs):
        return [self.params[i][j] for i, j in enumerate(idxs)]

    def next_predictor(self):
        if len(self.workset) > 0:
            self.last_params = self.workset.pop()
            translated_params = self._idxs_to_params(self.last_params)
            return translated_params, self.new_predictor(translated_params)
        self.last_params = None
        return None, None

    def set_score(self, score):
        assert self.last_params is not None
        assert self.last_params in self.added_to_workset
        score_rounded = int(score / self.round_to) * self.round_to
        if score_rounded < self.best_rounded:
            return False
        if self.best == float('-inf') or score > self.best * (1.0 + self.min_separation):
            self.best, self.best_rounded = score, score_rounded
            self.best_translated_params = self._idxs_to_params(self.last_params)
            self._add_neighbors_to_workset(self.last_params)
            return True
        return False

    def best_predictor(self):
        return self.best_translated_params, self.new_predictor(self.best_translated_params)


def random_forest_models(random_seed=33, round_to=1e-5, min_separation=0.05):
    # These perform OK but not as well as the extremely random trees
    def _gen(params):
        return RandomForestRegressor(n_estimators=params[0], max_depth=params[1],
                                     random_state=random_seed,
                                     max_features=2,
                                     oob_score=True, bootstrap=True)
    return lambda: ModelFamily(_gen, [range(5, 105, 10), range(2, 10)],
                               round_to, min_separation=min_separation)


def extra_trees_models(random_seed=33, round_to=1e-5, min_separation=0.05):
    # These perform quite well
    def _gen(params):
        return ExtraTreesRegressor(n_estimators=params[0], max_depth=params[1],
                                   random_state=random_seed,
                                   max_features=params[2],
                                   oob_score=True, bootstrap=True)
    return lambda: ModelFamily(_gen, [range(5, 85, 5), range(3, 16, 2), [0.25, 0.5, 1.0]],
                               round_to, min_separation=min_separation)


def model_family(args):
    if args['model_family'] == 'RandomForest':
        return random_forest_models(args['seed'], args['optimization_tolerance'])
    elif args['model_family'] == 'ExtraTrees':
        return extra_trees_models(args['seed'], args['optimization_tolerance'])
    else:
        raise RuntimeError('Bad value for --model-family: "%s"' % args['model_family'])


def mkdir_quiet(dr):
    # Create output directory if needed
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def output_top_incorrect(pred, odir, args):
    mkdir_quiet(odir)
    if args['write_top_incorrect'] or args['write_all']:
        df = pred.summarize_incorrect(1000)
        df.to_csv(os.path.join(odir, 'top1000_incorrect.tsv'), sep='\t', index=False)


def make_plots(pred, odir, args, prefix=''):
    mkdir_quiet(odir)
    fmat = args['plot_format']
    if args['plot_cum_incorrect'] or args['plot_all']:
        if len(pred.pcor) > 1000000:
            logging.warning(prefix + 'SKIPPING cumulative-incorrect plots because there were >1M predictions')
        else:
            logging.info(prefix + 'Making cumulative-incorrect plots')
            assert pred.correct is not None
            plot_drop_rate(pred.pcor, pred.correct, pcor2=pred.pcor_orig, log2ize=False).savefig(
                os.path.join(odir, 'drop_rate.' + fmat))
            plt.close()
            assert len(plt.get_fignums()) == 0
            plot_drop_rate(pred.pcor, pred.correct, pcor2=pred.pcor_orig, log2ize=True).savefig(
                os.path.join(odir, 'drop_rate_log2.' + fmat))
            plt.close()
            assert len(plt.get_fignums()) == 0
            plot_drop_rate_difference(pred.pcor, pred.pcor_orig, pred.correct, log2ize=False).savefig(
                os.path.join(odir, 'drop_rate_diff.' + fmat))
            plt.close()
            assert len(plt.get_fignums()) == 0
            plot_drop_rate_difference(pred.pcor, pred.pcor_orig, pred.correct, log2ize=True).savefig(
                os.path.join(odir, 'drop_rate_diff_log2.' + fmat))
            plt.close()
            assert len(plt.get_fignums()) == 0
    if args['plot_mapq_buckets'] or args['plot_all']:
        logging.info(prefix + 'Making MAPQ bucket plots')
        bucket_error_plot([pred.mapq, pred.mapq_orig], ['Predicted', 'Original'], ['b', 'g'], pred.correct).savefig(
            os.path.join(odir, 'mapq_buckets.' + fmat))
        plt.close()
        assert len(plt.get_fignums()) == 0


def go(args):
    remove_labels = None
    odir = args['output_directory']
    mkdir_quiet(odir)

    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG)

    if args['write_logs'] or args['write_all']:
        fn = os.path.join(odir, 'pred_logs.txt')
        fh = logging.FileHandler(fn)
        fh.setLevel(logging.DEBUG)
        logging.getLogger('').addHandler(fh)

    logging.info('Instantiating model family')
    fam = model_family(args)

    logging.info('Reading datasets')
    input_dir = args['input_directory']
    use_normalizers = not args['no_normalization']
    df_training = AlignmentTableReader(os.path.join(input_dir, 'training'),
                                       use_normalizers=use_normalizers, learn_normalizers=use_normalizers)
    df_test = AlignmentTableReader(os.path.join(input_dir, 'test'),
                                   use_normalizers=use_normalizers, normalizers=df_training.normalizers)

    if args['subsampling_series'] is not None:
        logging.info('Doing subsampling series')
        ss_odir = os.path.join(odir, 'subsampled')
        fractions = map(float, args['subsampling_series'].split(','))
        perf_dicts = {'test': [{'fraction': [],
                               'rank_err_diff_round_pct': [],
                               'auc_diff_round_pct': [],
                               'mse_diff_round_pct': [],
                               'mapq_avg': [],
                               'mapq_std': [],
                               'params': []} for _ in xrange(args['subsampling_replicates'])]}
        perf_dicts['training'] = copy.deepcopy(perf_dicts['test'])
        for fraction in fractions:
            logging.info('  Fraction=%0.3f' % fraction)
            for repl in xrange(1, args['subsampling_replicates']+1):
                my_seed = hash(str(args['seed'] + repl)) % 4294967296
                assert my_seed >= 0
                gc.collect()
                logging.info('    Replicate=%d' % repl)
                my_odir = os.path.join(ss_odir, '%0.3f' % fraction, str(repl))
                mkdir_quiet(my_odir)

                #
                # Model-related outputs
                #

                my_fit_fn = os.path.join(my_odir, 'fit.pkl')
                if os.path.exists(my_fit_fn) and not args['overwrite_fit']:
                    # Read model fit from a pickled file
                    logging.info('      Loading predictions from file')
                    with open(my_fit_fn, 'rb') as fh:
                        ss_fit = cPickle.load(fh)
                else:
                    # Fit model
                    logging.info('      Fitting')
                    ss_fit = MapqFit(df_training, fam, random_seed=my_seed, logger=logging.info,
                                     sample_fraction=fraction, remove_labels=remove_labels,
                                     include_ztzs=not args['ignore_ztzs'])
                    if args['serialize_fit']:
                        logging.info('      Serializing fit object')
                        with open(my_fit_fn, 'wb') as ofh:
                            cPickle.dump(ss_fit, ofh, 2)

                if args['write_feature_importances'] or args['write_all']:
                    logging.info('      Writing feature importances')
                    for ds, model in ss_fit.trained_models.iteritems():
                        my_fi_fn = os.path.join(my_odir, '%s_feature_importances.tsv' % ds)
                        with open(my_fi_fn, 'w') as fh:
                            importances = model.feature_importances_
                            ranks = np.argsort(importances)[::-1]
                            inv_ranks = [0] * len(ranks)
                            for i, r in enumerate(ranks):
                                inv_ranks[r] = i
                            i = 0
                            fh.write('feature\timportance\trank\n')
                            for im, r in zip(importances, inv_ranks):
                                fh.write('%s\t%0.4f\t%d\n' % (ss_fit.col_names[ds][i], im, r))
                                i += 1

                if args['write_oob_scores'] or args['write_all']:
                    logging.info('      Writing out-of-bag scores')
                    for ds, model in ss_fit.trained_models.iteritems():
                        my_fi_fn = os.path.join(my_odir, '%s_oob_score.txt' % ds)
                        with open(my_fi_fn, 'w') as fh:
                            fh.write('%0.5f\n' % model.oob_score_)

                if args['write_params'] or args['write_all']:
                    logging.info('      Writing parameters selected by cross-validation')
                    for ds, model in ss_fit.trained_models.iteritems():
                        my_fi_fn = os.path.join(my_odir, '%s_oob_score.txt' % ds)
                        with open(my_fi_fn, 'w') as fh:
                            fh.write('%0.5f\n' % model.oob_score_)

                #
                # Prediction-related outputs
                #

                # TODO: do this for training as well as test data
                for df, name in [(df_test, 'test'), (df_training, 'training')]:
                    pred_odir = os.path.join(my_odir, name)
                    mkdir_quiet(pred_odir)
                    logging.info('      Making %s predictions' % name)
                    pred_overall, _ = ss_fit.predict(
                        df, verbose=args['verbose'], keep_names=True, keep_data=True,
                        keep_per_category=True, remove_labels=remove_labels)
                    logging.info('        Outputting top incorrect alignments')
                    output_top_incorrect(pred_overall, pred_odir, args)
                    logging.info('        Making plots')
                    make_plots(pred_overall, pred_odir, args, prefix='        ')

                    perf_dicts[name][repl-1]['fraction'].append(fraction)
                    perf_dicts[name][repl-1]['rank_err_diff_round_pct'].append(pred_overall.rank_err_diff_round_pct)
                    perf_dicts[name][repl-1]['auc_diff_round_pct'].append(pred_overall.auc_diff_round_pct)
                    perf_dicts[name][repl-1]['mse_diff_round_pct'].append(pred_overall.mse_diff_round_pct)
                    perf_dicts[name][repl-1]['mapq_avg'].append(pred_overall.mapq_avg)
                    perf_dicts[name][repl-1]['mapq_std'].append(pred_overall.mapq_std)
                    perf_dicts[name][repl-1]['params'].append(str(ss_fit.trained_params))

                    # TODO: print columns giving SSE error for each distinct MAPQ
                    if args['write_roc_table'] or args['write_all']:
                        logging.info('      Writing ROC table')
                        my_roc_fn = os.path.join(pred_odir, 'roc_table.tsv')
                        df = roc_table(pred_overall.pcor, pred_overall.correct, rounded=True, mapqize=True)
                        df.to_csv(my_roc_fn, sep='\t', index=False)
                        my_roc_orig_fn = os.path.join(pred_odir, 'roc_table_orig.tsv')
                        df_orig = roc_table(pred_overall.pcor_orig, pred_overall.correct, rounded=True, mapqize=True)
                        df_orig.to_csv(my_roc_orig_fn, sep='\t', index=False)

                del ss_fit
                gc.collect()

            gc.collect()

        for name in ['test', 'training']:
            dfs = [pandas.DataFrame.from_dict(perf_dict) for perf_dict in perf_dicts[name]]
            for i, df in enumerate(dfs):
                df.to_csv(os.path.join(odir, 'subsampling_series_%s_%d.tsv' % (name, i+1)), sep='\t', index=False)
            plot_subsampling_series(dfs).savefig(os.path.join(odir, 'subsampling_series_%s.%s' % (name, args['plot_format'])))
            plt.close()
            assert len(plt.get_fignums()) == 0

    # if the fit already exists, use it unless --overwrite-fit is specified
    fit_fn = os.path.join(odir, 'fit.pkl')
    if os.path.exists(fit_fn) and not args['overwrite_fit']:
        logging.info('Loading fit from file')
        with open(fit_fn, 'rb') as fh:
            fit = cPickle.load(fh)
    else:
        logging.info('Fitting and making predictions')
        my_seed = hash(str(args['seed']) + '-1') % 4294967296
        assert my_seed >= 0
        fit = MapqFit(df_training, fam, random_seed=my_seed, logger=logging.info, include_ztzs=not args['ignore_ztzs'])
        if args['serialize_fit']:
            logging.info('Serializing fit object')
            with open(fit_fn, 'wb') as ofh:
                cPickle.dump(fit, ofh, 2)

    do_rewrite = True
    #if do_rewrite:
        #logging.info('Rewriting SAM')
        #fit.rewrite_sam(os.path.join(input_dir, 'input.sam'),
        #                os.path.join(input_dir, 'output.sam'))

    logging.info('Done')


def go_profile(args):
    pr = None
    if args['profile']:
        import cProfile
        import pstats
        import StringIO
        pr = cProfile.Profile()
        pr.enable()
    go(args)
    if args['profile']:
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'tottime'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(30)
        print s.getvalue()


def add_predict_args(parser):
    # Output-related options
    parser.add_argument('--input-directory', metavar='path', type=str, required=True,
                        help='Directory with output from tandem simulator')
    parser.add_argument('--output-directory', metavar='path', type=str, required=True,
                        help='Directory to write summaries and plots to')

    # Model-related options
    parser.add_argument('--model-family', metavar='family', type=str, required=False,
                        default='ExtraTrees', help='{RandomForest | ExtraTrees}')
    parser.add_argument('--optimization-tolerance', metavar='float', type=float, default=1e-3,
                        help='Tolerance when searching for best model parameters')
    parser.add_argument('--subsampling-fraction', metavar='float', type=float, default=1.0,
                        help='Subsample the training down to this fraction before fitting model')
    parser.add_argument('--no-normalization', action='store_const', const=True, default=False,
                        help='Don\'t normalize score differences with respect to minimum and maximum')
    parser.add_argument('--ignore-ztzs', action='store_const', const=True, default=False,
                        help='Don\'t include features specified in the ZT:Z extra flag')

    parser.add_argument('--overwrite-fit', action='store_const', const=True, default=False,
                        help='Re-fit the model even if a fit is already present in --output-directory')
    parser.add_argument('--serialize-fit', action='store_const', const=True, default=False,
                        help='Write fit model to a pickle file')
    parser.add_argument('--compression-effort', metavar='int', type=int, default=1,
                        help='How hard to try to compress the model when writing to .pkl file')

    # What to generate
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

    parser.add_argument('--write-top-incorrect', action='store_const', const=True, default=False,
                        help='Write information about the top 1000 misclassified alignments')
    parser.add_argument('--write-logs', action='store_const', const=True, default=False,
                        help='Write verbose prediction log to pred_log.txt in output directory')
    parser.add_argument('--write-feature-importances', action='store_const', const=True, default=False,
                        help='Write importance of each feature according to model to (cat)_feature_importances.tsv')
    parser.add_argument('--write-oob-scores', action='store_const', const=True, default=False,
                        help='Write out-of-bag scores to (cat)_oob_score.txt')
    parser.add_argument('--write-roc-table', action='store_const', const=True, default=False,
                        help='Write table with correct/incorrect stratified by MAPQ to roc_table.tsv in '
                             'output directory')
    parser.add_argument('--write-params', action='store_const', const=True, default=False,
                        help='Write hyper-parameters chosen with cross validation')
    parser.add_argument('--write-all', action='store_const', const=True, default=False,
                        help='Like specifying all the --write-* parameters')

    parser.add_argument('--rewrite-sam', action='store_const', const=True, default=False,
                        help='Like specifying all option beginning with --plot')


if __name__ == "__main__":
    import argparse

    _parser = argparse.ArgumentParser(description='Fit model, make predictions.')

    add_predict_args(_parser)

    # Other options
    _parser.add_argument('--seed', metavar='int', type=int, default=6277, help='Pseudo-random seed')
    _parser.add_argument('--test', action='store_const', const=True, default=False, help='Run unit tests')
    _parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Output profiling data')

    if '--version' in sys.argv:
        print 'Qsim predictor, version ' + VERSION
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

    myargs = _parser.parse_args()

    go_profile(vars(myargs))
