__author__ = 'langmead'

import math
import numpy as np
from mapq import pcor_to_mapq, pcor_to_mapq_np
from metrics import cum_squared_error, drop_rate_cum_sum, tally_cor_per, mseor
import matplotlib.pyplot as plt


def log2ize_p_minus(p, mx=10.0):
    """ Given a probability p, rescale it so that ps are mapped
        to values in [0, mx] and the right-hand end is magnified. """
    return abs(math.log(1 - (1 - (2 ** -mx)) * p, 2))


def unlog2ize_p_minus(p, mx=10.0):
    """ Inverse of log2ize_p_minus """
    return abs((2 ** -p) - 1)/(1 - (2 ** -mx))


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
    axes.fill_between(xs, ys, where=ys >= 0, interpolate=True, color='red')
    axes.fill_between(xs, ys, where=ys <= 0, interpolate=True, color='green')
    axes.grid(True)
    return fig


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
        partition_starts = [int(round(i*n/quantiles)) for i in range(quantiles+1)]
        estimated_lists.append([])
        actual_lists.append([])
        for i in range(quantiles):
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
        labs = [str(i+1) for i in range(len(seriess))]
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
    z = np.array(pcor_to_mapq_np(model.predict(xy), ''), dtype=np.float32).reshape((dim_x, dim_y))

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
