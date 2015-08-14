import numpy as np
from collections import defaultdict
from mapq import round_pcor_np, pcor_to_mapq


def tally_cor_per(level, cor):
    """ Some predictions from our model are the same; this helper function
        gathers all the correct/incorrect information for each group of equal
        predictions. """
    tally = defaultdict(lambda: [0, 0])
    for p, c in zip(level, cor):
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
    import pandas
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
            err += frac * sum(range(sofar, sofar + ntot))
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
    assert all(pcor[i] >= pcor[i+1] for i in range(len(pcor)-1))
    return np.cumsum((pcor - cor) ** 2)


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
        for i in range(ncor + nincor):
            cumsum += (float(nincor) / (ncor + nincor))
            cumsums.append(cumsum)
    return cumsums
