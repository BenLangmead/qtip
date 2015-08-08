"""
Given a directory with output from ts.py, predict new MAPQs.
"""

__author__ = 'langmead'

import pandas
import os
import sys
import math
import random
import numpy as np
import logging
import gc
import cPickle
import copy
from itertools import imap
from sklearn import cross_validation
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from collections import defaultdict, Counter
from mapq import pcor_to_mapq, mapq_to_pcor_np, pcor_to_mapq_np
from metrics import cum_squared_error, drop_rate_cum_sum, tally_cor_per, mseor, ranking_error, auc, roc_table


VERSION = '0.2.0'


class AlignmentTableReader(object):

    """ Reads a table of information describing alignments.  These are tables
        output by ts.py.  Tables might describe training or test alignments.

        The trickiest issue here is whether and how to normalize certain
        features, like alignment scores.  The goal would be to make scores
        more comparable across reads of different lengths.  A read of length
        1000 with alignment score 80 should probably not be considered "close
        to" a read of length 100 with alignment score 80.

        If we know the minimum and maximum possible alignment score (which
        depend on the alignment tool's scoring function, which might in turn
        depend on read length), we can standardize to that interval.  We can't
        generally expect to have that information, though.  Even if we had it,
        we might have to pass a different interval for each read, since they
        can vary in length.

        To avoid these issues, we simply add read length as a feature.  This
        has no effect for uniform-length data, since we also remove zero-
        variance features.
    """

    #            short name   suffix
    datasets = [('u',         '_unp.csv'),
                ('d',         '_disc.csv'),
                ('c',         '_conc.csv'),
                ('b',         '_bad_end.csv')]

    def __init__(self, prefix, chunksize=50000):
        self.prefix = prefix
        self.dfs = {}
        self.readers = {}
        for sn, suf in self.datasets:
            fn = self.prefix + suf
            if any(map(os.path.exists, [fn, fn + '.gz', fn + '.bz2'])):
                self.dfs[sn] = self._fn_to_iterator(fn, chunksize=chunksize)

                def _new_iter(_sn):  # a function that returns a new iterator
                    def _inner():
                        return iter(self.dfs[_sn]())
                    return _inner

                self.readers[sn] = _new_iter(sn)

    @staticmethod
    def _fn_to_iterator(fn, chunksize):
        if os.path.exists(fn):
            return pandas.io.parsers.read_csv(fn, quoting=2, chunksize=chunksize)
        elif os.path.exists(fn + '.gz'):
            return pandas.io.parsers.read_csv(fn + '.gz', quoting=2, chunksize=chunksize, compression='gzip')
        elif os.path.exists(fn + '.bz2'):
            return pandas.io.parsers.read_csv(fn + '.bz2', quoting=2, chunksize=chunksize, compression='bz2')
        else:
            raise RuntimeError('No such file: "%s"' % fn)

    @staticmethod
    def _postprocess_data_frame(df, sn):
        """ Changes 'correct' column to use 0/1 and replaces NAs in the score
            difference columns with small values. """

        def _fill_nas(_df, nm, best_nm, secbest_nm):
            _df[nm] = _df[best_nm].astype(float) - _df[secbest_nm]
            _df[nm] = _df[nm].fillna(np.nanmax(_df[nm])).fillna(0)
            assert not any([math.isnan(x) for x in _df[nm]])

        if df.shape[0] == 0:
            return

        # Turn the correct column into 0/1
        if df['correct'].count() == len(df['correct']):
            df['correct'] = df['correct'].map(lambda x: 1 if x == 'T' else 0)

        if sn in 'c':
            _fill_nas(df, 'diff_1', 'best1_1', 'best2_1')
            _fill_nas(df, 'diff_2', 'best1_2', 'best2_2')

            if math.isnan(df['best1conc'].sum()) or math.isnan(df['best2conc'].sum()):
                logging.warning('Difference of concordants not available, so using minimum of mates')
                df['diff_conc'] = df[['diff_1', 'diff_2']].min(axis=1)
            else:
                _fill_nas(df, 'diff_conc', 'best1conc', 'best2conc')

            for sc_lab in ['diff_1', 'best1_1', 'best2_1',
                           'diff_2', 'best1_2', 'best2_2']:
                assert not math.isnan(df[sc_lab])

        elif sn == 'd':
            _fill_nas(df, 'diff_1', 'best1_1', 'best2_1')
            for sc_lab in ['diff_1', 'best1_1', 'best2_1']:
                assert not math.isnan(df[sc_lab])

        else:
            assert sn in 'ub'
            _fill_nas(df, 'diff', 'best1', 'best2')
            for sc_lab in ['diff', 'best1', 'best2']:
                assert not math.isnan(df[sc_lab])

        return df

    def dataset_iter(self, sn):
        assert sn in self.readers
        return imap(lambda x: self._postprocess_data_frame(x, sn), self.readers[sn]())

    def __contains__(self, o):
        return o in self.readers


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
        return [x for x in range(len(self.correct)-1, -1, -1) if not self.correct[x]]

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
            for i in range(len(correct)-1, -1, -1):
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


"""
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
"""


class MapqFit:
    """
    """

    @staticmethod
    def _df_to_mat(data, shortname, include_ztzs=True):
        """ Convert a data frame read with read_dataset into a matrix
            suitable for use with scikit-learn, and parallel vectors
            giving the original MAPQ predictions and whether or not the
            alignments are correct. """
        if shortname == 'c':
            # extract relevant paired-end features from training data
            labs = ['best1_1',  # score of best alignment for mate
                    'best1_2',  # score of best alignment for opposite mate
                    'diff_1',  # difference for mate
                    'diff_2']  # difference for opposite mate
            if 'diff_conc' in data:
                labs.append('diff_conc')  # concordant difference
            labs.append('fraglen')
            data['rdlen_12'] = data['rdlen_1'] = data['rdlen_2']
            if data['rdlen_1'].nunique() > 1:
                labs.append('rdlen_1')
            if data['rdlen_12'].nunique() > 1:
                labs.append('rdlen_12')
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
                if colname.startswith('ztz') and data[colname].nunique() > 1:
                    labs.append(colname)
        for lab in labs:
            assert not np.isnan(data[lab]).any()
        data_mat = data[labs].values
        correct = np.array(map(lambda x: x == 1, data['correct']))
        assert not np.isinf(data_mat).any() and not np.isnan(data_mat).any()
        return data_mat, data[mapq_header], correct, labs

    @staticmethod
    def _subsample(x_train, mapq_orig_train, y_train, sample_fraction):
        """ Return a random subset of the data, MAPQs and labels.  Size of
            subset given by sample_fraction. """
        n_training_samples = x_train.shape[0]
        if sample_fraction < 1.0:
            sample_indexes = random.sample(range(n_training_samples), int(round(n_training_samples * sample_fraction)))
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
            better, much_better = mf.set_score(score)
            symbol = ''
            if much_better:
                symbol = '*'
            elif better:
                symbol = '+'
            logging.debug("%s, %s=%0.3f, %s%s" % (dataset_shortname, 'oob' if use_oob else 'score', score, str(params), symbol))
        best_params, best_pred = mf.best_predictor()
        logging.debug("BEST: %s, avg=%0.3f, %s" % (dataset_shortname, max(scores), str(best_params)))
        assert best_pred is not None
        return best_pred, best_params, max(scores)

    datasets = zip('dbcu', ['Discordant', 'Bad-end', 'Concordant', 'Unpaired'], [True, False, True, False])

    def _fit(self, dfs, logger=logging.info, sample_fraction=1.0, include_ztzs=True):

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
            x_train, mapq_orig_train, y_train, col_names = self._df_to_mat(train, ds, include_ztzs=include_ztzs)
            self.col_names[ds] = col_names
            # optionally subsample
            frac = sample_fraction
            if frac < 1.0:
                logger('Sampling %0.2f%% of %d rows of %s training data' % (100.0 * frac, train.shape[0], ds_long))
                x_train, mapq_orig_train, y_train = \
                    self._subsample(x_train, mapq_orig_train, y_train, frac)
                logger('  Now has %d rows' % x_train.shape[0])

            # use cross-validation to pick a model
            self.trained_models[ds], self.trained_params, self.crossval_avg[ds] = \
                self._crossval_fit(self.model_gen, x_train, y_train, ds)

            # fit all the training data with the model
            self.trained_models[ds].fit(x_train, y_train)

    def predict(self, dfs,
                keep_names=False, keep_data=False, keep_per_category=False, verbose=False,
                logger=logging.info):

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
                x_test, mapq_orig_test, y_test, col_names = self._df_to_mat(test_chunk, ds)
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

    """
    def predict_unpaired(self, al):
        " Given unpaided alignment, use fit model to predict a MAPQ for it "
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
        " Read SAM input line by line.  Output new lines "
        with open(fn_in) as ifh:
            with open(fn_out) as ofh:
                for ln in ifh:
                    if ln[0] != '@':
                        toks = ln.split()

                    ofh.write(ln)
    """

    def __init__(self,
                 dfs,  # dictionary of data frames, one per alignment type
                 model_gen,  # function that takes vector of hyperparameters, returns new model object
                 random_seed=628599,
                 logger=logging.info,
                 sample_fraction=1.0,  # fraction of training data to actually use
                 include_ztzs=True):

        assert random_seed >= 0
        self.model_gen = model_gen
        self.random_seed = random_seed
        self.trained_models = {}
        self.crossval_avg = {}
        self.crossval_std = {}
        self.col_names = {}
        self.trained_params = None
        self._fit(dfs, logger=logger, sample_fraction=sample_fraction, include_ztzs=include_ztzs)


class ModelFamily(object):
    """ Encapsulates a model family and a simple interface for search
        hyperparameter space. """

    def __init__(self, new_predictor, params, round_to, min_separation, start_in_middle=True):
        """
        new_predictor: function that takes set of params and returns
                       new predictor with those params
        params: list of lists
        """
        self.new_predictor = new_predictor  # returns new predictor given parameters
        self.params = params  # space of possible parameter choices
        self.last_params = None  # remember last set of params used for predictor
        self.round_to = round_to
        self.min_separation = min_separation  # have to improve by at least this much to declare new best
        self.best = self.best_base = self.best_rounded = float('-inf')
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
        if self.best == float('-inf') or score > self.best_base * (1.0 + self.min_separation):
            self.best, self.best_base, self.best_rounded = score, score, score_rounded
            self.best_translated_params = self._idxs_to_params(self.last_params)
            self._add_neighbors_to_workset(self.last_params)
            return True, True
        elif score > self.best:
            self.best, self.best_rounded = score, score_rounded
            self.best_translated_params = self._idxs_to_params(self.last_params)
            return True, False
        return False, False

    def best_predictor(self):
        return self.best_translated_params, self.new_predictor(self.best_translated_params)


def random_forest_models(random_seed=33, round_to=1e-5, min_separation=0.01):
    # These perform OK but not as well as the extremely random trees
    def _gen(params):
        return RandomForestRegressor(n_estimators=params[0], max_depth=params[1],
                                     random_state=random_seed,
                                     max_features=2,
                                     oob_score=True, bootstrap=True)
    return lambda: ModelFamily(_gen, [range(5, 105, 10), range(2, 10)],
                               round_to, min_separation=min_separation)


def extra_trees_models(random_seed=33, round_to=1e-5, min_separation=0.002):
    # These perform quite well
    def _gen(params):
        return ExtraTreesRegressor(n_estimators=params[0], max_depth=params[1],
                                   random_state=random_seed,
                                   max_features=0.5,
                                   oob_score=True, bootstrap=True)
    return lambda: ModelFamily(_gen, [range(5, 85, 2), range(3, 16, 1)],
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
    df_training = AlignmentTableReader(os.path.join(input_dir, 'training'))
    df_test = AlignmentTableReader(os.path.join(input_dir, 'test'))

    nrep = args['subsampling_replicates']
    if args['subsampling_series'] is not None:
        logging.info('Doing subsampling series')
        ss_odir = os.path.join(odir, 'subsampled')
        fractions = map(float, args['subsampling_series'].split(','))
        perf_dicts = {'test': [defaultdict(list) for _ in range(nrep)],
                      'training': [defaultdict(list) for _ in range(nrep)]}
        for fraction in fractions:
            logging.info('  Fraction=%0.3f' % fraction)
            for repl in range(1, nrep+1):
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
                                     sample_fraction=fraction, include_ztzs=not args['ignore_ztzs'])
                    if args['serialize_fit']:
                        logging.info('      Serializing fit object')
                        with open(my_fit_fn, 'wb') as ofh:
                            cPickle.dump(ss_fit, ofh, 2)

                feat_import_colnames = []
                feat_import_values = []
                logging.info('      Gathering and writing feature importances')
                write_importances = args['write_feature_importances'] or args['write_all']
                for ds, model in sorted(ss_fit.trained_models.items()):
                    my_fi_fn = os.path.join(my_odir, '%s_feature_importances.tsv' % ds)
                    with open(my_fi_fn if write_importances else os.devnull, 'w') as fh:
                        importances = model.feature_importances_
                        ranks = np.argsort(importances)[::-1]
                        inv_ranks = [0] * len(ranks)
                        for i, r in enumerate(ranks):
                            inv_ranks[r] = i
                        i = 0
                        fh.write('feature\timportance\trank\n')
                        for im, r in zip(importances, inv_ranks):
                            feat_import_colnames.append(ds + '_' + ss_fit.col_names[ds][i])
                            feat_import_values.append(im)
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

                for df, name in [(df_test, 'test'), (df_training, 'training')]:
                    pred_odir = os.path.join(my_odir, name)
                    mkdir_quiet(pred_odir)
                    logging.info('      Making %s predictions' % name)
                    pred_overall, _ = ss_fit.predict(
                        df, verbose=args['verbose'], keep_names=True, keep_data=True, keep_per_category=True)
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
                    for feat, val in zip(feat_import_colnames, feat_import_values):
                        perf_dicts[name][repl-1][feat].append(val)
                    # TODO: tease out ranking error due to conc, disc, bad_end

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
        ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
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
