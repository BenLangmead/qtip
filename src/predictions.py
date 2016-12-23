"""
Copyright 2016, Ben Langmead <langmea@cs.jhu.edu>

MapqPredictions class for storing and analyzing predictions.
"""

import os
import pandas
import logging
from collections import Counter
try:
    from itertools import izip
except ImportError:
    izip = zip

from roc import Roc
from metamat import MetaMat

# qtip imports
__author__ = 'langmead'


class MapqPredictions:
    """ Encapsulates mapq predictions for a dataset.  Sometimes the data has
        associated correctness information, in which case this class also
        encapsulates performance results. """

    def __init__(self, temp_man, name, calc_summaries=True, prediction_mem_limit=10000000):
        self.name = name
        self.temp_man = temp_man
        self.temp_group_name = 'MAPQ predictions %s' % name

        self.pred_fn_prefix = prefix = 'predictions_' + name
        pred_fn = '_'.join([prefix, '0.npy'])
        pred_meta_fn = '_'.join([prefix, '0.meta'])
        self.pred_fns = [temp_man.get_file(pred_fn, group=self.temp_group_name)]
        self.pred_fhs = [open(self.pred_fns[-1], 'wb')]
        self.pred_meta_fns = [temp_man.get_file(pred_meta_fn, group=self.temp_group_name)]
        self.pred_nrow = [0]
        self.pred_fh_new = True
        self.last_id = None

        self.mapq_precision = 3
        self.tally = Counter()
        self.tally_orig = Counter()
        self.tally_rounded = Counter()
        self.roc = None
        self.roc_orig = None
        self.roc_rounded = None
        self.df = None

        self.calc_summaries = calc_summaries
        self.prediction_mem_limit = prediction_mem_limit
        self.npredictions = 0
        self.has_correct = False
        self.auc_diff_pct = None
        self.auc_diff_round_pct = None
        self.mse_diff_pct = None
        self.mse_diff_round_pct = None

    def add(self, recs, first_id, last_id, mapq=None, mapq_orig=None, correct=None):
        """ Add a new batch of predictions. """
        if recs.shape[0] == 0:
            return

        # Open new csv if we got a discontiguous chunk. Requires merging later.
        if self.last_id is not None and first_id < self.last_id:
            pred_fn = '_'.join([self.pred_fn_prefix, str(len(self.pred_fns)) + '.npy'])
            pred_meta_fn = '_'.join([self.pred_fn_prefix, str(len(self.pred_fns)) + '.meta'])
            self.pred_fns.append(self.temp_man.get_file(pred_fn, group=self.temp_group_name))
            self.pred_meta_fns.append(self.temp_man.get_file(pred_meta_fn, group=self.temp_group_name))
            self.pred_fhs.append(open(self.pred_fns[-1], 'wb'))
            self.pred_nrow.append(0)

        # This is performance-critical
        recs.values.tofile(self.pred_fhs[-1], sep='')
        self.npredictions += recs.shape[0]
        self.pred_nrow[-1] += recs.shape[0]
        self.last_id = last_id

        # Update tallies if possible
        self.has_correct = mapq is not None
        if self.has_correct:
            assert mapq_orig is not None
            assert correct is not None
            self.tally.update(zip(mapq.round(decimals=self.mapq_precision), correct))
            self.tally_orig.update(zip(mapq_orig, correct))
            self.tally_rounded.update(zip(mapq.round(decimals=0), correct))

    def _load_predictions(self):
        """ Load all the predictions added with the 'add' member function into
            the appropriate fields of this object.  This should be called only
            with a reasonable number of predictions. """
        if self.npredictions > self.prediction_mem_limit:
            raise RuntimeError('Request to load %d predictions into memory exceeds limit (%d)' %
                               (self.npredictions, self.prediction_mem_limit))
        dfs = []
        for fn in self.pred_fns:
            assert fn.endswith('.npy')
            dfs.append(MetaMat(fn[:-4], chunk_size=-1).next())
        self.df = pandas.concat(dfs)

    def can_assess(self):
        """ Return true iff we have the data and the flags needed to do an
            accuracy assessment. """
        return self.calc_summaries and self.has_correct

    def incorrect_indexes(self):
        """ Return indexes of in correct alignments in order
            from highest to lowest predicted pcor """
        assert self.df.correct.max() > -1
        assert pandas.algos.is_monotonic_float64(self.df.mapq.values, True)
        return self.df.correct[self.df.correct == 0].index[::-1].tolist()

    def summarize_incorrect(self, n=50):
        """ Return a DataFrame summarizing information about """
        assert self.df.correct.max() > -1
        assert pandas.algos.is_monotonic_float64(self.df.mapq.values, True)
        cols = ['category', 'mapq', 'mapq_orig', 'data', 'correct']
        return self.df.loc[self.incorrect_indexes()[:n], cols]

    def write_rocs(self, roc_prefix):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        assert self.has_correct
        self.roc.tab.to_csv(roc_prefix + '.csv', sep=',', index=False, encoding='utf-8')
        self.roc_rounded.tab.to_csv(roc_prefix + '_round.csv', sep=',', index=False, encoding='utf-8')
        self.roc_orig.tab.to_csv(roc_prefix + '_orig.csv', sep=',', index=False, encoding='utf-8')

    def write_summary_measures(self, fn):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        auc_stats = [self.auc_diff_pct, self.auc_diff_round_pct]
        mse_stats = [self.mse_diff_pct, self.mse_diff_round_pct]
        with open(fn, 'wb') as fh:
            fh.write((','.join([self.name + '_auc_diff_pct',
                                self.name + '_auc_diff_pct_round',
                                self.name + '_mse_diff_pct',
                                self.name + '_mse_diff_pct_round']) + '\n').encode('utf-8'))
            fh.write((','.join(map(str, auc_stats + mse_stats)) + '\n').encode('utf-8'))

    def write_top_incorrect(self, fn, n=100):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        self.summarize_incorrect(n=n).to_csv(fn, sep=',', index=False, encoding='utf-8')

    def write_predictions(self, fn):
        """ Write all predictions, in order by the line number of the original
            alignment in the input SAM, to the provided filename. """
        # TODO: qtip-rewrite must now handle the merging
        return

        #sort_cmd = "sort -t',' -m -n -k 1,1 -S 20M"
        #cut_cmd = "cut -d',' -f1,2"
        #if len(self.pred_fns) == 1:
        #    ret = os.system(' '.join([cut_cmd, '<', self.pred_fns[0], '>', fn]))
        #    if ret != 0:
        #        raise RuntimeError('cut command returned %d' % ret)
        #else:
        #    ret = os.system(' '.join([sort_cmd] + self.pred_fns + ['|', cut_cmd, '>', fn]))
        #    if ret != 0:
        #        raise RuntimeError('sort & cut command returned %d' % ret)

    def purge_temporaries(self):
        self.temp_man.remove_group(self.temp_group_name)

    def finalize(self, log=logging):
        """ Close prediction file handle.  If we have the information and flags
            needed for accuracy assessment, then do that too. """

        # Finish writing numpy files
        for fh in self.pred_fhs:
            fh.close()

        # Write metadata for prediction files
        columns = ['ids', 'mapq', 'category', 'mapq_orig', 'correct']
        for fn, n_row in zip(self.pred_meta_fns, self.pred_nrow):
            with open(fn, 'wb') as fh:
                fh.write(b','.join(map(lambda x: x.encode(), columns)))
                fh.write(b',')
                fh.write(str(n_row).encode())

        log.info('  %d records written to %d files' % (self.npredictions, len(self.pred_fns)))

        # calculate error measures and other measures
        if self.can_assess():

            self.roc = Roc(self.tally)
            self.roc_orig = Roc(self.tally_orig)
            self.roc_rounded = Roc(self.tally_rounded)

            log.info('  Correctness information is present; loading predictions into memory')
            self._load_predictions()
            assert self.df is not None
            assert isinstance(self.df, pandas.DataFrame)
            assert 'mapq' in self.df
            assert 'mapq_orig' in self.df
            assert 'correct' in self.df
            self.df.mapq_orig = self.df.mapq_orig.astype(int)
            self.df.mapq_orig = self.df.mapq_orig.astype(int)
            self.df.correct = self.df.correct.astype(int)

            log.info('  Reordering')
            self.df.sort_values('mapq', ascending=False, inplace=True)

            log.info('  Calculating AUC')
            auc_orig = self.roc_orig.area_under_cumulative_incorrect()
            auc_raw = self.roc.area_under_cumulative_incorrect()
            auc_raw_round = self.roc_rounded.area_under_cumulative_incorrect()
            if auc_orig == 0.:
                if auc_raw > auc_orig:
                    self.auc_diff_pct = float('inf')
                else:
                    self.auc_diff_pct = 0.0
                if auc_raw_round > auc_orig > 0.:
                    self.auc_diff_round_pct = float('inf')
                else:
                    self.auc_diff_round_pct = 0
            else:
                self.auc_diff_pct = 100.0 * (auc_raw - auc_orig) / auc_orig
                self.auc_diff_round_pct = 100.0 * (auc_raw_round - auc_orig) / auc_orig
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.auc_diff_pct,
                                                               self.auc_diff_round_pct))

            log.info('  Calculating MSE')
            mse_orig = self.roc_orig.sum_of_squared_error()
            mse_raw = self.roc.sum_of_squared_error()
            mse_raw_round = self.roc_rounded.sum_of_squared_error()
            self.mse_diff_pct = 100.0 * (mse_raw - mse_orig) / mse_orig
            self.mse_diff_round_pct = 100.0 * (mse_raw_round - mse_orig) / mse_orig
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.mse_diff_pct, self.mse_diff_round_pct))
