import numpy as np
import pandas
import os
import logging
from itertools import repeat
from mapq import pcor_to_mapq, mapq_to_pcor
from metrics import mseor, ranking_error, auc, roc_table

__author__ = 'langmead'


class MapqPredictions:
    """ Encapsulates mapq predictions for a dataset.  Sometimes the data has
        associated correctness information, in which case this class also
        encapsulates performance results. """

    def __init__(self, temp_man, name, calc_summaries=True, prediction_mem_limit=10000000):
        self.temp_man = temp_man
        self.temp_group_name = 'MAPQ predictions %s' % name
        self.pred_fn = temp_man.get_file('predictions_' + name, group=self.temp_group_name)
        self.pred_fh = open(self.pred_fn, 'wb')
        self.pcor = []
        self.mapq = []
        self.category = []
        self.pcor_orig = None
        self.mapq_orig = None
        self.correct = None
        self.predictions_loaded = False

        self.ordered_by = "?"
        self.calc_summaries = calc_summaries
        self.prediction_mem_limit = prediction_mem_limit
        self.npredictions = 0
        self.correct_end, self.correct_run = 0, 0
        self.has_correct = None
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

    def add(self, pcor, ids, category, mapq_orig=None, data=None, correct=None):
        """ Add a new batch of predictions """
        mapq_orig_iter = iter(mapq_orig) if mapq_orig is not None else repeat([None])
        data_iter = iter(data) if data is not None else repeat([None])
        correct_iter = iter(correct) if correct is not None else repeat([None])
        self.has_correct = correct is not None
        for rec in zip(pcor, ids, category, mapq_orig_iter, data_iter, correct_iter):
            self.pred_fh.write(','.join(map(str, rec)) + '\n')
            self.npredictions += 1

    def _reset_mem_predictions(self):
        """ Erase in-memory copies of predictions """
        self.pcor = []
        self.mapq = []
        self.category = []
        self.pcor_orig = None
        self.mapq_orig = None
        self.correct = None

    def _load_predictions(self):
        """ Load all the predictions added with the 'add' member function into
            the appropriate fields of this object.  This should be called only
            with a reasonable number of predictions. """
        if self.npredictions > self.prediction_mem_limit:
            raise RuntimeError('Request to load %d predictions into memory exceeds limit (%d)' %
                               (self.npredictions, self.prediction_mem_limit))
        self._reset_mem_predictions()
        with open(self.pred_fn, 'rb') as ifh:
            for rec_ln in ifh:
                pc, _, ct, mq_orig, data, correct = rec_ln.rstrip().split(',')
                mq_orig = int(mq_orig)
                correct = int(correct)
                pc = float(pc)
                self.pcor.append(pc)
                self.mapq.append(pcor_to_mapq(pc))
                self.category.append(ct)
                if mq_orig != 'None':
                    if self.mapq_orig is None:
                        self.mapq_orig = []
                    mq_orig = int(mq_orig)
                    self.mapq_orig.append(mq_orig)
                    self.pcor_orig.append(mapq_to_pcor(mq_orig))
                if data != 'None':
                    if self.data is None:
                        self.data = []
                    self.data.append(data)
                if correct != 'None':
                    if self.correct is None:
                        self.correct = []
                    self.correct.append(correct)
        self.predictions_loaded = True

    def can_assess(self):
        """ Return true iff we have the data and the flags needed to do an
            accuracy assessment. """
        return self.calc_summaries and self.correct is not None and max(self.correct) > -1

    def incorrect_indexes(self):
        """ Return indexes of in correct alignments in order
            from highest to lowest predicted pcor """
        assert self.correct is not None
        assert self.ordered_by == "pcor"
        return [x for x in range(len(self.correct)-1, -1, -1) if not self.correct[x]]

    def summarize_incorrect(self, n=50):
        """ Return a DataFrame summarizing information about """
        assert self.correct is not None
        assert self.ordered_by == "pcor"
        incor_idx = self.incorrect_indexes()[:n]
        summ_dict = dict()
        summ_dict['category'] = [self.category[x] for x in incor_idx]
        summ_dict['mapq'] = [self.mapq[x] for x in incor_idx]
        summ_dict['mapq_orig'] = [self.mapq_orig[x] for x in incor_idx]
        if self.data is not None:
            summ_dict['data'] = map(lambda x: ','.join(map(lambda y: '%0.3f' % y, x)),
                                    [self.data[x] for x in incor_idx])
        summ_dict['correct'] = [self.correct[x] for x in incor_idx]
        return pandas.DataFrame.from_dict(summ_dict)

    def write_rocs(self, fn, fn_orig):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        assert self.correct is not None
        roc_table(self.pcor, self.correct, rounded=True, mapqize=True).to_csv(fn, sep=',', index=False)
        roc_table(self.pcor_orig, self.correct, rounded=True, mapqize=True).to_csv(fn_orig, sep=',', index=False)

    def write_summary_measures(self, fn):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        assert self.correct is not None
        rank_err_stats = [self.rank_err_diff_pct, self.rank_err_diff_round_pct, self.rank_err_raw, self.rank_err_orig]
        auc_stats = [self.auc_diff_pct, self.auc_diff_round_pct]
        mse_stats = [self.mse_diff_pct, self.mse_diff_round_pct]
        with open(fn, 'w') as fh:
            fh.write(','.join(['rank_err_diff_pct', 'rank_err_diff_pct_round', 'rank_err', 'rank_err_orig',
                               'auc_diff_pct', 'auc_diff_pct_round',
                               'mse_diff_pct', 'mse_diff_pct_round',
                               'correct_run']) + '\n')
            fh.write(','.join(map(str, rank_err_stats + auc_stats + mse_stats + [self.correct_run])) + '\n')

    def write_top_incorrect(self, fn, n=50):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        self.summarize_incorrect(n=n).to_csv(fn, sep=',', index=False)

    def write_predictions(self, fn):
        """ Write all predictions, in order by the line number of the original
            alignment in the input SAM, to the provided filename. """
        assert os.path.exists(self.pred_fn)
        with open(self.pred_fn, 'rb') as ifh:
            with open(fn, 'w') as ofh:
                for rec_ln in ifh:
                    pcor, ident, _, _, _, _ = rec_ln.rstrip().split(',')
                    ofh.write('%s,%0.3f\n' % (ident, pcor_to_mapq(float(pcor))))

    def _reorder_by(self, ls):
        """ Reordering helper function """
        if not self.predictions_loaded:
            raise RuntimeError('_reorder_by called without predictions loaded into memory')
        ordr = np.argsort(ls)
        self.pcor = self.pcor[ordr]
        self.mapq = self.mapq[ordr]
        self.ids = self.ids[ordr]
        self.pcor_orig = self.pcor_orig[ordr]
        self.mapq_orig = self.mapq_orig[ordr]
        self.category = [self.category[x] for x in ordr]
        if self.correct is not None:
            self.correct = [self.correct[x] for x in ordr]

    def order_by_ids(self, log=logging):
        """ Reorder in-memory predictions by id of the alignment (low to high) """
        if not self.predictions_loaded:
            raise RuntimeError('order_by_ids called without predictions loaded into memory')
        log.info('  Reordering by read id')
        self._reorder_by(self.ids)
        self.ordered_by = "id"

    def order_by_pcor(self, log=logging):
        """ Reorder in-memory predictions by pcor (high to low) """
        if not self.predictions_loaded:
            raise RuntimeError('order_by_pcor called without predictions loaded into memory')
        log.info('  Reordering by pcor')
        self._reorder_by(self.pcor)
        self.ordered_by = "pcor"

    def purge_temporaries(self):
        self.temp_man.remove_group(self.temp_group_name)

    def finalize(self, log=logging):
        # calculate error measures and other measures
        if self.can_assess():

            log.info('  Correctness information is present; loading predictions into memory')
            self._load_predictions()

            log.info('  Reordering')
            self.order_by_pcor(log=log)
            correct = self.correct

            log.info('  Compiling performance measures')

            # calculate # of highest pcors and max # pcors in a row that
            # correspond to correct alignments
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
            log.info('  Calculating rank error')
            self.rank_err_orig = ranking_error(self.pcor_orig, correct) + 1
            self.rank_err_raw = ranking_error(self.pcor, correct) + 1
            self.rank_err_raw_round = ranking_error(self.pcor, correct, rounded=True) + 1
            self.rank_err_diff = self.rank_err_raw - self.rank_err_orig
            self.rank_err_diff_pct = 100.0 * self.rank_err_diff / self.rank_err_orig
            self.rank_err_diff_round = self.rank_err_raw_round - self.rank_err_orig
            self.rank_err_diff_round_pct = 100.0 * self.rank_err_diff_round / self.rank_err_orig
            self.rank_err = self.rank_err_raw / self.rank_err_orig
            self.rank_err_round = self.rank_err_raw_round / self.rank_err_orig
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.rank_err_diff_pct,
                                                               self.rank_err_diff_round_pct))

            log.info('  Calculating AUC')
            self.auc_orig = auc(self.pcor_orig, correct)
            self.auc_raw = auc(self.pcor, correct)
            self.auc_raw_round = auc(self.pcor, correct, rounded=True)
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
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.auc_diff_pct, self.auc_diff_round_pct))

            log.info('  Calculating MSE')
            self.mse_orig = mseor(self.pcor_orig, correct)
            self.mse_raw = mseor(self.pcor, correct)
            self.mse_raw_round = mseor(self.pcor, correct, rounded=True)
            self.mse_diff = self.mse_raw - self.mse_orig
            self.mse_diff_pct = 100.0 * self.mse_diff / self.mse_orig
            self.mse_diff_round = self.mse_raw_round - self.mse_orig
            self.mse_diff_round_pct = 100.0 * self.mse_diff_round / self.mse_orig
            self.mse = self.mse_raw / self.mse_orig
            self.mse_round = mseor(self.pcor, correct, rounded=True) / self.mse_orig
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.mse_diff_pct, self.mse_diff_round_pct))
