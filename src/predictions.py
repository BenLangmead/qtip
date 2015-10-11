import numpy as np
import pandas
import logging
from collections import Counter
from metrics import mseor, ranking_error, auc, roc_table

__author__ = 'langmead'


class MapqPredictions:
    """ Encapsulates mapq predictions for a dataset.  Sometimes the data has
        associated correctness information, in which case this class also
        encapsulates performance results. """

    def __init__(self):
        # all these lists are parallel
        self.pcor = np.array([])  # predicted pcors
        self.ids = np.array([])  # line ids for original alignment records
        self.mapq_orig = np.array([])  # original mapping qualities
        self.category = []  # categories of alignments
        self.data = None  # data that gave rise to predictions
        self.correct = None  # whether or not alignment is correct
        self.pcor_hist = None
        self.ordered_by = "?"
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

    def add_pcors(self, pcor, ids, mapq_orig, category, data=None, correct=None):
        """ Add a new batch of predictions """
        self.pcor = np.append(self.pcor, pcor)
        self.ids = np.append(self.ids, ids)
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

    def has_correctness(self):
        """ Return true iff we have correctness values for the alignments """
        return self.correct is not None

    def incorrect_indexes(self):
        """ Return indexes of in correct alignments in order
            from highest to lowest predicted pcor """
        assert self.correct is not None
        assert self.ordered_by == "pcor"
        return [x for x in range(len(self.correct)-1, -1, -1) if not self.correct[x]]

    def top_incorrect(self, n=50):
        """ Get incorrect alignments with highest predicted MAPQ """
        assert self.data is not None
        assert self.ordered_by == "pcor"
        return [self.data[x] for x in self.incorrect_indexes()[:n]]

    def summarize_incorrect(self, n=50):
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
        rank_err_stats = [self.rank_err_diff_pct, self.rank_err_diff_round_pct]
        auc_stats = [self.auc_diff_pct, self.auc_diff_round_pct]
        mse_stats = [self.mse_diff_pct, self.mse_diff_round_pct]
        with open(fn, 'w') as fh:
            fh.write(','.join(['rank_err', 'rank_err_round',
                               'auc', 'auc_round',
                               'mse_err', 'mse_err_round']) + '\n')
            fh.write(','.join(map(str, rank_err_stats + auc_stats + mse_stats)) + '\n')

    def write_top_incorrect(self, fn, n=50):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        assert self.correct is not None
        assert self.ordered_by == "pcor", self.ordered_by
        self.summarize_incorrect(n=n).to_csv(fn, sep=',', index=False)

    def write_predictions(self, fn):
        """
        Write all predictions, in order by the line number of the original
        alignment in the input SAM, to the provided filename.
        """
        assert self.ordered_by == "id"
        with open(fn, 'w') as fh:
            for i in range(len(self.mapq)):
                fh.write('%d,%0.3f\n' % (self.ids[i], self.mapq[i]))

    def _reorder_by(self, ls):
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
        log.info('  Reordering by read id')
        self._reorder_by(self.ids)
        self.ordered_by = "id"

    def order_by_pcor(self, log=logging):
        log.info('  Reordering by pcor')
        self._reorder_by(self.pcor)
        self.ordered_by = "pcor"

    def finalize(self, log=logging):
        self.pcor_hist = Counter(self.pcor)

        pcor, mapq_orig = self.pcor, self.mapq_orig
        self.mapq = mapq = np.abs(-10.0 * np.log10(1.0 - pcor))
        self.pcor_orig = pcor_orig = 1.0 - 10.0 ** (-0.1 * mapq_orig)

        # calculate error measures and other measures
        if max(self.correct) > -1:
            log.info('  Correctness information is present; compiling error measures')
            self.order_by_pcor(log=log)
            correct = self.correct

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
            self.rank_err_orig = ranking_error(pcor_orig, correct) + 1
            self.rank_err_raw = ranking_error(pcor, correct) + 1
            self.rank_err_raw_round = ranking_error(pcor, correct, rounded=True) + 1
            self.rank_err_diff = self.rank_err_raw - self.rank_err_orig
            self.rank_err_diff_pct = 100.0 * self.rank_err_diff / self.rank_err_orig
            self.rank_err_diff_round = self.rank_err_raw_round - self.rank_err_orig
            self.rank_err_diff_round_pct = 100.0 * self.rank_err_diff_round / self.rank_err_orig
            self.rank_err = self.rank_err_raw / self.rank_err_orig
            self.rank_err_round = self.rank_err_raw_round / self.rank_err_orig
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.rank_err_diff_pct,
                                                               self.rank_err_diff_round_pct))

            log.info('  Calculating AUC')
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
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.auc_diff_pct, self.auc_diff_round_pct))

            log.info('  Calculating MSE')
            self.mse_orig = mseor(pcor_orig, correct)
            self.mse_raw = mseor(pcor, correct)
            self.mse_raw_round = mseor(pcor, correct, rounded=True)
            self.mse_diff = self.mse_raw - self.mse_orig
            self.mse_diff_pct = 100.0 * self.mse_diff / self.mse_orig
            self.mse_diff_round = self.mse_raw_round - self.mse_orig
            self.mse_diff_round_pct = 100.0 * self.mse_diff_round / self.mse_orig
            self.mse = self.mse_raw / self.mse_orig
            self.mse_round = mseor(pcor, correct, rounded=True) / self.mse_orig
            log.info('    Done: %+0.4f%%, %+0.4f%% rounded' % (self.mse_diff_pct, self.mse_diff_round_pct))

        log.info('  Calculating MAPQ summaries')
        self.mapq_avg, self.mapq_orig_avg = float(np.mean(mapq)), float(np.mean(mapq_orig))
        self.mapq_std, self.mapq_orig_std = float(np.std(mapq)), float(np.std(mapq_orig))
