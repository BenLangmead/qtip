import numpy as np
import pandas
import logging
from itertools import repeat
from mapq import pcor_to_mapq, mapq_to_pcor, round_pcor
from collections import defaultdict
from roc import Roc
try:
    from itertools import izip
except ImportError:
    izip = zip

__author__ = 'langmead'


class MapqPredictions:
    """ Encapsulates mapq predictions for a dataset.  Sometimes the data has
        associated correctness information, in which case this class also
        encapsulates performance results. """

    def __init__(self, temp_man, name, calc_summaries=True, prediction_mem_limit=10000000):
        self.name = name
        self.temp_man = temp_man
        self.temp_group_name = 'MAPQ predictions %s' % name

        self.pred_fns = []
        self.pred_fn_prefix = 'predictions_' + name
        self.temp_next_fn = '_'.join([self.pred_fn_prefix, str(len(self.pred_fns))])
        self.pred_fns.append(temp_man.get_file(self.temp_next_fn, group=self.temp_group_name))
        self.pred_fh = open(self.pred_fns[-1], 'wb')
        self.last_id = None

        self.pcor = []
        self.mapq = []
        self.category = []
        self.ids = []
        self.mapq_precision = 3
        self.tally = defaultdict(lambda: [0, 0])
        self.tally_orig = defaultdict(lambda: [0, 0])
        self.tally_rounded = defaultdict(lambda: [0, 0])
        self.roc = None
        self.roc_orig = None
        self.roc_rounded = None
        self.pcor_orig = None
        self.mapq_orig = None
        self.correct = None
        self.data = None
        self.predictions_loaded = False

        self.ordered_by = "?"
        self.calc_summaries = calc_summaries
        self.prediction_mem_limit = prediction_mem_limit
        self.npredictions = 0
        self.has_correct = False
        self.auc_diff_pct = None
        self.auc_diff_round_pct = None
        self.mse_diff_pct = None
        self.mse_diff_round_pct = None

    def add(self, pcor, ids, category, mapq_orig=None, data=None, correct=None):
        """ Add a new batch of predictions. """
        mapq_orig_iter = iter(mapq_orig) if mapq_orig is not None else repeat([None])
        data_iter = iter(data) if data is not None else repeat([None])
        correct_iter = iter(correct) if correct is not None else repeat([None])
        self.has_correct = correct is not None
        if self.last_id is not None and ids[0] < self.last_id:
            self.pred_fh.close()
            self.temp_next_fn = '_'.join([self.pred_fn_prefix, str(len(self.pred_fns))])
            self.pred_fns.append(self.temp_man.get_file(self.temp_next_fn, group=self.temp_group_name))
            self.pred_fh = open(self.pred_fns[-1], 'wb')
        for rec in izip(pcor, ids, repeat([category]), mapq_orig_iter, data_iter, correct_iter):
            self.last_id = rec[1]
            self.pred_fh.write(','.join(map(str, rec)) + '\n')
            self.npredictions += 1
            if rec[5] is not None:
                mapq = pcor_to_mapq(float(rec[0]))
                mapq = round(mapq, self.mapq_precision)
                mapq_orig = float(rec[3])
                self.tally[mapq][0 if rec[5] else 1] += 1
                self.tally_orig[mapq_orig][0 if rec[5] else 1] += 1
                self.tally_rounded[round(mapq)][0 if rec[5] else 1] += 1

    def _reset_mem_predictions(self):
        """ Erase in-memory copies of predictions """
        self.pcor = []
        self.mapq = []
        self.category = []
        self.ids = []
        self.pcor_orig = None
        self.mapq_orig = None
        self.correct = None
        self.data = None

    def _load_predictions(self):
        """ Load all the predictions added with the 'add' member function into
            the appropriate fields of this object.  This should be called only
            with a reasonable number of predictions. """
        if self.npredictions > self.prediction_mem_limit:
            raise RuntimeError('Request to load %d predictions into memory exceeds limit (%d)' %
                               (self.npredictions, self.prediction_mem_limit))
        self._reset_mem_predictions()
        for pred_fn in self.pred_fns:
            with open(pred_fn, 'rb') as ifh:
                for rec_ln in ifh:
                    pc, ident, ct, mq_orig, data, correct = rec_ln.rstrip().split(',')
                    mq_orig = int(mq_orig)
                    correct = int(correct)
                    pc = float(pc)
                    self.pcor.append(pc)
                    self.mapq.append(pcor_to_mapq(pc))
                    self.category.append(ct)
                    self.ids.append(ident)
                    if mq_orig != 'None':
                        if self.mapq_orig is None:
                            self.mapq_orig = []
                        if self.pcor_orig is None:
                            self.pcor_orig = []
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
        return self.calc_summaries and self.has_correct

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
            summ_dict['data'] = list(map(str, [self.data[x] for x in incor_idx]))
        summ_dict['correct'] = [self.correct[x] for x in incor_idx]
        return pandas.DataFrame.from_dict(summ_dict)

    def write_rocs(self, roc_prefix, cid_prefix, csed_prefix):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        assert self.correct is not None
        self.roc.tab.to_csv(roc_prefix + '.csv', sep=',', index=False)
        self.roc_rounded.tab.to_csv(roc_prefix + '_round.csv', sep=',', index=False)
        self.roc_orig.tab.to_csv(roc_prefix + '_orig.csv', sep=',', index=False)
        Roc.write_cum_incorrect_diff(self.roc, self.roc_orig, cid_prefix + '.csv')
        Roc.write_cum_incorrect_diff(self.roc_rounded, self.roc_orig, cid_prefix + '_round.csv')
        Roc.write_cum_squared_error(self.roc, self.roc_orig, csed_prefix + '.csv')
        Roc.write_cum_squared_error(self.roc_rounded, self.roc_orig, csed_prefix + '_round.csv')

    def write_summary_measures(self, fn):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        assert self.correct is not None
        auc_stats = [self.auc_diff_pct, self.auc_diff_round_pct]
        mse_stats = [self.mse_diff_pct, self.mse_diff_round_pct]
        with open(fn, 'w') as fh:
            fh.write(','.join([self.name + '_auc_diff_pct',
                               self.name + '_auc_diff_pct_round',
                               self.name + '_mse_diff_pct',
                               self.name + '_mse_diff_pct_round']) + '\n')
            fh.write(','.join(map(str, auc_stats + mse_stats)) + '\n')

    def write_top_incorrect(self, fn, n=50):
        """ Write a ROC table with # correct/# incorrect stratified by
            predicted MAPQ. """
        self.summarize_incorrect(n=n).to_csv(fn, sep=',', index=False)

    def write_predictions(self, fn):
        """ Write all predictions, in order by the line number of the original
            alignment in the input SAM, to the provided filename. """
        if len(self.pred_fns) == 1:
            last_id = -1
            with open(self.pred_fns[0], 'rb') as ifh:  # no merge needed
                with open(fn, 'w') as ofh:
                    for rec_ln in ifh:
                        pcor, ident, _, _, _, _ = rec_ln.rstrip().split(',')
                        int_ident = int(ident)
                        assert int_ident > last_id, "%d,%d:%s" % (int_ident, last_id, str(self.pred_fns))
                        last_id = int_ident
                        ofh.write('%s,%0.3f\n' % (ident, pcor_to_mapq(float(pcor))))
        else:
            # have to merge!  more complex
            pred_fhs = list(map(lambda x: open(x, 'rb'), self.pred_fns))
            recs = [None] * len(pred_fhs)
            done = [False] * len(pred_fhs)
            nmerged = 0
            last_min_rec = -1
            with open(fn, 'w') as ofh:
                while True:
                    min_rec = (None, float('inf'))
                    min_i = -1
                    for i, pred_fh in enumerate(pred_fhs):
                        if recs[i] is None and not done[i]:
                            ln = pred_fh.readline()
                            if len(ln) == 0:
                                recs[i] = None
                                done[i] = True
                            else:
                                pcor, ident, _, _, _, _ = ln.rstrip().split(',')
                                recs[i] = (pcor, int(ident))
                        if recs[i] is not None and recs[i][1] < min_rec[1]:
                            min_rec, min_i = recs[i], i
                    if min_i == -1:
                        assert all(done)
                        break
                    nmerged += 1
                    assert min_rec[1] > last_min_rec, "%d,%d:%s" % (min_rec[1], last_min_rec, str(self.pred_fns))
                    last_min_rec = min_rec[1]
                    ofh.write('%d,%0.3f\n' % (min_rec[1], pcor_to_mapq(float(min_rec[0]))))
                    recs[min_i] = None
            assert nmerged == self.npredictions, (nmerged, self.npredictions)
            for fh in pred_fhs:
                fh.close()

    def _reorder_by(self, ls):
        """ Reordering helper function """
        if not self.predictions_loaded:
            raise RuntimeError('_reorder_by called without predictions loaded into memory')
        ordr = np.argsort(ls)
        self.pcor = [self.pcor[i] for i in ordr]
        self.mapq = [self.mapq[i] for i in ordr]
        self.ids = [self.ids[i] for i in ordr]
        self.pcor_orig = [self.pcor_orig[i] for i in ordr]
        self.mapq_orig = [self.mapq_orig[i] for i in ordr]
        self.category = [self.category[i] for i in ordr]
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
        """ Close prediction file handle.  If we have the information and flags
            needed for accuracy assessment, then do that too. """

        self.pred_fh.close()
        log.info('  %d records written to %d files' % (self.npredictions, len(self.pred_fns)))

        # calculate error measures and other measures
        if self.can_assess():

            self.roc = Roc(self.tally)
            self.roc_orig = Roc(self.tally_orig)
            self.roc_rounded = Roc(self.tally_rounded)

            log.info('  Correctness information is present; loading predictions into memory')
            self._load_predictions()

            log.info('  Reordering')
            self.order_by_pcor(log=log)

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
