import os
import pandas
import logging
import csv
from mapq import pcor_to_mapq, mapq_to_pcor
from collections import Counter
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

    def add(self, recs):
        """ Add a new batch of predictions. """
        if recs.shape[0] == 0:
            return
        first_id = recs.ids[0]
        self.has_correct = recs.correct.max() > -1
        if self.last_id is not None and first_id < self.last_id:
            self.temp_next_fn = '_'.join([self.pred_fn_prefix, str(len(self.pred_fns))])
            self.pred_fns.append(self.temp_man.get_file(self.temp_next_fn, group=self.temp_group_name))
            self.pred_fh_new = True
        recs.to_csv(self.pred_fns[-1], header=self.pred_fh_new, sep=',', index=False,
                    quoting=csv.QUOTE_NONE, mode='wt' if self.pred_fh_new else 'at',
                    encoding='utf-8')
        self.pred_fh_new = False
        self.npredictions += len(recs)
        self.last_id = recs.ids.iloc[-1]
        if self.has_correct:
            self.tally.update(zip(recs.mapq.round(decimals=self.mapq_precision), recs.correct))
            self.tally_orig.update(zip(recs.mapq_orig, recs.correct))
            self.tally_rounded.update(zip(recs.mapq.round(decimals=0), recs.correct))

    def _load_predictions(self):
        """ Load all the predictions added with the 'add' member function into
            the appropriate fields of this object.  This should be called only
            with a reasonable number of predictions. """
        if self.npredictions > self.prediction_mem_limit:
            raise RuntimeError('Request to load %d predictions into memory exceeds limit (%d)' %
                               (self.npredictions, self.prediction_mem_limit))
        df = pandas.io.parsers.read_csv(self.pred_fns[0], header=0, quoting=3, encoding='utf-8')
        for pred_fn in self.pred_fns[1:]:
            df = df.append(pandas.io.parsers.read_csv(pred_fn, header=0, quoting=3, encoding='utf-8'), ignore_index=True)
        self.df = df

    def can_assess(self):
        """ Return true iff we have the data and the flags needed to do an
            accuracy assessment. """
        return self.calc_summaries and self.has_correct

    def incorrect_indexes(self):
        """ Return indexes of in correct alignments in order
            from highest to lowest predicted pcor """
        assert self.df.correct.max() > -1
        assert pandas.algos.is_monotonic_float64(self.df.pcor.values, True)
        return self.df.correct[self.df.correct == 0].index[::-1].tolist()

    def summarize_incorrect(self, n=50):
        """ Return a DataFrame summarizing information about """
        assert self.df.correct.max() > -1
        assert pandas.algos.is_monotonic_float64(self.df.pcor.values, True)
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
        if len(self.pred_fns) == 1:
            with open(fn, 'w') as ofh:
                for chunk in pandas.io.parsers.read_csv(
                        self.pred_fns[0], quoting=3, chunksize=100000, encoding='utf-8'):
                    chunk.to_csv(ofh, sep=',', index=False, columns=['ids', 'mapq'],
                                 header=False, encoding='utf-8')
        else:
            if True:
                # 4th column: ids, 5th column: mapqs
                sort_cmd = "sort -t',' -m -n -k 4,4 -S 20M"
                cut_cmd = "cut -d',' -f4,5"
                ret = os.system(' '.join([sort_cmd] + self.pred_fns + ['|', cut_cmd, '>', fn]))
                if ret != 0:
                    raise RuntimeError('sort & cut command returned %d' % ret)
            else:
                # have to merge!  more complex
                raise RuntimeError('not implemented yet')
                pred_fhs = list(map(lambda x: open(x, 'rb'), self.pred_fns))
                recs = [None] * len(pred_fhs)
                done = [False] * len(pred_fhs)
                nmerged = 0
                last_min_rec = -1
                with open(fn, 'wb') as ofh:
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
                                    pcor, ident = ln.rstrip().split(','.encode('utf-8'))[:2]
                                    recs[i] = (pcor, int(ident))
                            if recs[i] is not None and recs[i][1] < min_rec[1]:
                                min_rec, min_i = recs[i], i
                        if min_i == -1:
                            assert all(done)
                            break
                        nmerged += 1
                        assert min_rec[1] > last_min_rec, "%d,%d:%s" % (min_rec[1], last_min_rec, str(self.pred_fns))
                        last_min_rec = min_rec[1]
                        ofh.write(('%d,%0.3f\n' % (min_rec[1], pcor_to_mapq(float(min_rec[0])))).encode('utf-8'))
                        recs[min_i] = None
                assert nmerged == self.npredictions, (nmerged, self.npredictions)
                for fh in pred_fhs:
                    fh.close()

    def purge_temporaries(self):
        self.temp_man.remove_group(self.temp_group_name)

    def finalize(self, log=logging):
        """ Close prediction file handle.  If we have the information and flags
            needed for accuracy assessment, then do that too. """

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

            log.info('  Reordering')
            self.df.sort_values('pcor', ascending=False, inplace=True)

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
