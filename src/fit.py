from __future__ import print_function
import random
import logging
import pandas
import numpy as np
import itertools
import resource
import gc
import os
import sys
import multiprocessing
from predictions import MapqPredictions
from sklearn import cross_validation
from mapq import pcor_to_mapq_np
try:
    import itertools.izip as zip
except ImportError:
    pass

__author__ = 'langmead'


def _np_deduping_indexes(m):
    b = np.ascontiguousarray(m).view(np.dtype((np.void, m.dtype.itemsize * m.shape[1])))
    _, idx, inv = np.unique(b, return_index=True, return_inverse=True)
    return idx, inv


def _get_peak_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 * 1024.0)


def _clamp_predictions(pcor, min_pcor=0.0, max_pcor=0.999999):
    return np.maximum(np.minimum(pcor, max_pcor), min_pcor)


def postprocess_predictions(pcor_test, dataset_name, max_pcor=0.999999, log=logging):
    """ Deal with pcors equal to 1.0, which cause infinite MAPQs """
    mn, mx = min(pcor_test), max(pcor_test)
    if mx >= 1.0:
        if mn == mx:
            log.warning('All data points for %s are predicted correct; results unreliable' % dataset_name)
            pcor_test = [max_pcor] * len(pcor_test)
        max_noninf_pcor_test = max(filter(lambda x: x < 1.0, pcor_test))
        pcor_test = [max_noninf_pcor_test + 1e-6 if p >= 1.0 else p for p in pcor_test]
    return _clamp_predictions(pcor_test, 0.0, max_pcor)


def _df_to_mat(data, shortname, training, training_labs, log=logging, include_mapq=False):
    """ Convert a data frame read with read_dataset into a matrix suitable
        for use with scikit-learn, and parallel vectors giving the
        original MAPQ predictions, the ids for the alignments (i.e. their
        line of origin) and whether or not the alignments are correct. """
    labs = []
    exclude_cols = ['id', 'correct', 'rname']
    if not include_mapq:
        exclude_cols.append('mapq')
    if training:
        assert shortname not in training_labs
        log.info('  Removing duplicate columns')
        for col in data:
            if col not in exclude_cols and data[col].nunique() > 1:
                labs.append(col)
        to_remove = set()
        for x, y in itertools.combinations(labs, 2):
            if (data[x] == data[y]).all():
                to_remove.add(y)
        for lab in to_remove:
            labs.remove(lab)
        training_labs[shortname] = labs
        if len(labs) == 0:
            raise RuntimeError('Error: all training records were identical')
    else:
        assert shortname in training_labs, (shortname, str(training_labs.keys()))
        labs = training_labs[shortname]
        for lab in labs:
            assert lab in data, "Column %s in training data, but not in test (%s)" % (lab, shortname)
    for lab in labs:
        assert not np.isnan(data[lab]).any()
    data_mat = data[labs].values
    assert not np.isinf(data_mat).any() and not np.isnan(data_mat).any()
    return data_mat, np.array(data['id']), np.array(data['mapq']), np.array(data['correct']), labs


_prediction_worker_queue = multiprocessing.Queue()
_prediction_worker_trained_models = None
_prediction_worker_pred_overall = None
_prediction_worker_log = None


def _prediction_worker(my_test_chunk, training, training_labs, ds,
                       ds_long, pred_per_category, dedup,
                       keep_per_category=False, keep_data=False,
                       multiprocess=True, include_mapq=False):
    gc.collect()
    log = _prediction_worker_log
    trained_model = _prediction_worker_trained_models[ds]
    pred_overall = _prediction_worker_pred_overall
    log.info('  PID %d making predictions for %s %s chunk, %d rows (peak mem=%0.2fGB)' %
             (os.getpid(), 'training' if training else 'test', ds_long, my_test_chunk.shape[0], _get_peak_gb()))
    log.info('    Loading data')
    x_test, ids, mapq_orig_test, y_test, col_names = \
        _df_to_mat(my_test_chunk, ds, False, training_labs, log=log, include_mapq=include_mapq)
    del my_test_chunk
    gc.collect()
    if dedup:
        log.info('    Done loading data; collapsing and making predictions')
        idxs, invs = _np_deduping_indexes(x_test)
        log.info('    Collapsed %d rows to %d distinct rows (%0.2f%%)' %
                 (len(invs), len(idxs), 100.0 * len(idxs) / len(invs)))
        pcor = trained_model.predict(x_test[idxs])[invs]  # make predictions
    else:
        log.info('    Done loading data; making predictions')
        pcor = trained_model.predict(x_test)  # make predictions
    data = x_test.tolist() if keep_data else None
    del x_test
    gc.collect()
    log.info('    Done making predictions; about to postprocess (peak mem=%0.2fGB)' % _get_peak_gb())
    pcor = np.array(postprocess_predictions(pcor, ds_long))
    log.info('    Done postprocessing; adding to tally (peak mem=%0.2fGB)' % _get_peak_gb())
    if multiprocess:
        return pcor, ids, ds, mapq_orig_test, data, y_test
    else:
        for prd in [pred_overall, pred_per_category[ds]] if keep_per_category else [pred_overall]:
                prd.add(pcor, ids, ds, mapq_orig=mapq_orig_test, data=data, correct=y_test)
    log.info('    Done; peak mem usage so far = %0.2fGB' %
             (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 * 1024.0)))


def _prediction_worker_star(a_b):
    return _prediction_worker(*a_b)


def _prediction_worker_init(log):
    log.info('  Initializing worker process with PID %d' % (os.getpid()))


class MapqFit:
    """ Encapsulates an object that fits models and makes predictions """

    @staticmethod
    def _subsample(x_train, mapq_orig_train, y_train, sample_fraction):
        """ Return a random subset of the data, MAPQs and labels.  Size of
            subset given by sample_fraction. """
        n_training_samples = x_train.shape[0]
        assert x_train.shape[0] == y_train.shape[0]
        if sample_fraction < 1.0:
            sample_indexes = random.sample(range(n_training_samples), int(round(n_training_samples * sample_fraction)))
            x_train = x_train[sample_indexes, ]
            y_train = y_train[sample_indexes, ]
            mapq_orig_train = mapq_orig_train[sample_indexes, ]
        return x_train, mapq_orig_train, y_train

    def _fit_and_possibly_reweight_and_refit(self, predictor, x_train, y_train,
                                             reweight_ratio=1.0,
                                             reweight_mapq=False,
                                             reweight_mapq_offset=10.0):
        """ Fit, then, if request, use predictions to weigh samples and re-fit.
            Tends to force the model to fit the high-MAPQ points better so we have fewer
            incorrect alignments with high MAPQ. """
        predictor.fit(x_train, y_train)
        if reweight_ratio > 1.0:
            y_pred = predictor.predict(x_train)
            lower = 1.0 / reweight_ratio
            y_pred = lower + y_pred * (1.0 - lower)
            predictor.fit(x_train, y_train, y_pred)
        elif reweight_mapq:
            assert reweight_mapq_offset >= 0
            y_pred = pcor_to_mapq_np(_clamp_predictions(predictor.predict(x_train)))
            predictor.fit(x_train, y_train, y_pred + reweight_mapq_offset)

    def _crossval_fit(self, mf_gen, x_train, y_train, dataset_shortname, use_oob=True, log=logging,
                      reweight_ratio=1.0, reweight_mapq=False, reweight_mapq_offset=10.0):
        """ Use cross validation to pick the best model from a
            collection of possible models (model_family) """
        mf = mf_gen()
        self.model_fam_name = mf.name
        scores = []

        def _oob_score(pred_):
            assert x_train.shape[0] == y_train.shape[0]
            self._fit_and_possibly_reweight_and_refit(pred_, x_train, y_train,
                                                      reweight_ratio=reweight_ratio,
                                                      reweight_mapq=reweight_mapq,
                                                      reweight_mapq_offset=reweight_mapq_offset)
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
            log.debug("%s, %s=%0.3f, %s%s" % (dataset_shortname, 'oob' if use_oob else 'score', score, str(params), symbol))
        best_params, best_pred = mf.best_predictor()
        log.info("BEST: %s, avg=%0.3f, %s, using %s" % (dataset_shortname, max(scores), str(best_params),
                                                        'OOB' if use_oob else 'cross validation'))
        assert best_pred is not None
        return best_pred, best_params, max(scores)

    datasets = list(zip('dbcu', ['Discordant', 'Bad-end', 'Concordant', 'Unpaired'], [True, False, True, False]))

    def _fit(self, dfs, log=logging, frac=1.0, heap_profiler=None, include_mapq=False,
             reweight_ratio=1.0, reweight_mapq=False, reweight_mapq_offset=10.0, no_oob=False):
        """ Train one model per training table. Optionally subsample training
            data first. """
        for ds, ds_long, paired in self.datasets:
            if ds not in dfs:
                continue  # empty
            train = pandas.concat([x for x in dfs.dataset_iter(ds)])
            if train.shape[0] == 0:
                continue  # empty
            if train['correct'].nunique() == 1:
                logging.warning('Warning: All training data has correct=%d.  This might mean '
                                'the qtip software is making a mistake.  It could also '
                                'mean that, because of your data and reference genome, the aligner '
                                'can correctly resolve point of origin for all reads.  Treat '
                                'results circumspectly.' % train['correct'][0])
            # extract features, convert to matrix
            x_train, _, mapq_orig_train, y_train, self.col_names[ds] = \
                _df_to_mat(train, ds, True, self.training_labs, log=logging, include_mapq=False)
            assert x_train.shape[0] == y_train.shape[0]
            assert x_train.shape[1] > 0
            # optionally subsample
            if frac < 1.0:
                log.info('  Sampling %0.2f%% of %d rows of %s records' % (100.0 * frac, train.shape[0], ds_long))
                x_train, mapq_orig_train, y_train = \
                    self._subsample(x_train, mapq_orig_train, y_train, frac)
                log.info('  Now has %d rows' % x_train.shape[0])
            # use cross-validation to pick a model
            log.info('Fitting %d %s training records; %d features each' % (x_train.shape[0], ds_long, x_train.shape[1]))
            assert x_train.shape[0] == y_train.shape[0]
            self.trained_shape[ds] = x_train.shape
            self.trained_models[ds], self.trained_params[ds], self.model_score[ds] = \
                self._crossval_fit(self.model_gen, x_train, y_train, ds,
                                   use_oob=self.model_gen().calculates_oob() and not no_oob)
            log.info('    Chose parameters: %s' % str(self.trained_params[ds]))
            self._fit_and_possibly_reweight_and_refit(self.trained_models[ds], x_train, y_train,
                                                      reweight_ratio=reweight_ratio,
                                                      reweight_mapq=reweight_mapq,
                                                      reweight_mapq_offset=reweight_mapq_offset)
            del x_train
            del y_train
            gc.collect()
            log.info('    Done; peak mem usage so far = %0.2fGB' % _get_peak_gb())
            if heap_profiler is not None:
                print(heap_profiler.heap(), file=sys.stderr)

    def predict(self, dfs, temp_man,
                keep_data=False, keep_per_category=False, log=logging,
                dedup=False, training=False, calc_summaries=False, prediction_mem_limit=10000000,
                heap_profiler=None, include_mapq=False,
                multiprocess=False, n_multi=8):

        global _prediction_worker_trained_models
        global _prediction_worker_pred_overall
        global _prediction_worker_log

        name = '_'.join(['overall', 'training' if training else 'test'])
        pred_overall = MapqPredictions(temp_man, name, calc_summaries=calc_summaries,
                                       prediction_mem_limit=prediction_mem_limit)
        pred_per_category = {}
        log.info('  Created overall MapqPredictions (peak mem=%0.2fGB)' % _get_peak_gb())

        _prediction_worker_trained_models = self.trained_models
        _prediction_worker_pred_overall = pred_overall
        _prediction_worker_log = log

        p = None
        if multiprocess:
            assert n_multi is None or n_multi > 0
            p = multiprocessing.Pool(n_multi, _prediction_worker_init, (log,))

        for ds, ds_long, paired in self.datasets:  # outer loop over alignment types
            if ds not in dfs:
                continue
            if keep_per_category:
                name = '_'.join([ds_long, 'training' if training else 'test'])
                pred_per_category[ds] = MapqPredictions(temp_man, name, calc_summaries=calc_summaries,
                                                        prediction_mem_limit=prediction_mem_limit)
                log.info('  Created %s MapqPredictions (peak mem=%0.2fGB)' % (name, _get_peak_gb()))

            # assert type(training) is BooleanType
            # assert isinstance(self.trained_models[ds], object)
            # assert type(self.training_labs) is DictType
            # assert type(ds) is StringType
            # assert type(ds_long) is StringType
            # assert isinstance(pred_overall, MapqPredictions)
            # assert type(pred_per_category) is DictType
            # assert type(dedup) is BooleanType
            # assert isinstance(log, object)
            # assert type(keep_per_category) is BooleanType
            # assert type(keep_data) is BooleanType
            # assert type(include_mapq) is BooleanType

            if multiprocess:
                try:
                    r = itertools.repeat
                    for tup in p.imap(_prediction_worker_star,
                                      zip(dfs.dataset_iter(ds),
                                                     r(training),
                                                     r(self.training_labs),
                                                     r(ds),
                                                     r(ds_long),
                                                     r(pred_per_category),
                                                     r(dedup),
                                                     r(keep_per_category),
                                                     r(keep_data),
                                                     r(True),
                                                     r(include_mapq))):
                        pcor, ids, ds, mapq_orig_test, data, y_test = tup
                        for prd in [pred_overall, pred_per_category[ds]] if keep_per_category else [pred_overall]:
                            prd.add(pcor, ids, ds, mapq_orig=mapq_orig_test, data=data, correct=y_test)

                except KeyboardInterrupt:
                    p.terminate()
                    p.join()
                    raise
            else:
                for test_chunk in dfs.dataset_iter(ds):
                    _prediction_worker(test_chunk, training, self.training_labs,
                                       ds, ds_long, pred_per_category, dedup,
                                       keep_per_category=keep_per_category, keep_data=keep_data,
                                       multiprocess=False, include_mapq=include_mapq)

        log.info('Finalizing results for overall %s data (%d alignments)' %
                 ('training' if training else 'test', len(pred_overall.pcor)))
        pred_overall.finalize()
        if len(pred_per_category) > 1:
            for ds, pred in pred_per_category.items():
                log.info('Finalizing results for "%s" %s data (%d alignments)' %
                         (ds, 'training' if training else 'test', len(pred.pcor)))
                pred.finalize()
        log.info('Done')

        if keep_per_category:
            return pred_overall, pred_per_category
        else:
            return pred_overall

    def write_feature_importances(self, prefix):
        """
        Write feature importances for each model to an appropriately-named
        file with given prefix.
        """
        fi_colnames, fi_vals = [], []
        for ds, model in sorted(self.trained_models.items()):
            with open(prefix + '_' + ds + '.csv', 'w') as fh:
                assert len(self.col_names[ds]) == len(model.feature_importances_)
                ranks = np.argsort(model.feature_importances_)[::-1] + 1
                inv_ranks = [0] * len(ranks)
                for i, r in enumerate(ranks):
                    assert inv_ranks[r-1] == 0
                    inv_ranks[r-1] = i+1
                assert min(inv_ranks) == 1
                assert max(inv_ranks) == len(ranks)
                i = 0
                fh.write('feature,importance,rank\n')
                for im, r in zip(model.feature_importances_, inv_ranks):
                    fi_colnames.append(ds + '_' + self.col_names[ds][i])
                    fi_vals.append(im)
                    fh.write('%s,%0.4f,%d\n' % (self.col_names[ds][i], im, r))
                    i += 1

    def write_parameters(self, prefix):
        """
        Write out the hyperparameters for each model to an appropriately-named
        file with given prefix.
        """
        with open(prefix + '.csv', 'w') as fh:
            colnames = ['model_type', 'subsampling_fraction']
            for ds, _, _, in self.datasets:
                colnames.append(ds + '_model_params')
                colnames.append(ds + '_training_rows')
                colnames.append(ds + '_training_cols')
                colnames.append(ds + '_model_score')
            fh.write(','.join(colnames) + '\n')
            data = [self.model_fam_name, str(self.sample_fraction)]
            for ds, ds_long, paired in self.datasets:
                if ds in self.trained_shape:
                    data.append(':'.join(map(str, self.trained_params[ds])))
                    data.append(str(self.trained_shape[ds][0]))
                    data.append(str(self.trained_shape[ds][1]))
                    assert ds in self.trained_models
                    data.append(str(self.model_score[ds]))
                else:
                    data.extend(['NA', '0', '0', '0'])
            fh.write(','.join(data) + '\n')

    def __init__(self,
                 dfs,  # dictionary of data frames, one per alignment type
                 model_gen,  # function that takes vector of hyperparameters, returns new model object
                 log=logging,
                 sample_fraction=1.0,  # fraction of training data to actually use
                 heap_profiler=None,
                 include_mapq=False,
                 reweight_ratio=1.0,
                 reweight_mapq=False,
                 reweight_mapq_offset=10.0,
                 no_oob=False):
        self.model_gen = model_gen
        self.trained_models = {}
        self.crossval_std = {}
        self.col_names = {}
        self.trained_params = {}
        self.model_score = {}
        self.trained_shape = {}
        self.training_labs = {}
        self.model_fam_name = None
        self.sample_fraction = sample_fraction
        self.q = None
        self._fit(dfs, log=log, frac=sample_fraction, heap_profiler=heap_profiler, include_mapq=include_mapq,
                  reweight_ratio=reweight_ratio, reweight_mapq=reweight_mapq,
                  reweight_mapq_offset=reweight_mapq_offset, no_oob=no_oob)
