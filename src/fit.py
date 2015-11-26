import random
import logging
import pandas
import numpy as np
import itertools
from predictions import MapqPredictions
from sklearn import cross_validation

__author__ = 'langmead'


def _np_deduping_indexes(m):
    b = np.ascontiguousarray(m).view(np.dtype((np.void, m.dtype.itemsize * m.shape[1])))
    _, idx, inv = np.unique(b, return_index=True, return_inverse=True)
    return idx, inv


class MapqFit:
    """ Encapsulates an object that fits models and makes predictions """

    def _df_to_mat(self, data, shortname, training, log=logging):
        """ Convert a data frame read with read_dataset into a matrix suitable
            for use with scikit-learn, and parallel vectors giving the
            original MAPQ predictions, the ids for the alignments (i.e. their
            line of origin) and whether or not the alignments are correct. """
        labs = []
        if training:
            assert shortname not in self.training_labs
            log.info('  Removing duplicate columns')
            for col in data:
                if col not in ['id', 'mapq', 'correct'] and data[col].nunique() > 1:
                    labs.append(col)
            to_remove = set()
            for x, y in itertools.combinations(labs, 2):
                if (data[x] == data[y]).all():
                    to_remove.add(y)
            for lab in to_remove:
                labs.remove(lab)
            self.training_labs[shortname] = labs
            if len(labs) == 0:
                raise RuntimeError('Error: all training records were identical')
        else:
            assert shortname in self.training_labs, (shortname, str(self.training_labs.keys()))
            labs = self.training_labs[shortname]
            for lab in labs:
                assert lab in data, "Column %s in training data, but not in test (%s)" % (lab, shortname)
        for lab in labs:
            assert not np.isnan(data[lab]).any()
        data_mat = data[labs].values
        assert not np.isinf(data_mat).any() and not np.isnan(data_mat).any()
        return data_mat, np.array(data['id']), np.array(data['mapq']), np.array(data['correct']), labs

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

    @staticmethod
    def postprocess_predictions(pcor_test, dataset_name, max_pcor=0.999999, log=logging):
        """ Deal with pcors equal to 1.0, which cause infinite MAPQs """
        mn, mx = min(pcor_test), max(pcor_test)
        if mx >= 1.0:
            if mn == mx:
                log.warning('All data points for %s are predicted correct; results unreliable' % dataset_name)
                pcor_test = [max_pcor] * len(pcor_test)
            max_noninf_pcor_test = max(filter(lambda x: x < 1.0, pcor_test))
            pcor_test = [max_noninf_pcor_test + 1e-6 if p >= 1.0 else p for p in pcor_test]
        return np.maximum(np.minimum(pcor_test, max_pcor), 0.)

    def _crossval_fit(self, mf_gen, x_train, y_train, dataset_shortname, use_oob=True, log=logging):
        """ Use cross validation to pick the best model from a
            collection of possible models (model_family) """
        mf = mf_gen()
        self.model_fam_name = mf.name
        scores = []

        def _oob_score(pred_):
            assert x_train.shape[0] == y_train.shape[0]
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
            log.debug("%s, %s=%0.3f, %s%s" % (dataset_shortname, 'oob' if use_oob else 'score', score, str(params), symbol))
        best_params, best_pred = mf.best_predictor()
        log.debug("BEST: %s, avg=%0.3f, %s" % (dataset_shortname, max(scores), str(best_params)))
        assert best_pred is not None
        return best_pred, best_params, max(scores)

    datasets = zip('dbcu', ['Discordant', 'Bad-end', 'Concordant', 'Unpaired'], [True, False, True, False])

    def _fit(self, dfs, log=logging, frac=1.0):
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
                                'the qsim software is making a mistake.  It could also '
                                'mean that, because of your data and reference genome, the aligner '
                                'can correctly resolve point of origin for all reads.  Treat '
                                'results circumspectly.' % train['correct'][0])
            # extract features, convert to matrix
            x_train, _, mapq_orig_train, y_train, self.col_names[ds] = self._df_to_mat(train, ds, True, log=log)
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
            self.trained_models[ds], self.trained_params[ds], self.crossval_avg[ds] = \
                self._crossval_fit(self.model_gen, x_train, y_train, ds)
            # fit training data with the model
            self.trained_models[ds].fit(x_train, y_train)

    def predict(self, dfs, temp_man,
                keep_data=False, keep_per_category=False, log=logging,
                dedup=False, training=False, calc_summaries=False, prediction_mem_limit=10000000):

        name = '_'.join(['overall', 'training' if training else 'test'])
        pred_overall = MapqPredictions(temp_man, name, calc_summaries=calc_summaries,
                                       prediction_mem_limit=prediction_mem_limit)
        pred_per_category = {}

        for ds, ds_long, paired in self.datasets:  # outer loop over ailgnment types
            if ds not in dfs:
                continue
            if keep_per_category:
                name = '_'.join([ds_long, 'training' if training else 'test'])
                pred_per_category[ds] = MapqPredictions(temp_man, name, calc_summaries=calc_summaries,
                                                        prediction_mem_limit=prediction_mem_limit)
            nchunk = 0
            for test_chunk in dfs.dataset_iter(ds):  # inner loop over chunks of rows
                nchunk += 1
                log.info('  Getting ready to make predictions for %s %s chunk %d, %d rows' %
                         ('training' if training else 'test', ds_long, nchunk, test_chunk.shape[0]))
                x_test, ids, mapq_orig_test, y_test, col_names = self._df_to_mat(test_chunk, ds, False, log=log)
                log.info('  Making predictions')
                if dedup:
                    idxs, invs = _np_deduping_indexes(x_test)
                    log.info('    Collapsed %d rows to %d distinct rows (%0.2f%%)' %
                             (len(invs), len(idxs), 100.0 * len(idxs) / len(invs)))
                    pcor = self.trained_models[ds].predict(x_test[idxs])[invs]  # make predictions
                else:
                    pcor = self.trained_models[ds].predict(x_test)  # make predictions
                pcor = np.array(self.postprocess_predictions(pcor, ds_long))
                data = x_test.tolist() if keep_data else None
                for prd in [pred_overall, pred_per_category[ds]] if keep_per_category else [pred_overall]:
                    prd.add(pcor, ids, ds, mapq_orig=mapq_orig_test, data=data, correct=y_test)

        log.info('Finalizing results for overall %s data (%d alignments)' %
                 ('training' if training else 'test', len(pred_overall.pcor)))
        pred_overall.finalize()
        if len(pred_per_category) > 1:
            for ds, pred in pred_per_category.iteritems():
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
                colnames.append(ds + '_oob_score')
            fh.write(','.join(colnames) + '\n')
            data = [self.model_fam_name, str(self.sample_fraction)]
            for ds, ds_long, paired in self.datasets:
                if ds in self.trained_shape:
                    data.append(':'.join(map(str, self.trained_params[ds])))
                    data.append(str(self.trained_shape[ds][0]))
                    data.append(str(self.trained_shape[ds][1]))
                    assert ds in self.trained_models
                    data.append(str(self.trained_models[ds].oob_score_))
                else:
                    data.extend(['NA', 'NA', '0', '0', '0'])
            fh.write(','.join(data) + '\n')

    def __init__(self,
                 dfs,  # dictionary of data frames, one per alignment type
                 model_gen,  # function that takes vector of hyperparameters, returns new model object
                 log=logging,
                 sample_fraction=1.0):  # fraction of training data to actually use
        self.model_gen = model_gen
        self.trained_models = {}
        self.crossval_avg = {}
        self.crossval_std = {}
        self.col_names = {}
        self.trained_params = {}
        self.trained_shape = {}
        self.training_labs = {}
        self.model_fam_name = None
        self.sample_fraction = sample_fraction
        self._fit(dfs, log=log, frac=sample_fraction)
