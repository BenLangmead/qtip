"""
Given a directory with output from ts.py, predict new MAPQs.
"""

__author__ = 'langmead'

import pandas
import os
import sys
import logging
import gc
try:
    import cPickle as pickle
except ImportError:
    import pickle
try:
    import itertools.imap as map
except ImportError:
    pass
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
from plots import plot_drop_rate, plot_drop_rate_difference, plot_subsampling_series, bucket_error_plot
from metrics import mseor, ranking_error, auc, roc_table


VERSION = '0.2.0'


def mkdir_quiet(dr):
    """ Create directory if needed; don't complain """
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def output_top_incorrect(pred, odir, args):
    """ Output incorrect alignments with largest predicted MAPQ """
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
                        ss_fit = pickle.load(fh)
                else:
                    # Fit model
                    logging.info('      Fitting')
                    ss_fit = MapqFit(df_training, fam, random_seed=my_seed, logger=logging.info,
                                     sample_fraction=fraction, include_ztzs=not args['ignore_ztzs'])
                    if args['serialize_fit']:
                        logging.info('      Serializing fit object')
                        with open(my_fit_fn, 'wb') as ofh:
                            pickle.dump(ss_fit, ofh, 2)

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
            fit = pickle.load(fh)
    else:
        logging.info('Fitting and making predictions')
        my_seed = hash(str(args['seed']) + '-1') % 4294967296
        assert my_seed >= 0
        fit = MapqFit(df_training, fam, random_seed=my_seed, logger=logging.info, include_ztzs=not args['ignore_ztzs'])
        if args['serialize_fit']:
            logging.info('Serializing fit object')
            with open(fit_fn, 'wb') as ofh:
                pickle.dump(fit, ofh, 2)

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
