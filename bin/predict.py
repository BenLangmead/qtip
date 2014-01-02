"""
predict.py
"""

import os
import logging
from pandas import DataFrame
from pandas.io.parsers import read_csv
from sklearn.ensemble import RandomForestRegressor
from math import log10


def read_dataset(prefix):
    dfs = {}
    for name, suf, paired in [('u', '_unp.csv', False), ('m', '_mates.csv', False), ('c', '_conc.csv', True)]:
        fn = prefix + suf
        if os.path.exists(fn):
            dfs[name] = read_csv(fn, quoting=2)
        elif os.path.exists(fn + '.gz'):
            dfs[name] = read_csv(fn + '.gz', quoting=2, compression='gzip')
        elif os.path.exists(fn + '.bz2'):
            dfs[name] = read_csv(fn + '.bz2', quoting=2, compression='bz2')
        else:
            raise RuntimeError('No such file: "%s"' % fn)

    for df in dfs.itervalues():
        if df['correct'].count() == len(df['correct']):
            df['correct'] = df['correct'].map(lambda x: 1 if x == 'T' else 0)

    for df in [dfs['u'], dfs['m']]:
        diffv = df.maxv - df.minv
        df['bestnorm'] = (1.0 * df.best - df.minv) / diffv
        secbestnorm = (1.0 * df.secbest - df.minv) / diffv
        df['secbestnorm'] = secbestnorm.fillna(secbestnorm.min())
        df['diffnorm'] = df.bestnorm - df.secbestnorm

    conc = dfs['c']

    diffv_1 = conc.maxv1 - conc.minv1
    conc['bestnorm1'] = (1.0 * conc.best1 - conc.minv1) / diffv_1
    secbestnorm_1 = (1.0 * conc.secbest1 - conc.minv1) / diffv_1
    conc['secbestnorm1'] = secbestnorm_1.fillna(secbestnorm_1.min())
    conc['diffnorm1'] = conc['bestnorm1'] - conc['secbestnorm1']

    diffv_2 = conc.maxv2 - conc.minv2
    conc['bestnorm2'] = (1.0 * conc.best2 - conc.minv2) / diffv_2
    secbestnorm_2 = (1.0 * conc.secbest2 - conc.minv2) / diffv_2
    conc['secbestnorm2'] = secbestnorm_2.fillna(secbestnorm_2.min())
    conc['diffnorm2'] = conc['bestnorm2'] - conc['secbestnorm2']

    diff_conc = diffv_1 + diffv_2
    conc['bestnormconc'] = (1.0 * conc.best1 + conc.best2 - conc.minv1 - conc.minv2) / diff_conc
    secbest_conc = (1.0 * conc.secbest1 + conc.secbest2 - conc.minv1 - conc.minv2) / diff_conc
    conc['secbestnormconc'] = secbest_conc.fillna(secbest_conc.min())
    conc['diffnormconc'] = conc['bestnormconc'] - conc['secbestnormconc']

    return dfs


def mapqize(p):
    return map(lambda x: -10.0 * log10(1.0 - x), p)


def go(args):
    # Load training and test
    logging.info('Loading training data')
    train_dfs = read_dataset(args.training_prefix)
    logging.info('Loading test data')
    test_dfs = read_dataset(args.test_prefix)

    logging.info('Fitting models')
    datasets = [('Unpaired', train_dfs['u'], test_dfs['u'], False),
                ('Mate', train_dfs['m'], test_dfs['m'], False),
                ('Concordant', train_dfs['c'], test_dfs['c'], True)]
    models = [('RFR', RandomForestRegressor(n_estimators=10, max_depth=5))]
    for dataset_name, train, test, paired in datasets:
        if paired:
            x_train = train[['bestnorm1', 'diffnorm1', 'diffnormconc']].values
            x_test = test[['bestnorm1', 'diffnorm1', 'diffnormconc']].values
            mapq_orig_train = train['mapq1']
            mapq_orig_test = test['mapq1']
        else:
            x_train = train[['bestnorm', 'diffnorm']].values
            x_test = test[['bestnorm', 'diffnorm']].values
            mapq_orig_train = train['mapq']
            mapq_orig_test = test['mapq']
        y_train = map(lambda x: x == 1, train.correct)
        y_test = map(lambda x: x == 1, test.correct)
        for model_name, model in models:
            logging.info('  Fitting "%s" on test dataset "%s"' % (model_name, dataset_name))
            model.fit(x_train, y_train)
            pcor_test = model.predict(x_test)
            mapq_test = mapqize(pcor_test)
            logging.info('  Writing results')
            result_test_df = DataFrame.from_items([('pcor', pcor_test), ('mapq', mapq_test),
                                                   ('orig', mapq_orig_test), ('correct', y_test)])
            result_test_df.to_csv(args.test_prefix + '_' + model_name + '.csv', index=False)
            if args.training_results:
                logging.info('  Fitting "%s" on training dataset "%s"' % (model_name, dataset_name))
                pcor_train = model.predict(x_train)
                mapq_train = mapqize(pcor_train)
                logging.info('  Writing results')
                result_train_df = DataFrame.from_items([('pcor', pcor_train), ('mapq', mapq_train),
                                                        ('orig', mapq_orig_train), ('correct', y_train)])
                result_train_df.to_csv(args.training_prefix + '_' + model_name + '.csv', index=False)

    logging.info('Done')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Fit models, etc.')

    parser.add_argument('--training-prefix', metavar='path', type=str, default='training',
                        help='Prefix for files with training data')
    parser.add_argument('--test-prefix', metavar='path', type=str, default='test',
                        help='Prefix for files with test data')
    parser.add_argument('--training-results', action='store_const', const=True, default=False,
                        help='Use model to fit on training data and emit results')
    parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')

    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG if args.verbose else logging.INFO)

    if args.profile:
        import cProfile
        cProfile.run('go(args)')
    else:
        go(args)
