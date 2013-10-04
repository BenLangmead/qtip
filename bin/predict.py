import os
import logging
import math
import random
import abc
import itertools
import pandas as pd
import statsmodels.formula.api as sm
import numpy as np
import contextlib
from scipy.optimize import minimize, optimize
from sklearn.cross_validation import StratifiedKFold

# Modules that are part of the tandem simulator
from samples import UnpairedTuple, PairedTuple, Dataset

@contextlib.contextmanager
def noOutput():
    savestdout = sys.stdout
    savestderr = sys.stderr
    class Devnull(object):
        def write(self, _): pass
    sys.stdout = Devnull()
    sys.stderr = Devnull()
    yield
    sys.stdout = savestdout
    sys.stderr = savestderr

class Result(object):
    ''' Result of a set of predictions where correctness is known '''
    
    def __init__(self, pcors, cors, mapqize=False, pcorize=False, round=False, noise=True):
        if pcorize:
            # Caller gave us MAPQs, but we want pcors
            pcors = [ 1 - 10.0 ** (-0.1 * x) for x in pcors ]
        if mapqize:
            # Caller gave us pcors, but we want MAPQs
            pcors = [ 100.0 if p == 1.0 else -10 * math.log10(1.0 - p) for p in pcors ]
        if round:
            # Round pcors to nearest integer (only sensible if they're MAPQs)
            pcors = [ int(x) for x in pcors ]
        self.cors = cors
        self.pcors = pcors
        if noise:
            mn, mx = min(pcors), max(pcors)
            self.noisyPcors = np.minimum(mx, np.maximum(mn, np.add(np.random.uniform(0, 1e-10, len(pcors)), pcors)))
            self.items = zip(self.noisyPcors, cors)
        else:
            self.items = zip(self.pcors, cors)
        self.items.sort(key=lambda x: (x[0], -x[1]))
    
    def __iter__(self):
        return iter(self.items)
    
    #
    # Useful plots
    #
    
    def mapqBinPlot(self):
        pass

class RankingErrorObjective(object):
    ''' Return ranking error associated with result '''
    def score(self, result):
        i, tot = 1, 0
        for _, cor in iter(result):
            if not cor:
                tot += i
            i += 1
        return tot

class SmoothedOutcomeObjective(object):
    ''' An objective that decreases as the pcor predictions line up better
        with "actual" probabilities of being correct.  The "actual"
        probabilities are calculated by sorting the cor vector according to
        ascending pcor (with noise), then smoothing using running mean. '''
    
    def __init__(self, window=50):
        self.window = window # running mean window
    
    @classmethod
    def scoreHelper(cls, pcor, cor, window):
        smoothcor = pd.rolling_mean(np.array(cor, dtype=np.float32), window)
        ssd = np.sum((cor - np.mean(cor))**2)
        return sum([(0 if math.isnan(sc) else (pc-sc)**2) for pc, sc in itertools.izip(pcor, smoothcor)]) / ssd
    
    def score(self, result):
        ''' Compute and return a score equal to sum of squared differences
            between pcors and smoothed cor values. '''
        pcor, cor = zip(*result.items)
        return self.scoreHelper(pcor, cor, self.window)
    
    def plot(self, result, pdf):
        pcor, cor = zip(*result.items)
        smoothcor = self.smoothCor(cor)
        plt.figure(figsize=(7, 7))
        plt.title('Smooth outcome objective plot (score=%0.2f)' % self.score(result))
        ran = range(len(smoothcor))
        plt.plot(ran, pcor, 'bo', ran, smoothcor, 'g-', ran, cor, 'ko')
        plt.savefig(pdf, format='pdf')
        plt.close()

class DiffScale(object):
    def __init__(self, scalefact = 1.5):
        self.scalefact = scalefact
    
    def rescale(self, ds):
        ds['diff'] = np.minimum(self.scalefact, (ds['best1'] - ds['best2']).fillna(self.scalefact))

class LogitPredictor(object):
    ''' Use logistic regression model with two '''
    
    def __init__(self,
                 dscale=DiffScale(4.0),
                 logB=False,
                 logD=False,
                 ladd=0.5):
        self.dscale = dscale
        self.ladd = ladd
        self.res = None
        self.logB, self.logD = logB, logD
    
    def train(self, trainingData, objective, cv_fold=5):
        td = trainingData.copy()
        for col in [ 'correct', 'best1', 'best2' ]:
            assert col in td
        self.dscale.rescale(td)
        logging.info('Fitting logit model')
        bestStr, diffStr = 'best1', 'diff'
        if self.logB: bestStr = 'np.log(%f + best1)' % self.ladd
        if self.logD: diffStr = 'np.log(%f + diff)' % self.ladd
        form = 'correct ~ %s * %s' % (bestStr, diffStr)
        self.res = sm.logit(formula=form, data=td).fit()
        
        logging.debug('Doing stratified %d-fold cross validation' % cv_fold)
        skf = StratifiedKFold(td['correct'], cv_fold)
        foldErrs = []
        for traini, testi in skf:
            trainTd, testTd = td.ix[traini], td.ix[testi]
            try:
                mod = sm.logit(formula=form, data=trainTd).fit()
                res = Result(mod.predict(testTd), testTd['correct'])
                foldErrs.append(objective.score(res))
            except np.linalg.LinAlgError:
                logging.warning('Linear algebra error when fitting logit model!!!')
        return foldErrs
    
    def predict(self, testData):
        td = testData.copy()
        self.dscale.rescale(td)
        # Predict.  No need to apply log transformation to data -
        # self.res remembers that
        logging.info('Predicting with logit model')
        return self.res.predict(td, transform=True)

if __name__ == "__main__":
    import sys
    import unittest
    import argparse
    
    parser = argparse.ArgumentParser(\
        description='Fit model, make predictions.')
    
    # Output-related options
    parser.add_argument(\
        '--output-directory', metavar='path', type=str, required=True, help='Write outputs to this directory')
    parser.add_argument(\
        '--write-training-predictions', action='store_const', const=True,
        default=False, help='Write predictions + correctness info for training data to a csv file')
    parser.add_argument(\
        '--write-input-predictions', action='store_const', const=True,
        default=False, help='Write predictions + correctness info (if available) for input (test) data to a csv file')
    parser.add_argument(\
        '--write-all', action='store_const', const=True,
        default=False, help='Same as specifying all --write-* options')
    
    parser.add_argument(\
        '--seed', metavar='int', type=int, default=6277, help='Pseudo-random seed')
    parser.add_argument(\
        '--training', metavar='path', type=str, required=True, help='Input training data')
    parser.add_argument(\
        '--test', metavar='path', type=str, required=True, help='Input test data')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False, help='Be talkative')
    
    args = parser.parse_args()
    
    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG if args.verbose else logging.INFO)
    
    # Create output directory if needed
    outdir = args.output_directory
    if not os.path.isdir(outdir):
        try: os.makedirs(outdir)
        except OSError as exception:
            if exception.errno != errno.EEXIST: raise
    
    logging.info('Loading training data')
    trainingData = Dataset()
    trainingData.load(args.training)
    train = {}
    train['U'], train['M'], train['P'] = trainingData.toDataFrames()
    
    rankingObj = RankingErrorObjective()
    smoothObj = SmoothedOutcomeObjective()
    
    logging.info('Creating and fitting models')
    pred = {}
    for type in 'UMP':
        if not train[type].empty:
            logging.info('  "%s" model' % type)
            cv_fold = 2
            
            # Pick parameters based on an optimization
            def objf(x):
                random.seed(args.seed)
                np.random.seed(args.seed)
                pred = LogitPredictor(dscale=DiffScale(x))
                sc = np.mean(pred.train(train[type], smoothObj, cv_fold=cv_fold))
                logging.debug('  tried x=%0.6f, got mean cv=%0.6f' % (x, sc))
                return sc
            
            #minresult = minimize(objf, 4.0, method='Nelder-Mead', options={'maxiter':15})
            diffScale = optimize.brute(objf, ((1.0, 21.0, 2.0),), finish=None)
            logging.info('Best diffScale=%0.6f' % diffScale)
            random.seed(args.seed)
            np.random.seed(args.seed)
            pred[type] = LogitPredictor(dscale=DiffScale(diffScale))
            pred[type].train(train[type], smoothObj)
            assert pred[type] is not None
            train[type]['pcor'] = pred[type].predict(train[type])
            res = Result(train[type]['pcor'], train[type]['correct'])
            alignerRes = Result(train[type]['mapq'], train[type]['correct'], pcorize=True)
            myRank = rankingObj.score(res)
            alignerRank = rankingObj.score(alignerRes)
            logging.info('Score on training: ranking ratio=%0.4f, smooth=%0.4f' % (float(myRank)/alignerRank, smoothObj.score(res)))
            
            if args.write_all or args.write_training_predictions:
                train[type].to_csv(os.path.join(outdir, 'training_predictions_%s.csv' % type), index=False)
    
    logging.info('Loading test data')
    testData = Dataset()
    testData.load(args.test)
    test = {}
    test['U'], test['M'], test['P'] = testData.toDataFrames()
    
    logging.info('Predicting using models')
    for type in 'UMP':
        if not test[type].empty:
            assert not train[type].empty
            logging.info('  "%s" model' % type)
            assert 'pcor' not in test[type]
            # Predict pcor for all the input reads!
            test[type]['pcor'] = pred[type].predict(test[type])
            if any(map(lambda x: x is not None, test[type]['correct'])):
                logging.info('Input test data also has correctness information')
                res = Result(test[type]['pcor'], test[type]['correct'])
                alignerRes = Result(test[type]['mapq'], test[type]['correct'], pcorize=True)
                myRank = rankingObj.score(res)
                alignerRank = rankingObj.score(alignerRes)
                logging.info('Score on training: ranking ratio=%0.6f, smooth=%0.6f' % (float(myRank)/alignerRank, smoothObj.score(res)))
            
            if args.write_all or args.write_input_predictions:
                test[type].to_csv(os.path.join(outdir, 'input_predictions_%s.csv' % type), index=False)
    
    logging.info('DONE')
    
    class Test(unittest.TestCase):
        def test_rankingError(self):
            result = Result([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], [True, False, False, True, False, True, True, False], noise=False)
            r = RankingErrorObjective()
            self.assertEqual(2 + 3 + 5 + 8, r.score(result))

    unittest.main(argv=[sys.argv[0]])
