import abc
import pandas as pd
import statsmodels.formula.api as sm
from samples import UnpairedTuple, PairedTuple, Dataset
import logging

class Result(object):
    ''' Result of a set of predictions where correctness is known '''
    
    def __init__(self, pcors, cors):
        self.items = zip(pcors, cors)
        self.items.sort()
    
    def __iter__(self):
        return iter(self.items)

class Objective(object):
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def score(self, result):
        return 0.0

class RankingErrorObjective(Objective):
    ''' Return ranking error associated with result '''
    def score(self, result):
        i, tot = 0, 0
        for _, cor in iter(result):
            if not cor:
                tot += i
            i += 1
        return tot

class SmoothFitObjective(Objective):
    ''' Return error calculated by drawing a smooth curve through the cor data
        (sorted by pcor) and summing squared differences between the pcors and
        the curve. '''
    def __init__(self, smoothParam):
        self.smoothParam = smoothParam
    
    def score(self, result):
        pass

class DiffScale(object):
    def __init__(self, scalefact = 1.5):
        self.scalefact = scalefact
    
    def rescale(self, ds):
        ds['diff'] = (ds['best1'] - ds['best2']).fillna(self.scalefact)

class LogitPredictor(object):
    ''' Use logistic regression model with two '''
    
    def __init__(self, dscale=DiffScale(1.5), ladd=0.1):
        self.dscale = dscale
        self.ladd = ladd
    
    def train(self, trainingData):
        td = trainingData.copy()
        for col in [ 'correct', 'best1', 'best2' ]:
            assert col in td
        self.dscale.rescale(td)
        logging.info('Fitting logit model')
        #print td['correct']
        #print td['best1']
        res = sm.logit(formula='correct ~ best1 * diff', data=td).fit()
    
    def predict(self, testData):
        pass

class StratifiedPredictor(object):
    ''' Wrap another predictor and make it 'stratified' by read length. '''
    def __init__(self):
        pass

class KNN:
    
    def __init__(self, k=1000, leafSz=20, cacheSz=4096):
        self.leafSz = leafSz
        self.k = k
        self.kdt = None
        self.tups = []
        self.tupMap = {}
        self.hits, self.misses = 0, 0
        self.ncor, self.ntot = 0, 0
        
        # Very simple initial cache: remember last query and answer
        self.prevTup, self.prevProb = None, None
        
        # Cache is a doubly-linked list
        # Link layout:     [PREV, NEXT, KEY, RESULT]
        self.root = root = [None, None, None, None]
        self.cache = cache = {}
        last = root
        for i in range(cacheSz):
            key = object()
            cache[key] = last[1] = last = [last, root, key, None]
        root[0] = last
    
    def fit(self, tups, corrects):
        """ Add a collection of training tuples, each with associated
            boolean indicating whether tuple corresponds to a correct
            or incorrect alignment. """
        self.ntot = len(tups)
        for tup, correct in zip(tups, corrects):
            correcti = 1 if correct else 0
            self.ncor += correcti
            if tup in self.tupMap:
                self.tupMap[tup][0] += correcti # update # correct
                self.tupMap[tup][1] += 1        # update total
            else:
                self.tups.append(tup)
                self.tupMap[tup] = [correcti, 1]
        assert self.ntot > 0
        self.finalize()
    
    def finalize(self):
        """ Called when all training data tuples have been added. """
        assert self.kdt is None
        self.kdt = KDTree(self.tups, leafsize=self.leafSz)
    
    def __probCorrect(self, tup):
        """ Query k nearest neighbors of test tuple (tup); that is, the
            k training tuples most "like" the test tuple.  Return the
            fraction of the neighbors that are correct. """
        
        # Fast path: test point is on top of a stack of training points
        # that is already >= k points tall
        if tup in self.tupMap and self.tupMap[tup][1] >= self.k:
            ncor, ntot =  self.tupMap[tup]
            return float(ncor) / ntot
        
        radius = 50
        bestEst = None
        ntot = 0
        wtups = [] # weighted neighbor tuples
        maxDist, minDist = 0.0, 0.0
        while radius < 1000:
            neighbors = self.kdt.query_ball_point([tup], radius, p=2.0)
            for idx in neighbors[0]:
                ntup = self.tups[idx]
                assert ntup in self.tupMap
                ntot += self.tupMap[ntup][1]
                # Calculate distance 
                dist = numpy.linalg.norm(numpy.subtract(ntup, tup))
                # Calculate a weight that decreases with increasing distance
                maxDist = max(maxDist, dist)
                minDist = min(minDist, dist)
                wtups.append((self.tupMap[ntup][0], self.tupMap[ntup][1], dist))
            if ntot >= self.k:
                break
            radius *= 2.0
        if ntot == 0:
            print >> sys.stderr, "Tuple not within 1000 units of any other tuple: %s" % str(tup)
            return float(self.ncor) / self.ntot
        ncor = sum(map(lambda x: x[0] * (maxDist - x[2]) / (maxDist - minDist), wtups))
        ntot = sum(map(lambda x: x[1] * (maxDist - x[2]) / (maxDist - minDist), wtups))
        return float(ncor) / ntot
    
    def probCorrect(self, tup):
        """ Cacheing wrapper for self.__probCorrect. """
        assert self.kdt is not None
        assert isinstance(tup, collections.Hashable)
        if tup == self.prevTup:
            return self.prevProb
        self.prevTup = tup
        cache = self.cache
        root = self.root
        link = cache.get(tup)
        if link is not None:
            # Cache hit!
            link_prev, link_next, _, result = link
            link_prev[1] = link_next
            link_next[0] = link_prev
            last = root[0]
            last[1] = root[0] = link
            link[0] = last
            link[1] = root
            self.hits += 1
            self.prevProb = result
            return result
        # Cache miss
        result = self.__probCorrect(tup)
        root[2] = tup
        root[3] = result
        oldroot = root
        root = self.root = root[1]
        root[2], oldkey = None, root[2]
        root[3], oldvalue = None, root[3]
        del cache[oldkey]
        cache[tup] = oldroot
        self.misses += 1
        self.prevProb = result
        return result

def go(args):
    if not args.no_mapq:
        # Build KNN classifiers
        st = time.clock()
        training.fit(args.num_neighbors, scaleAs=args.as_scale, scaleDiff=args.diff_scale)
        print >> sys.stderr, "Finished fitting KNN classifiers on %d training tuples" % len(training)
        fitIval = time.clock() - st
        st = time.clock()
        
        mapqDiff = 0.0
        nrecs, npair, nunp = 0, 0, 0
        samTsvOfh = None
        
        exFields = set(["XQ:f"])
        if args.save_sam_tsv:
            # Figure out what all possible extra fields are so that we can write a
            # tsv version of the SAM output with a fixed number of columns per line 
            samIfh.seek(0)
            for samrec in samIfh:
                if samrec[0] == '@':
                    continue
                toks = string.split(samrec, '\t')
                for tok in toks[11:]:
                    idx1 = tok.find(':')
                    assert idx1 != -1
                    idx2 = tok.find(':', idx1+2)
                    assert idx2 != -1
                    exFields.add(tok[:idx2])
            samTsvOfh = open(args.save_sam_tsv, 'w')
            samTsvOfh.write("\t".join([
                "qname", "flag", "rname", "pos", "mapq", "cigar", "rnext",
                "pnext", "tlen", "seq", "qual"] + list(exFields)))
            samTsvOfh.write("\n")
            print >> sys.stderr, "Extra fields: " + str(exFields)
        
        print >> sys.stderr, "Assigning MAPQs"
        with open(args.S, 'w') as samOfh: # open output file
            
            def emitNewSam(samrec, al, probCorrect):
                probIncorrect = 1.0 - probCorrect # convert to probability incorrect
                mapq = args.max_mapq
                if probIncorrect > 0.0:
                    mapq = min(-10.0 * math.log10(probIncorrect), args.max_mapq)
                samrec = samrec.rstrip()
                xqIdx = samrec.find("\tXQ:f:")
                if xqIdx != -1:
                    samrec = samrec[:xqIdx]
                samOfh.write("\t".join([samrec, "XQ:f:" + str(mapq)]) + "\n")
                if samTsvOfh:
                    samToks = string.split(samrec, "\t")
                    samTsvOfh.write("\t".join(samToks[:11]))
                    for field in exFields:
                        samTsvOfh.write("\t")
                        if field == "XQ:f":
                            samTsvOfh.write(str(mapq))
                        else:
                            samTsvOfh.write(al.exFields[field] if field in al.exFields else "NA")
                    samTsvOfh.write("\n")
                return mapq - float(al.mapq)
            
            samIfh.seek(0)                # rewind temporary SAM file
            ival = 100
            while True:
                samrec = samIfh.readline()
                if len(samrec) == 0:
                    break
                if samrec[0] == '@':      # pass headers straight through
                    samOfh.write(samrec)
                    continue
                al = Alignment(samrec)    # parse alignment
                nrecs += 1
                if (nrecs % ival) == 0:
                    print >> sys.stderr, "  Assigned MAPQs to %d records (%d unpaired, %d paired)" % (nrecs, npair, nunp)
                if al.isAligned():
                    # Part of a concordantly-aligned paired-end read?
                    probCorrect = None
                    if al.concordant:
                        # Get the other mate
                        samrec2 = samIfh.readline()
                        assert len(samrec2) > 0
                        al2 = Alignment(samrec2)
                        al1, samrec1 = al, samrec
                        npair += 2
                        if al2.mate1:
                            al1, al2 = al2, al1
                            samrec1, samrec2 = samrec2, samrec1
                        assert al1.mate1 and not al1.mate2
                        assert al2.mate2 and not al2.mate1
                        probCorrect1, probCorrect2 = training.probCorrect(al1, al2)
                        mapqDiff += emitNewSam(samrec1.rstrip(), al1, probCorrect1)
                        mapqDiff += emitNewSam(samrec2.rstrip(), al2, probCorrect2)
                    else:
                        nunp += 1
                        probCorrect = training.probCorrect(al)
                        mapqDiff += emitNewSam(samrec.rstrip(), al, probCorrect)
                else:
                    samOfh.write(samrec)
                    if samTsvOfh:
                        samToks = string.split(samrec, "\t")
                        samTsvOfh.write("\t".join(samToks[:11]))
                        for field in exFields:
                            samTsvOfh.write("\t")
                            samTsvOfh.write(al.exFields[field] if field in al.exFields else "NA")
                        samTsvOfh.write("\n")
        
        if samTsvOfh: samTsvOfh.close()
        print >> sys.stderr, "Finished writing final SAM output (%d records) to '%s'" % (nrecs, args.S)
        print >> sys.stderr, "  concordant-pair SAM records: %d" % npair
        print >> sys.stderr, "  non-concordant-pair SAM records: %d" % nunp
        print >> sys.stderr, "  total mapping quality difference (new - old): %0.3f" % mapqDiff
        print >> sys.stderr, "  total KNN hits=%d, misses=%d" % (training.hits(), training.misses()) 

def makePredictor(args):
    return LogitPredictor()

if __name__ == "__main__":
    import sys
    import unittest
    import argparse
    
    parser = argparse.ArgumentParser(\
        description='Fit model, make predictions.')
    
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
    
    logging.info('Loading training data')
    trainingData = Dataset()
    trainingData.load(args.training)
    trainUnp, trainM, trainPair = trainingData.toDataFrames()
    
    print 'Training:'
    print 'Unp:'
    print trainUnp.head()
    print 'M:'
    print trainM.head()
    print 'Pair:'
    print trainPair.head()
    
    logging.info('Creating and fitting models')
    if not trainUnp.empty:
        logging.info('  unpaired model')
        pred = makePredictor(args)
        pred.train(trainUnp)
    if not trainM.empty:
        logging.info('  non-concordant paired-end model')
        pred = makePredictor(args)
        pred.train(trainM)
    if not trainPair.empty:
        logging.info('  concordant paired-end model')
        pred = makePredictor(args)
        pred.train(trainPair)
    
    logging.info('Loading test data')
    testData = Dataset()
    testData.load(args.test)
    testUnp, testM, testPair = testData.toDataFrames()
    
    print 'Test:'
    print 'Unp:'
    print testUnp.head()
    print 'M:'
    print testM.head()
    print 'Pair:'
    print testPair.head()
    
    class Test(unittest.TestCase):
        def test_rankingError(self):
            result = Result([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], [True, False, False, True, False, True, True, False])
            r = RankingErrorObjective()
            self.assertEqual(1 + 2 + 4 + 7, r.score(result))

    unittest.main(argv=[sys.argv[0]])
