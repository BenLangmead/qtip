import pandas as pd
import numpy as np

class UnpairedTuple(object):
    ''' Unpaired training/test tuple. '''
    def __init__(self, rdlen, minv, maxv, bestsc, best2sc, mapq):
        assert maxv >= minv
        assert bestsc >= minv and bestsc <= maxv, str([bestsc, minv, maxv])
        self.rdlen = rdlen              # read len
        self.minv = minv                # min valid score
        self.maxv = maxv                # max valid score
        self.bestsc = bestsc            # best
        self.best2sc = best2sc          # 2nd-best score
        self.mapq = mapq                # original mapq
        self.rescale()
    
    def rescale(self):
        ''' Rescale scores to be in [0, 1] where min/max are minv/maxv.  Put
            results in fields starting with "r". '''
        self.rbestsc = float(self.bestsc - self.minv) / (self.maxv - self.minv)
        self.rbest2sc = None
        if self.best2sc is not None:
            self.rbest2sc = float(self.best2sc - self.minv) / (self.maxv - self.minv)
    
    @classmethod
    def fromAlignment(cls, al):
        ''' Create unpaired training/test tuple from Alignment object '''
        secbest = max(al.secondBestScore, al.thirdBestScore)
        return cls(len(al), al.minValid, al.maxValid, al.bestScore, secbest, al.mapq)
    
    @classmethod
    def toDataFrame(cls, ptups, cor=None):
        ''' Convert the paired-end tuples to a pandas DataFrame '''
        rdlen, best1, best2, mapq = [], [], [], []
        for ptup in ptups:
            rdlen.append(ptup.rdlen)
            best1.append(ptup.rbestsc)
            best2.append(ptup.rbest2sc)
            mapq.append(ptup.mapq)
        df = pd.DataFrame.from_items([\
            ('rdlen', rdlen),
            ('best1', best1),
            ('best2', best2),
            ('mapq', mapq)])
        if cor is not None:
            df['correct'] = np.where(cor, 1, 0)
        return df

class PairedTuple(object):
    ''' Concordant paired-end training/test tuple.  One per mate alignment. '''
    def __init__(self, rdlen1, minv1, maxv1, bestsc1, best2sc1, mapq1,
                 rdlen2, minv2, maxv2, bestsc2, best2sc2, mapq2,
                 bestconcsc, best2concsc, fraglen):
        assert maxv1 >= minv1
        assert maxv2 >= minv2
        assert bestsc1 >= minv1 and bestsc1 <= maxv1
        assert bestsc2 >= minv2 and bestsc2 <= maxv2
        self.rdlen1 = rdlen1            # read len
        self.rdlen2 = rdlen2            # read len
        self.minv1 = minv1              # min valid score
        self.minv2 = minv2              # min valid score
        self.maxv1 = maxv1              # max valid score
        self.maxv2 = maxv2              # max valid score
        self.bestsc1 = bestsc1          # best
        self.bestsc2 = bestsc2          # best
        self.best2sc1 = best2sc1        # 2nd-best score
        self.best2sc2 = best2sc2        # 2nd-best score
        self.mapq1 = mapq1              # original mapq
        self.mapq2 = mapq2              # original mapq
        self.bestconcsc = bestconcsc    # best concordant
        self.best2concsc = best2concsc  # 2nd-best concordant
        self.fraglen = fraglen          # fragment length
        self.rescale()
    
    def rescale(self):
        ''' Rescale scores to be in [0, 1] where min/max are minv/maxv.  Put
            results in fields starting with "r". '''
        self.rbestsc1 = float(self.bestsc1 - self.minv) / (self.maxv - self.minv)
        self.rbest2sc1 = None
        if self.best2sc1 is not None:
            self.rbest2sc1 = float(self.best2sc1 - self.minv) / (self.maxv - self.minv)
        self.rbestsc2 = float(self.bestsc2 - self.minv) / (self.maxv - self.minv)
        self.rbest2sc2 = None
        if self.best2sc2 is not None:
            self.rbest2sc2 = float(self.best2sc2 - self.minv) / (self.maxv - self.minv)
    
    @classmethod
    def fromAlignments(cls, al1, al2, fraglen):
        ''' Create unpaired training/test tuple from pair of Alignments '''
        secbest1 = max(al1.secondBestScore, al1.thirdBestScore)
        secbest2 = max(al2.secondBestScore, al2.thirdBestScore)
        return cls(len(al1), al1.minValid, al1.maxValid, al1.bestScore,
                   secbest1, al1.mapq,
                   len(al2), al2.minValid, al2.maxValid, al2.bestScore,
                   secbest2, al2.mapq,
                   al1.bestConcordantScore,
                   al1.secondBestConcordantScore, al1.fragmentLength())
    
    @classmethod
    def toDataFrame(cls, ptups, cor=None):
        ''' Convert the paired-end tuples to a pandas DataFrame '''
        rdlen_1, rdlen_2 = [], []
        best1_1, best1_2 = [], []
        best2_1, best2_2 = [], []
        mapq_1, mapq_2 = [], []
        best1conc, best2conc = [], []
        fraglen = []
        for ptup in ptups:
            rdlen_1.append(ptup.rdlen1)
            rdlen_2.append(ptup.rdlen2)
            best1_1.append(ptup.rbestsc1)
            best1_2.append(ptup.rbestsc2)
            best2_1.append(ptup.rbest2sc1)
            best2_2.append(ptup.rbest2sc2)
            mapq_1.append(ptup.mapq1)
            mapq_2.append(ptup.mapq2)
            best1conc.append(ptup.bestconcsc)
            best2conc.append(ptup.best2concsc)
            fraglen.append(ptup.fraglen)
        df = pd.DataFrame.from_items([\
            ('rdlen_1', rdlen_1),
            ('rdlen_2', rdlen_2),
            ('best1_1', best1_1),
            ('best1_2', best1_2),
            ('best2_1', best2_1),
            ('best2_2', best2_2),
            ('mapq_1', mapq_1),
            ('mapq_2', mapq_2),
            ('best1conc', best1conc),
            ('best2conc', best2conc),
            ('fraglen', fraglen)])
        if cor is not None:
            df['correct'] = cor
        return df

class Dataset(object):
    
    ''' Encapsulates a collection of training or test data.  Training data is
        labeled, test data not. '''
    
    def __init__(self):
        # Data for individual reads and mates.  Tuples are (rdlen, minValid,
        # maxValid, bestSc, scDiff)
        self.dataUnp, self.labUnp = [], []
        self.dataM,   self.labM   = [], []
        # Data for concordant pairs.  Tuples are two tuples as described above,
        # one for each mate, plus the fragment length.  Label says whether the
        # first mate's alignment is correct.
        self.dataConc, self.labConc = [], []
        self.isTraining = False
    
    def __len__(self):
        ''' Return number of alignments added so far '''
        return len(self.dataUnp) + len(self.dataM) + len(self.dataConc)
    
    def addPaired(self, al1, al2, fraglen, correct1, correct2):
        ''' Add a concordant paired-end alignment to our dataset. '''
        if correct1 is not None or correct2 is not None:
            self.isTraining = True
        assert al1.concordant and al2.concordant
        rec1 = PairedTuple.fromAlignments(al1, al2, fraglen)
        rec2 = PairedTuple.fromAlignments(al2, al1, fraglen)
        for rec in [rec1, rec2]: self.dataConc.append(rec)
        self.labConc.extend([correct1, correct2])
    
    def addUnp(self, al, correct):
        ''' Add an alignment for a simulated unpaired read to our dataset. '''
        self.dataUnp.append(UnpairedTuple.fromAlignment(al))
        self.labUnp.append(correct)
    
    def addM(self, al, correct):
        ''' Add an alignment for a simulated mate that has been aligned in an
            unpaired fashion to our dataset. '''
        self.dataM.append(UnpairedTuple.fromAlignment(al))
        self.labM.append(correct)
    
    def save(self, fn, compress=True):
        ''' Save dataset to a (possibly) compressed pickle file. '''
        import cPickle
        save = (self.dataUnp,  self.labUnp,
                self.dataM,    self.labM,
                self.dataConc, self.labConc)
        if compress:
            import gzip
            fh = gzip.open(fn, 'wb')
        else:
            fh = open(fn, 'wb')
        cPickle.dump(save, fh, cPickle.HIGHEST_PROTOCOL)
        fh.close()
    
    def load(self, fn, compress=True):
        ''' Load dataset from a (possibly) compressed pickle file. '''
        import cPickle
        if compress:
            import gzip
            fh = gzip.open(fn, 'rb')
        else:
            fh = open(fn, 'rb')
        (self.dataUnp,  self.labUnp, \
         self.dataM,    self.labM, \
         self.dataConc, self.labConc) = cPickle.load(fh)
        fh.close()
    
    def toDataFrames(self):
        ''' Convert dataset to tuple of 3 pandas DataFrames. '''
        return (UnpairedTuple.toDataFrame(self.dataUnp,  self.labUnp),
                UnpairedTuple.toDataFrame(self.dataM,    self.labM),
                PairedTuple.toDataFrame(self.dataConc, self.labConc))
