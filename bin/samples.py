"""
samples.py

Some things that are easy to get from aligners:
1. Score of best alignment (AS:i) for an unpaired read or end
2. Score of second-best alignment (XS:i) for an unpaired read or end
3. Mapping quality (MAPQ) assigned by aligner
4. Implied fragment length for paired-end alignment (TLEN)

Some things that we can get from Bowtie 2 but are harder or impossible
to get from other aligners:
1. Score of second-best concordant alignment for a paired-end read
2. Maximum and minimum possible scores for alignments, for rescaling

This spurs us to the following decisions about our output format:
1. Alignment scores will be unscaled.  Maximum and minimum values will
   be NA by default, but non-NA when available.  Note that maximum and
   minimum are particularly valuable when putting scores on a common
   scale across various read lengths.
2. Second-best concordant alignment score will be NA by default,
   non-NA when available.  This is unfortunate since some models might
   depend on that value and will therefore have to include a special
   case for when it's not available.

Here's the output format.  It's CSV with numeric, integer and logical
data types.  It can be read in using read.csv in R.

For unpaired examples:
1. Score of best alignment (AS:i)
2. Score of second-best alignment (XS:i)
3. Minimum possible "valid" score (or NA if not available)
4. Maximum possible "valid" score (or NA if not available)
5. Length of read
6. Mapping quality (MAPQ)
7. Correct? (only non-NA if data is either training data or simulated
   input data)

For concordantly-aligned paired-end examples:
1. Mate 1: Score of best alignment (AS:i)
2. Mate 1: Score of second-best alignment (XS:i)
3. Mate 1: Minimum possible "valid" score (or NA if not available)
4. Mate 1: Maximum possible "valid" score (or NA if not available)
5. Mate 1: Length of read
6. Mate 1: Mapping quality (MAPQ)
7. Mate 2: Score of best alignment (AS:i)
8. Mate 2: Score of second-best alignment (XS:i)
9. Mate 2: Minimum possible "valid" score (or NA if not available)
10. Mate 2: Maximum possible "valid" score (or NA if not available)
11. Mate 2: Length of read
12. Mate 2: Mapping quality (MAPQ)
13. Score of best concordant paired-end alignment
14. Score of second-best concordant paired-end alignment
15. Fragment length
16. Mate 1: Correct? (only non-NA if data is either training data or
    simulated input data)
"""

import os
import csv

try:
    import numpypy as np
except ImportError:
    pass
import numpy as np
from itertools import izip


class UnpairedTuple(object):
    """ Unpaired training/test tuple. """
    def __init__(self, rdlen, minv, maxv, bestsc, best2sc, mapq):
        assert minv is None or minv <= bestsc <= maxv
        self.rdlen = rdlen              # read len
        self.minv = minv                # min valid score
        self.maxv = maxv                # max valid score
        self.bestsc = bestsc            # best
        self.best2sc = best2sc          # 2nd-best score
        self.mapq = mapq                # original mapq

    def __iter__(self):
        return iter([self.bestsc, self.best2sc, self.minv, self.maxv, self.rdlen, self.mapq])
    
    @classmethod
    def from_alignment(cls, al):
        """ Create unpaired training/test tuple from Alignment object """
        secbest = al.secondBestScore
        if hasattr(al, 'thirdBestScore'):
            secbest = max(secbest, al.thirdBestScore)
        min_valid, max_valid = None, None
        if hasattr(al, 'minValid'):
            assert hasattr(al, 'maxValid')
            min_valid, max_valid = al.minValid, al.maxValid
        return cls(len(al), min_valid, max_valid, al.bestScore, secbest, al.mapq)
    
    @classmethod
    def to_data_frame(cls, ptups, cor=None):
        """ Convert the paired-end tuples to a pandas DataFrame """
        from pandas import DataFrame
        rdlen, best1, best2, minv, maxv, mapq = [], [], [], [], [], []
        for ptup in ptups:
            best1.append(ptup.bestsc)
            best2.append(ptup.best2sc)
            minv.append(ptup.minv)
            maxv.append(ptup.maxv)
            rdlen.append(ptup.rdlen)
            mapq.append(ptup.mapq)
        df = DataFrame.from_items([('best1', best1),
                                   ('best2', best2),
                                   ('minv', minv),
                                   ('maxv', maxv),
                                   ('rdlen', rdlen),
                                   ('mapq', mapq)])
        if cor is not None:
            df['correct'] = np.where(cor, 1, 0)
        return df


class PairedTuple(object):
    """ Concordant paired-end training/test tuple.  One per mate alignment. """
    def __init__(self, rdlen1, minv1, maxv1, bestsc1, best2sc1, mapq1,
                 rdlen2, minv2, maxv2, bestsc2, best2sc2, mapq2,
                 bestconcsc, best2concsc, fraglen):
        assert minv1 is None or minv1 <= bestsc1 <= maxv1
        assert minv2 is None or minv2 <= bestsc2 <= maxv2
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

    def __iter__(self):
        return iter([self.bestsc1, self.best2sc1, self.minv1, self.maxv1, self.rdlen1, self.mapq1,
                     self.bestsc2, self.best2sc2, self.minv2, self.maxv2, self.rdlen2, self.mapq2,
                     self.bestconcsc, self.best2concsc, self.fraglen])
    
    @classmethod
    def from_alignments(cls, al1, al2):
        """ Create unpaired training/test tuple from pair of Alignments """
        secbest1, secbest2 = al1.secondBestScore, al2.secondBestScore
        if hasattr(al1, 'thirdBestScore'):
            assert hasattr(al2, 'thirdBestScore')
            secbest1 = max(secbest1, al1.thirdBestScore)
            secbest2 = max(secbest2, al2.thirdBestScore)
        min_valid1, min_valid2 = None, None
        max_valid1, max_valid2 = None, None
        if hasattr(al1, 'minValid'):
            min_valid1, min_valid2 = al1.minValid, al2.minValid
            max_valid1, max_valid2 = al1.maxValid, al2.maxValid
        best_concordant_score, second_best_concordant_score = None, None
        if hasattr(al1, 'bestConcordantScore'):
            assert hasattr(al2, 'bestConcordantScore')
            assert hasattr(al1, 'secondBestConcordantScore')
            assert hasattr(al2, 'secondBestConcordantScore')
            assert al1.bestConcordantScore == al2.bestConcordantScore
            assert al1.secondBestConcordantScore == al2.secondBestConcordantScore
            best_concordant_score, second_best_concordant_score = \
                al1.bestConcordantScore, al1.secondBestConcordantScore
        return cls(len(al1), min_valid1, max_valid1, al1.bestScore,
                   secbest1, al1.mapq,
                   len(al2), min_valid2, max_valid2, al2.bestScore,
                   secbest2, al2.mapq,
                   best_concordant_score, second_best_concordant_score,
                   al1.fragmentLength())

    @classmethod
    def columnize(cls, ptups):
        rdlen_1, rdlen_2 = [], []
        best1_1, best1_2 = [], []
        best2_1, best2_2 = [], []
        minv_1, maxv_1, minv_2, maxv_2 = [], [], [], []
        mapq_1, mapq_2 = [], []
        best1conc, best2conc = [], []
        fraglen = []
        for ptup in ptups:
            best1_1.append(ptup.bestsc1)
            best1_2.append(ptup.bestsc2)
            best2_1.append(ptup.best2sc1)
            best2_2.append(ptup.best2sc2)
            rdlen_1.append(ptup.rdlen1)
            rdlen_2.append(ptup.rdlen2)
            mapq_1.append(ptup.mapq1)
            mapq_2.append(ptup.mapq2)
            best1conc.append(ptup.bestconcsc)
            best2conc.append(ptup.best2concsc)
            minv_1.append(ptup.minv1)
            minv_2.append(ptup.minv2)
            maxv_1.append(ptup.maxv1)
            maxv_2.append(ptup.maxv2)
            fraglen.append(ptup.fraglen)
        return best1_1, best2_1, minv_1, maxv_1, rdlen_1, mapq_1, best1_2, best2_2, minv_2, maxv_2, rdlen_2, mapq_2,\
            best1conc, best2conc, fraglen

    @classmethod
    def to_data_frames(cls, ptups, cor=None):
        """ Convert the paired-end tuples to a pandas DataFrame """
        from pandas import DataFrame
        best1_1, best2_1, minv_1, maxv_1, rdlen_1, mapq_1, best1_2, best2_2, minv_2, maxv_2, rdlen_2, mapq_2,\
            best1conc, best2conc, fraglen = PairedTuple.columnize(ptups)
        df = DataFrame.from_items([('best1_1', best1_1),
                                   ('best2_1', best2_1),
                                   ('minv_1', minv_1),
                                   ('maxv_1', maxv_1),
                                   ('rdlen_1', rdlen_1),
                                   ('mapq_1', mapq_1),
                                   ('best1_2', best1_2),
                                   ('best2_2', best2_2),
                                   ('minv_2', minv_2),
                                   ('maxv_2', maxv_2),
                                   ('rdlen_2', rdlen_2),
                                   ('mapq_2', mapq_2),
                                   ('best1conc', best1conc),
                                   ('best2conc', best2conc),
                                   ('fraglen', fraglen)])
        if cor is not None:
            df['correct'] = cor
        return df


class Dataset(object):
    
    """ Encapsulates a collection of training or test data.  Training data is
        labeled, test data not.  Right now this is being stored row-wise.
        This works well for saving and loading the rows to/from a CSV file.
        But it doesn't work so well for other things we'd like to do, like
        rescaling. """
    
    def __init__(self):
        # Data for individual reads and mates.  Tuples are (rdlen, minValid,
        # maxValid, bestSc, scDiff)
        self.data_unp, self.lab_unp = [], []
        self.data_m, self.lab_m = [], []
        # Data for concordant pairs.  Tuples are two tuples as described above,
        # one for each mate, plus the fragment length.  Label says whether the
        # first mate's alignment is correct.
        self.data_conc, self.lab_conc = [], []

    def __len__(self):
        """ Return number of alignments added so far """
        return len(self.data_unp) + len(self.data_m) + len(self.data_conc)
    
    def add_paired(self, al1, al2, correct1, correct2):
        """ Add a concordant paired-end alignment to our dataset. """
        assert al1.concordant and al2.concordant
        rec1 = PairedTuple.from_alignments(al1, al2)
        rec2 = PairedTuple.from_alignments(al2, al1)
        for rec in [rec1, rec2]:
            self.data_conc.append(rec)
        self.lab_conc.extend([correct1, correct2])
    
    def add_unpaired(self, al, correct):
        """ Add an alignment for a simulated unpaired read to our dataset. """
        self.data_unp.append(UnpairedTuple.from_alignment(al))
        self.lab_unp.append(correct)
    
    def add_mate(self, al, correct):
        """ Add an alignment for a simulated mate that has been aligned in an
            unpaired fashion to our dataset. """
        self.data_m.append(UnpairedTuple.from_alignment(al))
        self.lab_m.append(correct)

    def save(self, fnprefix, compress=True):
        """ Save a file that we can load from R using read.csv with
            default arguments """
        fnmap = {'_unp': (self.data_unp, self.lab_unp),
                 '_mates': (self.data_m, self.lab_m),
                 '_conc': (self.data_conc, self.lab_conc)}
        for lab, p in fnmap.iteritems():
            data, corrects = p
            fn = fnprefix + lab + '.csv'
            if compress:
                import gzip
                fh = gzip.open(fn + '.gz', 'w')
            else:
                fh = open(fn, 'w')
            if lab == '_conc':
                fh.write(','.join(['best1', 'secbest1', 'minv1', 'maxv1', 'len1', 'mapq1',
                                   'best2', 'secbest2', 'minv2', 'maxv2', 'len2', 'mapq2',
                                   'bestconc', 'secbestconc', 'fraglen', 'correct']) + '\n')
            else:
                fh.write(','.join(['best', 'secbest', 'minv', 'maxv', 'len', 'mapq',
                                   'correct']) + '\n')
            for tup, correct in izip(data, corrects):
                correct_str = 'NA'
                if correct is not None:
                    correct_str = 'T' if correct else 'F'
                tup = map(lambda x: 'NA' if x is None else str(x), tup)
                tup.append(correct_str)
                fh.write(','.join(tup) + '\n')
            fh.close()

    def load(self, fnprefix):
        fnmap = {'_unp': (self.data_unp, self.lab_unp),
                 '_mates': (self.data_m, self.lab_m),
                 '_conc': (self.data_conc, self.lab_conc)}
        for lab, p in fnmap.iteritems():
            data, corrects = p
            fn = fnprefix + lab + '.csv'
            if os.path.exists(fn + '.gz'):
                import gzip
                fh = gzip.open(fn + '.gz')
            else:
                fh = open(fn)

            def int_or_none(s):
                return None if s == 'NA' else int(s)

            if lab == '_conc':
                for toks in csv.reader(fh):
                    assert 16 == len(toks)
                    if 'best1' == toks[0]:
                        continue  # skip header
                    # Note: pandas csv parser is much faster
                    best1, secbest1, minv1, maxv1, ln1, mapq1 = map(int_or_none, toks[0:6])
                    best2, secbest2, minv2, maxv2, ln2, mapq2 = map(int_or_none, toks[6:12])
                    bestconc, secbestconc, fraglen = map(int_or_none, toks[12:15])
                    data.append(PairedTuple(ln1, minv1, maxv1, best1, secbest1, mapq1,
                                            ln2, minv2, maxv2, best2, secbest2, mapq2,
                                            bestconc, secbestconc, fraglen))
                    corrects.append(toks[-1] == 'T')
            else:
                for toks in csv.reader(fh):
                    assert 7 == len(toks)
                    if 'best' == toks[0]:
                        continue  # skip header
                    # Note: pandas csv parser is much faster
                    best, secbest, minv, maxv, ln, mapq = map(int_or_none, toks[:6])
                    data.append(UnpairedTuple(ln, minv, maxv, best, secbest, mapq))
                    corrects.append(toks[-1] == 'T')
            fh.close()

    def to_data_frames(self):
        """ Convert dataset to tuple of 3 pandas DataFrames. """
        return (UnpairedTuple.to_data_frame(self.data_unp, self.lab_unp),
                UnpairedTuple.to_data_frame(self.data_m, self.lab_m),
                PairedTuple.to_data_frames(self.data_conc, self.lab_conc))
