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

import shutil
import os
from read import Alignment
from pandas import DataFrame, read_csv

try:
    import numpypy as np
except ImportError:
    pass


def _str_or_na(n):
    return 'NA' if n is None else str(n)


class UnpairedTuple(object):
    """ Unpaired training/test tuple.  An optional "other end length" field is
        0 for unpaired alignments and >0 for bad-end alignments. """

    def __init__(self, rdname, rdlen, minv, maxv, bestsc, best2sc, mapq, ordlen=0):
        assert minv is None or minv <= bestsc <= maxv
        self.rdname = rdname            # read name
        self.rdlen = rdlen              # read len
        self.ordlen = ordlen            # read len of opposite end
        self.minv = minv                # min valid score
        self.maxv = maxv                # max valid score
        self.bestsc = bestsc            # best
        self.best2sc = best2sc          # 2nd-best score
        self.mapq = mapq                # original mapq

    # header names
    csv_names = ['name', 'best1', 'best2', 'minv', 'maxv', 'rdlen', 'mapq', 'ordlen']
    
    @classmethod
    def from_alignment(cls, al, ordlen=0):
        """ Create unpaired training/test tuple from Alignment object """
        secbest = al.secondBestScore
        if hasattr(al, 'thirdBestScore'):
            secbest = max(secbest, al.thirdBestScore)
        min_valid, max_valid = None, None
        if hasattr(al, 'minValid'):
            assert hasattr(al, 'maxValid')
            min_valid, max_valid = al.minValid, al.maxValid
        return cls(al.name, len(al), min_valid, max_valid, al.bestScore, secbest, al.mapq, ordlen)

    @classmethod
    def append_csv_header(cls, fh):
        fh.write(','.join(cls.csv_names + ['correct']) + '\n')

    def append_csv(self, fh, correct=None):
        correct_str = 'NA'
        if correct is not None:
            correct_str = 'T' if correct else 'F'
        fh.write('%s,%d,%s,%s,%s,%d,%d,%s,%s\n' % (self.rdname, self.bestsc, _str_or_na(self.best2sc),
                                                   _str_or_na(self.minv), _str_or_na(self.maxv),
                                                   self.rdlen, self.mapq, _str_or_na(self.ordlen),
                                                   correct_str))


class PairedTuple(object):
    """ Concordant paired-end training/test tuple.  One per mate alignment. """
    def __init__(self, rdname1, rdlen1, minv1, maxv1,
                 bestsc1, best2sc1, mapq1,
                 rdname2, rdlen2, minv2, maxv2,
                 bestsc2, best2sc2, mapq2,
                 bestconcsc, best2concsc, fraglen):
        assert minv1 is None or minv1 <= bestsc1 <= maxv1
        assert minv2 is None or minv2 <= bestsc2 <= maxv2
        self.rdname1 = rdname1          # read name #1
        self.rdname2 = rdname2          # read name #2
        self.rdlen1 = rdlen1            # read len #1
        self.rdlen2 = rdlen2            # read len #2
        self.minv1 = minv1              # min valid score #1
        self.minv2 = minv2              # min valid score #2
        self.maxv1 = maxv1              # max valid score #1
        self.maxv2 = maxv2              # max valid score #2
        self.bestsc1 = bestsc1          # best #1
        self.bestsc2 = bestsc2          # best #2
        self.best2sc1 = best2sc1        # 2nd-best score #1
        self.best2sc2 = best2sc2        # 2nd-best score #2
        self.mapq1 = mapq1              # original mapq #1
        self.mapq2 = mapq2              # original mapq #2
        self.bestconcsc = bestconcsc    # best concordant
        self.best2concsc = best2concsc  # 2nd-best concordant
        self.fraglen = fraglen          # fragment length

    csv_names = ['name_1', 'best1_1', 'best2_1', 'minv_1', 'maxv_1', 'rdlen_1', 'mapq_1',
                 'name_2', 'best1_2', 'best2_2', 'minv_2', 'maxv_2', 'rdlen_2', 'mapq_2',
                 'best1conc', 'best2conc', 'fraglen']

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
        return cls(al1.name, len(al1), min_valid1, max_valid1, al1.bestScore,
                   secbest1, al1.mapq,
                   al2.name, len(al2), min_valid2, max_valid2, al2.bestScore,
                   secbest2, al2.mapq,
                   best_concordant_score, second_best_concordant_score,
                   Alignment.fragment_length(al1, al2))

    @classmethod
    def append_csv_header(cls, fh):
        fh.write(','.join(cls.csv_names + ['correct']) + '\n')

    def append_csv(self, fh, correct=None):
        correct_str = 'NA'
        if correct is not None:
            correct_str = 'T' if correct else 'F'
        fh.write('%s,%d,%s,%s,%s,%d,%d,%s,%d,%s,%s,%s,%d,%d,%s,%s,%d,%s\n' %
                 (self.rdname1, self.bestsc1, _str_or_na(self.best2sc1),
                  _str_or_na(self.minv1), _str_or_na(self.maxv1), self.rdlen1, self.mapq1,
                  self.rdname2, self.bestsc2, _str_or_na(self.best2sc2),
                  _str_or_na(self.minv2), _str_or_na(self.maxv2), self.rdlen2, self.mapq2,
                  _str_or_na(self.bestconcsc), _str_or_na(self.best2concsc), self.fraglen, correct_str))


class DatasetOnDisk(object):

    """ Encapsulates a collection of training or test data.  Training data is
        labeled, test data not.  Right now this is being stored row-wise.
        This works well for saving and loading the rows to/from a CSV file.
        But it doesn't work so well for other things we'd like to do, like
        rescaling. """

    def __init__(self, name, temp_man):
        # Data for individual reads and mates.  Tuples are (rdlen, minValid,
        # maxValid, bestSc, scDiff)
        self.data_unp, self.data_unp_fn = None, None
        # Data for concordant pairs.  Tuples are two tuples as described above,
        # one for each mate, plus the fragment length.  Label says whether the
        # first mate's alignment is correct.
        self.data_conc, self.data_conc_fn = None, None
        # Data for discordant pairs.
        self.data_disc, self.data_disc_fn = None, None
        # Data for bad ends
        self.data_bad_end, self.data_bad_end_fn = None, None
        self._len = 0  # total # alignments added
        self.name = name  # dataset name
        self.temp_man = temp_man  # temporary file allocator

    def __len__(self):
        """ Return # alignments added so far """
        return self._len

    def add_concordant(self, al1, al2, correct1, correct2):
        """ Add a concordant paired-end alignment to our dataset. """
        assert al1.concordant and al2.concordant
        rec1 = PairedTuple.from_alignments(al1, al2)
        rec2 = PairedTuple.from_alignments(al2, al1)
        if self.data_conc is None:
            self.data_conc_fn = self.temp_man.get_filename('%s_data_conc.csv' % self.name, 'dataset %s' % self.name)
            self.data_conc = open(self.data_conc_fn, 'w')
            PairedTuple.append_csv_header(self.data_conc)
        for rec, correct in zip([rec1, rec2], [correct1, correct2]):
            rec.append_csv(self.data_conc, correct)

    def add_discordant(self, al1, al2, correct1, correct2):
        """ Add a discordant paired-end alignment to our dataset. """
        assert al1.discordant and al2.discordant
        rec1 = PairedTuple.from_alignments(al1, al2)
        rec2 = PairedTuple.from_alignments(al2, al1)
        if self.data_disc is None:
            self.data_disc_fn = self.temp_man.get_filename('%s_data_disc.csv' % self.name, 'dataset %s' % self.name)
            self.data_disc = open(self.data_disc_fn, 'w')
            PairedTuple.append_csv_header(self.data_disc)
        for rec, correct in zip([rec1, rec2], [correct1, correct2]):
            rec.append_csv(self.data_disc, correct)

    def add_bad_end(self, al, unaligned, correct):
        """ Add a discordant paired-end alignment to our dataset. """
        assert al.paired
        if self.data_bad_end is None:
            self.data_bad_end_fn = self.temp_man.get_filename('%s_data_bad_end.csv' % self.name, 'dataset %s' % self.name)
            self.data_bad_end = open(self.data_bad_end_fn, 'w')
            UnpairedTuple.append_csv_header(self.data_bad_end)
        UnpairedTuple.from_alignment(al, len(unaligned.seq)).append_csv(self.data_bad_end, correct)

    def add_unpaired(self, al, correct):
        """ Add an alignment for a simulated unpaired read to our dataset. """
        assert not al.paired
        if self.data_unp is None:
            self.data_unp_fn = self.temp_man.get_filename('%s_data_unp.csv' % self.name, 'dataset %s' % self.name)
            self.data_unp = open(self.data_unp_fn, 'w')
            UnpairedTuple.append_csv_header(self.data_unp)
        UnpairedTuple.from_alignment(al).append_csv(self.data_unp, correct)

    def save(self, fnprefix):
        """ Save a file that we can load from R using read.csv with
            default arguments """
        if self.data_unp_fn is not None:
            shutil.copyfile(self.data_unp_fn, fnprefix + '_unp.csv')
        if self.data_bad_end_fn is not None:
            shutil.copyfile(self.data_bad_end_fn, fnprefix + '_bad_end.csv')
        if self.data_conc_fn is not None:
            shutil.copyfile(self.data_conc_fn, fnprefix + '_conc.csv')
        if self.data_disc_fn is not None:
            shutil.copyfile(self.data_disc_fn, fnprefix + '_disc.csv')

    @staticmethod
    def _from_file(fn, names):
        if fn is not None and os.path.exists(fn):
            return read_csv(fn, names=names)
        else:
            return DataFrame.from_dict(dict.fromkeys(names, []))

    def finalize(self):
        for fh in [self.data_conc, self.data_bad_end, self.data_disc, self.data_unp]:
            if fh is not None:
                fh.close()

    def to_data_frames(self):
        """ Convert dataset to tuple of 3 pandas DataFrames. """
        return (self._from_file(self.data_conc_fn, PairedTuple.csv_names),
                self._from_file(self.data_disc_fn, PairedTuple.csv_names),
                self._from_file(self.data_unp_fn, UnpairedTuple.csv_names),
                self._from_file(self.data_bad_end_fn, UnpairedTuple.csv_names))
