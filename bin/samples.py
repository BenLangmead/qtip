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

try:
    import numpypy as np
except ImportError:
    pass


def _str_or_na(n):
    return 'NA' if n is None else str(n)


class UnpairedTuple(object):
    """ Unpaired training/test tuple.  An optional "other end length" field is
        0 for unpaired alignments and >0 for bad-end alignments. """

    def __init__(self, rdname, rdlen, mapq, ztzs, ordlen=0):
        self.rdname = rdname            # read name
        self.mapq = mapq                # original mapq
        self.first = True

        # features for learning
        self.rdlen = rdlen              # read len
        self.ordlen = ordlen            # read len of opposite end
        self.ztzs = ztzs or []

    # header names
    csv_names = ['name', 'rdlen', 'mapq', 'ordlen']
    
    @classmethod
    def from_alignment(cls, al, ordlen=0):
        """ Create unpaired training/test tuple from Alignment object """
        return cls(al.name, len(al), al.mapq, al.ztzs, ordlen)

    @classmethod
    def append_csv_header(cls, fh, num_ztzs):
        ztz_colnames = ['ztz%d' % i for i in xrange(num_ztzs)]
        fh.write(','.join(cls.csv_names + ztz_colnames + ['correct']) + '\n')

    def append_csv(self, fh, correct=None):
        correct_str = 'NA'
        if correct is not None:
            correct_str = 'T' if correct else 'F'
        ls = [self.rdname, self.rdlen, self.mapq, _str_or_na(self.ordlen)] + self.ztzs + [correct_str]
        ls = map(str, ls)
        fh.write(','.join(ls) + '\n')


class PairedTuple(object):
    """ Concordant paired-end training/test tuple.  One per mate alignment. """
    def __init__(self, rdname1, rdlen1, mapq1, rdname2, rdlen2, mapq2, al1ztzs, fraglen):
        self.rdname1 = rdname1          # read name #1
        self.rdname2 = rdname2          # read name #2
        self.rdlen1 = rdlen1            # read len #1
        self.rdlen2 = rdlen2            # read len #2
        self.mapq1 = mapq1              # original mapq #1
        self.mapq2 = mapq2              # original mapq #2
        self.fraglen = fraglen          # fragment length
        self.ztzs1 = al1ztzs or []

    csv_names = ['name_1', 'rdlen_1', 'mapq_1', 'name_2', 'rdlen_2', 'mapq_2', 'fraglen']

    @classmethod
    def from_alignments(cls, al1, al2):
        """ Create unpaired training/test tuple from pair of Alignments """
        return cls(al1.name, len(al1), al1.mapq,
                   al2.name, len(al2), al2.mapq,
                   al1.ztzs, Alignment.fragment_length(al1, al2))

    @classmethod
    def append_csv_header(cls, fh, num_ztzs):
        ztz_colnames = ['ztz%d' % i for i in xrange(num_ztzs)]
        fh.write(','.join(cls.csv_names + ztz_colnames + ['correct']) + '\n')

    def append_csv(self, fh, correct=None):
        correct_str = 'NA'
        if correct is not None:
            correct_str = 'T' if correct else 'F'
        ls = [self.rdname1,
              self.rdlen1, self.mapq1, self.rdname2,
              self.rdlen2, self.mapq2, self.fraglen] + self.ztzs1 + [correct_str]
        ls = map(str, ls)
        fh.write(','.join(ls) + '\n')


class DatasetOnDisk(object):

    """ Encapsulates a collection of training or test data.  Training data is
        labeled, test data not.  Right now this is being stored row-wise.
        This works well for saving and loading the rows to/from a CSV file.
        But it doesn't work so well for other things we'd like to do, like
        rescaling. """

    def __init__(self, name, temp_man):
        # Data for individual reads and mates.
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
        self._len += 2
        if self.data_conc is None:
            self.data_conc_fn = self.temp_man.get_file('%s_conc.csv' % self.name, 'dataset %s' % self.name)
            self.data_conc = open(self.data_conc_fn, 'w')
            assert len(al1.ztzs or []) == len(al2.ztzs or [])
            PairedTuple.append_csv_header(self.data_conc, len(al1.ztzs or []))
        for rec, correct in zip([rec1, rec2], [correct1, correct2]):
            rec.append_csv(self.data_conc, correct)

    def add_discordant(self, al1, al2, correct1, correct2):
        """ Add a discordant paired-end alignment to our dataset. """
        assert al1.discordant and al2.discordant
        rec1 = PairedTuple.from_alignments(al1, al2)
        rec2 = PairedTuple.from_alignments(al2, al1)
        self._len += 2
        if self.data_disc is None:
            self.data_disc_fn = self.temp_man.get_file('%s_disc.csv' % self.name, 'dataset %s' % self.name)
            self.data_disc = open(self.data_disc_fn, 'w')
            assert len(al1.ztzs or []) == len(al2.ztzs or [])
            PairedTuple.append_csv_header(self.data_disc, len(al1.ztzs or []))
        for rec, correct in zip([rec1, rec2], [correct1, correct2]):
            rec.append_csv(self.data_disc, correct)

    def add_bad_end(self, al, unaligned, correct):
        """ Add a discordant paired-end alignment to our dataset. """
        assert al.paired
        self._len += 1
        if self.data_bad_end is None:
            self.data_bad_end_fn = self.temp_man.get_file('%s_bad_end.csv' % self.name, 'dataset %s' % self.name)
            self.data_bad_end = open(self.data_bad_end_fn, 'w')
            UnpairedTuple.append_csv_header(self.data_bad_end, len(al.ztzs or []))
        UnpairedTuple.from_alignment(al, len(unaligned.seq)).append_csv(self.data_bad_end, correct)

    def add_unpaired(self, al, correct):
        """ Add an alignment for a simulated unpaired read to our dataset. """
        assert not al.paired
        self._len += 1
        if self.data_unp is None:
            self.data_unp_fn = self.temp_man.get_file('%s_unp.csv' % self.name, 'dataset %s' % self.name)
            self.data_unp = open(self.data_unp_fn, 'w')
            UnpairedTuple.append_csv_header(self.data_unp, len(al.ztzs or []))
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
        import pandas
        if fn is not None and os.path.exists(fn):
            return pandas.read_csv(fn, names=names)
        else:
            return pandas.DataFrame.from_dict(dict.fromkeys(names, []))

    def finalize(self):
        for fh in [self.data_conc, self.data_bad_end, self.data_disc, self.data_unp]:
            if fh is not None:
                fh.close()

    def purge(self):
        self.temp_man.remove_group('dataset %s' % self.name)

    def to_data_frames(self):
        """ Convert dataset to tuple of 3 pandas DataFrames. """
        return (self._from_file(self.data_conc_fn, PairedTuple.csv_names),
                self._from_file(self.data_disc_fn, PairedTuple.csv_names),
                self._from_file(self.data_unp_fn, UnpairedTuple.csv_names),
                self._from_file(self.data_bad_end_fn, UnpairedTuple.csv_names))
