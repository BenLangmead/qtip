__author__ = 'langmead'

import os
import math
from collections import defaultdict
import csv


class AlignabilityStratifiedIndexed():
    """
    Like AlignabilityStratified but uses an index of the .wig input
    file(s) so that we don't have to load the entire .wig into memory.
    The index has per-wig file portions (ending in .index, .freq) and
    per-chromosome portions (ending in .csv).

    A final idea for speeding this up is to use pandas for csv file
    reading/writing.  the csv module and my custom loader are slow.
    """

    def _score_to_stratum(self, score):
        return min(int(math.ceil(abs(-math.log(score, 2)))), self.max_strata-1)

    def _ref_id_to_filename(self, ala_fn, ref_id):
        assert ref_id is not None
        dir_name = os.path.dirname(ala_fn)
        if dir_name != '':
            dir_name += '/'
        fn = '%s.%s.max_strata%d.csv' % (dir_name, ref_id, self.max_strata)
        self.ref_id_to_csv[ref_id] = fn
        return fn

    def _wig_fn_to_index_fn(self, wig_fn):
        dir_name = os.path.dirname(wig_fn)
        if dir_name != '':
            dir_name += '/'
        return '%s.%s.max_strata%d.index' % (dir_name, os.path.basename(wig_fn), self.max_strata)

    def _wig_fn_to_frequency_fn(self, wig_fn):
        dir_name = os.path.dirname(wig_fn)
        if dir_name != '':
            dir_name += '/'
        return '%s.%s.max_strata%d.freq' % (dir_name, os.path.basename(wig_fn), self.max_strata)

    @staticmethod
    def _dump_recs(recs, fn):
        with open(fn, 'wb') as ofh:
            writer = csv.writer(ofh)
            writer.writerows(recs)

    @staticmethod
    def _load_recs(fn):
        recs = []
        with open(fn, 'rb') as ofh:
            for ln in ofh:
                toks = ln.split(',')
                st, en, stratum = int(toks[0]), int(toks[1]), int(toks[2])
                recs.append([st, en, stratum])
        return recs

    def _index_wig(self, fn):
        # check if an index already exists for this wig file
        idx_fn, freq_fn = self._wig_fn_to_index_fn(fn), self._wig_fn_to_frequency_fn(fn)
        if os.path.exists(idx_fn):
            # load per-chromosome index information
            with open(idx_fn) as fh:
                for ln in fh:
                    ref_id, csv_fn = ln.split()
                    assert ref_id not in self.ref_id_to_csv
                    dir_name = os.path.dirname(fn)
                    if dir_name != '' and not dir_name.endswith('/'):
                        dir_name += '/'
                    self.ref_id_to_csv[ref_id] = dir_name + csv_fn
                    if not os.path.exists(self.ref_id_to_csv[ref_id]):
                        raise RuntimeError('No such per-chromosome index file: "%s"' % self.ref_id_to_csv[ref_id])
            # load stratum frequency information
            with open(freq_fn) as fh:
                for i, freq in enumerate(map(int, fh.read().split())):
                    self.frequencies[i] += freq
            return  # skip index building
        # index doesn't already exist; build it now
        index, frequencies = [], [0] * self.max_strata
        with open(fn, 'rb') as fh:
            last_ref_id, recs = None, []
            for ln in fh:
                if ln[0] == '#':
                    continue
                toks = ln.split()
                assert len(toks) >= 4
                ref_id, st, en, score = toks[0], int(toks[1]), int(toks[2]), float(toks[3])
                if last_ref_id is None or last_ref_id != ref_id:
                    if last_ref_id is not None:
                        # write index for this chromosome
                        csv_fn = self._ref_id_to_filename(fn, last_ref_id)
                        self._dump_recs(recs, csv_fn)
                        index.append((last_ref_id, os.path.basename(csv_fn)))
                        recs = []
                    assert ref_id not in self.ref_id_to_csv
                    last_ref_id = ref_id
                    assert len(recs) == 0
                stratum = self._score_to_stratum(score)
                if len(recs) > 0 and recs[-1][2] == stratum:
                    recs[-1][1] = en  # extend existing record
                else:
                    recs.append([st, en, stratum])
                frequencies[stratum] += (en - st)
            if len(recs) > 0:
                # write index for this chromosome
                csv_fn = self._ref_id_to_filename(fn, last_ref_id)
                self._dump_recs(recs, csv_fn)
                index.append((last_ref_id, os.path.basename(csv_fn)))
        # write index for this wig file
        with open(idx_fn, 'w') as ofh:
            for k, v in index:
                ofh.write('%s\t%s\n' % (k, v))
        # write stratum frequency information for this file
        with open(freq_fn, 'w') as ofh:
            ofh.write('\t'.join(map(str, [frequencies[i] for i in xrange(self.max_strata)])))
        # add stratum frequency information to tally
        for i in xrange(self.max_strata):
            self.frequencies[i] += frequencies[i]

    def __init__(self, wig_filenames, max_strata=6):
        assert max_strata > 0
        self.max_strata = max_strata
        self.frequencies = defaultdict(int)
        self.ref_id_to_csv = {}
        self.ref_id = None
        self.records = []
        for wig_fn in wig_filenames:
            self._index_wig(wig_fn)

    def has_ref_id(self, ref_id):
        return ref_id in self.ref_id_to_csv

    def current_ref_id(self):
        return self.ref_id

    def load_ref_id(self, ref_id):
        """ Return iterator over alignability records for given
            reference id. """
        if ref_id not in self.ref_id_to_csv:
            raise RuntimeError('Reference id "%s" not mentioned in alignability files' % ref_id)
        self.records = self._load_recs(self.ref_id_to_csv[ref_id])

if __name__ == "__main__":

    import sys
    import unittest
    import argparse
    from tempfile import mkdtemp
    from shutil import rmtree

    parser = argparse.ArgumentParser()
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')

    args = parser.parse_args()

    if args.test:
        import unittest

        class TestAlignabilityStratified(unittest.TestCase):

            def setUp(self):
                self.tmpdir = mkdtemp()
                self.wig_fn_1 = os.path.join(self.tmpdir, 'tmp1.wig')
                with open(self.wig_fn_1, 'w') as fh:
                    fh.write('''#bedGraph section chr1:10000-81873
chr1    10000   10014   0.00277778
chr1    10014   10015   0.333333
chr1    10015   10026   0.5
chr1    10026   10031   1
chr1    10031   10036   0.5
chr1    10036   10037   0.333333
chr1    10037   10048   0.2
chr1    10048   10050   0.166667
chr1    10050   10051   0.333333
chr1    10051   10061   1
chr1    10061   10064   0.5
chr1    10064   10533   1
chr1    10533   10603   0.5
chr1    10603   10616   0.25
chr1    10616   10624   0.333333
chr2    20000   20014   0.00277778
chr2    20014   20015   0.00277778
''')
            # [5, 2, 1, 0, 1, 2, 3, 3, 2, 0, 1, 0, 1, 2, 2]

            def tearDown(self):
                rmtree(self.tmpdir)

            def test_iter_records_1(self):
                al2 = AlignabilityStratifiedIndexed([self.wig_fn_1])
                al2.load_ref_id('chr1')
                al2_runs = [x for x in al2.records]
                self.assertEqual(13, len(al2_runs))

                # check strata
                self.assertEqual(5, al2_runs[0][2])
                self.assertEqual(2, al2_runs[1][2])
                self.assertEqual(1, al2_runs[2][2])
                self.assertEqual(0, al2_runs[3][2])
                self.assertEqual(1, al2_runs[4][2])
                self.assertEqual(2, al2_runs[5][2])
                self.assertEqual(3, al2_runs[6][2])
                self.assertEqual(2, al2_runs[7][2])
                self.assertEqual(0, al2_runs[8][2])
                self.assertEqual(1, al2_runs[9][2])
                self.assertEqual(0, al2_runs[10][2])
                self.assertEqual(1, al2_runs[11][2])
                self.assertEqual(2, al2_runs[12][2])

                # check starting positions
                self.assertEqual(10000, al2_runs[0][0])
                self.assertEqual(10014, al2_runs[1][0])
                self.assertEqual(10015, al2_runs[2][0])
                self.assertEqual(10026, al2_runs[3][0])
                self.assertEqual(10031, al2_runs[4][0])
                self.assertEqual(10036, al2_runs[5][0])
                self.assertEqual(10037, al2_runs[6][0])
                self.assertEqual(10050, al2_runs[7][0])
                self.assertEqual(10051, al2_runs[8][0])
                self.assertEqual(10061, al2_runs[9][0])
                self.assertEqual(10064, al2_runs[10][0])
                self.assertEqual(10533, al2_runs[11][0])
                self.assertEqual(10603, al2_runs[12][0])

                al2.load_ref_id('chr2')
                al2_runs_chr2 = [x for x in al2.records]
                self.assertEqual(1, len(al2_runs_chr2))
                self.assertEqual(5, al2_runs_chr2[0][2])
                self.assertEqual(20000, al2_runs_chr2[0][0])

            def test_iter_chunks_1(self):
                if False:
                    al2 = AlignabilityStratifiedIndexed([self.wig_fn_1])
                    al2_chunks = [x for x in al2.iter_chunks('chr1', chunk_size=100)]

                    for i, chunk in enumerate(al2_chunks):
                        self.assertEqual(100 * i, chunk[0])
                        if i < 100:
                            self.assertEqual(100 * [None], chunk[1])

                    self.assertEqual(100, len(al2_chunks[100][1]))
                    self.assertEqual(([5] * 14) +
                                     ([2] * 1) +
                                     ([1] * 11) +
                                     ([0] * 5) +
                                     ([1] * 5) +
                                     ([2] * 1) +
                                     ([3] * 13) +
                                     ([2] * 1) +
                                     ([0] * 10) +
                                     ([1] * 3) +
                                     ([0] * 36), al2_chunks[100][1])
                    self.assertEqual(100, len(al2_chunks[101][1]))
                    self.assertEqual([0] * 100, al2_chunks[101][1])
                    self.assertEqual([0] * 100, al2_chunks[102][1])
                    self.assertEqual([0] * 100, al2_chunks[103][1])
                    self.assertEqual([0] * 100, al2_chunks[104][1])
                    self.assertEqual(([0] * 33) + ([1] * 67), al2_chunks[105][1])
                    self.assertEqual(([1] * 3) + ([2] * 21), al2_chunks[106][1])

            def test_calculate_stratum_frequencies_1(self):
                for i in xrange(2):
                    al2 = AlignabilityStratifiedIndexed([self.wig_fn_1])
                    self.assertEqual(5 + 10 + 36 + 400 + 33, al2.frequencies[0])
                    self.assertEqual(11 + 5 + 3 + 67 + 3, al2.frequencies[1])
                    self.assertEqual(1 + 1 + 1 + 21, al2.frequencies[2])
                    self.assertEqual(13, al2.frequencies[3])
                    self.assertEqual(0, al2.frequencies[4])
                    self.assertEqual(29, al2.frequencies[5])

        unittest.main(argv=[sys.argv[0]])
