from __future__ import print_function
from __future__ import division
import numpy
import pandas
import os


class MetaMat(object):
    """
    Iterator that returns a large matrix of floats in chunks of rows, where the
    number of rows in a chunk is a parameter passed to the constructor.
    Assumes all elements are double-precision 8-byte floating-point numbers.
    """

    def __init__(self, prefix, chunk_size=1000000):
        """ Parse metadata, check that files exist and initialize members """
        self.prefix = prefix
        self.chunk_size = chunk_size
        self.fh = None
        self.cur = 0
        self.done = False

        meta_fn = prefix + '.meta'
        if not os.path.exists(meta_fn):
            raise RuntimeError('Metadata file does not exist: "%s"' % meta_fn)

        self.data_fn = prefix + '.npy'
        if not os.path.exists(self.data_fn):
            raise RuntimeError('Data file does not exist: "%s"' % self.data_fn)

        with open(meta_fn) as fh:
            fields = fh.readline().rstrip().split(',')
            self.nrow = int(fields[-1])
            self.cols = fields[:-1]

        # Start at first chunk
        self.fh = open(self.data_fn, 'rb')
        self.cur = 0
        self.done = False

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        """ Return next chunk """
        if self.done:
            self.fh.close()
            raise StopIteration
        if self.chunk_size > 0:
            row_i, row_f = self.cur, min(self.cur + self.chunk_size, self.nrow)
        else:
            assert self.cur == 0
            row_i, row_f = 0, self.nrow
        self.done = row_f == self.nrow
        self.cur = row_f
        nelt = (row_f - row_i) * len(self.cols)
        m = numpy.fromfile(self.fh, dtype='float64', count=nelt, sep='')
        assert m.size == nelt, (row_i, row_f, len(self.cols), m.size, nelt)
        m = m.reshape((row_f - row_i, len(self.cols)))
        return pandas.DataFrame(data=m, columns=self.cols)

    def reset(self):
        if self.fh is not None:
            self.fh.close()
        self.fh = open(self.data_fn, 'rb')
        self.cur = 0
        self.done = False

    @staticmethod
    def write_metamat(prefix, col_names, floats=None, append=False):
        """ Caution: this doesn't use numpy tofile, and so is slow for writing
            lots of floats """
        if not append:
            with open(prefix + '.meta', 'wb') as ofh:
                n_col = len(col_names)
                assert len(floats) % n_col == 0
                n_row = len(floats) // n_col
                ofh.write(b','.join(col_names + [str(n_row).encode()]))

        if floats is not None:
            with open(prefix + '.npy', 'ab' if append else 'wb') as ofh:
                for f in floats:
                    ofh.write(struct.pack('d', f))


if __name__ == "__main__":

    import sys
    import unittest
    import struct


    class TestCases(unittest.TestCase):

        def setUp(self):
            """ Create a few matrices to be tested later """
            self.prefixes = ['.testmat_a', '.testmat_b']
            self.float_list = list(map(lambda i: float(i)/1.234534, range(-10000, 10000, 1)))

            # 500 x 2
            MetaMat.write_metamat(self.prefixes[0], [b'alpha', b'bravo'], self.float_list[:1000])

            # 100 x 7
            MetaMat.write_metamat(self.prefixes[1], [b'alpha', b'bravo', b'charlie',
                                                     b'delta', b'echo', b'foxtrot', b'golf'],
                                  self.float_list[:700])

        def test_roc_1a(self):
            m = MetaMat(self.prefixes[0], 7)
            df = next(m)
            self.assertEqual(7, df.shape[0])
            self.assertEqual(2, df.shape[1])
            self.assertAlmostEqual(self.float_list[0], df.alpha[0], places=3)
            self.assertAlmostEqual(self.float_list[1], df.bravo[0], places=3)
            self.assertAlmostEqual(self.float_list[2], df.alpha[1], places=3)
            self.assertAlmostEqual(self.float_list[3], df.bravo[1], places=3)

            for _ in range((500 + 6) // 7 - 2):
                df = next(m)
                self.assertEqual(7, df.shape[0])
                self.assertEqual(2, df.shape[1])

            df = next(m)
            self.assertEqual(500 % 7, df.shape[0])
            self.assertEqual(2, df.shape[1])
            self.assertAlmostEqual(self.float_list[999], df.bravo.iloc[-1], places=3)
            self.assertAlmostEqual(self.float_list[998], df.alpha.iloc[-1], places=3)
            self.assertAlmostEqual(self.float_list[997], df.bravo.iloc[-2], places=3)
            self.assertAlmostEqual(self.float_list[996], df.alpha.iloc[-2], places=3)

            try:
                next(m)
                self.assertTrue(False)
            except StopIteration:
                pass

        def test_roc_1b(self):
            n_row_chunk = 13
            n_col = 7
            m = MetaMat(self.prefixes[1], n_row_chunk)
            df = next(m)
            self.assertEqual(n_row_chunk, df.shape[0])
            self.assertEqual(n_col, df.shape[1])
            for i in range(13):
                self.assertAlmostEqual(self.float_list[0 + i*7], df.alpha[i], places=3)
                self.assertAlmostEqual(self.float_list[1 + i*7], df.bravo[i], places=3)
                self.assertAlmostEqual(self.float_list[2 + i*7], df.charlie[i], places=3)
                self.assertAlmostEqual(self.float_list[3 + i*7], df.delta[i], places=3)
                self.assertAlmostEqual(self.float_list[4 + i*7], df.echo[i], places=3)
                self.assertAlmostEqual(self.float_list[5 + i*7], df.foxtrot[i], places=3)
                self.assertAlmostEqual(self.float_list[6 + i*7], df.golf[i], places=3)

            n_row = 100
            for _ in range((n_row + n_row_chunk - 1) // n_row_chunk - 2):
                df = next(m)
                self.assertEqual(n_row_chunk, df.shape[0])
                self.assertEqual(n_col, df.shape[1])

            df = next(m)
            self.assertEqual(n_row % 13, df.shape[0])
            self.assertEqual(n_col, df.shape[1])

            try:
                next(m)
                self.assertTrue(False)
            except StopIteration:
                pass

        def tearDown(self):
            for prefix in self.prefixes:
                os.remove(prefix + '.meta')
                os.remove(prefix + '.npy')

    unittest.main(argv=[sys.argv[0]])
    sys.exit()
