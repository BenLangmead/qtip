import os
import pandas
import warnings
import math
import numpy as np
try:
    from itertools import imap
except ImportError:
    imap = map


__author__ = 'langmead'


def _exists_and_nonempty(fn):
    return os.path.exists(fn) and os.stat(fn).st_size > 0


class FeatureTableReader(object):

    """ Reads a table of information describing alignments.  These are tables
        output by qsim.  Tables might describe training or test alignments.

        The trickiest issue here is whether and how to normalize certain
        features, like alignment scores.  The goal would be to make scores
        more comparable across reads of different lengths.  A read of length
        1000 with alignment score 80 should probably not be considered "close
        to" a read of length 100 with alignment score 80.

        If we know the minimum and maximum possible alignment score (which
        depend on the alignment tool's scoring function, which might in turn
        depend on read length), we can standardize to that interval.  We can't
        generally expect to have that information, though.  Even if we had it,
        we might have to pass a different interval for each read, since they
        can vary in length.

        To avoid these issues, we simply add read length as a feature.  This
        has no effect for uniform-length data, since we also remove zero-
        variance features.
    """

    #            short name   suffix
    datasets = [('u',         '_rec_u.csv'),
                ('d',         '_rec_d.csv'),
                ('c',         '_rec_c.csv'),
                ('b',         '_rec_b.csv')]

    def __init__(self, prefix, chunksize=50000):
        self.prefix = prefix
        self.dfs = {}
        self.readers = {}
        nonempty = False
        fns = []
        for sn, suf in self.datasets:
            fn = self.prefix + suf
            fns.append(fn)
            for fnex in [fn, fn + '.gz', fn + '.bz2']:
                if _exists_and_nonempty(fnex):
                    nonempty = True

                    def gen_new_iterator(_fn, _chunksize):
                        def _new_iterator():
                            if os.path.exists(_fn):
                                return pandas.io.parsers.read_csv(_fn, quoting=2, chunksize=_chunksize)
                            elif os.path.exists(_fn + '.gz'):
                                return pandas.io.parsers.read_csv(_fn + '.gz', quoting=2, chunksize=_chunksize, compression='gzip')
                            elif os.path.exists(_fn + '.bz2'):
                                return pandas.io.parsers.read_csv(_fn + '.bz2', quoting=2, chunksize=_chunksize, compression='bz2')
                            else:
                                raise RuntimeError('No such file: "%s"' % _fn)

                        return _new_iterator

                    self.readers[sn] = gen_new_iterator(fn, chunksize)
                    break

        if not nonempty:
            raise RuntimeError('No non-empty input files with names like: ' + str(fns))


    @staticmethod
    def _postprocess_data_frame(df):
        """ Changes 'correct' column to use 0/1 and replaces NAs in the score
            difference columns with small values. """

        def _fill_nas(_df, nm):
            with warnings.catch_warnings(record=True) as w:
                _df[nm] = _df[nm].fillna(np.nanmax(_df[nm])+1).fillna(0)
                if len(w) > 0:
                    assert len(w) == 1
                    assert issubclass(w[0].category, RuntimeWarning)
                    _df[nm] = _df[nm].fillna(0)
            assert not any([math.isnan(x) for x in _df[nm]])

        if df.shape[0] == 0:
            return

        for col in df:
            if df[col].dtype != 'object':
                _fill_nas(df, col)
                assert not math.isnan(df[col].sum())

        return df

    def dataset_iter(self, sn):
        """ Return an iterator over chunks of rows from the data frame. """
        assert sn in self.readers
        return imap(lambda x: self._postprocess_data_frame(x), self.readers[sn]())

    def __contains__(self, o):
        return o in self.readers
