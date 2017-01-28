"""
Copyright 2016, Ben Langmead <langmea@cs.jhu.edu>

FeatureTableReader
"""

import os
import warnings
import math
import numpy as np
try:
    from itertools import imap
except ImportError:
    imap = map

# qtip imports
from metamat import MetaMat


__author__ = 'langmead'


class FeatureTableReader(object):
    """ Reads a table of information describing alignments.  These are tables
        output by qtip.  Tables might describe training or test alignments. """

    #            short
    #            name   suffix
    datasets = [('u',   '_rec_u'),
                ('d',   '_rec_d'),
                ('c',   '_rec_c'),
                ('b',   '_rec_b')]

    def __init__(self, prefix, chunksize=100000):
        self.prefix = prefix
        self.dfs = {}
        self.readers = {}
        nonempty = False
        fns = []
        for sn, suf in self.datasets:
            fn = self.prefix + suf
            fns.append(fn)
            if os.path.exists(fn + '.npy') and os.stat(fn + '.npy').st_size > 0:
                nonempty = True
                self.readers[sn] = MetaMat(fn, chunksize)

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
        self.readers[sn].reset()
        return imap(lambda x: self._postprocess_data_frame(x), self.readers[sn])

    def __contains__(self, o):
        return o in self.readers
