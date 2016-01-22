import pandas
import numpy
from mapq import mapq_to_pcor_np, pcor_to_mapq_np


class Roc(object):
    """
    Class encapsulating # correct/incorrect information stratified by MAPQ.
    Easy to derive various accuracy measures and diagnostic plots from this.
    """

    def __init__(self, tally, mapq_strata=True):
        mapqs, tups = zip(*sorted(tally.items(), reverse=True))
        mapqs = numpy.array(mapqs)
        if mapq_strata:
            pcor = mapq_to_pcor_np(mapqs)
        else:
            pcor = mapqs
            mapqs = pcor_to_mapq_np(mapqs)
        cors, incors = zip(*tups)
        self.tab = pandas.DataFrame.from_dict({'mapq': mapqs,
                                               'pcor': pcor,
                                               'cor': cors,
                                               'incor': incors,
                                               'cum_cor': numpy.cumsum(cors),
                                               'cum_incor': numpy.cumsum(incors)})
        self.tab['n'] = self.tab['cor'] + self.tab['incor']
        self.tab['cum'] = self.tab['cum_cor'] + self.tab['cum_incor']
        self.tab['se'] = self.tab['incor'] * self.tab['pcor'] * self.tab['pcor'] + \
                         self.tab['cor'] * (1.0 - self.tab['pcor']) * (1.0 - self.tab['pcor'])
        self.tab['cum_se'] = self.tab['se'].cumsum()
        self.tot = self.tab['n'].sum()

    def cum_incorrect_and_error(self):
        """
        Take two ROC tables.  Return vectors corresponding to the CID and CSED
        curves.
        """
        ci, ce = [0], [0]
        for idx, row in self.tab.iterrows():
            incor_incr = row['incor'] / float(row['n'])
            se_incr = row['se'] / float(row['n'])
            for i in range(int(row['n'])):
                ci.append(ci[-1] + incor_incr)
                ce.append(ce[-1] + se_incr)
        return ci, ce

    @staticmethod
    def write_cum_incorrect_diff(roc1, roc2, fn):
        with open(fn, 'wb') as fh:
            ci1l, _ = roc1.cum_incorrect_and_error()
            ci2l, _ = roc2.cum_incorrect_and_error()
            for ci1, ci2 in zip(ci1l, ci2l):
                fh.write(str(ci1 - ci2) + '\n')

    @staticmethod
    def write_cum_squared_error(roc1, roc2, fn):
        with open(fn, 'wb') as fh:
            _, csel1 = roc1.cum_incorrect_and_error()
            _, csel2 = roc2.cum_incorrect_and_error()
            for cse1, cse2 in zip(csel1, csel2):
                fh.write(str(cse1 - cse2) + '\n')

    def area_under_cumulative_incorrect(self):
        """
        Return area under the cumulative incorrect curve, accumulated from
        high to low mapping quality.
        """
        cum, auc = 0, 0.0
        for idx, row in self.tab.iterrows():
            ic, n = int(row['incor']), int(row['n'])
            height = cum + ic/2.0
            width = n
            auc += height * width
            cum += ic
        return auc

    def sum_of_squared_error(self):
        """
        Return sum of squared errors.
        """
        return self.tab['cum_se'].iloc[-1]


if __name__ == "__main__":

    import sys
    import unittest

    class TestCases(unittest.TestCase):

        def test_roc_1(self):
            roc = Roc({2: [1, 1],
                       1: [1, 2],
                       0: [1, 1]})
            self.assertEqual(list(roc.tab['cum']), [2, 5, 7])
            self.assertEqual(list(roc.tab['cum_incor']), [1, 3, 4])
            self.assertEqual(list(roc.tab['cum_cor']), [1, 2, 3])
            self.assertEqual(list(roc.tab['n']), [2, 3, 2])

        def test_roc_2(self):
            roc = Roc({0.0: [1, 1],
                       1.0: [1, 2]}, mapq_strata=False)
            self.assertEqual(list(roc.tab['cum']), [3, 5])
            self.assertEqual(list(roc.tab['se']), [2.0, 1.0])
            self.assertEqual(list(roc.tab['cum_se']), [2.0, 3.0])
            self.assertEqual(list(roc.tab['n']), [3, 2])

        def test_cum_inc_and_err_1(self):
            roc = Roc({0.0: [1, 1],
                       1.0: [1, 2]}, mapq_strata=False)
            ci, ce = roc.cum_incorrect_and_error()
            self.assertEqual(ci, [0, 2.0/3, 4.0/3, 6.0/3, 2.5, 3.0])
            self.assertEqual(ce, [0, 2.0/3, 4.0/3, 6.0/3, 2.5, 3.0])

        def test_cum_inc_and_err_2(self):
            roc = Roc({0.0: [1, 1],
                       0.1: [0, 1],
                       0.9: [1, 0],
                       1.0: [1, 2]}, mapq_strata=False)
            ci, ce = roc.cum_incorrect_and_error()
            self.assertEqual(ci, [0, 2.0/3, 4.0/3, 6.0/3, 2.0, 3.0, 3.5, 4.0])
            ex = [0, 2.0/3, 4.0/3, 6.0/3, 2.01, 2.02, 2.52, 3.02]
            for x, y in zip(ce, ex):
                self.assertAlmostEqual(x, y, places=5)

        def test_auc_1(self):
            roc = Roc({2: [1, 1],
                       1: [1, 2],
                       0: [1, 1]})
            self.assertEqual(0.5 * 2.0 + 2.0 * 3.0 + 3.5 * 2.0, roc.area_under_cumulative_incorrect())

        def test_sse_1(self):
            roc = Roc({0.0: [1, 1],
                       0.1: [0, 1],
                       0.9: [1, 0],
                       1.0: [1, 2]}, mapq_strata=False)
            self.assertAlmostEqual(3.02, roc.sum_of_squared_error())

    unittest.main(argv=[sys.argv[0]])
    sys.exit()
