__author__ = 'langmead'

import random
from read import Alignment
from collections import defaultdict
from randutil import ReservoirSampler


class ReadTemplate(object):
    """ Info about an aligned read we can use to generate others like it """
    def __init__(self, sc, fw, qual, rd_aln, rf_aln, rf_len, mate1, olen):
        self.sc = sc
        self.fw = fw
        self.qual = qual
        self.rd_aln = rd_aln
        self.rf_aln = rf_aln
        self.rf_len = rf_len
        self.mate1 = mate1
        self.ordlen = olen

    @property
    def fraglen(self):
        return self.rf_len


class PairTemplate(object):

    """ Info about an aligned pair we can use to generate others like it """
    def __init__(self, rt1, rt2, fraglen, upstream1):
        self.rt1 = rt1
        self.rt2 = rt2
        self._fraglen = fraglen
        self.upstream1 = upstream1

    @property
    def fraglen(self):
        return self._fraglen


class ScoreDist(object):

    """ Like CollapsedScoreDist but with fraction_even=1.0 """

    def __init__(self, small_k=100, big_k=10000):
        self.sample = ReservoirSampler(big_k)
        self.max_fraglen = 0
        self.avg_fraglen = None
        self.finalized = False
        self.num_added = 0
        self.tot_len = 0
        self.num_drawn = 0
        self.has_correctness_info = False

    def draw(self):
        assert self.finalized
        assert not self.empty()
        self.num_drawn += 1
        score, fw, qual, rd_aln, rf_aln, rf_len, mate1, olen = self.sample.draw()
        return ReadTemplate(score, fw, qual, rd_aln, rf_aln, rf_len, mate1, olen)

    def add(self, al, correct, ref, ordlen=0, use_ref_for_edit_distance=False):
        pos = self.sample.add_step_1()
        if pos is not None:
            sc = al.bestScore
            rd_aln, rf_aln, rd_len, rf_len =\
                al.stacked_alignment(use_ref_for_edit_distance=use_ref_for_edit_distance, ref=ref)
            self.sample.add_step_2(pos, (sc, al.fw, al.qual, rd_aln, rf_aln, rf_len, al.mate1, ordlen))
            self.max_fraglen = max(self.max_fraglen, rf_len)
            self.tot_len += rf_len
        else:
            self.max_fraglen = max(self.max_fraglen, len(al.seq))
            self.tot_len += len(al.seq)
        self.num_added += 1

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_added == 0

    def finalize(self):
        """ Sort the samples in preparation for draws that might be biased
            toward high or (more likely) low scores. """
        self.finalized = True
        if not self.empty():
            self.avg_fraglen = float(self.tot_len) / self.num_added

class CollapsedScoreDist(object):

    def __init__(self, small_k=100, big_k=10000, fraction_even=0.5, bias=1.0):
        self.score_to_sample = {}
        self.score_to_fraction_correct = defaultdict(lambda: [0, 0])
        self.scores = None
        self.sample = ReservoirSampler(big_k)
        self.max_fraglen = 0
        self.avg_fraglen = None
        self.finalized = False
        self.num_added = 0
        self.tot_len = 0
        self.num_drawn = 0
        self.k = small_k
        self.correct_mass = 0
        self.has_correctness_info = False
        self.fraction_even = fraction_even
        self.bias = bias

    def draw(self):
        assert self.finalized
        assert not self.empty()
        if self.scores is None:
            self.scores = sorted(self.score_to_sample.iterkeys())
        self.num_drawn += 1
        if random.random() > self.fraction_even:
            if self.bias > 1.0:
                rand_i = random.uniform(0, 1) / random.uniform(1.0, self.bias)
                score = self.scores[int(rand_i * len(self.scores))]
            else:
                score = random.choice(self.score_to_sample.keys())
            tup = random.choice(self.score_to_sample[score].r)
            fw, qual, rd_aln, rf_aln, rf_len, mate1, olen = tup
        else:
            score, fw, qual, rd_aln, rf_aln, rf_len, mate1, olen = self.sample.draw()
        if self.has_correctness_info:
            p_correct = float(self.score_to_fraction_correct[score][0]) / self.score_to_fraction_correct[score][1]
            self.correct_mass += p_correct
        return ReadTemplate(score, fw, qual, rd_aln, rf_aln, rf_len, mate1, olen)

    def add(self, al, correct, ref, ordlen=0, use_ref_for_edit_distance=False):
        sc = al.bestScore
        # TODO: don't call stacked_alignment unless we have to -- some calls
        # to sample.add will not add to any reservoirs
        rd_aln, rf_aln, rd_len, rf_len =\
            al.stacked_alignment(use_ref_for_edit_distance=use_ref_for_edit_distance, ref=ref)
        self.max_fraglen = max(self.max_fraglen, rf_len)
        self.tot_len += rf_len
        if sc not in self.score_to_sample:
            self.score_to_sample[sc] = ReservoirSampler(self.k)
        self.score_to_sample[sc].add((al.fw, al.qual, rd_aln, rf_aln, rf_len, al.mate1, ordlen))
        self.sample.add((sc, al.fw, al.qual, rd_aln, rf_aln, rf_len, al.mate1, ordlen))
        if correct is not None:
            self.score_to_fraction_correct[sc][0] += 1 if correct else 0
            self.score_to_fraction_correct[sc][1] += 1
            self.has_correctness_info = True
        self.num_added += 1

    def frac_correct(self):
        return float(self.correct_mass) / self.num_drawn

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_added == 0

    def finalize(self):
        """ Sort the samples in preparation for draws that might be biased
            toward high or (more likely) low scores. """
        self.finalized = True
        if not self.empty():
            self.avg_fraglen = float(self.tot_len) / self.num_added

class ScorePairDist(object):

    def __init__(self, small_k=30, big_k=10000, max_allowed_fraglen=100000):
        """ Make a reservoir sampler for holding the tuples. """
        self.sample = ReservoirSampler(big_k)
        self.max_allowed_fraglen = max_allowed_fraglen
        self.max_fraglen = 0  # maximum observed fragment length
        self.avg_fraglen = None
        self.finalized = False
        self.tot_len = 0
        self.k = small_k
        self.num_drawn = 0
        self.num_added = 0
        self.has_correctness_info = False

    def draw(self):
        """ Draw from the reservoir """
        assert self.finalized
        assert not self.empty()
        self.num_drawn += 1
        score, tup1, tup2, fraglen, upstream1 = self.sample.draw()
        fw_1, qual_1, rd_aln_1, rf_aln_1, sc_1, rf_len_1, mate1_1, _ = tup1
        fw_2, qual_2, rd_aln_2, rf_aln_2, sc_2, rf_len_2, mate1_2, _ = tup2
        return PairTemplate(ReadTemplate(sc_1, fw_1, qual_1, rd_aln_1, rf_aln_1, rf_len_1, mate1_1, rf_len_2),
                            ReadTemplate(sc_2, fw_2, qual_2, rd_aln_2, rf_aln_2, rf_len_2, mate1_2, rf_len_1),
                            fraglen, upstream1)

    def add(self, al1, al2, correct1, correct2, ref, use_ref_for_edit_distance=False):
        """ Convert given alignment pair to a tuple and add it to the
            reservoir sampler. """
        fraglen = Alignment.fragment_length(al1, al2)
        fraglen = min(fraglen, self.max_allowed_fraglen)
        self.max_fraglen = max(self.max_fraglen, fraglen)
        pos = self.sample.add_step_1()
        if pos is not None:
            sc1, sc2 = al1.bestScore, al2.bestScore
            # Make note of which end is upstream
            upstream1 = al1.pos < al2.pos
            # Get stacked alignment
            rd_aln_1, rf_aln_1, rd_len_1, rf_len_1 =\
                al1.stacked_alignment(use_ref_for_edit_distance=use_ref_for_edit_distance, ref=ref)
            rd_aln_2, rf_aln_2, rd_len_2, rf_len_2 =\
                al2.stacked_alignment(use_ref_for_edit_distance=use_ref_for_edit_distance, ref=ref)
            assert fraglen == 0 or max(rf_len_1, rf_len_2) <= fraglen
            self.sample.add_step_2(pos, (sc1 + sc2,
                                         (al1.fw, al1.qual, rd_aln_1, rf_aln_1, sc1, rf_len_1, True, rf_len_2),
                                         (al2.fw, al2.qual, rd_aln_2, rf_aln_2, sc2, rf_len_2, False, rf_len_1),
                                         fraglen, upstream1))
        self.num_added += 1
        self.tot_len += fraglen

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_added == 0

    def finalize(self):
        """ Sort the samples in preparation for draws that might be biased
            toward high or (more likely) low scores. """
        self.finalized = True
        if not self.empty():
            self.avg_fraglen = float(self.tot_len) / self.num_added


class CollapsedScorePairDist(object):

    def __init__(self, small_k=30, big_k=10000, max_allowed_fraglen=100000, fraction_even=0.5, bias=1.0):
        """ Make a reservoir sampler for holding the tuples. """
        self.score_to_sample = {}
        self.sample = ReservoirSampler(big_k)
        self.score_to_fraction_correct1 = defaultdict(lambda: [0, 0])
        self.score_to_fraction_correct2 = defaultdict(lambda: [0, 0])
        self.scores = None
        self.max_allowed_fraglen = max_allowed_fraglen
        self.max_fraglen = 0  # maximum observed fragment length
        self.avg_fraglen = None
        self.finalized = False
        self.tot_len = 0
        self.k = small_k
        self.num_drawn = 0
        self.num_added = 0
        self.correct_mass1 = 0
        self.correct_mass2 = 0
        self.has_correctness_info = False
        self.fraction_even = fraction_even
        self.bias = bias

    def draw(self):
        """ Draw from the reservoir """
        assert self.finalized
        assert not self.empty()
        if self.scores is None:
            self.scores = sorted(self.score_to_sample.keys())
        self.num_drawn += 1
        if random.random() > self.fraction_even:
            if self.bias > 1.0:
                rand_i = random.uniform(0, 1) / random.uniform(1, self.bias)
                score = self.scores[int(rand_i * len(self.scores))]
            else:
                score = random.choice(self.score_to_sample.keys())
            tup1, tup2, fraglen, upstream1 = random.choice(self.score_to_sample[score].r)
        else:
            score, tup1, tup2, fraglen, upstream1 = self.sample.draw()
        fw_1, qual_1, rd_aln_1, rf_aln_1, sc_1, rf_len_1, mate1_1, _ = tup1
        fw_2, qual_2, rd_aln_2, rf_aln_2, sc_2, rf_len_2, mate1_2, _ = tup2
        if self.has_correctness_info:
            p_correct1 = float(self.score_to_fraction_correct1[score][0]) / self.score_to_fraction_correct1[score][1]
            self.correct_mass1 += p_correct1
            p_correct2 = float(self.score_to_fraction_correct2[score][0]) / self.score_to_fraction_correct2[score][1]
            self.correct_mass2 += p_correct2
        return PairTemplate(ReadTemplate(sc_1, fw_1, qual_1, rd_aln_1, rf_aln_1, rf_len_1, mate1_1, rf_len_2),
                            ReadTemplate(sc_2, fw_2, qual_2, rd_aln_2, rf_aln_2, rf_len_2, mate1_2, rf_len_1),
                            fraglen, upstream1)

    def add(self, al1, al2, correct1, correct2, ref, use_ref_for_edit_distance=False):
        """ Convert given alignment pair to a tuple and add it to the
            reservoir sampler. """
        # TODO: don't call stacked_alignment unless we have to -- some calls
        # to sample.add will not add to any reservoirs
        sc1, sc2 = al1.bestScore, al2.bestScore
        # Make note of fragment length
        fraglen = Alignment.fragment_length(al1, al2)
        fraglen = min(fraglen, self.max_allowed_fraglen)
        self.max_fraglen = max(self.max_fraglen, fraglen)
        # Make note of which end is upstream
        upstream1 = al1.pos < al2.pos
        # Get stacked alignment
        rd_aln_1, rf_aln_1, rd_len_1, rf_len_1 =\
            al1.stacked_alignment(use_ref_for_edit_distance=use_ref_for_edit_distance, ref=ref)
        rd_aln_2, rf_aln_2, rd_len_2, rf_len_2 =\
            al2.stacked_alignment(use_ref_for_edit_distance=use_ref_for_edit_distance, ref=ref)
        assert fraglen == 0 or max(rf_len_1, rf_len_2) <= fraglen
        score = sc1 + sc2
        if score not in self.score_to_sample:
            self.score_to_sample[score] = ReservoirSampler(self.k)
        if correct1 is not None:
            self.score_to_fraction_correct1[score][0] += 1 if correct1 else 0
            self.score_to_fraction_correct1[score][1] += 1
            self.score_to_fraction_correct2[score][0] += 1 if correct2 else 0
            self.score_to_fraction_correct2[score][1] += 1
            self.has_correctness_info = True
        self.score_to_sample[score].add(((al1.fw, al1.qual, rd_aln_1, rf_aln_1, sc1, rf_len_1, True, rf_len_2),
                                         (al2.fw, al2.qual, rd_aln_2, rf_aln_2, sc2, rf_len_2, False, rf_len_1),
                                         fraglen, upstream1))
        self.sample.add((score,
                         (al1.fw, al1.qual, rd_aln_1, rf_aln_1, sc1, rf_len_1, True, rf_len_2),
                         (al2.fw, al2.qual, rd_aln_2, rf_aln_2, sc2, rf_len_2, False, rf_len_1),
                         fraglen, upstream1))
        self.num_added += 1
        self.tot_len += fraglen

    def frac_correct1(self):
        return float(self.correct_mass1) / self.num_drawn

    def frac_correct2(self):
        return float(self.correct_mass2) / self.num_drawn

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_added == 0

    def finalize(self):
        """ Sort the samples in preparation for draws that might be biased
            toward high or (more likely) low scores. """
        self.finalized = True
        if not self.empty():
            self.avg_fraglen = float(self.tot_len) / self.num_added
