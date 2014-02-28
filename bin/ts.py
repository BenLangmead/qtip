#!/usr/bin/env python

"""
ts.py

A "tandem simulator," which wraps an alignment tool as it runs, eavesdrops on
the input and output, and builds a model that can be used to improve the
quality values calculated for aligned reads.

In a typical run, the user supplies an index and some read files to ts.py,
ts.py runs Bowtie 2 and produces alignments.  The reads and alignments are
used to create a model of the input data.  The model is then used to simulate

Output files encapsulate:
1. Input data model
2. Simulated reads
3. Alignments for simulated reads
4. (3) converted into training-data records
5. Trained models
6. Results of running the trained models on the training data

Paired-end reads come with several additional issues
- Input has mix of alignment types: UP, DP, and CP, where probably only CP is
  well represented.
  + Do we actively try to simulate UP and DP?  How do you simulate a DP?  It
    might be discordant because the fragment spans a breakpoint, but we don't
    really know why.
- We're eventually going to build models that have to be good enough to
  estimate MAPQs for all three kinds

Things we learn from reads
==========================

- Read length distribution
- Quality values

Things we learn from alignments
===============================

- Alignment type (aligned, unaligned, concordant, discordant)
- Fragment length distribution
- Number and placement of mismatches and gaps

SAM extra fields used
=====================

 Normal: AS:i, XS:i, MD:Z
 
 Bowtie-specific: Xs:i, YT:Z, YS:i, Zp:i, YN:i, Yn:i
 (and we potentially don't need YS:i or YT:Z?)

TODO:
- Handle large fragment lengths gracefully.  We have a 10K nt cutoff now.

"""

__author__ = "Ben Langmead"
__email__ = "langmea@cs.jhu.edu"

# depends on bowtie2.py, bwamem.py, randutil.py, sam.py, samples.py, simualtors.py

import os
import sys
from string import maketrans
import random
import time
import logging
import gzip
import errno
import numpy as np
from collections import defaultdict
from traceback import print_tb
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full  # python 3.x
from threading import Thread

# Modules that are part of the tandem simulator
import simulators
from samples import Dataset
from randutil import ReservoirSampler, WeightedRandomGenerator
from simplesim import FragmentSimSerial2, mutate
from read import Read, Alignment
from bowtie2 import AlignmentBowtie2, Bowtie2
from bwamem import AlignmentBwaMem, BwaMem
from reference import ReferenceIndexed, ReferenceOOB

_revcomp_trans = maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s[::-1].translate(_revcomp_trans)


class ScoreDist(object):
    """ Capture a list of tuples, where each tuples represents the following
        traits observed in an alignment: (a) orientation, (b) quality string,
        (c) read side of stacked alignment, (d) ref side of stacked alignment,
        (e) score.  Keeping these tuples allows us to generate new reads that
        mimic observed reads in these key ways. """
    
    def __init__(self, k=10000):
        """ Make a reservoir sampler for holding the tuples. """
        self.res = ReservoirSampler(k)
        self.max_fraglen = 0
        self.avg_fraglen = None

    def finalize(self):
        self.avg_fraglen = np.mean([x[5] for x in self.res])

    def draw(self):
        """ Draw from the reservoir """
        assert not self.empty()
        return self.res.draw()
    
    def add(self, al, ref, ordlen=0):
        """ Convert given alignment to a tuple and add it to the reservoir
            sampler. """
        sc = al.bestScore
        # Get stacked alignment
        rd_aln, rf_aln = al.stacked_alignment(alignSoftClipped=True, ref=ref)
        rl = len(rf_aln) - rf_aln.count('-')
        self.max_fraglen = max(self.max_fraglen, rl)
        self.res.add((al.fw, al.qual, rd_aln, rf_aln, sc, rl, al.mate1, ordlen))
    
    def num_added(self):
        """ Return the number of tuples that have been added, which could be
            much greater than the number sampled """
        return self.res.num_added()

    def num_sampled(self):
        """ Return number of tuples that have been sampled """
        return self.res.num_sampled()
    
    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_sampled() == 0


class ScoreDistEven(object):
    """ Capture a list of tuples, where each tuples represents the following
        traits observed in an alignment: (a) orientation, (b) quality string,
        (c) read side of stacked alignment, (d) ref side of stacked alignment,
        (e) score.  Keeping these tuples allows us to generate new reads that
        mimic observed reads in these key ways.

        The "even" version of this class attempts to make the sample
        relatively uniform across scores.  Without this, our sample will be
        hugely biased toward very high-scoring alignments, which are typical
        for high-quality datasets. """

    def __init__(self, k=10000, num_strata=10, train_from=1000):
        """ Make a reservoir sampler for holding the tuples. """
        self.ress = [ReservoirSampler(k // num_strata) for _ in xrange(num_strata)]
        self.staged_for_training = []
        self.num_strata = num_strata
        self.num_added = 0
        self.train_from = train_from
        self.trained = False
        self.max_fraglen = 0
        self.min_sc, self.max_sc = float('inf'), float('-inf')
        self.avg_fraglen = None

    def _bucketize_score(self, sc):
        """ Convert an alignment score to a bucket index """
        standardized = (float(sc) - self.min_sc) / (self.max_sc - self.min_sc)
        standardized = max(min(standardized, 1.0), 0.0)
        return int(round(standardized * (self.num_strata-1)))

    def _train(self):
        """ Train using data staged so far. """
        for al in self.staged_for_training:
            self.min_sc = min(self.min_sc, al[4])
            self.max_sc = max(self.max_sc, al[4])
        for al in self.staged_for_training:
            self._add_post_training(al)
        self.staged_for_training = []
        self.trained = True

    def finalize(self):
        """ Train using data staged so far. """
        self.avg_fraglen = np.mean([x[5] for x in self.res])
        if not self.trained:
            self._train()

    def draw(self):
        """ Draw from the reservoir """
        if self.empty():
            raise RuntimeError('Call add before draw')
        if not self.trained and len(self) > 0:
            raise RuntimeError('Call finalize before draw')
        resi = random.randint(0, len(self.ress)-1)
        while self.ress[resi].empty():
            resi = random.randint(0, len(self.ress)-1)
        return self.ress[resi].draw()

    def _tupleize_alignment(self, al, ordlen, ref):
        """ Turn alignment object into relevant tuple: (a) orientation,
            (b) quality string, (c) read half of stacked alignment, (d)
            reference half of stacked alignment, (e) score. """
        sc = al.bestScore
        rd_aln, rf_aln = al.stacked_alignment(alignSoftClipped=True, ref=ref)
        rl = len(rf_aln) - rf_aln.count('-')
        self.max_fraglen = max(self.max_fraglen, rl)
        return al.fw, al.qual, rd_aln, rf_aln, sc, rl, al.mate1, ordlen

    def add(self, al, ref, ordlen=0):
        """ Convert given alignment to a tuple and add it to the reservoir
            sampler. """
        tup = self._tupleize_alignment(al, ordlen, ref)
        if self.trained:
            self._add_post_training(tup)
        else:
            self.staged_for_training.append(tup)
            assert len(self.staged_for_training) <= self.train_from
            if len(self.staged_for_training) == self.train_from:
                self.finalize()
        self.num_added += 1

    def _add_post_training(self, tup):
        fw, qual, rd_aln, rf_aln, sc = tup
        bucket = self._bucketize_score(sc)
        assert 0 <= bucket < len(self.ress)
        self.ress[bucket].add(tup)

    def num_added(self):
        """ Return the number of tuples that have been added, which could be
            much greater than the number sampled """
        return self.res.num_added()

    def num_sampled(self):
        """ Return number of tuples that have been sampled """
        return self.res.num_sampled()

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_sampled() == 0


class ScorePairDist(object):
    """ Capture a list of tuple pairs, where each tuple pair represents the
        following in an alignment of a paired-end reads where both ends
        aligned: (a) orientation, (b) quality string, (c) read side of stacked
        alignment, (d) ref side of stacked alignment, (e) score, (f) fragment
        length.  We record pairs of tuples where the first element of the pair
        corresponds to mate 1 and the second to mate 2. """
    
    def __init__(self, k=10000):
        """ Make a reservoir sampler for holding the tuples. """
        self.res = ReservoirSampler(k)
        self.max_allowed_fraglen = 50000  # maximum allowed fragment length
        self.max_fraglen = 0  # maximum observed fragment length
        self.avg_fraglen = None

    def finalize(self):
        self.avg_fraglen = np.mean([x[2] for x in self.res])

    def draw(self):
        """ Draw from the reservoir """
        assert not self.empty()
        return self.res.draw()
    
    def add(self, al1, al2, ref):
        """ Convert given alignment pair to a tuple and add it to the
            reservoir sampler. """
        sc1, sc2 = al1.bestScore, al2.bestScore
        # Make note of fragment length
        fraglen = Alignment.fragment_length(al1, al2)
        self.max_fraglen = max(self.max_fraglen, fraglen)
        # Make note of which end is upstream
        upstream1 = al1.pos < al2.pos
        # Get stacked alignment
        rd_aln_1, rf_aln_1 = al1.stacked_alignment(alignSoftClipped=True, ref=ref)
        rd_aln_2, rf_aln_2 = al2.stacked_alignment(alignSoftClipped=True, ref=ref)
        rlen1 = len(rf_aln_1) - rf_aln_1.count('-')
        rlen2 = len(rf_aln_2) - rf_aln_2.count('-')
        self.res.add(((al1.fw, al1.qual, rd_aln_1, rf_aln_1, sc1, rlen1, True, rlen2),
                      (al2.fw, al2.qual, rd_aln_2, rf_aln_2, sc2, rlen2, False, rlen1), fraglen, upstream1))
    
    def num_added(self):
        """ Return the number of tuples that have been added, which could be
            much greater than the number sampled """
        return self.res.num_added()

    def num_sampled(self):
        """ Return number of tuples that have been sampled """
        return self.res.num_sampled()

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_sampled() == 0


class ScorePairDistEven(object):
    """ Capture a list of tuple pairs, where each tuple pair represents the
        following in an alignment of a paired-end reads where both ends
        aligned: (a) orientation, (b) quality string, (c) read side of stacked
        alignment, (d) ref side of stacked alignment, (e) score, (f) fragment
        length.  We record pairs of tuples where the first element of the pair
        corresponds to mate 1 and the second to mate 2.

        The "even" version of this class attempts to make the sample
        relatively uniform across scores.  Without this, our sample will be
        hugely biased toward very high-scoring alignments, which are typical
        for high-quality datasets. """

    def __init__(self, k=10000, num_strata=10, train_from=1000):
        """ Make a reservoir sampler for holding the tuples. """
        self.ress = [ReservoirSampler(k // num_strata) for _ in xrange(num_strata)]
        self.staged_for_training = []
        self.num_strata = num_strata
        self.num_added = 0
        self.train_from = train_from
        self.trained = False
        self.max_allowed_fraglen = 50000  # maximum allowed fragment length
        self.max_fraglen = 0  # maximum observed fragment length
        self.min_sc, self.max_sc = float('inf'), float('-inf')
        self.avg_fraglen = None

    def _bucketize_score(self, sc):
        """ Convert an alignment score to a bucket index """
        standardized = (float(sc) - self.min_sc) / (self.max_sc - self.min_sc)
        standardized = max(min(standardized, 1.0), 0.0)
        return int(round(standardized * (self.num_strata-1)))

    def _train(self):
        """ Train using data staged so far. """
        for al1, al2, _ in self.staged_for_training:
            sc = al1[4] + al2[4]
            self.min_sc = min(self.min_sc, sc)
            self.max_sc = max(self.max_sc, sc)
        for al1, al2, fraglen, upstream1 in self.staged_for_training:
            self._add_post_training(al1, al2, fraglen, upstream1)
        self.staged_for_training = []
        self.trained = True

    def finalize(self):
        """ Train using data staged so far. """
        self.avg_fraglen = np.mean([x[2] for x in self.res])
        if not self.trained:
            self._train()

    def draw(self):
        """ Draw from the reservoir """
        if self.empty():
            raise RuntimeError('Call add before draw')
        if not self.trained and len(self) > 0:
            raise RuntimeError('Call finalize before draw')
        resi = random.randint(0, len(self.ress)-1)
        while self.ress[resi].empty():
            resi = random.randint(0, len(self.ress)-1)
        return self.ress[resi].draw()

    def _tupleize_alignments(self, al1, al2, ref):
        """ Turn alignment object into relevant tuple: (a) orientation,
            (b) quality string, (c) read half of stacked alignment, (d)
            reference half of stacked alignment, (e) score. """
        sc1, sc2 = al1.bestScore, al2.bestScore
        rd_aln_1, rf_aln_1 = al1.stacked_alignment(alignSoftClipped=True, ref=ref)
        rd_aln_2, rf_aln_2 = al2.stacked_alignment(alignSoftClipped=True, ref=ref)
        # Record fragment length
        fraglen = Alignment.fragment_length(al1, al2)
        # Record whether mate 1 is upstream
        upstream1 = al1.pos < al2.pos
        rlen1 = len(rf_aln_1) - rf_aln_1.count('-')
        rlen2 = len(rf_aln_2) - rf_aln_2.count('-')
        return (al1.fw, al1.qual, rd_aln_1, rf_aln_1, sc1, rlen1, True, rlen2), \
               (al2.fw, al2.qual, rd_aln_2, rf_aln_2, sc2, rlen2, False, rlen1), fraglen, upstream1

    def add(self, al1, al2, ref):
        """ Convert given alignment pair to a tuple and add it to the
            reservoir sampler. """
        al1tup, al2tup, fraglen, upstream1 = self._tupleize_alignments(al1, al2, ref)
        self.max_fraglen = max(self.max_fraglen, fraglen)
        if self.trained:
            self._add_post_training(al1tup, al2tup, fraglen, upstream1)
        else:
            self.staged_for_training.append((al1tup, al2tup, fraglen, upstream1))
            assert len(self.staged_for_training) <= self.train_from
            if len(self.staged_for_training) == self.train_from:
                self.finalize()
        self.num_added += 1

    def _add_post_training(self, al1tup, al2tup, fraglen, upstream1):
        _, _, _, _, sc1 = al1tup
        _, _, _, _, sc2 = al2tup
        bucket = self._bucketize_score(sc1 + sc2)
        assert 0 <= bucket < len(self.ress)
        self.ress[bucket].add((al1tup, al2tup, fraglen, upstream1))

    def num_added(self):
        """ Return the number of tuples that have been added, which could be
            much greater than the number sampled """
        return self.res.num_added()

    def num_sampled(self):
        """ Return number of tuples that have been sampled """
        return self.res.num_sampled()

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_sampled() == 0


class Dist(object):
    """ Basically a histogram. """
    
    def __init__(self):
        self.hist = defaultdict(int)
        self.tot = 0
        self.changed = True
        self.gen = None
        self.max, self.min = None, None
    
    def __str__(self):
        """ Return string representation of histogram """
        return str(self.hist)
    
    def __len__(self):
        """ Return the number of items added to the histogram """
        return self.tot
    
    def empty(self):
        """ Return true iff there are no elements in the histogram """
        return len(self) == 0
    
    def draw(self):
        """ Draw an element from the histogram """
        if len(self) == 0:
            raise RuntimeError("Attempt to draw from empty empirical distribution")
        if self.changed:
            self.gen = WeightedRandomGenerator(self.hist)
            self.changed = False
        return self.gen.next()
    
    def add(self, key):
        """ Add new element to the histogram """
        self.hist[key] += 1
        self.tot += 1
        self.changed = True
        # keep track of min and max
        if self.max is None or key > self.max:
            self.max = key
        if self.min is None or key < self.min:
            self.min = key


class SimulatorWrapper(object):

    # This doesn't quite work because we're returning iterators.  The user
    # can't mix and match calls to these various methods.  Have to think
    # more about the interface between this and its caller.  Another reason
    # we might want this class to pick which kind of read gets generated for
    # each draw: otherwise we might be biased against the last one we pick
    # if we run out of reference.

    def __init__(self, frag_iter, dists):
        self.frag_iter = frag_iter
        self.dists = dists

    def _sim_fragment(self):
        try:
            # Try to simulate a fragment
            refid, refoff, seq = self.frag_iter.next()
        except StopIteration:
            return None, None, None  # already went through the whole reference
        return refid, refoff, seq

    def _next_unpaired(self, sc_dist, typ):
        """ Simulate an unpaired read """
        while True:
            sc_draw = sc_dist.draw()
            # sc_draw is tuple: (fw, qual, rd_aln, rf_aln, sc)
            fw, _, _, rf_aln, sc = sc_draw
            rl = len(rf_aln) - rf_aln.count('-')
            refid, refoff, seq = self._sim_fragment()
            if seq is None:
                break  # exhausted reference
            assert len(seq) >= rl
            seq = seq[:rl] if fw else seq[-rl:]
            read = Read.from_simulator(seq, None, refid, refoff, fw, sc, typ)
            mutate(read, fw, sc_draw)  # mutate unpaired read
            assert read.qual is not None
            yield read

    def next_unpaired(self):
        for x in self._next_unpaired(self.dists.sc_dist_unp, 'Unp'):
            yield x

    def next_bad_end(self):
        for x in self._next_unpaired(self.dists.sc_dist_bad_end, 'BadEnd'):
            yield x

    @staticmethod
    def _draw_pair(sc_dist):
        m1fw, m2fw = None, None
        fl, rl1, rl2 = None, None, None
        sc1, sc2, upstream1 = None, None, None
        sc1_draw, sc2_draw = None, None
        # Draw a long enough fragment to accomodate both ends
        while fl < rl1 or fl < rl2:
            sc1_draw, sc2_draw, fl, upstream1 = sc_dist.draw()
            m1fw, _, _, rf_aln_1, sc1 = sc1_draw
            m2fw, _, _, rf_aln_2, sc2 = sc2_draw
            rl1 = len(rf_aln_1) - rf_aln_1.count('-')
            rl2 = len(rf_aln_2) - rf_aln_2.count('-')
        return fl, rl1, rl2, sc1, sc2, m1fw, m2fw, upstream1, sc1_draw, sc2_draw

    def _next_pair(self, sc_dist, typ):
        """ Simulate a paired-end read """
        while True:
            fl, rl1, rl2, sc1, sc2, m1fw, m2fw, upstream1,\
                sc1_draw, sc2_draw = self._draw_pair(sc_dist)
            refid, refoff, seq = self._sim_fragment()
            if seq is None:
                break  # exhausted reference
            assert len(seq) >= fl
            if upstream1:
                seq1, seq2 = seq[:rl1], seq[-rl2:]
            else:
                seq2, seq1 = seq[:rl1], seq[-rl2:]
            if not m1fw:
                seq1 = revcomp(seq1)
            if not m2fw:
                seq2 = revcomp(seq2)
            # Now we have the Watson offset for one mate or the other,
            # depending on which mate is to the left w/r/t Watson.
            if upstream1:
                refoff1 = refoff
                refoff2 = refoff + fl - rl2
            else:
                refoff1 = refoff + fl - rl1
                refoff2 = refoff
            rdp1 = Read.from_simulator(seq1, None, refid, refoff1, m1fw, sc1, typ)
            rdp2 = Read.from_simulator(seq2, None, refid, refoff2, m2fw, sc2, typ)
            mutate(rdp1, m1fw, sc1_draw)
            mutate(rdp2, m2fw, sc2_draw)
            yield rdp1, rdp2

    def next_concordant_pair(self):
        for x in self._next_pair(self.dists.sc_dist_conc, 'Conc'):
            yield x

    def next_discordant_pair(self):
        for x in self._next_pair(self.dists.sc_dist_disc, 'Disc'):
            yield x


class SimulatorWrapperOld(object):
    
    """ Wrapper that sends requests to the Simulator but uses information
        gathered during alignment so far to select such parameters as read
        length, concordant/discordant fragment length, etc.
        
        With each call to next(), we either simulate a paired-end read or an
        unpaired read.
        """
    
    def __init__(self, sim, dists):
        self.sim = sim      # sequence simulator
        self.dists = dists  # empirical distributions
    
    def next_unpaired(self):
        """ Simulate an unpaired read """
        sc_draw = self.dists.sc_dist_unp.draw()
        # sc_draw is tuple: (fw, qual, rd_aln, rf_aln, sc)
        fw, _, _, rf_aln, sc = sc_draw
        rl = len(rf_aln) - rf_aln.count('-')
        refid, refoff, seq = self.sim.sim(rl)  # simulate it
        assert rl == len(seq)
        read = Read.from_simulator(seq, None, refid, refoff, fw, sc, "Unp")
        mutate(read, fw, sc_draw)  # mutate unpaired read
        assert read.qual is not None
        return read
    
    def next_concordant_pair(self):
        """ Simulate a paired-end read """
        fl, rl1, rl2 = None, None, None
        sc1_draw, sc2_draw, upstream1 = None, None, None
        sc1, sc2, m1fw, m2fw = None, None, None, None
        while fl < rl1 or fl < rl2:
            sc1_draw, sc2_draw, fl, upstream1 = self.dists.sc_dist_conc.draw()
            m1fw, _, _, rf_aln_1, sc1 = sc1_draw
            m2fw, _, _, rf_aln_2, sc2 = sc2_draw
            rl1 = len(rf_aln_1) - rf_aln_1.count('-')
            rl2 = len(rf_aln_2) - rf_aln_2.count('-')
        refid, refoff, fw, seq = self.sim.sim(fl)  # simulate fragment
        assert len(seq) == fl
        if upstream1:
            seq1, seq2 = seq[:rl1], seq[-rl2:]
        else:
            seq2, seq1 = seq[:rl1], seq[-rl2:]
        if not m1fw:
            seq1 = revcomp(seq1)
        if not m2fw:
            seq2 = revcomp(seq2)
        # Now we have the Watson offset for one mate or the other,
        # depending on which mate is to the left w/r/t Watson.
        if upstream1:
            refoff1 = refoff
            refoff2 = refoff + fl - rl2
        else:
            refoff1 = refoff + fl - rl1
            refoff2 = refoff
        rdp1 = Read.from_simulator(seq1, None, refid, refoff1, m1fw, sc1, "Conc")
        rdp2 = Read.from_simulator(seq2, None, refid, refoff2, m2fw, sc2, "Conc")
        mutate(rdp1, m1fw, sc1_draw)
        mutate(rdp2, m2fw, sc2_draw)
        return rdp1, rdp2


class Input(object):
    """ Class that parses reads from input files and yields the reads/pairs
        produced using generators """

    @staticmethod
    def fasta_parse(fh):
        """ Parse a single FASTA-format read from given filehandle.  Return
            None if input is exhausted. """
        lns = [fh.readline().rstrip() for _ in xrange(2)]
        orig = '\n'.join(lns) + '\n'
        if len(lns[0]) == 0:
            return None
        return Read(lns[0][1:], lns[1], None, orig)

    @staticmethod
    def fastq_parse(fh):
        """ Parse a single FASTQ-format read from given filehandle.  Return
            None if input is exhausted. """
        lns = [fh.readline().rstrip() for _ in xrange(4)]
        orig = '\n'.join(lns) + '\n'
        if len(lns[0]) == 0:
            return None
        return Read(lns[0][1:], lns[1], lns[3], orig)
    
    def __init__(self, file_format="fastq", unp_fns=None, m1_fns=None, m2_fns=None):
        self.format = file_format
        if file_format == "fastq":
            self.parse = self.fastq_parse
        elif file_format == "fasta":
            self.parse = self.fasta_parse
        else:
            raise RuntimeError("Bad input format: '%s'" % file_format)
        self.unp_fns, self.m1_fns, self.m2_fns = unp_fns, m1_fns, m2_fns
    
    def __iter__(self):
        """ Generator for all the reads.  Yields pairs of Read objects, where
            second element is None for unpaired reads. """
        # Yield all the unpaired reads first
        if self.unp_fns is not None:
            for unpFn in self.unp_fns:
                with gzip.open(unpFn, 'r') if unpFn.endswith('.gz') else open(unpFn, 'r') as unpFh:
                    while True:
                        rd = self.parse(unpFh)
                        if rd is not None:
                            yield (rd, None)
                        else:
                            break  # next file
        # Yield all the paired-end reads
        if self.m1_fns is not None:
            assert self.m2_fns is not None
            for (m1Fn, m2Fn) in zip(self.m1_fns, self.m2_fns):
                with gzip.open(m1Fn, 'r') if m1Fn.endswith('.gz') else open(m1Fn, 'r') as m1Fh:
                    with gzip.open(m2Fn, 'r') if m2Fn.endswith('.gz') else open(m2Fn, 'r') as m2Fh:
                        while True:
                            rd1, rd2 = self.parse(m1Fh), self.parse(m2Fh)
                            if rd1 is not None:
                                yield (rd1, rd2)
                            else:
                                break  # next pair of files


class Dists(object):
    
    """ Encapsulates the distributions that we capture from real data.
        We collect random subsets of qualities/edits or unpaired reads
        and same for paired-end mate 1s and mate 2s.  We also collect
        data on concordant and discordantly aligned pairs, such as
        their fragment length and strands. """

    def __init__(self, even=False):
        self.even = even
        if even:
            self.sc_dist_unp = ScoreDistEven()
            self.sc_dist_bad_end = ScoreDistEven()
            self.sc_dist_conc = ScorePairDistEven()
            self.sc_dist_disc = ScorePairDistEven()
        else:
            self.sc_dist_unp = ScoreDist()
            self.sc_dist_bad_end = ScoreDist()
            self.sc_dist_conc = ScorePairDist()
            self.sc_dist_disc = ScorePairDist()

    def finalize(self):
        self.sc_dist_unp.finalize()
        self.sc_dist_bad_end.finalize()
        self.sc_dist_conc.finalize()
        self.sc_dist_disc.finalize()

    def add_concordant_pair(self, al1, al2, ref):
        """ Add concordant paired-end read alignment to the model """
        self.sc_dist_conc.add(al1, al2, ref)

    def add_discordant_pair(self, al1, al2, ref):
        """ Add discordant paired-end read alignment to the model """
        self.sc_dist_disc.add(al1, al2, ref)

    def add_unpaired_read(self, al, ref):
        """ Add unpaired read alignment to the model """
        self.sc_dist_unp.add(al, ref)

    def add_bad_end_read(self, al, ordlen, ref):
        """ Add bad-end read alignment to the model """
        self.sc_dist_bad_end.add(al, ref, ordlen)

    def has_pairs(self):
        return self.has_concordant_pairs() or self.has_discordant_pairs() or \
            self.has_bad_end_reads()

    def has_concordant_pairs(self):
        """ Return true iff at least one concordant paired-end read was
            added. """
        return not self.sc_dist_conc.empty()

    def has_discordant_pairs(self):
        """ Return true iff at least one concordant paired-end read was
            added. """
        return not self.sc_dist_disc.empty()

    def has_bad_end_reads(self):
        """ Return true iff at least one bad-end was added. """
        return not self.sc_dist_bad_end.empty()

    def has_unpaired_reads(self):
        """ Return true iff at least one unpaired read was added. """
        return not self.sc_dist_unp.empty()

    def avg_concordant_fraglen(self):
        return self.sc_dist_conc.avg_fraglen

    def avg_discordant_fraglen(self):
        return self.sc_dist_disc.avg_fraglen

    def avg_unpaired_readlen(self):
        return self.sc_dist_unp.avg_fraglen

    def longest_fragment(self):
        """ Return length of longest fragment we'll have to sample
            from genome. """
        return max(self.sc_dist_conc.max_fraglen, self.sc_dist_disc.max_fraglen)

    def longest_unpaired(self):
        """ Return length of longest substring we'll have to sample
            from genome in order to simulate an unpaired read using
            this distribution. """
        return self.sc_dist_unp.max_fraglen


def is_correct(al, args):
    """ Return true if the given alignment is both simulated and "correct" --
        i.e. in or very near it's true point of origin """
    if args.correct_chromosomes is not None:
        return (al.refid in args.correct_chromosomes), 0
    elif simulators.isExtendedWgsim(al.name):
        return simulators.correctExtendedWgsim(al, args.wiggle)
    else:
        return None, 0


class AlignmentReader(Thread):
    
    """ Aligner output reader.  Reads SAM output from aligner, updates
        empirical distributions (our "model" for what the input reads look
        like), and looks for records that correspond to simulated reads.  If a
        record corresponds to a simulated read, its correctness is checked and
        an output tuple written """
    
    def __init__(self, args, sam_q, dataset, dists, ref,
                 alignment_class, dist_hist_cor, dist_hist_incor,
                 result_q, sam_ofh=None, ival=10000):
        Thread.__init__(self)
        self.args = args
        self.sam_q = sam_q
        self.sam_ofh = sam_ofh
        self.dataset = dataset
        self.dists = dists
        self.ref = ref
        self.alignment_class = alignment_class
        self.sc_diffs = defaultdict(int)
        self.ival = ival
        self.dist_hist_cor = dist_hist_cor
        self.dist_hist_incor = dist_hist_incor
        self.result_q = result_q
        self.deamon = True
        self.typ_hist = defaultdict(int)
    
    def run(self):
        """ Collect SAM output """
        args = self.args
        try:
            last_al, last_correct = None, None
            nal, nunp, npair = 0, 0, 0
            n_mate_first, n_mate_second = 0, 0
            # Following loop involves maintaining 'last_al' across
            # iterations so that we can match up the two ends of a pair
            correct = None
            while True:
                try:
                    ln = self.sam_q.get(block=True, timeout=0.5)
                    if ln is None:
                        break
                except Empty:
                    continue  # keep trying
                
                # Send SAM to SAM output filehandle
                if self.sam_ofh is not None:
                    self.sam_ofh.write(ln)
                
                # Skip headers
                if ln[0] == '@':
                    continue
                
                # Parse SAM record
                al = self.alignment_class()
                al.parse(ln)
                nal += 1
                if al.paired:
                    npair += 1
                else:
                    nunp += 1
                if (nal % self.ival) == 0:
                    logging.info('      # alignments parsed: %d (%d paired, %d unpaired)' % (nal, npair, nunp))

                if al.paired and last_al is not None:
                    assert al.concordant == last_al.concordant
                    assert al.discordant == last_al.discordant
                    mate1, mate2 = al, last_al
                    assert mate1.mate1 != mate2.mate1
                    if mate2.mate1:
                        mate1, mate2 = mate2, mate1
                    # if only one end aligned, possibly swap so that 'al' is
                    # the one that aligned
                    if not al.is_aligned() and last_al.is_aligned():
                        al, last_al = last_al, al
                        correct, last_correct = last_correct, correct
                else:
                    mate1, mate2 = None, None

                # Add the alignment to Dataset
                if self.dataset is not None:
                    is_training = (al.name[0] == '!' and al.name.startswith('!!ts!!'))
                    if al.is_aligned() and is_training:
                        _, refid, fw, refoff, sc, training_nm = al.name.split('!!ts-sep!!')
                        sc, refoff = int(sc), int(refoff)
                        self.sc_diffs[sc - al.bestScore] += 1
                        correct = False
                        # Check reference id, orientation
                        if refid == al.refid and fw == al.orientation():
                            # Check offset
                            correct = abs(refoff - al.pos) < args.wiggle
                        assert training_nm in ['unp', 'conc', 'disc', 'bad_end', 'bad_end2']
                        # Add to training dataset
                        if training_nm == 'unp':
                            self.typ_hist['unp'] += 1
                            self.dataset.add_unpaired(al, correct)
                        elif mate1 is not None:
                            n_mate_second += 1
                            correct1, correct2 = correct, last_correct
                            if last_al.mate1:
                                correct1, correct2 = correct2, correct1
                            if training_nm == 'conc':
                                self.typ_hist['conc'] += 1
                                if mate1.concordant:
                                    assert mate2.concordant
                                    self.dataset.add_concordant(mate1, mate2, correct1, correct2)
                            elif training_nm == 'disc':
                                self.typ_hist['disc'] += 1
                                if mate1.discordant:
                                    assert mate2.discordant
                                    self.dataset.add_discordant(mate1, mate2, correct1, correct2)
                            elif training_nm == 'bad_end':
                                self.typ_hist['bad_end'] += 1
                                self.dataset.add_bad_end(al, last_al, correct)
                            elif training_nm != 'bad_end2':
                                raise RuntimeError('Bad training data type: "%s"' % training_nm)
                        elif al.paired:
                            n_mate_first += 1

                    elif is_training:
                        pass

                    elif al.is_aligned():
                        # Test data
                        if mate1 is not None:
                            n_mate_second += 1
                            if not al.concordant and not al.discordant:
                                # bad-end
                                correct, dist = is_correct(al, args)
                                if not correct and dist < 20 * args.wiggle:
                                    self.dist_hist_incor[dist] += 1
                                elif correct:
                                    self.dist_hist_cor[dist] += 1
                                self.dataset.add_bad_end(al, last_al, correct)
                            else:
                                correct1, dist1 = is_correct(mate1, args)
                                correct2, dist2 = is_correct(mate2, args)
                                if not correct1 and dist1 < 20 * args.wiggle:
                                    self.dist_hist_incor[dist1] += 1
                                elif correct1:
                                    self.dist_hist_cor[dist1] += 1
                                if not correct2 and dist2 < 20 * args.wiggle:
                                    self.dist_hist_incor[dist2] += 1
                                elif correct2:
                                    self.dist_hist_cor[dist2] += 1
                                if al.concordant:
                                    self.dataset.add_concordant(mate1, mate2, correct1, correct2)
                                elif al.discordant:
                                    self.dataset.add_discordant(mate1, mate2, correct1, correct2)
                        elif not al.paired:
                            # For unpaired reads
                            correct, dist = is_correct(al, args)
                            if not correct and dist < 20 * args.wiggle:
                                self.dist_hist_incor[dist] += 1
                            elif correct:
                                self.dist_hist_cor[dist] += 1
                            self.dataset.add_unpaired(al, correct)
                        else:
                            n_mate_first += 1

                # Add to our input-read summaries
                if self.dists is not None and al.is_aligned():
                    try:
                        if mate1 is not None:
                            if al.concordant:
                                self.dists.add_concordant_pair(mate1, mate2, self.ref)
                            elif al.discordant:
                                self.dists.add_discordant_pair(mate1, mate2, self.ref)
                            else:
                                self.dists.add_bad_end_read(al, len(last_al.seq), self.ref)
                        elif not al.paired:
                            self.dists.add_unpaired_read(al, self.ref)
                    except ReferenceOOB:
                        pass
                
                if mate1 is None and al.paired:
                    last_al, last_correct = al, correct
                else:
                    last_al, last_correct = None, None

            logging.debug('Output monitoring thread returning successfully')
            assert n_mate_first == n_mate_second, (n_mate_first, n_mate_second)
            self.result_q.put(True)
        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print_tb(exc_traceback)
            logging.error(str(e))
            logging.debug('Output monitoring thread returning error')
            self.result_q.put(False)


def create_input(args):
    """ Return an Input object that reads all user-provided input reads """
    assert args.fasta or args.fastq
    return Input(file_format="fastq" if args.fastq else "fasta", unp_fns=args.U,
                 m1_fns=args.m1, m2_fns=args.m2)


def go(args, aligner_args):
    """ Main driver for tandem simulator """

    random.seed(args.seed)
    np.random.seed(args.seed)
    
    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG if args.verbose else logging.INFO)

    # Create output directory if needed
    if not os.path.isdir(args.output_directory):
        try:
            os.makedirs(args.output_directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    
    # Construct command to invoke aligner
    aligner_class, alignment_class = Bowtie2, AlignmentBowtie2
    align_cmd = None
    if args.aligner == 'bowtie2':
        align_cmd = 'bowtie2 '
        if args.bt2_exe is not None:
            align_cmd = args.bt2_exe + " "
        align_cmd += ' '.join(aligner_args)
        align_cmd += ' --reorder --sam-no-qname-trunc --mapq-extra'
    elif args.aligner == 'bwa-mem':
        align_cmd = 'bwa mem '
        if args.bwa_exe is not None:
            align_cmd = args.bwa_exe + ' mem '
        align_cmd += ' '.join(aligner_args)
        aligner_class, alignment_class = BwaMem, AlignmentBwaMem
    elif args.aligner is not None:
        raise RuntimeError('Aligner not supported: "%s"' % args.aligner)

    st = time.clock()
    
    logging.info('Loading reference data')
    with ReferenceIndexed(args.ref) as ref:
    
        # ##################################################
        # ALIGN REAL DATA (or read in alignments from SAM)
        # ##################################################
        
        aligner = aligner_class(align_cmd, args.index)
        sam_fn = os.path.join(args.output_directory, 'input.sam')
        input_sam_fh = open(sam_fn, 'w') if (args.write_input_sam or args.write_all) else None
        test_data = Dataset()
        dists = Dists(even=args.even)
        cor_dist, incor_dist = defaultdict(int), defaultdict(int)
        
        logging.info('Real-data aligner command: "%s"' % align_cmd)
        
        # Create the thread that eavesdrops on output from aligner
        result_test_q = Queue()
        othread = AlignmentReader(
            args,
            aligner.outQ,         # SAM queue
            test_data,            # Dataset to gather alignments into
            dists,                # empirical dists
            ref,                  # reference genome
            alignment_class,       #
            cor_dist,             #
            incor_dist,           #
            result_test_q,        # result queue
            sam_ofh=input_sam_fh)  # SAM output filehandle
        othread.daemon = True
        othread.start()
        
        # Stop timing setup
        setup_ival = time.clock() - st
        st = time.clock()
        
        logging.info('Initializing threads, queues and FIFOs')
        
        # Read through all the input read files and direct all reads to the
        # appropriate queue
        upto = args.upto or sys.maxint
        num_reads, num_unp, num_pair = 0, 0, 0
        for (rd1, rd2) in iter(create_input(args)):
            aligner.put(rd1, rd2)
            num_reads += 1
            if rd2 is None:
                num_unp += 1
            else:
                num_pair += 1
            if num_reads >= upto:
                break
            if not othread.is_alive():
                logging.error('Output thread died! - leaving input loop early.')
                break
        aligner.done()
        
        logging.info('Fed %d input reads (%d unpaired, %d pairs) to aligner' % (num_reads, num_unp, num_pair))
        logging.debug('Waiting for aligner to finish')
        while aligner.pipe.poll() is None:
            time.sleep(0.5)
        while othread.is_alive():
            othread.join(0.5)
        logging.debug('Aligner process and output thread finished')
        othread_result = result_test_q.get()
        if not othread_result:
            raise RuntimeError('Aligner output monitoring thread encountered error')
        
        logging.info('Finished aligning input reads')
        if input_sam_fh is not None:
            logging.info('  Input read alignments written to "%s"' % sam_fn)
        
        # Writing test dataset
        if args.write_test_data or args.write_all:
            test_csv_fn_prefix = os.path.join(args.output_directory, 'test')
            test_data.save(test_csv_fn_prefix, compress=args.compress_output)
            logging.info('Test data (CSV format) written to "%s*"' % test_csv_fn_prefix)
        
        # Writing correct/incorrect distances
        if args.write_test_distances or args.write_all:
            for short_name, long_name, hist in [('cor', 'Correct', cor_dist), ('incor', 'Incorrect', incor_dist)]:
                test_dist_fn = os.path.join(args.output_directory, 'test_%s_dists.csv' % short_name)
                with open(test_dist_fn, 'w') as fh:
                    for k, v in sorted(hist.iteritems()):
                        fh.write('%d,%d\n' % (k, v))
                logging.info('%s-alignment distances written to "%s"' % (long_name, test_dist_fn))
        
        # Stop timing interval for alignment phase 1
        al1_ival = time.clock() - st
        st = time.clock()
        
        # ##################################################
        # ALIGN SIMULATED DATA
        # ##################################################

        # Construct sequence and quality simulators
        logging.info('  Finalizing distributions')
        dists.finalize()

        # If the training data is all unpaired, or if the aligner
        # allows us to provide a mix a paired-end and unpaired data,
        # then we always align the training data with one aligner
        # process.  Otherwise, we need two aligner processes; one for
        # the paired-end and one for the unpaired data.
        
        iters = (1 if dists.has_unpaired_reads() else 0) + (1 if dists.has_pairs() else 0)
        if iters == 2 and aligner.supportsMix():
            iters = 1
        if iters == 2:
            logging.info('Aligner does not accept unpaired/paired mix; training will have 2 rounds')

        training_data = Dataset()
        for paired, lab in [(True, 'paired-end'), (False, 'unpaired')]:
            both = False
            if paired and not dists.has_pairs():
                continue
            if not paired and not dists.has_unpaired_reads():
                continue
            if aligner.supportsMix() and dists.has_pairs() and dists.has_unpaired_reads():
                # Do both unpaired and paired simualted reads in one round
                both, lab = True, 'both paired-end and unpaired'
            
            simw = FragmentSimSerial2(args.ref, dists)

            logging.info('Opening aligner process')
            aligner = aligner_class(align_cmd, args.index, pairsOnly=paired)
            sam_fn = os.path.join(args.output_directory, 'training.sam')
            training_sam_fh = open(sam_fn, 'w') if (args.write_training_sam or args.write_all) else None
            cor_dist, incor_dist = defaultdict(int), defaultdict(int)
            
            # Create thread that eavesdrops on output from aligner with simulated input
            logging.info('  Opening output-parsing thread (%s)' % lab)
            result_training_q = Queue()
            othread = AlignmentReader(
                args,
                aligner.outQ,
                training_data,
                None,
                ref,
                alignment_class,
                cor_dist,
                incor_dist,
                result_training_q,
                sam_ofh=training_sam_fh)
            othread.daemon = True
            othread.start()
            assert othread.is_alive()
            
            # Simulate concordant paired-end reads
            write_training_reads = args.write_training_reads or args.write_all
            training_out_fn, training_out_fh = {}, {}
            if write_training_reads:
                types = []
                if paired or both:
                    types.extend(['conc', 'disc', 'bad_end'])
                if not paired or both:
                    types.append('unp')
                for typ in types:
                    training_out_fn[typ] = fn = os.path.join(args.output_directory, 'training_%s.tab6' % typ)
                    training_out_fh[typ] = open(fn, 'w')

            logging.info('  Simulating reads')
            n_simread, typ_count = 0, defaultdict(int)
            for typ, rdp1, rdp2 in simw.simulate_batch(args.simulation_fraction, args.simulation_minimum):
                if write_training_reads:
                    training_out_fh[typ].write(Read.to_tab6(rdp1, rdp2) + '\n')
                aligner.put(rdp1, rdp2)
                typ_count[typ] += 1
                n_simread += 1
                if (n_simread % 1000) == 0:
                    logging.info('    simulated %d reads (%d conc, %d disc, %d bad-end, %d unp)' %
                                 (n_simread, typ_count['conc'], typ_count['disc'], typ_count['bad_end'], typ_count['unp']))

            logging.info('  Finished simulating reads (%d conc, %d disc, %d bad_end, %d unp)' %
                         (typ_count['conc'], typ_count['disc'], typ_count['bad_end'], typ_count['unp']))
            aligner.done()
            
            logging.debug('Waiting for aligner (%s) to finish' % lab)
            while aligner.pipe.poll() is None:
                time.sleep(0.5)
            while othread.is_alive():
                othread.join(0.5)
            logging.debug('aligner process finished')
            othread_result = result_training_q.get()
            if not othread_result:
                raise RuntimeError('Aligner output monitoring thread encountered error')
            logging.info('Finished simulating and aligning training reads')

            # Check the fraction of simualted reads that were aligned and
            # where we got an alignment of the expected type
            logging.info('Tally of how many of each simualted type aligned as that type:')
            for typ, cnt in typ_count.iteritems():
                if cnt == 0:
                    continue
                logging.info('  %s: simulated:%d aligned:%d (%0.2f%%)' %
                             (typ, cnt, othread.typ_hist[typ], 100.0 * othread.typ_hist[typ] / cnt))

            for typ in training_out_fh.iterkeys():
                training_out_fh[typ].close()
                logging.info('  Training reads written to "%s"' % training_out_fn[typ])

            if both:
                break
        
        # Stop timing interval for alignment phase 2
        al2_ival = time.clock() - st
    
        logging.info('Score difference (expected - actual) histogram:')
        for k, v in sorted(othread.sc_diffs.iteritems()):
            logging.info('  %d: %d' % (k, v))
        
        if args.write_training_data or args.write_all:
            # Writing training data
            training_csv_fn_prefix = os.path.join(args.output_directory, 'training')
            training_data.save(training_csv_fn_prefix, compress=args.compress_output)
            logging.info('Training data (CSV format) written to "%s*"' % training_csv_fn_prefix)
        
        if setup_ival is not None:
            logging.info('Setup running time: %0.3f secs' % setup_ival)
        if al1_ival is not None:
            logging.info('Alignment (real reads): %0.3f secs' % al1_ival)
        if al2_ival is not None:
            logging.info('Alignment (simulated reads): %0.3f secs' % al2_ival)

if __name__ == "__main__":
    
    import argparse

    parser = argparse.ArgumentParser(
        description='Evaluate the sensitivity of an alignment tool w/r/t given genome and set of alignment parameters.')

    # Inputs
    parser.add_argument('--ref', metavar='path', type=str, nargs='+', required=True,
                        help='FASTA file(s) containing reference genome sequences')
    parser.add_argument('--pickle-ref', metavar='path', type=str,
                        help='Pickle FASTA input for speed, or use pickled version if it exists already.  Pickled '
                             'version is stored at given path')
    parser.add_argument('--alignability-wig', metavar='path', type=str, nargs='+',
                        help='.wig files with alignability info')
    parser.add_argument('--U', metavar='path', type=str, nargs='+', help='Unpaired read files')
    parser.add_argument('--m1', metavar='path', type=str, nargs='+',
                        help='Mate 1 files; must be specified in same order as --m2')
    parser.add_argument('--m2', metavar='path', type=str, nargs='+',
                        help='Mate 2 files; must be specified in same order as --m1')
    parser.add_argument('--fastq', action='store_const', const=True, default=True, help='Input reads are FASTQ')
    parser.add_argument('--fasta', action='store_const', const=True, default=False, help='Input reads are FASTA')
    parser.add_argument('--index', metavar='path', type=str, help='Index file to use (usually a prefix).')

    parser.add_argument('--seed', metavar='int', type=int, default=99099, required=False,
                        help='Integer to initialize pseudo-random generator')
    parser.add_argument('--simulation-fraction', metavar='fraction', type=float, default=0.1, required=False,
                        help='When determining the number of simulated reads to generate for each type of '
                             'alignment (concordant, discordant, unpaired), let it be no less than this '
                             'fraction times the number of alignment of that type in the input data.')
    parser.add_argument('--simulation-minimum', metavar='int', type=int, default=1000, required=False,
                        help='When determining the number of simulated reads to generate for each type of '
                             'alignment (concordant, discordant, unpaired) that appears in the input, let '
                             'it be no less than this number.')
    parser.add_argument('--upto', metavar='int', type=int, default=None, required=False,
                        help='Stop after this many input reads')
    parser.add_argument('--wiggle', metavar='int', type=int, default=30, required=False,
                        help='Wiggle room to allow in starting position when determining whether alignment is correct')

    parser.add_argument('--sam-input', metavar='path', type=str,
                        help='Input SAM file to apply training data to.  Use with --training-input.')
    parser.add_argument('--bt2-exe', metavar='path', type=str, help='Path to Bowtie 2 exe')
    parser.add_argument('--bwa-exe', metavar='path', type=str, help='Path to BWA exe')
    parser.add_argument('--aligner', metavar='name', default='bowtie2', type=str, help='bowtie2 or bwa-mem')

    # Some basic flags
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')
    parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
    parser.add_argument('--version', action='store_const', const=True, default=False, help='Print version and quit')

    # For when input is itself simulated, so we can output a Dataset with the
    # 'correct' column filled in properly 
    parser.add_argument('--input-from-mason', action='store_const', const=True, default=False,
                        help='Input reads were simulated from Mason')
    parser.add_argument('--input-from-wgsim', action='store_const', const=True, default=False,
                        help='Input reads were simulated from wgsim')
    parser.add_argument('--input-from-grinder', action='store_const', const=True, default=False,
                        help='Input reads were simulated from Grinder')
    parser.add_argument('--correct-chromosomes', metavar='list', type=str, nargs='+',
                        help='Label test data originating from any of these chromosomes as "correct."  Useful for '
                             'tests on real-world data where it is known that the data came from a parituclar '
                             'chromosome.')

    # Sampling-related arguments
    parser.add_argument('--even', action='store_const', const=True, default=False,
                        help="Don't try to even out sampling.  This usually leads to a worse fit since not as many "
                             "alignment scores and degrees of repetitiveness are well represented in the sample.")

    # Output file-related arguments
    parser.add_argument('--output-directory', metavar='path', type=str, required=True,
                        help='Write outputs to this directory')
    parser.add_argument('--write-input-sam', action='store_const', const=True, default=False,
                        help='Write SAM alignments for the real input reads to "input.sam" in output directory')
    parser.add_argument('--write-training-reads', action='store_const', const=True, default=False,
                        help='Write FASTQ for the training reads to "training.fastq" in output directory')
    parser.add_argument('--write-test-data', action='store_const', const=True, default=False,
                        help='Write Dataset object for training data.  "Correct" column set to all None\'s.')
    parser.add_argument('--write-training-sam', action='store_const', const=True, default=False,
                        help='Write SAM alignments for the training reads to "training.sam" in output directory')
    parser.add_argument('--write-training-data', action='store_const', const=True, default=False,
                        help='Write Dataset object for training data.')
    parser.add_argument('--write-test-distances', action='store_const', const=True, default=False,
                        help='Write distances between true/actual alignments.')
    parser.add_argument('--write-training-distances', action='store_const', const=True, default=False,
                        help='Write distances between true/actual alignments.')
    parser.add_argument('--write-all', action='store_const', const=True, default=False,
                        help='Same as specifying all --write-* options')
    parser.add_argument('--compress-output', action='store_const', const=True, default=False,
                        help='gzip all output files')

    argv = sys.argv
    outer_aligner_args = []
    in_args = False
    for i in xrange(1, len(sys.argv)):
        if in_args:
            outer_aligner_args.append(sys.argv[i])
        if sys.argv[i] == '--':
            argv = sys.argv[:i]
            in_args = True
    
    outer_args = parser.parse_args(argv[1:])
    
    if outer_args.test:
        
        import unittest
        
        class Test(unittest.TestCase):
            def test_1(self):
                dst = Dist()
                dst.add("CP")
                dst.add("CP")
                dst.add("DP")
                dr = dst.draw()
                self.assertTrue(dr == "CP" or dr == "DP")

        unittest.main(argv=[sys.argv[0]])
    elif outer_args.profile:
        import cProfile
        cProfile.run('go(outer_args, outer_aligner_args)')
    else:
        go(outer_args, outer_aligner_args)
