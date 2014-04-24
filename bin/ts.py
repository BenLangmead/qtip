#!/usr/bin/env python

"""
ts.py

A "tandem simulator," which wraps an alignment tool as it runs, eavesdrops on
input and output, simulates a new dataset similar to the input data, aligns
it, uses those alignments as training data to build a model to predict MAPQ,
then re-calcualtes MAPQs for the original input using that predictor.

Output files encapsulate:
1. Input data model
2. Simulated reads
3. Alignments for simulated reads
4. (3) converted into training-data records
5. Trained models
6. Results of running the trained models on the training data

Things we learn from reads
==========================

- Read length distribution

Things we learn from alignments
===============================

- Alignment type (aligned, unaligned, concordant, discordant)
- Fragment length distribution
- Number and placement of mismatches and gaps and corresponding quality values

SAM extra fields used
=====================

 Normal: AS:i, XS:i, MD:Z
 
 Bowtie-specific: Xs:i, YT:Z, YS:i, Zp:i, YN:i, Yn:i
 (and we potentially don't need YS:i or YT:Z?)

"""

__author__ = "Ben Langmead"
__email__ = "langmea@cs.jhu.edu"

import os
import sys
from string import maketrans
import random
import time
import logging
import gzip
import errno
import tempfile
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
from simplesim import FragmentSimSerial2
from read import Read, Alignment
from bowtie2 import AlignmentBowtie2, Bowtie2
from bwamem import AlignmentBwaMem, BwaMem
from reference import ReferenceIndexed, ReferenceOOB

VERSION = '0.1.0'

_revcomp_trans = maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s[::-1].translate(_revcomp_trans)


class ReadTemplate(object):
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
    def __init__(self, rt1, rt2, fraglen, upstream1):
        self.rt1 = rt1
        self.rt2 = rt2
        self._fraglen = fraglen
        self.upstream1 = upstream1

    @property
    def fraglen(self):
        return self._fraglen


class ScoreDist(object):
    """ Capture a list of tuples, where each tuples represents the following
        traits observed in an alignment: (a) orientation, (b) quality string,
        (c) read side of stacked alignment, (d) ref side of stacked alignment,
        (e) score.  Keeping these tuples allows us to generate new reads that
        mimic observed reads in these key ways. """
    
    def __init__(self, k=50000):
        """ Make a reservoir sampler for holding the tuples. """
        self.res = ReservoirSampler(k)
        self.max_fraglen = 0
        self.avg_fraglen = None
        self.finalized = False
        self.samp = None

    def draw(self, bias=2):
        """ Draw from the reservoir """
        assert self.finalized
        assert not self.empty()
        rnd_i = int(len(self.samp) * random.uniform(0, 1) / random.uniform(1, bias))
        assert 0 <= rnd_i < len(self.samp)
        sc, _, fw, qual, rd_aln, rf_aln, rf_len, mate1, olen = self.samp[rnd_i]
        return ReadTemplate(sc, fw, qual, rd_aln, rf_aln, rf_len, mate1, olen)
    
    def add(self, al, ref, ordlen=0):
        """ Convert given alignment to a tuple and add it to the reservoir
            sampler. """
        sc = al.bestScore
        # Get stacked alignment
        rd_aln, rf_aln, rd_len, rf_len = al.stacked_alignment(align_soft_clipped=True, ref=ref)
        self.max_fraglen = max(self.max_fraglen, rf_len)
        self.res.add((sc, random.uniform(0, 1),
                      al.fw, al.qual, rd_aln, rf_aln, rf_len, al.mate1, ordlen))

    @property
    def num_added(self):
        """ Return the number of tuples that have been added, which could be
            much greater than the number sampled """
        return self.res.num_added()

    @property
    def num_sampled(self):
        """ Return number of tuples that have been sampled """
        return self.res.num_sampled()
    
    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_sampled == 0

    def finalize(self):
        """ Sort the samples in preparation for draws that might be biased
            toward high or (more likely) low scores. """
        self.finalized = True
        if not self.empty():
            self.samp = sorted(self.res.r)
            self.avg_fraglen = np.mean([x[6] for x in self.res])


class ScorePairDist(object):
    """ Capture a list of tuple pairs, where each tuple pair represents the
        following in an alignment of a paired-end reads where both ends
        aligned: (a) orientation, (b) quality string, (c) read side of stacked
        alignment, (d) ref side of stacked alignment, (e) score, (f) fragment
        length.  We record pairs of tuples where the first element of the pair
        corresponds to mate 1 and the second to mate 2. """
    
    def __init__(self, k=50000, max_allowed_fraglen=100000):
        """ Make a reservoir sampler for holding the tuples. """
        self.res = ReservoirSampler(k)
        self.max_allowed_fraglen = max_allowed_fraglen
        self.max_fraglen = 0  # maximum observed fragment length
        self.avg_fraglen = None
        self.finalized = False
        self.samp = None

    def draw(self, bias=2):
        """ Draw from the reservoir """
        assert self.finalized
        assert not self.empty()
        rnd_i = int(len(self.samp) * random.uniform(0, 1) / random.uniform(1, bias))
        assert 0 <= rnd_i < len(self.samp)
        tup1, tup2, tup3, fraglen, upstream1 = self.samp[rnd_i]
        fw_1, qual_1, rd_aln_1, rf_aln_1, sc_1, rf_len_1, mate1_1, _ = tup2
        fw_2, qual_2, rd_aln_2, rf_aln_2, sc_2, rf_len_2, mate1_2, _ = tup3
        return PairTemplate(ReadTemplate(sc_1, fw_1, qual_1, rd_aln_1, rf_aln_1, rf_len_1, mate1_1, rf_len_2),
                            ReadTemplate(sc_2, fw_2, qual_2, rd_aln_2, rf_aln_2, rf_len_2, mate1_2, rf_len_1),
                            fraglen, upstream1)

    def add(self, al1, al2, ref):
        """ Convert given alignment pair to a tuple and add it to the
            reservoir sampler. """
        sc1, sc2 = al1.bestScore, al2.bestScore
        # Make note of fragment length
        fraglen = Alignment.fragment_length(al1, al2)
        fraglen = min(fraglen, self.max_allowed_fraglen)
        self.max_fraglen = max(self.max_fraglen, fraglen)
        # Make note of which end is upstream
        upstream1 = al1.pos < al2.pos
        # Get stacked alignment
        rd_aln_1, rf_aln_1, rd_len_1, rf_len_1 = al1.stacked_alignment(align_soft_clipped=True, ref=ref)
        rd_aln_2, rf_aln_2, rd_len_2, rf_len_2 = al2.stacked_alignment(align_soft_clipped=True, ref=ref)
        assert fraglen == 0 or max(rf_len_1, rf_len_2) <= fraglen
        self.res.add(((sc1 + sc2, random.uniform(0, 1)),
                      (al1.fw, al1.qual, rd_aln_1, rf_aln_1, sc1, rf_len_1, True, rf_len_2),
                      (al2.fw, al2.qual, rd_aln_2, rf_aln_2, sc2, rf_len_2, False, rf_len_1),
                      fraglen,
                      upstream1))

    @property
    def num_added(self):
        """ Return the number of tuples that have been added, which could be
            much greater than the number sampled """
        return self.res.num_added()

    @property
    def num_sampled(self):
        """ Return number of tuples that have been sampled """
        return self.res.num_sampled()

    def empty(self):
        """ Return true iff no tuples have been added """
        return self.num_sampled == 0

    def finalize(self):
        """ Sort the samples in preparation for draws that might be biased
            toward high or (more likely) low scores. """
        self.finalized = True
        if not self.empty():
            self.samp = sorted(self.res.r)
            self.avg_fraglen = np.mean([x[3] for x in self.res])


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


class Dists(object):
    
    """ Encapsulates the distributions that we capture from real data.
        We collect random subsets of qualities/edits or unpaired reads
        and same for paired-end mate 1s and mate 2s.  We also collect
        data on concordant and discordantly aligned pairs, such as
        their fragment length and strands. """

    def __init__(self, k=50000, max_allowed_fraglen=100000):
        self.sc_dist_unp = ScoreDist(k=k)
        self.sc_dist_bad_end = ScoreDist(k=k)
        self.sc_dist_conc = ScorePairDist(k=k, max_allowed_fraglen=max_allowed_fraglen)
        self.sc_dist_disc = ScorePairDist(k=k, max_allowed_fraglen=max_allowed_fraglen)

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
        al_refid_trimmed = al.refid.split()[0]
        return (al_refid_trimmed in args.correct_chromosomes), 0
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
                 result_q, sam_ofh=None, ival=20000):
        Thread.__init__(self)
        self.args = args
        self.sam_q = sam_q
        self.sam_ofh = sam_ofh
        self.dataset = dataset
        self.dists = dists
        self.ref = ref  # TODO: what kind of reference object?  Needs to be random access?  Needs to be fast?
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
        has_readline = hasattr(self.sam_q, 'read')
        try:
            last_al, last_correct = None, None
            nal, nunp, nignored, npair = 0, 0, 0, 0
            #n_mate_first, n_mate_second = 0, 0
            # Following loop involves maintaining 'last_al' across
            # iterations so that we can match up the two ends of a pair
            correct = None
            while True:
                try:
                    if has_readline:
                        ln = self.sam_q.readline()
                        if len(ln) == 0:
                            break
                    else:
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
                if al.flags >= 2048:
                    nignored += 1
                    continue
                elif al.paired:
                    npair += 1
                else:
                    nunp += 1

                if (nal % self.ival) == 0:
                    logging.info('      # alignments parsed: %d (%d paired, %d unpaired, %d ignored)' %
                                 (nal, npair, nunp, nignored))

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
                            #n_mate_second += 1
                            #assert n_mate_second == n_mate_first, (n_mate_first, n_mate_second)
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
                        #elif al.paired:
                        #    n_mate_first += 1

                    elif is_training:
                        pass

                    elif al.is_aligned():
                        # Test data
                        if mate1 is not None:
                            #n_mate_second += 1
                            #assert n_mate_second == n_mate_first, (n_mate_first, n_mate_second)
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
                        #else:
                        #    n_mate_first += 1

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
            #assert n_mate_first == n_mate_second, (n_mate_first, n_mate_second)
            self.result_q.put(True)
        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print_tb(exc_traceback)
            logging.error(str(e))
            logging.debug('Output monitoring thread returning error')
            self.result_q.put(False)


class TemporaryFileManager():

    def __init__(self):
        self.dir = tempfile.mkdtemp()
        self.files = set()
        self.groups = defaultdict(list)
        self.peak_size = 0

    def get_filename(self, fn_basename, group=''):
        """ Return filename for new temporary file in temp dir """
        if fn_basename in self.files:
            raise RuntimeError('Temporary file with name "%s" already exists' % fn_basename)
        self.groups[group].append(fn_basename)
        self.files.add(fn_basename)
        return os.path.join(self.dir, fn_basename)

    def remove_group(self, group):
        """ Remove all the temporary files belonging to the named group """
        for fn_basename in self.groups[group]:
            self.files.remove(fn_basename)
            os.remove(os.path.join(self.dir, fn_basename))
        del self.groups[group]

    def size(self):
        """ Return total size of all the files in the temp dir """
        return sum(os.path.getsize(os.path.join(self.dir, f)) for f in self.files)

    def update_peak(self):
        """ Update peak size of temporary files """
        self.peak_size = max(self.peak_size, self.size())


class Timing(object):

    def __init__(self):
        self.labs = []
        self.timers = dict()

    def start_timer(self, lab):
        self.labs.append(lab)
        self.timers[lab] = time.time()

    def end_timer(self, lab):
        self.timers[lab] = time.time() - self.timers[lab]

    def __str__(self):
        ret = []
        for lab in self.labs:
            ret.append('\t'.join([lab, str(self.timers[lab])]))
        return '\n'.join(ret) + '\n'


def go(args, aligner_args):
    """ Main driver for tandem simulator """

    random.seed(args.seed)
    np.random.seed(args.seed)

    tim = Timing()
    tim.start_timer('Overall')

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

    if args.U is not None and args.m1 is not None:
        raise RuntimeError('Input must consist of only unpaired or only paired-end reads')
    
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

    temp_man = TemporaryFileManager()

    logging.info('Loading reference data')
    with ReferenceIndexed(args.ref) as ref:
    
        # ##################################################
        # ALIGN REAL DATA (or read in alignments from SAM)
        # ##################################################

        unpaired_arg = args.U
        paired_arg = None
        if args.m1 is not None:
            paired_arg = zip(args.m1, args.m2)

        input_sam_fh = None
        write_input_sam = args.write_input_sam or args.write_all
        sam_arg = None
        sam_fn = os.path.join(args.output_directory, 'input.sam')
        if not write_input_sam:
            sam_fn = temp_man.get_filename('input.sam', 'input sam')

        if args.use_temporary_files:
            sam_arg = sam_fn
        else:
            input_sam_fh = open(sam_fn, 'w') if write_input_sam else None

        if args.use_temporary_files:
            tim.start_timer('Aligning input reads')

        logging.info('Command for aligning input data: "%s"' % align_cmd)
        aligner = aligner_class(align_cmd, args.index,
                                unpaired=unpaired_arg, paired=paired_arg,
                                sam=sam_arg)

        if args.use_temporary_files:
            logging.debug('Waiting for aligner to finish')
            while aligner.pipe.poll() is None:
                time.sleep(0.5)
            logging.debug('Aligner finished')
            tim.end_timer('Aligning input reads')

        test_data = Dataset()
        dists = Dists(args.max_allowed_fraglen)
        cor_dist, incor_dist = defaultdict(int), defaultdict(int)
        
        result_test_q = Queue()
        if args.use_temporary_files:
            tim.start_timer('Parsing input read alignments')
            with open(sam_fn, 'r') as sam_fh:
                othread = AlignmentReader(
                    args,
                    sam_fh,                # SAM file
                    test_data,             # Dataset to gather alignments into
                    dists,                 # empirical dists
                    ref,                   # reference genome
                    alignment_class,       # class to instantiate for alignment
                    cor_dist,              # dist. of correct alignment deviations
                    incor_dist,            # dist. of incorrect alignment deviations
                    result_test_q,         # result queue
                    sam_ofh=input_sam_fh)  # SAM output filehandle
                othread.run()

            temp_man.update_peak()
            temp_man.remove_group('input sam')
            tim.end_timer('Parsing input read alignments')
        else:
            # Create the thread that eavesdrops on output from aligner
            othread = AlignmentReader(
                args,
                aligner.outQ,          # SAM queue
                test_data,             # Dataset to gather alignments into
                dists,                 # empirical dists
                ref,                   # reference genome
                alignment_class,       # class to instantiate for alignment
                cor_dist,              # dist. of correct alignment deviations
                incor_dist,            # dist. of incorrect alignment deviations
                result_test_q,         # result queue
                sam_ofh=input_sam_fh)  # SAM output filehandle
            othread.daemon = True
            othread.start()
            logging.info('Initializing threads, queues and FIFOs')

            logging.debug('Waiting for wrapped aligner to finish')
            while aligner.pipe.poll() is None:
                time.sleep(0.5)
            while othread.is_alive():
                othread.join(0.5)
            logging.debug('Aligner process and output thread finished')

        othread_result = result_test_q.get()
        if not othread_result:
            raise RuntimeError('Aligner output monitoring thread encountered error')

        if input_sam_fh is not None:
            logging.info('  Input alignments written to "%s"' % sam_fn)
        
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

        # ##################################################
        # ALIGN SIMULATED DATA
        # ##################################################

        # Construct sequence and quality simulators
        logging.info('  Finalizing distributions')
        dists.finalize()
        logging.info('    Longest unpaired=%d, fragment=%d' % (dists.longest_unpaired(), dists.longest_fragment()))
        # TODO: print something about average length and average alignment score

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

            def simulate(simw, aligner=None, format='tab6'):
                """ Need to re-think this to accomodate the case where we're
                    writing to a file instead of directly to aligner """
                write_training_reads = args.write_training_reads or args.write_all
                training_out_fn, training_out_fh = {}, {}
                if write_training_reads or args.use_temporary_files:
                    types = []
                    if paired or both:
                        types.extend(['conc', 'disc', 'bad_end'])
                    if not paired or both:
                        types.append('unp')
                    for t in types:
                        fn_base = 'training_%s.%s' % (t, format)
                        fn = os.path.join(args.output_directory, fn_base)
                        if not write_training_reads:
                            fn = temp_man.get_filename(fn_base, 'tandem reads')
                        training_out_fn[t] = fn
                        training_out_fh[t] = open(fn, 'w')

                logging.info('  Simulating reads')
                n_simread, typ_count = 0, defaultdict(int)
                for t, rdp1, rdp2 in simw.simulate_batch(args.sim_fraction, args.sim_unp_min,
                                                         args.sim_conc_min, args.sim_disc_min,
                                                         args.sim_bad_end_min,
                                                         bias=args.low_score_bias):
                    if t in training_out_fh:
                        if format == 'tab6':
                            training_out_fh[t].write(Read.to_tab6(rdp1, rdp2) + '\n')
                        elif format == 'interleaved_fastq':
                            training_out_fh[t].write(Read.to_interleaved_fastq(rdp1, rdp2) + '\n')
                        else:
                            raise RuntimeError('Bad training read output format "%s"' % format)
                    if aligner is not None:
                        aligner.put(rdp1, rdp2)
                    typ_count[t] += 1
                    n_simread += 1
                    if (n_simread % 20000) == 0:
                        logging.info('    simulated %d reads (%d conc, %d disc, %d bad-end, %d unp)' %
                                     (n_simread, typ_count['conc'], typ_count['disc'],
                                      typ_count['bad_end'], typ_count['unp']))

                for t in training_out_fh.iterkeys():
                    training_out_fh[t].close()
                    logging.info('  Training reads written to "%s"' % training_out_fn[t])

                logging.info('  Finished simulating reads (%d conc, %d disc, %d bad_end, %d unp)' %
                             (typ_count['conc'], typ_count['disc'], typ_count['bad_end'], typ_count['unp']))

                return typ_count, training_out_fn

            simw = FragmentSimSerial2(args.ref, dists)

            if args.use_temporary_files:
                tim.start_timer('Simulating tandem reads')
                logging.info('Simulating tandem reads (%s)' % lab)
                preferred_format = aligner.preferred_paired_format() if (paired or both) \
                    else aligner.preferred_unpaired_format()
                typ_sim_count, training_out_fn = simulate(simw, format=preferred_format)
                logging.info('Finished simulating tandem reads')
                tim.end_timer('Simulating tandem reads')

                logging.info('Aligning tandem reads (%s)' % lab)
                sam_fn_base = 'training.sam'
                if paired and not both:
                    sam_fn_base = 'training_paired.sam'
                elif not both:
                    sam_fn_base = 'training_unpaired.sam'

                if args.write_training_sam:
                    sam_fn = os.path.join(args.output_directory, sam_fn_base)
                else:
                    sam_fn = temp_man.get_filename(sam_fn_base, 'tandem sam')

                unpaired_arg = [training_out_fn['unp']] if 'unp' in training_out_fn else None
                paired_combined_arg = None
                if 'conc' in training_out_fn or 'disc' in training_out_fn:
                    paired_combined_arg = []
                    for t in ['conc', 'disc']:
                        if t in training_out_fn:
                            paired_combined_arg.append(training_out_fn[t])

                aligner = aligner_class(align_cmd, args.index,
                                        unpaired=unpaired_arg, paired_combined=paired_combined_arg,
                                        sam=sam_fn, format=preferred_format)
                cor_dist, incor_dist = defaultdict(int), defaultdict(int)

                while aligner.pipe.poll() is None:
                    time.sleep(0.5)
                while othread.is_alive():
                    othread.join(0.5)
                logging.debug('Finished aligning tandem reads')

                # remove temporary reads
                temp_man.update_peak()
                temp_man.remove_group('tandem reads')

                logging.info('Parsing tandem alignments (%s)' % lab)
                tim.start_timer('Parsing tandem alignments')
                with open(sam_fn, 'r') as sam_fh:
                    result_training_q = Queue()
                    reader = AlignmentReader(
                        args,
                        sam_fh,
                        training_data,
                        None,
                        ref,
                        alignment_class,
                        cor_dist,
                        incor_dist,
                        result_training_q)
                    reader.run()
                    typ_align_count, sc_diffs = reader.typ_hist, reader.sc_diffs

                # remove temporary alignments
                temp_man.update_peak()
                temp_man.remove_group('tandem sam')

                othread_result = result_training_q.get()
                if not othread_result:
                    raise RuntimeError('Tandem alignment parser encountered error')
                logging.info('Finished parsing tandem alignments')
                tim.end_timer('Parsing tandem alignments')
            else:

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

                typ_sim_count, training_out_fn = simulate(simw, aligner, format='tab6')
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

                typ_align_count, sc_diffs = othread.typ_hist, othread.sc_diffs

            # Check the fraction of simualted reads that were aligned and
            # where we got an alignment of the expected type
            logging.info('Tally of how many of each simualted type aligned as that type:')
            for typ, cnt in typ_sim_count.iteritems():
                if cnt == 0:
                    continue
                if typ not in typ_align_count:
                    logging.warning('  %s: simulated:%d but ALIGNED NONE' % (typ, cnt))
                else:
                    func = logging.warning if (cnt / float(typ_align_count[typ])) < 0.3 else logging.info
                    func('  %s: simulated:%d aligned:%d (%0.2f%%)' %
                         (typ, cnt, typ_align_count[typ], 100.0 * typ_align_count[typ] / cnt))

            if both:
                break
        
    logging.info('Score difference (expected - actual) histogram:')
    for k, v in sorted(sc_diffs.iteritems()):
        logging.info('  %d: %d' % (k, v))

    tim.start_timer('Writing training data')
    if args.write_training_data or args.write_all:
        # Writing training data
        training_csv_fn_prefix = os.path.join(args.output_directory, 'training')
        training_data.save(training_csv_fn_prefix, compress=args.compress_output)
        logging.info('Training data (CSV format) written to "%s*"' % training_csv_fn_prefix)
    tim.end_timer('Writing training data')

    logging.info('Peak temporary file size %0.2fMB' % (temp_man.peak_size / (1024.0 * 1024)))

    tim.end_timer('Overall')
    for ln in str(tim).split('\n'):
        logging.info(ln)
    if args.write_timings or args.write_all:
        with open(os.path.join(args.output_directory, 'timing.tsv'), 'w') as fh:
            fh.write(str(tim))

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

    parser.add_argument('--sim-fraction', metavar='fraction', type=float, default=0.01, required=False,
                        help='When determining the number of simulated reads to generate for each type of '
                             'alignment (concordant, discordant, bad-end, unpaired), let it be no less '
                             'than this fraction times the number of alignment of that type in the input '
                             'data.')
    parser.add_argument('--sim-unp-min', metavar='int', type=int, default=30000, required=False,
                        help='Number of simulated unpaired reads will be no less than this number.')
    parser.add_argument('--sim-conc-min', metavar='int', type=int, default=30000, required=False,
                        help='Number of simulated concordant pairs will be no less than this number.')
    parser.add_argument('--sim-disc-min', metavar='int', type=int, default=10000, required=False,
                        help='Number of simulated discordant pairs will be no less than this number.')
    parser.add_argument('--sim-bad-end-min', metavar='int', type=int, default=10000, required=False,
                        help='Number of simulated pairs with-one-bad-end will be no less than this number.')

    parser.add_argument('--upto', metavar='int', type=int, default=None, required=False,
                        help='Stop after this many input reads')
    parser.add_argument('--max-allowed-fraglen', metavar='int', type=int, default=100000, required=False,
                        help='When simulating fragments, observed fragments longer than this will be'
                             'truncated to this length')
    parser.add_argument('--low-score-bias', metavar='float', type=int, default=7, required=False,
                        help='When simulating reads, we randomly select a real read\'s alignment profile'
                             'as a template.  A higher value for this parameter makes it more likely'
                             'we\'ll choose a low-scoring alignment.  If set to 1, all templates are'
                             'equally likely.')

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
    parser.add_argument('--write-timings', action='store_const', const=True, default=False,
                        help='Write timing info to "timing.tsv".')
    parser.add_argument('--write-all', action='store_const', const=True, default=False,
                        help='Same as specifying all --write-* options')
    parser.add_argument('--compress-output', action='store_const', const=True, default=False,
                        help='gzip all output files')

    # Resource usage
    parser.add_argument('--use-concurrency', action='store_const', const=True, default=False,
                        help='Use pipes instead of temporary files.  Reduces disk usage '
                             'but requires more memory.')
    parser.add_argument('--use-temporary-files', action='store_const', const=True, default=True,
                        help='Use temporary files instead of pipes.  Reduces peak memory footprint '
                             'but requires more disk space.')

    if '--version' in sys.argv:
        print 'Tandem simulator, version ' + VERSION
        sys.exit(0)

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
