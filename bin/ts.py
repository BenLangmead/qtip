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
from samples import DatasetOnDisk
from simplesim import FragmentSimSerial2
from read import Read
from bowtie2 import AlignmentBowtie2, Bowtie2
from bwamem import AlignmentBwaMem, BwaMem
from mosaik import AlignmentMosaik, Mosaik
from snap import AlignmentSnap, SnapAligner
from reference import ReferenceIndexed, ReferenceOOB
from tempman import TemporaryFileManager
from score_dists import CollapsedScoreDist, CollapsedScorePairDist

VERSION = '0.1.0'

_revcomp_trans = maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s[::-1].translate(_revcomp_trans)


class Dists(object):
    
    """ Encapsulates the distributions that we capture from real data.
        We collect random subsets of qualities/edits or unpaired reads
        and same for paired-end mate 1s and mate 2s.  We also collect
        data on concordant and discordantly aligned pairs, such as
        their fragment length and strands. """

    def __init__(self, max_allowed_fraglen=100000, fraction_even=0.5, bias=1.0):
        self.sc_dist_unp = CollapsedScoreDist(fraction_even=fraction_even, bias=bias)
        self.sc_dist_bad_end = CollapsedScoreDist(fraction_even=fraction_even, bias=bias)
        self.sc_dist_conc = CollapsedScorePairDist(max_allowed_fraglen=max_allowed_fraglen,
                                                   fraction_even=fraction_even, bias=bias)
        self.sc_dist_disc = CollapsedScorePairDist(max_allowed_fraglen=max_allowed_fraglen,
                                                   fraction_even=fraction_even, bias=bias)

    def finalize(self):
        self.sc_dist_unp.finalize()
        self.sc_dist_bad_end.finalize()
        self.sc_dist_conc.finalize()
        self.sc_dist_disc.finalize()

    def add_concordant_pair(self, al1, al2, correct1, correct2, ref):
        """ Add concordant paired-end read alignment to the model """
        self.sc_dist_conc.add(al1, al2, correct1, correct2, ref)

    def add_discordant_pair(self, al1, al2, correct1, correct2, ref):
        """ Add discordant paired-end read alignment to the model """
        self.sc_dist_disc.add(al1, al2, correct1, correct2, ref)

    def add_unpaired_read(self, al, correct, ref):
        """ Add unpaired read alignment to the model """
        self.sc_dist_unp.add(al, correct, ref)

    def add_bad_end_read(self, al, correct, ordlen, ref):
        """ Add bad-end read alignment to the model """
        self.sc_dist_bad_end.add(al, correct, ref, ordlen)

    def has_pairs(self):
        return not self.sc_dist_conc.empty() or not self.sc_dist_disc.empty() or not self.sc_dist_bad_end.empty()

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
    if args['correct_chromosomes'] is not None:
        al_refid_trimmed = al.refid.split()[0]
        return (al_refid_trimmed in args['correct_chromosomes']), 0
    elif simulators.isExtendedWgsim(al.name):
        return simulators.correctExtendedWgsim(al, args['wiggle'])
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
                        name_split = al.name.split('!!ts-sep!!')
                        if len(name_split) == 6:
                            # unpaired
                            _, refid, fw, refoff, sc, training_nm = name_split
                        else:
                            # Paired.  Both mates must have same name.
                            assert len(name_split) == 10
                            assert al.paired
                            _, refid1, fw1, refoff1, sc1, refid2, fw2, refoff2, sc2, training_nm = name_split
                            if al.mate1:
                                refid, fw, refoff, sc = refid1, fw1, refoff1, sc1
                            else:
                                refid, fw, refoff, sc = refid2, fw2, refoff2, sc2
                        sc, refoff = int(sc), int(refoff)
                        self.sc_diffs[sc - al.bestScore] += 1
                        correct = False
                        # Check reference id, orientation
                        if refid == al.refid and fw == al.orientation():
                            # Check offset
                            correct = abs(refoff - al.pos) < args['wiggle']
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
                                if not correct and dist < 20 * args['wiggle']:
                                    self.dist_hist_incor[dist] += 1
                                elif correct:
                                    self.dist_hist_cor[dist] += 1
                                self.dataset.add_bad_end(al, last_al, correct)
                            else:
                                correct1, dist1 = is_correct(mate1, args)
                                correct2, dist2 = is_correct(mate2, args)
                                if not correct1 and dist1 < 20 * args['wiggle']:
                                    self.dist_hist_incor[dist1] += 1
                                elif correct1:
                                    self.dist_hist_cor[dist1] += 1
                                if not correct2 and dist2 < 20 * args['wiggle']:
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
                            if not correct and dist < 20 * args['wiggle']:
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
                                self.dists.add_concordant_pair(mate1, mate2, correct1, correct2, self.ref)
                            elif al.discordant:
                                self.dists.add_discordant_pair(mate1, mate2, correct1, correct2, self.ref)
                            else:
                                self.dists.add_bad_end_read(al, correct, len(last_al.seq), self.ref)
                        elif not al.paired:
                            self.dists.add_unpaired_read(al, correct, self.ref)
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

    random.seed(args['seed'])
    np.random.seed(args['seed'])

    tim = Timing()
    tim.start_timer('Overall')

    # Create output directory if needed
    if not os.path.isdir(args['output_directory']):
        try:
            os.makedirs(args['output_directory'])
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s', datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG if args['verbose'] else logging.INFO)

    if args['write_logs'] or args['write_all']:
        fn = os.path.join(args['output_directory'], 'ts_logs.txt')
        fh = logging.FileHandler(fn)
        fh.setLevel(logging.DEBUG)
        logging.getLogger('').addHandler(fh)

    if args['U'] is not None and args['m1'] is not None:
        raise RuntimeError('Input must consist of only unpaired or only paired-end reads')
    
    # Construct command to invoke aligner - right now we only support Bowtie 2 and BWA-MEM
    aligner_class, alignment_class = Bowtie2, AlignmentBowtie2
    align_cmd = None
    if args['aligner'] == 'bowtie2':
        align_cmd = 'bowtie2 '
        if args['bt2_exe'] is not None:
            align_cmd = args['bt2_exe'] + " "
        align_cmd += ' '.join(aligner_args)
        align_cmd += ' --reorder --sam-no-qname-trunc --mapq-extra'
    elif args['aligner'] == 'bwa-mem':
        align_cmd = 'bwa mem '
        if args['bwa_exe'] is not None:
            align_cmd = args['bwa_exe'] + ' mem '
        align_cmd += ' '.join(aligner_args)
        aligner_class, alignment_class = BwaMem, AlignmentBwaMem
    elif args['aligner'] == 'mosaik':
        if args['use_concurrency']:
            raise RuntimeError('--use-concurrency cannott be combined with --aligner mosaik; '
                               'MOSAIK writes BAM directly to a file')
        align_cmd = 'MosaikAlign '
        if args['mosaik_align_exe'] is not None:
            align_cmd = args['mosaik_align_exe'] + ' '
        align_cmd += ' '.join(aligner_args)
        aligner_class, alignment_class = Mosaik, AlignmentMosaik
    elif args['aligner'] == 'snap':
        align_cmd = 'snap-aligner '
        if args['snap_exe'] is not None:
            align_cmd = args['snap_exe'] + ' '
        align_cmd += ' '.join(aligner_args)
        aligner_class, alignment_class = SnapAligner, AlignmentSnap
    elif args['aligner'] is not None:
        raise RuntimeError('Aligner not supported: "%s"' % args['aligner'])

    # for storing temp files and keep track of how big they get
    temp_man = TemporaryFileManager(args['temp_directory'])

    logging.info('Loading reference data')
    with ReferenceIndexed(args['ref']) as ref:
    
        # ##################################################
        # ALIGN REAL DATA (or read in alignments from SAM)
        # ##################################################

        unpaired_arg = args['U']
        paired_arg = None if args['m1'] is None else zip(args['m1'], args['m2'])

        input_sam_fh = None
        sam_arg = None

        # we're going to write alignment SAM somewhere, either to the output
        # directory (if requested) or to a temporary file
        sam_fn = os.path.join(args['output_directory'], 'input.sam')

        if not args['use_concurrency']:
            sam_arg = sam_fn
        else:
            input_sam_fh = open(sam_fn, 'w')

        if not args['use_concurrency']:
            tim.start_timer('Aligning input reads')

        logging.info('Command for aligning input data: "%s"' % align_cmd)
        aligner = aligner_class(align_cmd, args['index'],
                                unpaired=unpaired_arg, paired=paired_arg,
                                sam=sam_arg)

        if not args['use_concurrency']:
            logging.debug('Waiting for aligner to finish')
            while aligner.pipe.poll() is None:
                time.sleep(0.5)
            logging.debug('Aligner finished')
            tim.end_timer('Aligning input reads')

        test_data = DatasetOnDisk('test_data', temp_man)
        dists = Dists(args['max_allowed_fraglen'], fraction_even=args['fraction_even'], bias=args['low_score_bias'])
        cor_dist, incor_dist = defaultdict(int), defaultdict(int)
        
        result_test_q = Queue()
        if not args['use_concurrency']:
            tim.start_timer('Parsing input read alignments')
            with open(sam_fn, 'r') as sam_fh:
                # TODO: support parsing of BAM, so we can use MOSAIK
                othread = AlignmentReader(
                    args,
                    sam_fh,                # SAM/BAM file
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

        test_data.finalize()

        # Writing test dataset
        if args['write_test_data'] or args['write_all']:
            test_csv_fn_prefix = os.path.join(args['output_directory'], 'test')
            test_data.save(test_csv_fn_prefix)
            logging.info('Test data (CSV format) written to "%s*"' % test_csv_fn_prefix)

        test_data.purge()

        # Writing correct/incorrect distances
        if args['write_test_distances'] or args['write_all']:
            for short_name, long_name, hist in [('cor', 'Correct', cor_dist), ('incor', 'Incorrect', incor_dist)]:
                test_dist_fn = os.path.join(args['output_directory'], 'test_%s_dists.csv' % short_name)
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
        if iters == 2 and aligner.supports_mix():
            iters = 1
        if iters == 2:
            logging.info('Aligner does not accept unpaired/paired mix; training will have 2 rounds')

        training_data = DatasetOnDisk('training_data', temp_man)
        for paired, lab in [(True, 'paired-end'), (False, 'unpaired')]:
            both = False
            if paired and not dists.has_pairs():
                continue
            if not paired and not dists.has_unpaired_reads():
                continue
            if aligner.supports_mix() and dists.has_pairs() and dists.has_unpaired_reads():
                # Do both unpaired and paired simualted reads in one round
                both, lab = True, 'both paired-end and unpaired'

            def simulate(simw, unpaired_format, paired_format, aligner=None):
                """ Simulates reads.  Either pipes them directly to the aligner
                    (when not write_training_reads and args['use_concurrency']
                    or writes them to various files and returns those. """
                type_to_format = {'conc': paired_format,
                                  'disc': paired_format,
                                  'bad_end': paired_format,
                                  'unp': unpaired_format}
                write_training_reads = args['write_training_reads'] or args['write_all']
                training_out_fn, training_out_fh = {}, {}
                if write_training_reads or not args['use_concurrency']:
                    types = []
                    if paired or both:
                        types.extend(zip(['conc', 'disc', 'bad_end'], [paired_format] * 3))
                    if not paired or both:
                        types.append(('unp', unpaired_format))
                    for t, frmt in types:
                        fn_base = 'training_%s.%s' % (t, frmt)
                        fn = os.path.join(args['output_directory'], fn_base)
                        if not write_training_reads:
                            fn = temp_man.get_filename(fn_base, 'tandem reads')
                        training_out_fn[t] = fn
                        training_out_fh[t] = open(fn, 'w')

                logging.info('  Simulating reads')
                n_simread, typ_count = 0, defaultdict(int)
                for t, rdp1, rdp2 in simw.simulate_batch(args['sim_fraction'], args['sim_unp_min'],
                                                         args['sim_conc_min'], args['sim_disc_min'],
                                                         args['sim_bad_end_min']):
                    frmt = type_to_format[t]
                    if t in training_out_fh:
                        # read is going to a file
                        if frmt == 'tab6':
                            # preferred format for Bowtie 2
                            training_out_fh[t].write(Read.to_tab6(rdp1, rdp2) + '\n')
                        elif frmt == 'interleaved_fastq':
                            # preferred paired-end format for BWA & SNAP
                            training_out_fh[t].write(Read.to_interleaved_fastq(rdp1, rdp2) + '\n')
                        elif frmt == 'fastq':
                            # preferred unpaired format for BWA & SNAP
                            assert rdp2 is None
                            training_out_fh[t].write(Read.to_fastq(rdp1) + '\n')
                        else:
                            raise RuntimeError('Bad training read output format "%s"' % frmt)
                    if aligner is not None:
                        # here, there's no guarantee about the order in which
                        # reads are being fed to the aligner, so the aligner
                        # had better be ready to accept a mixed stream of
                        # unpaired and paired
                        assert aligner.supports_mix()
                        aligner.put(rdp1, rdp2)  # read is going directly to the aligner
                    typ_count[t] += 1
                    n_simread += 1
                    if (n_simread % 20000) == 0:
                        logging.info('    simulated %d reads (%d conc, %d disc, %d bad-end, %d unp)' %
                                     (n_simread, typ_count['conc'], typ_count['disc'],
                                      typ_count['bad_end'], typ_count['unp']))

                for t in training_out_fh.keys():
                    training_out_fh[t].close()
                    logging.info('  Training reads written to "%s"' % training_out_fn[t])

                logging.info('  Finished simulating reads (%d conc, %d disc, %d bad_end, %d unp)' %
                             (typ_count['conc'], typ_count['disc'], typ_count['bad_end'], typ_count['unp']))

                return typ_count, training_out_fn

            # TODO: optional random-access simulator
            simw = FragmentSimSerial2(args['ref'], dists)

            if not args['use_concurrency']:

                #
                # Simulate
                #

                tim.start_timer('Simulating tandem reads')
                logging.info('Simulating tandem reads (%s)' % lab)
                typ_sim_count, training_out_fn = simulate(simw,
                                                          aligner.preferred_unpaired_format(),
                                                          aligner.preferred_paired_format())

                unpaired_arg = None
                if 'unp' in training_out_fn or 'bad_end' in training_out_fn:
                    unpaired_arg = []
                    for t in ['unp', 'bad_end']:
                        if t in training_out_fn:
                            unpaired_arg.append(training_out_fn[t])
                paired_combined_arg = None
                if 'conc' in training_out_fn or 'disc' in training_out_fn:
                    paired_combined_arg = []
                    for t in ['conc', 'disc']:
                        if t in training_out_fn:
                            paired_combined_arg.append(training_out_fn[t])
                    if len(paired_combined_arg) > 1:
                        # new file
                        fn_base = 'training_concdisc.%s' % aligner.preferred_paired_format()
                        fn = temp_man.get_filename(fn_base, 'tandem reads')
                        with open(fn, 'w') as fh:
                            for ifn in paired_combined_arg:
                                with open(ifn) as ifh:
                                    for ln in ifh:
                                        fh.write(ln)
                        paired_combined_arg = [fn]

                assert unpaired_arg is not None or paired_combined_arg is not None

                logging.info('Finished simulating tandem reads')
                tim.end_timer('Simulating tandem reads')

                #
                # Align
                #

                logging.info('Aligning tandem reads (%s)' % lab)
                tim.start_timer('Aligning tandem reads')

                def _wait_for_aligner(_al):
                    while _al.pipe.poll() is None:
                        time.sleep(0.5)

                sam_fn = temp_man.get_filename('training.sam', 'tandem sam')
                if aligner.supports_mix():
                    aligner = aligner_class(align_cmd, args['index'],  # no concurrency
                                            unpaired=unpaired_arg, paired_combined=paired_combined_arg,
                                            sam=sam_fn, input_format=aligner.preferred_paired_format())
                    # the aligner_class gets to decide what order to do unpaired/paired
                    _wait_for_aligner(aligner)
                    logging.debug('Finished aligning unpaired and paired-end tandem reads')
                else:
                    paired_sam, unpaired_sam = None, None
                    if unpaired_arg is not None:
                        unpaired_sam = temp_man.get_filename('training_unpaired.sam', 'tandem sam')
                        aligner = aligner_class(align_cmd, args['index'],  # no concurrency
                                                unpaired=unpaired_arg, paired_combined=None,
                                                sam=unpaired_sam, input_format=aligner.preferred_unpaired_format())
                        _wait_for_aligner(aligner)
                        logging.debug('Finished aligning unpaired tandem reads')

                    if paired_combined_arg is not None:
                        paired_sam = temp_man.get_filename('training_paired.sam', 'tandem sam')
                        aligner = aligner_class(align_cmd, args['index'],  # no concurrency
                                                unpaired=None, paired_combined=paired_combined_arg,
                                                sam=paired_sam, input_format=aligner.preferred_paired_format())
                        _wait_for_aligner(aligner)
                        logging.debug('Finished aligning paired-end tandem reads')

                    logging.debug('Concatenating unpaired and paired-end files')
                    with open(sam_fn, 'w') as ofh:
                        for fn in [paired_sam, unpaired_sam]:
                            if fn is not None:
                                with open(fn) as fh:
                                    for ln in fh:
                                        ofh.write(ln)

                # remove temporary reads
                temp_man.update_peak()
                temp_man.remove_group('tandem reads')

                cor_dist, incor_dist = defaultdict(int), defaultdict(int)

                logging.info('Parsing tandem alignments (%s)' % lab)
                tim.end_timer('Aligning tandem reads')
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

                logging.info('Opening aligner process using concurrency')
                # Does the aligner always get unpaireds before paireds or vice versa?
                aligner = aligner_class(align_cmd, args['index'], pairs_only=paired)
                sam_fn = os.path.join(args['output_directory'], 'training.sam')
                training_sam_fh = open(sam_fn, 'w') if (args['write_training_sam'] or args['write_all']) else None
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

                # uses aligner.put to pump simulated reads directly into aligner
                # always simulates
                typ_sim_count, training_out_fn = simulate(simw, 'tab6', 'tab6', aligner)
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

            if not dists.sc_dist_unp.empty() and dists.sc_dist_unp.has_correctness_info:
                logging.info('    %d unpaired draws, %0.3f%% correct' % (dists.sc_dist_unp.num_drawn,
                                                                         100*dists.sc_dist_unp.frac_correct()))
            if not dists.sc_dist_conc.empty() and dists.sc_dist_conc.has_correctness_info:
                logging.info('    %d concordant draws, %0.3f%%/%0.3f%% correct' % (dists.sc_dist_conc.num_drawn,
                                                                                   100*dists.sc_dist_conc.frac_correct1(),
                                                                                   100*dists.sc_dist_conc.frac_correct2()))
            if not dists.sc_dist_disc.empty() and dists.sc_dist_disc.has_correctness_info:
                logging.info('    %d discordant draws, %0.3f%%/%0.3f%% correct' % (dists.sc_dist_disc.num_drawn,
                                                                                   100*dists.sc_dist_disc.frac_correct1(),
                                                                                   100*dists.sc_dist_disc.frac_correct2()))
            if not dists.sc_dist_bad_end.empty() and dists.sc_dist_bad_end.has_correctness_info:
                logging.info('    %d bad-end draws, %0.3f%% correct' % (dists.sc_dist_bad_end.num_drawn,
                                                                        100*dists.sc_dist_bad_end.frac_correct()))

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

    training_data.finalize()

    tim.start_timer('Writing training data')

    # Writing training data
    training_csv_fn_prefix = os.path.join(args['output_directory'], 'training')
    training_data.save(training_csv_fn_prefix)
    logging.info('Training data (CSV format) written to "%s*"' % training_csv_fn_prefix)
    training_data.purge()

    tim.end_timer('Writing training data')

    temp_man.purge()
    logging.info('Peak temporary file size %0.2fMB' % (temp_man.peak_size / (1024.0 * 1024)))

    tim.end_timer('Overall')
    for ln in str(tim).split('\n'):
        if len(ln) > 0:
            logging.info(ln)
    if args['write_timings'] or args['write_all']:
        with open(os.path.join(args['output_directory'], 'timing.tsv'), 'w') as fh:
            fh.write(str(tim))


def add_args(parser):
    # Inputs
    parser.add_argument('--ref', metavar='path', type=str, nargs='+', required=True,
                        help='FASTA file(s) containing reference genome sequences')
    parser.add_argument('--pickle-ref', metavar='path', type=str,
                        help='Pickle FASTA input for speed, or use pickled version if it exists already.  Pickled '
                             'version is stored at given path')
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

    parser.add_argument('--max-allowed-fraglen', metavar='int', type=int, default=100000, required=False,
                        help='When simulating fragments, observed fragments longer than this will be'
                             'truncated to this length')
    parser.add_argument('--low-score-bias', metavar='float', type=float, default=1.0, required=False,
                        help='When simulating reads, we randomly select a real read\'s alignment profile'
                             'as a template.  A higher value for this parameter makes it more likely'
                             'we\'ll choose a low-scoring alignment.  If set to 1, all templates are'
                             'equally likely.')
    parser.add_argument('--fraction-even', metavar='float', type=float, default=1.0, required=False,
                        help='Fraction of the time to sample templates from the unstratified input '
                             'sample versus the stratified sample.')

    parser.add_argument('--wiggle', metavar='int', type=int, default=30, required=False,
                        help='Wiggle room to allow in starting position when determining whether alignment is correct')

    parser.add_argument('--sam-input', metavar='path', type=str,
                        help='Input SAM file to apply training data to.  Use with --training-input.')
    parser.add_argument('--bt2-exe', metavar='path', type=str, help='Path to Bowtie 2 exe')
    parser.add_argument('--bwa-exe', metavar='path', type=str, help='Path to BWA exe')
    parser.add_argument('--mosaik-align-exe', metavar='path', type=str, help='Path to MosaikAlign exe')
    parser.add_argument('--mosaik-build-exe', metavar='path', type=str, help='Path to MosaikBuild exe')
    parser.add_argument('--snap-exe', metavar='path', type=str, help='Path to snap-aligner exe')
    parser.add_argument('--aligner', metavar='name', default='bowtie2', type=str,
                        help='bowtie2 | bwa-mem | mosaik | snap')

    # For when input is itself simulated, so we can output a Dataset with the
    # 'correct' column filled in properly
    parser.add_argument('--correct-chromosomes', metavar='list', type=str, nargs='+',
                        help='Label test data originating from any of these chromosomes as "correct."  Useful for '
                             'tests on real-world data where it is known that the data came from a parituclar '
                             'chromosome.')

    # Output file-related arguments
    parser.add_argument('--temp-directory', metavar='path', type=str, required=False,
                        help='Write temporary files to this directory; default: uses environment variables '
                             'like TMPDIR, TEMP, etc')
    parser.add_argument('--output-directory', metavar='path', type=str, required=True,
                        help='Write outputs to this directory')
    parser.add_argument('--write-training-reads', action='store_const', const=True, default=False,
                        help='Write FASTQ for the training reads to "training.fastq" in output directory')
    parser.add_argument('--write-test-data', action='store_const', const=True, default=False,
                        help='Write Dataset object for training data.  "Correct" column set to all None\'s.')
    parser.add_argument('--write-training-sam', action='store_const', const=True, default=False,
                        help='Write SAM alignments for the training reads to "training.sam" in output directory')
    parser.add_argument('--write-test-distances', action='store_const', const=True, default=False,
                        help='Write distances between true/actual alignments.')
    parser.add_argument('--write-training-distances', action='store_const', const=True, default=False,
                        help='Write distances between true/actual alignments.')
    parser.add_argument('--write-timings', action='store_const', const=True, default=False,
                        help='Write timing info to "timing.tsv".')
    parser.add_argument('--write-logs', action='store_const', const=True, default=False,
                        help='Write logs to "ts_log.txt" in the output directory.')
    parser.add_argument('--write-all', action='store_const', const=True, default=False,
                        help='Same as specifying all --write-* options')
    parser.add_argument('--compress-output', action='store_const', const=True, default=False,
                        help='gzip all output files')

    # Resource usage
    parser.add_argument('--use-concurrency', action='store_const', const=True, default=False,
                        help='Use pipes instead of temporary files.  Reduces disk usage '
                             'but requires more memory.  Doesn\'t work with all aligners.')


def go_profile(args, aligner_args):
    if args['profile']:
        import cProfile
        cProfile.run('go(args, aligner_args)')
    else:
        go(args, aligner_args)


def parse_aligner_parameters_from_argv(_argv):
    argv = _argv[:]
    aligner_args = []
    in_args = False
    new_argv = None
    for i, ag in enumerate(argv):
        if in_args:
            aligner_args.append(ag)
        elif ag == '--':
            new_argv = argv[:i]
            in_args = True
    if new_argv is None:
        new_argv = argv
    return new_argv, aligner_args


if __name__ == "__main__":
    
    import argparse

    _parser = argparse.ArgumentParser(
        description='Align a collection of input reads, simulate a tandem'
                    'dataset, align the tandem dataset, and emit both the'
                    'input read alignments and the training data derived from'
                    'the tandem read alignments.')

    if '--version' in sys.argv:
        print 'Tandem simulator, version ' + VERSION
        sys.exit(0)

    add_args(_parser)

    # Some basic flags
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    _parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
    _parser.add_argument('--version', action='store_const', const=True, default=False, help='Print version and quit')

    _argv, _aligner_args = parse_aligner_parameters_from_argv(sys.argv)
    _args = _parser.parse_args(_argv[1:])

    go_profile(vars(_args), _aligner_args)
