#!/usr/bin/env python

from __future__ import print_function

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
import random
import time
import logging
import errno
import csv
from collections import defaultdict
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full  # python 3.x

# Modules that are part of the tandem simulator
from simplesim import FragmentSimSerial2
from read import Read
from bowtie2 import AlignmentBowtie2, Bowtie2
from bwamem import AlignmentBwaMem, BwaMem
from snap import AlignmentSnap, SnapAligner
from reference import ReferenceIndexed, ReferenceOOB
from tempman import TemporaryFileManager
from score_dists import CollapsedScoreDist, ScoreDist, CollapsedScorePairDist, ScorePairDist

bin_dir = os.path.dirname(os.path.realpath(__file__))

VERSION = '0.2.0'


class Dists(object):
    
    """ Encapsulates the distributions that we capture from real data.
        We collect random subsets of qualities/edits or unpaired reads
        and same for paired-end mate 1s and mate 2s.  We also collect
        data on concordant and discordantly aligned pairs, such as
        their fragment length and strands. """

    def __init__(self,
                 max_allowed_fraglen=100000,
                 fraction_even=0.5,
                 bias=1.0,
                 use_ref_for_edit_distance=False,
                 reference=None,
                 reservoir_size=10000):
        self.use_ref_for_edit_distance = use_ref_for_edit_distance
        self.reference = reference
        if fraction_even >= 1.0:
            self.sc_dist_unp = ScoreDist(reference=self.reference, big_k=reservoir_size)
            self.sc_dist_bad_end = ScoreDist(reference=self.reference, big_k=reservoir_size)
            self.sc_dist_conc = ScorePairDist(reference=self.reference, big_k=reservoir_size,
                                              max_allowed_fraglen=max_allowed_fraglen)
            self.sc_dist_disc = ScorePairDist(reference=self.reference, big_k=reservoir_size,
                                              max_allowed_fraglen=max_allowed_fraglen)
        else:
            self.sc_dist_unp = CollapsedScoreDist(reference=self.reference, big_k=reservoir_size,
                                                  fraction_even=fraction_even, bias=bias)
            self.sc_dist_bad_end = CollapsedScoreDist(reference=self.reference, big_k=reservoir_size,
                                                      fraction_even=fraction_even, bias=bias)
            self.sc_dist_conc = CollapsedScorePairDist(reference=self.reference,
                                                       big_k=reservoir_size, max_allowed_fraglen=max_allowed_fraglen,
                                                       fraction_even=fraction_even, bias=bias)
            self.sc_dist_disc = CollapsedScorePairDist(reference=self.reference,
                                                       big_k=reservoir_size, max_allowed_fraglen=max_allowed_fraglen,
                                                       fraction_even=fraction_even, bias=bias)

    def finalize(self):
        self.sc_dist_unp.finalize()
        self.sc_dist_bad_end.finalize()
        self.sc_dist_conc.finalize()
        self.sc_dist_disc.finalize()

    def add_concordant_pair(self, al1, al2, correct1, correct2):
        """ Add concordant paired-end read alignment to the model """
        self.sc_dist_conc.add(al1, al2, correct1, correct2, use_ref_for_edit_distance=self.use_ref_for_edit_distance)

    def add_discordant_pair(self, al1, al2, correct1, correct2):
        """ Add discordant paired-end read alignment to the model """
        self.sc_dist_disc.add(al1, al2, correct1, correct2, use_ref_for_edit_distance=self.use_ref_for_edit_distance)

    def add_unpaired_read(self, al, correct):
        """ Add unpaired read alignment to the model """
        self.sc_dist_unp.add(al, correct, use_ref_for_edit_distance=self.use_ref_for_edit_distance)

    def add_bad_end_read(self, al, correct, ordlen):
        """ Add bad-end read alignment to the model """
        self.sc_dist_bad_end.add(al, correct, ordlen=ordlen, use_ref_for_edit_distance=self.use_ref_for_edit_distance)

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

    def from_tables(self):
        pass


"""
Pass 1: parsing SAM resulting from alignment of input reads.
======

The fast C++ code does the parsing and emits tables.

One question is: can the output of the C++ code allow us to completely bypass
the process of writing via the DatasetOnDisk object?  I believe so.  I think
we can completely get rid of samples.py and just modify the code in predict.py
accordingly.  Note that predict.py does not use DatasetOnDisk; it uses pandas
read_csv instead.

Can we get rid of the AlignmentReader object?  Yes.  There's no longer any
reason for the Python code to go through the SAM line by line.  Instead, the
fast C++ code emits two tables which we parse separately from Python.  One
is the input model, which ts.py parses (maybe using pandas).  And the other
is the record table, which predict.py parses (maybe using pandas).

Another question is: can we add parsing so that, in the event that the read
is simulated, we can assess whether the alignment is correct.  We can do this,
it just involves porting some Python code over to C++.

Pass 2: parsing SAM resulting from alignment of tandem reads.
======

The fast C++ code does the parsing and emits tables.  One additional
complication is that we want it to also assess correctness of the alignments,
which involves parsing some crud that we put in the read names.

Again, seems like we can avoid DatasetOnDisk entirely.

"""


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


def sanity_check_binary(exe):
    if not os.path.exists(exe):
        raise RuntimeError('Binary "%s" does not exist' % exe)
    else:
        if not os.access(exe, os.X_OK):
            raise RuntimeError('Binary "%s" exists but is not executable' % exe)


def go(args, aligner_args, aligner_unpaired_args, aligner_paired_args):
    """ Main driver for tandem simulator """

    random.seed(args['seed'])
    import numpy
    numpy.random.seed(args['seed'])

    tim = Timing()
    tim.start_timer('Overall')

    # Create output directory if needed
    if not os.path.isdir(args['output_directory']):
        try:
            os.makedirs(args['output_directory'])
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    if not args['input_reads_simulated'] and args['write_test_distances']:
        raise RuntimeError('if --write-test-distances is set, --input-reads-simulated must also be')

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
        aligner_args.extend(['--mapq-extra'])  # TODO: do we want --reorder?
    elif args['aligner'] == 'bwa-mem':
        align_cmd = 'bwa mem '
        if args['bwa_exe'] is not None:
            align_cmd = args['bwa_exe'] + ' mem '
        aligner_class, alignment_class = BwaMem, AlignmentBwaMem
    elif args['aligner'] == 'snap':
        align_cmd = 'snap-aligner '
        if args['snap_exe'] is not None:
            align_cmd = args['snap_exe'] + ' '
        aligner_class, alignment_class = SnapAligner, AlignmentSnap
    elif args['aligner'] is not None:
        raise RuntimeError('Aligner not supported: "%s"' % args['aligner'])

    # for storing temp files and keep track of how big they get
    temp_man = TemporaryFileManager(args['temp_directory'])

    def _get_pass1_file_prefix():
        if args['write_all'] or args['write_input_intermediates']:
            def _purge():
                pass
            return os.path.join(args['output_directory'], 'input_intermediates'), _purge
        else:
            dr = temp_man.get_dir('input_intermediates', 'input_intermediates')

            def _purge():
                dr.remove_group('input_intermediates')
            return os.path.join(dr, 'tmp'), _purge

    # Note: this is only needed when either MD:Z information is missing from
    # alignments or when --ref-soft-clipping is specified
    logging.info('Loading reference data')
    with ReferenceIndexed(args['ref']) as ref:
    
        # ##################################################
        # ALIGN REAL DATA (or read in alignments from SAM)
        # ##################################################

        tim.start_timer('Aligning input reads')
        sam_fn = os.path.join(args['output_directory'], 'input.sam')
        logging.info('Command for aligning input data: "%s"' % align_cmd)
        aligner = aligner_class(
            align_cmd,
            aligner_args,
            aligner_unpaired_args,
            aligner_paired_args,
            args['index'],
            unpaired=args['U'],
            paired=None if args['m1'] is None else zip(args['m1'], args['m2']),
            sam=sam_fn)

        logging.debug('  waiting for aligner to finish...')
        while aligner.pipe.poll() is None:
            time.sleep(0.5)
        logging.debug('  aligner finished; results in "%s"' % sam_fn)
        tim.end_timer('Aligning input reads')

        tim.start_timer('Parsing and generating models/records from alignments')

        tim.start_timer('Parsing input reads')
        parse_input_exe = "%s/qsim-parse-input" % bin_dir
        sanity_check_binary(parse_input_exe)
        pass1_prefix, pass1_cleanup = _get_pass1_file_prefix()
        os.system("%s %s %s" % (parse_input_exe, sam_fn, pass1_prefix))

        logging.debug('  parsing finished; results in "%s.*"' % pass1_prefix)

        # TODO: Set up table files
        dists = Dists(
            args['max_allowed_fraglen'],
            fraction_even=args['fraction_even'],
            bias=args['low_score_bias'],
            use_ref_for_edit_distance=args['ref_soft_clipping'],
            reservoir_size=args['input_model_size'])
        dists.parse_from_tables('some way of referring to all the tables')

        tim.start_timer('Parsing and generating models/records from alignments')

        # ##################################################
        # SIMULATE TANDEM READS
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

        for paired, lab in [(True, 'paired-end'), (False, 'unpaired')]:
            both = False
            if paired and not dists.has_pairs():
                logging.debug('No paired-end reads in the input model')
                continue
            if not paired and not dists.has_unpaired_reads():
                logging.debug('No unpaired reads in the input model')
                continue
            if aligner.supports_mix() and dists.has_pairs() and dists.has_unpaired_reads():
                # Do both unpaired and paired simualted reads in one round
                both, lab = True, 'both paired-end and unpaired'

            def simulate(simw, unpaired_format, paired_format, aligner=None):
                type_to_format = {'conc': paired_format,
                                  'disc': paired_format,
                                  'bad_end': paired_format,
                                  'unp': unpaired_format}
                write_training_reads = args['write_training_reads'] or args['write_all']
                training_out_fn, training_out_fh = {}, {}
                types = []
                if paired or both:
                    types.extend(zip(['conc', 'disc', 'bad_end'], [paired_format] * 3))
                if not paired or both:
                    types.extend(zip(['unp'], [unpaired_format]))
                for t, frmt in types:
                    fn_base = 'training_%s.%s' % (t, frmt)
                    fn = os.path.join(args['output_directory'], fn_base)
                    if not write_training_reads:
                        fn = temp_man.get_file(fn_base, 'tandem reads')
                    training_out_fn[t] = fn
                    training_out_fh[t] = open(fn, 'w')

                logging.info('  Simulating reads')

                _finish_sim_profiling = lambda: None
                if args['profile_simulating']:
                    import cProfile
                    import pstats
                    import StringIO
                    pr = cProfile.Profile()
                    pr.enable()

                    def _finish_sim_profiling():
                        pr.disable()
                        s = StringIO.StringIO()
                        ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
                        ps.print_stats(30)
                        print(s.getvalue())

                n_simread, typ_count = 0, defaultdict(int)
                for t, rdp1, rdp2 in simw.simulate_batch(args['sim_fraction'], args['sim_unp_min'],
                                                         args['sim_conc_min'], args['sim_disc_min'],
                                                         args['sim_bad_end_min']):
                    assert t in type_to_format
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

                _finish_sim_profiling()

                logging.info('  Finished simulating reads (%d conc, %d disc, %d bad_end, %d unp)' %
                             (typ_count['conc'], typ_count['disc'], typ_count['bad_end'], typ_count['unp']))

                return typ_count, training_out_fn

            # TODO: optional random-access simulator
            simw = FragmentSimSerial2(args['ref'], dists)

            #
            # Simulate
            #

            tim.start_timer('Simulating tandem reads')
            logging.info('Simulating tandem reads (%s)' % lab)
            typ_sim_count, training_out_fn = simulate(simw,
                                                      aligner.preferred_unpaired_format(),
                                                      aligner.preferred_paired_format())

            # TODO: why is this correct?  the simulated bad_end reads seem
            # to be paired-end, and should be aligned in paired-end mode
            unpaired_arg = None
            if False:
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
                        fn_base = 'training_paired.%s' % aligner.preferred_paired_format()
                        fn = temp_man.get_file(fn_base, 'tandem reads')
                        with open(fn, 'w') as fh:
                            for ifn in paired_combined_arg:
                                with open(ifn) as ifh:
                                    for ln in ifh:
                                        fh.write(ln)
                        paired_combined_arg = [fn]
            else:
                if 'unp' in training_out_fn:
                    unpaired_arg = []
                    for t in ['unp']:
                        if t in training_out_fn:
                            unpaired_arg.append(training_out_fn[t])
                paired_combined_arg = None
                if 'conc' in training_out_fn or 'disc' in training_out_fn or 'bad_end' in training_out_fn:
                    paired_combined_arg = []
                    for t in ['conc', 'disc', 'bad_end']:
                        if t in training_out_fn:
                            paired_combined_arg.append(training_out_fn[t])
                    if len(paired_combined_arg) > 1:
                        # new file
                        fn_base = 'training_paired.%s' % aligner.preferred_paired_format()
                        fn = temp_man.get_file(fn_base, 'tandem reads')
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

            tim.start_timer('Aligning tandem reads')

            def _wait_for_aligner(_al):
                while _al.pipe.poll() is None:
                    time.sleep(0.5)

            if args['write_training_sam'] or args['write_all']:
                sam_fn = os.path.join(args['output_directory'], 'training.sam')
            else:
                sam_fn = temp_man.get_file('training.sam', 'tandem sam')

            if aligner.supports_mix():
                logging.info('Aligning tandem reads (%s, mix)' % lab)
                aligner = aligner_class(align_cmd,
                                        aligner_args, aligner_unpaired_args, aligner_paired_args,
                                        args['index'],
                                        unpaired=unpaired_arg, paired_combined=paired_combined_arg,
                                        sam=sam_fn, input_format=aligner.preferred_paired_format())
                # the aligner_class gets to decide what order to do unpaired/paired
                _wait_for_aligner(aligner)
                logging.debug('Finished aligning unpaired and paired-end tandem reads')
            else:
                paired_sam, unpaired_sam = None, None
                if unpaired_arg is not None:
                    logging.info('Aligning tandem reads (%s, unpaired)' % lab)
                    unpaired_sam = temp_man.get_file('training_unpaired.sam', 'tandem sam')
                    aligner = aligner_class(align_cmd,
                                            aligner_args, aligner_unpaired_args, aligner_paired_args,
                                            args['index'],
                                            unpaired=unpaired_arg, paired_combined=None,
                                            sam=unpaired_sam, input_format=aligner.preferred_unpaired_format())
                    _wait_for_aligner(aligner)
                    logging.debug('Finished aligning unpaired tandem reads')

                if paired_combined_arg is not None:
                    logging.info('Aligning tandem reads (%s, paired)' % lab)
                    paired_sam = temp_man.get_file('training_paired.sam', 'tandem sam')
                    aligner = aligner_class(align_cmd,
                                            aligner_args, aligner_unpaired_args, aligner_paired_args,
                                            args['index'],
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
                sam_reader = csv.reader(sam_fh, delimiter='\t', quotechar=None)
                reader = AlignmentReader(
                    args,
                    sam_reader,
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
    parser.add_argument('--input-model-size', metavar='int', type=int, default=10000, required=False,
                        help='Number of templates to keep when building input model.')
    parser.add_argument('--ref-soft-clipping', action='store_const', const=True, default=False,
                        help='Use bases from reference (instead of random bases) to re-align soft clipped bases.')
    parser.add_argument('--input-reads-simulated', action='store_const', const=True, default=False,
                        help='Input reads are simulated and should be evaluated for True/False.')

    parser.add_argument('--wiggle', metavar='int', type=int, default=30, required=False,
                        help='Wiggle room to allow in starting position when determining whether alignment is correct')

    parser.add_argument('--sam-input', metavar='path', type=str,
                        help='Input SAM file to apply training data to.  Use with --training-input.')
    parser.add_argument('--bt2-exe', metavar='path', type=str, help='Path to Bowtie 2 exe')
    parser.add_argument('--bwa-exe', metavar='path', type=str, help='Path to BWA exe')
    parser.add_argument('--snap-exe', metavar='path', type=str, help='Path to snap-aligner exe')
    parser.add_argument('--aligner', metavar='name', default='bowtie2', type=str,
                        help='bowtie2 | bwa-mem | snap')

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
    parser.add_argument('--write-timings', action='store_const', const=True, default=False,
                        help='Write timing info to "timing.tsv".')
    parser.add_argument('--write-logs', action='store_const', const=True, default=False,
                        help='Write logs to "ts_log.txt" in the output directory.')
    parser.add_argument('--write-all', action='store_const', const=True, default=False,
                        help='Same as specifying all --write-* options')
    parser.add_argument('--compress-output', action='store_const', const=True, default=False,
                        help='gzip all output files')
    parser.add_argument('--profile-parsing', action='store_const', const=True, default=False,
                        help='output profiling info related to parsing')
    parser.add_argument('--profile-simulating', action='store_const', const=True, default=False,
                        help='output profiling info related to simulating reads')


def go_profile(args, aligner_args, aligner_unpaired_args, aligner_paired_args):
    if args['profile']:
        import cProfile
        cProfile.run('go(args, aligner_args, aligner_unpaired_args, aligner_paired_args)')
    else:
        go(args, aligner_args, aligner_unpaired_args, aligner_paired_args)


def parse_aligner_parameters_from_argv(_argv):
    argv = _argv[:]
    sections = [[]]
    for arg in argv:
        if arg == '--':
            sections.append([])
        else:
            sections[-1].append(arg)
    new_argv = sections[0]
    aligner_args = [] if len(sections) < 2 else sections[1]
    aligner_unpaired_args = [] if len(sections) < 3 else sections[2]
    aligner_paired_args = [] if len(sections) < 4 else sections[3]
    return new_argv, aligner_args, aligner_unpaired_args, aligner_paired_args


if __name__ == "__main__":
    
    import argparse

    _parser = argparse.ArgumentParser(
        description='Align a collection of input reads, simulate a tandem'
                    'dataset, align the tandem dataset, and emit both the'
                    'input read alignments and the training data derived from'
                    'the tandem read alignments.')

    if '--version' in sys.argv:
        print('Tandem simulator, version ' + VERSION)
        sys.exit(0)

    add_args(_parser)

    # Some basic flags
    _parser.add_argument('--profile', action='store_const', const=True, default=False, help='Print profiling info')
    _parser.add_argument('--verbose', action='store_const', const=True, default=False, help='Be talkative')
    _parser.add_argument('--version', action='store_const', const=True, default=False, help='Print version and quit')

    _argv, _aligner_args, _aligner_unpaired_args, _aligner_paired_args = parse_aligner_parameters_from_argv(sys.argv)
    _args = _parser.parse_args(_argv[1:])

    go_profile(vars(_args), _aligner_args, _aligner_unpaired_args, _aligner_paired_args)
