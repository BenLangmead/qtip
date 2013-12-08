#!/usr/bin/env python

'''
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
- Handle local alignments

To load training and test data, copy samples.py to dir, then:
import cPickle
import gzip
from samples import Dataset
testData, trainingData = Dataset(), Dataset()
testData.load('test.pickle.gz')
trainingData.load('training.pickle.gz')
testU, testM, testC = testData.toDataFrames()
trainingU, trainingM, trainingC = trainingData.toDataFrames()
'''

__author__ = "Ben Langmead"
__email__ = "langmea@cs.jhu.edu"

# depends on bowtie2.py, bwamem.py, randutil.py, sam.py, samples.py, simualtors.py

import os
import sys
import math
import string
import random
import bisect
import subprocess
import tempfile
import signal
import traceback
import time
import logging
import gzip
from collections import defaultdict
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full  # python 3.x
from threading import Thread

# Modules that are part of the tandem simulator
import simulators
from sam import cigarToList, mdzToList, cigarMdzToStacked
from samples import UnpairedTuple, PairedTuple, Dataset
from randutil import ReservoirSampler, WeightedRandomGenerator
from simplesim import SequenceSimulator, mutate
from read import Read, Alignment
from bowtie2 import AlignmentBowtie2, Bowtie2
from reference import ReferenceIndexed, ReferencePicklable, ReferenceOOB

def quit_handler(signum,frame): traceback.print_stack()
signal.signal(signal.SIGQUIT,quit_handler)

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

class ScoreDist(object):
    ''' Capture a list of tuples, where each tuples represents the following
        traits observed in an alignment: (a) orientation, (b) quality string,
        (c) read side of stacked alignment, (d) ref side of stacked alignment,
        (e) score.  Keeping these tuples allows us to generate new reads that
        mimic observed reads in these key ways. '''
    
    def __init__(self, k=10000):
        ''' Make a reservoir sampler for holding the tuples. '''
        self.res = ReservoirSampler(k)
    
    def draw(self):
        ''' Draw from the reservoir '''
        assert len(self.res) > 0
        return self.res.draw()
    
    def add(self, al, ref):
        ''' Convert given alignment to a tuple and add it to the reservoir
            sampler. '''
        assert al.cigar is not None and al.mdz is not None
        cigarList = cigarToList(al.cigar)
        mdzList = mdzToList(al.mdz)
        sc = al.bestScore
        # Get stacked alignment
        rdAln, rfAln = al.stackedAlignment(alignSoftClipped=True, ref=ref)
        self.res.add((al.fw, al.qual, rdAln, rfAln, sc))
    
    def __len__(self):
        ''' Return number of tuples that have been added '''
        return len(self.res)
    
    def empty(self):
        ''' Return true iff no tuples have been added '''
        return len(self) == 0

class ScorePairDist(object):
    ''' Capture a list of tuple pairs, where each tuple pair represents the
        following in an alignment of a paired-end reads where both ends
        aligned: (a) orientation, (b) quality string, (c) read side of stacked
        alignment, (d) ref side of stacked alignment, (e) score.  We record
        pairs of tuples where the first element of the pair corresponds to
        mate 1 and the second to mate 2. '''
    
    def __init__(self, k=10000):
        ''' Make a reservoir sampler for holding the tuples. '''
        self.res = ReservoirSampler(k)
    
    def draw(self):
        ''' Draw from the reservoir '''
        assert len(self.res) > 0
        return self.res.draw()
    
    def add(self, al1, al2, ref):
        ''' Convert given alignment pair to a tuple and add it to the
            reservoir sampler. '''
        # Extract quality string
        assert al1.cigar is not None and al1.mdz is not None
        assert al2.cigar is not None and al2.mdz is not None
        # Extract CIGAR, MD:Z
        cigarList1, cigarList2 = cigarToList(al1.cigar), cigarToList(al2.cigar)
        mdzList1, mdzList2 = mdzToList(al1.mdz), mdzToList(al2.mdz)
        sc1, sc2 = al1.bestScore, al2.bestScore
        # Get stacked alignment
        rdAln1, rfAln1 = al1.stackedAlignment(alignSoftClipped=True, ref=ref)
        rdAln2, rfAln2 = al2.stackedAlignment(alignSoftClipped=True, ref=ref)
        self.res.add(((al1.fw, al1.qual, rdAln1, rfAln1, sc1),
                      (al2.fw, al2.qual, rdAln2, rfAln2, sc2)))
    
    def __len__(self):
        ''' Return number of tuples that have been added '''
        return len(self.res)
    
    def empty(self):
        ''' Return true iff no tuples have been added '''
        return len(self) == 0

class Dist(object):
    ''' Basically a histogram. '''
    
    def __init__(self):
        self.hist = defaultdict(int)
        self.tot = 0
        self.changed = True
        self.gen = None
    
    def __str__(self):
        ''' Return string representation of histogram '''
        return str(self.hist)
    
    def __len__(self):
        ''' Return the number of items added to the histogram '''
        return self.tot
    
    def empty(self):
        ''' Return true iff there are no elements in the histogram '''
        return len(self) == 0
    
    def draw(self):
        ''' Draw an element from the histogram '''
        if len(self) == 0:
            raise RuntimeError("Attempt to draw from empty empirical distribution")
        if self.changed:
            self.gen = WeightedRandomGenerator(self.hist)
            self.changed = False
        return self.gen.next()
    
    def add(self, key):
        ''' Add new element to the histogram '''
        self.hist[key] += 1
        self.tot += 1
        self.changed = True

class SimulatorWrapper(object):
    
    ''' Wrapper that sends requests to the Simulator but uses information
        gathered during alignment so far to select such parameters as read
        length, concordant/discordant fragment length, etc.
        
        With each call to next(), we either simulate a paired-end read or an
        unpaired read.
        '''
    
    def __init__(self, sim, m1fw, m2fw, dists):
        self.sim  = sim    # sequence simulator
        self.m1fw = m1fw   # whether mate 1 is revcomped w/r/t fragment
        self.m2fw = m2fw   # whether mate 2 is revcomped w/r/t fragment
        self.dists = dists # empirical distributions
    
    def nextUnpaired(self):
        ''' Simulate an unpaired read '''
        scDraw = self.dists.scDistUnp.draw()
        _, _, _, rfAln, sc = scDraw
        rl = len(rfAln) - rfAln.count('-')
        refid, refoff, fw, seq = self.sim.sim(rl) # simulate it
        assert rl == len(seq)
        read = Read.fromSimulator(seq, None, refid, refoff, fw, sc, "Unp")
        mutate(read, fw, scDraw) # mutate unpaired read
        assert read.qual is not None
        return read
    
    def nextPair(self):
        ''' Simulate a paired-end read '''
        fl, rl1, rl2 = 0, 1, 1
        sc1Draw, sc2Draw = None, None
        while fl < rl1 or fl < rl2:
            fl = self.dists.fragDist.draw()
            sc1Draw, sc2Draw = self.dists.scDistM12.draw()
            _, _, _, rfAln1, _ = sc1Draw
            _, _, _, rfAln2, _ = sc2Draw
            rl1 = len(rfAln1) - rfAln1.count('-')
            rl2 = len(rfAln2) - rfAln2.count('-')
        refid, refoff, fw, seq = self.sim.sim(fl) # simulate fragment
        assert len(seq) == fl
        # get mates from fragment
        seq1, seq2 = seq[:rl1], seq[-rl2:]
        if not self.m1fw: seq1 = revcomp(seq1)
        if not self.m2fw: seq2 = revcomp(seq2)
        # Now we have the Watson offset for one mate or the other,
        # depending on which mate is to the left w/r/t Watson.
        refoff1, refoff2 = refoff, refoff
        if fw: refoff2 = refoff + fl - rl2
        else:  refoff1 = refoff + fl - rl1
        _, _, _, _, sc1 = sc1Draw
        _, _, _, _, sc2 = sc2Draw
        rdp1 = Read.fromSimulator(seq1, None, refid, refoff1, fw == self.m1fw, sc1, "Conc")
        rdp2 = Read.fromSimulator(seq2, None, refid, refoff2, fw == self.m2fw, sc2, "Conc")
        rdm1 = Read.fromSimulator(seq1, None, refid, refoff1, fw == self.m1fw, sc1, "M1")
        rdm2 = Read.fromSimulator(seq2, None, refid, refoff2, fw == self.m2fw, sc2, "M2")
        mutate(rdp1, fw == self.m1fw, sc1Draw)
        mutate(rdp2, fw == self.m2fw, sc2Draw)
        mutate(rdm1, fw == self.m1fw, sc1Draw)
        mutate(rdm2, fw == self.m2fw, sc2Draw)
        return rdp1, rdp2, rdm1, rdm2

class Input(object):
    ''' Class that parses reads from input files and yields the reads/pairs
        produced using generators '''
    
    @staticmethod
    def fastaParse(fh):
        ''' Parse a single FASTA-format read from given filehandle.  Return
            None if input is exhausted. '''
        lns = [ fh.readline().rstrip() for x in xrange(0, 2) ]
        orig = '\n'.join(lns) + '\n'
        if len(lns[0]) == 0: return None
        return Read(lns[0][1:], lns[1], None, orig)

    @staticmethod
    def fastqParse(fh):
        ''' Parse a single FASTQ-format read from given filehandle.  Return
            None if input is exhausted. '''
        lns = [ fh.readline().rstrip() for x in xrange(0, 4) ]
        orig = '\n'.join(lns) + '\n'
        if len(lns[0]) == 0: return None
        return Read(lns[0][1:], lns[1], lns[3], orig)
    
    def __init__(self, format="fastq", unpFns=None, m1Fns=None, m2Fns=None):
        self.format = format
        if format == "fastq":
            self.parse = self.fastqParse
        elif format == "fasta":
            self.parse = self.fastaParse
        else: raise RuntimeError("Bad input format: '%s'" % format)
        self.unpFns, self.m1Fns, self.m2Fns = unpFns, m1Fns, m2Fns
    
    def __iter__(self):
        ''' Generator for all the reads.  Yields pairs of Read objects, where
            second element is None for unpaired reads. '''
        # Yield all the unpaired reads first
        if self.unpFns is not None:
            for unpFn in self.unpFns:
                with gzip.open(unpFn, 'r') if unpFn.endswith('.gz') else open(unpFn, 'r') as unpFh:
                    while True:
                        rd = self.parse(unpFh)
                        if rd is not None: yield (rd, None)
                        else: break # next file
        # Yield all the paired-end reads
        if self.m1Fns is not None:
            assert self.m2Fns is not None
            for (m1Fn, m2Fn) in zip(self.m1Fns, self.m2Fns):
                with gzip.open(m1Fn, 'r') if m1Fn.endswith('.gz') else open(m1Fn, 'r') as m1Fh:
                    with gzip.open(m2Fn, 'r') if m2Fn.endswith('.gz') else open(m2Fn, 'r') as m2Fh:
                        while True:
                            rd1, rd2 = self.parse(m1Fh), self.parse(m2Fh)
                            if rd1 is not None: yield (rd1, rd2)
                            else: break # next pair of files

class Dists(object):
    
    ''' Encapsulates all distributions that we train with real data.
        We collect random subsets of qualities/edits or unpaired reads
        and same for paired-end mate 1s and mate 2s.  We also collect
        a fragment length distribution. '''
    
    def __init__(self):
        self.fragDist  = Dist()      # fragment length distribution
        self.scDistUnp = ScoreDist() # qualities/edits for unpaired
        self.scDistM12 = ScorePairDist() # qualities/edits for mate #1s
    
    def addFraglen(self, fraglen):
        ''' Add given fragment length to fragment length distribution '''
        if fraglen < 10000:
            self.fragDist.add(fraglen) # remember fragment length
    
    def addPair(self, al1, al2, ref):
        ''' Add concordant paired-end read alignment to the model '''
        self.scDistM12.add(al1, al2, ref)
    
    def addRead(self, al, ref):
        ''' Add unpaired read alignment to the model '''
        self.scDistUnp.add(al, ref)
    
    def hasPaired(self):
        ''' Return true iff at least one concordant paired-end read was
            added. '''
        return not self.scDistM12.empty()
    
    def hasUnpaired(self):
        ''' Return true iff at least one unpaired read was added. '''
        return not self.scDistUnp.empty()

def isCorrect(al):
    ''' Return true if the given alignment is both simulated and "correct" --
        i.e. in or very near it's true point of origin '''
    print simulators.isExtendedWgsim(al.name)
    if args.correct_chromosomes is not None:
        return al.refid in args.correct_chromosomes
    elif simulators.isExtendedWgsim(al.name):
        return simulators.correctExtendedWgsim(al, args.wiggle)
    else: return None

class Output(Thread):
    
    ''' Aligner output reader.  Reads SAM output from aligner, updates
        empirical distributions (our "model" for what the input reads look
        like), and looks for records that correspond to simulated reads.  If a
        record corresponds to a simulated read, its correctness is checked and
        an output tuple written '''
    
    def __init__(self, samQ, dataset, dists, ref, resultQ, samOfh=None, ival=10000):
        Thread.__init__(self)
        self.samQ = samQ             # SAM records come from here
        self.samOfh = samOfh         # normal (non-simulated) SAM recs go here
        self.dataset = dataset       # dataset
        self.dists = dists           # distributions
        self.ref = ref               # reference sequences
        self.scDiffs = defaultdict(int) # histogram of expected-versus-observed score differences
        self.ival = ival
        self.resultQ = resultQ
        self.deamon = True
    
    def run(self):
        ''' Collect SAM output '''
        try:
            lastAl = None
            nal, nunp, npair = 0, 0, 0
            # Following loop involves maintaining 'lastAl' across
            # iterations so that we can match up the two ends of a pair
            while True:
                try:
                    ln = self.samQ.get(block=True, timeout=0.5)
                    if ln is None: break
                except Empty:
                    continue # keep trying
                
                # Send SAM to SAM output filehandle
                if self.samOfh is not None:
                    self.samOfh.write(ln)
                
                # Headers are skipped
                if ln[0] == '@': continue
                
                # Alignment parses the SAM line, but assumes that output is
                # from Bowtie 2
                al = AlignmentBowtie2(ln)
                nal += 1
                if al.paired: npair += 1
                else: nunp += 1
                if (nal % self.ival) == 0:
                    logging.info('      # alignments parsed: %d (%d paired, %d unpaired)' % (nal, npair, nunp))
                
                # If this is one mate from a concordantly-aligned pair,
                # match this mate up with its opposite (if we've seen it).
                # Note: Depends on Bowtie2-only 'YT:Z' flag
                mate1, mate2 = None, None
                if al.isAligned() and lastAl is not None and lastAl.alType == "CP":
                    mate1, mate2 = al, lastAl
                    assert mate1.mate1 != mate2.mate1
                    if mate2.mate1:
                        mate1, mate2 = mate2, mate1
                
                # Add the alignment to a Dataset data structure.  If it's an
                # aligned training read, check whether alignment is correct. 
                correct = None
                if self.dataset is not None:
                    isTraining = (al.name[0] == '!' and al.name.startswith('!!ts!!'))
                    if al.isAligned() and isTraining:
                        _, refid, fw, refoff, sc, trainingNm = string.split(al.name, '!!ts-sep!!')
                        sc, refoff = int(sc), int(refoff)
                        scDiff = sc - al.bestScore
                        self.scDiffs[scDiff] += 1
                        correct = False
                        # Check reference id, orientation
                        if refid == al.refid and fw == al.orientation():
                            # Check offset
                            correct = abs(refoff - al.pos) < args.wiggle
                        assert trainingNm in ['Unp', 'M1', 'M2', 'Conc']
                        # Add to training dataset
                        if trainingNm == "Unp":
                            self.dataset.addUnp(al, correct)
                        elif trainingNm[0] == "M":
                            self.dataset.addM(al, correct)
                        if mate1 is not None:
                            assert trainingNm == "Conc"
                            correct1, correct2 = correct, lastCorrect
                            if (lastAl.flags & 64) != 0:
                                correct1, correct2 = correct2, correct1
                            self.dataset.addPaired(mate1, mate2, abs(mate1.tlen), correct1, correct2)
                    elif al.isAligned():
                        # Test data
                        if mate1 is not None:
                            assert al.alType == 'CP'
                            correct1, correct2 = isCorrect(mate1), isCorrect(mate2)
                            self.dataset.addPaired(mate1, mate2, abs(mate1.tlen), correct1, correct2)
                        elif al.alType == 'UU':
                            self.dataset.addUnp(al, isCorrect(al))
                        elif al.alType in ['UP', 'DP']:
                            print ln
                            self.dataset.addM(al, isCorrect(al))
                
                # Add to our input-read summaries
                if self.dists is not None:
                    if al.isAligned():
                        if al.rnext == "=" and al.mateMapped():
                            # Can't easily get read length for opposite mate, so
                            # just use length of this mate as a surrogate
                            fraglen = abs(al.pnext - al.pos) + len(al.seq)
                            self.dists.addFraglen(fraglen)
                        try:
                            if mate1 is not None:
                                self.dists.addPair(mate1, mate2, self.ref)
                            else:
                                self.dists.addRead(al, self.ref)
                        except ReferenceOOB: pass
                
                if mate1 is None:
                    lastAl = al
                    lastCorrect = correct
                else:
                    lastAl = None
                    lastCorrect = None
                
                self.resultQ.put(True)
        except:
            self.resultQ.put(False)

def createInput():
    ''' Return an Input object that reads all user-provided input reads '''
    assert args.fasta or args.fastq
    return Input(format="fastq" if args.fastq else "fasta", unpFns=args.U,
                 m1Fns=args.m1, m2Fns=args.m2)

def go(args, alignerArgs):
    ''' Main driver for tandem simulator '''
    
    import tempfile
    
    random.seed(args.seed)
    
    # Set up logger
    logging.basicConfig(\
        format='%(asctime)s:%(levelname)s:%(message)s',
        datefmt='%m/%d/%y-%H:%M:%S',
        level=logging.DEBUG if args.verbose else logging.INFO)
    
    # Create output directory if needed
    if not os.path.isdir(args.output_directory):
        try: os.makedirs(args.output_directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    
    # Construct command to invoke aligner
    if args.aligner == 'bowtie2':
        align_cmd = 'bowtie2 '
        if args.bt2_exe is not None:
            align_cmd = args.bt2_exe + " "
        align_cmd += ' '.join(alignerArgs)
        align_cmd += ' --reorder --sam-no-qname-trunc --mapq-extra'
    elif args.aligner == 'bwa-mem':
        align_cmd = 'bwa mem '
        if args.bwa_exe is not None:
            align_cmd = args.bwa_exe + " "
        align_cmd += ' '.join(alignerArgs)
        align_cmd += ' --reorder --sam-no-qname-trunc --mapq-extra'
    elif args.aligner == 'bwa-aln':
        align_cmd = 'bwa aln '
        if args.bwa_exe is not None:
            align_cmd = args.bwa_exe + " "
        align_cmd += ' '.join(alignerArgs)
        align_cmd += ' --reorder --sam-no-qname-trunc --mapq-extra'
    
    # Timing measurements
    setupIval, al1Ival, al2Ival = None, None, None
    
    st = time.clock()
    
    logging.info('Loading reference data')
    indexed = True
    refClass = ReferenceIndexed if indexed else ReferencePicklable
    with refClass(args.ref) as ref:
    
        # ##################################################
        # ALIGN REAL DATA (or read in alignments from SAM)
        # ##################################################
        
        bt2 = Bowtie2(align_cmd)
        samFn = os.path.join(args.output_directory, 'input.sam')
        inputSamFh = open(samFn, 'w') if (args.write_input_sam or args.write_all) else None
        testData = Dataset()
        dists = Dists()
        
        logging.info('Real-data Bowtie 2 command: "%s"' % align_cmd)
        
        # Create the thread that eavesdrops on output from bowtie2
        resultQ = Queue()
        othread = Output(
            bt2.outQ,          # SAM queue
            testData,          # Dataset to gather alignments into
            dists,             # empirical dists
            ref,               # reference genome
            resultQ,           # result queue
            samOfh=inputSamFh) # SAM output filehandle
        othread.daemon = True
        othread.start()
        
        # Stop timing setup
        setupIval = time.clock() - st
        st = time.clock()
        
        logging.info('Initializing threads, queues and FIFOs')
        
        # Read through all the input read files and direct all reads to the
        # appropriate queue
        upto = args.upto or sys.maxint
        numReads = 0
        for (rd1, rd2) in iter(createInput()):
            bt2.inQ.put(Read.toTab6(rd1, rd2, truncateName=True) + '\n')
            numReads += 1
            if numReads >= upto: break
        bt2.inQ.put(None)
        
        logging.debug('Waiting for bt2 to finish')
        while bt2.pipe.poll() is None:
            time.sleep(0.5)
        while othread.is_alive():
            othread.join(0.5)
        logging.debug('bt2 process and output thread finished')
        othreadResult = resultQ.get()
        if not othreadResult:
            raise RuntimeError('Aligner output monitoring thread encountered error; aborting...')
        
        logging.info('Finished aligning input reads')
        if inputSamFh is not None:
            logging.info('  Input read alignments written to "%s"' % samFn)
        
        # Writing test dataset
        if args.write_test_data or args.write_all:
            if False:
                testPickleFn = os.path.join(args.output_directory, 'test.pickle' + ('.gz' if args.compress_output else ''))
                testData.save(testPickleFn, compress=args.compress_output)
                logging.info('Test data written to "%s"' % testPickleFn)
            
            testCsvFnPrefix = os.path.join(args.output_directory, 'test')
            testData.saveCsvs(testCsvFnPrefix, compress=args.compress_output)
            logging.info('Test data (CSV format) written to "%s*"' % testCsvFnPrefix)
        
        # Stop timing interval for alignment phase 1
        al1Ival = time.clock() - st
        st = time.clock()
        
        # ##################################################
        # ALIGN SIMULATED DATA
        # ##################################################
        #
        # Open Bowtie 2 w/ same arguments.  Give it only simulated reads.
        
        logging.info('Opening Bowtie 2 process for training reads')
        
        # We're potentially going to send four different types of simulated
        # reads to Bowtie 2, depending on what type of training data the read
        # pertains to.
        bt2 = Bowtie2(align_cmd)
        samFn = os.path.join(args.output_directory, 'training.sam')
        trainingSamFh = open(samFn, 'w') if (args.write_training_sam or args.write_all) else None
        trainingData = Dataset()
        
        # Construct sequence and quality simulators
        logging.info('  Creating sequence simulator')
        sim = SequenceSimulator(ref)
        simw = SimulatorWrapper(\
            sim,             # sequence simulator
            not args.m1rc,   # whether m1 is revcomped w/r/t fragment
            args.m2fw,       # whether m2 is revcomped w/r/t fragment
            dists)           # empirical distribution
        
        # Create thread that eavesdrops on output from bowtie2 with simulated input
        logging.info('  Opening output-parsing thread')
        resultQ = Queue()
        othread = Output(\
            bt2.outQ,        # SAM input filehandle
            trainingData,    # Training data
            None,            # qualities/edits for unpaired
            ref,             # reference genome
            resultQ,         # result queue
            samOfh=trainingSamFh)
        othread.daemon = True
        othread.start()
        
        # Simulate concordant paired-end reads from empirical distributions
        trainingPairFn, trainingUnpFn = None, None
        trainingPairFh, trainingUnpFh = None, None
        writeTrainingReads = args.write_training_reads or args.write_all
        if writeTrainingReads:
            if dists.hasPaired():
                trainingPairFn = os.path.join(args.output_directory, 'training_pairs.tab6')
                trainingPairFh = open(trainingPairFn, 'w')
            if dists.hasUnpaired():
                trainingUnpFn = os.path.join(args.output_directory, 'training_unpaired.tab6')
                trainingUnpFh = open(trainingUnpFn, 'w')
        
        logging.info('  Simulating %d paired-end reads' % (args.num_reads/2))
        if dists.hasPaired():
            for i in xrange(0, args.num_reads/2):
                if (i+1 % 100) == 0:
                    logging.info('    Simulating paired-end read %d' % i)
                rdp1, rdp2, rdm1, rdm2 = simw.nextPair()
                if writeTrainingReads:
                    trainingPairFh.write(Read.toTab6(rdp1, rdp2))
                    trainingPairFh.write('\n')
                bt2.inQ.put(Read.toTab6(rdp1, rdp2, truncateName=True) + '\n')
                bt2.inQ.put(Read.toTab6(rdm1, truncateName=True) + '\n')
                bt2.inQ.put(Read.toTab6(rdm2, truncateName=True) + '\n')
            # Signal that we're done supplying paired-end reads
            if writeTrainingReads:
                trainingPairFh.close()
        
        # Simulate unpaired reads from empirical distributions
        logging.info('  Simulating %d unpaired reads' % args.num_reads)
        if dists.hasUnpaired():
            for i in xrange(0, args.num_reads):
                if (i+1 % 100) == 0:
                    logging.info('    Simulating unpaired read %d' % i)
                rd = simw.nextUnpaired()
                if writeTrainingReads:
                    trainingUnpFh.write(Read.toTab6(rd))
                    trainingUnpFh.write('\n')
                bt2.inQ.put(Read.toTab6(rd, truncateName=True) + '\n')
            if writeTrainingReads:
                trainingUnpFh.close()
        
        bt2.inQ.put(None)
        
        # Signal that we're done supplying unpaired reads
        logging.debug('Waiting for bt2 to finish')
        while bt2.pipe.poll() is None:
            time.sleep(0.5)
        while othread.is_alive():
            othread.join(0.5)
        logging.debug('bt2 process finished')
        othreadResult = resultQ.get()
        if not othreadResult:
            raise RuntimeError('Aligner output monitoring thread encountered error; aborting...')
        logging.info('Finished simulating and aligning training reads')
        if trainingPairFh is not None:
            logging.info('  Paired training reads written to "%s"' % trainingPairFn)
        if trainingUnpFh is not None:
            logging.info('  Unpaired training reads written to "%s"' % trainingUnpFn)
        
        # Stop timing interval for alignment phase 2
        al2Ival = time.clock() - st
    
        logging.info('Score difference (expected - actual) histogram:')
        for k, v in sorted(othread.scDiffs.iteritems()):
            logging.info('  %d: %d' % (k, v))
        
        if args.write_training_data or args.write_all:
            # Writing training data
            if False:
                trainingPickleFn = os.path.join(args.output_directory, 'training.pickle' + ('.gz' if args.compress_output else ''))
                testData.save(testFn, compress=args.compress_output)
                trainingData.save(trainingPickleFn, compress=args.compress_output)
                logging.info('Training data written to "%s"' % trainingPickleFn)
            
            trainingCsvFnPrefix = os.path.join(args.output_directory, 'training')
            trainingData.saveCsvs(trainingCsvFnPrefix, compress=args.compress_output)
            logging.info('Training data (CSV format) written to "%s*"' % trainingCsvFnPrefix)
        
        if setupIval is not None:
            logging.info('Setup running time: %0.3f secs' % setupIval)
        if al1Ival is not None:
            logging.info('Alignment (real reads): %0.3f secs' % al1Ival)
        if al2Ival is not None:
            logging.info('Alignment (simulated reads): %0.3f secs' % al2Ival)

if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser(\
        description='Evaluate the sensitivity of an alignment tool w/r/t given'
                    'genome and set of alignment parameters.')
    
    parser.add_argument(\
        '--ref', metavar='path', type=str, nargs='+', required=True,
        help='FASTA file(s) containing reference genome sequences')
    parser.add_argument(\
        '--pickle-ref', metavar='path', type=str,
        help='Pickle FASTA input for speed, or use pickled version if it '
             'exists already.  Pickled version is stored at given path')
    parser.add_argument(\
        '--U', metavar='path', type=str, nargs='+', help='Unpaired read files')
    parser.add_argument(\
        '--m1', metavar='path', type=str, nargs='+', help='Mate 1 files; must be specified in same order as --m2')
    parser.add_argument(\
        '--m2', metavar='path', type=str, nargs='+', help='Mate 2 files; must be specified in same order as --m1')
    parser.add_argument(\
        '--m1rc', action='store_const', const=True, default=False,
        help='Set if mate 1 is reverse-complemented w/r/t the fragment')
    parser.add_argument(\
        '--m2fw', action='store_const', const=True, default=False,
        help='Set if mate 2 is oriented forward w/r/t the fragment')
    parser.add_argument(\
        '--fasta', action='store_const', const=True, default=False, help='Input reads are FASTA')
    parser.add_argument(\
        '--fastq', action='store_const', const=True, default=True, help='Input reads are FASTQ')
    parser.add_argument(\
        '--seed', metavar='int', type=int, default=99099,
        required=False, help='Integer to initialize pseudo-random generator')
    parser.add_argument(\
        '--num-reads', metavar='int', type=int, default=1000,
        required=False, help='Number of reads to simulate')
    parser.add_argument(\
        '--upto', metavar='int', type=int, default=None,
        required=False, help='Stop after this many input reads')
    parser.add_argument(\
        '--wiggle', metavar='int', type=int, default=30,
        required=False, help='Wiggle room to allow in starting position when determining whether alignment is correct')
    
    parser.add_argument(\
        '--sam-input', metavar='path', type=str,
        help='Input SAM file to apply training data to.  Use with --training-input.')
    parser.add_argument(\
        '--bt2-exe', metavar='path', type=str, help='Path to Bowtie 2 exe')
    parser.add_argument(\
        '--bwa-exe', metavar='path', type=str, help='Path to BWA exe')
    parser.add_argument(\
        '--aligner', metavar='name', default='bowtie2', type=str, help='bowtie2 or bwa-mem')
    
    # Some basic flags
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False, help='Do unit tests')
    parser.add_argument(\
        '--profile', action='store_const', const=True, default=False, help='Print profiling info')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False, help='Be talkative')
    parser.add_argument(\
        '--version', action='store_const', const=True, default=False, help='Print version and quit')
    
    # For when input is itself simulated, so we can output a Dataset with the
    # 'correct' column filled in properly 
    parser.add_argument(\
        '--input-from-mason', action='store_const', const=True, default=False, help='Input reads were simulated from Mason')
    parser.add_argument(\
        '--input-from-wgsim', action='store_const', const=True, default=False, help='Input reads were simulated from wgsim')
    parser.add_argument(\
        '--input-from-grinder', action='store_const', const=True, default=False, help='Input reads were simulated from Grinder')
    parser.add_argument(\
        '--correct-chromosomes', metavar='list', type=str, nargs='+',
        help='Label test data originating from any of these '
             'chromosomes as "correct."  Useful for tests on '
             'real-world data where it is known that the data came '
             'from a parituclar chromosome.')

    # Output file-related arguments
    parser.add_argument(\
        '--output-directory', metavar='path', type=str, required=True, help='Write outputs to this directory')
    parser.add_argument(\
        '--write-input-sam', action='store_const', const=True,
        default=False, help='Write SAM alignments for the real input reads to "input.sam" in output directory')
    parser.add_argument(\
        '--write-training-reads', action='store_const', const=True,
        default=False, help='Write FASTQ for the training reads to "training.fastq" in output directory')
    parser.add_argument(\
        '--write-test-data', action='store_const', const=True,
        default=False, help='Write Dataset object for training data.  "Correct" column set to all None\'s.')
    parser.add_argument(\
        '--write-training-sam', action='store_const', const=True,
        default=False, help='Write SAM alignments for the training reads to "training.sam" in output directory')
    parser.add_argument(\
        '--write-training-data', action='store_const', const=True,
        default=False, help='Write Dataset object for training data.')
    parser.add_argument(\
        '--write-all', action='store_const', const=True,
        default=False, help='Same as specifying all --write-* options')
    parser.add_argument(\
        '--compress-output', action='store_const', const=True,
        default=False, help='gzip all output files')
    
    argv = sys.argv
    alignerArgs = []
    in_args = False
    for i in xrange(1, len(sys.argv)):
        if in_args:
            alignerArgs.append(sys.argv[i])
        if sys.argv[i] == '--':
            argv = sys.argv[:i]
            in_args = True
    
    args = parser.parse_args(argv[1:])
    
    if args.test:
        
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
    elif args.profile:
        import cProfile
        cProfile.run('go(args, alignerArgs)')
    else:
        go(args, alignerArgs)
