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

'''

__author__ = "Ben Langmead"
__email__ = "langmea@cs.jhu.edu"

import os
import sys
import re
import threading
import math
import string
import random
import bisect
import subprocess
import tempfile
import signal
import traceback
import cPickle
import time
import logging
from collections import defaultdict
from Queue import Queue
from sam import cigarToList, mdzToList, cigarMdzToStacked

def quit_handler(signum,frame):
    traceback.print_stack()

signal.signal(signal.SIGQUIT,quit_handler)

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

class Read(object):
    ''' Encapsulates a read '''
    
    def __init__(self, name, seq, qual, orig=None):
        ''' Initialize new read given name, sequence, quality '''
        self.name = name
        self.seq = seq
        self.qual = qual
        self.orig = orig
        assert self.repOk()
    
    @classmethod
    def fromSimulator(cls, seq, qual, refid, refoff, fw, sc, trainingNm):
        ''' Construct appropriate read object (with appropriate name) given
            some simulated properties of the read. '''
        rdname = "!!ts-sep!!".join(["!!ts!!", refid, "+" if fw else "-", str(refoff), str(sc), trainingNm])
        return cls(rdname, seq, qual)
    
    @classmethod
    def toTab6(cls, rd1, rd2=None):
        ''' Convert either an unpaired read or a pair of reads to tab6
            format '''
        if rd2 is not None:
            return "\t".join([rd1.name, rd1.seq, rd1.qual, rd2.name, rd2.seq, rd2.qual])
        return "\t".join([rd1.name, rd1.seq, rd1.qual])
    
    def __len__(self):
        ''' Return number of nucleotides in read '''
        return len(self.seq)
    
    def __str__(self):
        ''' Return string representation '''
        if self.orig is not None:
            return self.orig # original string preferred
        elif self.qual is not None:
            assert len(self.seq) == len(self.qual)
            return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
        else:
            return ">%s\n%s\n" % (self.name, self.seq)
    
    def repOk(self):
        ''' Check that read is internally consistent '''
        if self.qual is not None:
            assert len(self.seq) == len(self.qual)
        return True

class Alignment(object):
    ''' Encapsulates an alignment record for a single aligned read.  Has
        facilities for parsing certain important SAM extra fields output by
        Bowtie 2. '''
    
    __asRe  = re.compile('AS:i:([-]?[0-9]+)') # best score
    __xsRe  = re.compile('XS:i:([-]?[0-9]+)') # second-best score
    __ytRe  = re.compile('YT:Z:([A-Z]+)')     # alignment type
    __ysRe  = re.compile('YS:i:([-]?[0-9]+)') # score of opposite
    __mdRe  = re.compile('MD:Z:([^\s]+)')     # MD:Z string
    __xlsRe = re.compile('Xs:i:([-]?[0-9]+)') # 3rd best
    __zupRe = re.compile('ZP:i:([-]?[0-9]+)') # best concordant
    __zlpRe = re.compile('Zp:i:([-]?[0-9]+)') # 2nd best concordant
    __kNRe  = re.compile('YN:i:([-]?[0-9]+)') # min valid score
    __knRe  = re.compile('Yn:i:([-]?[0-9]+)') # max valid score
    
    def __init__(self, ln):
        ''' Parse ln, which is a line of SAM output from Bowtie 2.  The line
            must correspond to an aligned read. '''
        self.name, self.flags, self.refid, self.pos, self.mapq, self.cigar, \
        self.rnext, self.pnext, self.tlen, self.seq, self.qual, self.extra = \
        string.split(ln, '\t', 11)
        assert self.flags != "*"
        assert self.pos != "*"
        assert self.mapq != "*"
        assert self.tlen != "*"
        self.flags = int(self.flags)
        self.pos = int(self.pos)
        self.mapq = int(self.mapq)
        self.tlen = int(self.tlen)
        self.pnext = int(self.pnext)
        self.fw = (self.flags & 16) == 0
        self.mate1 = (self.flags & 64) != 0
        self.mate2 = (self.flags & 128) != 0
        self.paired = self.mate1 or self.mate2
        self.exFields = {}
        for tok in string.split(self.extra, '\t'):
            idx1 = tok.find(':')
            assert idx1 != -1
            idx2 = tok.find(':', idx1+2)
            assert idx2 != -1
            self.exFields[tok[:idx2]] = tok[idx2+1:].rstrip()
        # Parse AS:i
        se = Alignment.__asRe.search(self.extra)
        self.bestScore = None
        if se is not None:
            self.bestScore = int(se.group(1))
        # Parse XS:i
        se = Alignment.__xsRe.search(self.extra)
        self.secondBestScore = None
        if se is not None:
            self.secondBestScore = int(se.group(1))
        # Parse Xs:i
        se = Alignment.__xlsRe.search(self.extra)
        self.thirdBestScore = None
        if se is not None:
            self.thirdBestScore = int(se.group(1))
        # Parse ZP:i
        se = Alignment.__zupRe.search(self.extra)
        self.bestConcordantScore = None
        if se is not None:
            self.bestConcordantScore = int(se.group(1))
        # Parse Zp:i
        se = Alignment.__zlpRe.search(self.extra)
        self.secondBestConcordantScore = None
        if se is not None:
            self.secondBestConcordantScore = int(se.group(1))
        # Parse YT:Z
        se = Alignment.__ytRe.search(self.extra)
        self.alType = None
        if se is not None:
            self.alType = se.group(1)
        # Parse YS:i
        se = Alignment.__ysRe.search(self.extra)
        self.mateBest = None
        if se is not None:
            self.mateBest = int(se.group(1))
        # Parse MD:Z
        self.mdz = None
        se = Alignment.__mdRe.search(self.extra)
        if se is not None:
            self.mdz = se.group(1)
        # Parse YN:i
        self.minValid = None
        se = Alignment.__kNRe.search(self.extra)
        if se is not None:
            self.minValid = se.group(1)
        # Parse Yn:i
        self.maxValid = None
        se = Alignment.__knRe.search(self.extra)
        if se is not None:
            self.maxValid = se.group(1)
        self.concordant = self.alType == "CP"
        self.discordant = self.alType == "DP"
        self.unpPair = self.alType == "UP"
        self.unp = self.alType == "UU"
        assert self.repOk()
    
    def isAligned(self):
        ''' Return true iff read aligned '''
        return (self.flags & 4) == 0
    
    def orientation(self):
        ''' Return orientation as + or - '''
        if (self.flags & 16) != 0:
            return "-"
        else:
            return "+"
    
    def mateMapped(self):
        ''' Return true iff opposite mate aligned '''
        return (self.flags & 8) == 0
    
    def fragmentLength(self):
        ''' Return fragment length '''
        return abs(self.tlen)
    
    def __len__(self):
        ''' Return read length '''
        return len(self.seq)
    
    def repOk(self):
        ''' Check alignment for internal consistency '''
        assert self.alType is not None
        assert self.paired or self.fragmentLength() == 0
        assert not self.isAligned() or self.bestScore is not None
        assert self.alType in ["CP", "DP", "UP", "UU"]
        return True

class WeightedRandomGenerator(object):
    
    ''' Given an ordered list of weights, generate with each call to next() an
        offset into the list of the weights with probability equal to the
        fraction of the total weight. '''
    
    def __init__(self, weights):
        self.totals = []
        running_total = 0
        for w in iter(weights):
            running_total += w
            self.totals.append(running_total)
    
    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

class ReservoirSampler(object):
    ''' Simple reservoir sampler '''
    
    def __init__(self, k):
        ''' Initialize given k, the size of the reservoir '''
        self.k = k
        self.r = []
        self.n = 0
    
    def add(self, obj):
        ''' Add object to sampling domain '''
        if self.n < self.k:
            self.r.append(obj)
        else:
            j = random.randint(0, self.n)
            if j < self.k:
                self.r[j] = obj
        self.n += 1
    
    def draw(self):
        ''' Make a draw from the reservoir '''
        return random.choice(self.r)
    
    def __len__(self):
        ''' Return number of items added to the domain '''
        return self.n
    
    def empty(self):
        ''' Return true iff no items have been added '''
        return len(self) == 0

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
    
    def add(self, al):
        ''' Convert given alignment to a tuple and add it to the reservoir
            sampler. '''
        # Extract quality string
        assert al.cigar is not None
        assert al.mdz is not None
        # Extract CIGAR, MD:Z
        cigarList = cigarToList(al.cigar)
        mdzList = mdzToList(al.mdz)
        sc = al.bestScore
        # Get stacked alignment
        try:
            rdAln, rfAln = cigarMdzToStacked(al.seq, cigarList, mdzList)
        except AssertionError:
            logger.error("cigar=%s, mdz=%s" % (al.cigar, al.mdz))
            raise AssertionError
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
    
    def add(self, al1, al2):
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
        rdAln1, rfAln1 = cigarMdzToStacked(al1.seq, cigarList1, mdzList1)
        rdAln2, rfAln2 = cigarMdzToStacked(al2.seq, cigarList2, mdzList2)
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
            self.gen = WeightedRandomGenerator(self.hist.itervalues())
            self.changed = False
        return self.hist.keys()[self.gen.next()]
    
    def add(self, key):
        ''' Add new element to the histogram '''
        self.hist[key] += 1
        self.tot += 1
        self.changed = True

class SequenceSimulator(object):
    ''' Class that, given a collection of FASTA files, samples intervals of
        specified length from the strings contained in them.  Simulated reads
        are not permitted to overlap a non-A/C/G/T character in the reference.
    '''
    
    def __init__(self, fafns, pickleFn, idx_fn=None, verbose=False):
        self.__re = re.compile('[^ACGTacgt]')
        self.refs, self.names, self.lens = {}, [], []
        pickleExists = False
        if pickleFn is not None:
            pickleExists = os.path.exists(pickleFn)
        if pickleFn is not None and pickleExists:
            self.refs, self.names, self.lens = cPickle.load(open(pickleFn, 'rb'))
        else:
            totlen = 0
            pt_sz = 50000000
            last_pt = 0
            max_bases = sys.maxint
            if args.max_bases is not None:
                max_bases = args.max_bases
            abort = False
            for fafn in fafns:
                fafh = open(fafn, 'r')
                name = None
                for line in fafh:
                    line = line.rstrip()
                    ln = len(line)
                    if ln > 0 and line[0] == '>':
                        ind = line.find(" ")
                        if ind == -1: ind = len(line)
                        line = line[1:ind]
                        name = line
                        self.refs[name] = []
                        self.names.append(name)
                        self.lens.append(0)
                    else:
                        assert name is not None
                        self.refs[name].append(line)
                        self.lens[-1] += ln
                        totlen += ln
                        if verbose:
                            pt = totlen / pt_sz
                            if pt > last_pt:
                                logging.info("Read %d FASTA bytes..." % totlen)
                            last_pt = pt
                        if totlen > max_bases:
                            abort = True
                            break
                fafh.close()
                if abort:
                    break
            for k in self.refs.iterkeys():
                self.refs[k] = ''.join(self.refs[k])
        if pickleFn is not None and not pickleExists:
            cPickle.dump((self.refs, self.names, self.lens), open(pickleFn, 'wb'), cPickle.HIGHEST_PROTOCOL)
        self.rnd = WeightedRandomGenerator(self.lens)
    
    def sim(self, ln, verbose=False):
        ''' Simulate a read '''
        # Pick a reference sequence in a weighted random fashon
        refi = self.rnd.next()
        assert refi < len(self.names)
        fw = True
        nm = self.names[refi] # reference name
        refoff = random.randint(0, self.lens[refi] - ln) # pick offset
        seq = self.refs[nm][refoff:refoff+ln] # extract substring
        # Simulated read can't overlap non-A-C-G-T character in reference
        while self.__re.search(seq):
            refoff = random.randint(0, self.lens[refi] - ln) # pick new offset
            seq = self.refs[nm][refoff:refoff+ln] # extract substring again
        seq = seq.upper()
        if random.random() > 0.5: # possibly reverse-complement
            fw = False
            seq = revcomp(seq) # reverse complement
        assert not "N" in seq
        return (nm, refoff, fw, seq) # return ref id, ref offset, orientation, sequence

def mutate(rd, rdfw, scDistDraw):
    ''' Given a read that already has the appropriate length (i.e. equal to #
        characters on the reference side of the alignment), take the alignment
        information contained in the scDistDraw object and modify rd to
        contain the same pattern of edits.  Modifies rd in place. '''
    fw, qual, rdAln, rfAln, sc = scDistDraw
    assert 'N' not in rd.seq
    assert len(rd.seq) == len(rfAln) - rfAln.count('-')
    if rdfw != fw:
        qual, rdAln, rfAln = qual[::-1], rdAln[::-1], rfAln[::-1]
    rd.qual = qual # Install qualities
    # Walk along stacked alignment, making corresponding changes to
    # read sequence
    seq = []
    i, rdi, rfi = 0, 0, 0
    while i < len(rdAln):
        if rdAln[i] == '-':
            rfi += 1
        elif rfAln[i] == '-':
            seq.append(rdAln[i])
            rdi += 1
        elif rdAln[i] != rfAln[i] and rdAln[i] == 'N':
            seq.append('N')
            rfi += 1; rdi += 1
        elif rdAln[i] != rfAln[i]:
            assert rfi < len(rd.seq)
            oldc = rd.seq[rfi].upper()
            cs = ['A', 'C', 'G', 'T']
            assert oldc in cs, "oldc was %s" % oldc
            cs.remove(oldc)
            newc = random.choice(cs)
            seq.append(newc)
            rfi += 1; rdi += 1
        else:
            assert rfi < len(rd.seq)
            seq.append(rd.seq[rfi])
            rfi += 1; rdi += 1
        i += 1
    rd.seq = ''.join(seq)
    assert len(rd.seq) == len(rd.qual)

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
                with open(unpFn, 'r') as unpFh:
                    while True:
                        rd = self.parse(unpFh)
                        if rd is not None: yield (rd, None)
                        else: break # next file
        # Yield all the paired-end reads
        if self.m1Fns is not None:
            assert self.m2Fns is not None
            for (m1Fn, m2Fn) in zip(self.m1Fns, self.m2Fns):
                with open(m1Fn, 'r') as m1Fh:
                    with open(m2Fn, 'r') as m2Fh:
                        while True:
                            rd1, rd2 = self.parse(m1Fh), self.parse(m2Fh)
                            if rd1 is not None: yield (rd1, rd2)
                            else: break # next pair of files

class UnpairedTuple(object):
    ''' Unpaired training/test tuple. '''
    def __init__(self, rdlen, minv, maxv, bestsc, best2sc, mapq):
        self.rdlen = rdlen              # read len
        self.minv = minv                # min valid score
        self.maxv = maxv                # max valid score
        self.bestsc = bestsc            # best
        self.best2sc = best2sc          # 2nd-best score
        self.mapq = mapq                # original mapq
    
    @classmethod
    def fromAlignment(cls, al):
        ''' Create unpaired training/test tuple from Alignment object '''
        secbest = al.secondBestScore or al.thirdBestScore
        return cls(len(al), al.minValid, al.maxValid, al.bestScore, secbest, al.mapq)

class PairedTuple(object):
    ''' Concordant paired-end training/test tuple.  One per mate alignment. '''
    def __init__(self, rdlen1, minv1, maxv1, bestsc1, best2sc1, mapq1,
                 rdlen2, minv2, maxv2, bestsc2, best2sc2, mapq2,
                 bestconcsc, best2concsc, fraglen):
        self.rdlen1 = rdlen1            # read len
        self.rdlen2 = rdlen2            # read len
        self.minv1 = minv1              # min valid score
        self.minv2 = minv2              # min valid score
        self.maxv1 = maxv1              # max valid score
        self.maxv2 = maxv2              # max valid score
        self.bestsc1 = bestsc1          # best
        self.bestsc2 = bestsc2          # best
        self.best2sc1 = best2sc1        # 2nd-best score
        self.best2sc2 = best2sc2        # 2nd-best score
        self.mapq1 = mapq1              # original mapq
        self.mapq2 = mapq2              # original mapq
        self.bestconcsc = bestconcsc    # best concordant
        self.best2concsc = best2concsc  # 2nd-best concordant
        self.fraglen = fraglen          # fragment length
    
    @classmethod
    def fromAlignments(cls, al1, al2, fraglen):
        ''' Create unpaired training/test tuple from pair of Alignments '''
        secbest1 = al1.secondBestScore or al1.thirdBestScore
        secbest2 = al2.secondBestScore or al2.thirdBestScore
        return cls(len(al1), al1.minValid, al1.maxValid, al1.bestScore,
                   secbest1, al1.mapq,
                   len(al2), al2.minValid, al2.maxValid, al2.bestScore,
                   secbest2, al2.mapq,
                   al1.bestConcordantScore,
                   al1.secondBestConcordantScore, al1.fragmentLength())

class Training(object):
    
    ''' Encapsulates all 'training data' and classifiers.  Training data =
        alignments for simulated reads. '''
    
    def __init__(self):
        # Training data for individual reads and mates.  Training tuples are
        # (rdlen, minValid, maxValid, bestSc, scDiff)
        self.trainUnp, self.labUnp = [], []
        self.trainM,   self.labM   = [], []
        # Training data for concordant pairs.  Training tuples are two tuples
        # like the one described above, one for each mate, plus the fragment
        # length.  Label says whether the first mate's alignment is correct.
        self.trainConc, self.labConc = [], []
    
    def __len__(self):
        ''' Return number of pieces of training data added so far '''
        return len(self.trainUnp) + len(self.trainM) + len(self.trainConc)
    
    def addPaired(self, al1, al2, fraglen, correct1, correct2):
        ''' Add a concordant paired-end alignment to our collection of
            training data. '''
        assert al1.concordant and al2.concordant
        rec1 = PairedTuple.fromAlignments(al1, al2, fraglen)
        rec2 = PairedTuple.fromAlignments(al2, al1, fraglen)
        for rec in [rec1, rec2]: self.trainConc.append(rec)
        self.labConc.extend([correct1, correct2])
    
    def addUnp(self, al, correct):
        ''' Add an alignment for a simulated unpaired read to our collection of
            training data. '''
        self.trainUnp.append(UnpairedTuple.fromAlignment(al))
        self.labUnp.append(correct)
    
    def addM(self, al, correct):
        ''' Add an alignment for a simulated mate that has been aligned in an
            unpaired fashion to our collection of training data. '''
        self.trainM.append(UnpairedTuple.fromAlignment(al))
        self.labM.append(correct)
    
    def save(self, fn, compress=True):
        ''' Save training data to a (possibly) compressed pickle file. '''
        import cPickle
        save = (self.trainUnp,  self.labUnp,
                self.trainM,    self.labM,
                self.trainConc, self.labConc)
        if compress:
            import gzip
            fh = gzip.open(fn, 'wb')
        else:
            fh = open(fn, 'wb')
        cPickle.dump(save, fh, cPickle.HIGHEST_PROTOCOL)
        fh.close()
    
    def load(self, fn, compress=True):
        ''' Load training data from a (possibly) compressed pickle file. '''
        import cPickle
        if compress:
            import gzip
            fh = gzip.open(fn, 'rb')
        else:
            fh = open(fn, 'rb')
        (self.trainUnp,  self.labUnp, \
         self.trainM,    self.labM, \
         self.trainConc, self.labConc) = cPickle.load(fh)
        fh.close()

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
    
    def addPair(self, al1, al2):
        ''' Add concordant paired-end read alignment to the model '''
        self.scDistM12.add(al1, al2)
    
    def addRead(self, al):
        ''' Add unpaired read alignment to the model '''
        self.scDistUnp.add(al)
    
    def hasPaired(self):
        ''' Return true iff at least one concordant paired-end read was
            added. '''
        return not self.scDistM12.empty()
    
    def hasUnpaired(self):
        ''' Return true iff at least one unpaired read was added. '''
        return not self.scDistUnp.empty()

class Output(threading.Thread):
    
    ''' Aligner output reader.  Reads SAM output from aligner, updates
        empirical distributions (our "model" for what the input reads look
        like), and looks for records that correspond to simulated reads.  If a
        record corresponds to a simulated read, its correctness is checked and
        an output tuple written '''
    
    def __init__(self, samIfh, samOfh, unalFh, training, dists, ival=1000):
        threading.Thread.__init__(self)
        self.samIfh = samIfh         # SAM records come from here
        self.samOfh = samOfh         # normal (non-simulated) SAM recs go here
        self.training = training     # training data store
        self.dists = dists           # distributions
        self.unalFh = unalFh         # write unaligned reads here
        self.scDiffs = {}            # histogram of expected-versus-observed score differences
        self.ival = ival
    
    def run(self):
        ''' Run the Bowtie 2 output collection thread '''
        lastAl = None
        nal, nunp, npair = 0, 0, 0
        for ln in self.samIfh:
            if ln[0] == '@':
                if self.samOfh is not None:
                    self.samOfh.write(ln) # header line
                continue
            al = Alignment(ln)
            nal += 1
            if al.paired: npair += 1
            else: nunp += 1
            if (nal % self.ival) == 0:
                logging.info('      # Bowtie 2 alignments parsed: %d (%d paired, %d unpaired)' % (nal, npair, nunp))
            nm, flags, refid, pos, _, _, _, _, _, seq, qual, _ = string.split(ln, '\t', 11)
            
            # If this is one mate from a concordantly-aligned pair,
            # match this mate up with its opposite (if we've seen it)
            mate1, mate2 = None, None
            if al.isAligned() and lastAl is not None and lastAl.alType == "CP":
                mate1, mate2 = al, lastAl
                assert mate1.mate1 != mate2.mate1
                if mate2.mate1:
                    mate1, mate2 = mate2, mate1
            
            correct = None
            if al.name[0] == '!' and al.name.startswith('!!ts!!'):
                # Simulated read
                assert self.training is not None
                if al.isAligned():
                    _, refid, fw, refoff, sc, trainingNm = string.split(al.name, '!!ts-sep!!')
                    sc, refoff = int(sc), int(refoff)
                    scDiff = sc - al.bestScore
                    self.scDiffs[scDiff] = self.scDiffs.get(scDiff, 0) + 1
                    correct = False
                    # Check reference id, orientation
                    if refid == al.refid and fw == al.orientation():
                        # Check offset
                        correct = abs(refoff - al.pos) < args.wiggle
                    assert trainingNm in ['Unp', 'M1', 'M2', 'Conc']
                    if trainingNm == "Unp":
                        self.training.addUnp(al, correct)
                    elif trainingNm[0] == "M":
                        self.training.addM(al, correct)
                    if mate1 is not None:
                        assert trainingNm == "Conc"
                        correct1, correct2 = correct, lastCorrect
                        if (lastAl.flags & 64) != 0:
                            correct1, correct2 = correct2, correct1
                        self.training.addPaired(mate1, mate2, abs(mate1.tlen), correct1, correct2)
                elif self.unalFh is not None:
                    # Perhaps write unaligned simulated read to file
                    self.unalFh.write("@%s\n%s\n+\n%s\n" % (nm, seq, qual))
            else:
                # Take alignment info into account
                if al.isAligned():
                    if al.rnext == "=" and al.mateMapped():
                        # Can't easily get read length for opposite mate, so
                        # just use length of this mate as a surrogate
                        fraglen = abs(al.pnext - al.pos) + len(al.seq)
                        self.dists.addFraglen(fraglen)
                    if mate1 is not None:
                        self.dists.addPair(mate1, mate2)
                    else:
                        self.dists.addRead(al)
            # Send SAM to SAM output filehandle
            if self.samOfh is not None:
                self.samOfh.write(ln)
            
            if mate1 is None:
                lastAl = al
                lastCorrect = correct
            else:
                lastAl = None
                lastCorrect = None

class AsyncWriter(threading.Thread):
    
    ''' Thread that takes items off a queue and writes them to a named file,
        which could be a FIFO '''
    
    def __init__(self, fn, q):
        threading.Thread.__init__(self)
        self.fn = fn
        self.q = q
    
    def run(self):
        ''' Write reads to the FIFO.  When None appears on the queue, we're
            done receiving reads and we close the FIFO filehandle. '''
        i = 0
        fh = open(self.fn, 'w')
        while True:
            item = self.q.get()
            if item is None:
                # None signals the final write
                self.q.task_done()
                fh.close()
                break
            fh.write(item)
            self.q.task_done()
            i += 1
        fh.close() # close FIFO filehandle

def createInput():
    ''' Return an Input object that reads all user-provided input reads '''
    assert args.fasta or args.fastq
    return Input(format="fastq" if args.fastq else "fasta", unpFns=args.U,
                 m1Fns=args.m1, m2Fns=args.m2)

class Bowtie2(object):
    
    ''' Encapsulates a Bowtie 2 process '''
    
    def __init__(self, bt2_cmd, unpaired, paired):
        # FIFO will become the input to Bowtie 2
        tmpdir = tempfile.mkdtemp()
        self.fn = os.path.join(tmpdir, 'fifo')
        os.mkfifo(self.fn)
        bt2_cmd += (" --tab6 " + self.fn)
        logging.info('Bowtie 2 command: ' + bt2_cmd)
        self.pipe = subprocess.Popen(bt2_cmd, shell=True, stdout=subprocess.PIPE, bufsize=-1)
        self.Q = Queue()
        # Takes reads from queue and writes them to FIFO
        self.W = AsyncWriter(self.fn, self.Q)
        self.W.start()
    
    def close(self):
        ''' Write None, indicating no more reads '''
        self.Q.put(None)
    
    def join(self):
        ''' Join the AyncWriters '''
        logging.debug('Joining Q')
        self.Q.join()
        logging.debug('Joining W')
        self.W.join()
        os.unlink(self.fn)
    
    def put(self, rd1, rd2=None):
        ''' Put a read on the appropriate input queue '''
        self.Q.put(Read.toTab6(rd1, rd2) + "\n")
    
    def stdout(self):
        ''' Get standard-out filehandle for Bowtie 2 process '''
        return self.pipe.stdout

def go(args, bowtieArgs):
    ''' Main driver for tandem simulator '''
    
    import tempfile
    
    if args.ref is None:
        raise RuntimeError("Must specify --ref")
    
    random.seed(args.seed)
    
    # Set up logger
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',
                        datefmt='%m/%d/%y-%H:%M:%S',
                        level=logging.DEBUG if args.verbose else logging.INFO)
    
    # Create output directory if needed
    if not os.path.isdir(args.output_directory):
        try: os.makedirs(args.output_directory)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    
    # Construct command to invoke bowtie2
    bt2_cmd = "bowtie2 "
    if args.bt2_exe is not None:
        bt2_cmd = args.bt2_exe + " "
    bt2_cmd += ' '.join(bowtieArgs)
    bt2_cmd += " --reorder --sam-no-qname-trunc -q --mapq-extra"
    
    # Timing measurements
    setupIval, al1Ival, al2Ival = None, None, None
    
    st = time.clock()
    
    # ##################################################
    # ALIGN REAL DATA
    # ##################################################
    
    bt2 = Bowtie2(bt2_cmd, args.U is not None, args.m1 is not None)
    samFn = os.path.join(args.output_directory, 'input.sam')
    inputSamFh = open(samFn, 'w') if (args.write_input_sam or args.write_all) else None
    dists = Dists()
    
    logging.info('Real-data Bowtie 2 command: "%s"' % bt2_cmd)
    
    # Create the thread that eavesdrops on output from bowtie2
    othread = Output(
        bt2.stdout(),    # SAM input filehandle
        inputSamFh,      # SAM output filehandle
        None,            # Write unaligned reads here
        None,            # Training data
        dists)           # empirical dists
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
        bt2.put(rd1, rd2)
        numReads += 1
        if numReads >= upto: break
    
    logging.debug('Closing bt2'); bt2.close();
    logging.debug('Joining bt2'); bt2.join();
    logging.debug('Joining othread'); othread.join()
    
    logging.info('Finished aligning input reads')
    if inputSamFh is not None:
        logging.info('  Input read alignments written to "%s"' % samFn)
    
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
    bt2 = Bowtie2(bt2_cmd, True, dists.hasPaired())
    samFn = os.path.join(args.output_directory, 'training.sam')
    trainingSamFh = open(samFn, 'w') if (args.write_training_sam or args.write_all) else None
    training = Training()
    
    # Create thread that eavesdrops on output from bowtie2 with simulated input
    logging.info('  Opening output-parsing thread')
    othread = Output(\
        bt2.stdout(),    # SAM input filehandle
        trainingSamFh,   # SAM output filehandle
        None,            # Write unaligned reads here
        training,        # Training data
        None)            # qualities/edits for unpaired
    othread.start()
    
    # Construct sequence and quality simulators
    logging.info('  Creating sequence simulator')
    sim = SequenceSimulator(args.ref, args.pickle_ref, idx_fn=args.ref_idx)
    simw = SimulatorWrapper(\
        sim,             # sequence simulator
        not args.m1rc,   # whether m1 is revcomped w/r/t fragment
        args.m2fw,       # whether m2 is revcomped w/r/t fragment
        dists)           # empirical distribution
    
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
            bt2.put(rdp1, rdp2)
            bt2.put(rdm1); bt2.put(rdm2)
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
            bt2.put(rd)
        if writeTrainingReads:
            trainingUnpFh.close()
    
    # Signal that we're done supplying unpaired reads
    logging.debug('Closing bt2'); bt2.close();
    logging.debug('Joining bt2'); bt2.join();
    logging.debug('Joining othread'); othread.join()
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
    
    # Writing training data
    trainingFn = os.path.join(args.output_directory, 'training.pickle')
    if not trainingFn.endswith('.gz'):
        trainingFn += '.gz'
    training.save(trainingFn)
    logging.info('Training data written to "%s"' % trainingFn)
    
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
        '--ref', metavar='path', type=str, nargs='+', help='FASTA file(s) containing reference genome sequences')
    parser.add_argument(\
        '--pickle-ref', metavar='path', type=str,
        help='Pickle FASTA input for speed, or use pickled version if it '
             'exists already.  Pickled version is stored at given path')
    parser.add_argument(\
        '--ref-idx', metavar='path', type=str, help='File containing index for FASTA file(s)')
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
        '--max-ref-bases', metavar='int', dest='max_bases', type=int, default=None,
        required=False, help='Stop reading in FASTA once we exceed this many reference nucleotides')
    
    parser.add_argument(\
        '--sam-input', metavar='path', type=str,
        help='Input SAM file to apply training data to.  Use with --training-input.')
    parser.add_argument(\
        '--bt2-exe', metavar='path', dest='bt2_exe', type=str, help='Path to Bowtie 2 exe')
    
    parser.add_argument(\
        '--sanity', action='store_const', const=True, default=False, help='Do various sanity checks')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False, help='Do unit tests')
    parser.add_argument(\
        '--profile', action='store_const', const=True, default=False, help='Print profiling info')
    parser.add_argument(\
        '--verbose', action='store_const', const=True, default=False, help='Be talkative')
    parser.add_argument(\
        '--version', action='store_const', const=True, default=False, help='Print version and quit')
    
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
        '--write-training-sam', action='store_const', const=True,
        default=False, help='Write SAM alignments for the training reads to "training.sam" in output directory')
    parser.add_argument(\
        '--write-all', action='store_const', const=True,
        default=False, help='Same as specifying all --write-* options')
    parser.add_argument(\
        '--compress-output', action='store_const', const=True,
        default=False, help='gzip all output files')
    
    argv = sys.argv
    bowtieArgs = []
    in_args = False
    for i in xrange(1, len(sys.argv)):
        if in_args:
            bowtieArgs.append(sys.argv[i])
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
        cProfile.run('go(args, bowtieArgs)')
    else:
        go(args, bowtieArgs)
