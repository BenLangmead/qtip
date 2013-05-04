#!/usr/bin/env python

"""
ts.py

A "tandem simulator," which wraps an alignment tool as it runs, eavesdrops on
the input and output, and builds a model that can be used to improve the
quality values calculated for aligned reads.

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

TODO:
- Test whether paired-end is doing a good job

"""

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
import numpy
from Queue import Queue
from sklearn.neighbors import KNeighborsClassifier
from scipy.spatial import KDTree
import collections

def quit_handler(signum,frame):
    traceback.print_stack()

signal.signal(signal.SIGQUIT,quit_handler)

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

class Read(object):
    """ Encapsulates one read """
    def __init__(self, name, seq, qual, orig=None):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.orig = orig
        assert self.repOk()
    
    @classmethod
    def fromSimulator(cls, seq, qual, refid, refoff, fw, sc, trainingNm):
        # Construct appropriate name
        rdname = "!!ts-sep!!".join(["!!ts!!", refid, "+" if fw else "-", str(refoff), str(sc), trainingNm])
        return cls(rdname, seq, qual)
    
    def __len__(self):
        """ Return number of nucleotides in read """
        return len(self.seq)
    
    def __str__(self):
        """ Return string representation """
        if self.orig is not None:
            return self.orig # original string preferred
        elif self.qual is not None:
            assert len(self.seq) == len(self.qual)
            return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
        else:
            return ">%s\n%s\n" % (self.name, self.seq)
    
    def repOk(self):
        """ Check that read is internally consistent """
        if self.qual is not None:
            assert len(self.seq) == len(self.qual)
        return True

class Alignment(object):
    """ Encapsulates an alignment record for a single aligned read """

    __asRe  = re.compile('AS:i:([-]?[0-9]+)')
    __xsRe  = re.compile('XS:i:([-]?[0-9]+)')
    __ytRe  = re.compile('YT:Z:([A-Z]+)')
    __ysRe  = re.compile('YS:i:([-]?[0-9]+)')
    __mdRe  = re.compile('MD:Z:([^\s]+)')
    __xlsRe = re.compile('Xs:i:([-]?[0-9]+)')
    __kNRe  = re.compile('YN:i:([-]?[0-9]+)')
    __knRe  = re.compile('Yn:i:([-]?[0-9]+)')
    
    def __init__(self, ln):
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
        se = Alignment.__asRe.search(self.extra)
        self.bestScore = None
        if se is not None:
            self.bestScore = int(se.group(1))
        se = Alignment.__xsRe.search(self.extra)
        self.secondBestScore = None
        if se is not None:
            self.secondBestScore = int(se.group(1))
        se = Alignment.__xlsRe.search(self.extra)
        self.thirdBestScore = None
        if se is not None:
            self.thirdBestScore = int(se.group(1))
        se = Alignment.__ytRe.search(self.extra)
        self.alType = None
        if se is not None:
            self.alType = se.group(1)
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
        """ Return true iff read aligned """
        return (self.flags & 4) == 0
    
    def orientation(self):
        """ Return orientation as + or - """
        if (self.flags & 16) != 0:
            return "-"
        else:
            return "+"
    
    def fragmentLength(self):
        """ Return fragment length """
        return abs(self.tlen)
    
    def __len__(self):
        """ Return read length """
        return len(self.seq)
    
    def repOk(self):
        """ Check alignment for internal consistency """
        assert self.alType is not None
        assert self.paired or self.fragmentLength() == 0
        assert not self.isAligned() or self.bestScore is not None
        assert self.alType in ["CP", "DP", "UP", "UU"]
        return True

class WeightedRandomGenerator(object):
    
    """ Given an ordered list of weights, generate with each call to next() an
        offset into the list of the weights with probability equal to the
        fraction of the total weight. """
    
    def __init__(self, weights):
        self.totals = []
        running_total = 0
        for w in iter(weights):
            running_total += w
            self.totals.append(running_total)
    
    def next(self):
        rnd = random.random() * self.totals[-1]
        return bisect.bisect_right(self.totals, rnd)

def cigarToList(cigar):
    """ Parse CIGAR string into a list of CIGAR operations. """
    ret = [];
    i = 0
    # CIGAR operations: MIDNSHP
    op_map = {'M':0, '=':0, 'X':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6}
    while i < len(cigar):
        run = 0
        while i < len(cigar) and cigar[i].isdigit():
            run *= 10
            run += int(cigar[i])
            i += 1
        assert i < len(cigar)
        op = cigar[i]
        i += 1
        assert op in op_map
        ret.append([op_map[op], run])
    return ret

def mdzToList(md):
    """ Parse MD:Z string into a list of operations, where 0=match,
        1=read gap, 2=mismatch. """
    i = 0;
    ret = [] # list of (op, run, str) tuples
    while i < len(md):
        if md[i].isdigit(): # stretch of matches
            run = 0
            while i < len(md) and md[i].isdigit():
                run *= 10
                run += int(md[i])
                i += 1 # skip over digit
            if run > 0:
                ret.append([0, run, ""])
        elif md[i].isalpha(): # stretch of mismatches
            mmstr = ""
            while i < len(md) and md[i].isalpha():
                mmstr += md[i]
                i += 1
            assert len(mmstr) > 0
            ret.append([1, len(mmstr), mmstr])
        elif md[i] == "^": # read gap
            i += 1 # skip over ^
            refstr = ""
            while i < len(md) and md[i].isalpha():
                refstr += md[i]
                i += 1 # skip over inserted character
            assert len(refstr) > 0
            ret.append([2, len(refstr), refstr])
        else: assert False
    return ret

def cigarMdzToStacked(seq, cgp, mdp_orig):
    """ Takes parsed cigar and parsed MD:Z and generates a stacked alignment:
        a pair of strings with gap characters inserted (possibly) and where
        characters at at the same offsets are opposite each other in the
        alignment.  Returns tuple of 2 parallel strings: read string, ref
        string. """
    mdp = mdp_orig[:]
    rds, rfs = "", ""
    mdo, rdoff = 0, 0
    for c in cgp:
        assert mdo < len(mdp), "mdo=%d" % mdo
        op, run = c
        if op == 0:   # M
            # Look for block matches and mismatches in MD:Z string
            mdrun = 0
            runleft = run
            while runleft > 0 and mdo < len(mdp):
                op_m, run_m, st_m = mdp[mdo]
                run_comb = min(runleft, run_m)
                runleft -= run_comb
                assert op_m == 0 or op_m == 1
                rds += seq[rdoff:rdoff + run_comb]
                if op_m == 0:
                    rfs += seq[rdoff:rdoff + run_comb]
                else:
                    assert len(st_m) == run_comb
                    rfs += st_m
                    for i in xrange(0, run_comb):
                        assert st_m[i] != seq[rdoff+i:rdoff+i+1]
                mdrun += run_comb
                rdoff += run_comb
                # A stretch of matches in MD:Z could span M and I sections of
                # CIGAR
                if run_comb < run_m:
                    assert op_m == 0
                    mdp[mdo][1] -= run_comb
                else:
                    mdo += 1
        elif op == 1: # I
            rds += seq[rdoff:rdoff + run]
            rfs += "-" * run
            rdoff += run
        elif op == 2: # D
            op_m, run_m, st_m = mdp[mdo]
            assert op_m == 2
            assert run == run_m
            assert len(st_m) == run
            mdo += 1
            rds += "-" * run
            rfs += st_m
        elif op == 3: # N
            # If this is really long, this probably isn't the best way to do
            # this
            rds += "-" * run
            rfs += "-" * run
        elif op == 4: # S
            rdoff += run
        elif op == 5: # H
            pass
        elif op == 6: # P
            assert False
        elif op == 7: # =
            assert False
        elif op == 8: # X
            assert False
        else: assert False
    assert mdo == len(mdp)
    return rds, rfs

class ReservoirSampler(object):
    """ Simple reservoir sampler """
    def __init__(self, k):
        self.k = k
        self.r = []
        self.n = 0
    
    def add(self, obj):
        if self.n < self.k:
            self.r.append(obj)
        else:
            j = random.randint(0, self.n)
            if j < self.k:
                self.r[j] = obj
        self.n += 1
    
    def draw(self):
        return random.choice(self.r)
    
    def __len__(self):
        return self.n
    
    def empty(self):
        return len(self) == 0

class ScoreDist(object):
    """ Capture a list of tuples, where each tuples represents the following
        traits observed in an alignment: (a) orientation, (b) quality string,
        (c) read side of stacked alignment, (d) ref side of stacked alignment,
        score.  """
    
    def __init__(self, k=10000):
        self.res = ReservoirSampler(k)
    
    def draw(self):
        assert len(self.res) > 0
        return self.res.draw()
    
    def add(self, al):
        # Extract quality string
        assert al.cigar is not None
        assert al.mdz is not None
        # Extract CIGAR, MD:Z
        cigarList = cigarToList(al.cigar)
        mdzList = mdzToList(al.mdz)
        sc = al.bestScore
        # Get stacked alignment
        rdAln, rfAln = cigarMdzToStacked(al.seq, cigarList, mdzList)
        self.res.add((al.fw, al.qual, rdAln, rfAln, sc))
    
    def __len__(self):
        return len(self.res)
    
    def empty(self):
        return len(self) == 0

class Dist(object):
    """ Basically a histogram. """
    
    def __init__(self):
        self.hist = dict()
        self.tot = 0
        self.changed = True
        self.gen = None
        self.vals = None
    
    def __str__(self):
        return str(self.hist)
    
    def __len__(self):
        return self.tot
    
    def empty(self):
        return len(self) == 0
    
    def draw(self):
        if len(self) == 0:
            raise RuntimeError("Attempt to draw from empty empirical distribution")
        if self.changed:
            self.gen = WeightedRandomGenerator(self.hist.itervalues())
            self.changed = False
        return self.hist.keys()[self.gen.next()]
    
    def add(self, key):
        self.hist[key] = self.hist.get(key, 0) + 1
        self.tot += 1
        self.changed = True

class SequenceSimulator(object):
    """ Class that, given a collection of FASTA files, samples intervals of
        specified length from the strings contained in them. """
    
    def __init__(self, fafns, pickleFn, idx_fn=None, verbose=False):
        self.__re = re.compile('[^ACGTacgt]')
        self.refs = dict()
        self.names = []
        self.lens = []
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
                                print >> sys.stderr, "Read %d FASTA bytes..." % totlen
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
        if verbose: print >>sys.stderr, "sim called..."
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
        if verbose:
            print >>sys.stderr, "...done"
        return (nm, refoff, fw, seq) # return ref id, ref offset, orientation, sequence

def mutate(rd, rdfw, scDistDraw):
    """ Given a read that already has the appropriate length (i.e.
        equal to the number of characters on the reference side of
        the alignment) 
        
        Modifies 'rd' in place. """ 
    fw, qual, rdAln, rfAln, sc = scDistDraw
    assert "N" not in rd.seq
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
    
    """ Wrapper that sends requests to the Simulator but uses information
        gathered during alignment so far to select such parameters as read
        length, concordant/discordant fragment length, etc.
        
        With each call to next(), we either simulate a paired-end read or an
        unpaired read.
        """
    
    def __init__(self, sim, m1fw, m2fw, dists):
        self.sim  = sim    # sequence simulator
        self.m1fw = m1fw   # whether mate 1 is revcomped w/r/t fragment
        self.m2fw = m2fw   # whether mate 2 is revcomped w/r/t fragment
        self.dists = dists # empirical distributions
    
    def nextUnpaired(self):
        # Simulating unpaired read
        scDraw = self.dists.scDistUnp.draw()
        _, _, _, rfAln, _ = scDraw
        rl = len(rfAln) - rfAln.count('-')
        refid, refoff, fw, seq = self.sim.sim(rl) # simulate it
        assert rl == len(seq)
        _, _, _, _, sc = scDraw
        read = Read.fromSimulator(seq, None, refid, refoff, fw, sc, "Unp")
        mutate(read, fw, scDraw) # mutate unpaired read
        assert read.qual is not None
        return read
    
    def nextPair(self):
        # Simulating paired-end read
        fl, rl1, rl2 = 0, 1, 1
        sc1Draw, sc2Draw = None, None
        while fl < rl1 or fl < rl2:
            fl = self.dists.fragDist.draw()
            sc1Draw = self.dists.scDistM1.draw()
            sc2Draw = self.dists.scDistM2.draw()
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
    """ Class that parses reads from input files and yields the reads/pairs
        produced using generators """
    
    @staticmethod
    def fastaParse(fh):
        lns = [ fh.readline().rstrip() for x in xrange(0, 2) ]
        orig = '\n'.join(lns) + '\n'
        if len(lns[0]) == 0: return None
        return Read(lns[0][1:], lns[1], None, orig)

    @staticmethod
    def fastqParse(fh):
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
        else:
            raise RuntimeError("Bad input format: '%s'" % format)
        self.unpFns = unpFns
        self.m1Fns = m1Fns
        self.m2Fns = m2Fns
    
    def __iter__(self):
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

class KNN:
    
    def __init__(self, k=1000, leafSz=20, cacheSz=4096):
        self.leafSz = leafSz
        self.k = k
        self.kdt = None
        self.tups = []
        self.tupMap = {}
        self.hits, self.misses = 0, 0
        self.ncor, self.ntot = 0, 0
        
        # Very simple initial cache: remember last query and answer
        self.prevTup, self.prevProb = None, None
        
        # Cache is a doubly-linked list
        # Link layout:     [PREV, NEXT, KEY, RESULT]
        self.root = root = [None, None, None, None]
        self.cache = cache = {}
        last = root
        for i in range(cacheSz):
            key = object()
            cache[key] = last[1] = last = [last, root, key, None]
        root[0] = last
    
    def fit(self, tups, corrects):
        """ Add a collection of training tuples, each with associated
            boolean indicating whether tuple corresponds to a correct
            or incorrect alignment. """
        self.ntot = len(tups)
        for tup, correct in zip(tups, corrects):
            correcti = 1 if correct else 0
            self.ncor += correcti
            if tup in self.tupMap:
                self.tupMap[tup][0] += correcti # update # correct
                self.tupMap[tup][1] += 1        # update total
            else:
                self.tups.append(tup)
                self.tupMap[tup] = [correcti, 1]
        assert self.ntot > 0
        self.finalize()
    
    def finalize(self):
        """ Called when all training data tuples have been added. """
        assert self.kdt is None
        self.kdt = KDTree(self.tups, leafsize=self.leafSz)
    
    def __probCorrect(self, tup):
        """ Query k nearest neighbors of test tuple (tup); that is, the
            k training tuples most "like" the test tuple.  Return the
            fraction of the neighbors that are correct. """
        
        # Fast path: test point is on top of a stack of training points
        # that is already >= k points tall
        if tup in self.tupMap and self.tupMap[tup][1] >= self.k:
            ncor, ntot =  self.tupMap[tup]
            return float(ncor) / ntot
        
        radius = 50
        bestEst = None
        ntot = 0
        wtups = [] # weighted neighbor tuples
        maxDist, minDist = 0.0, 0.0
        while radius < 1000:
            neighbors = self.kdt.query_ball_point([tup], radius, p=2.0)
            for idx in neighbors[0]:
                ntup = self.tups[idx]
                assert ntup in self.tupMap
                ntot += self.tupMap[ntup][1]
                # Calculate distance 
                dist = numpy.linalg.norm(numpy.subtract(ntup, tup))
                # Calculate a weight that decreases with increasing distance
                maxDist = max(maxDist, dist)
                minDist = min(minDist, dist)
                wtups.append((self.tupMap[ntup][0], self.tupMap[ntup][1], dist))
            if ntot >= self.k:
                break
            radius *= 2.0
        if ntot == 0:
            print >> sys.stderr, "Tuple not within 1000 units of any other tuple: %s" % str(tup)
            return float(self.ncor) / self.ntot
        ncor = sum(map(lambda x: x[0] * (maxDist - x[2]) / (maxDist - minDist), wtups))
        ntot = sum(map(lambda x: x[1] * (maxDist - x[2]) / (maxDist - minDist), wtups))
        return float(ncor) / ntot
    
    def probCorrect(self, tup):
        """ Cacheing wrapper for self.__probCorrect. """
        assert self.kdt is not None
        assert isinstance(tup, collections.Hashable)
        if tup == self.prevTup:
            return self.prevProb
        self.prevTup = tup
        cache = self.cache
        root = self.root
        link = cache.get(tup)
        if link is not None:
            # Cache hit!
            link_prev, link_next, _, result = link
            link_prev[1] = link_next
            link_next[0] = link_prev
            last = root[0]
            last[1] = root[0] = link
            link[0] = last
            link[1] = root
            self.hits += 1
            self.prevProb = result
            return result
        # Cache miss
        result = self.__probCorrect(tup)
        root[2] = tup
        root[3] = result
        oldroot = root
        root = self.root = root[1]
        root[2], oldkey = None, root[2]
        root[3], oldvalue = None, root[3]
        del cache[oldkey]
        cache[tup] = oldroot
        self.misses += 1
        self.prevProb = result
        return result

class Training(object):
    
    """ Encapsulates all training data and classifiers. """
    
    def __init__(self):
        # Training data for individual reads and mates.  Training tuples are
        # (rdlen, minValid, maxValid, bestSc, scDiff)
        self.trainUnp,    self.labUnp,    self.classUnp    = [], [], None
        self.trainM,      self.labM,      self.classM      = [], [], None
        # Training data for concordant pairs.  Training tuples are two tuples
        # like the one described above, one for each mate, plus the fragment
        # length.  The label says whether the first mate's alignment is
        # correct. 
        self.trainConc,   self.labConc,   self.classConc   = [], [], None
        self.classPair = None
        # These scale factors are set for real when we fit
        self.scaleAs, self.scaleDiff = 1.0, 1.0
    
    def __len__(self):
        """ Return number of pieces of training data added so far """
        return len(self.trainUnp) + len(self.trainM) + len(self.trainConc)
    
    def addPaired(self, al1, al2, fraglen, correct1, correct2):
        """ Add a concordant paired-end alignment to our collection of
            training data. """
        rec1 = TrainingRecord.fromAlignment(al1)
        rec2 = TrainingRecord.fromAlignment(al2)
        self.trainConc.append(rec1.toListV2() + rec2.toListV2() + [fraglen])
        self.trainConc.append(rec2.toListV2() + rec1.toListV2() + [fraglen])
        self.labConc.extend([correct1, correct2])
    
    def addUnp(self, al, correct):
        """ Add an alignment for a simulated unpaired read to our collection of
            training data. """
        rec = TrainingRecord.fromAlignment(al)
        self.trainUnp.append(rec.toListV2()); self.labUnp.append(correct)

    def addM(self, al, correct):
        """ Add an alignment for a simulated mate 1 that has been aligned in an
            unpaired fashion to our collection of training data. """
        rec = TrainingRecord.fromAlignment(al)
        self.trainM.append(rec.toListV2()); self.labM.append(correct)
    
    def probCorrect(self, al1, al2=None):
        """ Return probability that given alignment or is correct.  Arguments
            might form a concordant paired-end alignment. """
        assert al1.isAligned()
        rec = TrainingRecord.fromAlignment(al1, self.scaleAs, self.scaleDiff)
        if al2 is not None and al1.concordant:
            assert al2.isAligned()
            assert al1.mate1 and not al1.mate2
            assert al2.mate2 and not al2.mate1
            assert self.classPair is not None
            rec1, rec2 = rec, TrainingRecord.fromAlignment(al2, self.scaleAs, self.scaleDiff)
            probCor1 = self.classM.probCorrect(rec1.toTupleV1())
            probCor2 = self.classM.probCorrect(rec2.toTupleV1())
            concRec1 = ( probCor1, probCor2, abs(al1.tlen) )
            concRec2 = ( probCor2, probCor1, abs(al1.tlen) )
            newProbCor1 = self.classPair.probCorrect(concRec1)
            newProbCor2 = self.classPair.probCorrect(concRec2)
            return newProbCor1, newProbCor2
        elif al1.mate1 or al1.mate2:
            return self.classM.probCorrect(rec.toTupleV1())
        else:
            return self.classUnp.probCorrect(rec.toTupleV1())
    
    def fit(self, num_neighbors, scaleAs=1.0, scaleDiff=1.0):
        """ Train our KNN classifiers """
        self.scaleAs = scaleAs
        self.scaleDiff = scaleDiff
        
        def __adjustScale(torig):
            # Scale fields as requested
            ts = []
            for t in torig:
                rdlen, _, _, bestSc, scDiff = t
                bestSc *= scaleAs
                scDiff *= scaleDiff
                ts.append((rdlen, bestSc, scDiff))
            return ts
        
        # Create and train each classifier
        if len(self.trainUnp) > 0:
            self.classUnp = KNN(k=num_neighbors)
            ts = __adjustScale(self.trainUnp)
            self.classUnp.fit(ts, self.labUnp)
        if len(self.trainM) > 0:
            self.classM = KNN(k=num_neighbors)
            ts = __adjustScale(self.trainM)
            self.classM.fit(ts, self.labM)
        
        # Create training data for the concordant-alignment classifier.
        # This depends on M classifier already having been fit.
        if len(self.trainConc) > 0:
            pairTrain, pairLab = [], []
            for i in xrange(0, len(self.trainConc)):
                rdlenA, _, _, bestScA, scDiffA, \
                rdlenB, _, _, bestScB, scDiffB, fraglen = self.trainConc[i]
                correctA = self.labConc[i]
                bestScA *= scaleAs; scDiffA *= scaleDiff
                bestScB *= scaleAs; scDiffB *= scaleDiff
                pA = self.classM.probCorrect((rdlenA, bestScA, scDiffA))
                pB = self.classM.probCorrect((rdlenB, bestScB, scDiffB))
                pairTrain.append((pA, pB, fraglen))
                pairLab.append(correctA)
            self.classPair = KNN(k=num_neighbors)
            self.classPair.fit(pairTrain, pairLab)
    
    def save(self, fn):
        """ Save all training data to a file """
        save = (\
        self.trainUnp,    self.labUnp, \
        self.trainM,      self.labM, \
        self.trainConc,   self.labConc, \
        self.scaleAs,     self.scaleDiff )
        with open(fn, 'wb') as trainOfh:
            cPickle.dump(save, trainOfh, cPickle.HIGHEST_PROTOCOL)
    
    def load(self, fn):
        """ Load all training data from a file """
        with open(fn, 'rb') as trainIfh:
            (self.trainUnp,  self.labUnp, \
             self.trainM,    self.labM, \
             self.trainConc, self.labConc, \
             self.scaleAs,   self.scaleDiff) = cPickle.load(trainIfh)
    
    def saveTsv(self, fnPrefix):
        """ Save all training data to file.  We generate up to four files: one
            for unpaired reads, one for unpaired alignments of mate1s, one for
            unpaired alignments of mate2s, and one for concordant paired-end
            alignments. """
        import csv
        tups = [ (fnPrefix + "unp.training.tsv",  self.trainUnp,  self.labUnp),
                 (fnPrefix + "m.training.tsv",    self.trainM,    self.labM),
                 (fnPrefix + "conc.training.tsv", self.trainConc, self.labConc) ]
        for i in xrange(0, len(tups)):
            ofn, l, c = tups[i]
            if len(l) > 0:
                with open(ofn, 'w') as ofh:
                    writer = csv.writer(ofh, delimiter='\t')
                    if i == len(tups)-1:
                        # Concordant alignment training tuples
                        writer.writerows([["lenA", "minValidA", "maxValidA",
                                           "bestA", "diffA", "lenB",
                                           "minValidB", "maxValidB", "bestB",
                                           "diffB", "fraglen", "correctA"]])
                        writer.writerows([ list(x) + [int(y)] for x, y in zip(l, c) ])
                    else:
                        # Unpaired training tuples
                        writer.writerows([["len", "minValid", "maxValid",
                                           "best", "diff", "correct"]])
                        writer.writerows([ list(x) + [int(y)] for x, y in zip(l, c) ])
    
    def hits(self):
        return self.classUnp.hits  if self.classUnp  is not None else 0 + \
               self.classM.hits    if self.classM    is not None else 0 + \
               self.classConc.hits if self.classConc is not None else 0 + \
               self.classPair.hits if self.classPair is not None else 0
    
    def misses(self):
        return self.classUnp.misses  if self.classUnp  is not None else 0 + \
               self.classM.misses    if self.classM    is not None else 0 + \
               self.classConc.misses if self.classConc is not None else 0 + \
               self.classPair.misses if self.classPair is not None else 0

class TrainingRecord(object):
    """ Per-read tuple of training data including:
        1. Minimum valid score (derived by Bowtie from length)
        2. Maximum valid score (derived by Bowtie from length)
        3. Alignment score of best alignment
        4. Difference in alignment score b/t best, second-best """
    def __init__(self, rdlen, minValid, maxValid, bestSc, secbestSc,
                 scaleAs=1.0, scaleDiff=1.0):
        self.rdlen = rdlen
        self.minValid = minValid
        self.maxValid = maxValid
        self.bestSc = float(bestSc) * scaleAs
        self.secbestSc = secbestSc
        if secbestSc is not None:
            self.scDiff = float(abs(self.bestSc - self.secbestSc))
        else:
            self.scDiff = 999999.0
        self.scDiff *= scaleDiff
        assert self.repOk()
    
    @classmethod
    def fromAlignment(cls, al, scaleAs=1.0, scaleDiff=1.0):
        """ Initialize training record with respect to an alignment and a
            boolean indicating whether it is correct. """
        secbest = al.secondBestScore
        if secbest is None or al.thirdBestScore > secbest:
            secbest = al.thirdBestScore
        return cls(len(al), al.minValid, al.maxValid, al.bestScore,
                   secbest, scaleAs=scaleAs, scaleDiff=scaleDiff)
    
    @classmethod
    def v2tov1(cls, v2):
        rdlen, minValid, maxValid, bestSc, scDiff = v2
        return rdlen, bestSc, scDiff
    
    def toTupleV1(self):
        """ Return simple tuple form """
        return (self.rdlen, self.bestSc, self.scDiff)
    
    def toTupleV2(self):
        """ Return simple tuple form """
        return (self.rdlen, self.minValid, self.maxValid, self.bestSc, self.scDiff)
    
    def toListV1(self):
        """ Return simple list form """
        return [self.rdlen, self.bestSc, self.scDiff]
    
    def toListV2(self):
        """ Return simple list form """
        return [self.rdlen, self.minValid, self.maxValid, self.bestSc, self.scDiff]
    
    def repOk(self):
        """ Check for internal consistency """
        assert self.rdlen is not None
        assert self.bestSc is not None
        assert self.scDiff is not None
        return True

class Dists(object):
    
    """ Encapsulates all distributions that we train with real data.
        We collect random subsets of qualities/edits or unpaired reads
        and same for paired-end mate 1s and mate 2s.  We also collect
        a fragment length distribution. """
    
    def __init__(self):
        self.fragDist  = Dist()      # fragment length distribution
        self.scDistUnp = ScoreDist() # qualities/edits for unpaired
        self.scDistM1  = ScoreDist() # qualities/edits for mate #1s
        self.scDistM2  = ScoreDist() # qualities/edits for mate #2s
    
    def addConcordantPair(self, mate1, mate2, fraglen):
        self.fragDist.add(abs(mate1.tlen)) # remember fragment length
    
    def addRead(self, al):
        # remember quality string and edit pattern
        if   (al.flags &  64) != 0: self.scDistM1.add(al)
        elif (al.flags & 128) != 0: self.scDistM2.add(al)
        else:                       self.scDistUnp.add(al)
    
    def hasPaired(self):
        """ Return true iff at least one paired-end read has been added
            to our empirical distributions. """
        return not self.scDistM1.empty()
    
    def hasUnpaired(self):
        """ Return true iff at least one unpaired read has been added
            to our empirical distributions. """
        return not self.scDistUnp.empty()

class Output(threading.Thread):
    
    """ Encapsulates the output reader.  Reads SAM output from the aligner,
        updates empirical distributions for e.g. fragment length, and looks for
        records that correspond to simulated reads.  If a record corresponds to
        a simulated read, its correctness will be checked and a tuple
        written """
    
    def __init__(self, samIfh, samOfh, unalFh, training, dists):
        threading.Thread.__init__(self)
        self.samIfh = samIfh         # SAM records come from here
        self.samOfh = samOfh         # normal (non-simulated) SAM recs go here
        self.training = training     # training data store
        self.dists = dists           # distributions
        self.unalFh = unalFh         # write unaligned reads here
        self.scDiffs = {}            # histogram of expected-versus-observed score differences
    
    def run(self):
        lastAl = None
        for ln in self.samIfh:
            if ln[0] == '@':
                if self.samOfh is not None:
                    self.samOfh.write(ln) # header line
                continue
            al = Alignment(ln)
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
                    if trainingNm == "Unp":
                        self.training.addUnp(al, correct)
                    elif trainingNm[0] == "M":
                        self.training.addM(al, correct)
                    else:
                        assert trainingNm == "Conc"
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
                    if mate1 is not None:
                        self.dists.addConcordantPair(mate1, mate2, abs(mate1.tlen))
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
    
    """ Thread that takes items off a queue and writes them to a named file,
        which could be a FIFO """
    
    def __init__(self, fn, q):
        threading.Thread.__init__(self)
        self.fn = fn
        self.q = q
    
    def run(self):
        i = 0
        fh = open(self.fn, 'w')
        while True:
            item = self.q.get()
            if item is None:
                # None signals the final write
                self.q.task_done()
                fh.flush()
                fh.close()
                break
            fh.write(item)
            fh.flush()
            self.q.task_done()
            i += 1
        fh.close()

def createInput():
    """ Return an Input object that reads all user-provided input reads """
    if args.fastq:
        return Input(format="fastq", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2)
    elif args.fasta:
        return Input(format="fasta", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2)

class Bowtie2(object):
    
    """ Encapsulates a Bowtie 2 process """
    
    def __makeFifos(self, unpaired, paired):
        fifoArgs = ''
        fifoFns = [None, None, None] # FIFO filenames
        # Make temporary directory to store FIFOs
        tmpdir = tempfile.mkdtemp()
        # Make the appropriate FIFOs given the types of input we have
        if unpaired:
            fifoFns[0] = os.path.join(tmpdir, 'U.fifo')
            os.mkfifo(fifoFns[0])
            fifoArgs += (" -U " + fifoFns[0])
        if paired:
            fifoFns[1] = os.path.join(tmpdir, 'm1.fifo')
            os.mkfifo(fifoFns[1])
            fifoArgs += (" -1 " + fifoFns[1])
            fifoFns[2] = os.path.join(tmpdir, 'm2.fifo')
            os.mkfifo(fifoFns[2])
            fifoArgs += (" -2 " + fifoFns[2])
        return fifoArgs, fifoFns
    
    def __init__(self, bt2_cmd, unpaired, paired):
        self.Qs  = [None, None, None] # FIFO read queues
        self.Ws  = [None, None, None] # FIFO worker threads
        args, self.fns = self.__makeFifos(unpaired, paired)
        bt2_cmd += args
        
        # Open the Bowtie 2 process, which is going to want to start reading from
        # one or more of the FIFOs
        self.pipe = subprocess.Popen(bt2_cmd, shell=True, stdout=subprocess.PIPE)
        
        # For each input type (unpaired, mate1, mate2), initialize a queue and a
        # thread that takes reads from the queue and passes each along to the
        # appropriate FIFO.  It's important to have a separate thread for each FIFO
        # or we get deadlocks.
        for i in xrange(0, 3):
            if self.fns[i] is not None:
                self.Qs[i] = Queue()
                self.Ws[i] = AsyncWriter(self.fns[i], self.Qs[i])
                self.Ws[i].start()
    
    def close(self):
        # Write None to all the FIFOs to inform them we're done giving them reads
        for i in xrange(0, 3):
            if self.fns[i] is not None:
                self.Qs[i].put(None)
        # Join all the FIFOs and input handling threads, close and delete the FIFOs 
        for i in xrange(0, 3):
            if self.fns[i] is not None:
                self.Qs[i].join()
                if self.Ws[i] is not None:
                    self.Ws[i].join()
                    self.Ws[i] = None
                os.unlink(self.fns[i])
    
    def put(self, rd1, rd2=None):
        """ Put a read on the appropriate input queue """
        if rd2 is not None:
            # Write to -1/-2 filehandles
            self.Qs[1].put(str(rd1))
            self.Qs[2].put(str(rd2))
        else:
            # Write to -U filehandle
            self.Qs[0].put(str(rd1))
    
    def stdout(self):
        """ Get standard-out filehandle associated with this Bowtie2
            process """
        return self.pipe.stdout

def go():
    """ Main driver for tandem simulator """
    
    import tempfile
    
    if args.ref is None and args.training_input is None:
        raise RuntimeError("Must specify --ref")
    
    random.seed(args.seed)
    
    # When do we start generating training data?  If we train as we go,
    # we can start refining MAPQ estimates during the first pass, but we have
    # less data to direct the simulator to the most appropriate portions of the
    # training space.  If we wait until all the reads have been aligned and
    # generate it then, we need a second pass to refine any/all of the MAPQ
    # estimates but we have a lot of data to direct the simulator to the most
    # relevant portions of the training space. 
    
    # Construct command to invoke bowtie2
    bt2_cmd = "bowtie2 "
    if args.bt2_exe is not None:
        bt2_cmd = args.bt2_exe + " "
    if args.bt2_args is not None:
        bt2_cmd += args.bt2_args
    bt2_cmd += " --reorder --sam-no-qname-trunc -q --mapq-extra"
    
    # Structure encapsulating training data and KNN classifiers
    training = Training()
    
    # Timing measurements
    setupIval, al1Ival, al2Ival, fitIval, samIval = None, None, None, None, None
    
    if args.training_input is None:
        
        # The training data is not being given as input, so we have to
        # generate it with two Bowtie 2 runs.  The first run collects
        # some informations about data distributions.  The second run
        # draws from those distributions to simulate reads.  Each
        # simulated read that aligns is used as a training tuple.
        
        st = time.clock()
        
        # ##################################################
        # ALIGN REAL DATA
        # ##################################################
        
        bt2 = Bowtie2(bt2_cmd, args.U is not None, args.m1 is not None)
        samIfh = tempfile.TemporaryFile()
        
        if args.verbose:
            print >> sys.stderr, "Real-data Bowtie 2 command: '%s'" % bt2_cmd
        
        dists = Dists()
        
        # Create the thread that eavesdrops on output from bowtie2
        othread = Output(
            bt2.stdout(),    # SAM input filehandle
            samIfh,          # SAM output filehandle
            None,            # Write unaligned reads here
            None,            # Training data
            dists)           # empirical dists
        othread.start()
        
        # Stop timing setup
        setupIval = time.clock() - st
        st = time.clock()
        
        print >> sys.stderr, "Initializing threads, queues and FIFOs"
        
        # Read through all the input read files and direct all reads to the
        # appropriate queue
        upto = args.upto or sys.maxint
        numReads = 0
        for (rd1, rd2) in iter(createInput()):
            bt2.put(rd1, rd2)
            numReads += 1
            if numReads >= upto: break
        
        print >> sys.stderr, "Finished reading real reads"
    
        bt2.close()
        othread.join() # join the thread that monitors aligner output
        print >> sys.stderr, "Finished closing data FIFOs and joining output-monitoring thread"
        
        # Stop timing interval for alignment phase 1
        al1Ival = time.clock() - st
        st = time.clock()
        
        # ##################################################
        # ALIGN SIMULATED DATA
        # ##################################################
        #
        # Now we re-open Bowtie 2 using the same arguments.  This time we only give
        # it simulated reads.
        #
        
        print >> sys.stderr, "Opening new Bowtie 2"
        
        # We're potentially going to send four different types of simulated
        # reads to Bowtie 2, depending on what type of training data the read
        # pertains to.
        bt2 = Bowtie2(bt2_cmd, True, dists.hasPaired())
        
        unalFh = None
        if args.un_sim is not None:
            unalFh = open(args.un_sim, 'w')
        
        # Create the thread that eavesdrops on output from bowtie2 with simulated
        # input
        othread = Output(\
            bt2.stdout(),    # SAM input filehandle
            None,            # SAM output filehandle
            unalFh,          # Write unaligned reads here
            training,        # Training data
            None)            # qualities/edits for unpaired
        othread.start()
        
        # Construct sequence and quality simulators
        sim = SequenceSimulator(args.ref, args.pickle_ref, idx_fn=args.ref_idx)
        simw = SimulatorWrapper(\
            sim,           # sequence simulator
            not args.m1rc, # whether m1 is revcomped w/r/t fragment
            args.m2fw,     # whether m2 is revcomped w/r/t fragment
            dists)         # empirical distribution
        
        # Simulate concordant paired-end reads from empirical distributions
        if dists.hasPaired():
            for i in xrange(0, args.num_reads):
                if (i+1 % 100) == 0:
                    print >> sys.stderr, "Simulating paired-end read %d" % i
                rdp1, rdp2, rdm1, rdm2 = simw.nextPair()
                bt2.put(rdp1, rdp2)
                bt2.put(rdm1); bt2.put(rdm2)
        
        # Simulate unpaired reads from empirical distributions
        if dists.hasUnpaired():
            for i in xrange(0, args.num_reads):
                if (i+1 % 100) == 0:
                    print >> sys.stderr, "Simulating unpaired read %d" % i
                bt2.put(simw.nextUnpaired())
        
        print >> sys.stderr, "Finished simulating reads"
        bt2.close()
        print >> sys.stderr, "Closed FIFOs"
        othread.join() # join the thread that monitors aligner output
        if unalFh is not None:
            unalFh.close()
        print >> sys.stderr, "Finished closing simulation FIFOs and joining output-monitoring thread"
        
        # Stop timing interval for alignment phase 2
        al2Ival = time.clock() - st
        
        print >> sys.stderr, "Score difference (expected - actual) histogram:"
        for k, v in sorted(othread.scDiffs.iteritems()):
            print >>sys.stderr, "  %d: %d" % (k, v)
        
        if args.save_training is not None:
            training.save(args.save_training)
    else:
        st = time.clock()
        
        # Training data is being given to us
        if args.sam_input is None:
            raise RuntimeError("Must specify --sam-input along with --training-input")
        training.load(args.training_input)
        print >> sys.stderr, "Read %d training tuples from '%s'" % (len(training), args.training_input)
        samIfh = open(args.sam_input, 'r')
        
        # Stop timing interval for setup phase
        setupIval = time.clock() - st
    
    if args.save_training_tsv: training.saveTsv(args.save_training_tsv + ".")
    
    # Build KNN classifiers
    st = time.clock()
    training.fit(args.num_neighbors, scaleAs=args.as_scale, scaleDiff=args.diff_scale)
    print >> sys.stderr, "Finished fitting KNN classifiers on %d training tuples" % len(training)
    fitIval = time.clock() - st
    st = time.clock()
    
    mapqDiff = 0.0
    nrecs, npair, nunp = 0, 0, 0
    samTsvOfh = None
    
    exFields = set(["XQ:f"])
    if args.save_sam_tsv:
        # Figure out what all possible extra fields are so that we can write a
        # tsv version of the SAM output with a fixed number of columns per line 
        samIfh.seek(0)
        for samrec in samIfh:
            if samrec[0] == '@':
                continue
            toks = string.split(samrec, '\t')
            for tok in toks[11:]:
                idx1 = tok.find(':')
                assert idx1 != -1
                idx2 = tok.find(':', idx1+2)
                assert idx2 != -1
                exFields.add(tok[:idx2])
        samTsvOfh = open(args.save_sam_tsv, 'w')
        samTsvOfh.write("\t".join([
            "qname", "flag", "rname", "pos", "mapq", "cigar", "rnext",
            "pnext", "tlen", "seq", "qual"] + list(exFields)))
        samTsvOfh.write("\n")
        print >>sys.stderr, "Extra fields: " + str(exFields)
    
    with open(args.S, 'w') as samOfh: # open output file
        
        def emitNewSam(samrec, al, probCorrect):
            probIncorrect = 1.0 - probCorrect # convert to probability incorrect
            mapq = args.max_mapq
            if probIncorrect > 0.0:
                mapq = min(-10.0 * math.log10(probIncorrect), args.max_mapq)
            samrec = samrec.rstrip()
            xqIdx = samrec.find("\tXQ:f:")
            if xqIdx != -1:
                samrec = samrec[:xqIdx]
            samOfh.write("\t".join([samrec, "XQ:f:" + str(mapq)]) + "\n")
            if samTsvOfh:
                samToks = string.split(samrec, "\t")
                samTsvOfh.write("\t".join(samToks[:11]))
                for field in exFields:
                    samTsvOfh.write("\t")
                    if field == "XQ:f":
                        samTsvOfh.write(str(mapq))
                    else:
                        samTsvOfh.write(al.exFields[field] if field in al.exFields else "NA")
                samTsvOfh.write("\n")
            return mapq - float(al.mapq)
        
        samIfh.seek(0)                # rewind temporary SAM file
        concMap = {}
        for samrec in samIfh:         # read each record
            if samrec[0] == '@':      # pass headers straight through
                samOfh.write(samrec)
                continue
            al = Alignment(samrec)    # parse alignment
            nrecs += 1
            if al.isAligned():
                # Is this alignment part of a concordantly-aligned paired-end read?
                probCorrect = None
                if al.concordant:
                    if al.name in concMap:
                        al1, samrec1 = al, samrec
                        al2, samrec2 = concMap[al.name]
                        npair += 2
                        del concMap[al.name]
                        if al2.mate1:
                            al1, al2 = al2, al1
                            samrec1, samrec2 = samrec2, samrec1
                        probCorrect1, probCorrect2 = training.probCorrect(al1, al2)
                        mapqDiff += emitNewSam(samrec1.rstrip(), al1, probCorrect1)
                        mapqDiff += emitNewSam(samrec2.rstrip(), al2, probCorrect2)
                    else:
                        concMap[al.name] = (al, samrec)
                else:
                    nunp += 1
                    probCorrect = training.probCorrect(al)
                    mapqDiff += emitNewSam(samrec.rstrip(), al, probCorrect)
            else:
                samOfh.write(samrec)
                if samTsvOfh:
                    samToks = string.split(samrec, "\t")
                    samTsvOfh.write("\t".join(samToks[:11]))
                    for field in exFields:
                        samTsvOfh.write("\t")
                        samTsvOfh.write(al.exFields[field] if field in al.exFields else "NA")
                    samTsvOfh.write("\n")
        assert len(concMap) == 0
    
    if samTsvOfh: samTsvOfh.close()
    
    samIval = time.clock() - st
    
    print >> sys.stderr, "Finished writing final SAM output (%d records) to '%s'" % (nrecs, args.S)
    print >> sys.stderr, "  concordant-pair SAM records: %d" % npair
    print >> sys.stderr, "  non-concordant-pair SAM records: %d" % nunp
    print >> sys.stderr, "  total mapping quality difference (new - old): %0.3f" % mapqDiff
    print >> sys.stderr, "  total KNN hits=%d, misses=%d" % (training.hits(), training.misses()) 
    
    if setupIval is not None:
        print >> sys.stderr, "Setup running time: %0.3f secs" % setupIval 
    if al1Ival is not None:
        print >> sys.stderr, "Alignment (real reads): %0.3f secs" % al1Ival 
    if al2Ival is not None:
        print >> sys.stderr, "Alignment (simulated reads): %0.3f secs" % al2Ival 
    if fitIval is not None:
        print >> sys.stderr, "Fitting: %0.3f secs" % fitIval
    if samIval is not None:
        print >> sys.stderr, "MAPQ calculation / SAM rewriting: %0.3f secs" % samIval
    
    # Print timing results
    
    # Close temporary SAM output file; it will be deleted immediately
    samIfh.close()

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
        '--S', metavar='path', type=str, required=True, help='Write SAM output here')
    parser.add_argument(\
        '--m1', metavar='path', type=str, nargs='+', help='Mate 1 files')
    parser.add_argument(\
        '--m2', metavar='path', type=str, nargs='+', help='Mate 2 files')
    parser.add_argument(\
        '--m1rc', action='store_const', const=True, default=False,
        help='Set if mate 1 is reverse-complemented w/r/t the fragment')
    parser.add_argument(\
        '--m2fw', action='store_const', const=True, default=False,
        help='Set if mate 2 is oriented forward w/r/t the fragment')
    parser.add_argument(\
        '--fasta', action='store_const', const=True, default=False, help='Reads are FASTA')
    parser.add_argument(\
        '--fastq', action='store_const', const=True, default=True, help='Reads are FASTQ')
    parser.add_argument(\
        '--max-mapq', metavar='float', type=float, default=100.0,
        required=False, help='Maximum MAPQ possible for an alignment')
    parser.add_argument(\
        '--seed', metavar='int', type=int, default=99099,
        required=False, help='Integer to initialize pseudo-random generator')
    parser.add_argument(\
        '--num-reads', metavar='int', type=int, default=500000,
        required=False, help='Number of reads to simulate')
    parser.add_argument(\
        '--num-neighbors', metavar='int', type=int, default=1000,
        required=False, help='Number of neighbors to use for k-nearest-neighbors')
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
        '--training-input', metavar='path', type=str,
        help='Training data to use to predict new MAPQs for a SAM file.  Use with --sam-input.')
    parser.add_argument(\
        '--save-training', metavar='path', type=str,
        help='Save training data to file.')
    parser.add_argument(\
        '--save-training-tsv', metavar='path-prefix', type=str,
        help='Save training data as a set of TSV files with given prefix.  '
             'File-specific suffixes with .tsv at the end will be added.')
    parser.add_argument(\
        '--save-sam-tsv', metavar='path', type=str,
        help='Save SAM output as a TSV file.')
    parser.add_argument(\
        '--sam-input', metavar='path', type=str,
        help='Input SAM file to apply training data to.  Use with --training-input.')
    parser.add_argument(\
        '--un-sim', metavar='path', type=str,
        help='Write unaligned simulated reads to this file.')
    parser.add_argument(\
        '--bt2-exe', metavar='path', dest='bt2_exe', type=str,
        help='Path to Bowtie 2 exe')
    parser.add_argument(\
        '--bt2-args', metavar='args', dest='bt2_args', type=str,
        help='Arguments to pass to Bowtie 2 (besides input an output)')
    parser.add_argument(\
        '--distance-weight', action='store_const', const=True, default=False,
        help='Do distance weighting when doing KNN')
    parser.add_argument(\
        '--scale', metavar='type', type=str, default="unit",
        help='"unit" for min/max scaling or "z" for standard units')
    parser.add_argument(\
        '--as-scale', metavar='float', type=float, default=1.0,
        help='Multiply AS:i scores by this before making them into training records')
    parser.add_argument(\
        '--diff-scale', metavar='float', type=float, default=3.0,
        help='Multiply AS:i - XS:i scores by this before making them into training records')
    parser.add_argument(\
        '--sanity', dest='sanity', action='store_const', const=True, default=False,
        help='Do various sanity checks')
    parser.add_argument(\
        '--test', dest='test', action='store_const', const=True, default=False,
        help='Do unit tests')
    parser.add_argument(\
        '--profile', action='store_const', const=True,
        default=False, help='Print profiling info')
    parser.add_argument(\
        '--verbose', dest='verbose', action='store_const', const=True,
        default=False, help='Be talkative')
    parser.add_argument(\
        '--version', dest='version', action='store_const', const=True,
        default=False, help='Print version and quit')
    
    args = parser.parse_args()

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
        cProfile.run('go()')
    else:
        go()
