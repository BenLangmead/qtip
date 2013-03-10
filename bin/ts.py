#!/usr/bin/env python

"""
ts.py

A "tandem simulator," which wraps an alignment tool as it runs, eavesdrops on
the input and output, and builds a model that can be used to improve the
quality values calculated for aligned reads.

Right now we have to load the reference FASTA separately, and we need the user
to tell us about the scoring params so we can simulate reads with known best
score.

Things we learn from reads
==========================

- Read length distribution
- Quality values

Things we learn from alignments
===============================

- Alignment type (aligned, unaligned, concordant, discordant)
- Fragment length distribution
- Number and placement of mismatches and gaps

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
import pickle
from Queue import Queue
from sklearn.neighbors import KNeighborsClassifier

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
    def fromSimulator(cls, seq, qual, refid, refoff, fw, sc):
        # Construct appropriate name
        rdname = "!!ts-sep!!".join(["!!ts!!", refid, "+" if fw else "-", str(refoff), str(sc)])
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

    __asRe = re.compile('AS:i:([-]?[0-9]+)')
    __xsRe = re.compile('XS:i:([-]?[0-9]+)')
    __ytRe = re.compile('YT:Z:([A-Z]+)')
    __ysRe = re.compile('YS:i:([-]?[0-9]+)')
    __mdRe = re.compile('MD:Z:([^\s]+)')
    
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
        se = Alignment.__asRe.search(self.extra)
        self.bestScore = None
        if se is not None:
            self.bestScore = int(se.group(1))
        se = Alignment.__xsRe.search(self.extra)
        self.secondBestScore = None
        if se is not None:
            self.secondBestScore = int(se.group(1))
        se = Alignment.__ytRe.search(self.extra)
        self.alType = None
        if se is not None:
            self.alType = se.group(1)
        se = Alignment.__ysRe.search(self.extra)
        self.mateBest = None
        if se is not None:
            self.mateBest = int(se.group(1))
        self.mdz = None
        se = Alignment.__mdRe.search(self.extra)
        if se is not None:
            self.mdz = se.group(1)
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

class AlignmentPair(object):
    """ Encapsulates a pair of alignments for two ends of a paired-end read """
    def __init__(self, al1, al2):
        self.al1, self.al2 = al1, al2
    
    def fraglen(self):
        return abs(self.al1.tlen)

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

class ScoreDist(object):
    """ Capture a map from scores observed to patterns of mismatches
        and gaps observed.  We do this instead of trying to generate
        patterns of gaps and mismatches to achieve a desired score,
        because the latter has many issues.  Chief amongst them:
        1. Hard to  """
    
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

class Dist(object):
    """ Capture an empirical distribution.  Basically a histogram. """
    
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

class Scoring(object):
    """ Parameters governing how to score sequence differences between
        a read and the reference genome. """ 
    
    def __init__(self, ma, mmp, np, rdg, rfg):
        self.ma = ma             # match bonus: integer
        self.mmp = mmp           # mismatch penalty: (min, max) integer tuple
        self.np = np             # N penalty: integer
        self.rdg = rdg           # affine read gap penalty: (a, b)
        self.rfg = rfg           # affine reference gap penalty: (a, b)
    
    def rfgPenalty(self, length=1):
        return self.rfg[0] + length * self.rfg[1]
    
    def rdgPenalty(self, length=1):
        return self.rdg[0] + length * self.rdg[1]

class SimpleFunc(object):
    """ One-argument function with constant term and term that is proportional
        to some function of the argument.  Function can be f(x) = 0 (constant),
        f(x) = x (linear), f(x) = ln(x) (log), or f(x) = sqrt(x) (sqrt) """
    
    def __init__(self, type="const", I=None, X=None, C=0.0, L=0.0):
        self.type = type
        self.I = I
        self.X = X
        self.C = C
        self.L = L
        if I is None: I = float('-inf')
        if X is None: X = float('inf')
    
    def f(self, x):
        if self.type == "const":
            x = 0.0
        elif self.type == "sqrt":
            x = math.sqrt(x)
        elif self.type == "log":
            x = math.log(x)
        elif self.type != "linear":
            raise RuntimeError("Bad simple function type: '%s'" % self.type)
        return min(X, max(I, self.C + self.L * x))

class QualitySimulator(object):
    """ Class that, given read lengths and a distribution of observed quality
        values, returns a string of quality values of the desired length. """
    
    def __init__(self, qdistUnp, qdistM1, qdistM2):
        self.qdistUnp = qdistUnp
        self.qdistM1 = qdistM1
        self.qdistM2 = qdistM2
    
    def sim(self, ln):
        return 'I' * ln

class FastaSeqIdx(object):
    
    ''' Object encapsulating a FASTA file that has been indexed by sequence.
        That is, the index makes it easy to seek directly to the beginning of
        a specific sequence (but not to specific points within it). '''
    
    def __init__(self, fa_fns, truncate=True, idx_fn=None):
        self.fa_fns = fa_fns  # FASTA filenames to open/index
        self.idx_fn = idx_fn  # filename where index should be stored/retrieved
        self.idx = {}         # The index!!
        self.idxIval = 10 * 1024
        self.truncate = truncate
        self.fh = None
        self._lastc = ''
        if idx_fn is not None and os.path.exists(idx_fn):
            self.loadIdx(idx_fn)
        else:
            self.buildIdx()
        if idx_fn is not None and not os.path.exists(idx_fn):
            self.saveIdx(idx_fn)
    
    def __del__(self):
        if self.fh is not None: self.fh.close()
    
    def close(self):
        if self.fh is not None: self.fh.close()
    
    def nextc(self, n=1, skip=False):
        ''' Return a string consisting of the next n non-whitespace characters.
            The read is with respect to the current file cursor. '''
        self._lastc = ''
        if n == 1:
            c = self.fh.read(1)
            while c.isspace():
                c = self.fh.read(1)
            self._lastc = c
            return c # Empty string if EOF
        else:
            s = []
            slen = 0
            lastc = None
            while slen < n:
                c = self.fh.read(1)
                while c.isspace():
                    c = self.fh.read(1)
                if not skip:
                    s.append(c)
                lastc = c
                slen += 1
            if slen > 0 and lastc is not None:
                self._lastc = lastc
            return ''.join(s)
    
    def get(self, ref, off=0, ln=None):
        ''' Get the string at 0-based offset off, with length ln '''
        if off == 0 and ln is None:
            return self.getWhole(ref)
        if ln is None: ln = self.idx[ref][1]
        self.scanTo(ref, off)
        return self.nextc(n=ln, skip=False)
    
    def getpad(self, ref, off=0, ln=None, pad='x'):
        ''' Get the string at 0-based offset off, with length ln.  If this
            string falls off either end of the reference, fill those slots with
            the "pad" string. '''
        if off == 0 and ln is None:
            return self.getWhole(ref)
        if ln is None: ln = self.idx[ref][1]
        lf, rt = '', ''
        if off < 0: lf = pad * -off
        if off+ln > self.idx[ref][1]: rt = pad * (off+ln - self.idx[ref][1])
        self.scanTo(ref, off)
        return lf + self.nextc(n=ln, skip=False) + rt
    
    def getWhole(self, ref):
        """ Return an entire sequence from the FASTA file. """ 
        self.scanTo(ref, 0)
        lns = []
        while True:
            ln = self.fh.readline().rstrip()
            if len(ln) == 0 or ln[0] == '>': break
            lns.append(ln)
        return "".join(lns)
    
    def peek(self):
        ''' Take a peek at the next non-whitespace character and return it. '''
        off = self.fh.tell()
        c = self.fh.read(1)
        while c.isspace():
            c = self.fh.read(1)
        self.fh.seek(off)
        return c
    
    def lastc(self):
        ''' Return last character read '''
        return self._lastc
    
    def scanTo(self, ref, off):
        ''' Scan to the beginning of the FASTA sequence of the given name '''
        if self.fh is not None:
            self.fh.close()
            self.fh = None
        if self.truncate:
            idx = ref.find(' ')
            # Truncate name at first whitespace
            if idx >= 0:
                ref = ref[:idx]
        # Reference should be in the index
        assert ref in self.idx
        fn, seekoff = self.idx[ref]
        offlist = self.idx[ref][3]
        i = bisect.bisect_right(offlist, (off, sys.maxint))
        assert i <= len(offlist)
        self.fh = open(fn, 'rb')
        if i == 0:
            self.fh.seek(self.idx[ref][2], 0)
            self.nextc(off, skip=True)
        else:
            diff = offlist[i][0] - off
            self.fh.seek(offlist[i][1], 0)
            self.nextc(diff, skip=True)
    
    def loadIdx(self, idx_fn):
        ''' Load an already-calculated index from a file. '''
        with open(self.idx_fn, 'rb') as pkl_file:
            self.idx = pickle.load(pkl_file)
    
    def saveIdx(self, idx_fn):
        ''' Save an index to a file. '''
        with open(self.idx_fn, 'wb') as pkl_file:
            pickle.dump(self.idx, pkl_file)
    
    def length(self, nm):
        assert nm in self.idx
        return self.idx[nm][1]
    
    def buildIdx(self):
        ''' Build an index dictionary over for all the reference sequences in
            all the FASTA files. '''
        for fn in self.fa_fns:
            with open(fn, 'rb') as fh:
                cur = None
                while 1:
                    line = fh.readline()
                    if not line: break
                    line = line.rstrip()
                    off = fh.tell()
                    if len(line) == 0:
                        continue
                    if line[0] == '>': # Header
                        refid = line[1:]
                        if self.truncate:
                            idx = refid.find(' ')
                            if idx >= 0:
                                refid = refid[:idx]
                        assert refid not in self.idx, "Already saw refid '%s' from line '%s'" % (refid, line)
                        # Index record is:
                        # 1. Filename
                        # 2. Length of sequence
                        # 3. File offset
                        # 4. List of internal offsets
                        self.idx[refid] = [fn, 0, off, []]
                        cur = refid
                    elif cur is not None:
                        ln = self.idx[cur][1]
                        ln_old, ln_new = ln, ln + len(line)
                        if (ln_old % self.idxIval) != (ln_new % self.idxIval):
                            self.idx[cur][3].append((ln_new, off))
                        # See if we want to add an offset to the index
                        self.idx[cur][1] += len(line)

class SequenceSimulator(object):
    """ Class that, given a collection of FASTA files, samples intervals of
        specified length from the strings contained in them. """
    
    def __init__(self, fafns, idx_fn=None, verbose=False):
        self.refs = dict()
        self.names = []
        self.lens = []
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
                if line.startswith('>'):
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
        self.rnd = WeightedRandomGenerator(self.lens)
        self.__re = re.compile('[^ACGTacgt]')
    
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
        while self.__re.match(seq):
            refoff = random.randint(0, self.lens[refi] - ln) # pick new offset
            seq = self.refs[nm][refoff:refoff+ln] # extract substring again
        seq = seq.upper()
        if random.random() > 0.5: # possibly reverse-complement
            fw = False
            seq = revcomp(seq) # reverse complement
        if verbose:
            print >>sys.stderr, "...done"
        return (nm, refoff, fw, seq) # return ref id, ref offset, orientation, sequence

def addTo(total, summands=[1], memo = {}):
    """ Given a list of integers, perhaps with some integers repeated,
        calculate all ways to add those numbers to get the given total.
        Uses top-down dynamic programming.
        
        Problems with this idea:
        1. It doesn't necessarily generate 'realistic' combinations of gaps
           and mismatches
        2. The lower penalties are going to tend to be overrepresented.  For
           instance, what about N getting a penalty of -1?  It will show up in
           most combinations produced.
        3. Could take far longer if summands are mix of positive and negative
        
        What we might do instead is capture patterns of gaps and mismatches
        that occur in alignments.  We can't trust these 100% because they've
        been through the filter of the aligner.  So we might "smooth" them or
        otherwise try to include a bunch of "nearby" score combinations as
        well. """ 
    def addToHelper(total, totalSofar):
        solns = set()
        left = total - totalSofar
        if left in memo:
            return memo[left]
        for i in xrange(0, len(summands)):
            if totalSofar + summands[i] < total:
                newSolns = addToHelper(total, totalSofar + summands[i])
                for soln in newSolns:
                    solnList = list(soln)
                    solnList[i] += 1
                    solns.add(tuple(solnList))
            elif totalSofar + summands[i] == total:
                soln = [0] * len(summands)
                soln[i] = 1
                solns.add(tuple(soln))
        memo[left] = solns
        return solns
    return addToHelper(total, 0)

class MyMutableString(object):
    """ A string supporting efficient insertions and deletions """
    
    def __init__(self, s=None):
        self.slist = []
        if s is not None:
            self.append(s)

    def __str__(self):
        return ''.join(self.slist)
    
    def append(self, s):
        self.slist.extend(list(s))
    
    def delete(self, i, j=None):
        if j is None:
            j = i + 1
        # Set the appropriate elements to the empty string
        for k in xrange(i, j):
            self.slist[k] = ""
    
    def set(self, i, s):
        assert i < len(self.slist)
        self.slist[i] = s
    
    def get(self, i):
        return self.slist[i]
    
def mutate(rd, rdfw, scDistDraw):
    """ Given a read that already has the appropriate length (i.e.
        equal to the number of characters on the reference side of
        the alignment) 
        
        Modifies 'rd' in place. """ 
    fw, qual, rdAln, rfAln, sc = scDistDraw
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
            assert sc != 0 # score can't be perfect
            rfi += 1
        elif rfAln[i] == '-':
            assert sc != 0 # score can't be perfect
            seq.append(rdAln[i])
            rdi += 1
        elif rdAln[i] != rfAln[i] and rdAln[i] == 'N':
            assert sc != 0 # score can't be perfect
            seq.append('N')
            rfi += 1; rdi += 1
        elif rdAln[i] != rfAln[i]:
            assert sc != 0 # score can't be perfect
            assert rfi < len(rd.seq)
            oldc = rd.seq[rfi].upper()
            cs = ['A', 'C', 'G', 'T']
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
    
    """ Wrapper that sends requests to the Simualtor but uses information
        gathered during alignment so far to select such parameters as read
        length, concordant/discordant fragment length, etc. """
    
    def __init__(self, sim, qsim, tyd, scd, scd1, scd2, fld):
        self.sim  = sim  # sequence simulator
        self.qsim = qsim # quality simulator
        self.tyd  = tyd  # type distribution (UU/CP/DP/UP)
        self.scd  = scd  # score distribution for unpaired reads
        self.scd1 = scd1 # score distribution for mate 1
        self.scd2 = scd2 # score distribution for mate 2
        self.fld  = fld  # fragment length distribution
    
    def next(self):
        """ Simulate the next read/pair and associated quality values.  Return
            the simulated read along with information about where it
            originated. """
        ty = self.tyd.draw()
        if ty[1] == 'U':
            # Simulating unpaired read
            scDraw = self.scd.draw()
            _, _, _, rfAln, _ = scDraw
            rl = len(rfAln) - rfAln.count('-')
            refid, refoff, fw, seq = self.sim.sim(rl) # simulate it
            assert rl == len(seq)
            qual = self.qsim.sim(rl)
            _, _, _, _, sc = scDraw
            read = Read.fromSimulator(seq, qual, refid, refoff, fw, sc)
            mutate(read, fw, scDraw) # mutate unpaired read
            return read, None
        else:
            # Simulating paired-end read
            fl = self.fld.draw() # draw a fragment length
            sc1Draw, sc2Draw = self.scd1.draw(), self.scd2.draw()
            _, _, _, rfAln1, _ = sc1Draw
            _, _, _, rfAln2, _ = sc2Draw
            rl1 = len(rfAln1) - rfAln1.count('-')
            rl2 = len(rfAln2) - rfAln2.count('-')
            refid, refoff, fw, seq = self.sim.sim(fl) # simulate fragment
            assert len(seq) == fl
            # get mates from fragment
            seq1 = seq[:rl1]
            seq2 = seq[-rl2:]
            # TODO: possibly reverse-comp according to paired-end parameters
            if not m1fw: seq1 = revcomp(seq1)
            if not m2fw: seq2 = revcomp(seq2)
            refoff1, refoff2 = refoff, refoff
            if fw: refoff2 = refoff + fl - rl2
            else:  refoff1 = refoff + fl - rl1
            qual1, qual2 = self.qsim.sim(len(seq1)), self.qsim.sim(len(seq2))
            _, _, _, _, sc1 = sc1Draw
            _, _, _, _, sc2 = sc2Draw
            rd1 = Read.fromSimulator(seq1, qual1, refid, refoff1, fw, sc1)
            rd2 = Read.fromSimulator(seq2, qual2, refid, refoff2, fw, sc2)
            mutate(rd1, fw, sc1Draw) # mutate mate 1
            mutate(rd2, fw, sc2Draw) # mutate mate 2
            return rd1, rd2

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

class InputWrapper(object):
    """ Wraps the input reader so that we can eavesdrop on the reads and use
        them to build empirical distributions.  Also allows us to inject
        simulated reads now and again. """
    def __init__(self, inp, rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2):
        self.inp = inp
        self.rddistUnp, self.rddistM1, self.rddistM2 = rddistUnp, rddistM1, rddistM2
        self.qdistUnp, self.qdistM1, self.qdistM2 = qdistUnp, qdistM1, qdistM2
    
    def __iter__(self):
        for (rd1, rd2) in self.inp:
            assert rd1 is not None
            # ignore simulated reads
            if not rd1.name.startswith('!!ts!!'):
                if rd2 is not None:
                    # Paired-end case
                    # add read lengths to empirical distributions
                    self.rddist1.add(len(rd1))
                    self.rddist2.add(len(rd2))
                    # add quality values to empirical distributions
                    for q in map(lambda x: ord(x)-33, rd1.qual): self.qdist1.add(q)
                    for q in map(lambda x: ord(x)-33, rd2.qual): self.qdist2.add(q)
                else:
                    # Unpaired case
                    # add read length to empirical distribution
                    self.rddistUnp.add(len(rd1))
                    # add quality values to empirical distribution
                    for q in map(lambda x: ord(x)-33, rd1.qual): self.qdistUnp.add(q)
            yield (rd1, rd2)

class Training(object):
    
    """ Encapsulates all the training data and all the classifiers. """
    def __init__(self):
        # These tables are for collecting training data pertaining to
        # individual reads and mates
        self.trainUnp,    self.labUnp,    self.classUnp    = [], [], None
        self.trainM1Disc, self.labM1Disc, self.classM1Disc = [], [], None
        self.trainM2Disc, self.labM2Disc, self.classM2Disc = [], [], None
        self.trainM1Conc, self.labM1Conc, self.classM1Conc = [], [], None
        self.trainM2Conc, self.labM2Conc, self.classM2Conc = [], [], None
        # Following tables are for a second layer of training for
        # concordantly-aligned pairs
        self.trainConc,   self.labConcM1,   self.labConcM2 = [], [], []
        self.classConcM1, self.classConcM2 = [], []
        self.trainConcFraglen = []
        # These scale factors are set for real when we fit
        self.scaleAs, self.scaleDiff = 1.0, 1.0
    
    def __len__(self):
        """ Return number of pieces of training data added so far """
        return len(self.trainUnp) + len(self.trainM1Disc) + len(self.trainM2Disc) + \
            len(self.trainM1Conc) + len(self.trainM2Conc)
    
    def add(self, al, correct):
        """ Add an alignment for a simulated read to our collection of
            training data. """
        rec = TrainingRecord.fromAlignment(al)
        if al.concordant:
            if al.mate1:
                self.trainM1Conc.append(rec.toList())
                self.labM1Conc.append(correct)
            else:
                self.trainM2Conc.append(rec.toList())
                self.labM2Conc.append(correct)
        elif al.discordant or al.unpPair:
            if al.mate1:
                self.trainM1Disc.append(rec.toList())
                self.labM1Disc.append(correct)
            else:
                self.trainM2Disc.append(rec.toList())
                self.labM2Disc.append(correct)
        else:
            self.trainUnp.append(rec.toList())
            self.labUnp.append(correct)
    
    def probCorrect(self, al1, al2=None):
        """ Return probability that given alignment is """
        assert al1.isAligned()
        rec = TrainingRecord.fromAlignment(al1, self.scaleAs, self.scaleDiff)
        if al2 is not None:
            assert al2.isAligned()
            assert al1.concordant
            assert al2.concordant
            assert al1.mate1 and not al1.mate2
            assert al2.mate2 and not al2.mate1
            rec1, rec2 = rec, TrainingRecord.fromAlignment(al2, self.scaleAs, self.scaleDiff)
            probCor1 = self.classM1Conc.predict_proba([rec1.toList()])[0][-1]
            probCor2 = self.classM2Conc.predict_proba([rec2.toList()])[0][-1]
            concRec = [ probCor1, probCor2, abs(al1.tlen) ]
            newProbCor1 = self.classConcM1.predict_proba(concRec)[0][-1]
            newProbCor2 = self.classConcM2.predict_proba(concRec)[0][-1]
            return newProbCor1, newProbCor2
        elif al1.discordant or al1.unpPair:
            classifier = self.classM1Disc if al1.mate1 else self.classM2Disc
            return classifier.predict_proba([rec.toList()])[0][-1]
        else:
            return self.classUnp.predict_proba([rec.toList()])[0][-1]
    
    def fit(self, num_neighbors, weights, scaleAs=1.0, scaleDiff=1.0):
        """ Train our KNN classifiers """
        self.scaleAs = scaleAs
        self.scaleDiff = scaleDiff
        
        def __adjustScale(t):
            # Scale fields as requested
            if scaleAs != 1.0 or scaleDiff != 1.0:
                for i in xrange(0, len(t)):
                    t[i][1] *= scaleAs
                    t[i][2] *= scaleDiff
        
        # Create and train each classifier
        if len(self.trainUnp) > 0:
            self.classUnp = KNeighborsClassifier(\
                n_neighbors=num_neighbors, warn_on_equidistant=False, weights=weights)
            __adjustScale(self.trainUnp)
            self.classUnp.fit(self.trainUnp, self.labUnp)
        if len(self.trainM1Disc) > 0:
            self.classM1Disc = KNeighborsClassifier(\
                n_neighbors=num_neighbors, warn_on_equidistant=False, weights=weights)
            __adjustScale(self.trainM1Disc)
            self.classM1Disc.fit(self.trainM1Disc, self.labM1Disc)
        if len(self.trainM2Disc) > 0:
            self.classM2Disc = KNeighborsClassifier(\
                n_neighbors=num_neighbors, warn_on_equidistant=False, weights=weights)
            __adjustScale(self.trainM2Disc)
            self.classM2Disc.fit(self.trainM2Disc, self.labM2Disc)
        if len(self.trainM1Conc) > 0:
            self.classM1Conc = KNeighborsClassifier(\
                n_neighbors=num_neighbors, warn_on_equidistant=False, weights=weights)
            __adjustScale(self.trainM1Conc)
            self.classM1Conc.fit(self.trainM1Conc, self.labM1Conc)
        if len(self.trainM2Conc) > 0:
            self.classM2Conc = KNeighborsClassifier(\
                n_neighbors=num_neighbors, warn_on_equidistant=False, weights=weights)
            __adjustScale(self.trainM2Conc)
            self.classM2Conc.fit(self.trainM2Conc, self.labM2Conc)
        
        assert len(self.trainM1Conc) == len(self.trainM2Conc)
        assert len(self.trainM1Conc) == len(self.trainConcFraglen)
        assert len(self.trainM1Disc) == len(self.trainM2Disc)
        
        # Create training data for the concordant-alignment classifier.
        # This depends on M1Conc and M2Conc classifiers already having
        # been fit.
        if len(self.trainM1Conc) > 0:
            for i in xrange(0, len(self.trainM1Conc)):
                m1c, m2c = self.trainM1Conc[i], self.trainM2Conc[i]
                p1 = self.classM1Conc.predict_proba([m1c.toList()])[0][-1]
                p2 = self.classM2Conc.predict_proba([m2c.toList()])[0][-1]
                self.trainConc.append((p1, p2, self.trainConcFraglen[i]))
                self.labConcM1.append(self.labM1Conc[i])
                self.labConcM2.append(self.labM2Conc[i])
        
            self.classConcM1 = KNeighborsClassifier(\
                n_neighbors=num_neighbors, warn_on_equidistant=False, weights=weights)
            self.classConcM2 = KNeighborsClassifier(\
                n_neighbors=num_neighbors, warn_on_equidistant=False, weights=weights)
            
            self.classConcM1.fit(self.trainConc, self.labM1Conc)
            self.classConcM2.fit(self.trainConc, self.labM2Conc)
    
    def save(self, fn):
        """ Save all training data to a file """
        save = (\
        self.trainUnp,    self.labUnp, \
        self.trainM1Disc, self.labM1Disc, \
        self.trainM2Disc, self.labM2Disc, \
        self.trainM1Conc, self.labM1Conc, \
        self.trainM2Conc, self.labM2Conc, \
        self.trainConc,   self.labConcM1,   self.labConcM2, \
        self.classConcM1, self.classConcM2, \
        self.trainConcFraglen)
        with open(fn, 'wb') as trainOfh: pickle.dump(save, trainOfh)
    
    def load(self, fn):
        """ Load all training data from a file """
        with open(fn, 'rb') as trainIfh:
            (\
            self.trainUnp,    self.labUnp, \
            self.trainM1Disc, self.labM1Disc, \
            self.trainM2Disc, self.labM2Disc, \
            self.trainM1Conc, self.labM1Conc, \
            self.trainM2Conc, self.labM2Conc, \
            self.trainConc,   self.labConcM1,   self.labConcM2, \
            self.classConcM1, self.classConcM2, \
            self.trainConcFraglen) = pickle.load(trainIfh)

class TrainingRecord(object):
    """ A single tuple of per-read training data including:
        1. Read length
        . Read's best score
        34. Read's second-best score """
    def __init__(self, rdlen, bestSc, secbestSc, scaleAs=1.0, scaleDiff=1.0):
        self.rdlen = rdlen
        self.bestSc = float(bestSc) * scaleAs
        self.secbestSc = secbestSc
        if secbestSc is not None:
            self.scDiff = float(abs(self.bestSc - self.secbestSc))
        else:
            self.scDiff = 10000.0
        self.scDiff *= scaleDiff
        assert self.repOk()
    
    @classmethod
    def fromAlignment(cls, al, scaleAs=1.0, scaleDiff=1.0):
        """ Initialize training record with respect to an alignment and a
            boolean indicating whether it is correct. """
        return cls(len(al), al.bestScore, al.secondBestScore,
                   scaleAs=scaleAs, scaleDiff=scaleDiff)
    
    def toList(self):
        """ Return simple list form """
        return [ self.rdlen, self.bestSc, self.scDiff ]
    
    def repOk(self):
        """ Check for internal consistency """
        assert self.rdlen is not None
        assert self.bestSc is not None
        assert self.scDiff is not None
        return True

class Output(threading.Thread):
    """ Encapsulates the output reader.  Reads SAM output from the aligner,
        updates empirical distributions for e.g. fragment length, and looks for
        records that correspond to simulated reads.  If a record corresponds to
        a simulated read, its correctness will be checked and a tuple
        written """
    
    def __init__(self, samIfh, samOfh, unalFh, trainSink,
                 scDistUnp, scDistM1, scDistM2, fragDist, typeDist):
        threading.Thread.__init__(self)
        self.samIfh = samIfh       # SAM records come from here
        self.samOfh = samOfh       # normal (non-simulated) SAM records go here
        self.trainSink = trainSink # write training data here
        self.fragDist = fragDist   # fragment length distribution
        self.scDistUnp = scDistUnp # score distribution for unpaired
        self.scDistM1 = scDistM1   # score distribution for mate #1s
        self.scDistM2 = scDistM2   # score distribution for mate #2s
        self.typeDist = typeDist   # alignment type (UU/CP/DP/UP) distribution
        self.unalFh = unalFh       # write unaligned reads here
        self.scDiffs = {}
    
    def run(self):
        lastAl = None
        for ln in self.samIfh:
            if ln.startswith('@'):
                if self.samOfh is not None:
                    self.samOfh.write(ln) # header line
                continue
            al = Alignment(ln)
            nm, flags, refid, pos, _, _, _, _, _, seq, qual, _ = string.split(ln, '\t', 11)
            if al.name.startswith('!!ts!!'):
                # this is a simulated read
                _, refid, fw, refoff, sc = string.split(al.name, '!!ts-sep!!')
                sc = int(sc)
                refoff = int(refoff)
                if al.isAligned():
                    scDiff = sc - al.bestScore
                    self.scDiffs[scDiff] = self.scDiffs.get(scDiff, 0) + 1
                    correct = False
                    # Check reference id, orientation
                    if refid == al.refid and fw == al.orientation():
                        # Check offset
                        correct = abs(refoff - al.pos) < args.wiggle
                    if self.trainSink is not None:
                        self.trainSink(al, correct)
                elif self.unalFh is not None:
                    # Perhaps write unaligned simulated read to file
                    self.unalFh.write("@%s\n%s\n+\n%s\n" % (nm, seq, qual))
            else:
                # Take alignment info into account
                if lastAl is not None and al.name == lastAl.name and al.alType == "CP":
                    assert lastAl.alType == "CP"
                    mate1, mate2 = al, lastAl
                    if (lastAl.flags & 64) != 0:
                        mate1, mate2 = mate2, mate1
                    alPair = AlignmentPair(mate1, mate2)
                    self.fragDist.add(alPair.fraglen()) # alignment pair
                if al.isAligned():
                    if (al.flags & 64) != 0:
                        self.scDistM1.add(al)
                    elif (al.flags & 128) != 0:
                        self.scDistM2.add(al)
                    else:
                        self.scDistUnp.add(al)
                    self.typeDist.add(al.alType)
                # Send SAM to SAM output filehandle
                if self.samOfh is not None:
                    self.samOfh.write(ln)
            lastAl = al

class AsyncWriter(threading.Thread):
    def __init__(self, fh, q, name=""):
        threading.Thread.__init__(self)
        self.fh = fh
        self.q = q
        self.name = name
    
    def run(self):
        i = 0
        while True:
            item = self.q.get()
            if item is None:
                self.q.task_done()
                self.fh.flush()
                self.fh.close()
                break
            self.fh.write(item)
            self.fh.flush()
            self.q.task_done()
            i += 1

class Aligner(object):
    """ Encapsulates the aligner """
    def __init__(self):
        pass

def createInput(rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2):
    """ Return an Input object that reads all user-provided input reads """
    if args.fastq:
        return InputWrapper(\
            Input(format="fastq", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2),
            rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2)
    elif args.fasta:
        return InputWrapper(\
            Input(format="fasta", unpFns=args.U, m1Fns=args.m1, m2Fns=args.m2),
            rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2)

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
    bt2_cmd += " --reorder --sam-no-qname-trunc -q "
    
    fifoArgs = ""
    fifoFns = [None, None, None] # FIFO filenames
    fifoFhs = [None, None, None] # FIFO filehandles
    fifoQs  = [None, None, None] # FIFO read queues
    fifoWs  = [None, None, None] # FIFO worker threads
    
    def makeFifos():
        fifoArgs = ''
        # Make temporary directory to store FIFOs
        tmpdir = tempfile.mkdtemp()
        # Make the appropriate FIFOs given the types of input we have
        if args.U is not None:
            fifoFns[0] = os.path.join(tmpdir, 'U.fifo')
            os.mkfifo(fifoFns[0])
            fifoArgs += (" -U " + fifoFns[0])
        if args.m1 is not None:
            fifoFns[1] = os.path.join(tmpdir, 'm1.fifo')
            os.mkfifo(fifoFns[1])
            fifoArgs += (" -1 " + fifoFns[1])
        if args.m2 is not None:
            fifoFns[2] = os.path.join(tmpdir, 'm2.fifo')
            os.mkfifo(fifoFns[2])
            fifoArgs += (" -2 " + fifoFns[2])
        return fifoArgs
    
    def openBowtie(bt2_cmd):
        bt2_cmd += makeFifos()
        # Open the Bowtie 2 process, which is going to want to start reading from
        # one or more of the FIFOs
        pipe = subprocess.Popen(bt2_cmd, shell=True, stdout=subprocess.PIPE)
        if args.U  is not None:
            fifoFhs[0] = open(fifoFns[0], 'w')
        if args.m1 is not None:
            fifoFhs[1] = open(fifoFns[1], 'w')
        if args.m2 is not None:
            fifoFhs[2] = open(fifoFns[2], 'w')
        # For each input type (unpaired, mate1, mate2), initialize a queue and a
        # thread that takes reads from the queue and passes each along to the
        # appropriate FIFO.  It's important to have a separate thread for each FIFO
        # or we get deadlocks.
        for i in xrange(0, 3):
            if fifoFhs[i] is not None:
                fifoQs[i] = Queue()
                fifoWs[i] = AsyncWriter(fifoFhs[i], fifoQs[i], "Thread %d" % i)
                fifoWs[i].start()
        return pipe
    
    def closeFifos():
        # Write None to all the FIFOs to inform them we're done giving them reads
        for i in xrange(0, 3):
            if fifoFhs[i] is not None:
                fifoQs[i].put(None)
        # Join all the FIFOs and input handling threads, close and delete the FIFOs 
        for i in xrange(0, 3):
            if fifoFhs[i] is not None:
                fifoQs[i].join()
                if fifoWs[i] is not None:
                    fifoWs[i].join()
                    fifoWs[i] = None
                if fifoFhs[i] is not None:
                    fifoFhs[i].close()
                    fifoFhs[i] = None
                os.unlink(fifoFns[i])
    
    # Structure encapsulating training data and KNN classifiers
    training = Training()
    
    if args.training_input is None:
        
        # The training data is not being given as input, so we have to
        # generate it with two Bowtie 2 runs.  The first run collects
        # some informations about data distributions.  The second run
        # draws from those distributions to simulate reads.  Each
        # simulated read that aligns is used as a training tuple.
        
        pipe = openBowtie(bt2_cmd)
        samIfh = tempfile.TemporaryFile()
        
        if args.verbose:
            print >> sys.stderr, "Real-data Bowtie 2 command: '%s'" % bt2_cmd
        
        scDistUnp, scDistM1, scDistM2, typeDist, fragDist = \
            ScoreDist(), ScoreDist(), ScoreDist(), Dist(), Dist()
        
        # Create the thread that eavesdrops on output from bowtie2
        othread = Output(
            pipe.stdout,     # SAM input filehandle
            samIfh,          # SAM output filehandle
            None,            # Write unaligned reads here
            None,            # Training record sink
            scDistUnp,       # Score dist for unpaired
            scDistM1,        # Score dist for mate1
            scDistM2,        # Score dist for mate2
            fragDist,        # Fragment dist
            typeDist)        # Alignment type dist
        othread.start()
        rddistUnp, rddistM1, rddistM2 = Dist(), Dist(), Dist()
        qdistUnp, qdistM1, qdistM2 = Dist(), Dist(), Dist()
        sctok = map(int, string.split(args.scoring, ','))
        
        # Create object that captures the scoring scheme used by the aligner
        scoring = Scoring(sctok[0], (sctok[1], sctok[2]), sctok[3], (sctok[4], sctok[5]), (sctok[6], sctok[7]))
        
        # Construct sequence and quality simulators
        sim, qsim = SequenceSimulator(args.ref, idx_fn=args.ref_idx), QualitySimulator(qdistUnp, qdistM1, qdistM2)
        simw = SimulatorWrapper(\
            sim,       # sequence simulator
            qsim,      # quality-value simulator
            typeDist,  # alignment type distribution
            scDistUnp, # score distribution for unpaired
            scDistM1,  # score distribution for mate 1
            scDistM2,  # score distribution for mate 2
            fragDist)  # fragment-length distribution
        
        print >> sys.stderr, "Initializing threads, queues and FIFOs"
        
        # Read through all the input read files and direct all reads to the
        # appropriate queue
        upto = args.upto or sys.maxint
        numReads = 0
        for (rd1, rd2) in iter(createInput(rddistUnp, rddistM1, rddistM2, qdistUnp, qdistM1, qdistM2)):
            numReads += 1
            if rd2 is not None:
                # Write to -1/-2 filehandles
                assert fifoFhs[1] is not None
                assert fifoFhs[2] is not None
                fifoQs[1].put(str(rd1))
                fifoQs[2].put(str(rd2))
            else:
                # Write to -U filehandle
                assert fifoFhs[0] is not None
                fifoQs[0].put(str(rd1))
            if numReads >= upto:
                break
        
        print >> sys.stderr, "Finished reading real reads"
    
        closeFifos()
        othread.join() # join the thread that monitors aligner output
        print >> sys.stderr, "Finished closing data FIFOs and joining output-monitoring thread"
        
        #
        # Now we re-open Bowtie 2 using the same arguments.  This time we only give
        # it simulated reads.
        #
        
        print >> sys.stderr, "Opening new Bowtie 2"
        pipe = openBowtie(bt2_cmd)
        
        # Function gets called with each new piece of training data
        def trainingSink(al, correct): training.add(al, correct)
        
        unalFh = None
        if args.un_sim is not None:
            unalFh = open(args.un_sim, 'w')
        
        # Create the thread that eavesdrops on output from bowtie2 with simulated
        # input
        othread = Output(\
            pipe.stdout,     # SAM input filehandle
            None,            # SAM output filehandle
            unalFh,          # Write unaligned reads here
            trainingSink,    # Training record sink
            None,            # Score dist for unpaired
            None,            # Score dist for mate1
            None,            # Score dist for mate2
            None,            # Fragment dist
            None)            # Alignment type dist
        othread.start()
        
        # Simulate reads from empirical distributions
        for i in xrange(0, args.num_reads):
            if (i+1 % 1000) == 0: print "Generating read %d" % i
            rd1, rd2 = simw.next()
            if rd2 is not None:
                # Paired-end simulated read
                fifoQs[1].put(str(rd1))
                fifoQs[2].put(str(rd2))
            else:
                # Unpaired simulated read
                fifoQs[0].put(str(rd1))
        
        print >> sys.stderr, "Finished simulating reads"
        closeFifos()
        print >> sys.stderr, "Closed FIFOs"
        othread.join() # join the thread that monitors aligner output
        if unalFh is not None:
            unalFh.close()
        print >> sys.stderr, "Finished closing simulation FIFOs and joining output-monitoring thread"
        
        print >> sys.stderr, "Score difference (expected - actual) histogram:"
        for k, v in sorted(othread.scDiffs.iteritems()):
            print >>sys.stderr, "  %d: %d" % (k, v)
        
        if args.save_training is not None: training.save(args.save_training)
    else:
        # Training data is being given to us
        if args.sam_input is None:
            raise RuntimeError("Must specify --sam-input along with --training-input")
        training.load(args.training_input)
        print >> sys.stderr, "Read %d training tuples from '%s'" % (len(training), args.training_input)
        samIfh = open(args.sam_input, 'r')
    
    # Build KNN classifiers
    weights = 'distance' if args.distance_weight else 'uniform'
    training.fit(args.num_neighbors, weights, scaleAs=args.as_scale, scaleDiff=args.diff_scale)
    print >> sys.stderr, "Finished fitting KNN classifiers on %d training tuples" % len(training)
    
    # TODO: Need to pair up alignments so we can potentially use the
    # concordant-alignment model to calculate their MAPQs
    nrecs = 0
    with open(args.S, 'w') as samOfh:
        samIfh.seek(0)
        for samrec in samIfh:
            if samrec.startswith('@'):
                samOfh.write(samrec)
                continue
            al = Alignment(samrec)
            nrecs += 1
            if al.isAligned():
                probCorrect = training.probCorrect(al)
                probIncorrect = 1.0 - probCorrect # convert to probability incorrect
                mapq = args.max_mapq
                if probIncorrect > 0.0:
                    mapq = min(-10.0 * math.log10(probIncorrect), args.max_mapq)
                if args.verbose:
                    if True or rec.bestSc == rec.secbestSc:
                        print "=="
                        print rec.toList()
                        print "--"
                        print mapq
                samrec = samrec.rstrip()
                xqIdx = samrec.find("\tXQ:f:")
                if xqIdx != -1:
                    samrec = samrec[:xqIdx]
                samOfh.write("\t".join([samrec, "XQ:f:" + str(mapq)]) + "\n")
            else:
                samOfh.write(samrec)
    
    print >> sys.stderr, "Finished writing final SAM output (%d records) to '%s'" % (nrecs, args.S)
    
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
        '--fasta', action='store_const', const=True, default=False, help='Reads are FASTA')
    parser.add_argument(\
        '--fastq', action='store_const', const=True, default=True, help='Reads are FASTQ')
    parser.add_argument(\
        '--scoring', metavar='str', type=str, required=False,
        default='1,2,6,1,5,3,5,3',
        help='MatchBonus,MismatchMinPen,MismatchMaxPen,NPen,ReadGapConst,ReadGapLinear,RefGapConst,RefGapLinear')
    parser.add_argument(\
        '--max-mapq', metavar='float', type=float, default=100.0,
        required=False, help='Maximum MAPQ possible for an alignment')
    parser.add_argument(\
        '--seed', metavar='int', type=int, default=99099,
        required=False, help='Integer to initialize pseudo-random generator')
    parser.add_argument(\
        '--num-reads', metavar='int', type=int, default=100000,
        required=False, help='Number of reads to simulate')
    parser.add_argument(\
        '--num-neighbors', metavar='int', type=int, default=500,
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
