import os
import re
import random
import cPickle
import logging
import sys
import string
from randutil import WeightedRandomGenerator

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

class SequenceSimulator(object):
    ''' Class that, given a collection of FASTA files, samples intervals of
        specified length from the strings contained in them.  Simulated reads
        are not permitted to overlap a non-A/C/G/T character in the reference.
    '''
    
    def __init__(self, ref):
        self.__re = re.compile('[^ACGTacgt]')
        self.ref = ref
        lens = {}
        for name in ref.names():
            lens[name] = ref.length(name)
        self.rnd = WeightedRandomGenerator(lens)
    
    def sim(self, ln, verbose=False):
        ''' Simulate a read '''
        # Pick a reference sequence in a weighted random fashion
        attempts = 0
        while True:
            refi = self.rnd.next()
            assert self.ref.hasName(refi)
            if self.ref.length(refi) >= ln:
                break
            if attempts > 5:
                raise RuntimeError('Tried more than %d times to find reference with length at least %d' % (attempts, ln))
            attempts += 1
        fw = True
        refoff = random.randint(0, self.ref.length(refi) - ln) # pick offset
        seq = self.ref.get(refi, refoff, ln) # extract substring
        # Simulated read can't overlap non-A-C-G-T character in reference
        while self.__re.search(seq):
            refoff = random.randint(0, self.ref.length(refi) - ln) # pick new offset
            seq = self.ref.get(refi, refoff, ln) # extract substring again
        seq = seq.upper()
        if random.random() > 0.5: # possibly reverse-complement
            fw = False
            seq = revcomp(seq) # reverse complement
        assert not "N" in seq
        return (refi, refoff, fw, seq) # return ref id, ref offset, orientation, sequence

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
    assert len(rd.seq) == len(rd.qual), "%s\n%s\n%s\n%s" % (rdAln, rfAln, rd.seq, rd.qual)
