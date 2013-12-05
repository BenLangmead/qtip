import os
import re
import random
import cPickle
import logging
import sys
import string
from randutil import WeightedRandomGenerator
from read import Read

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

# TODO: get rid of this version
def mutateOld(rd, rdfw, scDistDraw):
    ''' Given a read that already has the appropriate length (i.e. equal to #
        characters on the reference side of the alignment), take the alignment
        information contained in the scDistDraw object and modify rd to
        contain the same pattern of edits.  Modifies rd in place. '''
    fw, qual, rdAln, rfAln, sc = scDistDraw
    #assert 'N' not in rd.seq
    #assert len(rd.seq) == len(rfAln) - rfAln.count('-')
    if rdfw != fw:
        qual, rdAln, rfAln = qual[::-1], rdAln[::-1], rfAln[::-1]
    rd.qual = qual # Install qualities
    # Walk along stacked alignment, making corresponding changes to
    # read sequence
    seq = []
    i, rdi, rfi = 0, 0, 0
    while i < len(rdAln):
        if rdAln[i] == rfAln[i] and rfAln[i] != 'N':
            seq.append(rd.seq[rfi])
            rfi += 1; rdi += 1
        elif rdAln[i] == '-':
            rfi += 1
        elif rfAln[i] == '-':
            seq.append(rdAln[i])
            rdi += 1
        elif rdAln[i] != rfAln[i] or rdAln[i] == 'N':
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
            raise RuntimeError('Should not get here')
        i += 1
    rd.seq = ''.join(seq)
    assert len(rd.seq) == len(rd.qual), "%s\n%s\n%s\n%s" % (rdAln, rfAln, rd.seq, rd.qual)

def mutate(rd, rdfw, scDistDraw):
    ''' Given a read that already has the appropriate length (i.e. equal to #
        characters on the reference side of the alignment), take the alignment
        information contained in the scDistDraw object and modify rd to
        contain the same pattern of edits.  Modifies rd in place. '''
    fw, qual, rdAln, rfAln, sc = scDistDraw
    if rdfw != fw:
        qual, rdAln, rfAln = qual[::-1], rdAln[::-1], rfAln[::-1]
    rd.qual = qual # Install qualities
    # Walk along stacked alignment, making corresponding changes to
    # read sequence
    seq = []
    i, rdi, rfi = 0, 0, 0
    while i < len(rdAln):
        if rdAln[i] == rfAln[i] and rfAln[i] != 'N':
            seq.append(rd.seq[rfi])
            rfi += 1; rdi += 1
        elif rdAln[i] == '-':
            rfi += 1
        elif rfAln[i] == '-':
            seq.append(rdAln[i])
            rdi += 1
        elif rdAln[i] != rfAln[i] and (rdAln[i] == 'N'):
            seq.append('N')
            rfi += 1; rdi += 1
        elif rfAln[i] == 'N':
            seq.append(random.choice('ACGT'))
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
            raise RuntimeError('Should not get here')
        i += 1
    rd.seq = ''.join(seq)
    assert len(rd.seq) == len(rd.qual), "%s\n%s\n%s\n%s" % (rdAln, rfAln, rd.seq, rd.qual)

if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False, help='Do unit tests')
    
    args = parser.parse_args()
    
    if args.test:
        import unittest
        
        random.seed(100)
        
        class TestMutate(unittest.TestCase):
            
            def test_1(self):
                # No edits
                rd = Read('r1', 'ACGT', '!!!!')
                mutate(rd, True, (True, 'IIII', 'TTTT', 'TTTT', 10))
                self.assertEqual(rd.seq, 'ACGT')
                self.assertEqual(rd.qual, 'IIII')
            
            def test_2(self):
                # 1 mismatch
                rd = Read('r2', 'ACGT', '!!!!')
                mutate(rd, True, (True, 'IIII', 'TTTT', 'TTTA', 10))
                self.assertNotEqual(rd.seq, 'ACGT')
                self.assertEqual(rd.qual, 'IIII')
            
            def test_3(self):
                # 1 read gap
                rd = Read('r3', 'A' * 12, '!' * 12)
                # rd: TTTTT--TTTTT
                # rf: TTTTTTTTTTTT
                mutate(rd, True, (True, 'I' * 10, 'TTTTT--TTTTT', 'T' * 12, 10))
                self.assertEqual(rd.seq, 'A' * 10)
                self.assertEqual(rd.qual, 'I' * 10)
            
            def test_4(self):
                # 1 reference gap
                rd = Read('r4', 'C' * 10, '#' * 10)
                # rd: TTTTTTTTTTTT
                # rf: TTTTT--TTTTT
                mutate(rd, True, (True, 'J' * 12, 'TTTTTTTTTTTT', 'TTTTT--TTTTT', 10))
                self.assertEqual(rd.seq, 'CCCCCTTCCCCC')
                self.assertEqual(rd.qual, 'J' * 12)
            
            def test_5(self):
                # N in the reference
                rd = Read('r5', 'CCCCC', '#####')
                # rd: TTTTT
                # rf: TTNTT
                mutate(rd, True, (True, 'KKKKK', 'TTTTT', 'TTNTT', 10))
                self.assertEqual(rd.seq[:2], 'CC')
                self.assertEqual(rd.seq[-2:], 'CC')
                self.assertNotEqual(rd.seq[2], 'N')
                self.assertEqual(rd.qual, 'KKKKK')
            
            def test_6(self):
                # N in the read
                rd = Read('r6', 'CCCCC', '#####')
                # rd: TTNTT
                # rf: TTTTT
                mutate(rd, True, (True, 'KKKKK', 'TTNTT', 'TTTTT', 10))
                self.assertEqual(rd.seq, 'CCNCC')
                self.assertEqual(rd.qual, 'KKKKK')
        
        unittest.main(argv=[sys.argv[0]])
