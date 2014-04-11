#!/usr/bin/env python

"""
fastq_sample.py

Author: Ben Langmead
Date: 7/3/2013
Contact: langmea@cs.jhu.edu

Sample --n reads from FASTQ inputs specified with --in.  Write to --out.
"""

import sys
import gzip
import random
import argparse
    
parser = argparse.ArgumentParser(description='Sample FASTQ reads from file')
parser.add_argument(\
    '--in', dest='input', metavar='path', type=str, nargs='+', required=True, help='FASTQ file input.')
parser.add_argument(\
    '--out', metavar='path', type=str, required=False, help='Output.')
parser.add_argument(\
    '--seed', metavar='integer', type=int, default=92039, help='Pseudo-random seed')
parser.add_argument(\
    '--n', metavar='integer', type=int, required=True, help='Number of reads to sample')
args = parser.parse_args()


class ReservoirSampler(object):
    """ Simple reservoir sampler """
    def __init__(self, k):
        self.k = k  # # elts to collect
        self.r = [] # the elts
        self.n = 0  # # elts scanned
    
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
    
    def shuffle(self):
        random.shuffle(self.r)

random.seed(args.seed)
sampler = ReservoirSampler(args.n)
for ifn in args.input:
    fh = gzip.GzipFile(ifn, 'r') if ifn.endswith('.gz') else open(ifn, 'r')
    print >> sys.stderr, "Processing '%s' ..." % ifn
    while True:
        name1 = fh.readline()
        if len(name1) == 0: break
        seq = fh.readline()
        name2 = fh.readline()
        qual = fh.readline()
        sampler.add((name1, seq, name2, qual))
    fh.close()

print >> sys.stderr, "Shuffling ..."
sampler.shuffle()

print >> sys.stderr, "Writing output ..."
ofn = args.out
ofh = sys.stdout if ofn is None else open(ofn, 'w')
for rd in sampler.r:
    ofh.write(''.join(rd))
ofh.close()
