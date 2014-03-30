"""
check_mason.py

Given a FASTQ file simulated with Mason and the collection of FASTA files
simulated from, check that the simulated reads seem to match the reference in
the proper spot.
"""

import sys
import gzip
import string
import re
from align import editDistance
from reference import ReferenceIndexed

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

if len(sys.argv) < 3:
    raise RuntimeError('Not enough arguments; need 1 FASTQ file followed by FASTA file(s)')

fastqFn = sys.argv[1]
ref = ReferenceIndexed(sys.argv[2:])

_mason_orig_beg = re.compile('orig_begin=([0-9]*)')
_mason_orig_end = re.compile('orig_end=([0-9]*)')
_mason_contig = re.compile('contig=([^\s]*)')
_mason_strand = re.compile('strand=([^\s]*)')


def parse_mason(nm):
    be = _mason_orig_beg.search(nm)
    en = _mason_orig_end.search(nm)
    assert be is not None and en is not None
    
    # Mason's offsets are 0-based
    left, right = int(be.group(1)), int(en.group(1))
    
    rr = _mason_contig.search(nm)
    assert rr is not None
    refid = rr.group(1)
    
    sr = _mason_strand.search(nm)
    assert sr is not None
    strand = sr.group(1)
    
    return left, right, refid, strand == 'forward'

with gzip.open(fastqFn) if fastqFn.endswith('.gz') else open(fastqFn) as fh:
    while True:
        l1_1 = fh.readline().rstrip()
        l2_1 = fh.readline().rstrip()
        l3_1 = fh.readline()
        l4_1 = fh.readline()
        if len(l4_1) == 0:
            break
        left, right, refid, strand = parse_mason(l1_1)
        refstr = ref.get(refid, left, right - left)
        if not strand:
            refstr = revcomp(refstr)
        print 'rd: ' + l2_1
        print 'rf: ' + refstr
        print 'ed: ' + str(editDistance(l2_1.upper(), refstr.upper()))
