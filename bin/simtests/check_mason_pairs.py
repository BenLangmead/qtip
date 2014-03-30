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
from collections import defaultdict

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

if len(sys.argv) < 3:
    raise RuntimeError('Not enough arguments; need 2 FASTQ files')

fastq_fn_1 = sys.argv[1]
fastq_fn_2 = sys.argv[2]

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

n_tot, n_ref_id_mismatch, n_strands_match, n_strands_not_compat = 0, 0, 0, 0
n_fraglen = defaultdict(int)
with gzip.open(fastq_fn_1) if fastq_fn_1.endswith('.gz') else open(fastq_fn_1) as fh1:
    with gzip.open(fastq_fn_2) if fastq_fn_2.endswith('.gz') else open(fastq_fn_2) as fh2:
        while True:
            l1_1 = fh1.readline().rstrip()
            l2_1 = fh1.readline().rstrip()
            l3_1 = fh1.readline()
            l4_1 = fh1.readline()
            l1_2 = fh2.readline().rstrip()
            l2_2 = fh2.readline().rstrip()
            l3_2 = fh2.readline()
            l4_2 = fh2.readline()
            if len(l4_1) == 0:
                break
            left_1, right_1, refid_1, strand_1 = parse_mason(l1_1)
            left_2, right_2, refid_2, strand_2 = parse_mason(l1_2)
            n_tot += 1
            if refid_1 != refid_2:
                n_ref_id_mismatch += 1
                continue
            elif strand_1 == strand_2:
                n_strands_match += 1
                continue
            elif (left_1 < left_2) and strand_2 or (left_2 < left_1) and strand_1:
                print l1_1
                print l1_2
                print '---'
                n_strands_not_compat += 1
                continue
            fraglen = max(right_1, right_2) - min(left_1, left_2)
            n_fraglen[fraglen] += 1

print 'ref id mismatch: %d (%0.4f%%)' % (n_ref_id_mismatch, 100.0 * n_ref_id_mismatch / n_tot)
print 'strands match: %d (%0.4f%%)' % (n_strands_match, 100.0 * n_strands_match / n_tot)
print 'strands not compat: %d (%0.4f%%)' % (n_strands_not_compat, 100.0 * n_strands_not_compat / n_tot)
print 'fraglens:'
for fraglen, cnt in sorted(n_fraglen.iteritems()):
    print '  %d: %d' % (fraglen, cnt)
