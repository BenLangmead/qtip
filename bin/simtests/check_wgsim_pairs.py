"""
check_mason.py

Given a FASTQ file simulated with Mason and the collection of FASTA files
simulated from, check that the simulated reads seem to match the reference in
the proper spot.
"""

import sys
import gzip
import re
from collections import defaultdict

_wgsimex_re = re.compile('(.+)_([^_]+)_([^_]+)_([^:]+):([^:]+):([^_]+)_([^:]+):([^:]+):([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^/]+).*')
#                           1   2       3       4       5       6       7       8       9       10      11      12      13


def is_extended_wgsim(nm):
    ret = _wgsimex_re.match(nm) is not None
    return ret


def parse_extended_wgsim(name, mate1):
    """ Note: this is my extended version of wgsim's output """
    res = _wgsimex_re.match(name)
    refid, fragst1, fragen1 = res.group(1), int(res.group(2))-1, int(res.group(3))-1
    len1, len2 = int(res.group(10)), int(res.group(11))
    flip = res.group(12) == '1'
    ln = len1 if mate1 else len2
    if (not flip) == mate1:
        return fragst1, refid, True, fragen1 - fragst1 + 1
    else:
        return fragen1 - (ln-1), refid, False, fragen1 - fragst1 + 1

if len(sys.argv) < 3:
    raise RuntimeError('Not enough arguments; need 2 FASTQ files')

fastq_fn_1 = sys.argv[1]
fastq_fn_2 = sys.argv[2]

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
            left_1, refid_1, strand_1, fraglen = parse_extended_wgsim(l1_1, True)
            left_2, refid_2, strand_2, _ = parse_extended_wgsim(l1_2, False)
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
            n_fraglen[fraglen] += 1

print 'ref id mismatch: %d (%0.4f%%)' % (n_ref_id_mismatch, 100.0 * n_ref_id_mismatch / n_tot)
print 'strands match: %d (%0.4f%%)' % (n_strands_match, 100.0 * n_strands_match / n_tot)
print 'strands not compat: %d (%0.4f%%)' % (n_strands_not_compat, 100.0 * n_strands_not_compat / n_tot)
print 'fraglens:'
for fraglen, cnt in sorted(n_fraglen.iteritems()):
    print '  %d: %d' % (fraglen, cnt)
