'''
check_wgsim.py

Given a FASTQ file simulated with wgsim and the collection of FASTA files
simulated from, check that the simulated reads seem to match the reference in
the proper spot.
'''

import sys
import re
import gzip
import string
from align import editDistance
from reference import ReferenceIndexed

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

if len(sys.argv) < 4:
    raise RuntimeError('Not enough arguments; need 2 wgsim files followed by FASTA file(s)')

wgsimFn1, wgsimFn2 = sys.argv[1], sys.argv[2]
ref = ReferenceIndexed(sys.argv[3:])

nm_re = re.compile('@([^_]+)_([^_]+)_([^_]+)_([^:]+):([^:]+):([^_]+)_([^:]+):([^:]+):([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^/]+)/([12])')

with gzip.open(wgsimFn1) if wgsimFn1.endswith('.gz') else open(wgsimFn1) as fh1:
    with gzip.open(wgsimFn2) if wgsimFn2.endswith('.gz') else open(wgsimFn2) as fh2:
        while True:
            l1_1, l1_2 = fh1.readline().rstrip(), fh2.readline().rstrip()
            l2_1, l2_2 = fh1.readline().rstrip(), fh2.readline().rstrip()
            l3_1, l3_2 = fh1.readline(), fh2.readline()
            l4_1, l4_2 = fh1.readline(), fh2.readline()
            if len(l4_1) == 0:
                break
            # @chr1_64253297_64253619_3:0:0_2:1:0_0_0/1
            #                         ***** ***** !
            # ***** # edits
            # ! whether it's "flipped"
            res1 = nm_re.match(l1_1)
            refid, fragst1, fragen1 = res1.group(1), int(res1.group(2)), int(res1.group(3))
            flip = res1.group(12) == '1'
            assert fragen1 > fragst1
            ref1 = ref.get(refid, fragst1-1, len(l2_1))
            ref2 = ref.get(refid, fragen1-len(l2_2), len(l2_2))
            if flip:
                ref1, ref2 = ref2, ref1
                ref1 = revcomp(ref1)
            else:
                ref2 = revcomp(ref2)
            print 'flipped: ' + str(flip)
            print 'rd1: ' + l2_1
            print 'rf1: ' + ref1
            print 'ed1: ' + str(editDistance(l2_1.upper(), ref1.upper()))
            print 'rd2: ' + l2_2
            print 'rf2: ' + ref2
            print 'ed2: ' + str(editDistance(l2_2.upper(), ref2.upper()))
