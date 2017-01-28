'''
fastq_interleave.py

Given two (possibly compressed) FASTQ files, print interleaved FASTQ
to stdout.
'''

import sys

def openex(fn, mode="rb"):
    if fn.endswith(".gz"):
        import gzip
        return gzip.open(fn, mode)
    elif fn.endswith(".bz2"):
        import bz2
        return bz2.BZ2File(fn, mode)
    else:
        return open(fn, mode)

if len(sys.argv) != 3:
    raise RuntimeError('Must specify 2 input FASTQ files')

ofh = sys.stdout
with openex(sys.argv[1], 'r') as fh1:
    with openex(sys.argv[2], 'r') as fh2:
        while True:
            ln1_1 = fh1.readline().rstrip()
            ln1_2 = fh2.readline().rstrip()
            ln2_1 = fh1.readline()
            ln2_2 = fh2.readline()
            ln3_1 = fh1.readline()
            ln3_2 = fh2.readline()
            ln4_1 = fh1.readline()
            ln4_2 = fh2.readline()
            if len(ln4_1) == 0:
                assert len(ln4_2) == 0
                assert len(ln2_1) == 0
                assert len(ln2_2) == 0
                break
            if not ln1_1.endswith('/1'):
                ln1_1 += '/1'
            if not ln1_2.endswith('/2'):
                ln1_2 += '/2'
            ofh.write(''.join([ln1_1, '\n', ln2_1, ln3_1, ln4_1,
                               ln1_2, '\n', ln2_2, ln3_2, ln4_2]))
