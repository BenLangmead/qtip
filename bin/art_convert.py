'''
art_convert.py

Take FASTQ and SAM files output by Art.  Emit a new FASTQ file that
uses our modified wgsim-like read name encoding.

Usage:
python art_convert.py --in1 *.fastq [--in2 *.fastq] --sam *.sam \
                      --out1 out1.fastq [--out2 out2.fastq]
'''

import sys
import re
import gzip
import argparse

parser = argparse.ArgumentParser(description='Convert Mason FASTQ files to our modified wgsim-like encoding')
parser.add_argument(\
    '--in1', metavar='path', type=str, nargs='+', required=True, help='Mate #1s or unpaired reads from art_illumina.')
parser.add_argument(\
    '--in2', metavar='path', type=str, nargs='+', required=False, help='Mate #2s from art_illumina.')
parser.add_argument(\
    '--sam', metavar='path', type=str, nargs='+', required=True, help='SAM output from art_illumina -sam.')
parser.add_argument(\
    '--out1', metavar='path', type=str, required=True, help='Output for mate #1s or unpaired reads.')
parser.add_argument(\
    '--out2', metavar='path', type=str, required=False, help='Output for mate #2s.')
args = parser.parse_args()

def makeWgsim(refid, posl, posr, len1, len2, flipped, id, mate):
    return '%s_%d_%d_%d:%d:%d_%d:%d:%d_%d_%d_%d_%d/%d' % (refid, posl, posr, 0, 0, 0, 0, 0, 0, len1, len2, flipped, id, mate)

id = 1
if args.in2 is not None:
    # Paired-end case
    for infn1, infn2, samfn in zip(args.in1, args.in2, args.sam):
        with gzip.open(infn1) if infn1.endswith('.gz') else open(infn1) as infh1:
            with gzip.open(infn2) if infn2.endswith('.gz') else open(infn2) as infh2:
                with gzip.open(samfn) if samfn.endswith('.gz') else open(samfn) as samfh:
                    with open(args.out1, 'w') as ofh1:
                        with open(args.out2, 'w') as ofh2:
                            while True:
                                # Read FASTQ records
                                l1_1 = infh1.readline().rstrip()[1:]
                                l2_1 = infh1.readline()
                                l3_1 = infh1.readline()
                                l4_1 = infh1.readline()
                                l1_2 = infh2.readline().rstrip()[1:]
                                l2_2 = infh2.readline()
                                l3_2 = infh2.readline()
                                l4_2 = infh2.readline()
                                if len(l4_1) == 0: break
                                
                                # Read pair of SAM records
                                sam1 = '@'
                                while sam1[0] == '@':
                                    sam1 = samfh.readline()
                                    if len(sam1) == 0:
                                        raise RuntimeError('Ran out of SAM')
                                sam2 = samfh.readline()
                                if len(sam2) == 0:
                                    raise RuntimeError('Ran out of SAM')
                                
                                # Parse SAM
                                nm_1, flag_1, refid_1, off_1, _, _, _, _, _, seq_1, qual_1 = sam1.split('\t')[:12]
                                nm_2, flag_2, refid_2, off_2, _, _, _, _, _, seq_2, qual_2 = sam2.split('\t')[:12]
                                flag_1, flag_2 = int(flag_1), int(flag_2)
                                off_1, off_2 = int(off_1), int(off_2)
                                assert refid_1 == refid_2, ("\n%s\n%s" % (sam1, sam2))
                                
                                flipped = (flag_1 & 16) != 0
                                
                                refid = refid_1
                                posl = min(off_1, off_2)
                                posr = max(off_1 + len(seq_1), off_2 + len(seq_2)) - 1
                                
                                nm_1 = makeWgsim(refid, posl, posr, len(seq_1), len(seq_2), flipped, id, 1)
                                nm_2 = makeWgsim(refid, posl, posr, len(seq_1), len(seq_2), flipped, id, 2)
                                ofh1.write("@%s\n" % nm_1)
                                ofh2.write("@%s\n" % nm_2)
                                for ln in [l2_1, l3_1, l4_1]: ofh1.write(ln)
                                for ln in [l2_2, l3_2, l4_2]: ofh2.write(ln)
                                id += 1
else:
    # Unpaired case
    for infn, samfn in zip(args.in1, args.sam):
        with gzip.open(infn) if infn.endswith('.gz') else open(infn) as infh:
            with gzip.open(samfn) if samfn.endswith('.gz') else open(samfn) as samfh:
                with open(args.out1, 'w') as ofh:
                    while True:
                        # Read FASTQ records
                        l1 = infh.readline().rstrip()[1:]
                        l2 = infh.readline()
                        l3 = infh.readline()
                        l4 = infh.readline()
                        if len(l4) == 0: break
                        
                        # Read pair of SAM records
                        sam = '@'
                        while sam[0] == '@':
                            sam = samfh.readline()
                            if len(sam) == 0:
                                raise RuntimeError('Ran out of SAM')
                        
                        # Parse SAM
                        nm, flag, refid, off, _, _, _, _, _, seq, qual = sam.split('\t')[:12]
                        flag = int(flag)
                        off = int(off)
                        flipped = (flag & 16) != 0
                        
                        posl = off
                        posr = off + len(seq) - 1
                        nm = makeWgsim(refid, posl, posr, len(seq), len(seq), flipped, id, 1)
                        ofh.write("@%s\n" % nm)
                        for ln in [l2, l3, l4]: ofh.write(ln)
                        id += 1
