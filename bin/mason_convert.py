"""
mason_convert.py

Take FASTQ files output by Mason.  Emit a new FASTQ file that uses our
modified wgsim-like read name encoding.

Usage:
python mason_convert.py --in1 *.fastq [--in2 *.fastq] \
                        --out1 out1.fastq [--out2 out2.fastq]
"""

from __future__ import print_function
import re
import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description='Convert Mason FASTQ files to our modified wgsim-like encoding')
parser.add_argument(
    '--in1', metavar='path', type=str, nargs='+', required=True, help='Mate #1s or unpaired reads.')
parser.add_argument(
    '--in2', metavar='path', type=str, nargs='+', required=False, help='Mate #2s.')
parser.add_argument(
    '--out1', metavar='path', type=str, required=True, help='Output for mate #1s or unpaired reads.')
parser.add_argument(
    '--out2', metavar='path', type=str, required=False, help='Output for mate #2s.')
args = parser.parse_args()

_mason_orig_beg = re.compile('orig_begin=([0-9]*)')
_mason_orig_end = re.compile('orig_end=([0-9]*)')
_mason_contig = re.compile('contig=([^\s]*)')
_mason_strand = re.compile('strand=([^\s]*)')


def parse_mason(nm):
    be = _mason_orig_beg.search(nm)
    en = _mason_orig_end.search(nm)
    assert be is not None and en is not None, nm
    
    # Mason's offsets are 0-based
    left, right = int(be.group(1)), int(en.group(1))
    
    rr = _mason_contig.search(nm)
    assert rr is not None
    refid = rr.group(1)
    
    sr = _mason_strand.search(nm)
    assert sr is not None
    strand = sr.group(1)
    
    # Convert from 0-based to 1-based
    return left+1, right, refid, strand == 'forward'


def make_wgsim(refid, posl, posr, len1, len2, flipped, idx, mate):
    return '%s_%d_%d_%d:%d:%d_%d:%d:%d_%d_%d_%d_%d/%d' %\
           (refid, posl, posr, 0, 0, 0, 0, 0, 0, len1, len2, flipped, idx, mate)


def go():
    idx = 1
    n_tot, n_ref_id_mismatch, n_strands_match, n_strands_not_compat = 0, 0, 0, 0
    n_discarded = 0
    if args.in2 is not None:
        for infn1, infn2 in zip(args.in1, args.in2):
            with gzip.open(infn1) if infn1.endswith('.gz') else open(infn1) as infh1:
                with gzip.open(infn2) if infn2.endswith('.gz') else open(infn2) as infh2:
                    with open(args.out1, 'w') as ofh1:
                        with open(args.out2, 'w') as ofh2:
                            while True:
                                ln_1 = infh1.readline().rstrip()[1:]
                                seq_1 = infh1.readline()
                                ln3_1 = infh1.readline()
                                ln4_1 = infh1.readline()
                                if len(ln4_1) == 0:
                                    break
                                ln_2 = infh2.readline().rstrip()[1:]
                                seq_2 = infh2.readline()
                                ln3_2 = infh2.readline()
                                ln4_2 = infh2.readline()
                                posl_1, posr_1, refid_1, strand_1 = parse_mason(ln_1)
                                posl_2, posr_2, refid_2, strand_2 = parse_mason(ln_2)
                                n_tot += 1
                                if refid_1 != refid_2:
                                    n_ref_id_mismatch += 1
                                    n_discarded += 1
                                    continue
                                elif strand_1 == strand_2:
                                    n_strands_match += 1
                                    n_discarded += 1
                                    continue
                                elif (posl_1 < posl_2) and strand_2 or (posr_2 < posr_1) and strand_1:
                                    n_strands_not_compat += 1
                                    n_discarded += 1
                                    continue
                                assert refid_1 == refid_2
                                nm_1 = make_wgsim(refid_1, min(posl_1, posl_2), max(posr_1, posr_2), len(seq_1) - 1,
                                                  len(seq_2) - 1, strand_2, idx, 1)
                                nm_2 = make_wgsim(refid_1, min(posl_1, posl_2), max(posr_1, posr_2), len(seq_1) - 1,
                                                  len(seq_2) - 1, strand_2, idx, 2)
                                ofh1.write("@%s\n" % nm_1)
                                ofh2.write("@%s\n" % nm_2)
                                for ln in [seq_1, ln3_1, ln4_1]:
                                    ofh1.write(ln)
                                for ln in [seq_2, ln3_2, ln4_2]:
                                    ofh2.write(ln)
                                idx += 1
    else:
        for infn1 in args.in1:
            with gzip.open(infn1) if infn1.endswith('.gz') else open(infn1) as infh1:
                with open(args.out1, 'w') as ofh1:
                    while True:
                        ln = infh1.readline().rstrip()[1:]
                        seq = infh1.readline()
                        ln3 = infh1.readline()
                        ln4 = infh1.readline()
                        if len(ln4) == 0:
                            break
                        n_tot += 1
                        posl, posr, refid, strand = parse_mason(ln)
                        nm = make_wgsim(refid, posl, posr, len(seq) - 1, len(seq) - 1, not strand, idx, 1)
                        ofh1.write("@%s\n" % nm)
                        for ln in [seq, ln3, ln4]:
                            ofh1.write(ln)
                        idx += 1
    print('ref id mismatch: %d (%0.4f%%)' %
          (n_ref_id_mismatch, 100.0 * n_ref_id_mismatch / n_tot), file=sys.stderr)
    print('strands match: %d (%0.4f%%)' %
          (n_strands_match, 100.0 * n_strands_match / n_tot), file=sys.stderr)
    print('strands not compat: %d (%0.4f%%)' %
          (n_strands_not_compat, 100.0 * n_strands_not_compat / n_tot), file=sys.stderr)
    print('reads converted: %d (%0.4f%%)' %
          (n_tot - n_discarded, 100.0 * (n_tot - n_discarded) / n_tot), file=sys.stderr)

go()
