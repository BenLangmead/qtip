'''
check_art.py

Given a FASTQ file and associated SAM file simulated with Art, and the
collection of FASTA files simulated from, check that the simulated
reads seem to match the reference in the proper spot.
'''

import sys
import gzip
import string
from align import editDistance
from reference import ReferenceIndexed

_revcomp_trans = string.maketrans("ACGTacgt", "TGCAtgca")
def revcomp(s):
    return s[::-1].translate(_revcomp_trans)

if len(sys.argv) < 4:
    raise RuntimeError('Not enough arguments; need 1 FASTQ file followed by 1 SAM file, then FASTA file(s)')

fqfn = sys.argv[1]
samfn = sys.argv[2]
ref = ReferenceIndexed(sys.argv[3:])

with gzip.open(fqfn) if fqfn.endswith('.gz') else open(fqfn) as fqfh:
    with gzip.open(samfn) if samfn.endswith('.gz') else open(samfn) as samfh:
        while True:
            l1 = fqfh.readline().rstrip()
            l2 = fqfh.readline().rstrip()
            l3 = fqfh.readline()
            l4 = fqfh.readline()
            if len(l4) == 0: break
            # From mate 1 file:
            #
            # @chr1-24747/1
            # TGGTTAATGTAGTTGCTTCATTGTGTCATTGATCTGTCTATTTCAGTGTGTTTTTGTAGTGGCTCGTAAAGGTTTTTCCTTTTCCTTATATCACCAGGAG
            # +
            # ?C??FDFF>HHHA<I?GF#HJ#JHIG@IHGG#JGFJGGIJ#J<IIBIJIGJ.DICA>I?D5J#:#IH@CCI#))FCF<DCHD=E#C##8(?D##BDD5C@
            #
            # From mate 2 file:
            #
            # @chr1-24747/2
            # TAACAACAGACCTCTTAACAGAAACTCTACAAGCAAGACGAGATTGGGGGCCAATAATCTCAAAAAAAAAAAAAAAAGACTTTGTAACCCAGAATTTTAC
            # +
            # =?B+4FFFHEGCFI3#J##I#GFG3#2GAEH#FE#I0J#@CDBIJ?JIHI#GBHFG#GDI#IGF#;;#BE#C#D#ED@H#F93#B#?;4#FD>?9D>#D'
            #
            # Pair of SAM records:
            # chr1-24747      83      chr1    188739972       99      6=2X2=1X1=2X1=1X19=1X23=1X8=1X9=1X21=   =       188739775       -297    CTCCTGGTGATATAAGGAAAAGGAAAAACCTTTACGAGCCACTACAAAAACACACTGAAATAGACAGATCAATGACACAATGAAGCAACTACATTAACCA    @C5DDB##D?(8##C#E=DHCD<FCF))#ICC@HI#:#J5D?I>ACID.JGIJIBII<J#JIGGJFGJ#GGHI@GIHJ#JH#FG?I<AHHH>FFDF??C?
            # chr1-24747      163     chr1    188739775       99      15=1X1=1X7=1X8=1X3=1X17=1X3=1X18=1X3=1X1=1X13=1X        =       188739972       297  TAACAACAGACCTCTTAACAGAAACTCTACAAGCAAGACGAGATTGGGGGCCAATAATCTCAAAAAAAAAAAAAAAAGACTTTGTAACCCAGAATTTTAC    =?B+4FFFHEGCFI3#J##I#GFG3#2GAEH#FE#I0J#@CDBIJ?JIHI#GBHFG#GDI#IGF#;;#BE#C#D#ED@H#F93#B#?;4#FD>?9D>#D'
            #
            mate1 = l1.endswith('1')
            while True:
                sam1 = samfh.readline()
                if len(sam1) == 0:
                    raise RuntimeError('Ran out of SAM')
                if sam1[0] == '@':
                    continue # skip header
                break
            sam2 = samfh.readline()
            sam = sam1 if mate1 else sam2
            toks = sam.split('\t')
            nm, flag, refid, off, _, _, _, _, _, seq, qual = toks[:12]
            flag, off = int(flag), int(off)
            assert nm == l1[1:-2]
            fw = (flag & 16) == 0
            assert seq == (l2 if fw else revcomp(l2))
            refstr = ref.get(refid, off-1, len(seq))
            if not fw: refstr = revcomp(refstr)
            dist = editDistance(l2.upper(), refstr.upper())
            out = sys.stdout if dist < 12 else sys.stderr
            print >> out, 'rd: ' + l2
            print >> out, 'rf: ' + refstr
            print >> out, 'ed: ' + str(dist)
