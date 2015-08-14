import re
import random
import logging
from sam import Cigar, Mdz, cigar_mdz_to_stacked, cigar_ref_to_stacked
from abc import ABCMeta, abstractmethod
from align import editDistance
from reference import ReferenceOOB


class Read(object):
    """ A read, along with some helper functions for parsing and
        composing a few differnt read formats. """
    
    def __init__(self, name, seq, qual, orig=None):
        """ Initialize new read given name, sequence, quality """
        self.name = name
        self.seq = seq
        self.qual = qual
        self.orig = orig
        assert self.rep_ok()
    
    @classmethod
    def from_simulator(cls, seq, qual, refid, refoff, fw, sc, training_nm):
        """ Construct appropriate read object (with appropriate name) given
            some simulated properties of the read. """
        rdname = "!!ts-sep!!".join(["!!ts!!", refid, "+" if fw else "-", str(refoff), str(sc), training_nm])
        return cls(rdname, seq, qual)

    @classmethod
    def pair_from_simulator(cls, seq1, qual1, refid1, refoff1, fw1, sc1,
                            seq2, qual2, refid2, refoff2, fw2, sc2, training_nm):
        """ Construct appropriate read object (with appropriate name) given
            some simulated properties of the read. """
        rdname = "!!ts-sep!!".join(["!!ts!!",
                                    refid1, "+" if fw1 else "-", str(refoff1), str(sc1),
                                    refid2, "+" if fw2 else "-", str(refoff2), str(sc2), training_nm])
        return cls(rdname, seq1, qual1), cls(rdname, seq2, qual2)

    @classmethod
    def to_tab6(cls, rd1, rd2=None, truncate_name=False):
        """ Convert either an unpaired read or a pair of reads to tab6
            format """
        name1 = rd1.name
        if truncate_name:
            name1 = name1.split()[0]
        if rd2 is not None:
            name2 = rd2.name
            if truncate_name:
                name2 = name2.split()[0]
            return "\t".join([name1, rd1.seq, rd1.qual, name2, rd2.seq, rd2.qual])
        return "\t".join([name1, rd1.seq, rd1.qual])
    
    @classmethod
    def from_tab6(cls, ln):
        toks = ln.rstrip().split('\t')
        if len(toks) == 3:
            return Read(toks[0], toks[1], toks[2]), None
        else:
            return Read(toks[0], toks[1], toks[2]), Read(toks[3], toks[4], toks[5])
    
    @classmethod
    def to_fastq(cls, rd, truncate_name=False):
        """ Convert a single read to FASTQ format """
        name = rd.name
        if truncate_name:
            name = name.split()[0]
        return '\n'.join(['@' + name, rd.seq, '+', rd.qual])
    
    @classmethod
    def to_interleaved_fastq(cls, rd1, rd2=None, truncate_name=False):
        """ Convert a single read to interleaved FASTQ format """
        name1 = rd1.name
        if truncate_name:
            name1 = name1.split()[0]
        if rd2 is not None:
            name2 = rd2.name
            if truncate_name:
                name2 = name2.split()[0]
            if not name1.endswith('/1'):
                name1 += '/1'
            if not name2.endswith('/2'):
                name2 += '/2'
            return '\n'.join(['@' + name1, rd1.seq, '+', rd1.qual, '@' + name2, rd2.seq, '+', rd2.qual])
        return '\n'.join(['@' + name1, rd1.seq, '+', rd1.qual])
    
    def __len__(self):
        """ Return number of nucleotides in read """
        return len(self.seq)
    
    def __str__(self):
        """ Return string representation """
        if self.orig is not None:
            return self.orig  # original string preferred
        elif self.qual is not None:
            assert len(self.seq) == len(self.qual)
            return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
        else:
            return ">%s\n%s\n" % (self.name, self.seq)
    
    def rep_ok(self):
        """ Check that read is internally consistent """
        if self.qual is not None:
            assert len(self.seq) == len(self.qual), "seq=%s\nqual=%s" % (self.seq, self.qual)
        return True


class Alignment(object):
    """ Abstract base class encapsulating an alignment record for a
        single aligned read.  Concrete subclasses consider the various
        tools' SAM dialects.
        
        Fields:
        
        refid: name of the reference sequence aligned to
        pos: 0-based reference offset of leftmost (w/r/t Watson) base
             involved in alignment (ignore soft-clipped bases)
        mapq: aligner-estimated mapping quality, -10 log10 (p) where p
              = probability alignment is incorrect
        """
    
    __nonAcgt = re.compile('[^ACGT]')
    __cigarLclip = re.compile('^([0-9]+)S.*')
    __cigarRclip = re.compile('.*[^0-9]([0-9]+)S$')
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def parse(self, ln):
        pass
    
    def __init__(self):
        pass

    def is_aligned(self):
        """ Return true iff read aligned """
        return (self.flags & 4) == 0
    
    def orientation(self):
        """ Return orientation as + or - """
        return '-' if ((self.flags & 16) != 0) else '+'
    
    def mate_mapped(self):
        """ Return true iff opposite mate aligned """
        return (self.flags & 8) == 0

    @staticmethod
    def fragment_length(al1, al2):
        """ Return fragment length """
        cigar1, cigar2 = Cigar(al1.cigar), Cigar(al2.cigar)
        if abs(al1.tlen) == 0:
            return 0  # no implied fragment
        al1_left = al1.pos - al1.soft_clipped_left()
        al1_rght = al1.pos + cigar1.reference_length() + al1.soft_clipped_right()
        al2_left = al2.pos - al2.soft_clipped_left()
        al2_rght = al2.pos + cigar2.reference_length() + al2.soft_clipped_right()
        return max(al1_rght, al2_rght) - min(al1_left, al2_left)

    def __len__(self):
        """ Return read length """
        return len(self.seq)
    
    def soft_clipped_left(self):
        """ Return amt soft-clipped from LHS """
        res = self.__cigarLclip.match(self.cigar)
        return 0 if res is None else int(res.group(1))

    def soft_clipped_right(self):
        """ Return amt soft-clipped from RHS """
        res = self.__cigarRclip.match(self.cigar)
        return 0 if res is None else int(res.group(1))

    def stacked_alignment(self, use_ref_for_edit_distance, ref=None):
        """ Return a stacked alignment corresponding to this
            alignment.  Optionally re-align the soft-clipped portions
            of the read so that the stacked alignment includes all
            characters from read. """
        cigar = Cigar(self.cigar)
        if self.mdz is not None:
            rd_aln, rf_aln = cigar_mdz_to_stacked(self.seq, cigar, Mdz(self.mdz))
        else:
            if ref is None:
                raise RuntimeError('Reference must be specified when MD:Z is absent')
            rd_aln, rf_aln = cigar_ref_to_stacked(self.seq, cigar, ref, self.refid, self.pos)
        rd_aln, rf_aln = [rd_aln], [rf_aln]
        cl = cigar.cigar_list
        if cl[0][0] == 4 or cl[-1][0] == 4:  # soft clipping?
            nrf = sum(map(lambda x: x == '-', rf_aln))
            for i in [0, -1]:
                if cl[i][0] == 4:  # soft clipping?
                    # Check for best alignment of soft-clipped portion
                    fudge = 10
                    unal_ln = cigar.cigar_list[i][1]
                    read_ln = len(self)
                    # 0-based offsets of leftmost and rightmost
                    # reference characters involved in alignment
                    posl_rf, posr_rf = self.pos, self.pos + nrf - 1
                    if i == 0:
                        readl, readr = 0, unal_ln
                        refl, refr = posl_rf - unal_ln, posl_rf
                    else:
                        readl, readr = read_ln - unal_ln, read_ln
                        refl, refr = posr_rf, posr_rf + unal_ln
                    if use_ref_for_edit_distance:
                        assert ref is not None
                        if refl < 0 or refr > ref.length(self.refid):
                            raise ReferenceOOB('[%d, %d) fell outside bounds for "%s": [0, %d)' %
                                               (refl, refr, self.refid, ref.length(self.refid)))
                    # Align read to ref using edit distance
                    rdstr = self.seq[readl:readr].upper()
                    rdstr = re.sub(self.__nonAcgt, 'N', rdstr)
                    assert refr - refl <= unal_ln + fudge
                    if use_ref_for_edit_distance:
                        logging.debug('GETTING BASES FROM REFERENCE')
                        refstr = ref.get(self.refid, refl, refr - refl).upper()
                        refstr = re.sub(self.__nonAcgt, 'N', refstr)
                    else:
                        logging.debug('GETTING RANDOM BASES')
                        refstr = ''.join([random.choice('ACGT') for _ in range(refr - refl)])
                    assert len(rdstr) == len(refstr)

                    eddist, stack = editDistance(rdstr, refstr, stacked=True)
                    rd_aln_new, rf_aln_new = stack
                    rd_aln_new, rf_aln_new = rd_aln_new.lower(), rf_aln_new.lower()
                    # Add traceback to stacked alignment
                    if i == 0:
                        rd_aln.insert(0, rd_aln_new)
                        rf_aln.insert(0, rf_aln_new)
                    else:
                        rd_aln.append(rd_aln_new)
                        rf_aln.append(rf_aln_new)


        rd_aln, rf_aln = ''.join(rd_aln), ''.join(rf_aln)
        rd_len = len(rd_aln) - rd_aln.count('-')
        rf_len = len(rf_aln) - rf_aln.count('-')
        return rd_aln, rf_aln, rd_len, rf_len
    
    def rep_ok(self):
        """ Check alignment for internal consistency """
        #assert self.paired or self.fragment_length() == 0
        assert not self.is_aligned() or self.bestScore is not None
        return True
