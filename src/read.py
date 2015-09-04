import re
from sam import Cigar
from abc import ABCMeta, abstractmethod


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
        if len(rdname) >= 255:
            raise RuntimeError('Read name too long for SAM spec: "%s"' % rdname)
        return cls(rdname, seq, qual)

    @classmethod
    def pair_from_simulator(cls, seq1, qual1, refid1, refoff1, fw1, sc1,
                            seq2, qual2, refid2, refoff2, fw2, sc2, training_nm):
        """ Construct appropriate read object (with appropriate name) given
            some simulated properties of the read. """
        rdname = "!!ts-sep!!".join(["!!ts!!",
                                    refid1, "+" if fw1 else "-", str(refoff1), str(sc1),
                                    refid2, "+" if fw2 else "-", str(refoff2), str(sc2), training_nm])
        if len(rdname) >= 255:
            raise RuntimeError('Read name too long for SAM spec: "%s"' % rdname)
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

    __metaclass__ = ABCMeta
    
    @abstractmethod
    def parse(self, ln):
        pass
    
    def __init__(self):
        self.cigar = None
        self.cigar_obj = None
        self.left_clip = self.right_clip = None

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

    def parse_cigar(self):
        """ Lazily parse self.cigar into self.cigar_obj """
        if self.cigar_obj is None:
            self.cigar_obj = Cigar(self.cigar)

    def soft_clipped_left(self):
        """ Return amt soft-clipped from LHS """
        if self.left_clip is not None:
            return self.left_clip
        elif self.cigar_obj is not None:
            self.left_clip = self.cigar_obj.cigar_list[0][1] if (self.cigar_obj.cigar_list[0][0] == 4) else 0
            return self.left_clip
        for i, c in enumerate(self.cigar):
            if not str.isdigit(c):
                if c != 'S':
                    self.left_clip = 0;
                    return 0
                self.left_clip = int(self.cigar[:i])
                return self.left_clip

    def soft_clipped_right(self):
        """ Return amt soft-clipped from RHS """
        if self.right_clip is not None:
            return self.right_clip
        elif self.cigar_obj is not None:
            self.right_clip = self.cigar_obj.cigar_list[-1][1] if (self.cigar_obj.cigar_list[-1][0] == 4) else 0
            return self.right_clip
        if self.cigar[-1] != 'S':
            self.right_clip = 0
            return 0
        self.parse_cigar()
        self.right_clip = self.cigar_obj.cigar_list[-1][1] if (self.cigar_obj.cigar_list[-1][0] == 4) else 0
        return self.right_clip

    def rep_ok(self):
        """ Check alignment for internal consistency """
        #assert self.paired or self.fragment_length() == 0
        assert not self.is_aligned() or self.bestScore is not None
        return True
