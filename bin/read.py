from sam import cigarToList, mdzToList, cigarMdzToStacked
from abc import ABCMeta, abstractmethod

class Read(object):
    ''' A read, along with some helper functions for parsing and
        composing a few differnt read formats. '''
    
    def __init__(self, name, seq, qual, orig=None):
        ''' Initialize new read given name, sequence, quality '''
        self.name = name
        self.seq = seq
        self.qual = qual
        self.orig = orig
        assert self.repOk()
    
    @classmethod
    def fromSimulator(cls, seq, qual, refid, refoff, fw, sc, trainingNm):
        ''' Construct appropriate read object (with appropriate name) given
            some simulated properties of the read. '''
        rdname = "!!ts-sep!!".join(["!!ts!!", refid, "+" if fw else "-", str(refoff), str(sc), trainingNm])
        return cls(rdname, seq, qual)
    
    @classmethod
    def toTab6(cls, rd1, rd2=None, truncateName=False):
        ''' Convert either an unpaired read or a pair of reads to tab6
            format '''
        name1 = rd1.name
        if truncateName: name1 = name1.split()[0]
        if rd2 is not None:
            name2 = rd2.name
            if truncateName: name2 = name2.split()[0]
            return "\t".join([name1, rd1.seq, rd1.qual, name2, rd2.seq, rd2.qual])
        return "\t".join([name1, rd1.seq, rd1.qual])
    
    @classmethod
    def toFastq(cls, rd, truncateName=False):
        ''' Convert a single read to FASTQ format '''
        name = rd.name
        if truncateName: name = name.split()[0]
        return '\n'.join(['@' + name, rd.seq, '+', rd.qual])
    
    def __len__(self):
        ''' Return number of nucleotides in read '''
        return len(self.seq)
    
    def __str__(self):
        ''' Return string representation '''
        if self.orig is not None:
            return self.orig # original string preferred
        elif self.qual is not None:
            assert len(self.seq) == len(self.qual)
            return "@%s\n%s\n+\n%s\n" % (self.name, self.seq, self.qual)
        else:
            return ">%s\n%s\n" % (self.name, self.seq)
    
    def repOk(self):
        ''' Check that read is internally consistent '''
        if self.qual is not None:
            assert len(self.seq) == len(self.qual), "seq=%s\nqual=%s" % (self.seq, self.qual)
        return True

class Alignment(object):
    ''' Abstract base class encapsulating an alignment record for a
        single aligned read.  Concrete subclasses consider the various
        tools' SAM dialects.
        
        Fields:
        
        refid: name of the reference sequence aligned to
        pos: 0-based reference offset of leftmost (w/r/t Watson) base
             involved in alignment (ignore soft-clipped bases)
        mapq: aligner-estimated mapping quality, -10 log10 (p) where p
              = probability alignment is incorrect
        '''
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def parse(self, ln):
        pass
    
    def __init__(self, ln):
        self.parse(ln)
    
    def isAligned(self):
        ''' Return true iff read aligned '''
        return (self.flags & 4) == 0
    
    def orientation(self):
        ''' Return orientation as + or - '''
        return '-' if ((self.flags & 16) != 0) else '+'
    
    def mateMapped(self):
        ''' Return true iff opposite mate aligned '''
        return (self.flags & 8) == 0
    
    def fragmentLength(self):
        ''' Return fragment length '''
        return abs(self.tlen)
    
    def __len__(self):
        ''' Return read length '''
        return len(self.seq)
    
    def stackedAlignment(self, alignSoftClipped=False, ref=None):
        ''' Return a stacked alignment corresponding to this
            alignment.  Optionally re-align the soft-clipped portions
            of the read so that the stacked alignment includes all
            characters from read. '''
        cigarList = cigarToList(self.cigar)
        mdzList = mdzToList(self.mdz)
        # Get stacked alignment
        rdAln, rfAln = cigarMdzToStacked(self.seq, cigarList, mdzList)
        rdAln, rfAln = [rdAln], [rfAln]
        if alignSoftClipped:
            assert ref is not None
            nrd, nrf = 0, 0
            for i in xrange(len(rdAln)):
                if rdAln[i] != '-': nrd += 1
                if rfAln[i] != '-': nrf += 1
            for i in [0, -1]:
                if cigarList[i][0] == 4: # soft clipping?
                    # Check for best alignment of soft-clipped portion
                    fudge = 10
                    unalLn = cigarList[i][1]
                    readLn = len(self)
                    # 0-based offsets of leftmost and rightmost
                    # reference characters involved in alignment
                    posl_rf, posr_rf = self.pos, self.pos + nrf - 1
                    if i == 0:
                        readl, readr = 0, unalLn
                        refl, refr = posl_rf - unalLn, posl_rf
                    else:
                        readl, readr = readLn - unalLn, readLn
                        refl, refr = posr_rf, posr_rf + unalLn
                    if refl < 0:
                        diff = -refl
                        readl += diff
                        refl += diff
                    if refr > ref.length(self.refid):
                        diff = refr - ref.length(self.refid)
                        readr -= diff
                        refr -= diff
                    # Align read to ref using reasonable scoring
                    # function
                    rdstr = self.seq[readl:readr]
                    assert refr - refl <= unalLn + fudge
                    refstr = ref.get(self.refid, refl, refr - refl)
                    assert len(rdstr) == len(refstr), "\n%s\n%s\nrefl=%d,refr=%d,readl=%d,readr=%d" % (rdstr, refstr, refl, refr, readl, readr)
                    # Add traceback to stacked alignment
                    rdAlnNew, rfAlnNew = rdstr, refstr
                    if i == 0:
                        rdAln.insert(0, rdAlnNew)
                        rfAln.insert(0, rfAlnNew)
                    else:
                        rdAln.append(rdAlnNew)
                        rfAln.append(rfAlnNew)
        return ''.join(rdAln), ''.join(rfAln)
    
    def repOk(self):
        ''' Check alignment for internal consistency '''
        assert self.alType is not None
        assert self.paired or self.fragmentLength() == 0
        assert not self.isAligned() or self.bestScore is not None
        return True
