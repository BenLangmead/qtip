import os
import re
from abc import ABCMeta, abstractmethod
from cPickle import load, dump

class ReferenceOOB(Exception):
    ''' Out of bounds exception for reference sequences '''
    
    def __init__(self, value):
        self.value = value
    
    def __str__(self):
        return repr(self.value)

class Reference(object):
    ''' Abstract base class for concrete subclasses implementing
        different ways of getting at substrings in collections of
        FASTA files. '''
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def hasName(self, refid):
        ''' Return True iff our FASTA files have a sequence named
            'refid' '''
        return False
    
    @abstractmethod
    def length(self, refid):
        ''' Return the length of the sequence named 'refid' '''
        return 0
    
    @abstractmethod
    def get(self, refid, start, ln):
        ''' Return the length-'ln' substring beginning at offset
            'start' in the sequence named 'refid' '''
        return ''
    
    @abstractmethod
    def names(self):
        ''' Return iterator over names of sequences '''
        return None

class ReferencePicklable(Reference):
    ''' Encapsulates a collection of FASTA files.  If a pickle filename is
        specified and the pickle file exists, we read from there.
        Otherwise we read from the FASTA files.  If a pickle file name
        is specified but doesn't exist at first, we create it at the
        end. '''
    
    def __init__(self, fafns, pickleFn=None, verbose=False):
        self.refs, self.lens = {}, {}
        pickleExists = False
        if pickleFn is not None:
            pickleExists = os.path.exists(pickleFn)
        if pickleFn is not None and pickleExists:
            self.load(pickleFn)
        else:
            last_pt = 0
            abort = False
            for fafn in fafns:
                with open(fafn, 'r') as fafh:
                    name = None
                    for line in fafh:
                        line = line.rstrip()
                        ln = len(line)
                        if ln > 0 and line[0] == '>':
                            ind = line.find(" ")
                            if ind == -1: ind = len(line)
                            line = line[1:ind]
                            name = line
                            self.refs[name] = []
                            self.lens[name] = 0
                        else:
                            assert name is not None
                            self.refs[name].append(line)
                            self.lens[name] += ln
                if abort: break
            for k in self.refs.iterkeys():
                self.refs[k] = ''.join(self.refs[k])
            if pickleFn is not None and not pickleExists:
                self.save(pickleFn)
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        pass
    
    def hasName(self, refid):
        return refid in self.refs
    
    def names(self):
        return self.refs.iterkeys()
    
    def length(self, refid):
        return self.lens[refid]
    
    def save(self, fn):
        #dump((self.refs, self.lens), open(fn, 'wb'), cPickle.HIGHEST_PROTOCOL)
        dump((self.refs, self.lens), open(fn, 'wb'), 2)
    
    def load(self, fn):
        self.refs, self.lens = load(open(fn, 'rb'))
    
    def get(self, refid, pos, ln):
        ''' Return the specified substring of the reference '''
        assert refid in self.refs
        if pos + ln > self.lens[refid]:
            raise ReferenceOOB('"%s" has length %d; tried to get [%d, %d)' % (refid, self.lens[refid], pos, pos+ln))
        return self.refs[refid][pos:pos+ln]

class ReferenceIndexed(Reference):
    ''' Like Reference but uses .fai index files to avoid ever loading
        entire sequences into memory.  Use in Python 'with' block so
        that FASTA filehandles are closed appropriately. '''
    
    __removeWs = re.compile(r'\s+')
    
    def __init__(self, fafns):
        self.fafhs = {}
        self.faidxs = {}
        self.chr2fh = {}
        self.offset = {}
        self.lens = {}
        self.charsPerLine = {}
        self.bytesPerLine = {}
        
        for fafn in fafns:
            self.fafhs[fafn] = fh = open(fafn, 'r')
            # Parse the index files
            with open(fafn + '.fai') as idxfh:
                for ln in idxfh:
                    toks = ln.rstrip().split()
                    if len(toks) == 0:
                        continue
                    assert len(toks) == 5
                    chr, ln, offset, charsPerLine, bytesPerLine = toks
                    self.chr2fh[chr] = fh
                    self.offset[chr] = int(offset) # 0-based
                    self.lens[chr] = int(ln)
                    self.charsPerLine[chr] = int(charsPerLine)
                    self.bytesPerLine[chr] = int(bytesPerLine)
    
    def __enter__(self):
        return self
    
    def __exit__(self, type, value, traceback):
        # Close all the open FASTA files
        for fafh in self.fafhs.itervalues():
            fafh.close()
    
    def hasName(self, refid):
        return refid in self.offset
    
    def names(self):
        return self.offset.iterkeys()
    
    def length(self, refid):
        return self.lens[refid]
    
    def get(self, refid, start, ln):
        ''' Return the specified substring of the reference. '''
        assert refid in self.offset
        if start + ln > self.lens[refid]:
            raise ReferenceOOB('"%s" has length %d; tried to get [%d, %d)' % (refid, self.lens[refid], start, start + ln))
        fh, offset, charsPerLine, bytesPerLine = \
            self.chr2fh[refid], self.offset[refid], \
            self.charsPerLine[refid], self.bytesPerLine[refid]
        assert bytesPerLine > charsPerLine
        byteOff = offset
        byteOff += (start // charsPerLine) * bytesPerLine
        into = start % charsPerLine
        byteOff += into
        fh.seek(byteOff)
        left = charsPerLine - into
        # Count the number of line breaks interrupting the rest of the
        # string we're trying to read
        if ln < left:
            return fh.read(ln)
        else:
            nbreaks = 1 + (ln - left) // charsPerLine
            res = fh.read(ln + nbreaks * (bytesPerLine - charsPerLine))
            res = re.sub(self.__removeWs, '', res)
        assert len(res) == ln, 'len(%s) != %d' % (res, ln)
        return res
