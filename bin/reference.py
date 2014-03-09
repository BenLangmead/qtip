import os
import re
from abc import ABCMeta, abstractmethod
from cPickle import load, dump


def iter_fasta_chunks(fasta_filenames, chunk_size=500000):
    """
    Generate every substring from each of the named fasta files.
    Extracts substrings of length 'chunk_size'.
    """
    for fa in fasta_filenames:
        with open(fa) as fh:
            rtrim = 1
            short_name = None
            buf, buf_length, chunk_offset = [], 0, 0
            for ln in fh:
                if ln[0] == '>':
                    while buf_length > 0:
                        buf_string = ''.join(buf)
                        if short_name is not None:
                            yield short_name, chunk_offset, buf_string[:chunk_size]
                        buf = [buf_string[chunk_size:]]
                        buf_length = len(buf[0])
                        chunk_offset += chunk_size
                    short_name = ln[1:].split()[0]
                    assert len(short_name) > 0
                    buf, buf_length, chunk_offset = [], 0, 0
                else:
                    assert short_name is not None
                    if ln[-1] != '\n':
                        rtrim = 0
                    buf.append(ln[:-rtrim])
                    buf_length += (len(ln) - rtrim)
                    while buf_length >= chunk_size:
                        buf_string = ''.join(buf)
                        yield short_name, chunk_offset, buf_string[:chunk_size]
                        buf = [buf_string[chunk_size:]]
                        buf_length = len(buf[0])
                        chunk_offset += chunk_size
            while buf_length > 0:
                assert short_name is not None
                buf_string = ''.join(buf)
                yield short_name, chunk_offset, buf_string[:chunk_size]
                buf = [buf_string[chunk_size:]]
                buf_length = len(buf[0])
                chunk_offset += chunk_size


def iter_fasta_chunks_fast(fasta_filenames, chunk_size=500000):
    """
    Generate every substring from each of the named fasta files.
    Extracts substrings of length approximately 'chunk_size' (but
    might deviate).
    """
    import io
    buf, buf2 = bytearray(chunk_size), bytearray(chunk_size)
    for fa in fasta_filenames:
        with io.open(fa, 'rb') as fh:
            seq_name, seq_offset, i, end_i = None, 0, 0, 0
            while True:
                if i == end_i:
                    i = 0
                    if len(buf) != chunk_size:
                        buf = bytearray(chunk_size)
                    end_i = fh.readinto(buf)
                    if end_i == 0:
                        break
                if buf[i] == ord('>'):
                    ws_idx = buf.find('\n', i+1)
                    while ws_idx == -1:
                        # entire name line didn't make it into the buffer
                        end_i += fh.readinto(buf2)
                        buf.extend(buf2)
                        ws_idx = buf.find('\n', i+1)
                    seq_name = str(buf[i+1:ws_idx].split()[0])
                    seq_offset = 0
                    i = ws_idx+1
                # another sequence starts in the middle of the buffer?
                carat_off = buf.find('>', i)
                if carat_off == -1 or carat_off > end_i:
                    carat_off = end_i
                if carat_off > i:
                    # the split is the slowest part
                    seq = str(bytearray().join(buf[i:carat_off].split()))
                    assert seq_name is not None
                    yield seq_name, seq_offset, seq
                    seq_offset += len(seq)
                i = carat_off


# TODO: parse 2-bit files


class ReferenceOOB(Exception):
    """ Out of bounds exception for reference sequences """
    
    def __init__(self, value):
        self.value = value
    
    def __str__(self):
        return repr(self.value)


class Reference(object):
    """ Abstract base class for concrete subclasses implementing
        different ways of getting at substrings in collections of
        FASTA files. """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def hasName(self, refid):
        """ Return True iff our FASTA files have a sequence named
            'refid' """
        return False
    
    @abstractmethod
    def length(self, refid):
        """ Return the length of the sequence named 'refid' """
        return 0
    
    @abstractmethod
    def get(self, refid, start, ln):
        """ Return the length-'ln' substring beginning at offset
            'start' in the sequence named 'refid' """
        return ''
    
    @abstractmethod
    def names(self):
        """ Return iterator over names of sequences """
        return None


class ReferenceSimple(Reference):
    """ Loads a group of FASTA files into memory for easy random access. """

    def __init__(self, fafns):
        self.refs, self.lens = {}, {}
        abort = False
        for fafn in fafns:
            with open(fafn, 'r') as fafh:
                name = None
                for line in fafh:
                    line = line.rstrip()
                    ln = len(line)
                    if ln > 0 and line[0] == '>':
                        ind = line.find(" ")
                        if ind == -1:
                            ind = len(line)
                        line = line[1:ind]
                        name = line
                        self.refs[name] = []
                        self.lens[name] = 0
                    else:
                        assert name is not None
                        self.refs[name].append(line)
                        self.lens[name] += ln
            if abort:
                break
        for k in self.refs.iterkeys():
            self.refs[k] = ''.join(self.refs[k])

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

    def get(self, refid, pos, ln):
        """ Return the specified substring of the reference """
        assert refid in self.refs
        if pos + ln > self.lens[refid]:
            raise ReferenceOOB('"%s" has length %d; tried to get [%d, %d)' % (refid, self.lens[refid], pos, pos+ln))
        return self.refs[refid][pos:pos+ln]


class ReferencePicklable(Reference):
    """ Encapsulates a collection of FASTA files.  If a pickle filename is
        specified and the pickle file exists, we read from there.
        Otherwise we read from the FASTA files.  If a pickle file name
        is specified but doesn't exist at first, we create it at the
        end. """
    
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
        """ Return the specified substring of the reference """
        assert refid in self.refs
        if pos + ln > self.lens[refid]:
            raise ReferenceOOB('"%s" has length %d; tried to get [%d, %d)' % (refid, self.lens[refid], pos, pos+ln))
        return self.refs[refid][pos:pos+ln]


class ReferenceIndexed(Reference):
    """ Like Reference but uses .fai index files to avoid ever loading
        entire sequences into memory.  Use in Python 'with' block so
        that FASTA filehandles are closed appropriately. """
    
    __removeWs = re.compile(r'\s+')
    
    def __init__(self, fafns):
        self.fafhs = {}
        self.faidxs = {}
        self.chr2fh = {}
        self.offset = {}
        self.lens = {}
        self.chars_per_line = {}
        self.bytes_per_line = {}
        
        for fafn in fafns:
            self.fafhs[fafn] = fh = open(fafn, 'r')
            # Parse the index files
            with open(fafn + '.fai') as idxfh:
                for ln in idxfh:
                    toks = ln.rstrip().split()
                    if len(toks) == 0:
                        continue
                    assert len(toks) == 5
                    ref_id, ln, offset, chars_per_line, bytes_per_line = toks
                    self.chr2fh[ref_id] = fh
                    self.offset[ref_id] = int(offset)  # 0-based
                    self.lens[ref_id] = int(ln)
                    self.chars_per_line[ref_id] = int(chars_per_line)
                    self.bytes_per_line[ref_id] = int(bytes_per_line)
    
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
        """ Return the specified substring of the reference. """
        assert refid in self.offset
        if start + ln > self.lens[refid]:
            raise ReferenceOOB('"%s" has length %d; tried to get [%d, %d)' % (refid, self.lens[refid], start, start + ln))
        fh, offset, chars_per_line, bytes_per_line = \
            self.chr2fh[refid], self.offset[refid], \
            self.chars_per_line[refid], self.bytes_per_line[refid]
        assert bytes_per_line > chars_per_line
        byte_off = offset
        byte_off += (start // chars_per_line) * bytes_per_line
        into = start % chars_per_line
        byte_off += into
        fh.seek(byte_off)
        left = chars_per_line - into
        # Count the number of line breaks interrupting the rest of the
        # string we're trying to read
        if ln < left:
            return fh.read(ln)
        else:
            nbreaks = 1 + (ln - left) // chars_per_line
            res = fh.read(ln + nbreaks * (bytes_per_line - chars_per_line))
            res = re.sub(self.__removeWs, '', res)
        assert len(res) == ln, 'len(%s) != %d' % (res, ln)
        return res

if __name__ == "__main__":

    import sys
    import unittest
    import argparse
    from tempfile import mkdtemp
    from shutil import rmtree
    from collections import defaultdict

    parser = argparse.ArgumentParser()
    parser.add_argument('--test', action='store_const', const=True, default=False, help='Do unit tests')

    args = parser.parse_args()

    if args.test:
        import unittest

        class TestSerial(unittest.TestCase):

            def setUp(self):
                self.tmpdir = mkdtemp()
                self.fa_fn_1 = os.path.join(self.tmpdir, 'tmp1.fa')
                with open(self.fa_fn_1, 'w') as fh:
                    fh.write(""">short_name1 with some stuff after whitespace
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
A
>short_name2 with some stuff after whitespace
CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT
CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT
>short_name3 with some stuff after whitespace
CA
""")

            def tearDown(self):
                rmtree(self.tmpdir)

            def test_iter_fasta_chunks_1(self):
                chunks = [x for x in iter_fasta_chunks([self.fa_fn_1], chunk_size=15)]

                self.assertEqual(13, len(chunks))

                self.assertEqual(('short_name1',  0, 'ACGTACGTACGTACG'), chunks[0])
                self.assertEqual(('short_name1', 15, 'TACGTACGTACGTAC'), chunks[1])
                self.assertEqual(('short_name1', 30, 'GTACGTACGTACGTA'), chunks[2])
                self.assertEqual(('short_name1', 45, 'CGTACGTACGTACGT'), chunks[3])
                self.assertEqual(('short_name1', 60, 'ACGTACGTACGTACG'), chunks[4])
                self.assertEqual(('short_name1', 75, 'TACGTA'), chunks[5])

                self.assertEqual(('short_name2',  0, 'CAGTCAGTCAGTCAG'), chunks[6])
                self.assertEqual(('short_name2', 15, 'TCAGTCAGTCAGTCA'), chunks[7])
                self.assertEqual(('short_name2', 30, 'GTCAGTCAGTCAGTC'), chunks[8])
                self.assertEqual(('short_name2', 45, 'AGTCAGTCAGTCAGT'), chunks[9])
                self.assertEqual(('short_name2', 60, 'CAGTCAGTCAGTCAG'), chunks[10])
                self.assertEqual(('short_name2', 75, 'TCAGT'), chunks[11])

                self.assertEqual(('short_name3',  0, 'CA'), chunks[12])

            def test_iter_fasta_chunks_fast_1(self):
                for chunk_size in xrange(51, 200):
                    chunks = [x for x in iter_fasta_chunks_fast([self.fa_fn_1], chunk_size=chunk_size)]
                    seqs = defaultdict(str)
                    for chunk in chunks:
                        name, offset, seq = chunk
                        name, seq = str(name), str(seq)
                        seqs[name] += seq
                    self.assertEqual('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA', seqs['short_name1'])
                    self.assertEqual('CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT', seqs['short_name2'])
                    self.assertEqual('CA', seqs['short_name3'])

            def test_iter_fasta_chunks_2(self):
                if False:
                    # ultimately not how I decided to implement look_ahead
                    chunks = [x for x in iter_fasta_chunks([self.fa_fn_1], chunk_size=15, look_ahead=5)]

                    self.assertEqual(13, len(chunks))

                    self.assertEqual(('short_name1',  0, 'ACGTACGTACGTACGTACGT'), chunks[0])
                    self.assertEqual(('short_name1', 15, 'TACGTACGTACGTACGTACG'), chunks[1])
                    self.assertEqual(('short_name1', 30, 'GTACGTACGTACGTACGTAC'), chunks[2])
                    self.assertEqual(('short_name1', 45, 'CGTACGTACGTACGTACGTA'), chunks[3])
                    self.assertEqual(('short_name1', 60, 'ACGTACGTACGTACGTACGT'), chunks[4])
                    self.assertEqual(('short_name1', 75, 'TACGTA'), chunks[5])

                    self.assertEqual(('short_name2',  0, 'CAGTCAGTCAGTCAGTCAGT'), chunks[6])
                    self.assertEqual(('short_name2', 15, 'TCAGTCAGTCAGTCAGTCAG'), chunks[7])
                    self.assertEqual(('short_name2', 30, 'GTCAGTCAGTCAGTCAGTCA'), chunks[8])
                    self.assertEqual(('short_name2', 45, 'AGTCAGTCAGTCAGTCAGTC'), chunks[9])
                    self.assertEqual(('short_name2', 60, 'CAGTCAGTCAGTCAGTCAGT'), chunks[10])
                    self.assertEqual(('short_name2', 75, 'TCAGT'), chunks[11])

                    self.assertEqual(('short_name3',  0, 'CA'), chunks[12])

        unittest.main(argv=[sys.argv[0]])
