import os
import logging
import re
import string
import sys
from subprocess import Popen, PIPE
from read import Read, Alignment
from threading import Thread
from threadutil import fh2q, q2fh
from aligner import Aligner

try:
    from Queue import Queue
except ImportError:
    from queue import Queue  # python 3.x


class Mosaik(Aligner):
    
    """ Encapsulates a BWA-MEM process.  The input can be a FASTQ
        file, or a Queue onto which the caller enqueues reads.
        Similarly, output can be a SAM file, or a Queue from which the
        caller dequeues SAM records.  All records are textual; parsing
        is up to the user. """
    
    def __init__(self,
                 cmd,
                 index,
                 unpaired=None,
                 paired=None,
                 paired_combined=None,
                 pairsOnly=False,
                 sam=None,
                 quiet=False,
                 format=None):
        """ Create new process.
            
            Inputs:
            
            'unpaired' is an iterable over unpaired input filenames.
            'paired' is an iterable over pairs of paired-end input
            filenames.  If both are None, then input reads will be
            taken over the inQ.  If either are non-None, then a call
            to inQ will raise an exception.
            
            Outputs:
            
            'sam' is a filename where output SAM records will be
            stored.  If 'sam' is none, SAM records will be added to
            the outQ.
        """
        if paired is not None:
            raise RuntimeError('Mosaik does not accept separate mate1/mate2 files')
        if index is None:
            raise RuntimeError('Must specify --index when aligner is mosaik')
        options = []
        popen_stdin, popen_stdout, popen_stderr = None, None, None
        self.inQ, self.outQ = None, None
        # Compose input arguments
        if unpaired is not None and len(unpaired) > 1:
            raise RuntimeError('mosaik can\'t handle more than one input file at a time')
        if paired is not None and len(paired) > 1:
            raise RuntimeError('mosaik can\'t handle more than one input file at a time')
        if paired_combined is not None and len(paired_combined) > 1:
            raise RuntimeError('mosaik can\'t handle more than one input file at a time')
        if unpaired is not None and (paired is not None or paired_combined is not None):
            raise RuntimeError('bwa mem can\'t handle unpaired and paired-end inputs at the same time')
        self.input_is_queued = False
        self.output_is_queued = False
        input_args = []
        if unpaired is not None:
            input_args = ['-in', unpaired[0]]
        if paired_combined is not None:
            input_args = ['-in', paired_combined[0]]
        if unpaired is None and paired is None and paired_combined is None:
            input_args = ['-']
            popen_stdin = PIPE
            self.input_is_queued = True
        # Compose output arguments
        output_args = []
        if sam is not None:
            output_args.extend(['-out', sam])
        else:
            self.output_is_queued = True
            popen_stdout = PIPE
        # Tell bwa mem whether to expected paired-end interleaved input
        if pairsOnly:
            options.append('-p')
        # Put all the arguments together
        cmd += ' '
        cmd += ' '.join(options + ['-ia', index] + input_args + output_args)
        logging.info('MosaikAlign command: ' + cmd)
        if quiet:
            popen_stderr = open(os.devnull, 'w')
        self.pipe = Popen(cmd, shell=True,
                          stdin=popen_stdin, stdout=popen_stdout, stderr=popen_stderr,
                          bufsize=-1, close_fds='posix' in sys.builtin_module_names)
        # Create queue threads, if necessary
        timeout = 0.2
        if self.input_is_queued:
            self.inQ = Queue()
            self._inThread = Thread(target=q2fh, args=(self.inQ, self.pipe.stdin, timeout))
            self._inThread.daemon = True  # thread dies with the program
            self._inThread.start()
        if self.output_is_queued:
            self.outQ = Queue()
            self._outThread = Thread(target=fh2q, args=(self.pipe.stdout, self.outQ, timeout))
            self._outThread.daemon = True  # thread dies with the program
            self._outThread.start()

    @staticmethod
    def format_read(rd1, rd2=None, truncate_name=True):
        return Read.to_interleaved_fastq(rd1, rd2, truncate_name=truncate_name) + '\n'

    def put(self, rd1, rd2=None, truncate_name=True):
        assert self.input_is_queued
        self.inQ.put(self.format_read(rd1, rd2, truncate_name=truncate_name))

    @staticmethod
    def preferred_unpaired_format():
        return 'fastq'

    @staticmethod
    def preferred_paired_format():
        return 'interleaved_fastq'

    def done(self):
        assert self.input_is_queued
        self.inQ.put(None)
    
    def supports_mix(self):
        return False


class AlignmentMosaik(Alignment):
    """ Encapsulates a bwa mem SAM alignment record.  Parses certain
        important SAM extra fields output by bwa mem. """
    
    __asRe = re.compile('AS:i:([-]?[0-9]+)')  # best score
    __xsRe = re.compile('XS:i:([-]?[0-9]+)')  # second-best score
    __mdRe = re.compile('MD:Z:([^\s]+)')  # MD:Z string
    __zupRe = re.compile('ZP:i:([-]?[0-9]+)')  # best concordant
    __zlpRe = re.compile('Zp:i:([-]?[0-9]+)')  # 2nd best concordant
    __ztRe = re.compile('ZT:Z:([^\s]*)')  # extra features

    def __init__(self):
        super(AlignmentMosaik, self).__init__()
        self.name = None
        self.flags = None
        self.refid = None
        self.pos = None
        self.mapq = None
        self.cigar = None
        self.rnext = None
        self.pnext = None
        self.tlen = None
        self.seq = None
        self.qual = None
        self.extra = None
        self.fw = None
        self.mate1 = None
        self.mate2 = None
        self.paired = None
        self.concordant = None
        self.discordant = None
        self.bestScore = None
        self.secondBestScore = None
        self.bestConcordantScore = None
        self.secondBestConcordantScore = None
        self.ztzs = None
        self.mdz = None

    def parse(self, ln):
        """ Parse ln, which is a line of SAM output from bwa mem.  The line
            must correspond to an aligned read. """
        self.name, self.flags, self.refid, self.pos, self.mapq, self.cigar, \
            self.rnext, self.pnext, self.tlen, self.seq, self.qual, self.extra = \
            string.split(ln, '\t', 11)
        assert self.flags != "*"
        assert self.pos != "*"
        assert self.mapq != "*"
        assert self.tlen != "*"
        self.flags = flags = int(self.flags)
        self.pos = int(self.pos) - 1
        self.mapq = int(self.mapq)
        self.tlen = int(self.tlen)
        self.pnext = int(self.pnext)
        self.fw = (flags & 16) == 0
        self.mate1 = (flags & 64) != 0
        self.mate2 = (flags & 128) != 0
        self.paired = self.mate1 or self.mate2
        assert self.paired == ((flags & 1) != 0)
        self.concordant = ((flags & 2) != 0)
        self.discordant = ((flags & 2) == 0) and ((flags & 4) == 0) and ((flags & 8) == 0)
        # Parse AS:i
        se = self.__asRe.search(self.extra)
        self.bestScore = None
        if se is not None:
            self.bestScore = int(se.group(1))
        # Parse XS:i
        se = self.__xsRe.search(self.extra)
        self.secondBestScore = None
        if se is not None:
            self.secondBestScore = int(se.group(1))
        # Parse ZP:i
        se = self.__zupRe.search(self.extra)
        self.bestConcordantScore = None
        if se is not None:
            self.bestConcordantScore = int(se.group(1))
        # Parse Zp:i
        se = self.__zlpRe.search(self.extra)
        self.secondBestConcordantScore = None
        if se is not None:
            self.secondBestConcordantScore = int(se.group(1))
        # Parse ZT.Z
        self.ztzs = None
        se = self.__ztRe.search(self.extra)
        if se is not None:
            self.ztzs = se.group(1).split(',')
        # Parse MD:Z
        self.mdz = None
        se = self.__mdRe.search(self.extra)
        if se is not None:
            self.mdz = se.group(1)
        assert self.rep_ok()
    
    def rep_ok(self):
        # Call parent's repOk
        assert super(AlignmentMosaik, self).rep_ok()
        return True
