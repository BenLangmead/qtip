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

class BwaMem(Aligner):
    
    ''' Encapsulates a bwa mem process.  Note: I've disabled '''
    
    def __init__(self,
                 cmd,
                 index,
                 unpaired=None,
                 paired=None,
                 pairsOnly=False,
                 sam=None,
                 quiet=False):
        ''' Create new process.
            
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
        '''
        if index is None:
            raise RuntimeError('Must specify --index when aligner is bwa mem')
        cmdToks = cmd.split()
        options = []
        popenStdin, popenStdout, popenStderr = None, None, None
        self.inQ, self.outQ = None, None
        # Compose input arguments
        if unpaired is not None and len(unpaired) >= 1:
            raise RuntimeError('bwa mem can\'t handle more than one input file at a time')
        if paired is not None and len(paired) >= 1:
            raise RuntimeError('bwa mem can\'t handle more than one input file at a time')
        if unpaired is not None and paired is not None:
            raise RuntimeError('bwa mem can\'t handle unpaired and paired-end inputs at the same time')
        inputArgs = []
        if unpaired is not None:
            inputArgs = [unpaired[0]]
        if paired is not None:
            assert len(paired[0]) == 2
            inputArgs = [paired[0][0], paired[0][1]]
        if unpaired is None and paired is None:
            inputArgs = ['-']
            popenStdin = PIPE
        # Compose output arguments
        outputArgs = []
        if sam is not None:
            outputArgs.extend(['>', sam])
        else:
            popenStdout = PIPE
        # Tell bwa mem whether to expected paired-end interleaved input
        if pairsOnly:
            options.append('-p')
        # Put all the arguments together
        cmd += ' '
        cmd += ' '.join(options + [index] + inputArgs + outputArgs)
        logging.info('bwa mem command: ' + cmd)
        if quiet:
            popenStderr = open(os.devnull, 'w')
        ON_POSIX = 'posix' in sys.builtin_module_names
        self.pipe = Popen(\
            cmd, shell=True,
            stdin=popenStdin, stdout=popenStdout, stderr=popenStderr,
            bufsize=-1, close_fds=ON_POSIX)
        # Create queue threads, if necessary
        timeout = 0.2
        if unpaired is None and paired is None:
            self.inQ = Queue()
            self._inThread = Thread(target=q2fh, args=(self.inQ, self.pipe.stdin, timeout))
            self._inThread.daemon = True # thread dies with the program
            self._inThread.start()
        if sam is None:
            self.outQ = Queue()
            self._outThread = Thread(target=fh2q, args=(self.pipe.stdout, self.outQ, timeout))
            self._outThread.daemon = True # thread dies with the program
            self._outThread.start()
    
    def put(self, rd1, rd2=None):
        self.inQ.put(Read.to_interleaved_fastq(rd1, rd2, truncate_name=True) + '\n')
    
    def done(self):
        self.inQ.put(None)
    
    def supportsMix(self):
        return False


class AlignmentBwaMem(Alignment):
    """ Encapsulates a bwa mem SAM alignment record.  Parses certain
        important SAM extra fields output by bwa mem. """
    
    __asRe = re.compile('AS:i:([-]?[0-9]+)')  # best score
    __xsRe = re.compile('XS:i:([-]?[0-9]+)')  # second-best score
    
    def __init__(self):
        super(AlignmentBwaMem, self).__init__()
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
        # Parse MD:Z
        self.mdz = None
        assert self.rep_ok()
    
    def rep_ok(self):
        # Call parent's repOk
        assert super(AlignmentBwaMem, self).rep_ok()
        return True
