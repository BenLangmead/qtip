import os
import logging
import re
import tempfile
import string
import sys
from subprocess import Popen, PIPE
from read import Read, Alignment
from threading import Thread
from threadutil import fh2q, q2fh

try:
    from Queue import Queue
except ImportError:
    from queue import Queue  # python 3.x

class Bowtie2(object):
    
    ''' Encapsulates a Bowtie 2 process.  The input can be a FASTQ
        file, or a Queue onto which the caller enqueues reads.
        Similarly, output can be a SAM file, or a Queue from which the
        caller dequeues SAM records.  All records are textual; parsing
        is up to the user. '''
    
    def __init__(self, bt2_cmd, unpaired=None, paired=None, sam=None, quiet=False):
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
        cmdToks = bt2_cmd.split()
        popenStdin, popenStdout, popenStderr = None, None, None
        self.inQ, self.outQ = None, None
        # Make sure input arguments haven't been specified already
        for tok in ['-U', '-1', '-2']: assert tok not in cmdToks
        # Compose input arguments
        inputArgs = []
        if unpaired is not None:
            inputArgs.extend(['-U', ','.join(unpaired)])
        if paired is not None:
            inputArgs.extend(['-1', ','.join(map(lambda x: x[0], paired))])
            inputArgs.extend(['-2', ','.join(map(lambda x: x[1], paired))])
        if unpaired is None and paired is None:
            inputArgs.extend(['--tab6', '-'])
            popenStdin = PIPE
        # Make sure output arguments haven't been specified already
        assert '-S' not in cmdToks
        # Compose output arguments
        outputArgs = []
        if sam is not None:
            outputArgs.extend(['-S', sam])
        else:
            popenStdout = PIPE
        # Put all the arguments together
        bt2_cmd += ' '
        bt2_cmd += ' '.join(inputArgs + outputArgs)
        logging.info('Bowtie 2 command: ' + bt2_cmd)
        if quiet:
            popenStderr = open(os.devnull, 'w')
        ON_POSIX = 'posix' in sys.builtin_module_names
        # Create the concurrent subprocess
        # TODO: don't use shell=True?
        self.pipe = Popen(\
            bt2_cmd, shell=True,
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

class AlignmentBowtie2(Alignment):
    ''' Encapsulates a Bowtie 2 SAM alignment record.  Parses certain
        important SAM extra fields output by Bowtie 2. '''
    
    __asRe  = re.compile('AS:i:([-]?[0-9]+)') # best score
    __xsRe  = re.compile('XS:i:([-]?[0-9]+)') # second-best score
    __ytRe  = re.compile('YT:Z:([A-Z]+)')     # alignment type
    __ysRe  = re.compile('YS:i:([-]?[0-9]+)') # score of opposite
    __mdRe  = re.compile('MD:Z:([^\s]+)')     # MD:Z string
    __xlsRe = re.compile('Xs:i:([-]?[0-9]+)') # 3rd best
    __zupRe = re.compile('ZP:i:([-]?[0-9]+)') # best concordant
    __zlpRe = re.compile('Zp:i:([-]?[0-9]+)') # 2nd best concordant
    __kNRe  = re.compile('YN:i:([-]?[0-9]+)') # min valid score
    __knRe  = re.compile('Yn:i:([-]?[0-9]+)') # max valid score
    
    def parse(self, ln):
        ''' Parse ln, which is a line of SAM output from Bowtie 2.  The line
            must correspond to an aligned read. '''
        self.name, self.flags, self.refid, self.pos, self.mapq, self.cigar, \
        self.rnext, self.pnext, self.tlen, self.seq, self.qual, self.extra = \
        string.split(ln, '\t', 11)
        assert self.flags != "*"
        assert self.pos != "*"
        assert self.mapq != "*"
        assert self.tlen != "*"
        self.flags = int(self.flags)
        self.pos = int(self.pos) - 1
        self.mapq = int(self.mapq)
        self.tlen = int(self.tlen)
        self.pnext = int(self.pnext)
        self.fw = (self.flags & 16) == 0
        self.mate1 = (self.flags & 64) != 0
        self.mate2 = (self.flags & 128) != 0
        self.paired = self.mate1 or self.mate2
        self.exFields = {}
        for tok in string.split(self.extra, '\t'):
            idx1 = tok.find(':')
            assert idx1 != -1
            idx2 = tok.find(':', idx1+2)
            assert idx2 != -1
            self.exFields[tok[:idx2]] = tok[idx2+1:].rstrip()
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
        # Parse Xs:i
        se = self.__xlsRe.search(self.extra)
        self.thirdBestScore = None
        if se is not None:
            self.thirdBestScore = int(se.group(1))
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
        # Parse YT:Z
        se = self.__ytRe.search(self.extra)
        self.alType = None
        if se is not None:
            self.alType = se.group(1)
        # Parse YS:i
        se = self.__ysRe.search(self.extra)
        self.mateBest = None
        if se is not None:
            self.mateBest = int(se.group(1))
        # Parse MD:Z
        self.mdz = None
        se = self.__mdRe.search(self.extra)
        if se is not None:
            self.mdz = se.group(1)
        # Parse YN:i
        self.minValid = None
        se = self.__kNRe.search(self.extra)
        if se is not None:
            self.minValid = int(se.group(1))
        # Parse Yn:i
        self.maxValid = None
        se = self.__knRe.search(self.extra)
        if se is not None:
            self.maxValid = int(se.group(1))
        self.concordant = self.alType == "CP"
        self.discordant = self.alType == "DP"
        self.unpPair = self.alType == "UP"
        self.unp = self.alType == "UU"
        assert self.repOk()
    
    def repOk(self):
        # Call parent's repOk
        assert super(AlignmentBowtie2, self).repOk()
        assert self.alType in ["CP", "DP", "UP", "UU"]
        return True
