"""
Encapsulates the SNAP aligner.

An important message from the fine folks at SNAP:

You may process more than one alignment without restarting SNAP, and if possible without reloading
the index.  In order to do this, list on the command line all of the parameters for the first
alignment, followed by a comma (separated by a space from the other parameters) followed by the
parameters for the next alignment (including single or paired).  You may have as many of these
as you please.  If two consecutive alignments use the same index, it will not be reloaded.
So, for example, you could do

'snap single hg19-20 foo.fq -o foo.sam , paired hg19-20 end1.fq end2.fq -o paired.sam'

and it would not reload the index between the single and paired alignments.

And another important message:

When specifying an input or output file, you can simply list the filename, in which case
SNAP will infer the type of the file from the file extension (.sam or .bam for example),
or you can explicitly specify the file type by preceeding the filename with one of the
 following type specifiers (which are case sensitive):
    -fastq
    -compressedFastq
    -sam
    -bam
    -pairedFastq
    -pairedInterleavedFastq
    -pairedCompressedInterleavedFastq

Input and output may also be from/to stdin/stdout. To do that, use a - for the input or output file
name and give an explicit type specifier.  So, for example,
snap single myIndex -fastq - -o -sam -
would read FASTQ from stdin and write SAM to stdout.

A couple of "paired" mode parameters to know about:

  -s   min and max spacing to allow between paired ends (default: 50 1000).
  -fs  force spacing to lie between min and max.
"""

import os
import logging
import re
import string
import sys
import operator
from subprocess import Popen, PIPE
from read import Read, Alignment
from threading import Thread
from threadutil import fh2q, q2fh
from aligner import Aligner

try:
    from Queue import Queue
except ImportError:
    from queue import Queue  # python 3.x


class SnapAligner(Aligner):
    
    """ Encapsulates a snap-aligner process.  The input can be a FASTQ
        file, or a Queue onto which the caller enqueues reads.
        Similarly, output can be a SAM file, or a Queue from which the
        caller dequeues SAM records.  All records are textual; parsing
        is up to the user. """

    _input_args = ['-fastq',
                   '-compressedFastq',
                   '-sam',
                   '-bam',
                   '-pairedFastq',
                   '-pairedInterleavedFastq',
                   '-pairedCompressedInterleavedFastq']

    _output_args = ['-o', '-sam', '-bam']

    def __init__(self,
                 cmd,
                 index,
                 unpaired=None,  # -fastq
                 paired=None,  # -pairedFastq
                 paired_combined=None,  # -pairedInterleavedFastq
                 pairs_only=False,
                 sam=None,
                 quiet=False,
                 input_format=None):
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

        if index is None:
            raise RuntimeError('Must specify --index when aligner is SNAP')

        cmd_toks = cmd.split()
        popen_stdin, popen_stdout, popen_stderr = None, None, None
        self.inQ, self.outQ = None, None

        #
        # Compose input arguments
        #

        for tok in self._input_args:
            assert tok not in cmd_toks

        self.input_is_queued = False
        self.output_is_queued = False

        args_single, args_paired = ['single'], ['paired']

        if paired_combined is not None:
            assert paired is None
            args_paired.append('-pairedInterleavedFastq')
            args_paired.extend(paired_combined)
        elif paired is not None:
            assert paired_combined is None
            args_paired.append('-pairedFastq')
            args_paired.extend(list(reduce(operator.concat, paired)))

        if unpaired is not None:
            args_single.append('-fastq')
            args_single.extend(unpaired)

        if unpaired is None and paired is None and paired_combined is None:
            raise RuntimeError('Cannot instantiate SnapAligner without input file(s) specified')

        #
        # Compose output arguments
        #

        for tok in self._output_args:
            assert tok not in cmd_toks

        # Compose output arguments
        args_output = ['-o', '-sam']
        if sam is not None:
            args_output.append(sam)
        else:
            args_output.append('-')
            self.output_is_queued = True
            popen_stdout = PIPE

        # Put all the arguments together
        if len(args_single) > 1:
            cmd += ' '.join([index] + args_single)
        if len(args_paired) > 1:
            if len(cmd) > 0:
                cmd += ' , '
            cmd += ' '.join([index] + args_paired)
        cmd += ' ' + ' '.join(args_output)
        cmd += ' ' + ' '.join(cmd_toks)

        logging.info('SNAP command: ' + cmd)
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
        """
        Note: return value of true doesn't just mean it can take some unpaired
        and some paired in a given invocation; it also means that can be
        interleaved in the input.
        """
        return False

    def supports_concurrency(self):
        """ Can take input reads on a queue and write output alignments
            to a queue?  Otherwise, . """
        return False  # might be True

    def writes_bam(self):
        """ Writes BAM directly to a file (like MOSAIK)?  Otherwise, we
            assume it writes SAM to stdout. """
        return False


class AlignmentSnap(Alignment):
    """ Encapsulates a SNAP SAM alignment record.  Parses certain
        important SAM extra fields output by snap-aligner. """

    __asRe = re.compile('AS:i:([-]?[0-9]+)')  # best score
    __xsRe = re.compile('XS:i:([-]?[0-9]+)')  # second-best score
    __mdRe = re.compile('MD:Z:([^\s]+)')  # MD:Z string
    __zupRe = re.compile('ZP:i:([-]?[0-9]+)')  # best concordant
    __zlpRe = re.compile('Zp:i:([-]?[0-9]+)')  # 2nd best concordant
    __ztRe = re.compile('ZT:Z:([^\s]*)')  # extra features

    def __init__(self):
        super(AlignmentSnap, self).__init__()
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
        """ Parse ln, which is a line of SAM output from SNAP.  The line
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
        assert super(AlignmentSnap, self).rep_ok()
        return True
