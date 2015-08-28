import os
import logging
import re
import string
import sys
from operator import itemgetter
from subprocess import Popen, PIPE
from read import Read, Alignment
from threading import Thread
from threadutil import fh2q, q2fh
from aligner import Aligner

try:
    from Queue import Queue
except ImportError:
    from queue import Queue  # python 3.x


class Bowtie2(Aligner):
    
    """ Encapsulates a Bowtie 2 process.  The input can be a FASTQ
        file, or a Queue onto which the caller enqueues reads.
        Similarly, output can be a SAM file, or a Queue from which the
        caller dequeues SAM records.  All records are textual; parsing
        is up to the user. """
    
    def __init__(self,
                 cmd,
                 aligner_args,
                 aligner_unpaired_args,
                 aligner_paired_args,
                 index,
                 unpaired=None,
                 paired=None,
                 paired_combined=None,
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
            raise RuntimeError('Must specify --index when aligner is Bowtie 2')
        cmd_toks = cmd.split()
        popen_stdin, popen_stdout, popen_stderr = None, None, None
        self.inQ, self.outQ = None, None
        # Make sure input arguments haven't been specified already
        for tok in ['-U', '-1', '-2']:
            assert tok not in cmd_toks
        # Compose input arguments
        self.input_is_queued = False
        self.output_is_queued = False
        input_args = []
        # Some of the Bowtie 2's format-related parameters take an
        # argument (e.g. -U, -1, -2, --tab5, --tab6) and some don't
        # (-f, -q, -r, -c)
        if input_format in ['fastq', 'fasta', 'raw']:
            if input_format == 'fastq':
                input_args.append('-q')
            elif input_format == 'fastq':
                input_args.append('-f')
            elif input_format == 'raw':
                input_args.append('-r')
            input_format = None
        if unpaired is not None:
            input_args.append(('--%s' % input_format) if input_format is not None else '-U')
            input_args.append(','.join(unpaired))
            input_args.extend(aligner_unpaired_args)
        if paired is not None:
            assert input_format not in ['tab5', 'tab6', '12']
            input_args.extend(['-1', ','.join(map(itemgetter(0), paired))])
            input_args.extend(['-2', ','.join(map(itemgetter(1), paired))])
            input_args.extend(aligner_paired_args)
        if paired_combined is not None:
            assert input_format is not None
            input_args.extend(['--%s' % input_format, ','.join(paired_combined)])
            input_args.extend(aligner_paired_args)
        if unpaired is None and paired is None and paired_combined is None:
            assert input_format is not None
            input_args.extend(['--%s' % input_format, '-'])
            popen_stdin = PIPE
            self.input_is_queued = True
        # Make sure output arguments haven't been specified already
        assert '-S' not in cmd_toks
        # Compose output arguments
        output_args = []
        if sam is not None:
            output_args.extend(['-S', sam])
        else:
            self.output_is_queued = True
            popen_stdout = PIPE
        index_args = ['-x', index]
        # Put all the arguments together
        input_args.extend(aligner_args)
        cmd += ' '
        cmd += ' '.join(input_args + output_args + index_args)
        logging.info('Bowtie 2 command: ' + cmd)
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
        return Read.to_tab6(rd1, rd2, truncate_name=truncate_name) + '\n'

    def put(self, rd1, rd2=None, truncate_name=True):
        assert self.input_is_queued
        self.inQ.put(self.format_read(rd1, rd2, truncate_name=truncate_name))

    def preferred_unpaired_format(self):
        return 'tab6'

    def preferred_paired_format(self):
        return 'tab6'

    def done(self):
        assert self.input_is_queued
        self.inQ.put(None)
    
    def supports_mix(self):
        return True

    def supports_concurrency(self):
        """ Can take input reads on a queue and write output alignments
            to a queue?  Otherwise, . """
        return True

    def writes_bam(self):
        """ Writes BAM directly to a file (like MOSAIK)?  Otherwise, we
            assume it writes SAM to stdout. """
        return False


class AlignmentBowtie2(Alignment):
    """ Encapsulates a Bowtie 2 SAM alignment record.  Parses certain
        important SAM extra fields output by Bowtie 2. """
    
    __mdRe = re.compile('MD:Z:([^\s]+)')  # MD:Z string

    def __init__(self):
        super(AlignmentBowtie2, self).__init__()
        self.name = None
        self.flags = None
        self.refid = None
        self.pos = None
        self.mapq = None
        self.cigar = None
        self.tlen = None
        self.seq = None
        self.qual = None
        self.extra = None
        self.fw = None
        self.mate1 = None
        self.mate2 = None
        self.paired = None
        self.aligned = None
        self.concordant = None
        self.discordant = None
        self.al_type = None
        self.bestScore = None
        self.ztzs = None
        self.mdz = None

    def parse(self, ln):
        """ Parse ln, which is a line of SAM output from Bowtie 2 """
        self.name, self.flags, self.refid, self.pos, self.mapq, self.cigar, \
            _, _, self.tlen, self.seq, self.qual, self.extra = \
            string.split(ln, '\t', 11)
        assert self.flags != "*"
        assert self.pos != "*"
        assert self.mapq != "*"
        assert self.tlen != "*"
        self.flags = flags = int(self.flags)
        self.pos = int(self.pos) - 1
        self.mapq = int(self.mapq)
        self.tlen = int(self.tlen)
        self.fw = (flags & 16) == 0
        self.mate1 = (flags & 64) != 0
        self.mate2 = (flags & 128) != 0
        self.paired = self.mate1 or self.mate2
        self.aligned = (self.flags & 4) == 0
        assert self.paired == ((flags & 1) != 0)
        self.concordant = ((flags & 2) != 0)
        self.discordant = ((flags & 2) == 0) and ((flags & 4) == 0) and ((flags & 8) == 0)
        # Parse MD:Z
        self.mdz = None
        se = self.__mdRe.search(self.extra)
        if se is not None:
            self.mdz = se.group(1)
        # Parse ZT.Z
        if self.aligned:
            ztzoff = self.extra.rfind('ZT:Z:')
            assert ztzoff != -1, ln
            self.ztzs = self.extra[ztzoff+5:].split(',')
            self.bestScore = int(self.ztzs[0])

    def rep_ok(self):
        # Call parent's repOk
        assert super(AlignmentBowtie2, self).rep_ok()
        return True
