"""
Copyright 2016, Ben Langmead <langmea@cs.jhu.edu>

Concrete subclass for Bowtie 2 aligner.
"""

import os
import logging
import sys
from operator import itemgetter
from subprocess import Popen
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
            paired = list(paired)  # because we traverse it multiple times
            input_args.extend(['-1', ','.join(map(itemgetter(0), paired))])
            input_args.extend(['-2', ','.join(map(itemgetter(1), paired))])
            input_args.extend(aligner_paired_args)
        if paired_combined is not None:
            assert input_format is not None
            input_args.extend(['--%s' % input_format, ','.join(paired_combined)])
            input_args.extend(aligner_paired_args)
        if unpaired is None and paired is None and paired_combined is None:
            raise RuntimeError("Must specify one or more of: unpaired, paired, paired_combined")
        # Make sure output arguments haven't been specified already
        assert '-S' not in cmd_toks
        # Compose output arguments
        output_args = []
        if sam is not None:
            output_args.extend(['| samtools view -Sb - >', sam])
        else:
            raise RuntimeError("Must specify SAM output")
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

    @staticmethod
    def supports_mix():
        return True
