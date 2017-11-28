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
import sys
import operator
from subprocess import Popen, PIPE
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
                 aligner_args,
                 aligner_unpaired_args,
                 aligner_paired_args,
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

        args_single, args_paired = ['single', index], ['paired', index]

        if paired_combined is not None:
            assert paired is None
            compressed = paired_combined[0].endswith('.gz')
            args_paired.append('-pairedCompressedInterleavedFastq' if compressed else '-pairedInterleavedFastq')
            args_paired.extend(paired_combined)
        elif paired is not None:
            assert paired_combined is None
            compressed = paired[0][0].endswith('.gz')
            args_paired.append('-compressedFastq' if compressed else '-fastq')
            args_paired.extend(list(reduce(operator.concat, paired)))

        if unpaired is not None:
            compressed = unpaired[0].endswith('.gz')
            args_single.append('-compressedFastq' if compressed else '-fastq')
            args_single.extend(unpaired)

        if unpaired is None and paired is None and paired_combined is None:
            raise RuntimeError('Cannot instantiate SnapAligner without input file(s) specified')

        #
        # Compose output arguments
        #

        for tok in self._output_args:
            assert tok not in cmd_toks

        # Compose output arguments
        args_output = ['-o', '-bam']
        if sam is not None:
            args_output.append(sam)
        else:
            raise RuntimeError("Must specify SAM output")

        # Put all the arguments together
        cmd = ''
        if len(args_single) > 2:
            cmd += ' '.join(args_single)
            cmd += ' ' + ' '.join(args_output)
            if len(cmd_toks) > 1:
                cmd += ' ' + ' '.join(cmd_toks[1:])
            cmd += ' ' + ' '.join(aligner_unpaired_args)
            cmd += ' ' + ' '.join(aligner_args)

        if len(args_paired) > 2:
            if len(cmd) > 0:
                cmd += ' , '
            cmd += ' '.join(args_paired)
            cmd += ' ' + ' '.join(args_output)
            if len(cmd_toks) > 1:
                cmd += ' ' + ' '.join(cmd_toks[1:])
            cmd += ' ' + ' '.join(aligner_paired_args)
            cmd += ' ' + ' '.join(aligner_args)

        cmd = cmd_toks[0] + ' ' + cmd

        logging.info('SNAP command: ' + cmd)
        if quiet:
            popen_stderr = open(os.devnull, 'w')
        self.pipe = Popen(cmd, shell=True,
                          stdin=popen_stdin, stdout=popen_stdout, stderr=popen_stderr,
                          bufsize=-1, close_fds='posix' in sys.builtin_module_names)

    @staticmethod
    def supports_mix():
        """
        Note: return value of true doesn't just mean it can take some unpaired
        and some paired in a given invocation; it also means that can be
        interleaved in the input.
        """
        return False
