import os
import logging
import sys
from subprocess import Popen
from aligner import Aligner

try:
    from Queue import Queue
except ImportError:
    from queue import Queue  # python 3.x


class BwaMem(Aligner):
    
    """ Encapsulates a BWA-MEM process.  The input can be a FASTQ
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
            raise RuntimeError('Must specify --index when aligner is bwa mem')
        options = []
        popen_stdin, popen_stdout, popen_stderr = None, None, None
        self.inQ, self.outQ = None, None
        # Compose input arguments
        if unpaired is not None and len(unpaired) > 1:
            raise RuntimeError('bwa mem can\'t handle more than one input file at a time')
        if paired is not None and len(paired) > 1:
            raise RuntimeError('bwa mem can\'t handle more than one input file at a time')
        if paired_combined is not None and len(paired_combined) > 1:
            raise RuntimeError('bwa mem can\'t handle more than one input file at a time')
        if unpaired is not None and (paired is not None or paired_combined is not None):
            raise RuntimeError('bwa mem can\'t handle unpaired and paired-end inputs at the same time')
        input_args = []
        if unpaired is not None:
            input_args = [unpaired[0]]
            input_args.extend(aligner_unpaired_args)
        if paired is not None:
            assert len(paired[0]) == 2
            input_args = [paired[0][0], paired[0][1]]
            input_args.extend(aligner_paired_args)
        if paired_combined is not None:
            options.append('-p')
            input_args = [paired_combined[0]]
            input_args.extend(aligner_paired_args)
        if unpaired is None and paired is None and paired_combined is None:
            raise RuntimeError("Must specify one or more of: unpaired, paired, paired_combined")
        # Compose output arguments
        output_args = []
        if sam is not None:
            output_args.extend(['>', sam])
        else:
            raise RuntimeError("Must specify SAM output")
        # Tell bwa mem whether to expected paired-end interleaved input
        if pairs_only:
            options.append('-p')
        # Put all the arguments together
        options.extend(aligner_args)
        cmd += ' '
        cmd += ' '.join(options + [index] + input_args + output_args)
        logging.info('bwa mem command: ' + cmd)
        if quiet:
            popen_stderr = open(os.devnull, 'w')
        self.pipe = Popen(cmd, shell=True,
                          stdin=popen_stdin, stdout=popen_stdout, stderr=popen_stderr,
                          bufsize=-1, close_fds='posix' in sys.builtin_module_names)

    def supports_mix(self):
        return False
