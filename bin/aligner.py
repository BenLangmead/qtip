"""
aligner.py

Abstract parent class for an aligner.
"""

from abc import ABCMeta, abstractmethod


class Aligner(object):
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def put(self, rd1, rd2=None):
        """ Only relevant if concurrency is supported. """
        pass
    
    @abstractmethod
    def done(self):
        pass
    
    @abstractmethod
    def supports_mix(self):
        """ Can take a mix if unpaired and paired-end reads as input? """
        pass

    @abstractmethod
    def supports_concurrency(self):
        """ Can take input reads on a queue and write output alignments
            to a queue?  Otherwise, . """
        pass

    @abstractmethod
    def writes_bam(self):
        """ Writes BAM directly to a file (like MOSAIK)?  Otherwise, we
            assume it writes SAM to stdout. """
        pass

    @abstractmethod
    def preferred_unpaired_format(self):
        """ Preferred input format for unpaired reads """
        pass

    @abstractmethod
    def preferred_paired_format(self):
        """ Preferred input format for paired-end reads, assuming reads
            can be specified in a single file. """
        pass
