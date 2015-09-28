"""
aligner.py

Encapsulates an aligner.  Classes for specific aligners inherit from
this class and override the constructor and these three member
functions.
"""

from abc import ABCMeta, abstractmethod


class Aligner(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def supports_mix(self):
        """ Can take a mix if unpaired and paired-end reads as input? """
        pass
