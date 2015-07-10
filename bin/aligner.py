"""
aligner.py

Encapsulates an aligner.  Classes for specific aligners inherit from
this class and override the constructor and these three member
functions.  Some points about these overrides:

1. The constructor

2. Put is a function that, given a read, will
"""

from abc import ABCMeta, abstractmethod

class Aligner(object):
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def put(self, rd1, rd2=None):
        pass
    
    @abstractmethod
    def done(self):
        pass
    
    @abstractmethod
    def supports_mix(self):
        pass
