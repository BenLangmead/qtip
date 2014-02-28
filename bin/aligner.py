'''
aligner.py

Encapsulates an aligner.
'''

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
    def supportsMix(self):
        pass
