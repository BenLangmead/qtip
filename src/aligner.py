"""
aligner.py

Encapsulates an aligner.  Classes for specific aligners inherit from
this class and override the constructor and these three member
functions.
"""

from abc import ABCMeta


class Aligner(object):
    __metaclass__ = ABCMeta

    @staticmethod
    def supports_mix():
        return False
