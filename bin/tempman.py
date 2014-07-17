__author__ = 'langmead'

import os
import tempfile
from collections import defaultdict


class TemporaryFileManager(object):

    def __init__(self, dr=None):
        self.dir = tempfile.mkdtemp(dir=dr)
        self.files = set()
        self.groups = defaultdict(list)
        self.peak_size = 0

    def get_filename(self, fn_basename, group=''):
        """ Return filename for new temporary file in temp dir """
        if fn_basename in self.files:
            raise RuntimeError('Temporary file with name "%s" already exists' % fn_basename)
        self.groups[group].append(fn_basename)
        self.files.add(fn_basename)
        return os.path.join(self.dir, fn_basename)

    def remove_group(self, group):
        """ Remove all the temporary files belonging to the named group """
        for fn_basename in self.groups[group]:
            self.files.remove(fn_basename)
            os.remove(os.path.join(self.dir, fn_basename))
        del self.groups[group]

    def size(self):
        """ Return total size of all the files in the temp dir """
        return sum(os.path.getsize(os.path.join(self.dir, f)) for f in self.files)

    def update_peak(self):
        """ Update peak size of temporary files """
        self.peak_size = max(self.peak_size, self.size())
