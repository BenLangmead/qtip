"""
Copyright 2016, Ben Langmead <langmea@cs.jhu.edu>

TemporaryFileManager for maintaining and measuring a set of related
temporary files.
"""

import os
import tempfile
import errno
import shutil
import logging
from collections import defaultdict
from os.path import join, getsize


def _recursive_size(dr):
    tot = 0
    for root, dirs, files in os.walk(dr):
        tot += sum(getsize(join(root, name)) for name in files)
    return tot


class TemporaryFileManager(object):
    """
    Dishes out temporary files and directories, with ability to report total
    temporary-file footprint at any given point.
    """

    def __init__(self, dr=None):
        self.dir = tempfile.mkdtemp(dir=dr)
        self.files = set()
        self.dirs = set()
        self.groups = defaultdict(list)
        self.peak_size = 0

    def get_file(self, fn_basename, group=''):
        """ Return filename for new temporary file in temp dir """
        fullpath = join(self.dir, fn_basename)
        if fn_basename in self.files:
            #raise RuntimeError('Temporary file with name "%s" already exists' % fn_basename)
            return fullpath
        self.groups[group].append((fn_basename, False))
        self.files.add(fn_basename)
        return fullpath

    def get_dir(self, dir_basename, group=''):
        """ Return filename for new temporary subdir in temp dir """
        fullpath = join(self.dir, dir_basename)
        if dir_basename in self.dirs:
            #raise RuntimeError('Temporary directory with name "%s" already exists' % dir_basename)
            return fullpath
        if len(group) == 0:
            group = dir_basename
        self.groups[group].append((dir_basename, True))
        self.dirs.add(dir_basename)
        # Create output directory if needed
        try:
            os.makedirs(fullpath)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        return fullpath

    def remove_group(self, group):
        """ Remove all the temporary files belonging to the named group """
        self.update_peak()
        for base, is_dir in self.groups[group]:
            if is_dir:
                self.dirs.remove(base)
                shutil.rmtree(join(self.dir, base))
            else:
                self.files.remove(base)
                os.remove(join(self.dir, base))
        del self.groups[group]

    def purge(self, log=logging):
        """ Remove all temporary files created for caller """
        self.update_peak()
        for root, subdirs, files in os.walk(self.dir):
            for fn in files:
                log.warning("  still have file: %s/%s" % (root, fn))
                os.remove(join(root, fn))
            for dr in subdirs:
                log.warning("  still have subdir: %s/%s" % (root, dr))
                shutil.rmtree(join(root, dr))
        assert len(list(os.walk(self.dir))) == 1, str(list(os.walk(self.dir)))
        self.files = set()
        self.dirs = set()
        self.groups = defaultdict(list)

    def size(self):
        """ Return total size of all the files in the temp dir """
        return _recursive_size(self.dir)

    def update_peak(self):
        """ Update peak size of temporary files """
        self.peak_size = max(self.peak_size, self.size())
