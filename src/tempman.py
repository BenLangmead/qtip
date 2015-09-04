import os
import tempfile
import errno
import shutil
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
        self.purged = False

    def get_file(self, fn_basename, group=''):
        """ Return filename for new temporary file in temp dir """
        assert not self.purged
        if fn_basename in self.files:
            raise RuntimeError('Temporary file with name "%s" already exists' % fn_basename)
        self.groups[group].append((fn_basename, False))
        self.files.add(fn_basename)
        return join(self.dir, fn_basename)

    def get_dir(self, dir_basename, group=''):
        """ Return filename for new temporary file in temp dir """
        assert not self.purged
        if dir_basename in self.dirs:
            raise RuntimeError('Temporary directory with name "%s" already exists' % dir_basename)
        self.groups[group].append((dir_basename, True))
        self.dirs.add(dir_basename)
        fullpath = join(self.dir, dir_basename)
        # Create output directory if needed
        try:
            os.makedirs(fullpath)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        return fullpath

    def remove_group(self, group):
        """ Remove all the temporary files belonging to the named group """
        assert not self.purged
        for base, is_dir in self.groups[group]:
            if is_dir:
                self.dirs.remove(base)
                shutil.rmtree(base)
            else:
                self.files.remove(base)
                os.remove(join(self.dir, base))
        del self.groups[group]

    def purge(self):
        """ Remove all temporary files created for the instantiator """
        shutil.rmtree(self.dir)
        self.purged = True

    def size(self):
        """ Return total size of all the files in the temp dir """
        assert not self.purged
        return _recursive_size(self.dir)

    def update_peak(self):
        """ Update peak size of temporary files """
        assert not self.purged
        self.peak_size = max(self.peak_size, self.size())
