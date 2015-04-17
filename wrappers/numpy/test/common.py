import unittest as ut
from tempfile import mkstemp
import os

class Slicee(object):
    def __getitem__(self, index):
        return index

class TempFile(object):
    def __init__(self):
        self.fh, self.path = mkstemp(suffix='.bp', prefix='adios-test_')

    def __del__(self):
        try:
            if self.fh:
                os.close(self.fh)
                os.remove(self.path)
        except:
            pass

