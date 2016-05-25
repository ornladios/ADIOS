import unittest as ut
import adios as ad
import numpy as np
from common import TempFile
from common import Slicee

class AdiosTestCase(ut.TestCase):
    f    = None ## Adios File class
    temp = None ## TempFile

    def setUp(self):
        NX = 10
        size = 2
        t = np.array(range(NX*size), dtype=np.float64)
        tt = t.reshape((size, NX))

        self.temp = TempFile()

        ad.init_noxml()
        ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10);

        fw = ad.writer(self.temp.path)
        fw.declare_group('group', method='POSIX1')

        fw['NX'] = NX
        fw['size'] = size
        fw['temperature'] = tt
        fw.attrs['/temperature/description'] = "Global array written from 'size' processes"
        fw.attrs["/someSubGroup/anOtherGroup/anOtherAttribute"] = 99
        fw["/someSubGroup/anOtherGroup/anOtherVariable"] = 77
        fw.close()

        ad.finalize()

        self.f = ad.file(self.temp.path)

    def tearDown(self):
        try:
            if self.f:
                self.f.close()
        except:
            pass

    def test_adios_group(self):
        f = ad.file(self.temp.path)
        t = f["temperature"]
        # here t could have a member dictionary attr(s) again
        # which looks up all attributes starting with t.name

        # now match: f.attrs["/temperature/description"]
        self.assertEqual(t.attrs["description"], f.attrs["/temperature/description"])

        # the same should be possible for groups
        g = f["/someSubGroup/anOtherGroup"]

        # now match: f.attrs["/someSubGroup/anOtherGroup/anOtherAttribute"]
        self.assertEqual(g.attrs["anOtherAttribute"], f.attrs["/someSubGroup/anOtherGroup/anOtherAttribute"])

        # now match: f["/someSubGroup/anOtherGroup/anOtherVariable"]
        self.assertEqual(g["anOtherVariable"], f["/someSubGroup/anOtherGroup/anOtherVariable"])

    def test_adios_groupname(self):
        f = ad.file(self.temp.path)

        # Missing '/'
        g = f["someSubGroup/anOtherGroup"]

        # now match: f.attrs["/someSubGroup/anOtherGroup/anOtherAttribute"]
        self.assertEqual(g.attrs["anOtherAttribute"], f.attrs["/someSubGroup/anOtherGroup/anOtherAttribute"])

        # now match: f["/someSubGroup/anOtherGroup/anOtherVariable"]
        self.assertEqual(g["anOtherVariable"], f["/someSubGroup/anOtherGroup/anOtherVariable"])

if __name__ == '__main__':
    ut.main()
