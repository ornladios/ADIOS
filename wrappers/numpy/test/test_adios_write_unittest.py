import unittest as ut
import adios as ad
import numpy as np
from common import TempFile
from common import Slicee

class AdiosTestCase(ut.TestCase):

    def setUp(self):
        self.temp = TempFile()

    def tearDown(self):
        try:
            if self.temp:
                del self.temp
        except:
            pass

    def write_scalar(self, adtype, val, varname='val'):
        ad.init_noxml()
        ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10);
        g = ad.declare_group("group", "", ad.FLAG.YES)
        ad.define_var(g, varname, "", adtype, "", "", "")
        ad.select_method(g, "POSIX1", "", "")

        if adtype == ad.DATATYPE.string:
            dtype = ad.adios2npdtype(adtype, len(str(val)))
            npval = np.array(val, dtype=dtype)
        else:
            dtype = ad.adios2npdtype(adtype)
            npval = np.array(val, dtype=dtype)

        fd = ad.open("group", self.temp.path, "w")
        ad.set_group_size(fd, npval.nbytes)
        ad.write(fd, varname, val, dtype)
        ad.close(fd)
        ad.finalize()

        f = ad.file(self.temp.path)
        v = f.vars['val']
        self.assertEqual(v.read(), npval)

    def test_write_scalar_int(self):
        self.write_scalar(ad.DATATYPE.integer, 123)

    def test_write_scalar_double(self):
        self.write_scalar(ad.DATATYPE.double, 123)

    def test_write_scalar_string(self):
        self.write_scalar(ad.DATATYPE.string, 123)

if __name__ == '__main__':
    ut.main()
