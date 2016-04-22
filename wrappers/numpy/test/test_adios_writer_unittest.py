import unittest as ut
import adios as ad
import numpy as np
from common import TempFile
from common import Slicee

class AdiosTestCase(ut.TestCase):
    f    = None ## Adios File class
    temp = None ## TempFile

    def setUp(self):
        ad.init_noxml()
        ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10)

    def tearDown(self):
        ad.finalize()

    def test_writer_var(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(range(NX), dtype=np.int32)
        val2 = np.array(range(5), dtype='f8')

        fw = ad.writer(self.temp.path)
        fw.declare_group("group", method="POSIX1")

        fw.define_var("NX")
        fw.define_var("val1", "NX")
        fw.define_var("val2", val2.shape)

        fw['NX'] = NX
        fw['val1'] = val1
        fw['val2'] = val2
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['NX'][...], NX)
        self.assertTrue((f['val1'][:] == val1).all())
        self.assertTrue((f['val2'][:] == val2).all())

    def test_writer_attr(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(range(NX), dtype=np.int32)
        val2 = np.array(range(5), dtype='f8')

        single_string = "ABCD"
        three_string = ("AA","BBB","CCCC")
        single_int = 10
        five_int = np.array(range(5), dtype=np.int32)
        single_double = 1.1
        five_double = np.array(range(5), dtype='double')*1.1

        fw = ad.writer(self.temp.path)
        fw.declare_group("group", method="POSIX1")

        fw.define_attr("single_string")
        fw.define_attr("three_string")
        fw.define_attr("single_int")
        fw.define_attr("five_int")
        fw.define_attr("single_double")
        fw.define_attr("five_double")

        fw['single_string'] = single_string
        fw['three_string'] = three_string
        fw['single_int'] = single_int
        fw['five_int'] = five_int
        fw['single_double'] = single_double
        fw['five_double'] = five_double
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['single_string'].value, single_string)
        self.assertTrue((f['three_string'].value == three_string).all())
        self.assertEqual(f['single_int'].value, single_int)
        self.assertTrue((f['five_int'].value == five_int).all())
        self.assertEqual(f['single_double'].value, single_double)
        self.assertTrue((f['five_double'].value == five_double).all())

    def test_writer_undefined_var(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(range(NX), dtype=np.int32)
        val2 = np.array(range(5), dtype='f8')

        fw = ad.writer(self.temp.path)
        fw.declare_group("group", method="POSIX1")

        fw['NX'] = NX
        fw['val1'] = val1
        fw['val2'] = val2
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['NX'][...], NX)
        self.assertTrue((f['val1'][:] == val1).all())
        self.assertTrue((f['val2'][:] == val2).all())

    def test_writer_undefined_var2(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(range(NX), dtype=np.int32)
        val2 = np.array(range(5), dtype='f8')

        fw = ad.writer(self.temp.path)
        fw.declare_group("group", method="POSIX1")

        fw.vars['NX'] = NX
        fw.vars['val1'] = val1
        fw.vars['val2'] = val2
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['NX'][...], NX)
        self.assertTrue((f['val1'][:] == val1).all())
        self.assertTrue((f['val2'][:] == val2).all())

if __name__ == '__main__':
    ut.main()
