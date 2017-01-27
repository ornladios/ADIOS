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

    def tearDown(self):
        ad.finalize()

    def test_writer_var(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(list(range(NX)), dtype=np.int32)
        val2 = np.array(list(range(5)), dtype='f8')

        fw = ad.writer(self.temp.path, method="POSIX1")

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
        val1 = np.array(list(range(NX)), dtype=np.int32)
        val2 = np.array(list(range(5)), dtype='f8')

        single_string = "ABCD"
        three_string = ("AA","BBB","CCCC")
        single_int = 10
        five_int = np.array(list(range(5)), dtype=np.int32)
        single_double = 1.1
        five_double = np.array(list(range(5)), dtype='double')*1.1
        unicode_string = u"unicode"
        bytes_string = u"bytes"

        fw = ad.writer(self.temp.path, method="POSIX1")

        fw.define_attr("single_string")
        fw.define_attr("three_string")
        fw.define_attr("single_int")
        fw.define_attr("five_int")
        fw.define_attr("single_double")
        fw.define_attr("five_double")
        fw.define_attr("unicode_string")
        fw.define_attr("bytes_string")

        fw['single_string'] = single_string
        fw['three_string'] = three_string
        fw['single_int'] = single_int
        fw['five_int'] = five_int
        fw['single_double'] = single_double
        fw['five_double'] = five_double
        fw['unicode_string'] = unicode_string
        fw['bytes_string'] = bytes_string
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['single_string'].value, single_string.encode())
        ##self.assertTrue((f['three_string'].value == three_string).all())
        self.assertTrue(f['three_string'].value[0], three_string[0].encode())
        self.assertTrue(f['three_string'].value[1], three_string[1].encode())
        self.assertTrue(f['three_string'].value[2], three_string[2].encode())
        self.assertEqual(f['single_int'].value, single_int)
        self.assertTrue((f['five_int'].value == five_int).all())
        self.assertEqual(f['single_double'].value, single_double)
        self.assertTrue((f['five_double'].value == five_double).all())
        self.assertTrue(f['unicode_string'].value, unicode_string.encode())
        self.assertTrue(f['bytes_string'].value, bytes_string)

    def test_writer_undefined_var(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(list(range(NX)), dtype=np.int32)
        val2 = np.array(list(range(5)), dtype='f8')

        fw = ad.writer(self.temp.path, method="POSIX1")

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
        val1 = np.array(list(range(NX)), dtype=np.int32)
        val2 = np.array(list(range(5)), dtype='f8')

        fw = ad.writer(self.temp.path, method="POSIX1")

        fw.vars['NX'] = NX
        fw.vars['val1'] = val1
        fw.vars['val2'] = val2
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['NX'][...], NX)
        self.assertTrue((f['val1'][:] == val1).all())
        self.assertTrue((f['val2'][:] == val2).all())

    def test_writer_default_group(self):
        self.temp = TempFile()
        fw = ad.writer(self.temp.path, method="POSIX1")
        fw.close()

        f = ad.file(self.temp.path)
        f.close()

    def test_writer_varname(self):
        self.temp = TempFile()

        fw = ad.writer(self.temp.path, method="POSIX1")

        NVARS = 99
        fw.vars['nvars'] = NVARS
        fw.vars['/nvars'] = NVARS
        fw.vars['_nvars'] = NVARS
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f.nvars, 3)
        self.assertEqual(f._nvars[...], NVARS)

    def test_writer_transform(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(list(range(NX)), dtype=np.int32)
        val2 = np.array(list(range(5)), dtype='f8')

        fw = ad.writer(self.temp.path, method="POSIX1")

        fw['NX'] = NX
        fw['val1'] = val1
        fw['val2'] = val2
        fw['val1'].transform = 'identity'
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['NX'][...], NX)
        self.assertTrue((f['val1'][:] == val1).all())
        self.assertTrue((f['val2'][:] == val2).all())

    def test_writer_timeaggregation(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(list(range(NX)), dtype=np.int32)
        val2 = np.array(list(range(5)), dtype='f8')

        fw = ad.writer(self.temp.path, method="POSIX1")
        fw.set_time_aggregation(3200)

        fw['NX'] = NX
        fw['val1'] = val1
        fw['val2'] = val2
        fw.close()
        ad.finalize()

        f = ad.file(self.temp.path)
        self.assertEqual(f['NX'][...], NX)
        self.assertTrue((f['val1'][:] == val1).all())
        self.assertTrue((f['val2'][:] == val2).all())

    def test_writer_empty_define_(self):
        self.temp = TempFile()

        NX = 10
        val1 = np.array(list(range(NX)), dtype=np.int32)
        val2 = np.array(list(range(5)), dtype='f8')

        fw = ad.writer(self.temp.path, method="POSIX1")

        fw.define_var("NX")
        fw.define_var("val1", "NX")
        fw.define_var("val2", val2.shape)
        fw.define_var("extra")

        fw['NX'] = NX
        fw['val1'] = val1
        fw['val2'] = val2
        fw.close()

        f = ad.file(self.temp.path)
        self.assertEqual(f['NX'][...], NX)
        self.assertTrue((f['val1'][:] == val1).all())
        self.assertTrue((f['val2'][:] == val2).all())

if __name__ == '__main__':
    ut.main()
