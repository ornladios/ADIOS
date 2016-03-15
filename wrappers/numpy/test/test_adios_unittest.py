import unittest as ut
import adios as ad
import numpy as np
from common import TempFile
from common import Slicee

class AdiosTestCase(ut.TestCase):
    f    = None ## Adios File class
    temp = None ## TempFile

    def setUp(self):
        self.temp = TempFile()

        ad.init_noxml()

        ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10);
        g = ad.declare_group("temperature", "", ad.FLAG.YES)
        ad.define_var(g, "NX", "", ad.DATATYPE.integer, "", "", "")
        ad.define_var(g, "size", "", ad.DATATYPE.integer, "", "", "")
        ad.define_var(g, "temperature", "", ad.DATATYPE.double, "size,NX", "size,NX", "0,0")
        self.msg = "this is a test"
        self.unit = "C"
        ## attr must be <varpath> + / + something without / (line 857, common_read.c)
        ad.define_attribute(g, "temperature/desc", "", ad.DATATYPE.string, self.msg, "")
        ad.define_attribute(g, "temperature/unit", "", ad.DATATYPE.string, self.unit, "")
        ad.select_method(g, "POSIX1", "verbose=3", "")

        fd = ad.open("temperature", self.temp.path, "w")
        self.NX = 10
        self.size = 2
        groupsize =  4 + 4 + 8 * self.size * self.NX
        t = np.array(range(self.NX * self.size), dtype=np.float64)
        self.tt = t.reshape((self.size, self.NX))
        ad.set_group_size(fd, groupsize)
        ad.write_int(fd, "NX", self.NX)
        ad.write_int(fd, "size", self.size)
        ad.write(fd, "temperature", self.tt)
        ad.close(fd)

        ad.finalize()

        self.f = ad.file(self.temp.path)

    def tearDown(self):
        try:
            if self.f:
                self.f.close()
        except:
            pass

    def test_adios_file(self):
        self.assertEqual(self.f.nattrs, 2)
        self.assertEqual(self.f.nvars, 3)
        self.assertEqual(self.f.current_step, 0)
        self.assertEqual(self.f.last_step, 0)
        self.assertEqual(sorted(self.f.vars.keys()), \
                         sorted(['NX', 'size', 'temperature']))
        self.assertEqual(sorted(self.f.attrs.keys()), \
                         sorted(['temperature/unit', 'temperature/desc']))

    def test_adios_attr(self):
        self.assertEqual(self.f.attrs['temperature/desc'].value, self.msg)
        self.assertEqual(self.f.attrs['temperature/desc'].dtype, np.dtype('S14'))

    def test_adios_file_getitem(self):
        self.assertRaises(TypeError, self.f.__getitem__, Slicee()[1])
        self.assertRaises(KeyError, self.f.__getitem__, Slicee()[:,:])
        self.assertRaises(KeyError, self.f.__getitem__, Slicee()['NONE'])

    def test_adios_var_scalar(self):
        v = self.f['NX']
        self.assertEqual(v.ndim, 0)
        self.assertEqual(v.dims, ())
        self.assertEqual(v.nsteps, 1)

        val = v.read()
        self.assertEqual(val, v[...])
        self.assertEqual(val, self.NX)

    def test_adios_var_array(self):
        v = self.f['temperature']
        self.assertEqual(v.ndim, 2)
        self.assertEqual(v.dims, (2L, 10L))
        self.assertEqual(v.nsteps, 1)

        val = v.read()
        self.assertEqual(val.dtype, np.dtype('float64'))
        self.assertEqual(val.ndim, 2)
        self.assertEqual(val.shape, (2, 10))
        self.assertTrue((val == v[:]).all())
        self.assertTrue((val == v[:,:]).all())
        self.assertTrue((v.read(offset=(0,5), count=(2,5)) == v[:,5:]).all())
        self.assertTrue((v.read(offset=(0,5), count=(2,5)) == v[:,-5:]).all())

    def test_adios_var_getitem(self):
        v = self.f['temperature']
        val = v.read()
        #self.assertTrue((val == v[...,...]).all())
        self.assertTrue((val == v[:,...]).all())
        self.assertTrue((val == v[:,::1]).all())

        ## equivalent to v[::2]
        ##self.assertRaises(IndexError, v.__getitem__, Slicee()[::2])
        #self.assertRaises(IndexError, v.__getitem__, Slicee()[:,:,:])

    def test_adios_var_array_squeeze(self):
        v = self.f['temperature']
        val = v[:,1]
        self.assertEqual(val.shape, (2,))

        val = v[:,1:2]
        self.assertEqual(val.shape, (2,1))

        val = v[1,:]
        self.assertEqual(val.shape, (10,))

        val = v[1:2,:]
        self.assertEqual(val.shape, (1,10))

    def test_adios_var_getattr(self):
        v = self.f['temperature']
        self.assertEqual(v.attrs['unit'].value, self.unit)

    def test_adios_var_read_points(self):
        v = self.f['temperature']
        x1 = ((0,0),)
        x2 = ((0,0),(0,1),)
        x3 = ((0,0),(0,1),(0,2),)

        #import ipdb; ipdb.set_trace()
        self.assertEqual(len(v.read_points()), 0)
        self.assertTrue((v.read_points(x1) == self.tt[0,0:1]).all())
        self.assertTrue((v.read_points(x2) == self.tt[0,0:2]).all())
        self.assertTrue((v.read_points(x3) == self.tt[0,0:3]).all())

        xerr = ((0,0),(0,1),(0,2,0),)
        self.assertRaises(IndexError, v.read_points, xerr)

    def test_adios_var_read_fancy(self):
        v = self.f['temperature']
        m = np.arange(10) % 3 != 0

        self.assertTrue((v[:,m] == self.tt[:,m]).all())
        self.assertTrue((v[1,m] == self.tt[1,m]).all())
        self.assertEqual(v[1,m].shape, (6,))
        self.assertTrue((v[:1,m] == self.tt[:1,m]).all())
        self.assertRaises(ValueError, v.__getitem__, Slicee()[3,m])

if __name__ == '__main__':
    ut.main()
