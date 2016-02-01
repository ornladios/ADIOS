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
        ad.define_attribute(g, "desc", "", ad.DATATYPE.string, self.msg, "")
        ad.select_method(g, "POSIX1", "verbose=3", "")

        for i in range(5):
            mode = "a"
            if i == 0: mode = "w"
            fd = ad.open("temperature", self.temp.path, mode)
            self.NX = 10
            self.size = 2
            groupsize =  4 + 4 + 8 * self.size * self.NX
            t = np.array(range(self.NX * self.size), dtype=np.float64) + 100*i
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
        self.assertEqual(self.f.current_step, 0)
        self.assertEqual(self.f.last_step, 4)

    def test_adios_var_scalar(self):
        v = self.f['NX']
        self.assertEqual(v.ndim, 0)
        self.assertEqual(v.dims, ())
        self.assertEqual(v.nsteps, 5)

        val = v.read()
        self.assertEqual(val.dtype, np.dtype('int32'))
        self.assertEqual(val.ndim, 1)
        self.assertEqual(val.shape, (5,))

        self.assertTrue((val == v[:]).all())

    def test_adios_var_array(self):
        v = self.f['temperature']

        self.assertEqual(v.ndim, 2)
        self.assertEqual(v.dims, (2L, 10L))
        self.assertEqual(v.nsteps, 5)

        val = v.read()
        self.assertEqual(val.dtype, np.dtype('float64'))

        self.assertEqual(val.ndim, v.ndim + 1)
        self.assertEqual(val.shape, (5, self.size, self.NX))

        self.assertTrue((val == v[:]).all())
        self.assertTrue((val == v[:,:]).all())
        self.assertTrue((val == v[:,:,:]).all())
        #self.assertRaises(NotImplementedError, v.__getitem__, Slicee()[::2])
        self.assertRaises(TypeError, v.__getitem__, Slicee()[:,:,:,:])

        self.assertTrue((v.read(offset=(0,5), count=(2,5)) == v[:,:,5:]).all())
        self.assertTrue((v.read(offset=(0,5), count=(2,5)) == v[:,:,-5:]).all())

    def test_adios_var_array_squeeze(self):
        v = self.f['temperature']
        self.assertEqual(v[...].shape, (5,2,10))
        self.assertEqual(v[1,...].shape, (2,10))
        self.assertEqual(v[:,1,...].shape, (5,10))
        self.assertEqual(v[:,1].shape, (5,10))
        self.assertEqual(v[:,:,1].shape, (5,2))
        self.assertEqual(v[...,1].shape, (5,2))

        self.assertEqual(v[1:2,...].shape, (1,2,10))
        self.assertEqual(v[:,1:2,...].shape, (5,1,10))

    def test_adios_var_getitem(self):
        v = self.f['temperature']

        self.assertTrue((v[0,] == v.read(from_steps=0, nsteps=1)).all())
        self.assertTrue((v[1,] == v.read(from_steps=1, nsteps=1)).all())

        self.assertTrue((v[:2,] == v.read(from_steps=0, nsteps=2)).all())
        self.assertTrue((v[0:2,] == v.read(from_steps=0, nsteps=2)).all())

        self.assertTrue((v[:2,:1,:5] == v.read(offset=(0,0), count=(1,5), from_steps=0, nsteps=2)).all())

        self.assertTrue((v[-1,...] == v.read(from_steps=4, nsteps=1)).all())
        self.assertTrue((v[-2,...] == v.read(from_steps=3, nsteps=1)).all())

        #import ipdb; ipdb.set_trace()
        self.assertTrue((v[:,...,-1] == v.read(offset=(0,9), count=(2,1), scalar=(False,True))).all())
        self.assertTrue((v[:,...,-3:-1] == v.read(offset=(0,7), count=(2,2))).all())

        self.assertRaises(ValueError, v.__getitem__, Slicee()[:,:,-1:-2])

    def test_adios_var_read_points(self):
        v = self.f['temperature']
        x1 = ((0,0),)
        x2 = ((0,0),(0,1),)
        x3 = ((0,0),(0,1),(0,2),)

        self.assertEqual(len(v.read_points()), 0)
        ##self.assertTrue((v.read_points(x1) == self.tt[0,0:1]).all())
        ##self.assertTrue((v.read_points(x2) == self.tt[0,0:2]).all())
        ##self.assertTrue((v.read_points(x3) == self.tt[0,0:3]).all())

    def test_adios_var_read_fancy(self):
        v = self.f['temperature']
        m0 = np.arange(2) % 2 != 0
        m1 = np.arange(10) % 3 != 0

        self.assertEqual(v[:,:,m1].shape, (5,2,6))
        self.assertEqual(v[:,m0,m1].shape, (5,1,6))
        self.assertEqual(v[:,1,m1].shape, (5,6))
        self.assertEqual(v[1,1,m1].shape, (6,))

if __name__ == '__main__':
    ut.main()
