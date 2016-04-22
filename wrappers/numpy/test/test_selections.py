import unittest as ut
import adios as ad
import numpy as np
from adios._hl import selections as sel
from common import TempFile
from common import Slicee

class SimpleSelectionTestCase(ut.TestCase):
    def test_simpleselection(self):
        sgen = Slicee()
        shape = (10,20,30)
        x = sel.select(shape, sgen[5])
        self.assertEqual(x.shape, shape)
        self.assertEqual(x.mshape, (20,30))

        x = sel.select(shape, sgen[5,:])
        self.assertEqual(x.mshape, (20,30))

        x = sel.select(shape, sgen[5,...])
        self.assertEqual(x.mshape, (20,30))

        self.assertRaises(ValueError, sel.select, shape, sgen[20])
        self.assertRaises(ValueError, sel.select, shape, sgen[...,40])
        self.assertRaises(ValueError, sel.select, shape, sgen[...,30,20])

        shape = (10,)
        x = sel.select(shape, sgen[1])
        self.assertEqual(x.mshape, ())

        self.assertRaises(ValueError, sel.select, shape, sgen[10])

    def test_simpleselection_noninclusive(self):
        sgen = Slicee()
        shape = (10,20,30)
        x = sel.select(shape, sgen[:,:,15])
        self.assertEqual(x.mshape, (10,20))

        x = sel.select(shape, sgen[:,15,:])
        self.assertEqual(x.mshape, (10,30))

        x = sel.select(shape, sgen[:,15:16,:])
        self.assertEqual(x.mshape, (10,1,30))

        x = sel.select(shape, sgen[:,15:17,:])
        self.assertEqual(x.mshape, (10,2,30))

    def test_simpleselection_ellipsis(self):
        sgen = Slicee()
        shape = (10,20,30)
        x = sel.select(shape, sgen[:,1])
        self.assertEqual(x.mshape, (10,30))

        x = sel.select(shape, sgen[...,1])
        self.assertEqual(x.mshape, (10,20))

        x = sel.select(shape, sgen[...])
        self.assertEqual(x.mshape, (10,20,30))

        x = sel.select(shape, sgen[:,...])
        self.assertEqual(x.mshape, (10,20,30))

    def test_simpleselection_erros(self):
        sgen = Slicee()
        shape = (10,20,30)
        self.assertRaises(ValueError, sel.select, shape, sgen[...,...])

    def test_fancyselection(self):
        sgen = Slicee()
        shape = (10,20,30)
        m1 = np.arange(10) % 3 != 0
        m2 = np.arange(20) % 3 != 0

        x1 = sel.FancySelection(shape)
        x1[m1,:,:]
        self.assertEqual(x1.shape, shape)
        self.assertEqual(x1.mshape, np.zeros(shape)[m1,:,:].shape)

        x2 = sel.FancySelection(shape)
        self.assertRaises(ValueError, x2.__getitem__, sgen[m2,:])

    def test_unitfunction(self):
        self.assertEqual(sel._translate_indexlist([], 5), (slice(0, 5, None),))
        self.assertEqual(sel._translate_indexlist([1,2], 5), (slice(1, 3, None),))
        self.assertEqual(sel._translate_indexlist([0,1,3,4], 5), (slice(0, 2, None), slice(3, 5, None)))

        self.assertRaises(ValueError, sel._translate_indexlist, [5], 5)


if __name__ == '__main__':
    ut.main()
