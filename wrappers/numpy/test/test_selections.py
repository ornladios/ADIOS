import unittest as ut
import adios as ad
import numpy as np
from adios._hl import selections as sel
from common import TempFile
from common import Slicee

class SimpleSelectionTestCase(ut.TestCase):
    def test_simpleselection(self):
        s = Slicee()
        shape = (10,20,30)
        x = sel.select(shape, s[5])
        self.assertEqual(x.shape, shape)
        self.assertEqual(x.mshape, (20,30))

        x = sel.select(shape, s[5,:])
        self.assertEqual(x.mshape, (20,30))

        x = sel.select(shape, s[5,...])
        self.assertEqual(x.mshape, (20,30))

        self.assertRaises(ValueError, sel.select, shape, s[20])
        self.assertRaises(ValueError, sel.select, shape, s[...,40])
        self.assertRaises(ValueError, sel.select, shape, s[...,30,20])

        shape = (10,)
        x = sel.select(shape, s[1])
        self.assertEqual(x.mshape, ())

        self.assertRaises(ValueError, sel.select, shape, s[10])

    def test_simpleselection_noninclusive(self):
        s = Slicee()
        shape = (10,20,30)
        x = sel.select(shape, s[:,:,15])
        self.assertEqual(x.mshape, (10,20))

        x = sel.select(shape, s[:,15,:])
        self.assertEqual(x.mshape, (10,30))

        x = sel.select(shape, s[:,15:16,:])
        self.assertEqual(x.mshape, (10,1,30))

        x = sel.select(shape, s[:,15:17,:])
        self.assertEqual(x.mshape, (10,2,30))

    def test_simpleselection_ellipsis(self):
        s = Slicee()
        shape = (10,20,30)
        x = sel.select(shape, s[:,1])
        self.assertEqual(x.mshape, (10,30))

        x = sel.select(shape, s[...,1])
        self.assertEqual(x.mshape, (10,20))

        x = sel.select(shape, s[...])
        self.assertEqual(x.mshape, (10,20,30))

        x = sel.select(shape, s[:,...])
        self.assertEqual(x.mshape, (10,20,30))

    def test_simpleselection_erros(self):
        s = Slicee()
        shape = (10,20,30)
        self.assertRaises(ValueError, sel.select, shape, s[...,...])

if __name__ == '__main__':
    ut.main()
