#!/usr/bin/env python
"""
Example:

$ python ./test_adios_writer.py
"""

""" Import ADIOS Python/Numpy wrapper """
import adios as ad
import numpy as np

print("\n>>> Prepare ...\n")
fname = 'adios_test_writer.bp'
NX = 10
size = 2
tt = np.arange(NX*size, dtype=np.float64).reshape((size, NX))

""" Writing """
print("\n>>> Writing ...\n")
ad.init_noxml()
ad.set_max_buffer_size (10)

fw = ad.writer(fname, method='POSIX1')

fw['NX'] = NX
fw['size'] = size
fw['temperature'] = tt
fw['temperature'].transform = "zfp:accuracy=0.001"
fw.attrs['/temperature/description'] = "Global array written from 'size' processes"
fw.close()

""" Reading """
print("\n>>> Reading ...\n")

f = ad.file(fname)
for key, val in f.vars.items():
    print(key, '=', val.read())

for key, val in f.attrs.items():
    print(key, '=', val.value)

""" Finalizing """
print("\n>>> Done.\n")
ad.finalize()
