#!/usr/bin/env python
"""
Example:

$ python ./test_adios.py
"""

import adios
import numpy as np

### Writing
adios.init("config.xml")
fd = adios.open("temperature", "adios_test.bp", "w")

NX = 10
size = 2
groupsize =  4 + 4 + 8 * size * NX
t = np.array(range(NX*size), dtype=np.float64)
tt = t.reshape((size, NX))
adios.set_group_size(fd, groupsize)
adios.write_int(fd, "NX", NX)
adios.write_int(fd, "size", size)
adios.write(fd, "temperature", tt)
adios.close(fd)

adios.finalize()

### Reading
v = adios.readvar("adios_test.bp", "temperature")

assert ((tt == v).all())

print "Done."
