#!/usr/bin/env python
"""
Example:

$ python ./test_adios.py
"""

import adios as ad
import numpy as np
import matplotlib.pyplot as plt

## Writing
print "\n>>> Writing ...\n"

ad.init("config.xml")
fd = ad.open("temperature", "adios_test.bp", "w")

NX = 10
size = 2
groupsize =  4 + 4 + 8 * size * NX
t = np.array(range(NX*size), dtype=np.float64)
tt = t.reshape((size, NX))
ad.set_group_size(fd, groupsize)
ad.write_int(fd, "NX", NX)
ad.write_int(fd, "size", size)
ad.write(fd, "temperature", tt)
ad.close(fd)

ad.finalize()

## Reading
print "\n>>> Reading ...\n"

f = ad.file("adios_test.bp")
f.printself()

v = f.vars['temperature']
v.printself()

val = v.read()
print val
assert ((tt == val).all())
f.close()

## Testing
print "\n>>> Test utility functions ...\n"

print "bpls:\n", ad.bpls('adios_test.bp')
print "readvar:\n", ad.readvar("adios_test.bp", "temperature")

plt.imshow(val)
plt.colorbar()
plt.show()

print "\n>>> Done.\n"
