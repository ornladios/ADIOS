#!/usr/bin/env python
"""
Example:

$ python ./test_adios.py
"""

import adios as ad
import numpy as np

## Writing
print "\n>>> Writing ...\n"

ad.init("config.xml")

for i in range(10):
    if i == 0:
        fd = ad.open("temperature", "adios_test.bp", "w")
    else:
        fd = ad.open("temperature", "adios_test.bp", "a")

    NX = 10
    size = 2
    groupsize =  4 + 4 + 8 * size * NX
    t = np.array(range(NX*size), dtype=np.float64) + i*NX*size
    tt = t.reshape((size, NX))
    ad.set_group_size(fd, groupsize)
    ad.write_int(fd, "NX", NX)
    ad.write_int(fd, "size", size)
    ad.write(fd, "temperature", tt)
    ad.close(fd)

ad.finalize()

## Reading
print "\n>>> Reading step-by-step ...\n"

f = ad.file("adios_test.bp")
f.printself()

v = f.vars['temperature']
v.printself()

for i in range(10):
    val = v.read(from_steps=i, nsteps=1)
    print "step =", i
    print val

assert ((tt == val).all())

## Testing
print "\n>>> Test utility functions ...\n"

print "bpls:\n", ad.bpls('adios_test.bp')
print "readvar:\n", ad.readvar("adios_test.bp", "temperature")

print "\n>>> Done.\n"
