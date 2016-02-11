#!/usr/bin/env python
"""
Example:

$ python ./test_adios.py
"""

import adios as ad
import numpy as np
import getopt, sys
import os
import datetime

method = "POSIX1"
init = "verbose=3;"

## Init
ad.init_noxml()
ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10);
g = ad.declare_group("temperature", "", 1)
ad.define_var(g, "NX", "", ad.DATATYPE.integer, "", "", "")
ad.define_var(g, "size", "", ad.DATATYPE.integer, "", "", "")
ad.define_var(g, "temperature", "", ad.DATATYPE.double, "size,NX", "size,NX", "0,0")
msg = str(datetime.datetime.now())
ad.define_attribute(g, "datetime", "", ad.DATATYPE.string, msg, "")
ad.define_attribute(g, "size2", "", ad.DATATYPE.integer, "0", "size")

ad.define_attribute(g, "temperature/unit", "", ad.DATATYPE.string, "cm", "")
ad.define_attribute(g, "temperature/unit/desc", "", ad.DATATYPE.string, "desc", "")
ad.define_attribute(g, "temperature/type", "", ad.DATATYPE.string, "None", "")

ad.select_method(g, method, init, "")
print ">>> Method:", method

## Writing
for i in range(5):
    fd = ad.open("temperature", "temp.bp", "a")

    NX = 10
    size = 2
    groupsize =  4 + 4 + 8 * size * NX
    t = np.array(range(NX*size), dtype=np.float64) + 100*i
    tt = t.reshape((size, NX))
    ad.set_group_size(fd, groupsize)
    ad.write_int(fd, "NX", NX)
    ad.write_int(fd, "size", size)
    ad.write(fd, "temperature", tt)
    ad.close(fd)

ad.finalize()
print ">>> Done."

f = ad.file('temp.bp')
v = f['temperature']
print(v.attrs['unit'].value)
print(f['temperature/unit'].value)

print(v.attrs['unit/desc'].value)
print(f['temperature/unit/desc'].value)

print(v.attrs['type'].value)
print(f['temperature/type'].value)
