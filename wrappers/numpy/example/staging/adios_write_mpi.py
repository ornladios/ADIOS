#!/usr/bin/env python
"""
Example:

$ python ./adios_write_mpi.py [METHOD] [params]
"""

import adios_mpi as ad
from mpi4py import MPI
import numpy as np
import getopt, sys
import os
import datetime

method = "POSIX1"
init = "verbose=3;"

if len(sys.argv) > 1:
    method = sys.argv[1]

if len(sys.argv) > 2:
    init = sys.argv[2]

## Init
ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10);
ad.init_noxml(MPI.COMM_WORLD)
g = ad.declare_group("temperature", "", 1)
ad.define_var(g, "NX", "", ad.DATATYPE.integer, "", "", "")
ad.define_var(g, "size", "", ad.DATATYPE.integer, "", "", "")
ad.define_var(g, "temperature", "", ad.DATATYPE.double, "size,NX", "size,NX", "0,0")
msg = str(datetime.datetime.now())
ad.define_attribute(g, "datetime", "", ad.DATATYPE.string, msg, "")
print ">>> Method:", method, init
ad.select_method(g, method, init, "")

## Writing
for i in range(5):
    print ">>> step:", i
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
