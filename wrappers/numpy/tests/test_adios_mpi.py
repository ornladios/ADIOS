#!/usr/bin/env python
"""
Example:

$ mpiexec -n 4 python ./test_adios_mpi.py
"""

import adios_mpi as ad
import numpy as np
from mpi4py import MPI
import sys

## Init
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

## Writing
print "\n>>> Writing ... (rank = %d)\n" % rank

config = "config_mpi.xml"
if len(sys.argv) > 1:
    config = sys.argv[1]

ad.init(config, comm)
fd = ad.open("temperature", "adios_test_mpi.bp", "w", comm)

NX = 10
groupsize =  4 + 4 + 4 + 8 * 1 * NX
t = np.array(range(NX), dtype=np.float64) + rank*NX
ad.set_group_size(fd, groupsize)
ad.write_int(fd, "NX", NX)
ad.write_int(fd, "rank", rank)
ad.write_int(fd, "size", size)
ad.write(fd, "temperature", t)
ad.close(fd)

ad.finalize()

## Reading
if rank == 0:
    print "\n>>> Reading ...\n"

    f = ad.file("adios_test_mpi.bp", comm=MPI.COMM_SELF)
    f.printself()

    v = f.vars['temperature']
    v.printself()

    val = v.read()
    print val
    assert (int(np.sum(val)) == (size*NX-1)*(size*NX)/2)
    f.close()

print "\n>>> Done.\n"

## Testing
if rank == 0:
    print "\n>>> Test utility functions ...\n"

    print "bpls:\n", ad.bpls('adios_test_mpi.bp')
    print "readvar:\n", ad.readvar("adios_test_mpi.bp", "temperature")

    print "\n>>> Done.\n"
