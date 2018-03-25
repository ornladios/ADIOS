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
print("\n>>> Writing ... (rank = %d)\n" % rank)

ad.init_noxml()

g = ad.declare_group("temperature", "", ad.FLAG.YES)
ad.define_var(g, "NX", "", ad.DATATYPE.integer, "", "", "")
ad.define_var(g, "size", "", ad.DATATYPE.integer, "", "", "")
ad.define_var(g, "offset", "", ad.DATATYPE.integer, "", "", "")
ad.define_var(g, "temperature", "", ad.DATATYPE.double, "1,NX", "size,NX", "offset,0")
ad.select_method(g, "MPI", "verbose=3", "")

fd = ad.open("temperature", "adios_test_mpi.bp", "w", comm=comm)
NX = 10
ad.write(fd, "NX", NX)
ad.write(fd, "size", size*2)

t = np.array(list(range(NX)), dtype=np.float64) + rank*NX
ad.write(fd, "offset", rank)
ad.write(fd, "temperature", t)

t = np.array(list(range(NX)), dtype=np.float64) + (rank+size)*NX
ad.write(fd, "offset", rank+size)
ad.write(fd, "temperature", t)

ad.close(fd)
ad.finalize()

## Reading
if rank == 0:
    print("\n>>> Reading ...\n")

    f = ad.file("adios_test_mpi.bp", comm=MPI.COMM_SELF)
    f.printself()

    v = f.vars['temperature']
    v.printself()

    val = v.read()
    print(val)
    assert (int(np.sum(val)) == (2*size*NX-1)*(2*size*NX)/2)
    f.close()

print("\n>>> Done.\n")
