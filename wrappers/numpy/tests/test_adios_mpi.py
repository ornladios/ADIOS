#!/usr/bin/env python
"""
Example:

$ mpiexec -n 4 python ./test_adios_mpi.py
"""

import adios_mpi
import numpy as np
from mpi4py import MPI

## Init
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

## Writing
adios_mpi.init("config_mpi.xml", comm)
fd = adios_mpi.open("temperature", "adios_test_mpi.bp", "w", comm)

NX = 10
groupsize =  4 + 4 + 4 + 8 * 1 * NX
t = np.array(range(NX), dtype=np.float64) + rank*NX
adios_mpi.set_group_size(fd, groupsize)
adios_mpi.write_int(fd, "NX", NX)
adios_mpi.write_int(fd, "rank", rank)
adios_mpi.write_int(fd, "size", size)
adios_mpi.write(fd, "temperature", t)
adios_mpi.close(fd)

adios_mpi.finalize()

## Reading
v = adios_mpi.readvar("adios_test_mpi.bp", "temperature")

assert ((t == v[rank,]).all())

print "Done."
