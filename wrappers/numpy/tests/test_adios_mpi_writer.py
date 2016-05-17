#!/usr/bin/env python
"""
Example:

$ mpiexec -n 4 python ./test_adios_mpi_writer.py
"""

import adios_mpi as ad
import numpy as np
from mpi4py import MPI

## Init
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

## Prepare
print "\n>>> Prepare ... (rank = %d)\n" % rank
fname = 'adios_test_mpi_writer.bp'
NX = 10
t = np.array(range(NX*size), dtype=np.float64) + rank*NX
gdim = (size, NX)
offset = (rank, 0)

print "\n>>> Writing ... (rank = %d)\n" % rank
ad.init_noxml()
ad.allocate_buffer (ad.BUFFER_ALLOC_WHEN.NOW, 10);

fw = ad.writer(fname, comm=comm)
fw.declare_group('group', method='MPI')
fw.define_var('temperature', ldim=(1,NX), gdim=gdim, offset=offset)

fw['NX'] = NX
fw['size'] = size
fw['temperature'] = t
fw.attrs['/temperature/description'] = "Global array written from 'size' processes"
fw.close()

## Reading
if rank == 0:
    print "\n>>> Reading ...\n"

    f = ad.file(fname, comm=MPI.COMM_SELF)
    for key, val in f.vars.iteritems():
        print key, '=', val.read()

    for key, val in f.attrs.iteritems():
        print key, '=', val.value

    print "\n>>> Done.\n"
