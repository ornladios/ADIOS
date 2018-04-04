#!/usr/bin/sed -f
##
## Example:
## sed -f adios_mpi2serial.sed adios_mpi.pyx > adios.pyx
##
s/^cdef extern from \"mpi-compat.h\": pass/#cdef extern from \"mpi-compat.h\": pass/g
s/^import mpi4py.MPI/#import mpi4py.MPI/g
s/^cimport mpi4py.MPI/#cimport mpi4py.MPI/g
s/MPI./MPI_/g
s/comm.ob_mpi/comm/g
s/comm.Clone()/comm/g
/ctypedef struct MPI_Comm:$/{
  N
  s/ctypedef struct MPI_Comm:\n[ ]*pass/ctypedef int MPI_Comm\
    int MPI_COMM_WORLD\
    int MPI_COMM_NULL\
    int MPI_COMM_SELF/
}
s/adios_mpi/adios/g
s/self.comm.Get_rank()/0/g
s/self.comm.bcast(ftmp, root=0)/ftmp/g
