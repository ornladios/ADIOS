#!/usr/bin/sed -f
##
## Example:
## sed -f adios_mpi2serial.sed adios_mpi.pyx > adios.pyx
##
s/^import mpi4py.MPI/#import mpi4py.MPI/g
s/^cimport mpi4py.MPI/#cimport mpi4py.MPI/g
s/MPI./MPI_/g
s/comm.ob_mpi/comm/g
/ctypedef struct MPI_Comm:$/{
  N
  s/ctypedef struct MPI_Comm:\n[ ]*pass/ctypedef int MPI_Comm\
    int MPI_COMM_WORLD\
    int MPI_COMM_SELF/
}
