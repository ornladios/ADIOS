!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

! ADIOS Fortran Example: write a global array from N processors with gwrite
!
! How to run: mpirun -np <N> global_array_write_F
! Output: global_array_F.bp
! ADIOS config file: global_array_F.xml
!

program global_array_F 
    use adios_write_mod
    implicit none
    include 'mpif.h'
    character(len=256)      :: filename = "global_array_F.bp"
    integer                 :: rank, size, i, ierr
    integer, parameter      :: NX = 10
    real*8, dimension(NX)   :: t
    integer                 :: comm

    ! ADIOS variables declarations for matching gwrite_temperature.fh 
    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    do i = 1, NX
        t(i)  = 10.0*rank+i-1
    enddo

    call adios_init ("global_array_F.xml", comm, adios_err)

    call adios_open (adios_handle, "temperature", filename, "w", comm, adios_err)
#include "gwrite_temperature.fh"
    call adios_close (adios_handle, adios_err)

    call MPI_Barrier (comm, ierr)

    call adios_finalize (rank, adios_err)

    call MPI_Finalize (ierr)
end program
