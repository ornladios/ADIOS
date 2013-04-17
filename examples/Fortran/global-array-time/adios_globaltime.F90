!
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

! ADIOS Fortran Example: write a global array from N processors with gwrite
! and write several timesteps into one BP file
!
! How to run: mpirun -np <N> adios_globaltime
! Output: adios_globaltime.bp
! ADIOS config file: adios_globaltime.xml
!

program adios_global
    use adios_write_mod
    implicit none
    include 'mpif.h'
    character(len=256)      :: filename = "adios_globaltime.bp"
    integer                 :: rank, size, i, it, ierr
    integer, parameter      :: NX = 10
    ! NY = 1 for testing purpose
    integer, parameter      :: NY = 1
    real*8                  :: t(NX)
    real*8                  :: p(NY)
    integer                 :: comm

    ! ADIOS variables declarations for matching gwrite_temperature.fh
    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    call adios_init ("adios_globaltime.xml", comm, adios_err)

    do it = 1, 13
        do i = 1, NX
            t(i)  = 100.0*it + NX*rank + i - 1
        enddo

        do i = 1, NY
            p(i)  = 1000.0*it + NY*rank + i - 1
        enddo

        ! We need to create the file in the first round,
        ! then we need to append to it
        if (it == 1) then
            call adios_open (adios_handle, "restart", filename, "w", comm, adios_err)
        else
            call adios_open (adios_handle, "restart", filename, "a", comm, adios_err)
        endif

#include "gwrite_restart.fh"

        call adios_close (adios_handle, adios_err)

        call MPI_Barrier (comm, ierr)
    enddo

    call MPI_Barrier (comm, ierr)

    call adios_finalize (rank, adios_err)

    call MPI_Finalize (ierr)
end program
