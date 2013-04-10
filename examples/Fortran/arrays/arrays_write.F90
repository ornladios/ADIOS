!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!/*************************************************************/
!/*      Example of writing of simple arrays in ADIOS         */
!/*                                                           */
!/*     Similar example is ../../C/manual/2_adios_write.c     */
!/*************************************************************/
program arrays
    use adios_write_mod
    implicit none
    include 'mpif.h'

    character(len=25)   :: filename = "arrays.bp"
    integer             :: rank, size, i, j, ierr
    integer             :: comm

    ! ADIOS variables declarations for matching gwrite_temperature.fh 
    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle

    ! variables to write out 
    integer, parameter      :: NX = 10, NY = 100
    real*8                  :: t(NX,NY)
    integer                 :: p(NX)

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    do j = 1,NY 
        do i = 1,NX
            t(i,j) = rank*NY + (j-1) + i*(1.0/NX)
        enddo     
    enddo

    do i = 1,NX
        p(i) = rank*NX + i
    enddo

    call adios_init ("arrays.xml", comm, adios_err);
    call adios_open (adios_handle, "arrays", filename, "w", comm, adios_err);
#include "gwrite_arrays.fh"
    call adios_close (adios_handle, adios_err)

    call MPI_Barrier (comm, ierr);

    call adios_finalize (rank, adios_err);

    call MPI_Finalize (ierr);

end program

