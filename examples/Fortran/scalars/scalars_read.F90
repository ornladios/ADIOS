!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!/*************************************************************/
!/*   Example of reading various types of variable in ADIOS   */
!/*************************************************************/
program scalars_read
    implicit none
    include 'mpif.h'

    character(len=25)   :: filename = "scalars.bp"
    integer             :: rank, size, i, ierr
    integer             :: comm

    ! ADIOS variables declarations for matching gread_scalars.fh 
    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle, adios_buf_size

    ! scalar variables to write out (including a string)
    integer*1           :: v1 = 0
    integer*2           :: v2 = 0
    integer*4           :: v3 = 0
    integer*8           :: v4 = 0

    integer*1           :: v5 = 0
    integer*2           :: v6 = 0
    integer*4           :: v7 = 0
    integer*8           :: v8 = 0

    real*4              :: v9 = 0.0
    real*8              :: v10 = 0.0

    character(len=20)   :: v11 = "undefined"

    complex*8           :: v12 = (0.0, 0.0)
    complex*16          :: v13 = (0.0, 0.0)

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr);

    call adios_init ("scalars.xml", adios_err);
    call adios_open (adios_handle, "scalars", filename, "r", comm, adios_err);
#include "gread_scalars.fh"
    call adios_close (adios_handle, adios_err)

    call MPI_Barrier (comm, ierr);

    call adios_finalize (rank, adios_err);

    call MPI_Finalize (ierr);

    if (rank == 0) then
        write (*, '("int*1      v1  = ",i3)') v1
        write (*, '("int*2      v2  = ",i3)') v2
        write (*, '("int*4      v3  = ",i3)') v3
        write (*, '("int*8      v4  = ",i3)') v4

        write (*, '("int*1      v5  = ",i3)') v5
        write (*, '("int*2      v6  = ",i3)') v6
        write (*, '("int*4      v7  = ",i3)') v7
        write (*, '("int*8      v8  = ",i3)') v8

        write (*, '("real*4     v9  = ",f6.2)') v9
        write (*, '("real*8     v10 = ",f6.2)') v10

        write (*, '("string     v11 = ",a)') trim(v11)

        write (*, '("complex*8  v12 = (",f6.2,", ", f6.2,")")') v12
        write (*, '("complex*16 v13 = (",f6.2,", ", f6.2,")")') v13
    endif

end program

