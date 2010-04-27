!  
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!/**************************************************************/
!/*      Example of reading of simple arrays in ADIOS          */
!/*      which were written from the same number of processors */
!/*                                                            */
!/*     Similar example is ../../C/manual/2_adios_read.c       */
!/**************************************************************/
program arrays
    implicit none
    include 'mpif.h'

    character(len=25)   :: filename = "arrays.bp"
    integer             :: rank, size, i, j, ierr
    integer             :: comm

    ! ADIOS variables declarations for matching gwrite_temperature.fh 
    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle, adios_buf_size

    ! variables to read in 
    integer                 :: NX, NY
    real*8, dimension(:,:), allocatable :: t
    integer, dimension(:), allocatable  :: p

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    call adios_init ("arrays.xml", adios_err);
    call adios_open (adios_handle, "arrays", filename, "r", comm, adios_err);

    adios_groupsize = 0
    adios_totalsize = 0
    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

    ! First read in the scalars to calculate the size of the arrays
    adios_buf_size = 4
    call adios_read (adios_handle, "NX", NX, adios_buf_size, adios_err)
    adios_buf_size = 4
    call adios_read (adios_handle, "NY", NY, adios_buf_size, adios_err)
    
    call adios_close (adios_handle, adios_err)
    ! Note, we have to close to perform the reading of the variables above

    write (*,'("rank=",i0," NX=",i0," NY=",i0)') rank, NX, NY

    ! Allocate space for the arrays
    allocate (t(NX,NY))
    allocate (p(NX))

    ! Read the arrays
    call adios_open (adios_handle, "arrays", filename, "r", comm, adios_err);
    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)
    adios_buf_size = 8 * (NX) * (NY)
    call adios_read (adios_handle, "var_double_2Darray", t, adios_buf_size, adios_err)
    adios_buf_size = 4 * (NX)
    call adios_read (adios_handle, "var_int_1Darray", p, adios_buf_size, adios_err)
    call adios_close (adios_handle, adios_err)


    ! Print the results
    write (*,'("rank=",i0," p = (",i0,$)') rank, p(1)
    do i = 2,NX 
        write (*,'(", ",i0,$)') p(i)
    enddo
    write (*,'(")")')

    write (*,'("rank=",i0," t(:,5) = (",f6.2,$)') rank, t(1,5)
    do i = 2,NX 
        write (*,'(", ",f6.2,$)') t(i,5)
    enddo
    write (*,'(")")')


    call MPI_Barrier (comm, ierr);

    call adios_finalize (rank, adios_err);

    call MPI_Finalize (ierr);

end program

