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
    use adios_read_mod
    implicit none
    include 'mpif.h'

    character(len=25)   :: filename = "arrays.bp"
    integer             :: rank, size, i, j, ierr
    integer             :: comm

    ! ADIOS variables declarations for matching gwrite_temperature.fh 
    integer*8               :: f
    integer                 :: method = ADIOS_READ_METHOD_BP
    integer*8               :: sel

    ! variables to read in 
    integer                 :: NX, NY
    real*8, dimension(:,:), allocatable :: t
    integer, dimension(:), allocatable  :: p

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr)

    call adios_read_init_method (method, comm, "verbose=3", ierr);

    call adios_read_open (f, filename, method, comm, ADIOS_LOCKMODE_NONE, 1.0, ierr);

    ! select specific writer's data (using rank since each rank of the writer
    ! wrote one block of the data)
    call adios_selection_writeblock (sel, rank)

    ! First get the scalars to calculate the size of the arrays.
    ! Note that we cannot use adios_get_scalar here because that
    !   retrieves the same NX for everyone (from writer rank 0).
    call adios_schedule_read (f, sel, "NX", 0, 1, NX, ierr)
    call adios_schedule_read (f, sel, "NY", 0, 1, NY, ierr)
    call adios_perform_reads (f, ierr)
    write (*,'("rank=",i0," NX=",i0," NY=",i0)') rank, NX, NY

    ! Allocate space for the arrays
    allocate (t(NX,NY))
    allocate (p(NX))

    ! Read the arrays
    call adios_schedule_read (f, sel, "var_double_2Darray", 0, 1, t, ierr)
    call adios_schedule_read (f, sel, "var_int_1Darray", 0, 1, p, ierr)
    call adios_perform_reads (f, ierr)


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


    call adios_read_close (f, ierr)
    call MPI_Barrier (comm, ierr);
    call adios_read_finalize_method (method, ierr);
    call MPI_Finalize (ierr);

end program

