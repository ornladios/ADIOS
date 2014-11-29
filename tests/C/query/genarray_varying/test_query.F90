
!  ADIOS is freely available under the terms of the BSD license described
!  in the COPYING file in the top level directory of this source distribution.
!
!  Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
!

!
! 
! Evaluate some queries on the output of genarray_varying
!
! (c) Oak Ridge National Laboratory, 2009
! Author: Norbert Podhorszki
!

program test_query
    use adios_query_mod
    implicit none
    include 'mpif.h'
    integer :: rank, nproc, comm
    integer :: ierr, j
    type(ADIOS_QUERY) :: q1, q2, q

    character(len=256) :: filename, errmsg
    integer*8           :: fh
    integer*8, dimension(2) :: offset=0, readsize=1
    integer*8           :: boxsel, pointsel  ! ADIOS selection objects

    integer             :: gndx,gndy
    integer             :: ts=0 ! timestep
    real*8, dimension(:,:),   allocatable :: xy
    real*8, dimension(:),     allocatable :: xy1D
    integer*8           :: batchsize = 10000000



    call MPI_Init (ierr)
    comm = MPI_COMM_WORLD
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, nproc , ierr)

    !print *,"call adios_init "
    call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=3", ierr)

    !write(filename,'(a6,i2.2,a3)') 'writer',ts,'.bp'
    write(filename,'("g.bp")') 
    call adios_read_open_file (fh, filename, ADIOS_READ_METHOD_BP, comm, ierr)

    if (fh .ne. 0) then
        call adios_get_scalar(fh, "/dimensions/gndx", gndx, ierr)
        call adios_get_scalar(fh, "/dimensions/gndy", gndy, ierr)

        readsize(1) = gndx/ nproc
        readsize(2) = gndy
        offset(1)   = rank * readsize(1)
        offset(2)   = 0

        if (rank == nproc-1) then  ! last process should read all the rest of columns
            readsize(1) = gndx - readsize(1)*(nproc-1)
        endif

        allocate( xy  (readsize(1), readsize(2)) )
        allocate( xy1D  (readsize(1)*readsize(2)) )
        xy1D = 0

        ! Create a 2D selection for the subset
        call adios_selection_boundingbox (boxsel, 2, offset, readsize)

        call adios_query_create (fh, boxsel, "xy", ADIOS_GT, "1.0", q)
        !call adios_query_create (fh, boxsel, "xy", ADIOS_GT, "1.0", q1)
        !call adios_query_create (fh, boxsel, "xy", ADIOS_LT, "2.0", q2)
        !call adios_query_combine (q1, ADIOS_QUERY_OP_AND, q2, q)

        call adios_query_evaluate (q, boxsel, ts, batchsize, pointsel, ierr)

        write (*, '("Evaluate return value = ",i0)') ierr 
        write (*, '("Selection pointer of the result = ",i0)') pointsel

        call adios_schedule_read (fh, pointsel, "xy", 0, 1, xy1D, ierr)
        call adios_perform_reads (fh, ierr)

        write (*, '("Points = ",100g12.6)') ( xy1D(j), j=1,readsize(1)*readsize(2))
        write (*, '("Number of points = ",i0)') 

        call adios_query_free (q)
        call adios_selection_delete (boxsel)
        call adios_selection_delete (pointsel)
        call adios_read_close (fh, ierr)

    endif


    ! Terminate
    call MPI_Barrier (comm, ierr)
    call adios_read_finalize_method (ADIOS_READ_METHOD_BP, ierr)
    !call MPI_Barrier (MPI_COMM_WORLD, ierr)
    !print *,"Writer calls MPI_Finalize"
    call MPI_Finalize (ierr)
    !print *,"Exit writer code "
end program test_query


