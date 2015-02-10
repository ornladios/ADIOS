program reader
    use adios_read_mod
    implicit none
    include 'mpif.h'

    character(len=256) :: filename, errmsg
    integer :: timesteps      ! number of times to read data    
    integer :: nproc          ! number of processors
    
    real*8, dimension(:,:),   allocatable :: xy

    ! Offsets and sizes
    integer :: nx_global, ny_global
    integer*8, dimension(2) :: offset=0, readsize=1
    integer*8           :: sel  ! ADIOS selection object

    ! MPI variables
    integer :: group_comm
    integer :: rank
    integer :: ierr

    integer :: ts   ! actual timestep
    integer :: i,j

    integer :: ntsteps
    
    integer*8 :: fh ! ADIOS file handle



    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)
    call adios_read_init_method (ADIOS_READ_METHOD_BP, group_comm, "verbose=3", ierr)

    ntsteps=2

    do ts = 0, ntsteps-1
       write(filename,'(a6,i2.2,a3)') 'writer',ts,'.bp'
       call adios_read_open_file (fh, filename, ADIOS_READ_METHOD_BP, group_comm, ierr)
  
       if (rank==0) write(6,*) 'ts=',ts
       if (ts==0) then
          ! adios_get_scalar() gets the value from metadata in memory
          call adios_get_scalar(fh, "nx_global", nx_global, ierr)
          call adios_get_scalar(fh, "ny_global", ny_global, ierr)
          readsize(1) = nx_global / nproc
          readsize(2) = ny_global

          offset(1)   = rank * readsize(1)
          offset(2)   = 0

          if (rank == nproc-1) then  ! last process should read all the rest of columns
             readsize(1) = nx_global - readsize(1)*(nproc-1)
          endif
          
          allocate( xy  (readsize(1), readsize(2)) )

          ! Create a 2D selection for the subset
          call adios_selection_boundingbox (sel, 2, offset, readsize)

       end if
       
       ! Arrays are read by scheduling one or more of them
       ! and performing the reads at once
       call adios_schedule_read(fh, sel, "xy", 0, 1, xy, ierr)
       call adios_perform_reads (fh, ierr)

       do j=1,readsize(2)
          do i=1,readsize(1)
             write (100+rank, '(3i5,f8.1)') ts,i-1+offset(1),j-1+offset(2),xy(i,j)
          enddo
       enddo
       call adios_read_close (fh, ierr)
    end do
    ! Terminate
    call adios_selection_delete (sel)
    call adios_read_finalize_method (ADIOS_READ_METHOD_BP, ierr)
    call MPI_Finalize (ierr)
  end program reader  
