program coupling
    implicit none
    include 'mpif.h'

    character(len=256) :: filename
    integer :: timesteps      ! number of times to read data    
    integer :: nproc          ! number of processors
    
    real*8, dimension(:,:),   allocatable :: xy

    ! Offsets and sizes
    integer :: nx_global, ny_global
    integer*8, dimension(2) :: offset=0, readsize=1

    ! MPI variables
    integer :: group_comm
    integer :: rank
    integer :: ierr

    integer :: ts   ! actual timestep
    integer :: i,j

    integer :: ntsteps
    
    ! This example can read from 1-260 readers
    integer             :: gcnt, vcnt, acnt
    integer*8           :: fh, gh, read_bytes



    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    ntsteps=2

    do ts = 0, ntsteps-1
       write(filename,'(a4,i2.2,a3)') 'cpes',ts,'.bp'
       call adios_fopen (fh, filename, group_comm, gcnt, ierr)
       call adios_gopen (fh, gh, "writer2D", vcnt, acnt, ierr)

       if (ts==0) then
          call adios_read_var(gh, "nx_global", offset, readsize, nx_global, read_bytes)
          call adios_read_var(gh, "ny_global", offset, readsize, ny_global, read_bytes)
          readsize(1) = nx_global / nproc
          readsize(2) = ny_global

          offset(1)   = rank * readsize(1)
          offset(2)   = 0

          if (rank == nproc-1) then  ! last process should read all the rest of columns
             readsize(1) = nx_global - readsize(1)*(nproc-1)
          endif
          
          allocate( xy  (readsize(1), readsize(2)) )
       end if
       
       call adios_read_var(gh, "xy", offset, readsize, xy, read_bytes)

       do j=1,readsize(2)
          do i=1,readsize(1)
             write (100+rank, '(3i5,f8.1)') ts,i-1+offset(1),j-1+offset(2),xy(i,j)
          enddo
       enddo
       call adios_gclose (gh, ierr)
       call adios_fclose (fh, ierr)
    end do
    ! Terminate
    call MPI_Finalize (ierr)
  end program coupling  
