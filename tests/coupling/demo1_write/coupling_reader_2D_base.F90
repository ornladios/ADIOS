program coupling
    implicit none
    include 'mpif.h'

    character(len=256) :: filename
    integer :: timesteps      ! number of times to read data    
    integer :: nproc          ! number of processors
    
    real*8, dimension(:,:),   allocatable :: xy

    ! Offsets and sizes
    integer :: nx_global, ny_global
    integer*4, dimension(2) :: offset, readsize
    
    ! MPI variables
    integer :: group_comm
    integer :: rank
    integer :: ierr
 
    integer :: ts   ! actual timestep
    integer :: i,j

    integer :: ntsteps

    ! This example only works if the number of readers = number of writers
    !ADIOS integer             :: gcnt, vcnt, acnt
    !ADIOS integer*8           :: fh, gh, read_bytes

    integer :: posx=4

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
    call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size (group_comm, nproc , ierr)

    ntsteps=2

    do ts = 0, ntsteps-1
       write(filename,'(a4,i2.2,a4,i2.2)') 'cpes',ts,'.bn.',rank
       open(100,file=filename,status='OLD',form='unformatted',action='read')
       !ADIOS

       !ADIOS       if (ts==0) then
       read(100) nx_global,ny_global
       !ADIOS
       read(100) readsize(1), readsize(2)       
       !ADIOS

       offset(1) = mod(rank, posx) * readsize(1)
       offset(2) = rank/posx * readsize(2)

       !ADIOS if (rank == nproc-1) then  ! last process should read all the rest of columns
       !ADIOS      readsize(1) = nx_global - readsize(1)*(nproc-1)
       !ADIOS endif
       if (ts==0) then
          allocate( xy  (readsize(1), readsize(2)) )
       end if

       read(100) xy
 
       do j=1,readsize(2)
          do i=1,readsize(1)
             write (200+rank, '(3i5,f8.1)') ts,i-1+offset(1),j-1+offset(2),xy(i,j)
          enddo
       enddo       
       close(100)
       !ADIOS
    end do
    ! Terminate
    call MPI_Finalize (ierr)
  end program coupling  
