program writer
  use adios_write_mod
  implicit none
  include 'mpif.h'
  character(len=256) :: filename
  integer :: npx, npy   ! # of processors in x-y direction
  integer :: ndx, ndy   ! size of array per processor
  integer :: timesteps=2  ! number of timesteps to run  
  integer :: nproc      ! number of processors
  
  real*8, dimension(:,:),   allocatable :: xy
  ! Offsets and sizes
  integer :: offs_x, offs_y
  integer :: nx_local, ny_local
  integer :: nx_global, ny_global
  integer :: posx, posy   ! position index in the array
  integer :: i,j
  ! MPI variables
  integer :: group_comm
  integer :: rank
  integer :: ierr
  ! ADIOS variables
  character (len=200) :: group
  integer*8 :: adios_handle, adios_totalsize, adios_groupsize, adios_buf_size
  integer*8 :: m_adios_group, id1_xy, id2_xy
  character(len=32)       :: local, global, offset
  ! actual timestep
  integer   :: ts 
  
  call MPI_Init (ierr)
  call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
  call MPI_Comm_rank (group_comm, rank, ierr)
  call MPI_Comm_size (group_comm, nproc , ierr)

  ! call adios_init ("writer.xml", group_comm, ierr) -->
  call adios_init_noxml (group_comm, ierr)
  call adios_allocate_buffer (20, ierr)
  call adios_declare_group (m_adios_group, "writer", "iter", 1, ierr)
  call adios_select_method (m_adios_group, "MPI", "", "", ierr)

  ! will work with 12 cores, which are arranged by npx=4, npy=3 (4x3)
  npx = 4
  npy = 3
  ndx = 4 !65
  ndy = 3 !129
  ! will work with each core writing ndx = 65, ndy = 129, (65*4,129*3) global
  if (nproc.ne.npx*npy) then
     if (rank==0) write(6,*) 'ERROR with nproc=',nproc,'npx=',npx,'npy=',npy
     call adios_finalize (rank, ierr)
     call MPI_Finalize (ierr)
     call exit(1)
  end if
  !!  2D array with block,block decomposition
  posx = mod(rank, npx)     ! 1st dim easy: 0, npx, 2npx... are in the same X position
  posy = rank/npx           ! 2nd dim: npx processes belong into one dim
  offs_x = posx * ndx
  offs_y = posy * ndy
  nx_local   = ndx
  ny_local   = ndy
  nx_global  = npx * ndx * 2
  ny_global  = npy * ndy
  allocate( xy(1:ndx, 1:ndy) )
  
  ! Define a global array, get a handle in id_xy
  ! option 2: define with string of actual values
  !          -> we can write the array only 
  !          -> but have to create strings here
  write (local,  '(i0,",",i0)') nx_local, ny_local
  write (global, '(i0,",",i0)') nx_global, ny_global
  write (offset, '(i0,",",i0)') offs_x, offs_y
  !print '("local=",a)', local

  call adios_define_var (m_adios_group, "xy", "",  &
                         adios_double,  &
                         local, global, offset, id1_xy)

  ! second block
  offs_x = posx * ndx + npx * ndx
  offs_y = posy * ndy
  write (offset, '(i0,",",i0)') offs_x, offs_y

  call adios_define_var (m_adios_group, "xy", "",  &
                         adios_double,  &
                         local, global, offset, id2_xy)



  do ts=0,timesteps-1
     write(filename,'(a6,i2.2,a3)') 'writer',ts,'.bp'
     if (rank==0) write(6,*) 'ts=',ts
     xy = 1.0*rank + 1.0*ts
     call adios_open (adios_handle, "writer", trim(filename), "w", group_comm, ierr)     
!#include "gwrite_writer.fh" --> add write calls manually
     ! calculate how many bytes we are going to write
     !  twice a block of nx_local*ny_local of 8 bytes
     adios_groupsize = 6*4 &
                     + 2 *  8 * (nx_local) * (ny_local)

     ! write first data block
     call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, ierr)
     call adios_write_byid (adios_handle, id1_xy, xy,  ierr)
     ! write second data block (here the same array, xy, again)
     call adios_write_byid (adios_handle, id2_xy, xy,  ierr)
     
     call adios_close (adios_handle,  ierr)
     call MPI_Barrier (group_comm, ierr)
  enddo
  ! Terminate
  call adios_finalize (rank,  ierr)
  call MPI_Finalize (ierr)
end program writer
