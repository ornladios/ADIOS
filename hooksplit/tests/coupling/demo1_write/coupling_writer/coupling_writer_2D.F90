program coupling_writer
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
  integer   :: adios_err
  ! actual timestep
  integer   :: ts 
  
  call MPI_Init (ierr)
  call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
  call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size (group_comm, nproc , ierr)
  call adios_init ("coupling2D_writer.xml", ierr)
  ! will work with 12 cores, which are arranged by npx=4, npy=3 (4x3)
  npx = 4
  npy = 3
  ndx = 65
  ndy = 129
  ! will work with each core writing ndx = 65, ndy = 129, (65*4,129*3) global
  if (nproc.ne.npx*npy) then
     if (rank==0) write(6,*) 'ERROR with nproc=',nproc,'npx=',npx,'npy=',npy
     call adios_finalize (rank, adios_err)
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
  nx_global  = npx * ndx
  ny_global  = npy * ndy
  allocate( xy(1:ndx, 1:ndy) )
  
  do ts=0,timesteps-1
     write(filename,'(a4,i2.2,a3)') 'cpes',ts,'.bp'
     if (rank==0) write(6,*) 'ts=',ts
     xy = 1.0*rank + 1.0*ts
     call adios_open (adios_handle, "writer2D", trim(filename), "w", group_comm, adios_err)     
#include "gwrite_writer2D.fh"
     
     
     call adios_close (adios_handle, adios_err)
     call MPI_Barrier (MPI_COMM_WORLD, ierr)
  enddo
  ! Terminate
  call adios_finalize (rank, adios_err)
  call MPI_Finalize (ierr)
end program coupling_writer
