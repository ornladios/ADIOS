program writer
  use adios_write_mod
  implicit none
  include 'mpif.h'
  character(len=256) :: filename
  integer :: npx, npy   ! # of processors in x-y direction
  integer :: ndx, ndy   ! size of array per processor
  integer :: timesteps=2  ! number of timesteps to run  
  integer :: nproc      ! number of processors
  
  integer :: i,j
  real*8 :: hx,hy,PI, x, y
  real*8, dimension(:,:),   allocatable :: xy
  ! Offsets and sizes
  integer :: offs_x, offs_y
  integer :: nx_local, ny_local
  integer :: nx_global, ny_global
  integer :: posx, posy   ! position index in the array
  ! MPI variables
  integer :: group_comm
  integer :: rank
  integer :: ierr
  ! ADIOS variables (used in gwrite_writer.fh)
  character (len=200) :: group
  integer*8 :: adios_handle, adios_totalsize, adios_groupsize, adios_buf_size
  integer :: adios_err
  ! actual timestep
  integer   :: ts 

  
  call MPI_Init (ierr)
  call MPI_Comm_dup (MPI_COMM_WORLD, group_comm, ierr)
  call MPI_Comm_rank (group_comm, rank, ierr)
  call MPI_Comm_size (group_comm, nproc , ierr)
  call adios_init ("writer.xml", group_comm, adios_err)
  ! will work with 12 cores, which are arranged by npx=4, npy=3 (4x3)
  npx = 4
  npy = 3
  !ndx = 65
  !ndy = 129
  ndx = 650
  ndy = 1290
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
  PI = 4.0*atan(1.0)
  hx= 2.0 * PI / (nx_global-1) 
  hy= 2.0 * PI / (ny_global-1) 
  

  do ts=0,timesteps-1
     write(filename,'(a6,i2.2,a3)') 'writer',ts,'.bp'
     if (rank==0) write(6,*) 'ts=',ts
      do j = 1 , ny_local
      do i = 1 , nx_local
            x = posx*(nx_local *hx) - PI  + (i-1) * hx
            y = posy*(ny_local *hy) - PI + (j-1) * hy
            xy(i,j) = (ts+1) * (cos(x)+sin(y)) 
         end do
      end do
     !xy = rank
     call adios_open (adios_handle, "writer", trim(filename), "w", group_comm, adios_err)
#include "gwrite_writer.fh"
     call adios_close (adios_handle, adios_err)
     
  enddo

  ! Terminate
  call MPI_Barrier (group_comm, ierr)
  call adios_finalize (rank, ierr)
  call MPI_Finalize (ierr)
end program writer
