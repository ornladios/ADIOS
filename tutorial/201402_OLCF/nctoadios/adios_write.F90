
program simple_xy_par_wr
  use adios_write_mod
  implicit none
  include 'mpif.h'

  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "simple_xy_par.bp"

  ! We are writing 2D data.
  integer, parameter :: NDIMS = 2
  integer, parameter :: N = 5  ! number of elements written by one process

  ! ADIOS uses 64bit integer ids
  integer*8 :: adios_handle, adios_groupsize, adios_totalsize
  integer :: adios_err
  integer :: x_dimid, y_dimid

  ! These will tell where in the data file this processor should
  ! write.
  integer :: start(NDIMS), count(NDIMS)

  ! This is the data array we will write. It will just be filled with
  ! the rank of this processor.
  integer, allocatable :: data_out(:)

  ! MPI stuff: number of processors, rank of this processor, and error
  ! code.
  integer :: p, my_rank, ierr

  ! Loop indexes, and error handling.
  integer :: x, stat

  ! Initialize MPI, learn local rank and total number of processors.
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

  call adios_init ("adios_write.xml", MPI_COMM_WORLD, adios_err)

  ! Create some pretend data. We just need one row.
  allocate(data_out(N), stat = stat)
  if (stat .ne. 0) stop 3
  do x = 1, N
     data_out(x) = my_rank
  end do

  call adios_open (adios_handle, "writer", trim(FILE_NAME), "w", MPI_COMM_WORLD, adios_err)
#include "gwrite_writer.fh"
  call adios_close (adios_handle, adios_err)

  ! Free my local memory.
  deallocate(data_out)

  ! MPI library must be shut down.
  call MPI_Finalize(ierr)

  if (my_rank .eq. 0) print *, "*** SUCCESS writing example file ", FILE_NAME, "! "

end program simple_xy_par_wr

