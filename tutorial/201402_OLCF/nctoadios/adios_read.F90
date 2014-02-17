program simple_xy_par_rd
  use adios_read_mod
  implicit none
  include 'mpif.h'

  ! This is the name of the data file we will read.
  character (len = *), parameter :: FILE_NAME = "simple_xy_par.bp"

  ! These will tell where in the data file this processor should
  ! write.
  integer, parameter :: NDIMS = 2
  integer :: N
  integer*8 :: start(NDIMS), count(NDIMS)

  ! We will read data into this array.
  integer, allocatable :: data_in(:)

  ! ADIOS variables
  integer*8           :: fh
  integer*8           :: sel  ! ADIOS selection object


  ! MPI stuff: number of processors, rank of this processor, and error
  ! code.
  integer :: p, my_rank, ierr

  ! Loop indexes, and error handling.
  integer :: x, y, stat


  ! Initialize MPI, learn local rank and total number of processors.
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
          "verbose=3", ierr)

  ! Open the file as "file", i.e. see all timesteps at once.
  call adios_read_open_file (fh, FILE_NAME, ADIOS_READ_METHOD_BP, &
          MPI_COMM_WORLD, ierr)

  ! Get the dimensions
  call adios_get_scalar(fh, "x", p, ierr)
  call adios_get_scalar(fh, "y", N, ierr)

  print *, "Global dimensions:",p,"x",N

  ! Allocate space to read in data.
  allocate(data_in(N), stat = stat)
  if (stat .ne. 0) stop 3

  ! Create a 2D selection for the subset
  start = (/ 0, my_rank/)
  count = (/ N, 1 /)
  call adios_selection_boundingbox (sel, 2, start, count)

  ! Read the data.
  ! Arrays are read by scheduling one or more of them
  ! and performing the reads at once
  call adios_schedule_read (fh, sel, "data", 0, 1, data_in, ierr)
  call adios_perform_reads (fh, ierr)

  ! Check the data.
  do x = 1, N
     if (data_in(x) .ne. my_rank) then
        print *, "data_in(", x, ") = ", data_in(x)
        stop "Stopped"
     endif
  end do

  ! Close the file, freeing all resources.
  call adios_read_close (fh, ierr)

  ! Free my local memory.
  deallocate(data_in)
  call adios_selection_delete (sel)

  ! ADIOS must be finalized
  call adios_read_finalize_method (ADIOS_READ_METHOD_BP, ierr)

  ! MPI library must be shut down.
  call MPI_Finalize(ierr)

  if (my_rank .eq. 0) print *,"*** SUCCESS reading example file ", FILE_NAME, "! "

end program simple_xy_par_rd

