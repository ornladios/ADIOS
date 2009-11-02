! ADIOS Fortran Example: write a global array from N processors with gwrite
!
! How to run: mpirun -np <N> adiosf_global
! Output: adiosf_global.bp
! ADIOS config file: adiosf_global.xml
!

program adios_global 
    implicit none
    include 'mpif.h'
    character(len=256)      :: filename = "adiosf_global.bp"
    integer                 :: rank, size, i, ierr
    integer                 :: NX = 10
    real*8, dimension(NX)   :: t
    integer                 :: comm

    ! ADIOS variables declarations for matching gwrite_temperature.fh 
    integer                 :: adios_err
    integer*8               :: adios_groupsize, adios_totalsize
    integer*8               :: adios_handle

    call MPI_Init (ierr)
    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
    call MPI_Comm_rank (comm, rank, ierr)
    call MPI_Comm_size (comm, size, ierr);

    do i = 1, NX
        t(i)  = 10*rank+i;
    enddo

    call adiosf_init ("adiosf_global.xml", ierr);

    call adiosf_open (adios_handle, "temperature", filename, "w", ierr);
#include "gwrite_temperature.fh"
    call adiosf_close (adios_handle, ierr);

    call MPI_Barrier (comm, ierr)

    call adiosf_finalize (rank, ierr);

    call MPI_Finalize ();
end program
