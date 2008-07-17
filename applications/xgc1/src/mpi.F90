!MPI related function 
! for absoft compiler call external name is UPPERCASE and 
! all mpi name is lowercase


! mympif.h is the same as mpif.h except the all text is lowercase.

subroutine MY_MPI_INIT
  use sml_module
  implicit none
  include 'mpif.h'
  integer :: ierror

  call mpi_init(ierror)
  call mpi_comm_size(mpi_comm_world,sml_totalpe,ierror)
  call mpi_comm_rank(mpi_comm_world,sml_mype,ierror)

end subroutine MY_MPI_INIT


subroutine MY_MPI_FINALIZE
  use sml_module
  implicit none
  include 'mpif.h'
  integer :: ierror
  
  call mpi_finalize(mpi_comm_world,sml_mype,ierror)

end subroutine MY_MPI_FINALIZE

subroutine MY_MPI_REDUCE(a,b,n)
  implicit none
  include 'mpif.h'
  integer ,intent(in):: n
  real (kind=8), intent(inout) :: a(n), b(n)
  integer :: ierror
  
  call mpi_reduce(a,b,n,mpi_real8,mpi_sum,0,mpi_comm_world,ierror)

end subroutine MY_MPI_REDUCE


subroutine MY_MPI_ALLREDUCE(a,b,n)
  implicit none
  include 'mpif.h'
  integer ,intent(in):: n
  real (kind=8), intent(inout) :: a(n), b(n)
  integer :: ierror
  
  call mpi_allreduce(a,b,n,mpi_real8,mpi_sum,mpi_comm_world,ierror)

end subroutine MY_MPI_ALLREDUCE
