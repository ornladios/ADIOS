!-----------------------------------------------------------------------
!    Module:       parallel_module
!    Author:       S. W. Bruenn
!    Date:         2/25/03
!
!    Declares all variables needed for parallelization
!-----------------------------------------------------------------------

MODULE parallel_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!  MPI Variables
!-----------------------------------------------------------------------

INTEGER                          :: ierr         ! initialization variable for MPI
INTEGER                          :: myid         ! rank of each processor (MPI)          
INTEGER                          :: myid_y       ! rank of each processor when split wrt k_block         
INTEGER                          :: myid_z       ! rank of each processor when split wrt j_block 
INTEGER                          :: MPI_COMM_ROW ! new communicator handle
INTEGER                          :: MPI_COMM_COL ! new communicator handle

END MODULE parallel_module
