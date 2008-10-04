      subroutine ckHeaderSetup( )
!
! Author: James Sutherland (sutherland@crsim.utah.edu)
! Date:   Sept 9, 2001
!
! This routine initializes the ckstrt.h header variables for parallel runs.
! It must be called by all processors after CKINIT has been called by process 0
!
      use topology_m
      implicit none
      include 'ckstrt.h'

  ! if ckstrt.h is modified, correct the following calls to broadcast
  ! the commons accurately

!******************************************
! The original ckstrt.h used in S3D was a Chemkin II header file (I think).
! June 6, 2002  James upgraded everything to Chemkin III (not sure exactly what version)
! and modified this accordingly.
!      call MPI_Bcast(nmm,    104,     MPI_INTEGER,   0, gcomm, ierr)
!      call MPI_Bcast(LPERT,  1,       MPI_LOGICAL,   0, gcomm, ierr)
!******************************************

! Chemkin III ckstrt.h

!jcs bug fix 6-18-03 should be 109, not 108 (as previously)...
      call MPI_Bcast(nmm,    109,     MPI_INTEGER,   0, gcomm, ierr)
      call MPI_Bcast(LPERT,  1,       MPI_LOGICAL,   0, gcomm, ierr)

      end subroutine ckHeaderSetup
