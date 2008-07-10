SUBROUTINE read_pack_init( nreadp, c_init_data, i_init_data, nrst, nouttmp )
!-----------------------------------------------------------------------
!
!    File:         readst_init
!    Module:       readst_init
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/25/03
!
!    Purpose:
!        To read in the model header and cycle number.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp      : unit number from which to read
!
!    Output arguments:
!  c_init_data : character array of initial data containing header
!  nrst        : cycle number to start simulation
!  nouttmp     : unit number to get restart data if nrst /= 0 (obsolete)
!
!    Include files:
!      edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE edit_module, ONLY : nprint, nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp   ! unit number to read from

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER(len = 128), INTENT(out), DIMENSION(1) :: c_init_data  ! character array of initial data

INTEGER, INTENT(out), DIMENSION(2)              :: i_init_data  ! integer array of initial data

INTEGER, INTENT(out)             :: nrst     ! cycle number to start simulation
INTEGER, INTENT(out)             :: nouttmp  ! unit number to get restart data if nrst /= 0

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER(len = 128)             :: head     ! header

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (a128)
    2 FORMAT (1x,a128)
    3 FORMAT (/)
    5 FORMAT (10x,i10)
    7 FORMAT (' nrst=',i10)
    9 FORMAT (20x,i10)
   11 FORMAT (' nouttmp=',i7)
  101 FORMAT (' Init variables have been read and initialized')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Read title card
!-----------------------------------------------------------------------

READ (nreadp,1) head
WRITE (nprint,2) head

!-----------------------------------------------------------------------
!
!                 \\\\\ READ NRST AND NOUTTMP /////
!
!   If nrst = 0, parameters are initialized and the problem configuration
!    is read from unit number nrrstp = nrrst.
!
!   If nrst > 0 and ndim = 1 (ny = 1) problem is read by processor 0 and
!    broadcast to all processors from file
!
!    rst_tmp1.d                : nouttmp = 1
!    rst_tmp2.d                : nouttmp = 2
!    restart_nrst.d            : nouttmp = 3
!    restart_final.d           : nouttmp = 4
!
!   If nrst > 0 and ndim > 1 (ny > 1 amd/or nz > 1) problem keys are read
!    by processor 0 and broadcast to all processors from file
!
!    rst_tmp1_keys.d           : nouttmp = 1
!    rst_tmp2_keys.d           : nouttmp = 2
!    restart_keys_nrst.d       : nouttmp = 3
!    restart_final_keys.d      : nouttmp = 4
!
!   The nuclear abundance data and model configuration are read by all
!    processors from file
!
!    rst_tmp_file1             : nouttmp = 1
!    rst_tmp_file2             : nouttmp = 2
!    restart_model_nrst_myid,d : nouttmp = 3
!    restart_final_mod.d       : nouttmp = 4
!
!-----------------------------------------------------------------------

READ (nreadp,5) nrst
WRITE (nprint,7) nrst

IF ( nrst > 0 ) THEN
  READ (nreadp,9) nouttmp
  WRITE (nprint,11) nouttmp
END IF

!-----------------------------------------------------------------------
!
!                    \\\\\ PACK INIT DATA /////
!
!-----------------------------------------------------------------------

c_init_data(1)                    = head
i_init_data(1)                    = nrst
i_init_data(2)                    = nouttmp

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE read_pack_init
