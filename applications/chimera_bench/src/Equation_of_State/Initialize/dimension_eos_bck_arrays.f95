SUBROUTINE dimension_eos_bck_arrays( nx )
!-----------------------------------------------------------------------
!
!    File:         dimension_eos_bck_arrays
!    Module:       dimension_eos_bck_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to eos bck arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x-array dimensions
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, eos_bck_module, parallel_module
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE eos_bck_module
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! radial array dimension

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' EOS BCK arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_eos_bck_arrays')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ ALLOCATE EOS_BCK_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (uea(nx,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uea       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (una(nx,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'una       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uhata(nx,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uhata     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (thetaa(nx,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'thetaa    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zaa(nx,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zaa       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xaa(nx,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xaa       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dtrana(nx,2,2,2), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtrana    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ueaa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ueaa      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (unaa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'unaa      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uhataa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uhataa    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (thetaaa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'thetaaa   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (zaaa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'zaaa      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (xaaa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xaaa      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dtranaa(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dtranaa   '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!             \\\\\ INITIALIZE EOS_BCK_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

uea                       = zero
una                       = zero
uhata                     = zero
thetaa                    = zero
zaa                       = zero
xaa                       = zero
dtrana                    = zero

ueaa                      = zero
unaa                      = zero
uhataa                    = zero
thetaaa                   = zero
zaaa                      = zero
xaaa                      = zero
dtranaa                   = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_eos_bck_arrays
