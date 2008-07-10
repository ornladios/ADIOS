SUBROUTINE dimension_convect_arrays( nx )
!-----------------------------------------------------------------------
!
!    File:         dimension_convect_arrays
!    Module:       dimension_convect_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the convection arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!   nx        : x (radial) array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  convect_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE convect_module
USE edit_module, ONLY : nlog
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

  101 FORMAT (' ML convection arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_convect_arrays')

!-----------------------------------------------------------------------
!        Allocate convect_module arrays
!-----------------------------------------------------------------------

ALLOCATE (aledoux(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'aledoux   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psclht(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psclht    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (usound(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'usound    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (lmix(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'lmix      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ulmixl(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ulmixl    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ulcnvct(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ulcnvct   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dmcnvct(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dmcnvct   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (pcnvct(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'pcnvct    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (scnvct(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'scnvct    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ INITIALIZE CONVECT_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

aledoux                   = zero
psclht                    = zero
usound                    = zero
lmix                      = zero
ulmixl                    = zero
ulcnvct                   = zero
dmcnvct                   = zero
pcnvct                    = zero
scnvct                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_convect_arrays
