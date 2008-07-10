SUBROUTINE dimension_boundary_arrays( nx )
!-----------------------------------------------------------------------
!
!    File:         dimension_boundary_arrays
!    Module:       dimension_boundary_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the boundary arrays.
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
!  boundary_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE boundary_module
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

  101 FORMAT (' Boundary arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a11,' in dimension_boundary_arrays')

!-----------------------------------------------------------------------
!        Allocate boundary_module arrays
!-----------------------------------------------------------------------

ALLOCATE (uset(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uset      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!            \\\\\ INITIALIZE BOUNDARY_MODULE ARRAYS /////
!
!-----------------------------------------------------------------------

ipbnd                     = 0
iubcjmn                   = 0
iubcjmx                   = 0
iuset                     = 0

pbound                    = zero
ubcjmx                    = zero
ubcjmn                    = zero
uset                      = zero
r_u                       = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_boundary_arrays
