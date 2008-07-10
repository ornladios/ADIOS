SUBROUTINE dimension_tov_pot_arrays( nx )
!-----------------------------------------------------------------------
!
!    File:         dimension_tov_pot_arrays
!    Module:       dimension_tov_pot_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the tov_pot arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x (radial) array dimension
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, tov_potential_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE tov_potential_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)               :: var_name

INTEGER                          :: istat         ! allocation status

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' tov_pot arrays have been dimensioned')
 1001 FORMAT (' Allocation problem for array ',a10,' in dimension_tov_pot_arrays')

!-----------------------------------------------------------------------
!
!                \\\\\ ALLOCATE TOV_POT ARRAYS /////
!
!-----------------------------------------------------------------------

ALLOCATE (effpot_c(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'effpot_c  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (effpot_e(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'effpot_e  '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!              \\\\\ INITIALIZE TOV_POT ARRAYS /////
!
!-----------------------------------------------------------------------

effpot_c                  = zero
effpot_e                  = zero

WRITE (nlog,101)

RETURN
END SUBROUTINE dimension_tov_pot_arrays
                                                                                
