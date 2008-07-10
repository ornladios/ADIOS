SUBROUTINE initialize_rezone_arrays
!-----------------------------------------------------------------------
!
!    File:         initialize_rezone_arrays
!    Module:       initialize_rezone_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To allocate dimensions to the mgfld remap arrays.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  edit_module, parallel_module, rezone_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid
USE rezone_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Rezone variables have been initialized')

!-----------------------------------------------------------------------
!        Initialize initialize_rezone_arrays
!-----------------------------------------------------------------------

!.........Eulerian rezoning controls....................................

n_eulgrid                 = 2
n1zoom                    = 0
n2zoom                    = 0
n3zoom                    = 0

r_1                       = zero
r_2                       = zero
r_3                       = zero
zoome1                    = zero
zoome2                    = zero
zoome3                    = zero

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE initialize_rezone_arrays
                                                                                
