SUBROUTINE initialize_variables
!-----------------------------------------------------------------------
!
!    File:         initialize_variables
!    Module:       initialize_variables
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/17/04
!
!    Purpose:
!      To initialize variables before read in not initialized in the
!       dimension array subroutines.
!
!    Subprograms called:
!  initialize_global_var
!  initialize_cycle_arrays
!  initialize_it_tol_arrays
!  initialize_bomb_arrays
!  initialize_rezone_arrays
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!        none
!
!-----------------------------------------------------------------------

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize global arrays.
!-----------------------------------------------------------------------

CALL initialize_global_var

!-----------------------------------------------------------------------
!  Initialize cycle and tolerance arrays.
!-----------------------------------------------------------------------

CALL initialize_cycle_arrays
CALL initialize_it_tol_arrays

!-----------------------------------------------------------------------
!  Initialize bomb and rezone arrays.
!-----------------------------------------------------------------------

CALL initialize_bomb_arrays
CALL initialize_rezone_arrays

RETURN
END SUBROUTINE initialize_variables
