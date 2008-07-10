SUBROUTINE initialize_bomb_arrays
!-----------------------------------------------------------------------
!
!    File:         initialize_bomb_arrays
!    Module:       initialize_bomb_arrays
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
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  bomb_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero

USE bomb_module
USE edit_module, ONLY : nlog
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Bomb variables have been initialized')

!-----------------------------------------------------------------------
!        Initialize bomb_module arrays
!-----------------------------------------------------------------------

jexpl_min                 = 0
jexpl_max                 = 0

e_bomb                    = zero
bomb_time                 = 1.d+100
t_start_bomb              = 1.d+100

IF ( myid == 0 ) WRITE (nlog,101)

RETURN
END SUBROUTINE initialize_bomb_arrays
                                                                                
