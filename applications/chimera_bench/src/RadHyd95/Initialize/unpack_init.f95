SUBROUTINE unpack_init( c_init_data, i_init_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_init
!    Module:       unpack_init
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/25/03
!
!    Purpose:
!      To read in the model configuration and the problem
!       specification parameters
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp     : unit number from which to read
!
!    Output arguments:
!  c_eos_data : character array of initial data
!
!    Include files:
!  cycle_module, edit_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : nrst
USE edit_module, ONLY : head, nouttmp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER(len = 128), INTENT(in), DIMENSION(1) :: c_init_data  ! character array of initial data

INTEGER, INTENT(in), DIMENSION(2)              :: i_init_data  ! integer array of initial data

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                   \\\\\ UNPACK INIT KEYS /////
!
!-----------------------------------------------------------------------

head                    = c_init_data(1)
nrst                    = i_init_data(1)
nouttmp                 = i_init_data(2)

RETURN
END SUBROUTINE unpack_init
