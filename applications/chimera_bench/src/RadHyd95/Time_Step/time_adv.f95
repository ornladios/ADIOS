SUBROUTINE time_adv( dtnph, dtnmh, time )
!-----------------------------------------------------------------------
!
!    File:         time_adv
!    Module:       time_adv
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To update the time and switch the time step
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  dtnph : time step for current cycle
!  dtnmh : time step for preceding cycle
!  time  : elapsed time
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double

USE t_cntrl_module, ONLY: dtnmh_t=>dtnmh, time_t=>time

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout)    :: dtnph   ! time step for current cycle
REAL(KIND=double), INTENT(inout)    :: dtnmh   ! time step for preceding cycle
REAL(KIND=double), INTENT(inout)    :: time    ! elapsed time

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Update the time and switch the time step.......................

dtnmh          = dtnph
dtnmh_t        = dtnph
time           = time + dtnph
time_t         = time

RETURN

END SUBROUTINE time_adv
