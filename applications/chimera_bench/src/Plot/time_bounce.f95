SUBROUTINE time_bounce( rhobar, nx )
!-----------------------------------------------------------------------
!
!    File:         time_bounce
!    Module:       time_bounce
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/2/00
!
!    Purpose:
!      To determine the time of core bounce.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!  rhobar     : mean density (cm^{-3})
!  nx         : x_array extent
!
!    Output arguments:
!      none
!
!    Include files:
!      numerical_module
!      t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE t_cntrl_module, ONLY : t_bounce, time

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                             :: nx            ! x-array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx)    :: rhobar        ! mean density (MGFLD indexed( (cm^{-3})

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: i_b            ! radial zone index
INTEGER, PARAMETER                :: n_bounce = 100 ! number of cycles rho(j) < rhomax criterion

REAL(KIND=double)                 :: rhomax = 0.d0  ! maximum density up to current cycle
REAL(KIND=double)                 :: t_b    = 0.d0  ! time of maximum density up to current cycle

!-----------------------------------------------------------------------
!        Return if t_bounce ne 0.
!-----------------------------------------------------------------------

IF ( t_bounce /= zero ) RETURN

!-----------------------------------------------------------------------
!        Determine when core bounce occurs.
!
!        Bounce has occurred if rho(2) < rhomax for at least
!         n_bounce cycles.
!-----------------------------------------------------------------------

IF ( rhobar(2) > rhomax ) THEN
  rhomax           = rhobar(2)
  t_b              = time
  i_b              = 0
ELSE
  i_b              = i_b + 1
  IF ( i_b > n_bounce ) t_bounce = t_b
END IF ! rho(2) > rhomax

RETURN
END
