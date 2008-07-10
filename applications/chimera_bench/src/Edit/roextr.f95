SUBROUTINE roextr
!-----------------------------------------------------------------------
!
!    File:         roextr
!    Module:       roextr
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/2/00
!
!    Purpose:
!      To determine the maximum and minimum values attained by the
!        central density between calls to editc.
!
!    Subprograms called:
!      none
!
!    Input arguments:
!      none
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!
!  rho(j)     : density of zone j (g/cm3)
!
!    Include files:
!  kind_module, numerical_module
!  mdl_cnfg_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero
      
USE mdl_cnfg_module, ONLY : rho

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)                :: rhomax = 0.0d+00 ! running maximum density
REAL(KIND=double)                :: rhomin = 1.d+100 ! running minimum density
REAL(KIND=double)                :: rhomaxrd      ! maximum density since last called
REAL(KIND=double)                :: rhominrd      ! minimum density since last called

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( rho(2) > rhomax ) rhomax = rho(2)
IF ( rho(2) < rhomin ) rhomin = rho(2)

RETURN

ENTRY roextrrd(rhominrd,rhomaxrd)

rhominrd           = rhomin
rhomaxrd           = rhomax
rhomin             = 1.d+100
rhomax             = zero

RETURN
END
