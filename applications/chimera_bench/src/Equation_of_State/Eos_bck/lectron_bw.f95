SUBROUTINE lectron
!-----------------------------------------------------------------------
!
!    File:         lectron
!    Module:       lectron
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!	 To compute the electron-positron contribution to the equation of state,
!	  given rho, t, and ye. The electron chemical potential is iterated until
!     the appropriate Fermi-Dirac integrals, or approximations thereof, for
!     the electron - positron number converge to the electron positron number
!     given by rho*ye/mb.
!
!    Variables that must be passed through common:
!  dbck       :  density (baryons/fm**3)
!  tbck       :  temperature (MeV)
!  yebck      :  Y_e^- - Y_e^+  (ye)
!
!    Subprograms called:
!        glaquad
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!        kind_module, array_module, numerical_module, physcnst_module
!        edit_module, eos_bck_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module
USE numerical_module, ONLY : zero, third, half
USE physcnst_module, ONLY : me, pi, hbarc, wnm, ws, dmnp

USE edit_module, ONLY : nlog
USE eos_bck_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                            :: first = .true.

REAL(KIND=double), PARAMETER       :: gelec = 2.d0   ! electron spin degeneracy
REAL(KIND=double)                  :: c0             ! constant for computing the chemical potential
REAL(KIND=double)                  :: ek             ! electron kinetic energy/baryon (MeV)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( first ) THEN
  c0            = ( me/( 2.d0 * pi * hbarc**2 ) )**( 3.d0/2.d0 )
END IF

therm           = c0 * tbck * dsqrt(tbck) / dbck
pe              = dbck * yebck * tbck
ek              = 1.5d0 * pe
ee              = ek/dbck ! + yebck * me
ue              = tbck * DLOG( DABS( yebck/( therm * gelec ) ) )
se              = yebck * ( 2.5d0 - ue/tbck )
yeplus          = zero
rel             = zero

RETURN


END SUBROUTINE lectron
