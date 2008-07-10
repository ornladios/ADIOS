SUBROUTINE e_advct_evol
!-----------------------------------------------------------------------
!
!    File:         e_advct_evol
!    Module:       e_advct_evol
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/28/01
!
!    Purpose:
!      To update the energy zone boundaries and the psi0.
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
!    Input arguments (common):
!  xk(k)       : energy grid zone edge locations prior to Lagrangian update
!  dk(k)       : xk(k+1) - xk(k) prior to Lagrangian update
!  dvolk_0(l)  : energy space volume prior to Lagrangian update
!  psi(k)      : neutrino occupation number prior to Lagrangian update
!
!    Output arguments (common):
!  xk(k)       : energy grid zone edge locations after Lagrangian update
!  dk(k)       : xk(k+1) - xk(k) after Lagrangian update
!  dvolk(l)    : energy space volume after Lagrangian update
!  psi(k)      : neutrino occupation number after Lagrangian update
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : third

USE e_advct_module, ONLY : nmin, nmax, xk, vk, dk, dvolk, dvolk_0, psi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Energy grid position evolution
!-----------------------------------------------------------------------

xk(nmin-3:nmax+4)    = xk(nmin-3:nmax+4) + vk(nmin-3:nmax+4)
dk(nmin-4:nmax+5)    = xk(nmin-3:nmax+6) - xk(nmin-4:nmax+5)

!-----------------------------------------------------------------------
!  Calculate the energy space volumes
!-----------------------------------------------------------------------

CALL e_advct_vol

!-----------------------------------------------------------------------
!  Lagrangian evolution of psi0
!-----------------------------------------------------------------------

psi(nmin-3:nmax+3)   = psi(nmin-3:nmax+3) * ( dvolk_0(nmin-3:nmax+3)/dvolk(nmin-3:nmax+3) )

RETURN
END SUBROUTINE e_advct_evol
