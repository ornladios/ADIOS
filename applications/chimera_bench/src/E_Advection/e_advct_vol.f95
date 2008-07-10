SUBROUTINE e_advct_vol
!-----------------------------------------------------------------------
!
!    File:         e_advct_vol
!    Module:       e_advct_vol
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/28/01
!
!    Purpose:
!      To compute the energy space volume.
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
!  nmin        : number of first real zone (=1+6)
!  nmax        : number of last  real zone (=nnugp(n)+6)
!  xk(k)       : energy grid zone edge locations
!  dk(k)       : xk(k+1) - xk(k)
!  xk_0(k)     : energy grid zone edge locations prior to lagrangian
!                 update
!  dk_0(k)     : dk(k) prior to lagrangian update
!  xk_1(k)     : energy grid zone edge locations at time m+1
!  dk_1(k)     : dk(k) at time m+1
!
!    Output arguments (common):
!  dvolk(l)    : energy space volume
!  dvolk_0(l)  : energy space volume prior to lagrangian updat
!  dvolk_1(l)  : energy space volume at time m+1!
!
!    Include files:
!  numerical_module
!  e_advct_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : third
USE e_advct_module, ONLY : nmin, nmax, dvolk, dvolk_0, dvolk_1, dk, xk, &
& dk_0, xk_0, dk_1, xk_1

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: l             ! paddwed zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Calculate the volumes
!-----------------------------------------------------------------------

DO l = nmin-3, nmax+4
  dvolk  (l)      = dk  (l) * ( xk ( l) * ( xk  (l) + dk  (l) ) + dk  (l) * dk  (l) * third )
  dvolk_0(l)      = dk_0(l) * ( xk_0(l) * ( xk_0(l) + dk_0(l) ) + dk_0(l) * dk_0(l) * third )
  dvolk_1(l)      = dk_1(l) * ( xk_1(l) * ( xk_1(l) + dk_1(l) ) + dk_1(l) * dk_1(l) * third )
END DO

RETURN
END SUBROUTINE e_advct_vol
