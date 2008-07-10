SUBROUTINE gammaj_x( j, rho, t, ij_ray, ik_ray, gamma1, gamma2, gamma3 )
!-----------------------------------------------------------------------
!
!    File:         gammaj_x
!    Module:       gammaj_x
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the three adiabatic exponents, gamma1, gamma2, and gamma3,
!       for radial zone j. Here
!
!                       d(lnP)
!           gamma1   = --------    (constant s and ye)
!                      d(lnrho)
!
!           gamma2     d(lnP)
!         ---------- = ------      (constant s and ye)
!         gamma2 - 1   d(lnt)
!
!                       d(lnt)
!         gamma3 - 1 = --------    (constant s and ye)
!                      d(lnrho)
!
!    Input arguments:
!
!  j                : radial zone index.
!  rho              : shifted matter density array (g/cm**3).
!  t                : shifted matter matter temperature array (K).
!  ij_ray           : j-index of a radial ray
!  ik_ray           : k-index of a radial ray
!
!    Output arguments:
!
!  gamma1           : first adiabatic index.
!  gamma2           : second adiabatic index.
!  gamma3           : third adiabatic index.
!
!    Input arguments (common):
!
!  aesv(j,i,ij_ray,ik_ray)  : interpolated equation of state quantity i for radial zone j.
!  aesvd(j,i,ij_ray,ik_ray) : d(aesv(j,i,ij_ray,ik_ray))/d(rho).
!  aesvd(j,i,ij_ray,ik_ray) : d(aesv(j,i,ij_ray,ik_ray))/d(t).
!  rho(j)                   : density of radial zone j (g/cm**3).
!  t(j)                     : temperature of radial zone j (K).
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  eos_snc_x_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY: nx
USE numerical_module, ONLY : one, epsilon

USE eos_snc_x_module, ONLY : aesv, aesvd, aesvt

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j               ! radial zone index
INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nx) :: t   ! shifted matter matter temperature array (K)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: gamma1          ! Gamma_{1}
REAL(KIND=double), INTENT(out)   :: gamma2          ! Gamma_{2}
REAL(KIND=double), INTENT(out)   :: gamma3          ! Gamma_{3}

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

REAL(KIND=double)                :: gammat          ! temporary variable

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

gamma3       = aesvt(j,1,ij_ray,ik_ray)/( aesvt(j,2,ij_ray,ik_ray) * rho(j)  &
&            + epsilon ) + one
gamma1       = ( t(j) * aesvt(j,1,ij_ray,ik_ray) * ( gamma3 - one ) + rho(j) &
&            * aesvd(j,1,ij_ray,ik_ray) )/( aesv(j,1,ij_ray,ik_ray) + epsilon )
gammat       = ( rho(j) * aesvd(j,1,ij_ray,ik_ray)/( gamma3 - one ) + t(j)   &
&            * aesvt(j,1,ij_ray,ik_ray) )/( aesv(j,1,ij_ray,ik_ray) + epsilon )
gamma2       = gammat/( gammat - one )

RETURN
END SUBROUTINE gammaj_x
