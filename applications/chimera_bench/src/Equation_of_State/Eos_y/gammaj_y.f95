SUBROUTINE gammaj_y( j, rho, t, ji_ray, jk_ray, gamma1, gamma2, gamma3 )
!-----------------------------------------------------------------------
!
!    File:         gammaj_y
!    Module:       gammaj_y
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the three adiabatic exponents, gamma1, gamma2, and gamma3,
!       for angular zone j. Here
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
!  j                : y (angular) zone index.
!  rho              : shifted matter density array (g/cm**3).
!  t                : shifted matter matter temperature array (K).
!  ji_ray           : x (radial) index of a specific y (angular) ray
!  jk_ray           : z (azimuthal) index of a specific y (angular) ray
!
!    Output arguments:
!
!  gamma1           : first adiabatic index.
!  gamma2           : second adiabatic index.
!  gamma3           : third adiabatic index.
!
!    Input arguments (common):
!
!  aesv(j,i,ji_ray,jk_ray)  : interpolated equation of state quantity i for angular zone j.
!  aesvd(j,i,ji_ray,jk_ray) : d(aesv(j,i,ji_ray,jk_ray))/d(rho).
!  aesvd(j,i,ji_ray,jk_ray) : d(aesv(j,i,ji_ray,jk_ray))/d(t).
!  rho(j)           : density of angular zone j (g/cm**3).
!  t(j)             : temperature of angular zone j (K).
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  eos_snc_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY: ny
USE numerical_module, ONLY : one, epsilon

USE eos_snc_y_module, ONLY : aesv, aesvd, aesvt

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j               ! y (angular) zone index
INTEGER, INTENT(in)              :: ji_ray          ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray          ! z (azimuthal) index of a specific y (angular) ray

REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: t   ! shifted matter matter temperature array (K)

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

gamma3       = aesvt(j,1,ji_ray,jk_ray)/( aesvt(j,2,ji_ray,jk_ray) * rho(j) + epsilon ) + one
gamma1       = ( t(j) * aesvt(j,1,ji_ray,jk_ray) * ( gamma3 - one ) + rho(j) * aesvd(j,1,ji_ray,jk_ray) ) &
&            / ( aesv(j,1,ji_ray,jk_ray) + epsilon )
gammat       = ( rho(j) * aesvd(j,1,ji_ray,jk_ray)/( gamma3 - one ) + t(j) * aesvt(j,1,ji_ray,jk_ray) ) &
&            / ( aesv(j,1,ji_ray,jk_ray) + epsilon )
gamma2       = gammat/( gammat - one )

RETURN
END SUBROUTINE gammaj_y
