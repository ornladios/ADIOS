SUBROUTINE gammaj_z( k, rho, t, ki_ray, kj_ray, gamma1, gamma2, gamma3 )
!-----------------------------------------------------------------------
!
!    File:         gammaj_z
!    Module:       gammaj_z
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the three adiabatic exponents, gamma1, gamma2, and gamma3,
!       for angular zone k. Here
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
!  k                : z (azimuthal) zone index
!  rho              : azimuthal array of matter density (g/cm**3)
!  t                : azimuthal array of matter temperature (K)
!  ki_ray           : x (radial) index of a specific z (azimuthal) ray
!  kj_ray           : y (angular) index of a specific z (azimuthal) ray
!
!    Output arguments:
!
!  gamma1           : first adiabatic index.
!  gamma2           : second adiabatic index.
!  gamma3           : third adiabatic index.
!
!    Input arguments (common):
!
!  aesv(k,i,kj_ray,ki_ray)  : interpolated equation of state quantity i for angular zone k.
!  aesvd(k,i,kj_ray,ki_ray) : d(aesv(k,i,kj_ray,ki_ray))/d(rho).
!  aesvd(k,i,kj_ray,ki_ray) : d(aesv(k,i,kj_ray,ki_ray))/d(t).
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  eos_snc_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY: nz
USE numerical_module, ONLY : one, epsilon

USE eos_snc_z_module, ONLY : aesv, aesvd, aesvt

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: k               ! z (azimuthal) zone index
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray

REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rho ! azimuthal array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: t   ! azimuthal array of matter matter temperature array (K)

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

gamma3       = aesvt(k,1,kj_ray,ki_ray)/( aesvt(k,2,kj_ray,ki_ray) * rho(k) + epsilon ) + one
gamma1       = ( t(k) * aesvt(k,1,kj_ray,ki_ray) * ( gamma3 - one ) + rho(k) * aesvd(k,1,kj_ray,ki_ray) ) &
&            / ( aesv(k,1,kj_ray,ki_ray) + epsilon )
gammat       = ( rho(k) * aesvd(k,1,kj_ray,ki_ray)/( gamma3 - one ) + t(k) * aesvt(k,1,kj_ray,ki_ray) ) &
&            / ( aesv(k,1,kj_ray,ki_ray) + epsilon )
gamma2       = gammat/( gammat - one )

RETURN
END SUBROUTINE gammaj_z
