SUBROUTINE gammaz_z( kmin, kmax, rho, t, ki_ray, kj_ray )
!-----------------------------------------------------------------------
!
!    File:         gammaz_z
!    Module:       gammaz_z
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To compute the three adiabatic exponents, gam1(j,kj_ray,ki_ray),
!       gam2(j,kj_ray,ki_ray), and gam3(j,kj_ray,ki_ray), for all angular
!       zones j between kmin and kmax.
!       Here
!
!                       d(lnP)
!          gam1(j)   = --------    (constant s and ye)
!                      d(lnrho)
!
!          gam2(j)     d(lnP)
!        ----------- = ------      (constant s and ye)
!        gam2(j) - 1   d(lnt)
!
!                       d(lnt)
!        gam3(j) - 1 = --------    (constant s and ye)
!                      d(lnrho)
!
!    Input arguments:
!  kmin             : minimum z (azimuthal) zone index.
!  kmax             : maximum z (azimuthal) zone index.
!  rho              : azimuthal array of matter density (g/cm**3)
!  t                : azimuthal array of matter temperature (K)
!  ki_ray           : x (radial) index of a specific z (azimuthal) ray
!  kj_ray           : y (angular) index of a specific z (azimuthal) ray
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  aesv(j,i,kj_ray,ki_ray)  : interpolated equation of state quantity i for angular zone j.
!  aesvd(j,i,kj_ray,ki_ray) : d(aesv(j,i))/d(rho).
!  aesvd(j,i,kj_ray,ki_ray  : d(aesv(j,i))/d(t).
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module,
!  eos_snc_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nz
USE numerical_module, ONLY : one, epsilon

USE eos_snc_z_module, ONLY : aesv, aesvd, aesvt, gam1, gam2, gam3

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: kmin            ! minimum z (azimuthal) zone index
INTEGER, INTENT(in)              :: kmax            ! maximum z (azimuthal) zone index
INTEGER, INTENT(in)              :: ki_ray          ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray          ! y (angular) index of a specific z (azimuthal) ray

REAL(KIND=double), INTENT(in), DIMENSION(nz) :: rho ! azimuthal array of matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(nz) :: t   ! azimuthal array of matter matter temperature array (K)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j               ! angular zone index

REAL(KIND=double)                :: gammat          ! temporary variable

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO j = kmin,kmax
  gam3(j,kj_ray,ki_ray) = aesvt(j,1,kj_ray,ki_ray)                      &
&                       / ( aesvt(j,2,kj_ray,ki_ray) * rho(j) + epsilon ) + one
  gam1(j,kj_ray,ki_ray) = ( t(j) * aesvt(j,1,kj_ray,ki_ray)             &
&                       * ( gam3(j,kj_ray,ki_ray) - one ) + rho(j)      &
&                       * aesvd(j,1,kj_ray,ki_ray) ) &
&                       / ( aesv(j,1,kj_ray,ki_ray) + epsilon )
  gammat                = ( rho(j) * aesvd(j,1,kj_ray,ki_ray)           &
&                       / ( gam3(j,kj_ray,ki_ray) - one ) + t(j)        &
&                       * aesvt(j,1,kj_ray,ki_ray) ) &
&                       / ( aesv(j,1,kj_ray,ki_ray) + epsilon )
  gam2(j,kj_ray,ki_ray) = gammat/( gammat - one )
END DO

RETURN
END SUBROUTINE gammaz_z
