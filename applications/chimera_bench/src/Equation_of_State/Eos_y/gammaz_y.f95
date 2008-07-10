SUBROUTINE gammaz_y( jmin, jmax, rho, t, ji_ray, jk_ray )
!-----------------------------------------------------------------------
!
!    File:         gammaz_y
!    Module:       gammaz_y
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute the three adiabatic exponents, gam1(j,ji_ray,jk_ray), gam2(j,ji_ray,jk_ray),
!       and gam3(j,ji_ray,jk_ray), for all angular zones j between jmin and jmax.
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
!  jmin             : minimum y (angular) zone index.
!  jmax             : maximum y (angular) zone index.
!  rho              : shifted matter density array (g/cm**3).
!  t                : shifted matter matter temperature array (K).
!  ji_ray           : i (radial) index of a specific y (angular) ray
!  jk_ray           : k (azimuthal) index of a specific y (angular) ray
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  aesv(j,i,ji_ray,jk_ray)  : interpolated equation of state quantity i for angular zone j.
!  aesvd(j,i,ji_ray,jk_ray) : d(aesv(j,i))/d(rho).
!  aesvd(j,i,ji_ray,jk_ray  : d(aesv(j,i))/d(t).
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module,
!  eos_snc_y_module
!
!-----------------------------------------------------------------------
USE array_module, ONLY : ny

USE kind_module, ONLY : double
USE numerical_module, ONLY : one, epsilon

USE eos_snc_y_module, ONLY : aesv, aesvd, aesvt, gam1, gam2, gam3

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jmin          ! minimum y (angular) zone index
INTEGER, INTENT(in)              :: jmax          ! maximum y (angular) zone index
INTEGER, INTENT(in)              :: ji_ray        ! i (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray        ! k (azimuthal) index of a specific y (angular) ray

REAL(KIND=double), INTENT(in), DIMENSION(ny) :: rho ! shifted matter density array (g/cm**3)
REAL(KIND=double), INTENT(in), DIMENSION(ny) :: t   ! shifted matter matter temperature array (K)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! angular zone index

REAL(KIND=double)                :: gammat        ! temporary variable

DO j = jmin,jmax
  gam3(j,ji_ray,jk_ray) = aesvt(j,1,ji_ray,jk_ray)                      &
&                       / ( aesvt(j,2,ji_ray,jk_ray) * rho(j) + epsilon ) + one
  gam1(j,ji_ray,jk_ray) = ( t(j) * aesvt(j,1,ji_ray,jk_ray)             &
&                       * ( gam3(j,ji_ray,jk_ray) - one ) + rho(j)      &
&                       * aesvd(j,1,ji_ray,jk_ray) ) &
&                       / ( aesv(j,1,ji_ray,jk_ray) + epsilon )
  gammat                = ( rho(j) * aesvd(j,1,ji_ray,jk_ray)           &
&                       / ( gam3(j,ji_ray,jk_ray) - one ) + t(j)        &
&                       * aesvt(j,1,ji_ray,jk_ray) ) &
&                       / ( aesv(j,1,ji_ray,jk_ray) + epsilon )
  gam2(j,ji_ray,jk_ray) = gammat/( gammat - one )
END DO

RETURN
END SUBROUTINE gammaz_y
