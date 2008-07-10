SUBROUTINE solid_angle
!-----------------------------------------------------------------------
!
!    File:         solid_angle
!    Module:       solid_angle
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/12/06
!
!    Purpose:
!      To compute the solid angles, normalized to 4pi, subtended by the
!       radial rays, and the unnormalized sum.
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
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : frpi
USE physcnst_module, ONLY : pi

USE radial_ray_module, ONLY : jmin, jmax, kmin, kmax, y_ei, y_ci, z_ei, &
& z_ci, d_omega, omega, cos_theta, sin_theta, cos_phi, sin_phi, ndim

IMPLICIT none

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                    :: j              ! y-array index
INTEGER                                    :: k              ! z-array index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  If ndim = 1, set d_omega = omega = 4pi
!-----------------------------------------------------------------------

IF ( ndim == 1 ) THEN
  d_omega(1,1)                 = frpi
  omega                        = frpi

!-----------------------------------------------------------------------
!  If ndim = 2, compute d_omega(j,1)
!-----------------------------------------------------------------------

ELSE IF ( ndim == 2 ) THEN

   d_omega(jmin:jmax,1)        = DCOS(y_ei(jmin:jmax)) - DCOS(y_ei(jmin+1:jmax+1))
   omega                       = DCOS(y_ei(jmin)) - DCOS(y_ei(jmax+1))

!-----------------------------------------------------------------------
!  Compute cosines and sines of theta
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    cos_theta(j)               = DCOS(y_ci(j))
    sin_theta(j)               = DSIN(y_ci(j))
  END DO ! j

!-----------------------------------------------------------------------
!  If ndim = 3, compute d_omega(j,k)
!-----------------------------------------------------------------------

ELSE IF ( ndim == 3 ) THEN

!-----------------------------------------------------------------------
!  Compute solid angles subtended by radial rays, normalized to 4pi
!-----------------------------------------------------------------------

  DO k = kmin,kmax
    d_omega(jmin:jmax,k)       = ( DCOS(y_ei(jmin:jmax)) - DCOS(y_ei(jmin+1:jmax+1)) ) &
&                              * ( z_ei(k+1) - z_ei(k) )
  END DO ! k
  omega                        = ( DCOS(y_ei(jmin)) - DCOS(y_ei(jmax+1)) )             &
&                              * ( z_ei(kmax+1) - z_ei(kmin) )
  d_omega(jmin:jmax,kmin:kmax) = frpi * d_omega(jmin:jmax,kmin:kmax)/omega

!-----------------------------------------------------------------------
!  Compute cosines and sines of theta
!-----------------------------------------------------------------------

  DO j = jmin,jmax
    cos_theta(j)               = DCOS(y_ci(j))
    sin_theta(j)               = DSIN(y_ci(j))
  END DO ! j

!-----------------------------------------------------------------------
!  Compute cosines and sines of phi
!-----------------------------------------------------------------------

  DO k = kmin,kmax
    cos_phi(k)                 = DCOS(z_ci(k))
    sin_phi(k)                 = DSIN(z_ci(k))
  END DO ! k

END IF ! ndim


RETURN
END SUBROUTINE solid_angle
