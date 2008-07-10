 SUBROUTINE store_int_angular_var( ji_ray, jk_ray, ny )
!-----------------------------------------------------------------------
!
!    File:         store_int_angular_var
!    Module:       store_int_angular_var
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To store coordinate and state variables before the angular
!       Lagrangian hydro step.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ji_ray : x (radial) index of a specific y (angular) ray
!  jk_ray : z (azimuthal) index of a specific y (angular) ray
!  ny     : y-array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  angular_ray_module
!
!-----------------------------------------------------------------------

USE angular_ray_module, ONLY : rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, &
& rho_yi, t_yi, ye_yi, ei_yi, u_yi, v_yi, w_yi

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ji_ray          ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray          ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: ny              ! y-array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!         \\\\\ STORE INITIAL VARIABLES FOR THE Y-SWEEP /////
!
!-----------------------------------------------------------------------

rho_yi(1:ny,ji_ray,jk_ray) = rho_y(1:ny,ji_ray,jk_ray)
t_yi  (1:ny,ji_ray,jk_ray) = t_y  (1:ny,ji_ray,jk_ray)
ye_yi (1:ny,ji_ray,jk_ray) = ye_y (1:ny,ji_ray,jk_ray)
ei_yi (1:ny,ji_ray,jk_ray) = ei_y (1:ny,ji_ray,jk_ray)
u_yi  (1:ny,ji_ray,jk_ray) = u_y  (1:ny,ji_ray,jk_ray)
v_yi  (1:ny,ji_ray,jk_ray) = v_y  (1:ny,ji_ray,jk_ray)
w_yi  (1:ny,ji_ray,jk_ray) = w_y  (1:ny,ji_ray,jk_ray)

RETURN
END SUBROUTINE store_int_angular_var
