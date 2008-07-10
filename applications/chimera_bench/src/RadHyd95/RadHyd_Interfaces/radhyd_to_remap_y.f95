SUBROUTINE radhyd_to_remap_y( nx, ny, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
& i_radial, nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_remap_y
!    Module:       radhyd_to_remap_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To port angular_ray_module variables into, and updated variabbles
!       out of, the y-array remap modules via subroutine remap_y_inout.
!
!    Input arguments:
!  nx            : x-array extent
!  ny            : y-array extent
!  ji_ray        : x (radial) index of a specific y (angular) ray
!  jk_ray        : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim     : the number of radial zones on a processor after swapping with y
!  ik_ray_dim    : the number of z-zones on a processor before swapping with z
!  i_radial      : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez           : neutrino energy array extent
!  nnu           : neutrino flavor array extent
!  nnc           : composition array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  remap_y_inout : executes the remap along the y-direction
!
!    Include files:
!  radial_ray_module, angular_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : jmin, jmax, y_el, dy_cl, y_cl, y_ef, dy_cf, &
& y_cf, time, t_bounce, tb_dy_shift, rhobar
USE angular_ray_module, ONLY : rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, psi0_y, &
& xn_y, a_nuc_rep_y, z_nuc_rep_y, be_nuc_rep_y, flat_x_y

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ki_ray, kj_ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER                          :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! lower composition index
INTEGER                          :: le            ! upper composition index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                    = 1
le                    = nnc

!-----------------------------------------------------------------------
!  Transfer variables to remap arrays
!-----------------------------------------------------------------------

CALL remap_y_inout( nx, ny, jmin, jmax, ji_ray, jk_ray, j_ray_dim,      &
& ik_ray_dim, i_radial, nez, nnu, ls, le, nnc, y_el, dy_cl, y_cl, y_ef, &
& dy_cf, y_cf, rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, psi0_y, xn_y,     &
& a_nuc_rep_y, z_nuc_rep_y, be_nuc_rep_y, flat_x_y, time, t_bounce,     &
& tb_dy_shift, rhobar )

RETURN
END SUBROUTINE radhyd_to_remap_y
