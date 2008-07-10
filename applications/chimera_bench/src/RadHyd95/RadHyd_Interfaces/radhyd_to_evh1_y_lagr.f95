SUBROUTINE radhyd_to_evh1_y_lagr( nx, ny, nz, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, i_radial, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_evh1_y_lagr
!    Module:       radhyd_to_evh1_y_lagr
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the x-array Lagrangian hydro modules via subroutine
!       evh1_y_lagr_inout.
!
!    Input arguments:
!  nz                : x-array extent
!  ny                : y-array extent
!  nz                : z-array extent
!  ji_ray            : x (radial) index of a specific y (angular) ray
!  jk_ray            : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim         : the number of radial zones on a processor after swapping with y
!  ik_ray_dim        : the number of z-zones on a processor before swapping with z
!  i_radial          : the unshifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!  nez               : neutrino energy array extent
!  nnu               : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  evh1_y_lagr_inout : executes the Lagrangian y-hydrodynamics
!
!    Include files:
!  radial_ray_module, angular_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : jmin, jmax, x_ef, dx_cf, x_cf, y_ei, dy_ci, &
& y_ci, z_ei, dz_ci, z_ci, y_el, dy_cl, y_cl, dtime=>dtnph, rhobar
USE angular_ray_module, ONLY : rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y,     &
& nu_str_cy, nu_str_ey, dyphys_c, flat_y, agr_y, grav_y_cy, e_nu_y, dt_y, &
& jdt_y
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ji_ray, jk_ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 
CALL evh1_y_lagr_inout( jmin, jmax, nx, ny, nz, ji_ray, jk_ray, j_ray_dim,    &
& ik_ray_dim, i_radial, nez, nnu, x_cf, y_ei, dy_ci, y_ci, y_el, dy_cl, y_cl, &
& dyphys_c, rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, nu_str_cy, nu_str_ey,      &
& dtime, dt_y, jdt_y, rhobar, flat_y, agr_y, grav_y_cy, e_nu_y )

END SUBROUTINE radhyd_to_evh1_y_lagr
