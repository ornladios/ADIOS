SUBROUTINE radhyd_to_evh1_z_lagr( nx, ny, nz, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, i_radial, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_evh1_z_lagr
!    Module:       radhyd_to_evh1_z_lagr
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the x-array Lagrangian hydro modules via subroutine
!       evh1_x_lagr_inout.
!
!    Input arguments:
!  nz                : x-array extent
!  ny                : y-array extent
!  nz                : z-array extent
!  ki_ray            : x (radial) index of a specific z (azimuthal) ray
!  kj_ray            : y (azimuthal) index of a specific z (azimuthal) ray
!  ij_ray_dim        : the number of radial zones on a processor before swapping with y
!  k_ray_dim         : the number of z-zones on a processor after swapping with z
!  i_radial          : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez               : neutrino energy array extent
!  nnu               : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  evh1_z_lagr_inout : executes the Lagrangian z-hydrodynamics
!
!    Include files:
!  radial_ray_module, azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : kmin, kmax, x_cf, z_ei, dz_ci, z_ci, z_el, &
& dz_cl, z_cl, dtime=>dtnph, rhobar, cos_theta, sin_theta
USE azimuthal_ray_module, ONLY : rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z,   &
& nu_str_cz, nu_str_ez, dzphys_c, flat_z, agr_z, grav_z_cz, e_nu_z, dt_z, &
& jdt_z
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (azimuthal) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ki_ray, kj_ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
 
CALL evh1_z_lagr_inout( kmin, kmax, nx, ny, nz, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, i_radial, nez, nnu, x_cf, z_ei, dz_ci, z_ci, z_el, dz_cl,      &
& z_cl, dzphys_c, rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, nu_str_cz,         &
& nu_str_ez, dtime, dt_z, jdt_z, rhobar, flat_z, agr_z, grav_z_cz,          &
& cos_theta, sin_theta, e_nu_z )

END SUBROUTINE radhyd_to_evh1_z_lagr
