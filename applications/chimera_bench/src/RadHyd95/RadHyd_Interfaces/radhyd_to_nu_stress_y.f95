SUBROUTINE radhyd_to_nu_stress_y( ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
& i_radial, j_radial, nx, ny, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_stress_y
!    Module:       radhyd_to_nu_stress_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To port angular_ray_module variables into, and updated variabbles
!       out of, the nu_stress_y module via subroutine
!       evh1_x_lagr_inout.
!
!    Input arguments:
!  ji_ray     : i (radial) index of a specific angular ray
!  jk_ray     : k (azimuthal) index of a specific angular ray
!  j_ray_dim  : the number of radial zones on a processor after swapping with y
!  ik_ray_dim : the number of z-zones on a processor before swapping with z
!  i_radial   : the unshifted radial zone (angular ray) corresponding to j_ray
!  j_radial   : the shifted radial zone (angular ray) corresponding to j_ray
!  nx         : x-array extent
!  ny         : y-array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  nu_stress_x_inout
!
!    Include files:
!  radial_ray_module, angular_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : jmin, jmax, x_ef, rhobar, y_ei
USE angular_ray_module, ONLY : rho_y, psi0_y, nu_str_cy, nu_str_ey
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ji_ray        ! i (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray        ! k (azimuthal) index of a specific angular ray
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to j_ray
INTEGER, INTENT(in)              :: j_radial      ! the shifted radial zone corresponding to j_ray
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CALL nu_stress_y_inout( jmin, jmax, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
& i_radial, j_radial, nx, ny, nez, nnu, x_ef, rhobar, y_ei, rho_y, psi0_y, &
& nu_str_cy, nu_str_ey )

END SUBROUTINE radhyd_to_nu_stress_y
