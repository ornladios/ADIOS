SUBROUTINE radhyd_to_nu_stress_z( ki_ray, kj_ray, ij_ray_dim, k_ray_dim, &
& i_radial, j_radial, nx, ny, nz, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_stress_z
!    Module:       radhyd_to_nu_stress_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/22/07
!
!    Purpose:
!      To port angular_ray_module variables into, and updated variabbles
!       out of, the nu_stress_y module via subroutine
!       evh1_x_lagr_inout.
!
!    Input arguments: 
!  ki_ray     : x (radial) index of a specific z (azimuthal) ray
!  kj_ray     : y (azimuthal) index of a specific z (azimuthal) ray
!  ij_ray_dim : the number of y-zones on a processor before swapping with y
!  k_ray_dim  : the number of radial zones on a processor after swapping with z
!  i_radial   : the unshifted radial zone (angular ray) corresponding to ki_ray
!  j_radial   : the shifted radial zone (angular ray) corresponding to ki_ray
!  nx         : x-array extent
!  ny         : y-array extent
!  nz         : z-array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  nu_stress_z_inout : supervises the computation of the neutrino y_stress
!
!    Include files:
!  radial_ray_module, azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : kmin, kmax, x_ef, rhobar, y_ei, z_ei
USE azimuthal_ray_module, ONLY : rho_z, psi0_z, nu_str_cz, nu_str_ez
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (azimuthal) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to j_ray
INTEGER, INTENT(in)              :: j_radial      ! the shifted radial zone corresponding to j_ray
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! x-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CALL nu_stress_z_inout( kmin, kmax, ki_ray, kj_ray, ij_ray_dim, k_ray_dim,   &
& i_radial, j_radial, nx, ny, nz, nez, nnu, x_ef, rhobar, y_ei, z_ei, rho_z, &
& psi0_z, nu_str_cz, nu_str_ez )

END SUBROUTINE radhyd_to_nu_stress_z
