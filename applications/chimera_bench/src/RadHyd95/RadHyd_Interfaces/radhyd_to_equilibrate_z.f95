SUBROUTINE radhyd_to_equilibrate_z( nx, nz, nez, nnu, ki_ray, kj_ray, &
& ij_ray_dim, k_ray_dim, i_radial )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_equilibrate_z
!    Module:       radhyd_to_equilibrate_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/23/07
!
!    Purpose:
!      To transfer variables to equilibrate_z, which equilibrates matter
!       and neutrinos along an azimuthal ray.
!
!    Subprograms called:
!      equilibrate_y
!
!    Input arguments:
!  nx         : x-array extent
!  nz         : z-array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  ki_ray     : x (radial) index of a specific z (azimuthal) ray
!  kj_ray     : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim : the number of radial zones on a processor before swapping with y
!  k_ray_dim  : the number of z-zones on a processor after swapping with z
!  i_radial   : the unshifted radial zone (angular ray) corresponding to j_ray
!
!    Output arguments:
!        none
!
!    Include files:
!  radial_ray_module, azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : kmin, kmax, nprint, rhobar, rho_equilibrate, &
& dtnph
USE azimuthal_ray_module, ONLY : rho_z, t_z, ye_z, psi0_z, agr_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: k_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to j_ray

!-----------------------------------------------------------------------
!  Transfer variables to equilibrate_y
!-----------------------------------------------------------------------

CALL equilibrate_z( kmin, kmax, nx, nz, nez, nnu, ki_ray, kj_ray,       &
& ij_ray_dim, k_ray_dim, i_radial, rho_z, t_z, ye_z, psi0_z, rhobar,    &
& dtnph, rho_equilibrate, agr_z )

RETURN
END SUBROUTINE radhyd_to_equilibrate_z
