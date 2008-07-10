SUBROUTINE radhyd_to_equilibrate_y( nx, ny, nez, nnu, ji_ray, jk_ray, &
& j_ray_dim, ik_ray_dim, i_radial )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_equilibrate_y
!    Module:       radhyd_to_equilibrate_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/20/05
!
!    Purpose:
!      To transfer variables equilibrate_y, which equilibrates matter and
!       neutrinos along an angular ray.
!
!    Subprograms called:
!      equilibrate_y
!
!    Input arguments:
!  nx         : x-array extent
!  ny         : y-array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  ji_ray     : i (radial) index of a specific angular ray
!  jk_ray     : k (azimuthal) index of a specific angular ray
!  j_ray_dim  : the number of radial zones on a processor after swapping with y
!  ik_ray_dim : the number of z-zones on a processor before swapping with z
!  i_radial   : the unshifted radial zone (angular ray) corresponding to j_ray
!
!    Output arguments:
!        none
!
!    Include files:
!  radial_ray_module, angular_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : jmin, jmax, nprint, rhobar, rho_equilibrate, &
& dtnph
USE angular_ray_module, ONLY : rho_y, t_y, ye_y, psi0_y, agr_y

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ji_ray        ! i (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray        ! k (azimuthal) index of a specific angular ray
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to j_ray

!-----------------------------------------------------------------------
!  Transfer variables to equilibrate_y
!-----------------------------------------------------------------------

CALL equilibrate_y( jmin, jmax, nx, ny, nez, nnu, ji_ray, jk_ray,       &
& j_ray_dim, ik_ray_dim, i_radial, rho_y, t_y, ye_y, psi0_y, rhobar,    &
& dtnph, rho_equilibrate, agr_y  )

RETURN
END SUBROUTINE radhyd_to_equilibrate_y
