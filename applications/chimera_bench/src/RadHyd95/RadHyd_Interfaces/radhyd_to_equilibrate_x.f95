SUBROUTINE radhyd_to_equilibrate_x( nx, nez, nnu, ij_ray_dim, ik_ray_dim, &
& ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_equilibrate_x
!    Module:       radhyd_to_equilibrate_x
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
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  ij_ray_dim : number of y-zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  ij_ray     : index denoting the j-index of a specific radial ray
!  ik_ray     : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, rho_c, t_c, ye_c, psi0_c, &
& dtnph, rho_equilibrate

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray with y
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray qith z

!-----------------------------------------------------------------------
!  Transfer variables to equilibrate_y
!-----------------------------------------------------------------------

CALL equilibrate_x( imin, imax, nx, nez, nnu, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, rho_c, t_c, ye_c, psi0_c, dtnph, rho_equilibrate )

RETURN
END SUBROUTINE radhyd_to_equilibrate_x
