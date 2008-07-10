SUBROUTINE radhyd_to_nu_stress_x( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& nx, nez, nnu)
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_stress_x
!    Module:       radhyd_to_nu_stress_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To port the radial_ray_module variables necessary for the computation
!       of the x-component of the neutrino stress into subroutine nu_stress_x_inout,
!       and port the computed stresses back.
!
!    Input arguments:
!  ij_ray            : j-index of a radial ray
!  ik_ray            : k-index of a radial ray
!  ij_ray_dim        : number of y-zones on a processor before swapping
!  ik_ray_dim        : number of z-zones on a processor before swapping
!  nx                : x-array extent
!  nez               : neutrino energy array extent
!  nnu               : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  nu_stress_x_inout : executes the calculation of the x-component of
!                       the neutrino stress
!
!    Include files:
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, x_ei, rho_c, rhs1_c, dc_e, psi0_c, &
& nu_str_c, nu_str_e, rhobar
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CALL nu_stress_x_inout( imin, imax, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& nx, nez, nnu, x_ei, rho_c, rhobar, rhs1_c, dc_e, psi0_c, nu_str_c, nu_str_e )

END SUBROUTINE radhyd_to_nu_stress_x
