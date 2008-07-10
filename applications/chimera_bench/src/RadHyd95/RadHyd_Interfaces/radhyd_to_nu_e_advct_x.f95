SUBROUTINE radhyd_to_nu_e_advct_x( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
& ik_ray, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_nu_e_advct_x
!    Module:       radhyd_to_nu_e_advct_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the neutrino energy advection modules via subroutine
!       nu_energy_advct_inout_x.
!
!    Subprograms called:
!  nu_energy_advct_inout_x : transfers variables to and from neutrino x-advection modules
!
!    Input arguments:
!  nx         : x-array extent
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_dim : number of z-zones on a processor before swapping
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, nprint, ncycle, dtime=>dtnph, &
& time, rho_ci, rho_c, x_ei, x_el, t_c, ye_c, u_c, psi0_c, psi1_e, dt,  &
& jdt, agr_c

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!  Transfer variables to nu_energy_advct_inout_x
!-----------------------------------------------------------------------

CALL nu_energy_advct_inout_x( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nez, nnu, nprint, rho_ci, rho_c, x_ei, x_el, t_c, ye_c, u_c, &
& psi0_c, psi1_e, dtime, agr_c , dt, jdt)

RETURN
END SUBROUTINE radhyd_to_nu_e_advct_x
