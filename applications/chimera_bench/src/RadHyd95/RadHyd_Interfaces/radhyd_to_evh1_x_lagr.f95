SUBROUTINE radhyd_to_evh1_x_lagr( nx, ny, nz, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_evh1_x_lagr
!    Module:       radhyd_to_evh1_x_lagr
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the x-array Lagrangian hydro modules via subroutine
!       evh1_x_lagr_inout.
!
!    Input arguments:
!  nx                : x-array extent
!  ny                : y-array extent
!  nz                : z-array extent
!  ij_ray            : j-index of a radial ray
!  ik_ray            : k-index of a radial ray
!  ij_ray_dim        : number of y-zones on a processor before swapping
!  ik_ray_dim        : number of z-zones on a processor before swapping
!  nez               : neutrino energy array extent
!  nnu               : neutrino flavor array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  evh1_x_lagr_inout : executes the Lagrangian x-hydrodynamics
!
!    Include files:
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, x_ei, dx_ci, x_ci, x_el, dx_cl, &
& x_cl, rho_c, rho_l, t_c, ye_c, ei_c, e_v_c, u_c, u_e, v_c, w_c, u_l,    &
& v_l, w_l, nu_str_c, nu_str_e, dtime=>dtnph, dt, jdt, rhobar, flat_x,    &
& grav_x_c, grav_pot_c, grav_x_e, grav_pot_e, e_nu_c_bar, f_nu_e_bar,     &
& agr_c, agr_e, e_nu_c
     
IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

CALL evh1_x_lagr_inout( imin, imax, nx, ny, nz, ij_ray, ik_ray, ij_ray_dim,      &
& ik_ray_dim, nez, nnu, x_ei, dx_ci, x_ci, x_el, dx_cl, x_cl, rho_c, rho_l,      &
& t_c, ye_c, ei_c, e_v_c, u_c, u_e, v_c, w_c, u_l, v_l, w_l, nu_str_c, nu_str_e, &
& dtime, dt, jdt, rhobar, flat_x, grav_x_c, grav_pot_c, grav_x_e, grav_pot_e,    &
& e_nu_c_bar, f_nu_e_bar, agr_c, agr_e, e_nu_c )

END SUBROUTINE radhyd_to_evh1_x_lagr
