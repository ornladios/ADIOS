SUBROUTINE radhyd_to_remap_x_e( nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, &
& nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_remap_x_e
!    Module:       radhyd_to_remap_x_e
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the y-array remap modules via subroutine remap_x_inout.
!
!    Input arguments:
!  nx            : x-array extent
!  i_ray         : index denoting a specific radial ray
!  i_ray_dim     : number of radial rays assigned to a processor
!  nez           : neutrino energy array extent
!  nnu           : neutrino flavor array extent
!  nnc           : composition array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  remap_x_inout : executes the remap along the x-direction
!
!    Include files:
!  radial_ray_module, evh1_global
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, x_ci, x_el, dx_cl, x_cl, x_ef, &
& dx_cf, x_cf, rho_l, rho_c, t_c, ye_c, ei_c, e_v_c, u_l, v_l, w_l, u_c, &
& v_c, w_c, p_c, xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, e_bind_zn_c, &
& eb_c, fluxbe_c, grav_pot_c, grav_pot_c_i
USE evh1_global, ONLY : lagrangian

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)              :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer variables to remap_e arrays
!-----------------------------------------------------------------------

IF ( .not. lagrangian ) THEN
  CALL remap_x_e_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, &
&  ik_ray_dim, nnc, x_ci, x_el, dx_cl, x_cl, x_ef, dx_cf, x_cf, rho_l, &
&  rho_c, t_c, ye_c, ei_c, e_v_c, u_l, v_l, w_l, u_c, v_c, w_c, p_c, &
&  xn_c, a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, e_bind_zn_c, eb_c, &
&  fluxbe_c, grav_pot_c, grav_pot_c_i )
END IF

RETURN
END SUBROUTINE radhyd_to_remap_x_e
