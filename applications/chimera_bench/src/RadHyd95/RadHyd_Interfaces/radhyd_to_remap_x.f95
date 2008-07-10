SUBROUTINE radhyd_to_remap_x( nx, ij_ray_dim, ik_ray_dim, ij_ray, ik_ray, &
& nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_remap_x
!    Module:       radhyd_to_remap_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/14/07
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the x-array remap modules via subroutine remap_x_inout.
!
!    Input arguments:
!  nx            : x-array extent
!  ij_ray_dim    : number of y-zones on a processor before swapping with y
!  ik_ray_dim    : number of z-zones on a processor before swapping with z
!  ij_ray        : index denoting the j-index of a specific radial ray
!  ik_ray        : index denoting the k-index of a specific radial ray
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

USE radial_ray_module, ONLY : imin, imax, x_el, dx_cl, x_cl, x_ef, dx_cf, &
& x_cf, rho_c, t_c, ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, &
& a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, agr_c, e_bind_zn_c, eb_c, fluxbe_c
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
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER                          :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! inner logical composition zone
INTEGER                          :: le            ! outer logical composition zone

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                    = 1
le                    = nnc

!-----------------------------------------------------------------------
!  Transfer variables to remap arrays
!-----------------------------------------------------------------------

IF ( .not. lagrangian ) THEN
  CALL remap_x_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&  nez, nnu, ls, le, nnc, x_el, dx_cl, x_cl, x_ef, dx_cf, x_cf, rho_c, t_c, &
&  ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c, z_nuc_rep_c, &
&  be_nuc_rep_c, agr_c, e_bind_zn_c, eb_c, fluxbe_c )
END IF

RETURN
END SUBROUTINE radhyd_to_remap_x
