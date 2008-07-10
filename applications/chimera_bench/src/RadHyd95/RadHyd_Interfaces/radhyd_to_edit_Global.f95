SUBROUTINE radhyd_to_edit_Global( nx, nez, nnu, nnc, ij_ray_dim, ny, &
& ik_ray_dim, nz )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_edit_Global
!    Module:       radhyd_to_edit_Global
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/06/06
!
!    Purpose:
!      To port variables into edit_Global_inout for global edits.
!
!    Subprograms called:
!  edit_Global_inout : calculates quantities for the global edit dumps
!
!    Input arguments:
!  nx           : x_array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  nnc          : composition array extent
!  ij_ray_dim   : number of y-zones on a processor before swapping with y
!  ny           : y_array extent
!  ik_ray_dim   : number of z-zones on a processor before swapping with z
!  nz           : z_array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  radhyd_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, x_ef, &
& x_cf, dx_cf, y_ef, u_c, v_c, w_c, rho_c, t_c, ye_c, psi0_c, psi1_e, &
& unu_c, dunu_c, unue_e, dunue_e, dtnph, time, t_bounce, ncycle, xn_c, &
& nse_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, d_omega
     
IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc              ! composition array extent
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ny               ! y-array extent
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nz               ! z_array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ TRANSFER VARIABLES TO EDIT_2D /////
!
!-----------------------------------------------------------------------

CALL edit_Global_inout( imin, imax, nx, nez, nnu, jmin, jmax, ij_ray_dim, &
& ny, kmin, kmax, ik_ray_dim, nz, x_ef, x_cf, dx_cf, u_c, v_c, w_c, rho_c, &
& t_c, ye_c, psi0_c, psi1_e, unu_c, dunu_c, unue_e, dunue_e, dtnph, time, &
& t_bounce, ncycle, nnc, xn_c, nse_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, &
& d_omega )

RETURN
END SUBROUTINE radhyd_to_edit_Global
