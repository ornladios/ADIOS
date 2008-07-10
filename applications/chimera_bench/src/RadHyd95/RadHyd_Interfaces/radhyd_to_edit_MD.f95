SUBROUTINE radhyd_to_edit_MD( nx, nez, nnu, ij_ray_dim, ny, ik_ray_dim, &
& nz, nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_edit_MD
!    Module:       radhyd_to_edit_MD
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!      To port variables into edit_MD_inout for M-D edits.
!
!    Subprograms called:
!  editMD_inout : calculates quantities for the M-D edit dumps
!
!    Input arguments:
!  nx           : x_array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  ij_ray_dim   : number of y-zones on a processor before swapping with y
!  ny           : y (angular) array dimension
!  ik_ray_dim   : number of z-zones on a processor before swapping with z
!  nz           : z (azimuthal) array dimension
!  nnc          : composition array extent
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

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, x_ef,      &
& x_cf, y_cf, z_cf, u_c, v_c, w_c, rho_c, t_c, ye_c, rhobar, psi0_c, psi1_e, &
& e_nu_c, f_nu_e, unu_c, dunu_c, unue_e, dunue_e, dtnph, time, t_bounce,     &
& xn_c, nse_c, ncycle, grav_x_e, grav_x_c, grav_y_c, grav_z_c, gtot_pot_c
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino flavor array extent
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ny               ! y (angular) array dimension
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nz               ! z (azimuthal) array dimension
INTEGER, INTENT(in)              :: nnc              ! composition array extent

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             \\\\\ TRANSFER VARIABLES TO EDIT_2D /////
!
!-----------------------------------------------------------------------

CALL edit_MD_inout( imin, imax, nx, nez, nnu, jmin, jmax, ij_ray_dim, ny,   &
& kmin, kmax, ij_ray_dim, nz, x_ef, x_cf, y_cf, z_cf, u_c, v_c, w_c, rho_c, &
& t_c, ye_c, rhobar, psi0_c, psi1_e, e_nu_c, f_nu_e, unu_c, dunu_c, unue_e, &
& dunue_e, dtnph, time, t_bounce, xn_c, nse_c, nnc, ncycle, grav_x_e,       &
& grav_x_c, grav_y_c, grav_z_c, gtot_pot_c )

RETURN
END SUBROUTINE radhyd_to_edit_MD
