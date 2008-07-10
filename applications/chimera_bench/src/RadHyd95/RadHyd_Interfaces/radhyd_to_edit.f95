SUBROUTINE radhyd_to_edit( ij_ray_min, ij_ray_max, ij_ray_dim, ik_ray_min, &
& ik_ray_max, ik_ray_dim, i_editp, nx, ny, nz, nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_edit
!    Module:       radhyd_to_edit
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To load variables into 3-d arrays for porting into MGFLD edit
!         via subroutine edit_in.
!
!    Subprograms called:
!        edit_in
!
!    Input arguments:
!  ij_ray_min : minimum j-index of radial ray
!  ij_ray_max : maximum j-index of radial ray
!  ij_ray_dim : number of y-zones on a processor before swapping
!  ik_ray_min : minimum k-index of radial ray
!  ik_ray_max : maximum k-index of radial ray
!  ik_ray_dim : number of z-zones on a processor before swapping
!  i_editp    : edit parameter
!  nx         : x_array extent
!  ny         : y_array extent
!  nz         : z_array extent
!  nez        : neutrino energy array extent
!  nnu        : neutrino flavor array extent
!  nnc        : neutrino abundance array extent
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module
!  edit_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, rho_ci, &
& rho_c, t_ci, t_c, ye_ci, ye_c, u_c, v_c, w_c, x_ef, psi0_c, psi1_e, xn_c, &
& be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, nse_c, nedc, nedmi, nedma, nedh, &
& nedps, nedu, nedy, nedsc, nedn, nedng, dtnph, time, ncycle, t_bounce, &
& rhobar, nu_str_e, nu_str_c, agr_e, agr_c, grav_x_e, grav_x_c, grav_y_c, &
& grav_pot_c, mass_ns, vel_ns, d_omega
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: ij_ray_min       ! minimum j-index of radial ray
INTEGER, INTENT(in)              :: ij_ray_max       ! maximum j-index of radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_min       ! minimum k-index of radial ray
INTEGER, INTENT(in)              :: ik_ray_max       ! maximum k-index of radial ray
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: i_editp          ! edit parameter
INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: ny               ! y-array extent
INTEGER, INTENT(in)              :: nz               ! z-array extent
INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc              ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: ij_ray           ! j-index of a radial ray
INTEGER                          :: ik_ray           ! k-index of a radial ray

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!             ||||| BEGIN LOOP OVER RADIAL ARRAYS |||||
!
!-----------------------------------------------------------------------

DO ik_ray = ik_ray_min,ik_ray_max
  DO ij_ray = ij_ray_min,ij_ray_max

!-----------------------------------------------------------------------
!
!            \\\\\ TRANSFER VARIABLES TO MGFLD EDIT /////
!
!-----------------------------------------------------------------------

    CALL edit_in( imin, imax, jmin, jmax, kmin, kmax, nx, ny, nz, ij_ray,  &
&    ik_ray, ij_ray_dim, ik_ray_dim,  nez, nnu, nnc, rho_ci, rho_c, t_ci,  &
&    t_c, ye_ci, ye_c, rhobar, x_ef, u_c, v_c, w_c, psi0_c, psi1_e, dtnph, &
&    time, i_editp, ncycle, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c,  &
&    nu_str_e, nu_str_c, agr_e, agr_c, grav_x_e, grav_x_c, grav_y_c,       &
&    grav_pot_c, mass_ns, vel_ns, nse_c, nedc, nedmi, nedma, nedh, nedps,  &
&    nedu, nedy, nedsc, nedn, nedng, d_omega, first )

!-----------------------------------------------------------------------
!
!           \\\\\ TRANSFER VARIABLES FROM MGFLD EDIT /////
!
!-----------------------------------------------------------------------

    CALL edit_out( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, nedc, nedmi, &
&    nedma, nedh, nedps, nedu, nedy, nedsc, nedn, nedng, t_bounce )

!-----------------------------------------------------------------------
!
!              ||||| END LOOP OVER RADIAL ARRAYS |||||
!
!-----------------------------------------------------------------------

  END DO ! ij_ray = ij_ray_min,ij_ray_max
END DO ! ik_ray = ik_ray_min,ik_ray_max

first                  = .false.

RETURN
END SUBROUTINE radhyd_to_edit
