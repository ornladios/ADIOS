SUBROUTINE radhyd_to_eos_x_reset( nx, ij_ray_dim, ik_ray_dim, ij_ray, &
& ik_ray, nnc, reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_eos_x_reset
!    Module:       radhyd_to_eos_x_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/05
!
!    Purpose:
!      To update the x-array EOS tables.
!
!    Input arguments:
!  nx                : x-array extent
!  ij_ray_dim        : number of y-zones on a processor before swapping
!  ik_ray_dim        : number of z-zones on a processor before swapping
!  ij_ray            : j-index of a radial ray
!  ik_ray            : k-index of a radial ray
!  nnc               : abundance array extent
!  nnc               : abundance array extent
!  reset_comp_eos    : composition reset flag
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  reset_eos_x_inout : updates the equation of state variables
!
!    Include files:
!  kind_module
!  radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module

USE radial_ray_module, ONLY : imin, imax, nprint, idty, nse_c, rho_c, t_c, &
& ye_c, ei_c, xn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos ! composition reset flag

INTEGER, INTENT(in)              :: nx             ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim     ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim     ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)              :: ij_ray         ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray         ! k-index of a radial ray
INTEGER, INTENT(in)              :: nnc            ! abundance array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! minimum abundance array index
INTEGER                          :: le            ! maximum abundance array index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                 = 1
le                 = nnc

!-----------------------------------------------------------------------
!  Transfer variables to mgfld arrays and reset_tables
!-----------------------------------------------------------------------

CALL reset_eos_x_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, ls, le, nnc, nprint, idty, rho_c, t_c, ye_c, ei_c, xn_c, &
& be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, nse_c, reset_comp_eos )

RETURN
END SUBROUTINE radhyd_to_eos_x_reset
