SUBROUTINE radhyd_to_eos_y_reset( ny, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
& nnc, reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_eos_y_reset
!    Module:       radhyd_to_eos_y_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/20/05
!
!    Purpose:
!      To update the y-array EOS tablese.
!
!    Subprograms called:
!  reset_eos_y_inout : directs the recomputation, if necessary, of the y-array
!   EOS tables
!
!    Input arguments:
!  ny                : y-array extent
!  ji_ray            : x (radial) index of a specific y (angular) ray
!  jk_ray            : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim         : the number of radial zones on a processor after swapping with y
!  ik_ray_dim        : the number of z-zones on a processor before swapping with z
!  nnc               : abundance array extent
!  reset_comp_eos    : composition reset flag
!
!    Output arguments:
!        none
!
!    Include files:
!  radial_ray_module, angular_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : jmin, jmax, nprint
USE angular_ray_module, ONLY : idty_y, nse_y, rho_y, t_y, ye_y, ei_y, xn_y, &
& be_nuc_rep_y, a_nuc_rep_y, z_nuc_rep_y

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos ! composition reset flag

INTEGER, INTENT(in)              :: ny             ! x-array extent
INTEGER, INTENT(in)              :: j_ray_dim      ! number of rays assigned to a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim     ! the number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: ji_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)              :: jk_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)              :: nnc            ! abundance array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls             ! minimum abundance array index
INTEGER                          :: le             ! maximum abundance array index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                 = 1
le                 = nnc

!-----------------------------------------------------------------------
!  Transfer variables to reset_tables
!-----------------------------------------------------------------------

CALL reset_eos_y_inout( jmin, jmax, ny, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, ls, le, nnc, nprint, idty_y, nse_y, rho_y, t_y, ye_y, &
& ei_y, xn_y, be_nuc_rep_y, a_nuc_rep_y, z_nuc_rep_y, reset_comp_eos )

RETURN
END SUBROUTINE radhyd_to_eos_y_reset
