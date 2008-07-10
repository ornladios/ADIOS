SUBROUTINE radhyd_to_eos_z_reset( nz, ki_ray, kj_ray, ij_ray_dim, k_ray_dim, &
& nnc, reset_comp_eos )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_eos_z_reset
!    Module:       radhyd_to_eos_z_reset
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To update the y-array EOS tablese.
!
!    Subprograms called:
!  reset_eos_z_inout : directs the recomputation, if necessary, of the z-array
!   EOS tables
!
!    Input arguments:
!  nz                : y-array extent
!  ki_ray            : x (radial) index of a specific z (azimthal) ray
!  kj_ray            : y (angular) index of a specific z (azimthal) ray
!  ij_ray_dim        : the number of radial zones on a processor before swapping with y
!  k_ray_dim         : the number of z-zones on a processor after swapping with z
!  nnc               : abundance array extent
!  reset_comp_eos    : composition reset flag
!
!    Output arguments:
!        none
!
!    Include files:
!  radial_ray_module, azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : kmin, kmax, nprint
USE azimuthal_ray_module, ONLY : idty_z, nse_z, rho_z, t_z, ye_z, ei_z, xn_z, &
& be_nuc_rep_z, a_nuc_rep_z, z_nuc_rep_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

CHARACTER(len=2), INTENT(in)     :: reset_comp_eos ! composition reset flag

INTEGER, INTENT(in)              :: nz             ! x-array extent
INTEGER, INTENT(in)              :: ki_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)              :: kj_ray         ! index denoting a specific radial ray
INTEGER, INTENT(in)              :: ij_ray_dim     ! the number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim      ! the number of z-zones on a processor after swapping with z
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

CALL reset_eos_z_inout( kmin, kmax, nz, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, ls, le, nnc, nprint, idty_z, nse_z, rho_z, t_z, ye_z, &
& ei_z, xn_z, be_nuc_rep_z, a_nuc_rep_z, z_nuc_rep_z, reset_comp_eos )

RETURN
END SUBROUTINE radhyd_to_eos_z_reset
