SUBROUTINE radhyd_to_remap_z( nx, nz, ki_ray, kj_ray, ij_ray_dim, k_ray_dim, &
& i_radial, nez, nnu, nnc )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_remap_z
!    Module:       radhyd_to_remap_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/12/07
!
!    Purpose:
!      To port radial_ray_module variables into, and updated variabbles
!       out of, the y-array Lagrangian hydro modules via subroutine
!       evh1_y_lagr_inout.
!
!    Input arguments:
!  nx            : x-array extent
!  nz            : z-array extent
!  ki_ray        : x (radial) index of a specific z (azimuthal) ray
!  kj_ray        : y (azimuthal) index of a specific z (azimuthal) ray
!  ij_ray_dim    : the number of radial zones on a processor before swapping with y
!  k_ray_dim     : the number of z-zones on a processor after swapping with z
!  i_radial      : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez           : neutrino energy array extent
!  nnu           : neutrino flavor array extent
!  nnc           : composition array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  remap_z_inout : executes the remap along the z-direction
!
!    Include files:
!  radial_ray_module, azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : kmin, kmax, z_el, dz_cl, z_cl, z_ef, dz_cf, &
& z_cf, time, t_bounce, tb_dy_shift, rhobar
USE azimuthal_ray_module, ONLY : rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, psi0_z, &
& xn_z, a_nuc_rep_z, z_nuc_rep_z, be_nuc_rep_z, flat_x_z

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: nz            ! y-array extent
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (azimuthal) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ki_ray, kj_ray
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER                          :: nnc           ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! lower composition index
INTEGER                          :: le            ! upper composition index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                    = 1
le                    = nnc

!-----------------------------------------------------------------------
!  Transfer variables to remap arrays
!-----------------------------------------------------------------------

CALL remap_z_inout( nx, nz, kmin, kmax, ki_ray, kj_ray, ij_ray_dim,    &
& k_ray_dim, i_radial, nez, nnu, ls, le, nnc, z_el, dz_cl, z_cl, z_ef, &
& dz_cf, z_cf, rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, psi0_z, xn_z,    &
& a_nuc_rep_z, z_nuc_rep_z, be_nuc_rep_z, flat_x_z, time, t_bounce,    &
& tb_dy_shift, rhobar )

RETURN
END SUBROUTINE radhyd_to_remap_z
