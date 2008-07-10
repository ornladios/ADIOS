SUBROUTINE radhyd_to_transpose_z_x( nx, ij_ray_dim, ik_ray_dim, nz, k_ray_dim, &
& nez, nnu, nnc, n_proc, n_proc_z )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_transpose_z_x
!    Module:       radhyd_to_transpose_z_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To transpose certain of the y-array variables in azimuthal_ray_module
!       and load them into radial_ray_module.
!
!    Input arguments:
!  nx         : x-array extent
!  ij_ray_dim : number of radial zones on a processor before swapping with y
!  ik_ray_dim : number of z-zones on a processor before swapping with z
!  nz         : z-array extent
!  k_ray_dim  : number of radial zones on a processor after swapping with z
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  nnc        : composition array dimension
!  n_proc     : number of processors assigned to the run
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  transpose_z_x : performs the actual variable transpose
!
!    Include files:
!  radial_ray_module, azimuthal_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, kmin, kmax, rho_c, t_c, ye_c, ei_c, &
& x_cf, y_ef, dy_cf, y_cf, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c, &
& z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_z_x
USE azimuthal_ray_module, ONLY : rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, psi0_z, &
& psi1_z, xn_z, a_nuc_rep_z, z_nuc_rep_z, be_nuc_rep_z, nse_z, flat_z
USE shock_module, ONLY : pq_z, pqz_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc           ! composition array extent
INTEGER, INTENT(in)              :: n_proc        ! number of processors assigned to the run
INTEGER, INTENT(in)              :: n_proc_z      ! number of processors assigned to the z-indices of rays

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! inner logical composition zone
INTEGER                          :: le            ! outer logical composition zone
INTEGER                          :: my_id = 0     ! outer logical composition zone

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                    = 1
le                    = nnc

!-----------------------------------------------------------------------
!  Transpose variables.
!-----------------------------------------------------------------------

CALL transpose_z_x(imin, imax, nx, ij_ray_dim, ik_ray_dim, kmin, kmax, &
& nz, k_ray_dim, nez, nnu, ls, le, nnc, n_proc, n_proc_z, rho_c, t_c,  &
& ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c,        &
& z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_z_x, pq_z, rho_z, t_z, ye_z,  &
& ei_z, u_z, v_z, w_z, psi0_z, psi1_z, xn_z, a_nuc_rep_z, z_nuc_rep_z, &
& be_nuc_rep_z, nse_z, flat_z, pqz_x )

RETURN
END SUBROUTINE radhyd_to_transpose_z_x
