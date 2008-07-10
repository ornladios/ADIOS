SUBROUTINE radhyd_to_transpose_x_z( nx, ij_ray_dim, ik_ray_dim, nz, &
& k_ray_dim, ny, nez, nnu, nnc, n_proc, n_proc_y, n_proc_z )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_transpose_x_z
!    Module:       radhyd_to_transpose_x_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To transpose certain of the x-array variables in radial_ray_module
!       and load them into azimuthal_ray_module.
!
!    Input arguments:
!  nx            : x-array extent
!  ij_ray_dim    : number of y-zones on a processor before swapping with y
!  ik_ray_dim    : number of z-zones on a processor before swapping with z
!  nz            : z-array extent
!  k_ray_dim     : the number of radial zones on a processor after swapping with z
!  ny            : y-array extent
!  nez           : neutrino energy array dimension
!  nnu           : neutrino flavor array dimension
!  nnc           : composition array dimension
!  n_proc        : number of processors assigned to the run
!  n_proc_y      : number of processors assigned to the y-indices of rays
!  n_proc_z      : number of processors assigned to the z-indices of rays
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  transpose_x_z : transposes variables for the z-sweeps
!
!    Include files:
!  radial_ray_module, azimuthal_ray_module, parallel_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, kmin, kmax, rho_c, &
& t_c, ye_c, ei_c, x_cf, y_cf, dz_cf, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, &
& a_nuc_rep_c, z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_x, agr_c, grav_z_c, &
& sin_theta, e_nu_c
USE azimuthal_ray_module, ONLY : rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, &
& psi0_z, psi1_z, xn_z, a_nuc_rep_z, z_nuc_rep_z, be_nuc_rep_z, dzphys_c, &
& nse_z, flat_x_z, agr_z, grav_z_cz, e_nu_z
USE parallel_module, ONLY : myid_z, myid_y
USE edit_module, ONLY : nlog

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: k_ray_dim     ! the number of radial zones on a processor after swapping with zi_ray_dim
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)              :: nnc           ! composition array extent
INTEGER, INTENT(in)              :: n_proc        ! number of processors assigned to the run
INTEGER, INTENT(in)              :: n_proc_y      ! number of processors assigned to the y-indices of rays
INTEGER, INTENT(in)              :: n_proc_z      ! number of processors assigned to the z-indices of rays

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! inner logical composition zone
INTEGER                          :: le            ! outer logical composition zone
INTEGER                          :: j             ! y-array index
INTEGER                          :: k             ! z-array index
INTEGER                          :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                         = 1
le                         = nnc

!-----------------------------------------------------------------------
!  Transpose variables
!-----------------------------------------------------------------------

CALL transpose_x_z(imin, imax, nx, kmin, kmax, nz, ij_ray_dim, ik_ray_dim,  &
& k_ray_dim, ny, nez, nnu, ls, le, nnc, n_proc, n_proc_y, n_proc_z, rho_c,  &
& t_c, ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c,        &
& z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_x, agr_c, e_nu_c, grav_z_c, rho_z, &
& t_z, ye_z, ei_z, u_z, v_z, w_z, nse_z, psi0_z, psi1_z, xn_z, a_nuc_rep_z, &
& z_nuc_rep_z, be_nuc_rep_z, flat_x_z, agr_z, grav_z_cz, e_nu_z )

!-----------------------------------------------------------------------
!  Compute physical z-coordinate differences
!-----------------------------------------------------------------------

DO ki_ray = 1,k_ray_dim
  DO j = 1,ij_ray_dim
    DO k = kmin,kmax
      dzphys_c(k,j,ki_ray) = x_cf(k_ray_dim * myid_z + ki_ray) * dz_cf(k) &
&                          * sin_theta(ij_ray_dim * myid_y + j)
    END DO ! k = kmin,kmax
  END DO ! j = jmin,jmax
END DO ! ji_ray = 1,j_ray_dim

RETURN
END SUBROUTINE radhyd_to_transpose_x_z
