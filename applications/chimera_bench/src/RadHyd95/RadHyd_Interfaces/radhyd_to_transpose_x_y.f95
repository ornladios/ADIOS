SUBROUTINE radhyd_to_transpose_x_y( nx, ij_ray_dim, ik_ray_dim, ny, &
& j_ray_dim, nz, nez, nnu, nnc, n_proc, n_proc_y, n_proc_z )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_transpose_x_y
!    Module:       radhyd_to_transpose_x_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To transpose certain of the x-array variables in radial_ray_module
!       and load them into angular_ray_module.
!
!    Input arguments:
!  nx            : x-array extent
!  ij_ray_dim    : number of y-zones on a processor before swapping with y
!  ik_ray_dim    : number of z-zones on a processor before swapping with z
!  ny            : y-array extent
!  j_ray_dim     : the number of radial zones on a processor after swapping with y
!  nz            : z-array extent
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
!  transpose_x_y : transposes variables for the y-sweeps
!
!    Include files:
!  radial_ray_module, angular_ray_module, parallel_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, rho_c, t_c, ye_c, ei_c, &
& x_cf, dy_cf, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c, z_nuc_rep_c, &
& be_nuc_rep_c, nse_c, flat_x, agr_c, grav_y_c, e_nu_c
USE angular_ray_module, ONLY : rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, psi0_y, &
& psi1_y, xn_y, a_nuc_rep_y, z_nuc_rep_y, be_nuc_rep_y, dyphys_c, nse_y, &
& flat_x_y, agr_y, grav_y_cy, e_nu_y
USE parallel_module, ONLY : myid_y

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of z-zones on a processor before swapping with z
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: j_ray_dim     ! the number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: nz            ! z-array extent
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
INTEGER                          :: ji_ray        ! x (radial) index of a specific y (angular) ray

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                        = 1
le                        = nnc

!-----------------------------------------------------------------------
!  Transpose variables
!-----------------------------------------------------------------------

CALL transpose_x_y(imin, imax, nx, jmin, jmax, ny, ij_ray_dim, ik_ray_dim,  &
& j_ray_dim, nz, nez, nnu, ls, le, nnc, n_proc, n_proc_y, n_proc_z, rho_c,  &
& t_c, ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c,        &
& z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_x, agr_c, e_nu_c, grav_y_c, rho_y, &
& t_y, ye_y, ei_y, u_y, v_y, w_y, nse_y, psi0_y, psi1_y, xn_y, a_nuc_rep_y, &
& z_nuc_rep_y, be_nuc_rep_y, flat_x_y, agr_y, grav_y_cy, e_nu_y )

!-----------------------------------------------------------------------
!  Compute physical y-coordinate differences
!-----------------------------------------------------------------------

DO ji_ray = 1,j_ray_dim
  DO j = jmin,jmax
    dyphys_c(j,ji_ray,:)  = x_cf(j_ray_dim * myid_y + ji_ray) * dy_cf(j)
  END DO ! j = jmin,jmax
END DO ! ji_ray = 1,j_ray_dim

RETURN
END SUBROUTINE radhyd_to_transpose_x_y
