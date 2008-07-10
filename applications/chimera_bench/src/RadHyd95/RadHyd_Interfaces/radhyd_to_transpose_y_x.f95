SUBROUTINE radhyd_to_transpose_y_x( nx, ij_ray_dim, ik_ray_dim, ny, j_ray_dim, &
& nez, nnu, nnc, n_proc, n_proc_y )
!-----------------------------------------------------------------------
!
!    File:         radhyd_to_transpose_y_x
!    Module:       radhyd_to_transpose_y_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To transpose certain of the y-array variables in angular_ray_module
!       and load them into radial_ray_module.
!
!    Input arguments:
!  nx         : x-array extent
!  ij_ray_dim : number of x (radial) zones on a processor before swapping with y
!  ik_ray_dim : number of z (azimuthal) zones on a processor before swapping with z
!  ny         : y-array extent
!  y_ray_dim  : number of angular rays assigned to a processor
!  nez        : neutrino energy array dimension
!  nnu        : neutrino flavor array dimension
!  nnc        : composition array dimension
!  n_proc     : number of processors assigned to the run
!  n_proc_y   : number of processors assigned to the y-indices of rays
!
!    Output arguments:
!        none
!
!    Subprograms called:
!  transpose_y_x : does the actual variable transpose
!
!    Include files:
!  radial_ray_module, angular_ray_module
!
!-----------------------------------------------------------------------

USE radial_ray_module, ONLY : imin, imax, jmin, jmax, rho_c, t_c, ye_c, ei_c, &
& x_cf, y_ef, dy_cf, y_cf, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c, &
& z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_y_x
USE angular_ray_module, ONLY : rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, psi0_y, &
& psi1_y, xn_y, a_nuc_rep_y, z_nuc_rep_y, be_nuc_rep_y, nse_y, flat_y
USE shock_module, ONLY : pq_y, pqy_x

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: j_ray_dim     ! number of angular rays assigned to a processor
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent
INTEGER                          :: nnc           ! composition array extent
INTEGER, INTENT(in)              :: n_proc        ! number of processors assigned to the run
INTEGER, INTENT(in)              :: n_proc_y      ! number of processors assigned to the y-indices of rays

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: ls            ! inner logical composition zone
INTEGER                          :: le            ! outer logical composition zone

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ls                    = 1
le                    = nnc

!-----------------------------------------------------------------------
!  Transpose variables.
!-----------------------------------------------------------------------

CALL transpose_y_x(imin, imax, nx, ij_ray_dim, ik_ray_dim, jmin,       &
& jmax, ny, j_ray_dim, nez, nnu, ls, le, nnc, n_proc, n_proc_y, rho_c, &
& t_c, ye_c, ei_c, u_c, v_c, w_c, psi0_c, psi1_e, xn_c, a_nuc_rep_c,   &
& z_nuc_rep_c, be_nuc_rep_c, nse_c, flat_y_x, pq_y, rho_y, t_y, ye_y,  &
& ei_y, u_y, v_y, w_y, psi0_y, psi1_y, xn_y, a_nuc_rep_y, z_nuc_rep_y, &
& be_nuc_rep_y, nse_y, flat_y, pqy_x )

RETURN
END SUBROUTINE radhyd_to_transpose_y_x
