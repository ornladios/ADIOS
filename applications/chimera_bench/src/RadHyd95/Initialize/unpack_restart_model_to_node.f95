SUBROUTINE unpack_restart_model_to_node( nx, nez, nnu, nnc, ij_ray_dim,   &
&  ik_ray_dim, i_model_data, d_model_data1, d_model_data2, d_model_data3, &
&  d_psi_data1, d_psi_data2, d_psi_data3, d_psi_data4, d_psi_data5,       &
&  i_nuc_data, d_nuc_data )
!-----------------------------------------------------------------------
!
!    File:         unpack_restart_model_to_node
!    Module:       unpack_restart_model_to_node
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/01/08
!
!    Purpose:
!      To unpack the complete model configuration arrays and restore
!       the values to the appropriate variables in radial_ray_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx            : radial array dimension
!  nez           : neutrino energy array extent
!  nnu           : neutrino flavor array extent
!  nnc           : composition array extent
!  ij_ray_dim    : number of y-zones on a processor before swapping
!  ik_ray_dim    : number of z-zones on a processor before swapping
!  i_model_data  : integer array of nuclear keys
!  d_model_data1 : 64 bit real array of initial model data
!  d_model_data2 : 64 bit real array of initial model data
!  d_model_data3 : 64 bit real array of initial model data
!  d_psi_data1   : 64 bit real array of neutrino data
!  d_psi_data2   : 64 bit real array of neutrino data
!  d_psi_data3   : 64 bit real array of neutrino data
!  d_psi_data4   : 64 bit real array of neutrino data
!  d_psi_data5   : 64 bit real array of neutrino data
!  i_nuc_data    : integer array of edit keys
!  d_nuc_data    : 64 bit real array of nuclear keys
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module, nucbrn_module,
!  nu_dist_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : n_proc
USE numerical_module, ONLY : zero, half

USE edit_module, ONLY : nlog, nu_r, nu_rt, nu_rho, nu_rhot, psi0dat, psi1dat
USE eos_snc_x_module, ONLY : nse_e=>nse, duesrc
USE mdl_cnfg_module, ONLY : jr_min, jr_max
USE nucbrn_module, ONLY : nse_n=>nse
USE nu_dist_module, ONLY : dnurad, unukrad, unujrad, nnukrad, nnujrad, &
& e_rad, unurad, elec_rad, nnurad
USE radial_ray_module, ONLY : rho_c, rho_ci, t_c, t_ci, ye_c, ye_ci, u_c, &
& v_c, w_c, dx_ci, x_ei, x_ci, psi0_c, xn_c, a_nuc_rep_c, z_nuc_rep_c, &
& be_nuc_rep_c, uburn_c, nse_c, e_nu_c_bar, f_nu_e_bar

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nx            ! energy array dimension
INTEGER, INTENT(in)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(inout)            :: nnc            ! composition array extent
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping

INTEGER, INTENT(in), DIMENSION(2)                           :: i_model_data  ! integer array of edit keys
INTEGER, INTENT(in), DIMENSION(10+nx,ij_ray_dim,ik_ray_dim) :: i_nuc_data    ! integer array of edit keys

REAL(KIND=double), INTENT(in), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)                :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(2,nx)                                      :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(2,nx+1)                                    :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(8,nez,nnu,ij_ray_dim,ik_ray_dim)           :: d_psi_data1   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim)        :: d_psi_data2   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,nx,nnu,ij_ray_dim,ik_ray_dim)            :: d_psi_data3   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,nnu,ij_ray_dim,ik_ray_dim)               :: d_psi_data4   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(2,ij_ray_dim,ik_ray_dim)                   :: d_psi_data5   ! 64 bit real array of edit keys
REAL(KIND=double), INTENT(in), DIMENSION(10+4*nnc+(nnc+4)*nx,ij_ray_dim,ik_ray_dim) :: d_nuc_data    ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                          :: i             ! do index
INTEGER                          :: n             ! composition index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                      \\\\\ UNPACK DATA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Model configuration data
!-----------------------------------------------------------------------

jr_min                 = i_model_data(1) + 1
jr_max                 = i_model_data(2) + 1

rho_c (:,:,:)          = d_model_data1(1,:,:,:)
rho_ci(:,:,:)          = d_model_data1(1,:,:,:)
t_c   (:,:,:)          = d_model_data1(2,:,:,:)
t_ci  (:,:,:)          = d_model_data1(2,:,:,:)
ye_c  (:,:,:)          = d_model_data1(3,:,:,:)
ye_ci (:,:,:)          = d_model_data1(3,:,:,:)
u_c   (:,:,:)          = d_model_data1(4,:,:,:)
v_c   (:,:,:)          = d_model_data1(5,:,:,:)
w_c   (:,:,:)          = d_model_data1(6,:,:,:)
duesrc(:,:,:)          = d_model_data1(7,:,:,:)

dx_ci     (:)          = d_model_data2(1,:)
x_ei      (:)          = d_model_data3(1,:)
x_ci(1:nx)             = half * ( x_ei(1:nx) + x_ei(2:nx+1) )

!-----------------------------------------------------------------------
!  Neutrino data
!-----------------------------------------------------------------------

psi0_c(:,:,:,:,:)      = d_psi_data2(1,:,:,:,:,:)
dnurad(:,:,:,:,:)      = d_psi_data2(2,:,:,:,:,:)

nu_r   (:,:,:,:)       = d_psi_data1(1,:,:,:,:)
nu_rt  (:,:,:,:)       = d_psi_data1(2,:,:,:,:)
nu_rho (:,:,:,:)       = d_psi_data1(3,:,:,:,:)
nu_rhot(:,:,:,:)       = d_psi_data1(4,:,:,:,:)
psi0dat(:,:,:,:)       = d_psi_data1(5,:,:,:,:)
psi1dat(:,:,:,:)       = d_psi_data1(6,:,:,:,:)
unukrad(:,:,:,:)       = d_psi_data1(7,:,:,:,:)
nnukrad(:,:,:,:)       = d_psi_data1(8,:,:,:,:)

unujrad(:,:,:,:)       = d_psi_data3(1,:,:,:,:)
nnujrad(:,:,:,:)       = d_psi_data3(2,:,:,:,:)

unurad(:,:,:)          = d_psi_data4(1,:,:,:)
nnurad(:,:,:)          = d_psi_data4(2,:,:,:)

e_rad   (:,:)          = d_psi_data5(1,:,:)
elec_rad(:,:)          = d_psi_data5(2,:,:)

e_nu_c_bar(:)          = d_model_data2(2,:)
f_nu_e_bar(:)          = d_model_data3(2,:)

!-----------------------------------------------------------------------
!  Abundance data
!-----------------------------------------------------------------------

nse_c(1:nx,:,:)        = i_nuc_data(11:nx+10,:,:)
nse_n(2:nx)            = i_nuc_data(11:nx+ 9,1,1)
nse_e(2:nx,:,:)        = i_nuc_data(11:nx+ 9,:,:)

DO n = 1,nnc
  DO i = 1,nx
    xn_c(i,n,:,:)      = d_nuc_data(10+4*nnc+(n-1)*nx+i,:,:)
  END DO ! i = 1,nx
END DO ! 1,nnc

DO i = 1,nx
  a_nuc_rep_c(i,:,:)   = d_nuc_data(10+4*nnc+(nnc+0)*nx+i,:,:)
  z_nuc_rep_c(i,:,:)   = d_nuc_data(10+4*nnc+(nnc+1)*nx+i,:,:)
  be_nuc_rep_c(i,:,:)  = d_nuc_data(10+4*nnc+(nnc+2)*nx+i,:,:)
  uburn_c(i,:,:)       = d_nuc_data(10+4*nnc+(nnc+3)*nx+i,:,:)
END DO ! i

RETURN
END SUBROUTINE unpack_restart_model_to_node
