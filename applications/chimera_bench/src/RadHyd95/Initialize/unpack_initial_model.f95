SUBROUTINE unpack_initial_model( nx, nez, nnu, ij_ray_dim, ik_ray_dim, &
& i_model_data, d_model_data1, d_model_data2, d_model_data3, d_psi_data2 )
!-----------------------------------------------------------------------
!
!    File:         unpack_initial_model
!    Module:       unpack_initial_model
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To unpack the model configuration arrays and restore the values
!       to the appropriate variables in radial_ray_module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx            : radial array dimension
!  ij_ray_dim    : number of y-zones on a processor before swapping
!  ik_ray_dim    : number of z-zones on a processor before swapping
!  i_model_data  : integer array of nuclear keys
!  d_model_data1 : 64 bit real array of initial model data
!  d_model_data2 : 64 bit real array of initial model data
!  d_model_data3 : 64 bit real array of initial model data
!  d_psi_data2   : 64 bit real array of neutrino data
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  mdl_cnfg_module, radial_ray_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : half

USE eos_snc_x_module, ONLY : duesrc
USE mdl_cnfg_module, ONLY : jr_min, jr_max
USE radial_ray_module, ONLY : rho_c, rho_ci, t_c, t_ci, ye_c, ye_ci, u_c, &
& v_c, w_c, dx_ci, x_ei, x_ci, psi0_c, e_nu_c_bar, f_nu_e_bar

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: nx            ! energy array dimension
INTEGER, INTENT(in)               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)               :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)               :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)               :: ik_ray_dim    ! number of z-zones on a processor before swapping

INTEGER, INTENT(in), DIMENSION(2) :: i_model_data  ! integer array of edit keys

REAL(KIND=double), INTENT(in), DIMENSION(7,nx,ij_ray_dim,ik_ray_dim)         :: d_model_data1 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(2,nx)                               :: d_model_data2 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(2,nx+1)                             :: d_model_data3 ! 64 bit real array of initial model data
REAL(KIND=double), INTENT(in), DIMENSION(2,nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: d_psi_data2   ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                      \\\\\ UNPACK DATA /////
!
!-----------------------------------------------------------------------

jr_min              = i_model_data(1) + 1
jr_max              = i_model_data(2) + 1

rho_c (:,:,:)       = d_model_data1(1,:,:,:)
rho_ci(:,:,:)       = d_model_data1(1,:,:,:)
t_c   (:,:,:)       = d_model_data1(2,:,:,:)
t_ci  (:,:,:)       = d_model_data1(2,:,:,:)
ye_c  (:,:,:)       = d_model_data1(3,:,:,:)
ye_ci (:,:,:)       = d_model_data1(3,:,:,:)
u_c   (:,:,:)       = d_model_data1(4,:,:,:)
v_c   (:,:,:)       = d_model_data1(5,:,:,:)
w_c   (:,:,:)       = d_model_data1(6,:,:,:)
duesrc(:,:,:)       = d_model_data1(7,:,:,:)

dx_ci     (:)       = d_model_data2(1,:)
e_nu_c_bar(:)       = d_model_data2(2,:)

x_ei      (:)       = d_model_data3(1,:)
f_nu_e_bar(:)       = d_model_data3(2,:)

x_ci(1:nx)          = half * ( x_ei(1:nx) + x_ei(2:nx+1) )

psi0_c(:,:,:,:,:)   = d_psi_data2(1,:,:,:,:,:)

RETURN
END SUBROUTINE unpack_initial_model
