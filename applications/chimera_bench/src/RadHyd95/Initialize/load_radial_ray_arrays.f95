SUBROUTINE load_radial_ray_arrays( nx, nez, nnu, nnc, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         load_radial_ray_arrays
!    Module:       load_radial_ray_arrays
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/10/04
!
!    Purpose:
!      To load variables computed in setup into radial_ray_module module.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nx        : x-array extent
!  nez       : neutrino energy array extent
!  nnu       : neutrino flavor array extent
!  nnc       : number of species
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module
!  radial_ray_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module, nu_dist_module
!
!-----------------------------------------------------------------------

USE numerical_module, ONLY : zero, half, one, epsilon

USE radial_ray_module, ONLY : rho_c, t_c, ye_c, ei_c, p_c, gc_c, ge_c,  &
& u_c, v_c, w_c, xn_c, nse_c, uburn_c, be_nuc_rep_c, a_nuc_rep_c, z_nuc_rep_c, &
& idty, dgrid, tgrid, ygrid, rhoes, u_bcr, v_bcr, w_bcr, r_bcr, p_bcr, ei_bcr, &
& ye_bcr, temp_bcr, gc_bcr, ge_bcr, psi0_bcr, comp_bcr, x_ei, y_ei, z_ei, dx_ci, &
& dy_ci, dz_ci, x_ci, y_ci, z_ci, nprint, imax, rho_ci, t_ci, ye_ci, ei_ci, &
& u_ci, v_ci, w_ci, nu_str_e, psi1_e, unu_c, unub_c, dunu_c, unue_e, unube_e, &
& dunue_e, dc_e, rhs1_c

USE edit_module, ONLY : nprintp=>nprint, nlog
USE eos_snc_x_module, ONLY : aesv, idtyp=>idty, dgridp=>dgrid, tgridp=>tgrid, &
& ygridp=>ygrid, rhoesp=>rhoes
USE mdl_cnfg_module, ONLY : rho, t, ye, u, r
USE nu_dist_module, ONLY : stress_x, psi1,  unu, unub, dunu, unue, unube, &
& dunue, dc, rhs1

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)    :: nx            ! x-array extent
INTEGER, INTENT(in)    :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)    :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)    :: nnc           ! composition array dimension
INTEGER, INTENT(in)    :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)    :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                :: i             ! radial array index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  General edit parameter
!-----------------------------------------------------------------------

nprint                           = nprintp

!-----------------------------------------------------------------------
!
!                    \\\\\ STATE VARIABLES /////
!
!-----------------------------------------------------------------------

ei_c(1:nx-1,ij_ray,ik_ray)       = aesv(2:nx,2,ij_ray,ik_ray)
p_c (1:nx-1,ij_ray,ik_ray)       = aesv(2:nx,1,ij_ray,ik_ray)
gc_c(1:nx-1,ij_ray,ik_ray)       = aesv(2:nx,12,ij_ray,ik_ray)
u_c (1:nx-1,ij_ray,ik_ray)       = u   (2:nx)
ge_c(1:nx-1,ij_ray,ik_ray)       = one + aesv(2:nx,1,ij_ray,ik_ray)/( aesv(2:nx,2,ij_ray,ik_ray) * rho(2:nx) + epsilon )

nu_str_e(:,ij_ray,ik_ray)        = zero

DO i = 1,nx
  nu_str_e(i,ij_ray,ik_ray)      = SUM( stress_x(i,:,ij_ray,ik_ray) )
END DO ! i = 1,nx

nse_c(imax+1:nx,ij_ray,ik_ray)   = 0

ei_ci(1:nx-1,ij_ray,ik_ray)      = aesv(2:nx,2,ij_ray,ik_ray)
u_ci (1:nx-1,ij_ray,ik_ray)      = u_c (1:nx-1,ij_ray,ik_ray)
v_ci (1:nx-1,ij_ray,ik_ray)      = v_c (1:nx-1,ij_ray,ik_ray)
w_ci (1:nx-1,ij_ray,ik_ray)      = w_c (1:nx-1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!
!                  \\\\\ RADIATION VARIABLES /////
!
!-----------------------------------------------------------------------

psi1_e(:,:,:,ij_ray,ik_ray)      = psi1(:,:,:)
dc_e  (:,:,:,ij_ray,ik_ray)      = dc  (:,:,:)
rhs1_c(1:nx-1,:,:,ij_ray,ik_ray) = rhs1(2:nx,:,:)

unu_c  (1:nx-1,:,ij_ray,ik_ray)  = unu  (2:nx,:)
unub_c (1:nx-1,:,ij_ray,ik_ray)  = unub (2:nx,:)
dunu_c (1:nx-1,:,ij_ray,ik_ray)  = dunu (2:nx,:)
unue_e (:,:,ij_ray,ik_ray)       = unue (:,:)
unube_e(:,:,ij_ray,ik_ray)       = unube(:,:)
dunue_e(:,:,ij_ray,ik_ray)       = dunue(:,:)

!-----------------------------------------------------------------------
!
!                  \\\\\ BOUNDARY CONDITIONS /////
!
!-----------------------------------------------------------------------

r_bcr                            = rho (imax+1)
ei_bcr                           = aesv(imax+1,2,1,1)
ye_bcr                           = ye  (imax+1)
temp_bcr                         = t   (imax+1)
u_bcr                            = u   (imax+1)
v_bcr                            = zero
w_bcr                            = zero

comp_bcr(:)                      = xn_c(imax,:,1,1)

!-----------------------------------------------------------------------
!
!                \\\\\ EOS CUBE REGIME INDICES /////
!
!-----------------------------------------------------------------------

idty(1:nx-1,ij_ray,ik_ray)       = idtyp(2:nx,ij_ray,ik_ray)

dgrid(:)                         = dgridp(:)
tgrid(:)                         = tgridp(:)
ygrid(:)                         = ygridp(:)
rhoes(:)                         = rhoesp(:)

RETURN
END SUBROUTINE load_radial_ray_arrays
