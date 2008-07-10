SUBROUTINE nu_stress_x_inout( imin, imax, ij_ray, ik_ray, ij_ray_dim,  &
& ik_ray_dim, nx, nez, nnu, x_ei, rho_c, rhobar, rhs1_c, dc_e, psi0_c, &
& nu_str_c, nu_str_e )
!-----------------------------------------------------------------------
!
!    File:         nu_stress_x_inout
!    Module:       nu_stress_x_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive the arrays from radial_ray_module, compute the x-componsnt
!       of neutrino stress, and return to radial_ray_module the updated
!       stresses.
!
!    Input arguments:
!  imin         : lower unshifted x-array index
!  imax         : upper unshifted x-array index
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  nx           : x-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  x_ei         : initial x grid zone faces
!  rho_c        : density (g cm^{-3})
!  rhobar       : angularly averaged density [g cm^{-3}]
!  rhs1_c       : the right-hand side of transport equation, first moment
!  dc_e         : the diffusion coefficient for neutrinos
!  psi0_c       : zero moment of the neutrino occupation probability
!
!    Output arguments:
!  nu_str_c     : zone-centered x-component of neutrino stress (dynes g^{-1})
!  nu_str_e     : zone-edge x-component of neutrino stress (dynes g^{-1})
!
!    Subprograms called:
!  nu_stress_x
!      
!
!    Include files:
!  kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                            :: imin       ! lower unshifted x-array index
INTEGER, INTENT(in)                                            :: imax       ! upper unshifted x-array index
INTEGER, INTENT(in)                                            :: ij_ray     ! j-index of a radial ray
INTEGER, INTENT(in)                                            :: ik_ray     ! k-index of a radial ray
INTEGER, INTENT(in)                                            :: ij_ray_dim ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                                            :: ik_ray_dim ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)                                            :: nx         ! x-array extent
INTEGER, INTENT(in)                                            :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)                                            :: nnu        ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                 :: x_ei       ! x-coordinate zone edge
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rho_c      ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                   :: rhobar    ! angularly averaged density [g cm^{-3}]
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: rhs1_c     ! the right-hand side of transport equation, first moment
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: dc_e       ! the diffusion coefficient for neutrinos
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0_c     ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)        :: nu_str_c   ! zone-centered x-component of neutrino stress
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)        :: nu_str_e   ! zone-edge x-component of neutrino stress

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                  :: jr_min    ! minimum shifted radial zone index
INTEGER                                  :: jr_max    ! maximum shifted radial zone index

REAL(KIND=double), DIMENSION(nx)         :: r         ! radius (cm)
REAL(KIND=double), DIMENSION(nx)         :: rho       ! shifted density (g cm^{-3})
REAL(KIND=double), DIMENSION(nx,nez,nnu) :: rhs1      ! the right-hand side of transport equation, first moment
REAL(KIND=double), DIMENSION(nx,nez,nnu) :: dc        ! the diffusion coefficient for neutrinos
REAL(KIND=double), DIMENSION(nx,nez,nnu) :: psi0      ! zero moment of the neutrino occupation probability
REAL(KIND=double), DIMENSION(nx)         :: nu_strs_c ! x-component of neutrino stress
REAL(KIND=double), DIMENSION(nx)         :: nu_strs_e ! x-component of neutrino stress

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Add shift to minimum and maximum x-array extent
!-----------------------------------------------------------------------

jr_min                      = imin + 1
jr_max                      = imax + 1

!-----------------------------------------------------------------------
!  Put radiation variables into shifted arrays for passing into
!   nu_stress_x
!-----------------------------------------------------------------------

r   (imin:imax+1)           = x_ei  (imin:imax+1)

rho (jr_min:jr_max)         = rho_c (imin:imax,ij_ray,ik_ray)

psi0(jr_min:jr_max,:,:)     = psi0_c(imin:imax,  :,:,ij_ray,ik_ray)
rhs1(jr_min:jr_max,:,:)     = rhs1_c(imin:imax,  :,:,ij_ray,ik_ray)

dc  (imin:imax+1,  :,:)     = dc_e  (imin:imax+1,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Compute neutrino stresses
!-----------------------------------------------------------------------

CALL nu_stress_x( jr_min, jr_max, ij_ray, ik_ray, nx, nez, nnu, rho, rhobar, &
& r, rhs1, dc, psi0, nu_strs_c, nu_strs_e )
 
!-----------------------------------------------------------------------
!  Return neutrino stress (shifter nu_strs_c --> unshifted nu_str_c)
!-----------------------------------------------------------------------

nu_str_c(imin:imax,ij_ray,ik_ray)   = nu_strs_c(jr_min:jr_max)
nu_str_e(imin:imax+1,ij_ray,ik_ray) = nu_strs_e(imin:imax+1)

RETURN
END SUBROUTINE nu_stress_x_inout
