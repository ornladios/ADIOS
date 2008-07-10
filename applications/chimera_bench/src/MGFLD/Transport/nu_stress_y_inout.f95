SUBROUTINE nu_stress_y_inout( jmin, jmax, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
& i_radial,  j_radial, nx, ny, nez, nnu, x_ef, rhobar, y_ei, rho_y, psi0_y, &
& nu_str_cy, nu_str_ey)
!-----------------------------------------------------------------------
!
!    File:         nu_stress_y_inout
!    Module:       nu_stress_y_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive the arrays from radial_ray_module and angular_ray_module,
!       compute the x-componsnt of neutrino stress, and return to
!       radial_ray_module the updated stresses.
!
!    Input arguments:
!  jmin         : lower y-array index
!  jmax         : upper y-array index
!  ji_ray       : x (radial) index of a specific angular ray
!  jk_ray       : z (azimuthal) index of a specific angular ray
!  j_ray_dim    : the number of radial zones on a processor after swapping with y
!  ik_ray_dim   : the number of z-zones on a processor before swapping with z
!  i_radial     : the unshifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!  j_radial     : the shifted radial zone (angular ray) corresponding to ji_ray, jk_ray
!  nx           : x-array extent
!  nx           : y-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  x_ef         : initial x grid zone faces
!  rhobar       : mean density at a given radius
!  y_ei         : initial y grid zone faces
!  rho_y        : density (g cm^{-3})
!  rhs1_c       : the right-hand side of transport equation, first moment
!  psi0_y       : zero moment of the neutrino occupation probability
!
!    Output arguments:
!  nu_str_cy    : y-component of zone-centered neutrino stress (dynes g^{-1})
!  nu_str_ey    : y-component of zone-edged neutrino stress (dynes g^{-1})
!
!    Subprograms called:
!  nu_stress_y  : computes the neutrino y-stress
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

INTEGER, INTENT(in)              :: jmin          ! minimum angular zone index
INTEGER, INTENT(in)              :: jmax          ! maximum angular zone index
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific angular ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific angular ray
INTEGER, INTENT(in)              :: j_ray_dim     ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim    ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ji_ray, jk_ray
INTEGER, INTENT(in)              :: j_radial      ! the shifted radial zone corresponding to ji_ray, jk_ray

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! y-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                 :: x_ef      ! x-coordinate of zone interface
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)                 :: y_ei      ! y-coordinate of zone interface
REAL(KIND=double), INTENT(in), DIMENSION(nx)                   :: rhobar    ! mean density at a given radius (g cm^{-3}
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: rho_y     ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim) :: psi0_y    ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)        :: nu_str_cy  ! y-component of zone-centered neutrino stress
REAL(KIND=double), INTENT(out), DIMENSION(ny+1,j_ray_dim,ik_ray_dim)      :: nu_str_ey  ! y-component of zone-edged neutrino stress

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(ny)         :: rho        ! shifted density (g cm^{-3})
REAL(KIND=double), DIMENSION(ny,nez,nnu) :: psi0       ! zero moment of the neutrino occupation probability
REAL(KIND=double), DIMENSION(ny)         :: nu_strs_cy ! y-component of zone-centered neutrino stress
REAL(KIND=double), DIMENSION(ny+1)       :: nu_strs_ey ! y-component of zone-edged neutrino stress

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Put radiation variables into y-arrays for passing into nu_stress_y
!-----------------------------------------------------------------------

rho (jmin:jmax)     = rho_y (jmin:jmax,ji_ray,jk_ray)

psi0(jmin:jmax,:,:) = psi0_y(jmin:jmax,:,:,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Compute neutrino stresses
!-----------------------------------------------------------------------

CALL nu_stress_y( jmin, jmax, ji_ray, jk_ray, j_ray_dim, ik_ray_dim, &
& i_radial, j_radial, nx, ny, nez, nnu, rho, x_ef, rhobar, y_ei, psi0, &
& nu_strs_cy, nu_strs_ey )
 
!-----------------------------------------------------------------------
!  Return neutrino stress
!-----------------------------------------------------------------------

nu_str_cy(:,ji_ray,jk_ray) = nu_strs_cy(:)
nu_str_ey(:,ji_ray,jk_ray) = nu_strs_ey(:)

RETURN
END SUBROUTINE nu_stress_y_inout
