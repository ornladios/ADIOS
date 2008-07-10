SUBROUTINE nu_stress_z_inout( kmin, kmax, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, i_radial, j_radial, nx, ny, nz, nez, nnu, x_ef, rhobar, &
& y_ei, z_ei, rho_z, psi0_z, nu_str_cz, nu_str_ez)
!-----------------------------------------------------------------------
!
!    File:         nu_stress_z_inout
!    Module:       nu_stress_z_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/22/07
!
!    Purpose:
!      To receive the arrays from radial_ray_module and angular_ray_module,
!       compute the x-componsnt of neutrino stress, and return to
!       radial_ray_module the updated stresses.
!
!    Input arguments:
!  kmin         : lower y-array index
!  kmax         : upper y-array index
!  ki_ray       : x (radial) index of a specific z (azimuthal) ray
!  kj_ray       : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim   : the number of y-zones on a processor before swapping with y
!  k_ray_dim    : the number of radial zones on a processor after swapping with z
!  i_radial     : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  j_radial     : the shifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nx           : x-array extent
!  ny           : y-array extent
!  nz           : z-array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  x_ef         : initial x grid zone faces
!  rhobar       : mean density at a given radius
!  y_ei         : initial y grid zone faces
!  z_ei         : initial z grid zone faces
!  rho_z        : density (g cm^{-3})
!  rhs1_c       : the right-hand side of transport equation, first moment
!  psi0_z       : zero moment of the neutrino occupation probability
!
!    Output arguments:
!  nu_str_cz    : z-component of zone-centered neutrino stress (dynes g^{-1})
!  nu_str_ez    : z-component of zone-edged neutrino stress (dynes g^{-1})
!
!    Subprograms called:
!  nu_stress_z  : computes the neutrino z-stress
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

INTEGER, INTENT(in)              :: kmin          ! minimum angular zone index
INTEGER, INTENT(in)              :: kmax          ! maximum angular zone index
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim    ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: k_ray_dim     ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial      ! the unshifted radial zone corresponding to ki_ray, kj_ray
INTEGER, INTENT(in)              :: j_radial      ! the shifted radial zone corresponding to ki_ray, kj_ray

INTEGER, INTENT(in)              :: nx            ! x-array extent
INTEGER, INTENT(in)              :: ny            ! x-array extent
INTEGER, INTENT(in)              :: nz            ! z-array extent
INTEGER, INTENT(in)              :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                 :: x_ef      ! x-coordinate of zone interface
REAL(KIND=double), INTENT(in), DIMENSION(ny+1)                 :: y_ei      ! y-coordinate of zone interface
REAL(KIND=double), INTENT(in), DIMENSION(nz+1)                 :: z_ei      ! z-coordinate of zone interface
REAL(KIND=double), INTENT(in), DIMENSION(nx)                   :: rhobar    ! mean density at a given radius (g cm^{-3}
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: rho_z     ! density (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz,nez,nnu,ij_ray_dim,k_ray_dim) :: psi0_z    ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)        :: nu_str_cz  ! z-component of zone-centered neutrino stress
REAL(KIND=double), INTENT(out), DIMENSION(nz+1,ij_ray_dim,k_ray_dim)      :: nu_str_ez  ! z-component of zone-edged neutrino stress

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nz)         :: rho        ! shifted density (g cm^{-3})
REAL(KIND=double), DIMENSION(nz,nez,nnu) :: psi0       ! zero moment of the neutrino occupation probability
REAL(KIND=double), DIMENSION(nz)         :: nu_strs_cz ! z-component of zone-centered neutrino stress
REAL(KIND=double), DIMENSION(nz+1)       :: nu_strs_ez ! z-component of zone-edged neutrino stress

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Put radiation variables into y-arrays for passing into nu_stress_z
!-----------------------------------------------------------------------

rho (kmin:kmax)     = rho_z (kmin:kmax,kj_ray,ki_ray)

psi0(kmin:kmax,:,:) = psi0_z(kmin:kmax,:,:,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Compute neutrino stresses
!-----------------------------------------------------------------------

CALL nu_stress_z( kmin, kmax, ki_ray, kj_ray, ij_ray_dim, k_ray_dim, &
& i_radial, j_radial, nx, ny, nz, nez, nnu, rho, x_ef, rhobar, y_ei, &
& z_ei, psi0, nu_strs_cz, nu_strs_ez )
 
!-----------------------------------------------------------------------
!  Return neutrino stress
!-----------------------------------------------------------------------

nu_str_cz(:,kj_ray,ki_ray) = nu_strs_cz(:)
nu_str_ez(:,kj_ray,ki_ray) = nu_strs_ez(:)

RETURN
END SUBROUTINE nu_stress_z_inout
