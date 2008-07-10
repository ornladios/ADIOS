!-----------------------------------------------------------------------
!    Module:       scat_nn_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE scat_nn_module

USE kind_module

SAVE


!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!
!  idrnn(j,ij_ray,ik_ray), itrnn(j,ij_ray,ik_ray), and iyrnn(j,ij_ray,ik_ray)
!   are integers defining the location of log(rho), log(t), and ye for 
!   radial zone j on the grid, for the neutrino-nucleon inelastic scattering,
!   i.e.,
!
!     idrnn(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrnn(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itrnn(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrnn(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrnn(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrnn(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Neutrino-nucleon inelastic
!   scattering rates for radial zone j are stored at the corners of unit
!   cube j. Rates for zone j are interpolated from the rates stored at
!   the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrnn
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrnn
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrnn


!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic scattering functions at cube corners
!-----------------------------------------------------------------------
!
!  sctnn0(j,k,kp,id,it,iy,ij_ray,ik_ray): array of the zero moments of
!   the neutrino-nucleon inelastic scattering functions as a function of
!   the radial zone (j), incident neutrino energy (k), final neutrino 
!   energy (kp) at the unit cube corners id, it, and iy (id, it, iy = 1,2).
!   Interpolations of neutrino-nucleon inelastic scattering functions
!   are computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: sctnn0


!-----------------------------------------------------------------------
!  nterpolated neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------
!
!  scnnf(j,k,n) : Neutrino-nucleon inelastic scattering function i
!
!                i = 1: a0w, i = 2: b0w, i = 3: c0w,
!                i = 4: a1w, i = 5: b1w, i = 6: c1w.
!
!      scnnfd(i,j,k,n)     = d(scnnf(i,j,k,n))/d(density)
!      scnnft(i,j,k,n)     = d(scnnf(i,j,k,n))/d(temperature)
!      scnnfy(i,j,k,n)     = d(scnnf(i,j,k,n))/d(electron fraction)
!      scnnfp0(i,kp,j,k,n) = d(scnnf(i,j,k,n))/d(psi0(j,kp,n))
!      scnnfp1(i,kp,j,k,n) = d(scnnf(i,j,k,n))/d(psi1(j,kp,n))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnnf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnnfd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnnft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnnfy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scnnfp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scnnfp1

END module scat_nn_module
