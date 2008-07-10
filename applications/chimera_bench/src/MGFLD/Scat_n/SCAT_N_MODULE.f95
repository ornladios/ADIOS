!-----------------------------------------------------------------------
!    Module:       scat_n_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE scat_n_module

USE kind_module, ONLY : single, double

SAVE


!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!
!  idrn(j,ij_ray,ik_ray), itrn(j,ij_ray,ik_ray), and iyrn(j,ij_ray,ik_ray)
!  are integers defining the location of log(rho), log(t), and ye for
!  radial zone j on the grid, for the neutrino-nucleon elastic scattering,
!  i.e.,
!
!     idrn(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrn(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itnn(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrn(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrn(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrn(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Neutrino-nucleon elastic
!   scattering rates for radial zone j are stored at the corners of unit
!   cube j. Rates for zone j are interpolated from the rates stored at
!   the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrn
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrn
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrn

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering functions at cube corners
!-----------------------------------------------------------------------
!
!  sctn0(j,k,kp,id,it,iy,ij_ray,ik_ray),
!  sctn1(j,k,kp,id,it,iy,ij_ray,ik_ray) : arrays of the zero moments of
!   the neutrino-nucleon elastic scattering functions as a function of
!   the radial zone (j), incident neutrino energy (k), final neutrino
!   energy (kp) at the unit cube corners id, it, and iy (id, it, iy = 1,2). 
!   Interpolations of neutrino-nucleon elastic scattering functions are
!   computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: sctn0

!-----------------------------------------------------------------------
!  Antieutrino-nucleon elastic scattering functions at cube corners
!-----------------------------------------------------------------------
!
!  sctnb0(j,k,kp,id,it,iy,ij_ray,ik_ray),
!  sctnb1(j,k,kp,id,it,iy,ij_ray,ik_ray) : arrays of the zero moments of
!   the antineutrino-nucleon elastic scattering functions as a function
!   of the radial zone (j), incident antineutrino energy (k), final
!   antineutrino energy (kp) at the unit cube corners id, it, and iy
!   (id, it, iy = 1,2). 
!  Interpolations of antineutrino-nucleon elastic scattering functions
!   are computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: sctnb0


!-----------------------------------------------------------------------
!  Interpolated neutrino-nucleon inelastic scattering functions
!-----------------------------------------------------------------------
!
!  scnf(j,k,n) : Neutrino-nucleon elastic scattering function i
!
!                i = 1: a0w, i = 2: b0w, i = 3: c0w,
!                i = 4: a1w, i = 5: b1w, i = 6: c1w.
!
!      scnfd(i,j,k,n)     = d(scnf(i,j,k,n))/d(density)
!      scnft(i,j,k,n)     = d(scnf(i,j,k,n))/d(temperature)
!      scnfy(i,j,k,n)     = d(scnf(i,j,k,n))/d(electron fraction)
!      scnfp0(i,kp,j,k,n) = d(scnf(i,j,k,n))/d(psi0(j,kp,n))
!      scnfp1(i,kp,j,k,n) = d(scnf(i,j,k,n))/d(psi1(j,kp,n))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnfd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnfy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scnfp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scnfp1

END module scat_n_module
