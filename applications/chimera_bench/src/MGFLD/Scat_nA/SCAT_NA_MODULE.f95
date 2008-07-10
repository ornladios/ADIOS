!-----------------------------------------------------------------------
!    Module:       scat_nA_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE scat_nA_module

USE kind_module

SAVE


!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!
!  idrnA(j,ij_ray,ik_ray), itrnA(j,ij_ray,ik_ray), and iyrnA(j,ij_ray,ik_ray)
!   are integers defining the location of log(rho), log(t), and ye for
!   radial zone j on the grid, for the neutrino-nucleus inelastic scattering,
!   i.e.,
!
!     idrnA(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrnA(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itrnA(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrnA(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrnA(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrnA(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Neutrino-nucleus inelastic
!   scattering rates for radial zone j are stored at the corners of unit
!   cube j. Rates for zone j are interpolated from the rates stored at
!   the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrnA
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrnA
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrnA


!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering functions at cube corners
!-----------------------------------------------------------------------
!
!  sctnA0(j,k,kp,id,it,iy,f(ij_ray,ik_ray)): array of the zero moments
!   of the neutrino-nucleus inelastic scattering functions as a function
!   of the radial zone (j), incident neutrino energy (k), final neutrino 
!   energy (kp) at the unit cube corners id, it, and iy (id, it, iy = 1,2).
!   Interpolations of neutrino-nucleus inelastic scattering functions
!   are computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: sctnA0


!-----------------------------------------------------------------------
!  Interpolated neutrino-nucleus inelastic scattering functions
!-----------------------------------------------------------------------
!
!  scnAf(j,k,n) : Neutrino-nucleus inelastic scattering function i
!
!                i = 1: a0w, i = 2: b0w, i = 3: c0w,
!                i = 4: a1w, i = 5: b1w, i = 6: c1w.
!
!      scnAfd(i,j,k,n)     = d(scnAf(i,j,k,n))/d(density)
!      scnAft(i,j,k,n)     = d(scnAf(i,j,k,n))/d(temperature)
!      scnAfy(i,j,k,n)     = d(scnAf(i,j,k,n))/d(electron fraction)
!      scnAfp0(i,kp,j,k,n) = d(scnAf(i,j,k,n))/d(psi0(j,kp,n))
!      scnAfp1(i,kp,j,k,n) = d(scnAf(i,j,k,n))/d(psi1(j,kp,n))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnAf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnAfd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnAft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scnAfy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scnAfp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scnAfp1

END module scat_nA_module
