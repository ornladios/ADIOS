!-----------------------------------------------------------------------
!    Module:       brem_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE brem_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!  idrb(j,ij_ray,ik_ray), itrb(j,ij_ray,ik_ray), and iyrb(j,ij_ray,ik_ray)
!   are integers defining the location of log(rho), log(t), and ye for
!   radial zone j on the grid for n-neutrino bremsstrahlung pair
!   annihilation functions, i.e.,
!
!     idrb(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrb(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itrb(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrb(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrb(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrb(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to  as the unit cube j. Bremsstrahlung pair
!   annihilation rates for radial zone j are stored at the corners of unit
!   cube j. Rates for zone j are interpolated from the rates stored at the
!   corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrb
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrb
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrb


!-----------------------------------------------------------------------
!  Bremsstrahlung pair annihilation functions at cube corners
!-----------------------------------------------------------------------
!  brema0(j,k,kp,id,it,iy,f(ij_ray,ik_ray)) : arrays of the zero moments
!   of the  bremsstrahlung pair annihilation  functions for neutrinos of
!   energy k, antineutrinos of energy kp, in radial zone j at the unit
!   cube corners id, it, and iy (id, it, iy = 1,2). Interpolations of
!   the bremsstrahlung pair annihilation functions are computed from
!   this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: brema0


!-----------------------------------------------------------------------
!  Interpolated bremsstrahlung pair annihilation functions
!-----------------------------------------------------------------------
!  baf(j,k,n)        : Bremsstrahlung pair annihilation function i
!
!     i = 1: a0w, i = 2: b0w, i = 3: c0w,
!     i = 4: a1w, i = 5: b1w, i = 6: c1w.
!
!  bafd(i,j,k,n)     : d(baf(i,j,k,n))/d(density)
!  baft(i,j,k,n)     : d(baf(i,j,k,n))/d(temperature)
!  bafy(i,j,k,n)     : d(baf(i,j,k,n))/d(electron fraction)
!  bafp0(i,kp,j,k,n) : d(baf(i,j,k,n))/d(psi0(j,kp,n))
!  bafp1(i,kp,j,k,n) : d(baf(i,j,k,n))/d(psi1(j,kp,n))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: baf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: bafd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: baft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: bafy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: bafp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: bafp1

END module brem_module
