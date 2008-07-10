!-----------------------------------------------------------------------
!    Module:       pair_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE pair_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!  idrpp(j), itrpp(j), and iyrpp(j) are integers defining the location of
!   log(rho), log(t), and ye for radial zone j on the grid for the
!   electron-positron pair annihilation process, i.e.,
!
!     idrpp(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrpp(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itrpp(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrpp(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrpp(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrpp(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Electron-positron pair
!   annihilation rates for radial zone j are stored at the corners of
!   unit cube j. Rates for zone j are interpolated from the rates stored
!   at the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrpp
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrpp
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrpp

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation functions at cube corners
!-----------------------------------------------------------------------
!  paira0i(j,k,kp,id,it,iy,ij_ray,ik_ray),
!  paira0ii(j,k,kp,id,it,iy,ij_ray,ik_ray) : arrays of the zero
!   moments of the  electron-positron pair annihilation functions for
!   neutrinos of energy k, antineutrinos of energy kp, in radial zone
!   j at the unit cube corners id, it, and iy (id, it, iy = 1,2). 
!   Interpolations of the electron-positron pair annihilation functions
!   are computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: paira0i
REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: paira0ii

!-----------------------------------------------------------------------
!  Interpolated electron-positron pair annihilation functions
!-----------------------------------------------------------------------
!  paf(j,k,n)        : electron-positron pair annihilation function i
!
!     i = 1: a0w, i = 2: b0w, i = 3: c0w,
!     i = 4: a1w, i = 5: b1w, i = 6: c1w.
!
!  pafd(i,j,k,n)     : d(baf(i,j,k,n))/d(density)
!  paft(i,j,k,n)     : d(baf(i,j,k,n))/d(temperature)
!  pafy(i,j,k,n)     : d(baf(i,j,k,n))/d(electron fraction)
!  pafp0(i,kp,j,k,n) : d(baf(i,j,k,n))/d(psi0(j,kp,n))
!  pafp1(i,kp,j,k,n) : d(baf(i,j,k,n))/d(psi1(j,kp,n))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: paf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: pafd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: paft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: pafy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: pafp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: pafp1

END module pair_module
