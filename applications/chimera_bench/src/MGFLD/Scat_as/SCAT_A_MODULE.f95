!-----------------------------------------------------------------------
!    Module:       scat_a_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE scat_a_module

USE kind_module

SAVE


!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!
!  idrsa(j,ij_ray,ik_ray), itrsa(j,ij_ray,ik_ray), and iyrsa(j,ij_ray,ik_ray)
!   are integers defining the location of log(rho), log(t), and ye for
!   radial zone j on the grid, or n-neutrino nucleus inelastic scattering,
!   i.e.,
!
!           idrsa(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrsa(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!           itrsa(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrsa(j,ij_ray,ik_ray) + 1 )/tgrid
!
!          0.5 - iyrsa(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrsa(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Neutrino-nucleus inelastic
!   scattering rates for radial zone j are stored at the corners of unit
!   cube j. Rates for zone j are interpolated from the rates stored at the
!   corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrsa
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrsa
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrsa


!-----------------------------------------------------------------------
!  n-neutrino-nucleus non-isoenergetic scattering data
!-----------------------------------------------------------------------
!
!  rncnu0(k,kp): zero moment of the n-neutrino-nucleus inelastic scattering
!   functions (per iron-like nucleus) and per incident neutrino of energy
!   k and final neutrino energy kp as given by Haxton. n-neutrino-nucleus
!   inelastic scattering functions are computed from this array in
!   subroutine scatacal, and the results are stored in the array scta0
!   below.                                   c
!
!  rncnb0(k,kp): zero moment of the n-antineutrino-nucleus inelastic
!   scattering functions (per iron-like nucleus) and per incident neutrino
!   of energy k and final neutrino energy kp as given by Haxton.
!   n-antineutrino-nucleus inelastic scattering functions are computed
!   from this array in subroutine scatacal, and the results are stored
!   in the array scta0 below.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: rncnu0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: rncnb0


!-----------------------------------------------------------------------
!  n-neutrino-nucleus inelastic scattering functions at cube corners
!-----------------------------------------------------------------------
!
!  scta0(j,k,kp,id,it,iy,f(ij_ray,ik_ray)): arrays of the n-neutrino-nucleus
!   inelastic scattering functions as a function of the radial zone (j),
!   incident neutrino energy (k), final neutrino energy (kp) at the unit
!   cube corners id, it, and iy (id, it, iy = 1,2). Interpolations of
!   zero moments of the n-neutrino-nucleus and n-antineutrino-nucleus
!   inelastic scattering functions are computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: scta0_nue
REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: scta0_nueb
REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: scta0_nux
REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: scta0_nuxb


!-----------------------------------------------------------------------
!  Interpolated n-neutrino-nucleus inelastic scattering functions.
!-----------------------------------------------------------------------
!
!  scaf(j,k,n) : neutrino-nucleus inelastic scattering function i
!
!                i = 1: a0w, i = 2: b0w, i = 3: c0w,
!                i = 4: a1w, i = 5: b1w, i = 6: c1w.
!
!      scafd(i,j,k,n) = d(scaf(i,j,k,n))/d(density)
!      scaft(i,j,k,n) = d(scaf(i,j,k,n))/d(temperature)
!      scafy(i,j,k,n) = d(scaf(i,j,k,n))/d(electron fraction)
!      scafp0(i,kp,j,k,n) = d(scaf(i,j,k,n))/d(psi0(j,kp,n))
!      scafp1(i,kp,j,k,n) = d(scaf(i,j,k,n))/d(psi1(j,kp,n))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scaf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scafd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scaft
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)         :: scafy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scafp0
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:)       :: scafp1

END module scat_a_module
