!-----------------------------------------------------------------------
!    Module:       abem_z_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE abem_z_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!  drae(j), itrae(j), and iyrae(j) are integers defining the location of
!   log(rho), log(t), and ye for radial zone j on the grid, or n-neutrino
!   absorption and emission, i.e.,
!
!     idrae(j,kj_ray,ki_ray)/dgrid < log(rho(j)) < ( idrae(j,kj_ray,ki_ray) + 1 )/dgrid 
!
!     itrae(j,kj_ray,ki_ray)/tgrid <  log(t(j))  < ( itrae(j,kj_ray,ki_ray) + 1 )/tgrid
!
!     0.5 - iyrae(j,kj_ray,ki_ray)/ygrid < ye < 0.5 - ( iyrae(j,kj_ray,ki_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Absorption and emission
!   rates for radial  zone j are stored at the corners of unit cube j.
!   Rates for zone j are interpolated from the rates stored at the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: idrae
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: itrae
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                     :: iyrae

!-----------------------------------------------------------------------
!  Haxton's emission and absporption rates
!-----------------------------------------------------------------------
!  rncem(k,kp,kj_ray,ki_ray) stores the e-neutrino-nucleus absorption rate
!   (per iron-like nucleus) and per incident neutrino of energy k
!   and final electrons of energy kp as given by Haxton.
!   Absorption and emission rates are computed from this array in
!   subroutine abemhnc, and the results contribute to the arrays
!   ab and em below.
!
!  rncep(k,kp,kj_ray,ki_ray) stores the e-antineutrino-nucleus absorption rate
!   (per iron-like nucleus) and per incident neutrino of energy k
!   and final electrons of energy kp as given by Haxton.
!   Absorption and emission rates are computed from this array in
!   subroutine abemhnc, and the results contribute to the arrays
!   ab and em below.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: rncem
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)             :: rncep

!-----------------------------------------------------------------------
!  Totsl emission and absporption inverse mean free paths at cube corners
!-----------------------------------------------------------------------
!  em(j,k,n,id,it,iy,kj_ray,ki_ray) is an array of emission rates for
!   neutrinos as a function of neutrino type (n), energy (k) and angular
!   zone (j) at the unit cube corners id, it, and iy (id, it, iy = 1,2).
!   Interpolations of emission rates are computed from this table.
!
!  ab(j,k,n,id,it,iy,kj_ray,ki_ray) is an array of absorption rates for
!   neutrinos as a function of neutrino type (n), energy (k) and angular
!   zone (j) at the unit cube corners id, it, and iy (id, it, iy = 1,2).
!   Interpolaitons of absorption rates are computed from this table.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: em
REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:)   :: ab

!-----------------------------------------------------------------------
!  Emission inverse mean free paths
!-----------------------------------------------------------------------
!  emis(j,k,n)  : is the emission inverse mean free path (1/cm) for
!   angular zone j, energy zone k, neutrino type n
!
!  emisd(j,k,n) : d(emis(j,k,n))/d(rho(j))
!  emist(j,k,n) : d(emis(j,k,n))/d(t  (j))
!  emisy(j,k,n) : d(emis(j,k,n))/d(ye (j))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: emis
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: emisd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: emist
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: emisy

!-----------------------------------------------------------------------
!  Absorption inverse mean free paths
!-----------------------------------------------------------------------
!
!  absor(j,k,n)  : is the absorption inverse mean free path (1/cm)
!   for angular zone j, energy zone k, neutrino type n
!
!  absord(j,k,n) : d(absor(j,k,n))/d(rho(j))
!  absort(j,k,n) : d(absor(j,k,n))/d(t  (j))
!  absory(j,k,n) : d(absor(j,k,n))/d(ye (j))
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: absor
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: absord
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: absort
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)           :: absory

END module abem_z_module
