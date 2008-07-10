!-----------------------------------------------------------------------
!    Module:       eos_snc_z_module
!    Author:       S. W. Bruenn
!    Date:         1/04/07
!-----------------------------------------------------------------------

MODULE eos_snc_z_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Nuclear statistical equilibrium
!-----------------------------------------------------------------------
!  nse(k,kj_ray,ki_ray) : a nuclear statistical equilibrium flag for
!   radial zone k.
!
!     nse(k,kj_ray,ki_ray) = 0 : material not in nuclear statistical
!      equilibrium; nuclear reaction network must be turned on to evolve
!      the matter composition.
!     nse(k,kj_ray,ki_ray) = 1 : material in nuclear statistical
!      equilibrium; nuclear reaction network turned off.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                  :: nse

!-----------------------------------------------------------------------
!  Equation of state identifier
!-----------------------------------------------------------------------
!  eos_i - equation of state identifier
!
!     eos_i = 'L'  : Lattimer-Swesty equation of state is used.
!     eos_i = 'B'  : Cooperstein-BCK equation of state is used.
!-----------------------------------------------------------------------

CHARACTER (LEN=1)                                       :: eos_i

!-----------------------------------------------------------------------
!  Equation of state borders
!-----------------------------------------------------------------------
!  rhopnu : density above which neutrino pressure is computed
!   isotropically as a equilibrated relativistic gas.
!     eos_i = 'L' : Lattimer-Swesty equation of state is used.
!     eos_i = 'B' : Cooperstein-BCK equation of state is used.
!
!  eosrho : the border density between the LS EOS and the BCK EOS (/fm3).
!-----------------------------------------------------------------------

REAL(KIND=double)                                       :: rhopnu
REAL(KIND=double)                                       :: eosrho

!-----------------------------------------------------------------------
!  Interpolated equation of state variables
!-----------------------------------------------------------------------
!  aesv(k,i,kj_ray,ki_ray)  : equation of state dependent variable i at
!   radial zone k.
!  aesvd(k,i,kj_ray,ki_ray) : derivative with respect to the density of
!   equation of state dependent variable i at radial zone k.
!  aesvt(k,i,kj_ray,ki_ray) : derivative with respect to the temperature
!   of equation of state dependent variable i at radial zone k.
!  aesvy(k,i,kj_ray,ki_ray) : derivative with respect to the electron
!   fraction of equation of state dependent variable i at radial zone k.
!
!     i = 1   : pressure
!     i = 2   : energy
!     i = 3   : entropy
!     i = 4   : neutron chemical potential
!     i = 5   : proton chemical potential
!     i = 6   : electron chemical potential
!     i = 7   : free neutron mass fraction
!     i = 8   : free proton mass fraction
!     i = 9   : heavy nucleus mass fraction
!     i = 10  : heavy nucleus mass number
!     i = 11  : heavy nucleus charge number
!     i = 12  : gamma1
!
!  gam1(k,kj_ray,ki_ray)    : first adiabatic index for zone k.
!  gam2(k,kj_ray,ki_ray)    : second adiabatic index for zone k.
!  gam3(k,kj_ray,ki_ray)    : third adiabatic index for zone k.
!
!  duesrc(k,kj_ray,ki_ray)  : cumulative energy 'glitches' per unit mass
!   in radial zone k due to passing from one unit cube to an adjacent
!   unit cube in the equation of state, and from rezoning.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: aesv
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: aesvd
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: aesvt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)      :: aesvy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)        :: gam1
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)        :: gam2
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)        :: gam3
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)        :: duesrc

!-----------------------------------------------------------------------
!  Cube grid
!-----------------------------------------------------------------------
!  dgrid(i), tgrid(i), ygrid(i) : log(rho), log(t), ye space is overlain
!   with a uniform grid of 'dgrid(i)' divisions per unit change in
!   log(rho), 'tgrid(i)' divisions per unit change in log(t), and
!   'ygrid(i)' divisions per 0.5 change in ye. Equation of state, nuclear
!   reaction rate, and neutrino interaction variables at each radial zone
!   are interpolated from values at the nearest corners on this grid.
!
!  rhoes(k) : The variables dgrid, tgrid, and ygrid are each 3 element
!   arrays permitting different partitionings of log(rho), log(t), ye
!   space in different density regimes delimited by rhoes.  These
!   different regimes are
!
!     regime 1 :             rho < rhoes(1)
!     regime 2 :        rhoes(1) < rho < rhoes(2)
!     regime 3 :             rhoes(2) < rho
!
!  idty(k,kj_ray,ki_ray)  : the density regime (i.e., 1, 2, or 3) of radial
!   zone k as given by the above inequalities.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                  :: idty

REAL(KIND=double), DIMENSION(3)                         :: dgrid
REAL(KIND=double), DIMENSION(3)                         :: tgrid
REAL(KIND=double), DIMENSION(3)                         :: ygrid
REAL(KIND=double), DIMENSION(2)                         :: rhoes

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!  idr(k,kj_ray,ki_ray), itr(k,kj_ray,ki_ray), and iyr(k,kj_ray,ki_ray) :
!   integers defining the location of log(rho), log(t), and ye for radial
!   zone k on the grid, or n-neutrino bremsstrahlung, i.e.,
!
!     idr(k,kj_ray,ki_ray)/dgrid < log(rho(k)) < ( idr(k,kj_ray,ki_ray) + 1 )/dgrid 
!
!     itrb(k,kj_ray,ki_ray)/tgrid <  log(t(k))  < ( itrb(k,kj_ray,ki_ray) + 1 )/tgrid
!
!     0.5 - iyrb(k,kj_ray,ki_ray)/ygrid < ye < 0.5 - ( iyrb(k,kj_ray,ki_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone k are referred to as the unit cube k. Equation of state quantities
!   for radial zone k are stored at the corners of unit cube k. Rates for
!   the equation of state quantities zone k are interpolated from the rates
!   stored at the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: idr
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: itr
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: iyr

!-----------------------------------------------------------------------
!  EOS table
!-----------------------------------------------------------------------
!  estble(i,k,id,it,iy,kj_ray,ki_ray) : equation of state variable i of
!   radial zone k y-zone ja, z-xone ka, at the unit cube corners id, it,
!   and iy (id, it, iy = 1,2). This is the table of equation of state
!   variables from which interpolations are performed.
!
!  escnst(i,k,kj_ray,ki_ray) : constant that is added to equation of state
!   variable i of radial zone k before taking logarithms for interpolation,
!   and then subtracted after interpolation in order to avoid taking the
!   log of a negative number.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: estble
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)       :: escnst

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------
!  nuc_number    : number of nuclear species (not counting representative
!   heavy nucleus).
!
!  xn(k,i)       : mass fraction of the ith nucleus.
!
!  a_nuc_rep(k)  : mass number of the representative heavy nucleus.
!
!  z_nuc_rep(k)  : charge number of the representative heavy nucleus.
!
!  be_nuc_rep(k) : binding energy of the representative heavy nucleus
!   (MeV).
!
!  a_nuc(n)      : mass number of the nth nuclear species.
!
!  z_nuc(n)      : charge number of the nth nuclear species.
!
!  be_nuc(n)     : binding energy of the nth nuclear species (MeV).
!-----------------------------------------------------------------------

INTEGER                                                :: nuc_number

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)         :: xn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)           :: be_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)           :: a_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)           :: z_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)           :: a_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)           :: z_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)           :: be_nuc

END module eos_snc_z_module
