!-----------------------------------------------------------------------
!    Module:       nucbrn_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE nucbrn_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Nuclear network controls
!-----------------------------------------------------------------------
!  inuc : nuclear reaction network switch.
!
!     inuc = 0: nuclear reactions (for nse(j) = 0 zones) bypassed, 
!      tnse set to 0.44 for flashing
!     inuc = 1: nuclear reactions (for nse(j) = 0 zones) included,
!      tnse set to higher value to assure smooth transition to NSE
!     inuc = 2: no nuclear network, no flashing
!-----------------------------------------------------------------------

INTEGER                                                  :: inuc

!-----------------------------------------------------------------------
!  General edit parameters
!-----------------------------------------------------------------------
!  nprint : a parameter used in the call to an edit subroutine denoting
!   the unit number of the file to which the print file is to be sent.
!-----------------------------------------------------------------------

INTEGER                                                  :: nprint

!-----------------------------------------------------------------------
!  Computational cycles
!-----------------------------------------------------------------------
!  ncycle : the cycle number of the calculation.
!-----------------------------------------------------------------------

INTEGER                                                  :: ncycle

!-----------------------------------------------------------------------
!  Nuclear statistical equilibrium
!-----------------------------------------------------------------------
!  nse(j) : a nuclear statistical equilibrium flag for radial zone j.
!
!     nse(j) = 0 : material not in nuclear statistical equilibrium;
!      nuclear reaction network must be turned on to evolve the matter
!      composition.
!     nse(j) = 1 : material in nuclear statistical equilibrium; nuclear
!      reaction network turned off.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:)                       :: nse

!-----------------------------------------------------------------------
!  Time step controls
!-----------------------------------------------------------------------
!  dTmax          : maximum relative temperature change arising from
!   nuclear burning.
!
!  jdTmax         : radiail zone limiting the time step due the
!   temperature change arising from nuclear burning.
!
!  dynmax         : maximum relative composition change arising from
!   nuclear burning.
!
!  jdynmax        : radiail zone limiting the time step due the
!   composiiton change arising from nuclear burning.
!
!  dynmax, ynmin  : parameters used in determining the composition change
!   time step due change arising to nuclear reactions.
!
!     dynmax = max( abs(dyn(i))/( yn(i) + ynmin ) )
!
!  jdynmax        : the mass zone for which dynmax is maximum.
!
!  t_cntl_burn(1) : nuclear burn temperature change time step criterion,
!   i.e., the maximum permitted abs( dT_burn(j)/t(j) ), where dT_burn(j)
!   is the nuclear burn temperature change in radial zone j.
!
!  t_cntl_burn(2) : nuclear burn composition change time step criterion,
!   i.e., the maximum permitted abs( dyn(j,i)/yn(j,i) ), where dyn(j,i)
!   is the abundance change of specie i in radial zone j.
!
!  dtime_burn(1)  : the maximum time step permitted by the nuclear burn
!   temperature change.
!  dtime_burn(2)  : the maximum time step permitted by the nuclear burn
!   composition change.
!-----------------------------------------------------------------------

INTEGER                                                  :: jdTmax
INTEGER                                                  :: jdynmax

REAL(KIND=double)                                        :: dTmax
REAL(KIND=double)                                        :: dynmax
REAL(KIND=double), DIMENSION(2)                          :: t_cntl_burn
REAL(KIND=double), DIMENSION(2)                          :: dtime_burn

!-----------------------------------------------------------------------
!  Temperature
!-----------------------------------------------------------------------
!  t(j) : temperature in radial zone j after nuclear burn
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: T_burn

!-----------------------------------------------------------------------
!  Temperature change
!-----------------------------------------------------------------------
!  dT_burn(j) : temperature change in radial zone j due to nuclear
!   reactions.
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: dT_burn

!-----------------------------------------------------------------------
!  Abundance parameters
!-----------------------------------------------------------------------
!  a_name            : Mass and symbol of nucleus.
!
!  nuc_number        : number of nuclear species (not counting
!   representative heavy nucleus).
!
!  xn(j,i)           : mass fraction of the ith nucleus.
!
!  dudt_nuc(j,ij_ray,ik_ray) : energy generation rate in zone j by
!   nuclear reactions for the current time step (ergs/gm).
!
!  uburn(j,ij_ray,ik_ray)    : cumulative energy generated in zone j
!   by nuclear reactions (ergs/gm).
!
!  a_nuc_rep(j)      : mass number of the representative heavy nucleus.
!
!  z_nuc_rep(j)      : charge number of the representative heavy
!   nucleus.
!
!  be_nuc_rep(j)     : binding energy of the representative heavy
!   nucleus (MeV).
!
!  a_nuc(n)          : mass number of the nth nuclear species.
!
!  z_nuc(n)          : charge number of the nth nuclear species.
!
!  m_ex_nuc(n)       : mass excess of the nth nuclear species (MeV).
!
!  be_nuc(n)         : binding energy of the nth nuclear species (MeV).
!-----------------------------------------------------------------------

CHARACTER (len=5), ALLOCATABLE, DIMENSION(:)             :: a_name

INTEGER                                                  :: nuc_number

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)           :: xn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)         :: dudt_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)         :: uburn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: be_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: a_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: z_nuc_rep
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: a_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: z_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: m_ex_nuc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)             :: be_nuc

!-----------------------------------------------------------------------
!  Nuclear reaction rate tolerances
!-----------------------------------------------------------------------
!  itnuc   : maximum number of iterations to attempt in order to obtain a
!   convergent solution of the reaction rate equations for the ion
!   abundances (i.e., the xn(j,i)'s) and the temperature in zones not
!   assumed to be in nse. If itnuc = 1, the variables are assumed to have
!   converged after the first iteration attempt.
!
!  ttolnuc : temperature convergence parameter for nuclear reaction rate
!   equations. The criterion for temperature convergence is that
!
!     abs(dt/t) < ttolnuc .
!
!  ytolnuc : ion abundance convergence parameter for nuclear reaction
!   rate equations. The criteria for ion abundance convergence is that
!
!     abs(dy)/(y(i) + ynmin) < ytolnuc .
!-----------------------------------------------------------------------

INTEGER                                                  :: itnuc

REAL(KIND=double)                                        :: ttolnuc
REAL(KIND=double)                                        :: ytolnuc
REAL(KIND=double)                                        :: ynmin

!-----------------------------------------------------------------------
!  Reaction rates
!-----------------------------------------------------------------------
!  rdpg  : 2H(p,g)3He               1
!
!  rhegp : 3He(g,p)2H               2
!
!  r3a   : 3 alpha reaction         3
!
!  rg3a  : 12C(g,a)2a               4
!
!  rcag  : 12C(a,g)16O              5
!
!  roga  : 16O(g,a)12C              6
!
!  roag  : 16O(a,g)Ne20             7
!
!  rnega : 20Ne(g,a)16O             8
!
!  rneag : 20Ne(a,g)24Mg            9
!
!  rmgga : 24Mg(g,a)20Ne           10
!
!  rmgag : 24Mg(a,g)28Si           11
!
!  rsiga : 28Si(g,a)24Mg           12
!
!  rcaag : 40Ca(a,g)44Ti           13
!
!  rtiga : 44Ti(g,a)40Ca           14
!
!  r1212 : 12C + 12C               15
!
!  r1216 : 12C + 16O               16
!
!  r1616 : 16O + 16O               17
!-----------------------------------------------------------------------

REAL(KIND=double)                                        :: rdpg
REAL(KIND=double)                                        :: rhegp
REAL(KIND=double)                                        :: r3a
REAL(KIND=double)                                        :: rg3a
REAL(KIND=double)                                        :: rcag
REAL(KIND=double)                                        :: roga
REAL(KIND=double)                                        :: roag
REAL(KIND=double)                                        :: rnega
REAL(KIND=double)                                        :: rneag
REAL(KIND=double)                                        :: rmgga
REAL(KIND=double)                                        :: rmgag
REAL(KIND=double)                                        :: rsiga
REAL(KIND=double)                                        :: rcaag
REAL(KIND=double)                                        :: rtiga
REAL(KIND=double)                                        :: r1212
REAL(KIND=double)                                        :: r1216
REAL(KIND=double)                                        :: r1616

!-----------------------------------------------------------------------
!  Screening corrections
!-----------------------------------------------------------------------
!  fescrn(j,i) : electron screening correction at the (2,2,2) cube corner
!   for the ith nuclei pair in radial zone j.
!
!     i =  1: p-p
!     i =  2: 4He-4He
!     i =  3: 4He-12C
!     i =  4: 4He-16O
!     i =  5: 4He-20Ne
!     i =  6: 4He-24Mg
!     i =  7: 4He-40Ca
!     i =  8: 12C-12C
!     i =  9: 12C-16O
!     i = 10: 16O-16O
!
!  fascrn(j,i) : ion screening correction at the (2,2,2) cube corner for
!   the ith nuclei pair in radial zone j.
!
!     i =  1: p-p
!     i =  2: 4He-4He
!     i =  3: 4He-12C
!     i =  4: 4He-16O
!     i =  5: 4He-20Ne
!     i =  6: 4He-24Mg
!     i =  7: 4He-40Ca
!     i =  8: 12C-12C
!     i =  9: 12C-16O
!     i = 10: 16O-16O
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)       :: fescrn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)       :: fascrn

!-----------------------------------------------------------------------
!  Cube indices
!-----------------------------------------------------------------------
!  idrnc(j,ij_ray,ik_ray), itrnc(j,ij_ray,ik_ray), and iyrnc(j,ij_ray,ik_ray) :
!   integers defining the location oflog(rho), log(t), and ye for radial
!   zone j on the grid, or reaction rates, i.e.,
!
!     idrnc(j,ij_ray,ik_ray)/dgrid < log(rho(j)) < ( idrncb(j,ij_ray,ik_ray) + 1 )/dgrid 
!
!     itrnc(j,ij_ray,ik_ray)/tgrid <  log(t(j))  < ( itrnc(j,ij_ray,ik_ray) + 1 )/tgrid
!
!     0.5 - iyrnc(j,ij_ray,ik_ray)/ygrid < ye < 0.5 - ( iyrnc(j,ij_ray,ik_ray) + 1 )/ygrid
!
!  The eight grid points surrounding log(rho), log(t), and ye for radial
!   zone j are referred to as the unit cube j. Reaction rates for radial
!   zone j are stored at the corners of unit cube j. Reaction ratws for
!   zone j are interpolated from the rates stored at the corners.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: idrnc
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: itrnc
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)                   :: iyrnc

!-----------------------------------------------------------------------
!  Cube grid
!-----------------------------------------------------------------------
!  dgrid(i), tgrid(i), ygrid(i): log(rho), log(t), ye space is overlain
!   with a uniform grid of 'dgrid(i)' divisions per unit change in
!   log(rho), 'tgrid(i)' divisions per unit change in log(t), and
!   'ygrid(i)' divisions per 0.5 change in ye. Equation of state, nuclear
!   reaction rate, and neutrino interaction variables at each radial zone
!   are interpolated from values at the nearest corners on this grid.
!
!  rhoes(k) : The variables dgrid, tgrid, and ygrid are each 3 element
!   arrays permitting different partitionings of log(rho), log(t), ye
!   space in different density regimes delimited by rhoes. These different
!   regimes are
!
!     regime 1:             rho < rhoes(1)
!     regime 2:        rhoes(1) < rho < rhoes(2)
!     regime 3:             rhoes(2) < rho
!
!  idty(j)  : the density regime (i.e., 1, 2, or 3) of radial zone j as
!   given by the above inequalities.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:)                       :: idty

REAL(KIND=double), DIMENSION(3)                          :: dgrid
REAL(KIND=double), DIMENSION(3)                          :: tgrid
REAL(KIND=double), DIMENSION(3)                          :: ygrid
REAL(KIND=double), DIMENSION(2)                          :: rhoes

!-----------------------------------------------------------------------
!  Reaction rate table
!-----------------------------------------------------------------------
!  rrdata(i,j,id,it,iy,ij_ray,ik_ray) : reaction rate i of radial zone
!   j at the unit cube corners id, it, and iy (id, it, iy = 1,2). This
!   is the table of reaction rates from which interpolations are performed.
!   The reaction rate numbers are as given above for rdpg, rhegp, etc.
!-----------------------------------------------------------------------

REAL(KIND=single), ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) ::  rrdata

END module nucbrn_module
