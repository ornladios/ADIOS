SUBROUTINE read_pack_transport_keys( nreadp, nprint, iskip, nez, nezp1, nnu, &
& i_trans_data, d_trans_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_transport_keys
!    Module:       read_pack_transport_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To read in the transport keys defining the transport sources, solution
!       tolerances, and time step criteria, and to pack them into an integer
!       and a real*8 array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp       : unit number from which to read
!  nprint       : unit number from which to print
!  iskip        : echo data read flag
!  nez          : energy array dimension
!  nez+1        : nez + 1
!  nnu          : neutrino flavor dimension
!  nrst         : cycle number at start or restart
!
!    Output arguments:
!  i_trans_data : integer array of transport keys
!  d_trans_data : real*8 array of transport keys
!
!    Include files:
!  kind_module, numerical_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nez           ! energy array dimension
INTEGER, INTENT(in)              :: nezp1         ! nez + 1
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor dimension
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(40+2*nnu)                      :: i_trans_data  ! integer array of transport keys

REAL(KIND=double), INTENT(out), DIMENSION((110+3*nnu+3*nez+1)) :: d_trans_data  ! 64 bit real array of transport keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line

INTEGER                          :: i             ! do i ndex
INTEGER                          :: k             ! neutrino eneregy index
INTEGER                          :: n             ! neutrino flavor index

INTEGER                          :: int           ! integer data variable to read in an interger datum

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl1           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl2           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl3           ! 64 bit real variable read in a real datum
REAL(KIND=double)                :: rl4           ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (a128)
  105 FORMAT (1x,a128)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)
  115 FORMAT (28x,2a)
  117 FORMAT (1x,a6,22x,a2,42x,a16)
  121 FORMAT (10x,2i10)
  123 FORMAT (1x,a6,4x,2i10,42x,a16)
  131 FORMAT (20x,i10,5x,e15.8)
  133 FORMAT (1x,a6,14x,i10,5x,es15.8,22x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  171 FORMAT (10x,i10,2(5x,e15.8))
  173 FORMAT (1x,a6,4x,i10,2(5x,es15.8),12x,a16)
  181 FORMAT (10x,i10,4e10.2)
  183 FORMAT (1x,a6,4x,i10,4(es10.2),12x,a16)
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)

!-----------------------------------------------------------------------
!
!                   \\\\\ TRANSPORT CONTROLS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  inutrn : Transport switch.
!     inutrn = 0 : neutrino transport bypassed.
!     inutrn = 1 : neutrino transport included.
!-----------------------------------------------------------------------

INTEGER                                         :: inutrn

!-----------------------------------------------------------------------
!  iyenu : transport electron fraction change switch.
!
!     iyenu = 0: electron fraction change due to n-type neutrinos bypassed.
!     iyenu = 1: electron fraction change due to n-type neutrinos included.
!-----------------------------------------------------------------------

INTEGER                                         :: iyenu

!-----------------------------------------------------------------------
!  itnu(n) : transport temperature change switch.
!
!     itnu(n) = 0 : temperature change due to n-type neutrinos bypassed.
!     itnu(n) = 1 : temperature change due to n-type neutrinos included.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(nnu)                         :: itnu

!-----------------------------------------------------------------------
!  jnumin : inner radial zone for which neutrino transport is computed.
!  jnumax : outer radial zone for which neutrino transport is computed.
!-----------------------------------------------------------------------

INTEGER                                         :: jnumin
INTEGER                                         :: jnumax

!-----------------------------------------------------------------------
!  idiff : Transport boundary conditions toggle.
!
!     idiff = 0 : dc(j,k,n) computed normally at each iteration.
!     idiff = 1 : dc(j,k,n) computed normally only at the first iteration.
!     idiff = 2 : dc(j,k,n) = 0 for j = jnumax, dc(j,k,n) computed normally
!      at each iteration for  j ne jnumax.
!     idiff = 3 : dc(j,k,n) = 0 for j = jnumax, dc(j,k,n) computed normally
!      only at the first iteration for j ne jnumax.
!     idiff = 4 : dc(j,k,n) = 0 for all j.
!     idiff = 5 : dc(j,k,n) computed normally only at the first iteration.
!-----------------------------------------------------------------------

INTEGER                                         :: idiff

!-----------------------------------------------------------------------
!  ireltrns : Newtonian - GR transport toggle.
!-----------------------------------------------------------------------

INTEGER                                         :: ireltrns

!-----------------------------------------------------------------------
!
!            \\\\\ ITERATION AND TOLERANCE CRITERIA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  tolnut: Temperature convergence parameter for neutrino transport. 
!   The criterion fortemperature convergence is that
!
!     abs(dt/t) < tolnut .
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: tolnut

!-----------------------------------------------------------------------
!  tolnuye : Electron fraction convergence parameter for neutrino transport.
!   The criterion for electron fraction convergence is that
!
!     abs( dye/ye ) < tolnuye .
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: tolnuye

!-----------------------------------------------------------------------
!  tolnupsi, tolpsimin : zero-monent neutrino occupation distribution
!   convergence parameters for neutrino transport. The criteria for neutrino
!   occupation convergence is that
!
!     abs( dpsi0/max( psi0, tolpsimin ) < tolnupsi .
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: tolnupsi
REAL(KIND=double)                               :: tolpsimin

!-----------------------------------------------------------------------
!  iternu : maximum number of iterations to attempt in order to obtain
!   convergence of the neutrino  transport variables (i.e., t, ye, and
!   the psi0's). If iternu = 1, the variables are assumed to have converged
!   after the first iteration attempt.                                          c
!-----------------------------------------------------------------------

INTEGER                                         :: iternu

!-----------------------------------------------------------------------
!  itfail : 0 - iteration failure does not stop calculation.
!  itfail : 1 - iteration failure stops calculation.
!-----------------------------------------------------------------------

INTEGER                                         :: itfail

!-----------------------------------------------------------------------
!  a_prec : recompute neutrino rates if
!
!     dabs( agr(j) - agr_prev(j) )/agr(j) > a_prec
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: a_prec

!-----------------------------------------------------------------------
!
!              \\\\\ TIME AND TIME STEP CONTROLS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  dtnph_trans : the time step for absorption, emission, and scattering,
!   set by the neutrino occupation distributions and variables affected by
!   the neutrino-matter interactions, between time cycle m and time cycle m + 1.
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: dtnph_trans

!-----------------------------------------------------------------------
!  dtnmhm_trans : the coordinate time step for absorption, emission, andx
!   scattering between time cycle m - 1 and time cycle m.
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: dtnmhn_trans

!-----------------------------------------------------------------------
!  dpsisp(n), psimin(n) : parameters used in determining the psi0 change
!   time step due to absorption, emission and transport.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nnu)               :: psimin

!-----------------------------------------------------------------------
!  dpsisp(n), psipmin(n) : parameters used in determining the psi0 change
!   time step due to inelastic scattering and pair production.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(nnu)               :: psipmin

!-----------------------------------------------------------------------
!  tcntrl(i) : numerical criterion for the ith time step control.
!
!     tcntrl(10+n): n-neutrino net temperature change time step criterion,
!      i.e., the maximum permitted
!          abs( tf(j) - ti(j) )/ti(j),
!      due to n-neutrinos. (n=1,2,3,4).
!
!     tcntrl(15+n): n-neutrino net electron fraction change time step
!      criterion, i.e., the maximum permitted
!          abs( yef(j,n,i_ray) - yei(j,n,i_ray) )/yei(j) ),
!      where dye(j,1,i_ray) is the change of ye in zone j due to all neutrinos.
!
!     tcntrl(20+n): n-neutrino zero-moment change time step criterion due
!      to absorption, emission, scattering, production and transport of
!      n-neutrinos, i.e., the maximum permitted
!          abs( psi0f(j,k,n) - psi0i(j,k,n) )/max( psi0f(j,k,n), tolpsimin )
!
!     tcntrl(49): maximum increase in the 'neutrino transport' time step,
!      i.e., the maximum allowed value of dtnphn/dtnmhn.
!
!     tcntrl(50): time step is set equal to tcntrl(50) if jdt(50) = -1.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(50)                :: tcntrl

!-----------------------------------------------------------------------
!  rdtmax : dtnphn > dtnph: neutrino transport is bypassed until the
!   accumulated time since the last neutrino transport update could
!   exceed dtnphn by a specified amount in the subsequent time step.
!   Thus, if delt is the accumulated time from the last implementation of
!   neutrino transport to the beginning of the current time cycle, neutrino
!   transport is implemented in the current time cycle with dtnphn set equal
!   to the the accumulated time if
!
!     delt + rdtmax*dtnph > dtnphn .
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: rdtmax

!-----------------------------------------------------------------------
!  intnu_trns : hydro subcycling parameters for the preceding time step.
!   i.e., the number of hydro updates per neutrino process i update at
!   the last neutrino transport update.
!-----------------------------------------------------------------------

INTEGER                                        :: intnu_trns

!-----------------------------------------------------------------------
!  ncynu_trns : the number of hydro updates since the last neutrino process i update.
!-----------------------------------------------------------------------

INTEGER                                        :: ncynu_trns

!-----------------------------------------------------------------------
!  intnur_trns : hydro subcycling parameters for the preceding time step. i.e.,
!   the number of hydro updates per neutrino process i update at the last
!   neutrino transport update.
!-----------------------------------------------------------------------

INTEGER                                        :: intnur_trns

!-----------------------------------------------------------------------
!
!               \\\\\ NEUTRINO GROUP ENERGIES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  nnugp(n) : the number of energy zones for neutrinos of type n.
!
!  unumn    : the minimum value of the zone-centered neutrino energy (MeV). 
!  unumx    : the maximum value of the zone-centered neutrino energy (MeV).
!
!  unui(k)  : the value of the zone-centered neutrino energy at infinity of
!              energy zone k (MeV).
!  dunui(k) : the width of energy zone k at infinity (MeV),
!              dunui(k) = unubi(k+1) - unubi(k) 
!  unubi(k) : the value of the inner zone-edge neutrino energy at infinity
!              of energy zone k (MeV).
!
!  stwt(n)  : the statistical weight of neutrinos of type n, i.e., the
!   number of distinct multiplicities of neutrinos subsumed as neutrinos
!   of type n. (Typically, stwt(1) = 1., stwt(2) = 1., stwt(3) = 2 stwt(4) = 2.)
!-----------------------------------------------------------------------

INTEGER, DIMENSION(nnu)                        :: nnugp

REAL(KIND=double)                              :: unumn
REAL(KIND=double)                              :: unumx

REAL(KIND=double), DIMENSION(nez)              :: unui
REAL(KIND=double), DIMENSION(nez)              :: dunui
REAL(KIND=double), DIMENSION(nezp1)            :: unubi

REAL(KIND=double), DIMENSION(nnu)              :: stwt

!-----------------------------------------------------------------------
!
!                    \\\\\ SOURCE CONTROLS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Absorption and emission of e-neutrinos on nuclei (FFN)
!-----------------------------------------------------------------------
!  iaence : emission and absorption of e-neutrinos on nuclei as given by
!   the FFN prescription switch
!
!     iaence = 0 : emission and absorption of e-neutrinos on nuclei as
!      given by the FFN  prescription turned off.
!     iaence = 1 : emission and absorption of e-neutrinos on nuclei as
!      given by the FFN  prescription turned on.
!
!  edmpe :  difference between the mean excited state energy of the daughter
!   nucleus and that ofthe parent nucleus (MeV).
!-----------------------------------------------------------------------

INTEGER                                         :: iaence

REAL(KIND=double)                               :: edmpe

!-----------------------------------------------------------------------
!  Absorption and emission of e-antineutrinos on nuclei (FFN)
!-----------------------------------------------------------------------
!      iaenca: emission and absorption of e-antineutrinos on nuclei as
!       given by the FFN prescription switch
!
!          iaenca = 0: emission and absorption of e-antineutrinos on nuclei
!           as given by the FFN prescription turned off.
!          iaenca = 1: emission and absorption of e-antineutrinos on nuclei
!           as given by the FFN prescription turned on.
!
!      edmpa: difference between the mean excited state energy of the
!       daughter nucleus and that of the parent nucleus (MeV).
!-----------------------------------------------------------------------

INTEGER                                         :: iaenca

REAL(KIND=double)                               :: edmpa

!-----------------------------------------------------------------------
!  Absorption and emission of e-neutrinos and e-antineutrinos on nucleons
!-----------------------------------------------------------------------
!  iaefnp : emission and absorption of e-neutrinos and e-antineutrinos on
!   free neutrons and protons switch.
!
!     iaefnp = 0 : emission and absorption of e-neutrinos and e-antineutrinos
!      on free neutrons and protons turned off.
!     iaefnp = 1 : emission and absorption of e-neutrinos and e-antineutrinos
!      on free neutrons and protons turned on.
!
!      rhoaefnp : density above which emission and absorption of e-neutrinos
!       and e-antineutrinos on free neutrons and protons is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: iaefnp

REAL(KIND=double)                               :: rhoaefnp

!-----------------------------------------------------------------------
!  Absorption and emission of e-neutrinos and e-antineutrinos on nuclei
!   (Haxton)
!-----------------------------------------------------------------------
!  iaencnu : emission and absorption of e-neutrinos and e-antineutrinos on
!   nuclei as given by Haxton's prescription
!
!  iaencnu = 0 : emission and absorption of e-neutrinos and e-antineutrinos
!   on nuclei as given by Haxton's prescription turned off.
!  iaencnu = 1 : emission and absorption of e-neutrinos and e-antineutrinos
!   on nuclei as given by Haxton's prescription turned on.
!
!  roaencnu : density above which emission and absorption of e-neutrinos
!   and e-antineutrinos on nuclei as given by Haxton is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: iaencnu

REAL(KIND=double)                               :: roaencnu

!-----------------------------------------------------------------------
!  Absorption and emission of e-neutrinos on nuclei
!   (NSE-folded tabular data)
!-----------------------------------------------------------------------
!  iaenct: emission and absorption of e-neutrinos on nuclei as given by
!   the NSE-folded tabular data for Hix et al (2003) 
!
!     iaenct  = 0: emission and absorption of e-neutrinos on nuclei as
!      given by the tabular data is turned off.
!     iaenct  = 1: emission and absorption of e-neutrinos on nuclei as
!      given by the tabular data is turned on.
!
!  roaenct: density above which emission and absorption of e-neutrinos on
!   nuclei as given by the is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: iaenct

REAL(KIND=double)                               :: roaenct

!-----------------------------------------------------------------------
!  Inclusion of recoll, thermal motions, and nucleon blocking factors
!   switch
!-----------------------------------------------------------------------
!  i_aeps : emission and aborption of e-neutrinos and e-antineutrinos
!   on free neutrons and protons phase space switch.
!
!     i_aeps = 0: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons are computed
!      assuming no recoil, or thermal motions, and approximate
!      nucleon phase space blocking.
!     iaefnp = 1: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons are computed
!      taking into account recoil, thermal motions, and nucleon
!      phase space blocking.
!-----------------------------------------------------------------------

INTEGER                                         :: i_aeps

!-----------------------------------------------------------------------
!  Weak magnetism correction for charged current interactions
!-----------------------------------------------------------------------
!  icc_wm : charged current weak magnetism switch.
!
!     icc_wm = 0 : weak magnetism correction for charged current
!      interactions turned off.
!     icc_wm = 1 : weak magnetism correction for charged current
!      interactions turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: icc_wm

!-----------------------------------------------------------------------
!  General scattering switch
!-----------------------------------------------------------------------
!  iscat : general scattering switch.
!
!     iscat = 0 : all scattering processes turned off.
!     iscat = 1 : any scattering process is allowed to proceed if its key
!      is turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: iscat

!-----------------------------------------------------------------------
!  Neutrino nucleon isoenergetic scattering controls
!-----------------------------------------------------------------------
!  in : neutrino-neutron scattering switch
!
!     in = 0 : neutrino-neutron isoenergetic scattering turned off.
!     in = 1 : neutrino-neutron isoenergetic scattering turned on.
!
!  ip : neutrino-proton scattering switch
!
!     ip = 0 : neutrino-proton isoenergetic scattering turned off.
!     ip = 1 : neutrino-proton isoenergetic scattering turned on.
!
!  ietann : final state blocking switch
!
!     ietann = 0 : final state blocking due to neutron or proton final state
!      occupancy in neutrino-neutron and neutrino-proton isoenergetic
!      scattering turned off.
!     ietann = 1 : final state blocking due to neutron or proton final state
!      occupancy in neutrino-neutron and neutrino-proton isoenergetic
!      scattering turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: in
INTEGER                                         :: ip
INTEGER                                         :: ietann

!-----------------------------------------------------------------------
!  Neutrino nuclei isoenergetic scattering controls
!-----------------------------------------------------------------------
!  ihe : neutrino-helium isoenergetic scattering switch
!
!     ihe = 0 : neutrino-helium isoenergetic scattering turned off.
!     ihe = 1 : neutrino-helium isoenergetic scattering turned on.
!
!  iheavy : neutrino-heavy nucleus isoenergetic scattering switch
!
!     iheavy = 0 : neutrino-heavy nucleus isoenergetic scattering turned off.
!     iheavy = 1 : neutrino-heavy nucleus isoenergetic scattering turned on.
!
!  iicor : ion-ion correlation switch.
!
!     iicor = 0 : neutrino-heavy nucleus isoenergetic scattering correction
!      due to ion-ion correlation effects turned off.
!     iicor = 1 : neutrino-heavy nucleus isoenergetic scattering correction
!      due to ion-ion correlation effects turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: ihe
INTEGER                                         :: iheavy
INTEGER                                         :: iicor

!-----------------------------------------------------------------------
!  Weak magnetism correction for neutral current interactions
!-----------------------------------------------------------------------
!  inc_wm : neutral current weak magnetism switch.
!
!     inc_wm = 0 : weak magnetism correction for neutral current
!      interactions turned off.
!     inc_wm = 1 : weak magnetism correction for neutral current
!      interactions turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: inc_wm

!-----------------------------------------------------------------------
!  Neutrino-electron scattering controls
!-----------------------------------------------------------------------
!  nes : neutrino-electron scattering switch.
!
!     nes = 0 : neutrino-electron scattering turned off.
!     nes = 1 : neutrino-electron scattering turned on.
!
!  rhonesmn : density below which n-neutrino-electron scattering is turned off.
!
!  rhonesmx : density above which n-neutrino-electron scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: nes

REAL(KIND=double)                               :: rhonesmn
REAL(KIND=double)                               :: rhonesmx

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering controls
!-----------------------------------------------------------------------
!  nncs : neutrino-nucleus inelastic scattering switch.
!
!     nncs = 0 : neutrino-nucleus inelastic scattering turned off.
!     nncs = 1 : neutrino-nucleus inelastic scattering turned on.
!
!  rhonncsmn : density below which n-neutrino-electron scattering is turned off.
!
!  rhonncsmx : density above which n-neutrino-electron scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: nncs

REAL(KIND=double)                               :: rhonncsmn
REAL(KIND=double)                               :: rhonncsmx

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation controls
!-----------------------------------------------------------------------
!  ipair : electron-positron pair annihilation switch.
!
!     ipair = 0 : electron-positron pair annihilation process turned off.
!     ipair = 1 : electron-positron pair annihilation process turned on.
!
!  rhopairemn : density below which e-neutrino - e-antineutrino pair
!   production by electron-positron pair annihilation is turned off.
!
!  rhopairemx : density above which e-neutrino - e-antineutrino pair
!   production by electron-positron pair annihilation is turned off.
!
!  rhopairtmn : density below which x-neutrino - x-antineutrino pair
!   production by electron-positron pair annihilation is turned off.
!
!  rhopairtmx : density above which x-neutrino - x-antineutrino pair
!   production by electron-positron pair annihilation is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: ipair

REAL(KIND=double)                               :: rhopairemn
REAL(KIND=double)                               :: rhopairemx
REAL(KIND=double)                               :: rhopairtmn
REAL(KIND=double)                               :: rhopairtmx

!-----------------------------------------------------------------------
!  Bremsstrahlung pair annihilation controls
!-----------------------------------------------------------------------
!  ibrem : bremsstrahlung pair annihilation switch.
!
!     ibrem = 0 : bremsstrahlung pair annihilation process turned off.
!     ibrem = 1 : bremsstrahlung pair annihilation process turned on.
!
!  rhobrememn : density below which e-neutrino - e-antineutrino pair
!   production by bremsstrahlung pair annihilation is turned off.
!
!  rhobrememx : density above which e-neutrino - e-antineutrino pair
!   production by bremsstrahlung pair annihilation is turned off.
!
!  rhobremtmn : density below which x-neutrino - x-antineutrino pair
!   production by bremsstrahlung pair annihilation is turned off..
!
!  rhobremtmx : density above which x-neutrino - x-antineutrino pair
!   production by bremsstrahlung pair annihilation is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: ibrem

REAL(KIND=double)                               :: rhobrememn
REAL(KIND=double)                               :: rhobrememx
REAL(KIND=double)                               :: rhobremtmn
REAL(KIND=double)                               :: rhobremtmx

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering controls
!-----------------------------------------------------------------------
!  isctn: neutrino-nucleon inelastic scattering switch.
!
!     isctn = 0: neutrino-nucleon inelastic scattering turned off.
!     isctn = 1: neutrino-nucleon inelastic scattering turned on.
!
!  rhosctnemn: density below which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon  elastic scattering is turned off.
!
!  rhosctnemx: density above which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon elastic scattering is turned off.
!
!  rhosctntmn: density below which x-neutrino-neutrino-nucleon,
!   x-antineutrino-neutrino-nucleon  elastic scattering is turned off.
!
!  rhosctntmx: density above which x-neutrino-neutrino-nucleon,
!   x-antineutrino-neutrino-nucleon elastic scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: isctn

REAL(KIND=double)                               :: rhosctnemn
REAL(KIND=double)                               :: rhosctnemx
REAL(KIND=double)                               :: rhosctntmn
REAL(KIND=double)                               :: rhosctntmx

!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic scattering controls
!-----------------------------------------------------------------------
!  isctnn: neutrino-nucleon inelastic scattering switch.
!
!     isctnn = 0: neutrino-nucleon inelastic scattering turned off.
!     isctnn = 1: neutrino-nucleon inelastic scattering turned on.
!
!  rhosctnnemn: density below which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon inelastic scattering is turned off.
!
!  rhosctnnemx: density above which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon inelastic scattering is turned off.
!
!  rhosctnntmn: density below which x-neutrino-neutrino-nucleon,
!   x-antineutrino-neutrino-nucleon inelastic scattering is turned off.
!
!  rhosctnntmx: density above which x-neutrino-neutrino-nucleon,
!   x-antineutrino-neutrino-nucleon scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: isctnn

REAL(KIND=double)                               :: rhosctnnemn
REAL(KIND=double)                               :: rhosctnnemx
REAL(KIND=double)                               :: rhosctnntmn
REAL(KIND=double)                               :: rhosctnntmx

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering controls
!-----------------------------------------------------------------------
!  isctnn : neutrino-nucleus inelastic scattering switch.
!
!     isctnn = 0 : neutrino-nucleus inelastic scattering turned off.
!     isctnn = 1 : neutrino-nucleus inelastic scattering turned on.
!
!  rhosctnnemn : density below which e-neutrino-neutrino-nucleus,
!   e-antineutrino-neutrino-nucleus  inelastic scattering is turned off.
!
!  rhosctnnemx : density above which e-neutrino-neutrino-nucleus,
!   e-antineutrino-neutrino-nucleus inelastic scattering is turned off.
!
!  rhosctnntmn : density below which x-neutrino-neutrino-nucleus,
!   x-antineutrino-neutrino-nucleus  inelastic scattering is turned off.
!
!  rhosctnntmx : density above which x-neutrino-neutrino-nucleus,
!   x-antineutrino-neutrino-nucleus scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: isctnA

REAL(KIND=double)                               :: rhosctnAemn
REAL(KIND=double)                               :: rhosctnAemx
REAL(KIND=double)                               :: rhosctnAtmn
REAL(KIND=double)                               :: rhosctnAtmx

!-----------------------------------------------------------------------
!
!                \\\\\ EQUATION OF STATE BORDERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  rhopnu : density above which neutrino pressure is computed isotropically
!   as a equilibrated relativistic gas.
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: rhopnu


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!               \\\\\ INITIALIZE TRANSPORT KEYS /////
!
!-----------------------------------------------------------------------

iyenu                             = 0
itnu                              = 0
jnumin                            = 0
jnumax                            = 0
idiff                             = 0
inutrn                            = 0
ireltrns                          = 0
iternu                            = 30
itfail                            = 0

tolnut                            = 1.d-02
tolnuye                           = 1.d-02
tolnupsi                          = 1.d-01
tolpsimin                         = 1.d-01

a_prec                            = zero
dtnph_trans                       = 1.d+20
dtnmhn_trans                      = 1.d+20
psimin                            = 1.d-01
psipmin                           = 1.d-01

tcntrl                            = 1.d+20
rdtmax                            = 1.2d+00

intnu_trns                        = 1
ncynu_trns                        = 0
intnur_trns                       = 1

nnugp                             = 0

unumn                             = zero
unumx                             = zero
unui                              = zero
dunui                             = zero
unubi                             = zero
stwt                              = zero

iaefnp                            = 0
iaence                            = 0
iaenca                            = 0
iaencnu                           = 0
iaenct                            = 0
i_aeps                            = 0
icc_wm                            = 0

rhoaefnp                          = 1.d+15
edmpe                             = zero
edmpa                             = zero
roaencnu                          = 1.d+15
roaenct                           = 1.d+13

iscat                             = 1

in                                = 0
ip                                = 0
ietann                            = 0

ihe                               = 0
iheavy                            = 0
iicor                             = 0
inc_wm                            = 0

nes                               = 0
rhonesmn                          = 1.d+7
rhonesmx                          = 1.d+15

nncs                              = 0
rhonncsmn                         = 1.d+7
rhonncsmx                         = 1.d+15

ipair                             = 0
rhopairemn                        = 1.d+7
rhopairemx                        = 1.d+15
rhopairtmn                        = 1.d+7
rhopairtmx                        = 1.d+15

ibrem                             = 0
rhobrememn                        = 1.d+7
rhobrememx                        = 1.d+15
rhobremtmn                        = 1.d+7
rhobremtmx                        = 1.d+15

isctn                             = 0
rhosctnemn                        = 1.d+7
rhosctnemx                        = 1.d+15
rhosctntmn                        = 1.d+7
rhosctntmx                        = 1.d+15

isctnn                            = 0
rhosctnnemn                       = 1.d+7
rhosctnnemx                       = 1.d+15
rhosctnntmn                       = 1.d+7
rhosctnntmx                       = 1.d+15

isctnA                            = 0
rhosctnAemn                       = 1.d+7
rhosctnAemx                       = 1.d+15
rhosctnAtmn                       = 1.d+7
rhosctnAtmx                       = 1.d+15

rhopnu                            = 1.d+15


!-----------------------------------------------------------------------
!
!                 \\\\\ READ TRANSPORT KEYS /////
!
!-----------------------------------------------------------------------

REWIND (nreadp)

READ: DO

!-----------------------------------------------------------------------
!  Read line
!-----------------------------------------------------------------------

  READ (nreadp,101,END=5000) line

  type               = line(1:6)
  name               = line(73:88)

!-----------------------------------------------------------------------
!  Ignore comment cards in the data stream
!-----------------------------------------------------------------------

  IF ( type == 'cccccc'  .or.  type(1:1) == '!'  .or.  type == '      ') THEN
    IF ( iskip == 0 ) WRITE (nprint,105) line
    CYCLE
  END IF ! type = 'cccccc' or type = '!'

!-----------------------------------------------------------------------
!
!                   \\\\\ TRANSPORT CONTROLS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  iyenu
!-----------------------------------------------------------------------

  IF ( type == 'iyenu ' ) THEN
    READ (line ,111) iyenu
    IF ( iskip == 0 ) WRITE (nprint,113) type,iyenu,name
    CYCLE
  END IF ! type = 'iyenu '

!-----------------------------------------------------------------------
!  itnu(n)
!-----------------------------------------------------------------------

  IF ( type == 'itnu  ' ) THEN
    READ (line ,121) n,int
    itnu(n)          = int
    IF ( iskip == 0 ) WRITE (nprint,123) type,n,itnu(n),name
    CYCLE
  END IF ! type = 'itnu  '

!-----------------------------------------------------------------------
!  jnu
!-----------------------------------------------------------------------

  IF ( type == 'jnu   ' ) THEN
    READ (line ,111) int

    IF ( name == 'jnumin  ' ) THEN
      jnumin         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jnumin,name
      CYCLE
    END IF ! name = 'jnumin  '

    IF ( name == 'jnumax  ' ) THEN
      jnumax         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jnumax,name
      CYCLE
    END IF ! name = 'jnumax  '

  END IF ! type = 'jnu   '

!-----------------------------------------------------------------------
!  idiff
!-----------------------------------------------------------------------

  IF ( type == 'idiff ' ) THEN
    READ (line ,111) idiff
    IF ( iskip == 0 ) WRITE (nprint,113) type,idiff,name
     CYCLE
  END IF ! type = 'idiff '

!-----------------------------------------------------------------------
!  inutrn
!-----------------------------------------------------------------------

  IF ( type == 'inutrn' ) THEN
    READ (line ,111) inutrn
    IF ( iskip == 0 ) WRITE (nprint,113) type,inutrn,name
    CYCLE
  END IF ! type = 'inutrn'

!-----------------------------------------------------------------------
!  irel
!-----------------------------------------------------------------------

  IF ( type == 'irel  ' ) THEN
    READ (line ,111) int

    IF ( name == 'ireltrns' ) THEN
      ireltrns       = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ireltrns,name
      CYCLE
    END IF ! name = 'ireltrns'

  END IF ! type = 'irel  '

!-----------------------------------------------------------------------
!
!            \\\\\ ITERATION AND TOLERANCE CRITERIA /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  tolnu
!-----------------------------------------------------------------------

  IF ( type == 'tolnu ' ) THEN

    IF ( name == 'tolnut  ') THEN
      READ (line ,141) tolnut
      IF ( iskip == 0 ) WRITE (nprint,143) type,tolnut,name
      CYCLE
    END IF ! name = 'tolnut  '

    IF ( name == 'tolnuye ') THEN
      READ (line ,141) tolnuye
      IF ( iskip == 0 ) WRITE (nprint,143) type,tolnuye,name
      CYCLE
    END IF ! name = 'tolnuye '

    IF ( name == 'tolnupsi') THEN
      READ (line ,141) tolnupsi
      IF ( iskip == 0 ) WRITE (nprint,143) type,tolnupsi,name
      CYCLE
    END IF ! name = 'tolnupsi'

    IF ( name == 'tolpsimin') THEN
      READ (line ,141) tolpsimin
      IF ( iskip == 0 ) WRITE (nprint,143) type,tolpsimin,name
      CYCLE
    END IF ! name = 'tolpsimin'

  END IF ! type = 'tolnu '

!-----------------------------------------------------------------------
!  iternu
!-----------------------------------------------------------------------

  IF ( type == 'iternu' ) THEN

    IF ( name == 'iternu  ') THEN
      READ (line ,111) iternu
      IF ( iskip == 0 ) WRITE (nprint,113) type,iternu,name
      CYCLE
    END IF ! name = 'iternu  '

    IF ( name == 'itfail  ' ) THEN
      READ (line,111) itfail
      IF ( iskip == 0 ) WRITE (nprint,113) type,itfail,name
      CYCLE
    END IF ! name = 'itfail  '

  END IF ! type = 'iternu'

  IF ( type == 'a_prec' ) THEN
    READ (line ,141) rl
    a_prec           = rl
    IF ( iskip == 0 ) WRITE (nprint,143) type,a_prec,name
    CYCLE
  END IF ! type = 'a_prec'

!-----------------------------------------------------------------------
!
!              \\\\\ TIME AND TIME STEP CONTROLS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  sptime
!-----------------------------------------------------------------------

  IF ( type == 'sptime' ) THEN
    READ (line ,141) rl

    IF ( name == 'dtnph_trans') THEN
      dtnph_trans    = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtnph_trans,name
      CYCLE
    END IF ! name = 'dtnph_trans'

    IF ( name == 'dtnmhn_trans') THEN
      dtnmhn_trans    = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtnmhn_trans,name
      CYCLE
    END IF ! name = 'dtnmhn_trans'

  END IF ! type = 'sptime'

!-----------------------------------------------------------------------
!  psitol
!-----------------------------------------------------------------------

  IF ( type == 'psitol' ) THEN

    IF ( name == 'psimin') THEN
      READ (line ,131) n,rl
      psimin(n)  = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,psimin(n),name
      CYCLE
    END IF ! name = 'psimin'

    IF ( name == 'psipmin') THEN
      READ (line ,131) n,rl
      psipmin(n)  = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,psipmin(n),name
      CYCLE
    END IF ! name = 'psipmin'

  END IF ! type = 'psitol'

!-----------------------------------------------------------------------
!  tcntrl
!-----------------------------------------------------------------------

  IF ( type == 'tcntrl' ) THEN

    IF ( name == 'tcntrl_trans' ) THEN
      READ (line ,131) n,rl
      tcntrl(n)     = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,tcntrl(n),name
      CYCLE
    END IF ! name = 'tcntrl_trans'

    IF ( name == 'rdtmax  ' ) THEN
      READ (line ,141) rl
      rdtmax         = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,rdtmax,name
      CYCLE
    END IF ! name = 'rdtmax  '

  END IF ! type = 'tcntrl'

!-----------------------------------------------------------------------
!  Multiple time step controls
!-----------------------------------------------------------------------

  IF ( type == 'intnu ' ) THEN
    READ (line ,111) int

    IF ( name == 'intnu_trns') THEN
      intnu_trns = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,intnu_trns,name
      CYCLE
    END IF ! name = 'intnu_trns'

    IF ( name == 'ncynu_trns') THEN
      ncynu_trns = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ncynu_trns,name
      CYCLE
    END IF ! name = 'ncynu_sp'

    IF ( name == 'intnur_trns') THEN
      intnur_trns= int
      IF ( iskip == 0 ) WRITE (nprint,113) type,intnur_trns,name
      CYCLE
    END IF ! name = 'intnur_trns'

  END IF ! type = 'intnu '

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO DISTRIBUTION ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  unu
!-----------------------------------------------------------------------

  IF ( type == 'unu   ' ) THEN

    IF ( name == 'nnugp   ') THEN
      READ (line ,121) n,int
      nnugp(n)       = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,nnugp(n),name
      CYCLE
    END IF ! name = 'nnugp   '

    IF ( name == 'unumn   ') THEN
      READ (line ,141) unumn
      IF ( iskip == 0 ) WRITE (nprint,143) type,unumn,name
      CYCLE
    END IF ! name = 'unumn   '

    IF ( name == 'unumx   ') THEN
      READ (line ,141) unumx
      IF ( iskip == 0 ) WRITE (nprint,143) type,unumx,name
      CYCLE
    END IF ! name = 'unumx   '

    IF ( name == 'unu     ') THEN
      READ (line ,131) k,rl
      unui(k)        = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,k,unui(k),name
      CYCLE
    END IF ! name = 'unu     '

    IF ( name == 'dunu    ') THEN
      READ (line ,131) k,rl
      dunui(k)       = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,k,dunui(k),name
      CYCLE
    END IF ! name = 'dunu    '

    IF ( name == 'unub    ') THEN
      READ (line ,131) k,rl
      unubi(k)       = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,k,unubi(k),name
      CYCLE
    END IF ! name = 'unub    '

  END IF

!-----------------------------------------------------------------------
!  stwt
!-----------------------------------------------------------------------

  IF ( type == 'stwt  ' ) THEN
    READ (line ,131) n,rl
    stwt(n)          = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,n,stwt(n),name
    CYCLE
  END IF ! type = 'stwt  '

!-----------------------------------------------------------------------
!
!                    \\\\\ SOURCE CONTROLS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  aence
!-----------------------------------------------------------------------

  IF ( type == 'aence ' ) THEN
    READ (line ,131) int,rl
    iaence            = int
    edmpe             = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,iaence,edmpe,name
    CYCLE
  END IF ! type = 'aence '

!-----------------------------------------------------------------------
!  aenca
!-----------------------------------------------------------------------

  IF ( type == 'aenca ' ) THEN
    READ (line ,131) int,rl
    iaenca            = int
    edmpa             = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,iaenca,edmpa,name
    CYCLE
  END IF ! type = 'aenca '

!-----------------------------------------------------------------------
!  iaefnp
!-----------------------------------------------------------------------

  IF ( type == 'iaefnp' ) THEN
    READ (line ,131) int,rl
    iaefnp            = int
    rhoaefnp          = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,iaefnp,rhoaefnp,name
    CYCLE
  END IF ! type = 'iaefnp'

!-----------------------------------------------------------------------
!  aencnu
!-----------------------------------------------------------------------

  IF ( type == 'aencnu' ) THEN
    READ (line ,131) int,rl
    iaencnu           = int
    roaencnu          = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,iaencnu,roaencnu,name
    CYCLE
  END IF ! type = 'aencnu'

!-----------------------------------------------------------------------
!  aenct
!-----------------------------------------------------------------------

  IF ( type == 'aenct' ) THEN
    READ (line ,131) int,rl
    iaenct            = int
    roaenct           = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,iaenct,roaenct,name
    CYCLE
  END IF ! type = 'aenct'

!-----------------------------------------------------------------------
!  i_aeps
!-----------------------------------------------------------------------

  IF ( type == 'i_aeps' ) THEN
    READ (line ,111) i_aeps
    IF ( iskip == 0 ) WRITE (nprint,113) type,i_aeps,name
    CYCLE
  END IF ! type = 'i_aeps'

!-----------------------------------------------------------------------
!  icc_wm
!-----------------------------------------------------------------------

  IF ( type == 'icc_wm' ) THEN
    READ (line ,111) icc_wm
    IF ( iskip == 0 ) WRITE (nprint,113) type,icc_wm,name
    CYCLE
  END IF ! type = 'icc_wm'

!-----------------------------------------------------------------------
!  iscat
!-----------------------------------------------------------------------

  IF ( type == 'iscat ' ) THEN
    READ (line ,111) iscat
    IF ( iskip == 0 ) WRITE (nprint,113) type,iscat,name
    CYCLE
  END IF ! type = 'iscat '

!-----------------------------------------------------------------------
!  in
!-----------------------------------------------------------------------

  IF ( type == 'in    ' ) THEN
    READ (line ,111) in
    IF ( iskip == 0 ) WRITE (nprint,113) type,in,name
    CYCLE
  END IF ! type = 'in    '

!-----------------------------------------------------------------------
!  ip
!-----------------------------------------------------------------------

  IF ( type == 'ip    ' ) THEN
    READ (line ,111) ip
    IF ( iskip == 0 ) WRITE (nprint,113) type,ip,name
    CYCLE
  END IF ! type = 'ip    '

!-----------------------------------------------------------------------
!  ihe
!-----------------------------------------------------------------------

  IF ( type == 'ihe   ' ) THEN
    READ (line ,111) ihe
    IF ( iskip == 0 ) WRITE (nprint,113) type,ihe,name
    CYCLE
  END IF ! type = 'ihe   '

!-----------------------------------------------------------------------
!  iheavy
!-----------------------------------------------------------------------

  IF ( type == 'iheavy' ) THEN
    READ (line ,111) iheavy
    IF ( iskip == 0 ) WRITE (nprint,113) type,iheavy,name
    CYCLE
  END IF ! type == 'iheavy'

!-----------------------------------------------------------------------
!  iicor
!-----------------------------------------------------------------------

  IF ( type == 'iicor ' ) THEN
    READ (line ,111) iicor
    IF ( iskip == 0 ) WRITE (nprint,113) type,iicor,name
    CYCLE
  END IF ! type == 'iicor '

!-----------------------------------------------------------------------
!  ietann
!-----------------------------------------------------------------------

  IF ( type == 'ietann' ) THEN
    READ (line ,111) ietann
    IF ( iskip == 0 ) WRITE (nprint,113) type,ietann,name
    CYCLE
  END IF ! type = 'ietann'

!-----------------------------------------------------------------------
!  inc_wm
!-----------------------------------------------------------------------

  IF ( type == 'inc_wm' ) THEN
    READ (line ,111) inc_wm
    IF ( iskip == 0 ) WRITE (nprint,113) type,inc_wm,name
    CYCLE
  END IF ! type = 'inc_wm'

!-----------------------------------------------------------------------
!  nes
!-----------------------------------------------------------------------

  IF ( type == 'nes   ' ) THEN
    READ (line ,171) nes,rhonesmn,rhonesmx
    IF ( iskip == 0 ) WRITE (nprint,173) type,nes,rhonesmn,rhonesmx,name
    CYCLE
  END IF ! type = 'nes   '

!-----------------------------------------------------------------------
!  sncnu
!-----------------------------------------------------------------------

  IF ( type == 'sncnu ' ) THEN
    READ (line ,171) nncs,rhonncsmn,rhonncsmx
    IF ( iskip == 0 ) WRITE (nprint,173) type,nncs,rhonncsmn,rhonncsmx,name
    CYCLE
  END IF ! type = 'sncnu '

!-----------------------------------------------------------------------
!  ipair
!-----------------------------------------------------------------------

  IF ( type == 'ipair ' ) THEN
    READ (line ,181) int,rl1,rl2,rl3,rl4
    ipair            = int
    rhopairemn       = rl1
    rhopairemx       = rl2
    rhopairtmn       = rl3
    rhopairtmx       = rl4
    IF ( iskip == 0 ) WRITE (nprint,183) type,ipair,rhopairemn,rhopairemx, &
&    rhopairtmn,rhopairtmx,name
    CYCLE
  END IF ! type = 'ipair '

!-----------------------------------------------------------------------
!  ibrem
!-----------------------------------------------------------------------

  IF ( type == 'ibrem ' ) THEN
    READ (line ,181) int,rl1,rl2,rl3,rl4
    ibrem            = int
    rhobrememn       = rl1
    rhobrememx       = rl2
    rhobremtmn       = rl3
    rhobremtmx       = rl4
    IF ( iskip == 0 ) WRITE (nprint,183) type,ibrem,rhobrememn,rhobrememx, &
&    rhobremtmn,rhobremtmx,name
    CYCLE
  END IF ! type = 'ibrem '

!-----------------------------------------------------------------------
!  isctn
!-----------------------------------------------------------------------

  IF ( type == 'isctn ' ) THEN
    READ (line ,181) int,rl1,rl2,rl3,rl4
    isctn            = int
    rhosctnemn       = rl1
    rhosctnemx       = rl2
    rhosctntmn       = rl3
    rhosctntmx       = rl4
    IF ( iskip == 0 ) WRITE (nprint,183) type,isctn,rhosctnemn,rhosctnemx, &
&    rhosctntmn,rhosctntmx,name
    CYCLE
  END IF ! type = 'isctn '

!-----------------------------------------------------------------------
!  isctnn
!-----------------------------------------------------------------------

  IF ( type == 'isctnn' ) THEN
    READ (line ,181) int,rl1,rl2,rl3,rl4
    isctnn           = int
    rhosctnnemn      = rl1
    rhosctnnemx      = rl2
    rhosctnntmn      = rl3
    rhosctnntmx      = rl4
    IF ( iskip == 0 ) WRITE (nprint,183) type,isctnn,rhosctnnemn,rhosctnnemx, &
&    rhosctnntmn,rhosctnntmx,name
    CYCLE
  END IF ! type = 'ibrem '

!-----------------------------------------------------------------------
!  isctnA
!-----------------------------------------------------------------------

  IF ( type == 'isctnA' ) THEN
    READ (line ,181) int,rl1,rl2,rl3,rl4
    isctnA           = int
    rhosctnAemn      = rl1
    rhosctnAemx      = rl2
    rhosctnAtmn      = rl3
    rhosctnAtmx      = rl4
    IF ( iskip == 0 ) WRITE (nprint,183) type,isctnA,rhosctnAemn,rhosctnAemx, &
&    rhosctnAtmn,rhosctnAtmx,name
    CYCLE
  END IF ! type = 'ibrem '

!-----------------------------------------------------------------------
!  rhopnu
!-----------------------------------------------------------------------

  IF ( type == 'rhopnu' ) THEN
    READ (line ,141) rhopnu
    IF ( iskip == 0 ) WRITE (nprint,143) type,rhopnu,name
    CYCLE
  END IF ! type = 'rhopnu'

!-----------------------------------------------------------------------
!  Unrecognized card
!-----------------------------------------------------------------------

  IF ( nrst == 0 ) THEN
    WRITE (nprint,401)
    WRITE (nprint,403) line
  END IF

END DO READ

!-----------------------------------------------------------------------
!
!                 \\\\\ PACK TRANSPORT KEYS /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE

i_trans_data(1)                   = iyenu
i_trans_data(2)                   = jnumin
i_trans_data(3)                   = jnumax
i_trans_data(4)                   = idiff
i_trans_data(5)                   = inutrn
i_trans_data(6)                   = ireltrns
i_trans_data(7)                   = iternu
i_trans_data(8)                   = itfail
i_trans_data(9)                   = intnu_trns
i_trans_data(11)                  = ncynu_trns
i_trans_data(13)                  = intnur_trns
i_trans_data(15)                  = iaence
i_trans_data(16)                  = iaenca
i_trans_data(17)                  = iaefnp
i_trans_data(18)                  = iaencnu
i_trans_data(19)                  = i_aeps
i_trans_data(20)                  = icc_wm
i_trans_data(21)                  = iscat
i_trans_data(22)                  = in
i_trans_data(23)                  = ip
i_trans_data(24)                  = ietann
i_trans_data(25)                  = ihe
i_trans_data(26)                  = iheavy
i_trans_data(27)                  = iicor
i_trans_data(28)                  = inc_wm
i_trans_data(29)                  = nes
i_trans_data(30)                  = nncs
i_trans_data(31)                  = ipair
i_trans_data(32)                  = ibrem
i_trans_data(33)                  = isctn
i_trans_data(34)                  = isctnn
i_trans_data(35)                  = isctnA
i_trans_data(36)                  = iaenct

DO i = 1,nnu
  i_trans_data(40+0*nnu+i)        = itnu(i)
END DO

DO i = 1,nnu
  i_trans_data(40+1*nnu+i)        = nnugp(i)
END DO

d_trans_data(1)                   = tolnut
d_trans_data(2)                   = tolnuye
d_trans_data(3)                   = tolnupsi
d_trans_data(4)                   = tolpsimin
d_trans_data(5)                   = a_prec
d_trans_data(13)                  = dtnph_trans
d_trans_data(15)                  = dtnmhn_trans
d_trans_data(19)                  = rdtmax
d_trans_data(24)                  = unumn
d_trans_data(25)                  = unumx
d_trans_data(26)                  = edmpe
d_trans_data(27)                  = edmpa
d_trans_data(28)                  = rhoaefnp
d_trans_data(29)                  = roaencnu
d_trans_data(30)                  = rhonesmn
d_trans_data(31)                  = rhonesmx
d_trans_data(32)                  = rhonncsmn
d_trans_data(33)                  = rhonncsmx
d_trans_data(34)                  = rhopairemn
d_trans_data(35)                  = rhopairemx
d_trans_data(36)                  = rhopairtmn
d_trans_data(37)                  = rhopairtmx
d_trans_data(38)                  = rhosctnemn
d_trans_data(39)                  = rhosctnemx
d_trans_data(40)                  = rhosctntmn
d_trans_data(41)                  = rhosctntmx
d_trans_data(42)                  = rhosctnnemn
d_trans_data(43)                  = rhosctnnemx
d_trans_data(44)                  = rhosctnntmn
d_trans_data(45)                  = rhosctnntmx
d_trans_data(46)                  = rhosctnAemn
d_trans_data(47)                  = rhosctnAemx
d_trans_data(48)                  = rhosctnAtmn
d_trans_data(49)                  = rhosctnAtmx
d_trans_data(50)                  = rhopnu
d_trans_data(51)                  = roaenct

DO i = 11,40
  d_trans_data(50 +0*nnu+0*nez+i) = tcntrl(i)
END DO

d_trans_data(50 +0*nnu+0*nez+50)  = tcntrl(50)

DO i = 1,nnu
  d_trans_data(100+0*nnu+0*nez+i) = psimin(i)
END DO

DO i = 1,nnu
  d_trans_data(100+1*nnu+0*nez+i) = psipmin(i)
END DO

DO i = 1,nnu
  d_trans_data(100+2*nnu+0*nez+i) = stwt(i)
END DO

DO i = 1,nez
  d_trans_data(100+3*nnu+0*nez+i) = unui(i)
END DO

DO i = 1,nez
  d_trans_data(100+3*nnu+1*nez+i) = dunui(i)
END DO

DO i = 1,nez+1
  d_trans_data(100+3*nnu+2*nez+i) = unubi(i)
END DO


RETURN
END SUBROUTINE read_pack_transport_keys
