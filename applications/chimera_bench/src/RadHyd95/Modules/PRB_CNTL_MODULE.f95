!-----------------------------------------------------------------------
!    Module:       prb_cntl_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE prb_cntl_module

USE kind_module
SAVE

!-----------------------------------------------------------------------
!
!             \\\\\ ZERO TRANSVERSE VELOCITY OPTION /////
!
!-----------------------------------------------------------------------
!  v_trans_0      : zero transverse velocity above shock switch
!
!     v_trans_0 = ye : transverse velocities 0 above shock
!     v_trans_0 = no : transverse velocities 0 above shock computed
!-----------------------------------------------------------------------

CHARACTER (len=2)                               :: v_trans_0

!-----------------------------------------------------------------------
!  Hydrodynamic Controls
!-----------------------------------------------------------------------
!  ihydro   : Hydrodynamics switch.
!
!     ihydro = 0   : hydrodynamics bypassed.
!     ihydro = 1   : hydrodynamics included.
!
!  irelhy   : Newtonian - GR hydrodynamics toggle.
!
!     irelhy = 0   : Newtonian hydrodynamics.
!     irelhy = 1   : PN hydrodynamics.
!     irelhy = 2   : general relativistic hydrodynamics.
!
!  ilapsehy : Newtonian - GR lapse toggle.
!
!     ilapsehy = 0 : lapse = 1 in the hydrodynamics if irelhy = 1.
!     ilapsehy = 1 : lapse computed in the hydrodynamics if irelhy = 1.
!
!  lagr     : Eulerian - Lagrangian hydrodynamics toggle.
!
!     lagr = ye    : Lagrangian hydrodynamics.
!     lagr = no    : Eulerian hydrodynamics.
!
!  rezn     : Regridding toggle.
!
!     rezn = ye    : regrid on startup or restart.
!     rezn = no    : no regriding on startup or restart.
!-----------------------------------------------------------------------

INTEGER                                         :: ihydro
INTEGER                                         :: irelhy
INTEGER                                         :: ilapsehy

CHARACTER (len=2)                               :: lagr
CHARACTER (len=2)                               :: rezn

!-----------------------------------------------------------------------
!  Gravitation Options
!-----------------------------------------------------------------------
!  i_grav     : gravitation key
!
!     i_grav = 1 : spherical symmetric gravity
!     i_grav = 2 : nonspherical components computed from by a Newtonian
!                   Poisson solver
!
!-----------------------------------------------------------------------

INTEGER                                         :: i_grav

!-----------------------------------------------------------------------
!  Transport controls
!-----------------------------------------------------------------------
!  inutrn   : Transport switch.
!
!     inutrn = 0 : neutrino transport bypassed.
!     inutrn = 1 : neutrino transport included.
!
!  ireltrns : Newtonian - GR transport toggle.
!
!     ireltrns = 0 : Newtonian neutrino transport.
!     ireltrns = 1 : general relativistic neutrino transport.
!
!  idiff    : Transport boundary conditions toggle.
!
!     idiff = 0 : dc(j,k,n) computed normally at each iteration.
!     idiff = 1 : dc(j,k,n) is the average of the current and preceding
!      timestep value.
!     idiff = 2 : dc(j,k,n) = 0 for j = jnumax, dc(j,k,n) computed normally
!      at each iteration for 
!      j ne jnumax.
!     idiff = 3 : dc(j,k,n) = 0 for j = jr_max, dc(j,k,n) as in idiff = 1
!      for j /= jr_max
!     idiff = 4 : dc(j,k,n) = 0 for all j.
!     idiff = 5 : dc(j,k,n) is computed only on the first iteration
!     idiff = 6 : dc(j,k,n) = 0 for j = jr_max, computed only on the first
!      iteration for j /= jr_max
!
!  jnumin   : inner radial zone for which neutrino transport is computed.
!  jnumax   : outer radial zone for which neutrino transport is computed.
!
!  itnu(n)  : transport temperature change switch.
!
!     itnu(n) = 0 : temperature change due to n-type neutrinos bypassed.
!     itnu(n) = 1 : temperature change due to n-type neutrinos included.
!
!  jexect   : transport temperature change for zone jexect switch.
!
!     jexect = 0 : computation of the temperature is normal.
!     jexect = j : if itnu = 0, computation of temperature is normal for
!      j = jexect.
!
!  iyenu    : transport electron fraction change switch.
!
!     iyenu = 0 : electron fraction change due to n-type neutrinos bypassed.
!     iyenu = 1 : electron fraction change due to n-type neutrinos included.
!
!  jexecy   : transport electron fraction change for zone jexect switch.
!
!     jexecy = 0 : computation of the temperature is normal.
!     jexecy = j : if iyenu = 0, computation of temperature is normal for
!      j = jexecy.
!-----------------------------------------------------------------------

INTEGER                                         :: inutrn
INTEGER                                         :: ireltrns
INTEGER                                         :: idiff
INTEGER                                         :: jnumin
INTEGER                                         :: jnumax
INTEGER, ALLOCATABLE, DIMENSION(:)              :: itnu
INTEGER                                         :: jexect
INTEGER                                         :: iyenu
INTEGER                                         :: jexecy

!-----------------------------------------------------------------------
!  Absorption and emission controls
!-----------------------------------------------------------------------
!  iaefnp   : emission and absorption of e-neutrinos and e-antineutrinos
!   on free neutrons and protons switch.
!
!     iaefnp = 0: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons turned off.
!     iaefnp = 1: emission and absorption of e-neutrinos and
!      e-antineutrinos on free neutrons and protons turned on.
!
!  i_aeps   : emission and aborption of e-neutrinos and e-antineutrinos
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
!
!  rhoaefnp : density above which emission and absorption of e-neutrinos
!   and e-antineutrinos on free neutrons and protons is turned off.
!
!  iaence   : emission and absorption of e-neutrinos on nuclei as given
!   by the FFN prescription  switch
!
!     iaence = 0: emission and absorption of e-neutrinos on nuclei as
!      given by the FFN prescription turned off.
!     iaence = 1: emission and absorption of e-neutrinos on nuclei as
!      given by the FFN prescription turned on.
!
!  edmpe    : difference between the mean excited state energy of the
!   daughter nucleus and that of the parent nucleus (MeV).
!
!  iaenca   : emission and absorption of e-antineutrinos on nuclei as
!   given by the FFN prescription  switch
!
!     iaenca = 0: emission and absorption of e-antineutrinos on nuclei as
!      given by the FFN prescription turned off.
!     iaenca = 1: emission and absorption of e-antineutrinos on nuclei as
!      given by the FFN prescription turned on.
!
!  edmpa    : difference between the mean excited state energy of the
!   daughter nucleus and that of  the parent nucleus (MeV).
!
!  iaencnu  : emission and absorption of e-neutrinos and e-antineutrinos
!   on nuclei as given by Haxton's prescription
!
!     iaencnu = 0: emission and absorption of e-neutrinos and
!      e-antineutrinos on nuclei as given by Haxton's prescription turned
!      off.
!     iaencnu = 1: emission and absorption of e-neutrinos and
!      e-antineutrinos on nuclei as given by Haxton's prescription turned
!      on.
!
!  roaencnu : density above which emission and absorption of e-neutrinos
!   and e-antineutrinos on nuclei as given by Haxton is turned off.
!
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

INTEGER                                         :: iaefnp
INTEGER                                         :: i_aeps
INTEGER                                         :: iaence
INTEGER                                         :: iaenca
INTEGER                                         :: iaencnu
INTEGER                                         :: iaenct

REAL(KIND=double)                               :: rhoaefnp
REAL(KIND=double)                               :: edmpe
REAL(KIND=double)                               :: edmpa
REAL(KIND=double)                               :: roaencnu
REAL(KIND=double)                               :: roaenct

!-----------------------------------------------------------------------
!  Weak magnetism correction for charged current interactions
!-----------------------------------------------------------------------
!  icc_wm : charged current weak magnetism switch.
!
!     icc_wm = 0: weak magnetism correction for charged current
!      interactions turned off.
!     icc_wm = 1: weak magnetism correction for charged current
!      interactions turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: icc_wm

!-----------------------------------------------------------------------
!  General scattering switch
!-----------------------------------------------------------------------
!  iscat : general scattering switch.
!
!     iscat = 0 : all scattering processes turned off.
!     iscat = 1 : any scattering process is allowed to proceed if its
!      key is turned on.                                           
!-----------------------------------------------------------------------

INTEGER                                         :: iscat

!-----------------------------------------------------------------------
!  Neutrino nucleon scattering controls
!-----------------------------------------------------------------------
!  in     : neutrino-neutron scattering switch
!
!     in = 0 : neutrino-neutron isoenergetic scattering turned off.
!     in = 1 : neutrino-neutron isoenergetic scattering turned on.
!
!  ip     : neutrino-proton scattering switch
!
!     ip = 0: neutrino-proton isoenergetic scattering turned off.
!     ip = 1: neutrino-proton isoenergetic scattering turned on.
!
!  ietann : final state blocking switch
!
!     ietann = 0 : final state blocking due to neutron or proton final
!      state occupancy in neutrino-neutron and neutrino-proton
!      isoenergetic scattering turned off.
!     ietann = 1 : final state blocking due to neutron or proton final
!      state occupancy in neutrino-neutron and neutrino-proton
!      isoenergetic scattering turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: in
INTEGER                                         :: ip
INTEGER                                         :: ietann

!-----------------------------------------------------------------------
!  Neutrino nuclei scattering controls
!-----------------------------------------------------------------------
!  ihe    : neutrino-helium isoenergetic scattering switch
!
!     ihe = 0: neutrino-helium isoenergetic scattering turned off.
!     ihe = 1: neutrino-helium isoenergetic scattering turned on.
!
!  iheavy : neutrino-heavy nucleus isoenergetic scattering switch
!
!     iheavy = 0: neutrino-heavy nucleus isoenergetic scattering turned off.
!     iheavy = 1: neutrino-heavy nucleus isoenergetic scattering turned on.
!
!  iicor  : ion-ion correlation switch.
!
!     iicor = 0: neutrino-heavy nucleus isoenergetic scattering
!      correction due to ion-ion correlation effects turned off.
!     iicor = 1: neutrino-heavy nucleus isoenergetic scattering
!      correction due to ion-ion correlation effects turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: ihe
INTEGER                                         :: iheavy
INTEGER                                         :: iicor

!-----------------------------------------------------------------------
!  Weak magnetism correction for neutral current interactions
!-----------------------------------------------------------------------
!  inc_wm : neutral current weak magnetism switch.
!
!     inc_wm = 0: weak magnetism correction for neutral current
!      interactions turned off.
!     inc_wm = 1: weak magnetism correction for neutral current
!      interactions turned on.
!-----------------------------------------------------------------------

INTEGER                                         :: inc_wm

!-----------------------------------------------------------------------
!  Neutrino-electron scattering controls
!-----------------------------------------------------------------------
!  nes      : neutrino-electron scattering switch.
!
!     nes = 0 : neutrino-electron scattering turned off.
!     nes = 1 : neutrino-electron scattering turned on.
!
!  rhonesmn : density below which n-neutrino-electron scattering is
!   turned off.
!
!  rhonesmx : density above which n-neutrino-electron scattering is
!   turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: nes

REAL(KIND=double)                               :: rhonesmn
REAL(KIND=double)                               :: rhonesmx

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering controls
!-----------------------------------------------------------------------
!  nncs      : neutrino-nucleus inelastic scattering switch.
!
!     nncs = 0 : neutrino-nucleus inelastic scattering turned off.
!     nncs = 1 : neutrino-nucleus inelastic scattering turned on.
!
!  rhonncsmn : density below which n-neutrino-electron scattering is
!   turned off.
!
!  rhonncsmx : density above which n-neutrino-electron scattering is
!   turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: nncs

REAL(KIND=double)                               :: rhonncsmn
REAL(KIND=double)                               :: rhonncsmx

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering controls
!-----------------------------------------------------------------------
!  isctn      : neutrino-nucleon inelastic scattering switch.
!
!     isctn = 0: neutrino-nucleon inelastic scattering turned off.
!     isctn = 1: neutrino-nucleon inelastic scattering turned on.
!
!  rhosctnemn : density below which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon  elastic scattering is turned off.
!
!  rhosctnemx : density above which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon elastic scattering is turned off.
!
!  rhosctntmn : density below which x-neutrino-neutrino-nucleon,
!   x-antineutrino-neutrino-nucleon  elastic scattering is turned off.
!
!  rhosctntmx : density above which x-neutrino-neutrino-nucleon,
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
!  isctnn      : neutrino-nucleon inelastic scattering switch.
!
!     isctnn = 0 : neutrino-nucleon inelastic scattering turned off.
!     isctnn = 1 : neutrino-nucleon inelastic scattering turned on.
!
!  rhosctnnemn : density below which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon  inelastic scattering is turned off.
!
!  rhosctnnemx : density above which e-neutrino-neutrino-nucleon,
!   e-antineutrino-neutrino-nucleon inelastic scattering is turned off.
!
!  rhosctnntmn : density below which x-neutrino-neutrino-nucleon,
!   x-antineutrino-neutrino-nucleon  inelastic scattering is turned off.
!
!  rhosctnntmx : density above which x-neutrino-neutrino-nucleon,
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
!  isctnA      : neutrino-nuclei inelastic scattering switch.
!
!     isctnA = 0 : neutrino-nucleus inelastic scattering turned off.
!     isctnA = 1 : neutrino-nucleus inelastic scattering turned on.
!
!  rhosctnAemn : density below which e-neutrino-neutrino-nucleus,
!   e-antineutrino-neutrino-nucleus  inelastic scattering is turned off.
!
!  rhosctnAemx : density above which e-neutrino-neutrino-nucleus,
!   e-antineutrino-neutrino-nucleus inelastic scattering is turned off.
!
!  rhosctnAtmn : density below which x-neutrino-neutrino-nucleus,
!   x-antineutrino-neutrino-nucleus  inelastic scattering is turned off.
!
!  rhosctnAtmx : density above which x-neutrino-neutrino-nucleus,
!   x-antineutrino-neutrino-nucleus scattering is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: isctnA

REAL(KIND=double)                               :: rhosctnAemn
REAL(KIND=double)                               :: rhosctnAemx
REAL(KIND=double)                               :: rhosctnAtmn
REAL(KIND=double)                               :: rhosctnAtmx

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation controls
!-----------------------------------------------------------------------
!  ipair      : electron-positron pair annihilation switch.
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
!   production by electron-positron pair annihilation is turned off..
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
!  Neutrino-nuclear pair annihilation controls
!-----------------------------------------------------------------------
!  ipairA      : muclear excitation pair annihilation switch.
!
!     ipairA = 0 : nuclear excitation pair annihilation process turned off.
!     ipairA = 1 : nuclear excitation pair annihilation process turned on.
!
!  rhopairAemn : density below which e-neutrino - e-antineutrino pair
!   production by nuclear deexcitation is turned off.
!
!  rhopairAemx : density above which e-neutrino - e-antineutrino pair
!   production by nuclear deexcitation is turned off.
!
!  rhopairAtmn : density below which x-neutrino - x-antineutrino pair
!   production by nuclear deexcitation is turned off..
!
!  rhopairAtmx : density above which x-neutrino - x-antineutrino pair
!   production by nuclear deexcitation is turned off.
!-----------------------------------------------------------------------

INTEGER                                         :: ipairA

REAL(KIND=double)                               :: rhopairAemn
REAL(KIND=double)                               :: rhopairAemx
REAL(KIND=double)                               :: rhopairAtmn
REAL(KIND=double)                               :: rhopairAtmx

!-----------------------------------------------------------------------
!  Bremsstrahlung pair annihilation controls
!-----------------------------------------------------------------------
!  ibrem      : bremsstrahlung pair annihilation switch.
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
!  Convection controls
!-----------------------------------------------------------------------
!  ilcnvct : mixing length convection switch.
!
!     ilcnvct = 0 : no convection.
!     ilcnvct = 1 : convection computed.
!
!  iemix   : energy mixing switch.
!
!     iemix = 0 : no energy mixing.
!     iemix = 1 : energy mixing computed.
!
!  iyemix  : ye mixing switch.
!
!     iyemix = 0 : no ye mixing.
!     iyemix = 1 : ye mixing computed.
!
!  ipsimix : psi0 mixing switch.
!
!     ipsimix = 0 : no psi0 mixing.
!     ipsimix = 1 : psi0 mixing computed.
!
!  ipcnvct : turbulent pressure switch.
!
!     ipcnvct = 0 : no convective, turbulent pressure.
!     ipcnvct = 1 : convective, turbulent pressure computed.
!
!  iscnvct : convective energy dissipation switch.
!
!     iscnvct = 0 : no convective energy dissipation.
!     iscnvct = 1 : convective energy dissipation computed.
!-----------------------------------------------------------------------

INTEGER                                         :: ilcnvct
INTEGER                                         :: iemix
INTEGER                                         :: iyemix
INTEGER                                         :: ipsimix
INTEGER                                         :: ipcnvct
INTEGER                                         :: iscnvct

!-----------------------------------------------------------------------
!  NSE flashing and deflashing controls
!-----------------------------------------------------------------------
!  tnse  : temperature at which material is flashed to nse (K)
!
!  tdnse : temperature at which material is deflashed from nse (K)
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: tnse
REAL(KIND=double)                               :: tdnse

!-----------------------------------------------------------------------
!  Explosion simulation control
!-----------------------------------------------------------------------
!  i_bomb : bomb switch
!
!     i_bomb = 0 : no bomb
!     i_bomb = 1 : bomb
!-----------------------------------------------------------------------

INTEGER                                         :: i_bomb

END module prb_cntl_module
