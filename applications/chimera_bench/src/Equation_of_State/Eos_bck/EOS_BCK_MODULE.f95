!-----------------------------------------------------------------------
!    Module:       eos_bck_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE eos_bck_module

USE kind_module, ONLY : double

SAVE



!-----------------------------------------------------------------------
!  Quantities used in the BCK EOS
!-----------------------------------------------------------------------
!  bad    : auxilary variable.
!
!  fmbad  : auxilary variable.
!
!  itsav  : the number of iterations done in SAHA.
!
!  jshel  : zone index ; for debugging.
!
!  af     : auxilary variable.
!
!  ahbck  : A[heavy nucleus].
!
!  b      : binding energy (heavy nucleus).
!
!  con    : auxilary variable.
!
!  dbck   : density (baryons/fm**3).
!
!  d0     : saturation density of dense phase, or density of nuclear matter
!   (baryons/fm**3).
!
!  d00    : zero value of saturation density of dense phase, or density of
!   Fe56 (baryons/fm**3).
!
!  dbdlu  : auxilary variable.
!
!  dedt - dE / dT  [approximate!].
!
!  delbet : ( ue - uneu - uhat ) [affinity] : not functional.
!
!  dtran  : transition density to nuclear matter.
!
!  ed     : energy in drip (including heavy translation).
!
!  ee     : energy in electrons.
!
!  egy0   : zero value of energy of Fe56 (MeV).
!
!  eh     : energy in heavy nuclei (or nuclear matter) w/o tranlation.
!
!  enu    : energy in neutrinos ; not calculated.
!
!  erad   : photon energy.
!
!  etot   : total energy ; no neutrinos ; includes radiation.
!
!  fa     : auxilary variable.
!
!  fm     : auxilary variable.
!
!  pd     : drip pressure (including heavy translation).
!
!  pe     : electron pressure.
!
!  ph     : heavy pressure (no translation).
!
!  phi    : auxilary variable.
!
!  pnu    : neutrino pressure ; not calculated.
!
!  prad   : photon pressure
!
!  ptotbck: total pressure (includes radiation).
!
!  rel    : 3 * pressure / kinetic energy (from 1(rel) to 2(nonrel)).
!
!  sd     : entropy in drip (includes heavy translation).
!
!  se     : entropy in electrons.
!
!  sh     : entropy in heavy nuclei (intrinsic only).
!
!  size0  : zero value of size of Fe56.
!
!  sneu   : neutrino entropy (not calculated).
!
!  stotbck: total entropy ; includes radiation.
!
!  tbck   : temperature (MeV).
!
!  therm  : auxilary variable.
!
!  theta  : density of dense phase / saturation density of dense phase.
!
!  theta0 : zero value of Fe56 density / saturation density of Fe56.
!
!  tsi    : auxilary variable.
!
!  upack  : packing fraction of dense phase.
!
!  u1pack : 1 - upack.
!
!  ue     : electron chemical potential (mu_e).
!
!  uhat   :  mu_n - mu_p.
!
!  uhtra  :  auxiliary variable
!
!  un     :  neutron chemical potential (mu_n).
!
!  uneu   :  neutrino chemical potential (mu_nu) ; not calculated.
!
!  xabck  :  alpha particle mass fraction ; = 0 for nuclear matter.
!
!  xhbck  :  heavy nuclei mass fraction ; = 1 for nuclear matter.
!
!  xmstar :  m* / m for heavy nuclei or nuclear matter.
!
!  xnbck  :  free neutron mass fraction ; = 0 for nuclear matter.
!
!  xpbck  :  free proton mass fraction  ; = 0 for nuclear matter.
!
!  ya     :  auxilary variable.
!
!  yeFe   :  zero value of Y_e^- - Y_e^+ for Fe56.
!
!  yebck  :  Y_e^- - Y_e^+.
!
!  yeplus :  Y_e^+.
!
!  yh     :  auxilary variable.
!
!  ynbck  :  auxilary variable.
!
!  ynu    :  Y_nu ; not used.
!
!  ypbck  :  auxilary variable.
!
!  zabck  :  Z / A for heavy nucleus ; = ye for nuclear matter. 
!
!  uea(j,id,it,iy)    : initial guess of ue, the electron chemical
!   potential, for cube corner (id,it,iy) of radial zone j when calling
!   the Cooperstein-BCK equation of state. When possible, these
!   quantities are saved from results of previous calls.
!
!  una(j,id,it,iy)    : initial guess of un, the neutron chemical potential,
!   for cube corner (id,it,iy) of radial zone j when calling the
!   Cooperstein-BCK equation of state. When possible, these quantities
!   are saved from results of previous calls.
!
!  uhata(j,id,it,iy)  : initial guess of un, the neutron chemical - proton
!   potential, for cube corner (id,it,iy) of radial zone j when calling
!   the Cooperstein-BCK equation of state. When possible, these quantities
!   are saved from results of previous calls.
!
!  thetaa(j,id,it,iy) : initial guess of theta, density of
!          dense phase / saturation density of dense phase
!   for cube corner (id,it,iy) of radial zone j when calling the
!   Cooperstein-BCK equation of state. When possible, these quantities
!   are saved from results of previous calls.
!
!  zaa(j,id,it,iy)    : initial guess of zabck, Z / A for heavy nucleus
!   (= ye for nuclear matter), for cube corner (id,it,iy) of radial zone
!   j when calling the Cooperstein-BCK equation of state. When possible,
!   these quantities are saved from results of previous calls.
!
!  xaa(j,id,it,iy)    : initial guess of xabck, alpha particle mass
!   fraction, for cube corner (id,it,iy) of radial zone j when calling
!   the Cooperstein-BCK equation of state. When possible, these quantities
!   are saved from results of previous calls.
!
!  dtrana(j,id,it,iy) : initial guess of theta, transition density to
!   nuclear matter, for cube corner (id,it,iy) of radial zone j when
!   calling the Cooperstein-BCK equation of state. When possible,
!   these quantities are saved from results of previous calls.
!-----------------------------------------------------------------------

LOGICAL                                        :: bad
LOGICAL                                        :: fmbad

INTEGER                                        :: itsav
INTEGER                                        :: jshel

REAL(KIND=double)                              :: af
REAL(KIND=double)                              :: ahbck
REAL(KIND=double)                              :: b
REAL(KIND=double)                              :: con
REAL(KIND=double)                              :: dbck
REAL(KIND=double)                              :: d0
REAL(KIND=double)                              :: d00
REAL(KIND=double)                              :: dbdlu
REAL(KIND=double)                              :: dedt
REAL(KIND=double)                              :: delbet
REAL(KIND=double)                              :: dtran  = 0.0d0
REAL(KIND=double)                              :: ed
REAL(KIND=double)                              :: ee
REAL(KIND=double)                              :: egy0
REAL(KIND=double)                              :: eh
REAL(KIND=double)                              :: enu
REAL(KIND=double)                              :: erad
REAL(KIND=double)                              :: etot
REAL(KIND=double)                              :: fa
REAL(KIND=double)                              :: fm
REAL(KIND=double)                              :: pd
REAL(KIND=double)                              :: pe
REAL(KIND=double)                              :: ph
REAL(KIND=double)                              :: phi
REAL(KIND=double)                              :: pnu
REAL(KIND=double)                              :: prad
REAL(KIND=double)                              :: ptotbck
REAL(KIND=double)                              :: rel
REAL(KIND=double)                              :: sd
REAL(KIND=double)                              :: se
REAL(KIND=double)                              :: sh
REAL(KIND=double)                              :: size0
REAL(KIND=double)                              :: sneu
REAL(KIND=double)                              :: stotbck
REAL(KIND=double)                              :: tbck
REAL(KIND=double)                              :: therm
REAL(KIND=double)                              :: theta
REAL(KIND=double)                              :: theta0
REAL(KIND=double)                              :: tsi
REAL(KIND=double)                              :: upack
REAL(KIND=double)                              :: u1pack
REAL(KIND=double)                              :: ue
REAL(KIND=double)                              :: uhat
REAL(KIND=double)                              :: uhtra
REAL(KIND=double)                              :: un
REAL(KIND=double)                              :: uneu
REAL(KIND=double)                              :: xabck
REAL(KIND=double)                              :: xhbck
REAL(KIND=double)                              :: xmstar
REAL(KIND=double)                              :: xnbck
REAL(KIND=double)                              :: xpbck
REAL(KIND=double)                              :: ya
REAL(KIND=double)                              :: yeFe
REAL(KIND=double)                              :: yebck
REAL(KIND=double)                              :: yeplus
REAL(KIND=double)                              :: yh
REAL(KIND=double)                              :: ynbck
REAL(KIND=double)                              :: ynu
REAL(KIND=double)                              :: ypbck
REAL(KIND=double)                              :: zabck

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:) :: uea
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:) :: una
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:) :: uhata
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:) :: thetaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:) :: zaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:) :: xaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:) :: dtrana

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: ueaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: unaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: uhataa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: thetaaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: zaaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: xaaa
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dtranaa

END module eos_bck_module
