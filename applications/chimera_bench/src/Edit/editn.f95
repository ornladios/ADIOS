SUBROUTINE editn( n, jr_min, jr_max, ij_ray, ik_ray )
!-----------------------------------------------------------------------
!
!    File:         editn
!    Module:       editn
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/21/97
!
!    Purpose:
!      To edit neutrino data.
!
!    Subprograms called:
!  date_and_time_print : gets date and time
!  enutcp              : computes thermodynamic fitting parameters for neutrino distributions
!  abem_cal            : computes inverse mean free paths for neutrino absorption and emission on mucleons
!  abemnc              : computes inverse mean free paths for neutrino absorption on nuclei (FFN)
!  abemhnc             : computes inverse mean free paths for neutrino absorption on nuclei (Haxton)
!  scatical            : computes isoenergetic scattering on nucleons (when selected) and nuclei
!  scatical_ed         : computes isoenergetic scattering on nucleons inverse mean free paths
!  sctednes            : computes neutrino electron scattering inverse mean free paths
!  sctednns            : computes non-isoenergetic scattering on nucleons inverse mean free paths
!  sctednAs            : computes non-isoenergetic scattering on nuclei (FM) inverse mean free paths
!  scatiicr            : computes ion-ion correlations
!  paired              : computes electron-positron pair annihilation inverse mean free paths
!
!    Input arguments:
!
!  n          : neutrino type
!  jr_min     : inner radial zone of region for which configuration
!                edit is to be made.
!  jr_max     : outer radial zone of region for which configuration
!                edit is to be made.
!  ij_ray     : j-index of a radial ray
!  ik_ray     : k-index of a radial ray
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!
!  prnttest   : true  - test to see if printing criteria is satisfied.
!               false - bypass printing criteria test.
!  iprint     : 0    - do not print to print file.
!               ne 0 - print to print file.
!  nprint     : unit number of print file.
!  iplot      : 0    - do not print to plot file.
!               ne 0 - print to plot file.
!  nplot      : unit number of plot file.
!  nedn(i)    : editc counter for data set i.
!  intedn(i)  : number of cycles between edits of data set i.
!  idxedn(i)  : edit jr_min, jr_max, and every idxedc(i) radial zone
!                between them for data set i.
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  abem_module, edit_module, eos_snc_x_module, mdl_cnfg_module,
!  nu_dist_module, nu_energy_grid_module, scat_i_module
!
!-----------------------------------------------------------------------
!
!                  n-type neutrino transport data
!-----------------------------------------------------------------------
! enu(k)      'w'            Energy of n-neutrinos in group k (MeV).
! psi0eq(k)   'psieq'        Beta-equilibrium n-neutrino occupation number, zone j, group k.
! psi0(j,k,n) 'psi0'         Zero moment of the n-neutrino occupation number, zone j, group k.
! psi1(j,k,n) 'psi1'         First moment of the n-neutrino occupation number, zone j, group k
! psiteq(k)   'psiteq'       Thermal equilibrium n-neutrino occupation  number, i.e., the Fermi-Dirac n-neutrino
!                             distribution having the same total number and energy as the actual n-neutrino
!                             distribution, zone j, group k.
! reanpin(k) 'dndt ec p in'  Rate of n-neutrino production by electron capture on free protons per unit volume
!                             per unit energy (cm^{-3} s^{-1} MeV^{-1})
! reanpin(k) 'dndt ec p nt'  Net rate of n-neutrino production by electron capture on free protons and the
!                             inverse reaction per unit volume per unit energy (cm^{-3} s^{-1} MeV^{-1}).
! reancin(k) 'dndt ec a in'  Rate of n-neutrino production by electron capture on nuclei per unit volume per unit
!                             energy (FFN/B) (cm^{-3} s^{-1} MeV^{-1}).
! reancnet(k)'dndt ec a nt'  Net rate of n-neutrino production by electron capture on nuclei and the inverse
!                             reaction per unit volume per unit energy (FFN/B) (cm^{-3} s^{-1} MeV^{-1}).
! reahncin(k)'dndt ec a in'  Rate of n-neutrino production by electron capture on nuclei per unit volume per unit
!                             energy (Haxton) (cm^{-3} s^{-1} MeV^{-1}).
! reahncnet(k)'dndt ec a nt' Net rate of n-neutrino production by electron capture on nuclei and the inverse
!                             reaction per unit volume per unit energy (Haxton) (cm^{-3} s^{-1} MeV^{-1}).
! rpnesin(k) 'dndt nes in'   Rate of n-neutrino production by neutrino-electron "in" scattering per unit volume
!                             per unit energy (cm^{-3} s^{-1} MeV^{-1}).
! rpnesnet(k)'dndt nes net'  Net rate of n-neutrino production by neutrino-electron scattering per unit
!                             volume per unit energy ( cm^{-3} s^{-1} MeV^{-1}).
! rpain(k)  'dndt p-a in'    Rate of n-neutrino production by electron-positron pair annihilation per unit volume
!                             per unit energy (cm^{-3} s^{-1} MeV^{-1}).
! rpanet(k) 'dndt p-a net'   Net rate of n-neutrino production by electron-positron pair annihilation and
!                             the inverse reaction per unit volume per unit energy (cm^{-3} s^{-1} MeV^{-1}).
! dndt_v(j,k,n) 'dndt v'     Rate of n-neutrino "production" due to energy advection (cm^{-3} s^{-1} MeV^{-1})
! rlocal(k)   'lcl tot'      Rate of e-neutrino "production" due to all local processes (cm^{-3} s^{-1} MeV^{-1})
! rtr(k)      'dndt tr'      Rate of e-neutrino "production" due to transport (cm^{-3} s^{-1} MeV^{-1})
! ymefrnp(k)  'mfp-1 npe'    n-neutrino free nucleon emission inverse mean free path (cm^{-1}).
! ymenc(k)    'mfp-1 nce'    n-neutrino nuclei emission inverse mean free path (FFN/B) (cm^{-1}).
! ymehnc(k)   'mfp-1 nca'    n-neutrino nuclei emission inverse mean free path (Haxton) (cm^{-1}).
! taueanp(k)  'mfp-1 npa'    n-neutrino free nucleon absorption inverse mean free path (cm^{-1}).
! taueanc(k)  'mfp-1 nca'    n-neutrino nuclei absorption inverse mean
!                             free path (FFN/B) (cm^{-1}).
! taueahnc(k) 'mfp-1 nca'    n-neutrino nuclei absorption inverse mean free path (Haxton) (cm^{-1}).
! ymdnps(k)   'mfp-1 nns'    n-neutrino-proton scattering isoenergetic inverse mean free path (cm^{-1}).
! ymdnns(k)   'mfp-1 nns'    n-neutrino-neutron scattering isoenergetic inverse mean free path (cm^{-1}).
! ymdnhs(k)   'mfp-1 nhs'    n-neutrino-helium scattering isoenergetic inverse mean free path (cm^{-1}).
! ymdnas(k)   'mfp-1 nas'    n-neutrino-nuclei scattering isoenergetic inverse mean free path (cm^{-1}).
! ymdnes(k)   'mfp-1 nes'    n-neutrino-electron scattering inverse mean free path (cm^{-1}).
! ymdnts(k)   'mfp-1 nts'    n-neutrino pair annihilation and production inverse mean free path (cm^{-1}).
! ymdnbs(k)   'mfp-1 nbs'    n-neutrino pair bremsstrahlung annihilation and production inverse mean free path (cm^{-1}).
! ymdnn(k)    'mfp-1 nN '    n-neutrino-nucleon elastic scattering inverse mean free path (cm^{-1}).
! ymdnnn(k)   'mfp-1 nNN'    n-neutrino-nucleon inelastic scattering inverse mean free path (cm^{-1}).
! ymdnAn(k)   'mfp-1 nA'     n-neutrino-nucleus inelastic scattering inverse mean free path (cm^{-1}).
! tautrns(k)  'tautrnsp'     Total n-neutrino transport optical depth from surface.
! tauopt(k)   'tauopt'       Total n-neutrino outbeam optical depth from surface.
! rnes(k)     'rnes'         Rate of n-neutrino-electron scattering for n-neutrinos of energy group k (/sec).
! qnes(k)	  'qnes'         Fraction of n-neutrino energy transferred to electrons by nes per scattering event.
! rnncs(k)    'rnncs'        Rate of n-neutrino-nuclei inelastic scattering for n-neutrinos of energy group k (/sec).
! qnncs(k)	  'qnncs'        Fraction of n-neutrino energy transferred to electrons by nncs per scattering event.
! rns(k)      'rnns'         Rate of n-neutrino-nucleon elastic scattering for n-neutrinos of energy group k (/sec).
! qns(k)      'qnns'         Fraction of n-neutrino energy transferred to electrons by nns  per scattering event.
!                             to electrons by nncs per scattering event.
! rnns(k)     'rnnns'        Rate of n-neutrino-nucleon inelastic scattering for n-neutrinos of energy group k (/sec).
! qnns(k)     'qnnns'        Fraction of n-neutrino energy transferred to electrons by nnns per scattering event.
! rnAs(k)     'rnAs'         Rate of n-neutrino-nucleus inelastic scattering for n-neutrinos of energy group k (/sec).
! qnAs(k)     'qnAs'         Fraction of n-neutrino energy transferred to electrons by nnns per scattering event.
! fluxn(k)    'd(nlum)/de'   Number luminosity per unit energy of n-neutrinos of group k across zone j
!                             (s^{-1} MeV^{-1}).
! dnuradjn(k) 'enu cross'    Total number of e-neutrinos per unit energy that have crossed zone j (MeV^{-1})
! vdrift(k)   'drift v'      Drift velocity of e-neutrinos of group k (cm/sec).
! strsnujn(k) 'estress'      n-neutrino-matter force (dynes g^{-1} MeV^{-1})
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY : nx, nez, nnu
USE numerical_module, ONLY : zero, half, one, ecoef, epsilon, frpi
USE physcnst_module, ONLY : cvel, pi, h, rmu, msolar, ergfoe, kmev

USE abem_module, ONLY : emis, absor
USE e_advct_module, ONLY : flxf, Ef, unucrv, dndt_v
USE edit_module, ONLY : head, prnttest, nedn, neden, niedn, intedn, idxedn, &
& nprint, nlog
USE eos_snc_x_module, ONLY : aesv
USE mdl_cnfg_module, ONLY : rho, t, ye, dmrst, rstmss, r, dr
USE nu_dist_module, ONLY : unu, dunu, stwt, psi0, psi1, dnurad, strsnu_x, &
& unulum, unurad, unucr, unucrea, unucrnis, unucrpa, unucrba, unucrnns, &
& unucrt, nnurad, nnucr, j_sphere, r_sphere, d_sphere, t_sphere, m_sphere, &
& tmfp_j, unukrad, fluxlm, dc
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE scat_i_module, ONLY : rmdnhes0, rmdnhs0

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min        ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)              :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=24)               :: type
CHARACTER (len=26)               :: under
CHARACTER (len=2), DIMENSION(4)  :: nutype = (/'e ','eb','t ','tb'/)
CHARACTER (len=10)               :: var_name

LOGICAL                          :: first = .true.
LOGICAL                          :: itaut         ! emission surface logical key

INTEGER, DIMENSION(3)            :: kmin = (/1,11,21/)  ! segments of energy binning index
INTEGER, DIMENSION(3)            :: kmax = (/10,20,30/) ! segments of energy binning index
INTEGER, DIMENSION(40)           :: nedens        ! temporary storage for neden

INTEGER                          :: i             ! do index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! do index
INTEGER                          :: l             ! do index
INTEGER                          :: jv            ! do index
INTEGER                          :: jj            ! radial zone index of emission surface
INTEGER                          :: icnt          ! energy group counter
INTEGER                          :: n_time        ! used for data-and-time
INTEGER                          :: jedit         ! idxedn(n)
INTEGER                          :: jp1           ! j + 1
INTEGER                          :: istat         ! allocation status

REAL(KIND=double)                :: cd3           ! c/3
REAL(KIND=double)                :: tthird        ! 2/3
REAL(KIND=double)                :: fcoef         ! coefficient for computing neutrino number flux ( cm/s /Mev**3 cm**3 )
REAL(KIND=double)                :: nstate        ! coefficient for computing neutrino number density ( /Mev**3 cm**3 )
REAL(KIND=double)                :: bedot         ! coefficient for computing heating rate (MeV/MeV s)
REAL(KIND=double)                :: tfluxc        ! coefficient for computing neutrino luminosity

REAL(KIND=double)                :: drj           ! dr(j)
REAL(KIND=double)                :: aream         ! area at j-1
REAL(KIND=double)                :: areaj         ! area at j
REAL(KIND=double)                :: enut          ! temperature of an neutrino distribution with same 0, 1, 2 moments
REAL(KIND=double)                :: enucmp        ! chemical pot of an neutrino distribution with same 0, 1, 2 moments

REAL(KIND=double)                :: w2            ! neutrino energy squared
REAL(KIND=double)                :: cflxnm        ! coefficient for computing rate of e-neutrino "production" by transport
REAL(KIND=double)                :: cfluxn        ! coefficient for computing rate of e-neutrino "production" by transport
REAL(KIND=double)                :: expa          ! argument of an exponential
REAL(KIND=double)                :: ekt           ! neutrino energy/kT
REAL(KIND=double)                :: xhe           ! mass fraction of helium

REAL(KIND=double)                :: emisnp        ! emission invmfp for electron capture on free protons (cm^-1)
REAL(KIND=double)                :: absornp       ! absorption invmfp for electron capture on free protons (cm^-1)
REAL(KIND=double)                :: emisnc        ! emission invmfp for electron capture on nuclei (FFN) (cm^-1)
REAL(KIND=double)                :: absornc       ! absorption invmfp for electron capture on nuclei (FFN) (cm^-1)
REAL(KIND=double)                :: emishnc       ! emission invmfp for electron capture on nuclei (Haxton) (cm^-1)
REAL(KIND=double)                :: absorhnc      ! absorption invmfp for electron capture on nuclei (Haxton) (cm^-1)

REAL(KIND=double)                :: fexp          ! external exponential function

REAL(KIND=double)                :: rmdnns        ! mfp^-1 for neutrino-neutron scattering
REAL(KIND=double)                :: rmdnps        ! mfp^-1 for neutrino-proton scattering
REAL(KIND=double)                :: rmdnbns       ! mfp^-1 for antineutrino-neutron scattering
REAL(KIND=double)                :: rmdnbps       ! mfp^-1 for antineutrino-proton scattering
REAL(KIND=double)                :: rmdnhes       ! mfp^-1 for neutrino-helium scattering
REAL(KIND=double)                :: rmdnhs        ! mfp^-1 for neutrino-heavy nucleus scattering
REAL(KIND=double)                :: coh           ! sum of above for neutrinos
REAL(KIND=double)                :: cohb          ! sum of above for antineutrinos
REAL(KIND=double)                :: rmdnns_ed     ! mfp^-1 for neutrino-neutron scattering (for editing purposes)
REAL(KIND=double)                :: rmdnps_ed     ! mfp^-1 for neutrino-proton scattering (for editing purposes)

REAL(KIND=double)                :: xi_p_wm       ! weak magnetism correction for neutrino-proton scattering
REAL(KIND=double)                :: xi_n_wm       ! weak magnetism correction for neutrino-neutron scattering
REAL(KIND=double)                :: xib_p_wm      ! weak magnetism correction for antineutrino-proton scattering
REAL(KIND=double)                :: xib_n_wm      ! weak magnetism correction for antineutrino-neutron scattering

REAL(KIND=double)                :: rmdnes        ! neutrino inverse mean free path
REAL(KIND=double)                :: rmdnes0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                :: rmdnes1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                :: a0w           ! coefficient for computing change in psi0
REAL(KIND=double)                :: b0w           ! coefficient for computing change in psi0
REAL(KIND=double)                :: c0w           ! coefficient for computing change in psi0
REAL(KIND=double)                :: a1w           ! coefficient for computing change in psi1
REAL(KIND=double)                :: b1w           ! coefficient for computing change in psi1
REAL(KIND=double)                :: c1w           ! coefficient for computing change in psi1

REAL(KIND=double)                :: arscte        ! coefficient for computing NES scattering rate
REAL(KIND=double)                :: brscte        ! coefficient for computing NES scattering rate
REAL(KIND=double)                :: crscte        ! coefficient for computing NES scattering rate
REAL(KIND=double)                :: arsctu        ! coefficient for computing NES upscattering rate
REAL(KIND=double)                :: brsctu        ! coefficient for computing NES upscattering rate
REAL(KIND=double)                :: crsctu        ! coefficient for computing NES upscattering rate
REAL(KIND=double)                :: arscti        ! coefficient for computing NES isoscattering rate
REAL(KIND=double)                :: brscti        ! coefficient for computing NES isoscattering rate
REAL(KIND=double)                :: crscti        ! coefficient for computing NES isoscattering rate
REAL(KIND=double)                :: arsctd        ! coefficient for computing NES downscattering rate
REAL(KIND=double)                :: brsctd        ! coefficient for computing NES downscattering rate
REAL(KIND=double)                :: crsctd        ! coefficient for computing NES downscattering rate
REAL(KIND=double)                :: aesctu        ! coefficient for computing NES upscattering energy
REAL(KIND=double)                :: besctu        ! coefficient for computing NES upscattering energy
REAL(KIND=double)                :: cesctu        ! coefficient for computing NES upscattering energy
REAL(KIND=double)                :: aescti        ! coefficient for computing NES isoscattering energy
REAL(KIND=double)                :: bescti        ! coefficient for computing NES isoscattering energy
REAL(KIND=double)                :: cescti        ! coefficient for computing NES isoscattering energy
REAL(KIND=double)                :: aesctd        ! coefficient for computing NES downscattering energy
REAL(KIND=double)                :: besctd        ! coefficient for computing NES downscattering energy
REAL(KIND=double)                :: cesctd        ! coefficient for computing NES downscattering energy
REAL(KIND=double)                :: arnes         ! coefficient for computing NES rate
REAL(KIND=double)                :: brnes         ! coefficient for computing NES rate
REAL(KIND=double)                :: crnes         ! coefficient for computing NES rate
REAL(KIND=double)                :: aenes         ! coefficient for computing NES fractional eneergy transfer
REAL(KIND=double)                :: benes         ! coefficient for computing NES fractional eneergy transfer
REAL(KIND=double)                :: cenes         ! coefficient for computing NES fractional eneergy transfer

REAL(KIND=double)                :: dpsi0nesin    ! coefficient for computing rate of NES "in" scattering
REAL(KIND=double)                :: dpsi0nesnet   ! coefficient for computing rate of NES net scattering
REAL(KIND=double)                :: einitial      ! coefficient for computing n-neutrino energy transferred to electrons
REAL(KIND=double)                :: efinal        ! coefficient for computing n-neutrino energy transferred to electrons

REAL(KIND=double)                :: rmdnncs       ! neutrino inverse mean free path due to inelastic neutrino-nucleus scat (Haxton)
REAL(KIND=double)                :: rmdnncs0      ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                :: rmdnncs1      ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                :: arnncs        ! coefficient for computing NNCS (Haxton) rate
REAL(KIND=double)                :: brnncs        ! coefficient for computing NNCS (Haxton) rate
REAL(KIND=double)                :: crnncs        ! coefficient for computing NNCS (Haxton) rate
REAL(KIND=double)                :: aenncs        ! coefficient for computing NNCS (Haxton) fractional eneergy transfer
REAL(KIND=double)                :: benncs        ! coefficient for computing NNCS (Haxton) fractional eneergy transfer
REAL(KIND=double)                :: cenncs        ! coefficient for computing NNCS (Haxton) fractional eneergy transfer

REAL(KIND=double)                :: dpsi0nncsin   ! coefficient for computing rate of NNCS "in" scattering
REAL(KIND=double)                :: dpsi0nncsnet  ! coefficient for computing rate of NNCS net scattering

REAL(KIND=double)                :: rmdnn         ! neutrino inverse mean free path due to ielastic neutrino-nucleon scat (Haxton)
REAL(KIND=double)                :: rmdnn0        ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                :: rmdnn1        ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                :: dpsi0nnsin    ! coefficient for computing rate of NNS "in" scattering
REAL(KIND=double)                :: dpsi0nnsnet   ! coefficient for computing rate of NNS net scattering

REAL(KIND=double)                :: rmdnnn        ! neutrino inverse mean free path due to ielastic neutrino-nucleon scat (Haxton)
REAL(KIND=double)                :: rmdnnn0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                :: rmdnnn1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                :: dpsi0nnnsin   ! coefficient for computing rate of NNNS "in" scattering
REAL(KIND=double)                :: dpsi0nnnsnet  ! coefficient for computing rate of NNNS net scattering

REAL(KIND=double)                :: rmdnnA        ! neutrino inverse mean free path due to ielastic neutrino-nuclei scat (Haxton)
REAL(KIND=double)                :: rmdnnA0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                :: rmdnnA1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                :: dpsi0nnAsin   ! coefficient for computing rate of NAS "in" scattering
REAL(KIND=double)                :: dpsi0nnAsnet  ! coefficient for computing rate of NAS net scattering

REAL(KIND=double)                :: rmdnts        ! neutrino inverse mean free path due to electron-positron pair annihilation
REAL(KIND=double)                :: rmdnts0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                :: rmdnts1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                :: artpe         ! coefficient for computing pair production rate
REAL(KIND=double)                :: brtpe         ! coefficient for computing pair production rate
REAL(KIND=double)                :: crtpe         ! coefficient for computing pair production rate
REAL(KIND=double)                :: artae         ! coefficient for computing pair annihilation rate
REAL(KIND=double)                :: brtae         ! coefficient for computing pair annihilation rate
REAL(KIND=double)                :: crtae         ! coefficient for computing pair annihilation rate

REAL(KIND=double)                :: dpsi0pain     ! coefficient for computing rate of pair "in" production
REAL(KIND=double)                :: dpsi0panet    ! coefficient for computing rate of pair net production

REAL(KIND=double)                :: rmdnbs        ! neutrino inverse mean free path due to nucleon bremsstrahlung
REAL(KIND=double)                :: rmdnbs0       ! zero moment of the neutrino inverse mean free path
REAL(KIND=double)                :: rmdnbs1       ! first moment of the neutrino inverse mean free path

REAL(KIND=double)                :: dpsi0bain     ! coefficient for computing rate of bremsstrahlung "in" production
REAL(KIND=double)                :: dpsi0banet    ! coefficient for computing rate of bremsstrahlung net production

REAL(KIND=double)                :: wenulr        ! core n-neutrino luminosity (foes/sec)
REAL(KIND=double)                :: wenul         ! the cumulative energy emitted from the core in n-type neutrinos (foes)
REAL(KIND=double)                :: wenut         ! the total energy of n-type neutrinos currently residing in the core (foes)
REAL(KIND=double)                :: wenuta        ! total energy transferred to n-type neutrinos by emission and absorption (foes)
REAL(KIND=double)                :: wenuts        ! the total energy transferred to n-type neutrinos by NES (foes)
REAL(KIND=double)                :: wenutp        ! total energy transferred to n-type neutrinos by e-p pair annihilation (foes)
REAL(KIND=double)                :: wenutb        ! total energy transferred to n-type neutrinos by bremsstrahlung (foes)
REAL(KIND=double)                :: wenutnn       ! total energy transferred to n-type neutrinos by NNNS (foes)
REAL(KIND=double)                :: wenutt        ! total energy transferred to n-neutrinos by microphysics (foes)
REAL(KIND=double)                :: wenutv        ! total energy transferred to n-type neutrinos by advection (foes)
REAL(KIND=double)                :: wnnurad       ! the cumulative number of n-type neutrinos emitted by the core
REAL(KIND=double)                :: wnnucr        ! the net number of n-type neutrinos currently residing in the core

REAL(KIND=double)                :: elecn         ! total electron number
REAL(KIND=double)                :: totlpn        ! total lepton number

REAL(KIND=double)                :: taut          ! total Planck averaged optical depth to zone j
REAL(KIND=double)                :: tautp         ! total Planck averaged optical depth to zone j + 1
REAL(KIND=double)                :: ymdavt        ! parameter for calculating taut
REAL(KIND=double)                :: flux          ! parameter for calculating taut
REAL(KIND=double)                :: ymdave        ! Planck averaged optical depth of zone j
REAL(KIND=double)                :: remsf         ! radius of emission surface (cm)
REAL(KIND=double)                :: demsf         ! density of emission surface (g cm^-3)
REAL(KIND=double)                :: gemst         ! enclosed mass of emission surface (g)
REAL(KIND=double)                :: gemsf         ! enclosed mass of emission surface (msolar)
REAL(KIND=double)                :: temsf         ! temperature of emission surface (MeV)

REAL(KIND=double)                :: ymdavat       ! parameter for calculating thermal emission surface
REAL(KIND=double)                :: ymdavst       ! parameter for calculating thermal emission surface
REAL(KIND=double)                :: ymdava        ! Planck averaged optical depth of zone j
REAL(KIND=double)                :: ymdavs        ! Planck averaged optical depth of zone j
REAL(KIND=double)                :: tautherm      ! Planck averaged optical depth of zone j

REAL(KIND=double)                :: dum1          ! dummy variable
REAL(KIND=double)                :: dum2          ! dummy variable
REAL(KIND=double)                :: dum3          ! dummy variable

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: psi0l
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: psi1l
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: psiteq
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: psi10

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: tautrns
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: tauopt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: taueanp
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: taueanc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: taueahnc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxn

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: tauajk
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: taujk
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: tausjk

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dnuradjn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: strsnujn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: edotnp

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: psi0eq
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: reanpin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: reanpnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: reancin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: reancnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: reahncin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: reahncnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rlocal
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rtr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: tflux
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: vdrift
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ciicra

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnesin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnesnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: qnes
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rnes

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnncsin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnncsnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: qnncs
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rnncs

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnnsin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnnsnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: qns
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rns

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnnnsin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnnnsnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: qnns
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rnns

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnnAsin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpnnAsnet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: qnAs
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rnAs

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: wk_mag_p
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: wk_mag_n

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpain
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rpanet
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rbain
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rbanet

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymefrnp
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymenc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymehnc
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnps
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnns
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnhs
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnas
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnes
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnncs
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnnn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnAn
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnts
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnbs
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnps_ed
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ymdnns_ed

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: ecoefp
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: uenulf

EXTERNAL fexp

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

    1 FORMAT (1x)
    3 FORMAT (1x,a128)
    5 FORMAT (1x,128('*'))
    7 FORMAT (' n has a value other than 1 - 4 in editn')
  101 FORMAT (23x,a24)
  103 FORMAT (22x,a26/)
  105 FORMAT ('0',120('.')/)
  201 FORMAT (1x,i4,'  w       ',10(1pe11.3))
  203 FORMAT (4x,' dw       ',10(1pe11.3))
  205 FORMAT (4x,' psieq    ',10(1pe11.3))
  207 FORMAT (4x,' psi0     ',10(1pe11.3))
  209 FORMAT (4x,'logpsi0   ',10(1pe11.3))
  211 FORMAT (4x,' psi0p    ',10(1pe11.3))
  213 FORMAT (4x,' psi1     ',10(1pe11.3))
  215 FORMAT (4x,'logpsi1   ',10(1pe11.3))
  217 FORMAT (4x,'fluxfactor',10(1pe11.3))
  219 FORMAT (4x,'Edd factor',10(1pe11.3))
  221 FORMAT (2x,'   psiteq   ',10(1pe11.3))
  223 FORMAT (2x,'dndt ec p in',10(1pe11.3))
  225 FORMAT (1x,'dndt ec p net',10(1pe11.3))
  227 FORMAT (2x,'dndt ec a in',10(1pe11.3))
  229 FORMAT (1x,'dndt ec a net',10(1pe11.3))
  231 FORMAT (1x,'dndt ec ah in',10(1pe11.3))
  233 FORMAT (1x,'dndt ec ah nt',10(1pe11.3))
  235 FORMAT (2x,'dndt nes in ',10(1pe11.3))
  237 FORMAT (2x,'dndt nes net',10(1pe11.3))
  239 FORMAT (1x,'dndt nncs in ',10(1pe11.3))
  241 FORMAT (1x,'dndt nncs net',10(1pe11.3))
  243 FORMAT (2x,'dndt p-a in ',10(1pe11.3))
  245 FORMAT (2x,'dndt p-a net',10(1pe11.3))
  247 FORMAT (2x,'dndt b-a in ',10(1pe11.3))
  249 FORMAT (2x,'dndt b-a net',10(1pe11.3))
  251 FORMAT (2x,'dndt nNs in ',10(1pe11.3))
  253 FORMAT (2x,'dndt nNs net',10(1pe11.3))
  255 FORMAT (1x,'dndt nNNs in ',10(1pe11.3))
  257 FORMAT (1x,'dndt nNNs net',10(1pe11.3))
  259 FORMAT (2x,'dndt v      ',10(1pe11.3))
  261 FORMAT (2x,'dndt lcl tot',10(1pe11.3))
  263 FORMAT (2x,'dndt tr     ',10(1pe11.3))
  265 FORMAT (2x,'  edot np   ',10(1pe11.3))
  267 FORMAT (2x,'   rnes     ',10(1pe11.3))
  269 FORMAT (2x,'   qnes     ',10(1pe11.3))
  271 FORMAT (2x,'  rnncs     ',10(1pe11.3))
  273 FORMAT (2x,'  qnncs     ',10(1pe11.3))
  275 FORMAT (2x,'  rnns      ',10(1pe11.3))
  277 FORMAT (2x,'  qnns      ',10(1pe11.3))
  279 FORMAT (2x,'  rnnns     ',10(1pe11.3))
  281 FORMAT (2x,'  qnnns     ',10(1pe11.3))
  283 FORMAT (2x,'  rnAs      ',10(1pe11.3))
  285 FORMAT (2x,'  qnAs      ',10(1pe11.3))
  287 FORMAT (2x,'mfp-1 npe   ',10(1pe11.3))
  289 FORMAT (2x,'mfp-1 nce   ',10(1pe11.3))
  291 FORMAT (1x,'mfp-1 hnce   ',10(1pe11.3))
  293 FORMAT (2x,'mfp-1 npa   ',10(1pe11.3))
  295 FORMAT (2x,'mfp-1 nca   ',10(1pe11.3))
  297 FORMAT (1x,'mfp-1 hnca   ',10(1pe11.3))
  299 FORMAT (2x,'mfp-1 nps   ',10(1pe11.3))
  301 FORMAT (2x,'mfp-1 nps ed',10(1pe11.3))
  303 FORMAT (2x,'mfp-1 nns   ',10(1pe11.3))
  305 FORMAT (2x,'mfp-1 nns ed',10(1pe11.3))
  307 FORMAT (2x,'mfp-1 nhs   ',10(1pe11.3))
  309 FORMAT (2x,'mfp-1 nas   ',10(1pe11.3))
  311 FORMAT (2x,'mfp-1 nes   ',10(1pe11.3))
  313 FORMAT (1x,'mfp-1 nncs   ',10(1pe11.3))
  315 FORMAT (2x,'mfp-1 nts   ',10(1pe11.3))
  317 FORMAT (2x,'mfp-1 nbs   ',10(1pe11.3))
  319 FORMAT (2x,'mfp-1 nN    ',10(1pe11.3))
  321 FORMAT (2x,'mfp-1 nNN   ',10(1pe11.3))
  323 FORMAT (2x,'mfp-1 nA    ',10(1pe11.3))
  325 FORMAT (3x,'wk_mag_p   ',10(1pe11.3))
  327 FORMAT (3x,'wk_mag_n   ',10(1pe11.3))
  329 FORMAT (3x,'tautrnsp   ',10(1pe11.3))
  331 FORMAT (5x,'tauopt   ',10(1pe11.3))
  333 FORMAT (4x,'ion-ion   ',10(1pe11.3))
  335 FORMAT (1x,'d(nlum)/de   ',10(1pe11.3))
  337 FORMAT (2x,'enu cross   ',10(1pe11.3))
  339 FORMAT (4x,'drift v   ',10(1pe11.3))
  341 FORMAT (5x,'fluxlm   ',10(1pe11.3))
  343 FORMAT (5x,'tr mfp   ',10(1pe11.3))
  345 FORMAT (5x,'dc       ',10(1pe11.3))
  347 FORMAT (4x,'estress   ',10(1pe11.3)/)
  401 FORMAT ('   j           k=1        k=2        k=3        k=4        k=5        k=6        k=7        k=8&
&        k=9        k=10'/)
  403 FORMAT ('   j           k=11       k=12       k=13       k=14       k=15       k=16       k=17       k=18&
&       k=19       k=20'/)
  405 FORMAT ('   j           k=21       k=22       k=23       k=24       k=25       k=26       k=27       k=28&
&         k=29       k=30'/)
  421 FORMAT (4x,'  w    ',10(1pe11.3))
  423 FORMAT (4x,'em surf',3x,10(i3,8x))
  425 FORMAT (2x,'r em surf',10(1pe11.3))
  427 FORMAT (1x,'ro em surf',10(1pe11.3))
  429 FORMAT (2x,'t em surf',10(1pe11.3))
  431 FORMAT (2x,'m em surf',10(1pe11.3)/)
  433 FORMAT (41x,'Luminosity of each energy group (foes/sec)')
  435 FORMAT (11x,10(1pe11.3)/)
  437 FORMAT (42x,'Differential luminosity (foes/sec/mev)')
  439 FORMAT (41x,'Net energy emitted from each energy group (foes)')
  441 FORMAT (39x,'Differential energy radiated from surface (foes/mev)')
  443 FORMAT (1x,10(1pe11.3))
  445 FORMAT (' The ',a1,'-neutrino emission       surface is at r=',1pe10.3, &
&             ' cm, rho=',1pe10.3,' g/cm3, t=',1pe10.3,' k, m=',1pe10.3,' msolar')
  447 FORMAT (' The ',a1,'-neutrino thermalization surface is at r=',1pe10.3, &
&             ' cm, rho=',1pe10.3,' g/cm3, t=',1pe10.3,' k, m=',1pe10.3,' msolar')
  501 FORMAT (' Total ',a1,'-neutrino luminosity from the surface = ',46('.'),1pe10.3,' foes/sec')
  503 FORMAT (' Total energy radiated from core as ',a1,'-neutrinos = ',44('.'),1pe10.3,' foes')
  505 FORMAT (' Total energy of ',a1,'-neutrinos currently residing in the core = ',32('.'),1pe10.3,' foes')
  507 FORMAT (' Total energy transferred to ',a1,'-neutrinos by emission and absorption = ',24('.'),1pe10.3,' foes')
  509 FORMAT (' Total energy transferred to ',a1,'-neutrinos by scattering on electrons = ',24('.'),1pe10.3,' foes')
  511 FORMAT (' Total energy transferred to ',a1,'-neutrinos by the pair annihilation process = ',18('.'),1pe10.3,' foes')
  513 FORMAT (' Total energy transferred to ',a1,'-neutrinos by the bremsstrahlung pair annihilation process = ',3('.'),&
&  1pe10.3,' foes')
  515 FORMAT (' Total energy transferred to ',a1,'-neutrinos by inelastic scattering on nucleons = ',15('.'),1pe10.3,' foes')
  517 FORMAT (' Total energy transferred to ',a1,'-neutrinos by all produc and scat processes = ',18('.'),1pe10.3,' foes')
  519 FORMAT (' Total energy transferred to ',a1,'-neutrinos by matter compression or expansion = ',16('.'),1pe10.3,' foes')
  521 FORMAT (' Number of ',a1,'-neutrinos emitted from core = ',51('.'),1pe10.3)
  523 FORMAT (' Net number of ',a1,'-neutrinos residing in the core at this time = ',31('.'),1pe10.3)
  601 FORMAT (' Number of electrons in the core = ',59('.'),1pe10.3)
  603 FORMAT (' Number of leptons emitted by or residing in core = ',42('.'),1pe10.3/)
 1001 FORMAT (' Allocation problem for array',a10,' in editn')
 2001 FORMAT (' Deallocation problem for array',a10,' in editn')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) RETURN

!-----------------------------------------------------------------------
!  Update edit control parameters.
!   If prnttest = true,
!    (1) add 1 to nedn(n)
!    (2) if nedu(n) < intedu(n), return
!    (3) if nedu(n) >= intedu(n), add 1 to neden(i) and proceed to the edit
!   If prnttest = false,
!    (1) set neden(i) = niedn(i) if niedn(i) gt 0 (to print)
!    (2) set neden(i) = -1       if niedn(i) eq 0 (not to print)
!    (printing is only implemented if neden(i) = niedn(i)).
!-----------------------------------------------------------------------

IF ( prnttest ) THEN
  nedn(n)            = nedn(n) + 1
  IF ( nedn(n) < intedn(n) ) RETURN
  nedn(n)            = 0
  DO i = 1,40
    neden(i)         = neden(i) + 1
  END DO
ELSE
  DO i = 1,40
    nedens(i)        = neden(i)
    IF ( niedn(i) /= 0 ) THEN
      neden(i)       = niedn(i)
    ELSE
      neden(i)       = -1
    END IF ! niedn(i) ne 0
  END DO
END IF ! prnttest

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (psi0l(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0l     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi1l(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1l     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psiteq(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psiteq    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi10(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi10     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (tautrns(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tautrns   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tauopt(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tauopt    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (taueanp(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taueanp   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (taueanc(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taueanc   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (taueahnc(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taueahnc  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxn(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxn     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (tauajk(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tauajk    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (taujk(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taujk     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tausjk(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tausjk    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dnuradjn(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dnuradjn  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (strsnujn(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'strsnujn  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (edotnp(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'edotnp    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (psi0eq(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0eq    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (reanpin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reanpin   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (reanpnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reanpnet  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (reancin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reancin   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (reancnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reancnet  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (reahncin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reahncin  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (reahncnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reahncnet '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rlocal(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rlocal    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rtr(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rtr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (tflux(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tflux     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (vdrift(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vdrift    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ciicra(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ciicra    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rpnesin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnesin   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpnesnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnesnet  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (qnes(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnes      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rnes(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnes      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rpnncsin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnncsin  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpnncsnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnncsnet '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (qnncs(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnncs     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rnncs(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnncs     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rpnnsin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnsin   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpnnsnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnsnet  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (qns(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qns       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rns(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rns       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rpnnnsin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnnsin  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpnnnsnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnnsnet '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (qnns(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnns      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rnns(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnns      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rpnnAsin(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnAsin  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpnnAsnet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnAsnet '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (qnAs(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnAs      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rnAs(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnAs      '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (wk_mag_p(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wk_mag_p  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (wk_mag_n(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wk_mag_n  '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rpain(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpain     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rpanet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpanet    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rbain(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rbain     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (rbanet(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rbanet    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ymefrnp(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymefrnp   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymenc(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymenc     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymehnc(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymehnc    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnps(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnps    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnns(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnns    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnhs(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnhs    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnas(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnas    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnes(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnes    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnncs(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnncs   '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnn(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnn     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnnn(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnnn    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnAn(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnAn    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnts(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnts    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnbs(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnbs    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnps_ed(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnps_ed '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ymdnns_ed(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnns_ed '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ecoefp(nx,nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefp    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (uenulf(nez), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uenulf    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

psi0l                = zero
psi1l                = zero
psiteq               = zero
psi10                = zero

tautrns              = zero
tauopt               = zero
taueanp              = zero
taueanc              = zero
taueahnc             = zero
fluxn                = zero

tauajk               = zero
taujk                = zero
tausjk               = zero

dnuradjn             = zero
strsnujn             = zero
edotnp               = zero

psi0eq               = zero
reanpin              = zero
reanpnet             = zero
reancin              = zero
reancnet             = zero
reahncin             = zero
reahncnet            = zero
rlocal               = zero
rtr                  = zero
tflux                = zero
vdrift               = zero
ciicra               = zero

rpnesin              = zero
rpnesnet             = zero
qnes                 = zero
rnes                 = zero

rpnncsin             = zero
rpnncsnet            = zero
qnncs                = zero
rnncs                = zero

rpnnsin              = zero
rpnnsnet             = zero
qns                  = zero
rns                  = zero

rpnnnsin             = zero
rpnnnsnet            = zero
qnns                 = zero
rnns                 = zero

rpnnAsin             = zero
rpnnAsnet            = zero
qnAs                 = zero
rnAs                 = zero

rpain                = zero
rpanet               = zero
rbain                = zero
rbanet               = zero

ymefrnp              = zero
ymenc                = zero
ymehnc               = zero
ymdnps               = zero
ymdnns               = zero
ymdnhs               = zero
ymdnas               = zero
ymdnes               = zero
ymdnncs              = zero
ymdnn                = zero
ymdnnn               = zero
ymdnAn               = zero
ymdnts               = zero
ymdnbs               = zero
ymdnps_ed            = zero
ymdnns_ed            = zero

ecoefp               = zero
uenulf               = zero

!-----------------------------------------------------------------------
!        Initialize constants.
!
!  fcoef  :   coefficient for computing neutrino number flux ( cm/s /Mev**3 cm**3 )
!  nstate :   coefficient for computing neutrino number density ( /Mev**3 cm**3 )
!  bedot  :   coefficient for computing heating rate (MeV/MeV s)
!-----------------------------------------------------------------------

IF ( first ) THEN
  cd3                = cvel/3.d+00
  tthird             = 2.0d+00/3.0d+00
  fcoef              = ( 4.d0 * pi/3.d0 ) * cvel/( h * cvel )**3
  nstate             = 4.d0 * pi/( h * cvel )**3
  bedot              = cvel * rmu * nstate
  first              = .false.
END IF ! first

!-----------------------------------------------------------------------
!  Initialize text
!-----------------------------------------------------------------------

IF ( n == 1 ) THEN
  type               = 'e-type neutrino data'
  under              = '----------------------'
ELSE IF ( n == 2 ) THEN
  type               = 'e-type antineutrino data'
  under              = '--------------------------'
ELSE IF ( n == 3 ) THEN
  type               = 't-type neutrino data'
  under              = '----------------------'
ELSE IF ( n == 4 ) THEN
  type               = 't-type antineutrino data'
  under              = '--------------------------'
ELSE
  WRITE (nprint,7)
  WRITE (nlog,7)
  STOP
END IF ! n = 1

!-----------------------------------------------------------------------
!  Energy group edit counters
!-----------------------------------------------------------------------

icnt                 = INT( REAL(nnugp(n)-1)/10.e0 ) + 1
kmax(icnt)           = nnugp(n)

!-----------------------------------------------------------------------
!  Initialize variables
!-----------------------------------------------------------------------

DO k = 1,nnugpmx
  DO j = jr_min,jr_max
    ecoefp(j,k)      = ecoef * stwt(n) * unu(j,k)**3 * dunu(j,k)
    tautrns(k)       = zero
    tauopt(k)        = zero
  END DO
END DO

jedit                = idxedn(n)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                    ||||| Begin radial loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

radial: DO jv = jr_min,jr_max
  j                  = jr_max - jv + jr_min

  jp1                = j + 1
  drj                = dr(j)
  aream              = frpi * r(j-1)**2
  areaj              = frpi * r(j)**2

!-----------------------------------------------------------------------
!  T and mmu of an equilibrium neutrino distribution
!-----------------------------------------------------------------------

  CALL enutcp( n, j, enut, enucmp )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!        ||||| Begin loop for sets of 10 energy groups |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

energy_ct:  DO l = 1,icnt

!-----------------------------------------------------------------------
!  Write header
!-----------------------------------------------------------------------

    n_time           = nprint
    IF ( l == 1 ) THEN
      IF ( j == jr_max ) THEN
        WRITE (nprint,1)
        WRITE (nprint,3) head
        WRITE (nprint,5)
        CALL date_and_time_print( n_time )
        WRITE (nprint,101) type
        WRITE (nprint,103) under
      END IF ! j = jr_max
    END IF ! l = 1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                   ||||| Begin energy loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

energy: DO k = 1,nnugp(n)

      w2             = unu(j,k)**2
      cflxnm         = fcoef * w2 * aream
      cfluxn         = fcoef * w2 * areaj

!-----------------------------------------------------------------------
!  Equilibrium n-type neutrino occupation numbers
!-----------------------------------------------------------------------

      IF ( n /= 3 ) THEN
          psi0eq(k)  = emis(j,k,n,ij_ray,ik_ray)/( emis(j,k,n,ij_ray,ik_ray) &
&                    + absor(j,k,n,ij_ray,ik_ray) + epsilon )
      ELSE
        ekt          = unu(j,k)/( kmev * t(j) )
        psi0eq(k)    = one/( one + fexp(ekt) )
      END IF ! n ne 3

!-----------------------------------------------------------------------
!  Ratio of psi1 to psi0
!-----------------------------------------------------------------------

      psi10(k)       = psi1(j,k,n)/( psi0(j,k,n) + epsilon )

!-----------------------------------------------------------------------
!  Equilibrium n-type neutrino occupation numbers assuming separate
!   n-neutrino temperature and chemical potential.
!-----------------------------------------------------------------------

      IF ( enut == zero  .and.  enucmp == zero ) THEN
        psiteq(k)    = zero
      ELSE
        expa         = ( unu(j,k) - enucmp )/( enut + 1.d-5 )
        psiteq(k)    = one/( one + fexp(expa) )
      END IF ! enut = 0  and  enucmp = 0

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by electron (positron) captures on
!   free nucleons (cm^{-3} s^{-1} MeV^{-1}).
!-----------------------------------------------------------------------
      CALL abem_cal( n, unu(j,k), rho(j), t(j), aesv(j,7,ij_ray,ik_ray),      &
& aesv(j,8,ij_ray,ik_ray), aesv(j,9,ij_ray,ik_ray), aesv(j,10,ij_ray,ik_ray), &
& aesv(j,11,ij_ray,ik_ray), aesv(j,4,ij_ray,ik_ray), aesv(j,5,ij_ray,ik_ray), &
& aesv(j,6,ij_ray,ik_ray), absornp, emisnp, ye(j) )
      reanpin(k)     = nstate * cvel * emisnp * w2 * ( one - psi0(j,k,n) )
      reanpnet(k)    = nstate * cvel * w2 * ( emisnp * ( one - psi0(j,k,n) ) - absornp * psi0(j,k,n) )
      ymefrnp(k)     = emisnp

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by electron (positron) captures on
!   heavy nuclei (FFN/B) (cm^{-3} s^{-1} MeV^{-1}).
!-----------------------------------------------------------------------

      CALL abemnc( n, unu(j,k), rho(j), t(j), aesv(j,9,ij_ray,ik_ray),         &
& aesv(j,10,ij_ray,ik_ray), aesv(j,11,ij_ray,ik_ray), aesv(j,4,ij_ray,ik_ray), &
& aesv(j,5,ij_ray,ik_ray), aesv(j,6,ij_ray,ik_ray), absornc, emisnc )
      reancin(k)     = nstate * cvel * emisnc * w2 * ( one - psi0(j,k,n) )
      reancnet(k)    = nstate * cvel * w2 * ( emisnc * ( one - psi0(j,k,n) ) - absornc * psi0(j,k,n) )
      ymenc(k)       = emisnc

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by electron (positron) captures on
!   heavy nuclei (Haxton) (cm^{-3} s^{-1} MeV^{-1}).
!-----------------------------------------------------------------------

      CALL abemhnc( n, k, rho(j), t(j), aesv(j,9,ij_ray,ik_ray), aesv(j,10,ij_ray,ik_ray), &
&      aesv(j,4,ij_ray,ik_ray), aesv(j,5,ij_ray,ik_ray), aesv(j,6,ij_ray,ik_ray), absorhnc, emishnc)
      reahncin(k)    = nstate * cvel * emishnc * w2 * ( one - psi0(j,k,n) )
      reahncnet(k)   = nstate * cvel * w2 * ( emishnc * ( one - psi0(j,k,n) ) - absorhnc * psi0(j,k,n) )
      ymehnc(k)      = emishnc

!-----------------------------------------------------------------------
!  e-neutrino - free nucleon absorption and emission inverse mean free
!   paths (cm^{-1}).
!-----------------------------------------------------------------------

      taueanp(k)     = absornp

!-----------------------------------------------------------------------
!  e-neutrino - nuclei absorption and emission inverse mean free paths
!   (FFN/B) (cm^{-1}).
!-----------------------------------------------------------------------

      taueanc(k)     = absornc

!-----------------------------------------------------------------------
!  e-neutrino - nuclei absorption and emission inverse mean free paths
!   (Haxton) (cm^{-1}).
!-----------------------------------------------------------------------

      taueahnc(k)    = absorhnc

!-----------------------------------------------------------------------
!  Rate of energy transfer from neutrinos to matter per baryon per
!   unit energy (Mev/MeV s) by neutrino absorption  and emission on
!   free neutrons and protons.
!
!    de              3
!   ----- = -coef * w  [emisnp *(1 - psi0 ) - absrnp *psi0 ]
!   dtime            k        k          k          k     k
!-----------------------------------------------------------------------

      edotnp(k)      = ( bedot/rho(j) ) * unu(j,k)**3 * ( emisnp - ( emisnp + absornp ) * psi0(j,k,n) )

!-----------------------------------------------------------------------
!  Number luminosity / energy of n-type neutrinos across zone
!   j (s^{-1} MeV^{-1}). 
!-----------------------------------------------------------------------

      fluxn(k)       = cfluxn * psi1(j,k,n)

!-----------------------------------------------------------------------
!  Total number of e-type neutrinos per unit energy that have crossed
!   zone j (MeV^{-1}).
!-----------------------------------------------------------------------

      dnuradjn(k)    = dnurad(j,k,n,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Force per unit mass per unit energy exerted on the matter by
!   n-neutrinos in radial zone j and energy zone k 
!   (dynes g^{-1} MeV^{-1}).
!-----------------------------------------------------------------------

      strsnujn(k)    = strsnu_x(j,k,n)

!-----------------------------------------------------------------------
!  Rate of e-neutrino "production" by transport.
!-----------------------------------------------------------------------

      rtr(k)         = ( cflxnm * psi1(j-1,k,n) - cfluxn * psi1(j,k,n) )/( dmrst(j)/rho(j) )

!-----------------------------------------------------------------------
!  Logarithms of the e-neutrino occupation numbers
!-----------------------------------------------------------------------

      psi0l(k)       = DLOG10( DABS( psi0(j,k,n) ) + epsilon )
      psi1l(k)       = DLOG10( DABS( psi1(j,k,n) ) + epsilon )

!-----------------------------------------------------------------------
!  n-type isoenergetic neutrino scattering inverse mean free paths (cm^{-1}).
!-----------------------------------------------------------------------

      xhe            = DMAX1( one - aesv(j,7,ij_ray,ik_ray) - aesv(j,8,ij_ray,ik_ray) &
&                    - aesv(j,9,ij_ray,ik_ray), zero )

      CALL scatical( rho(j), t(j), unu(j,k), aesv(j,7,ij_ray,ik_ray),      &
&                   aesv(j,8,ij_ray,ik_ray), xhe, aesv(j,9,ij_ray,ik_ray), &
&                   aesv(j,10,ij_ray,ik_ray), aesv(j,11,ij_ray,ik_ray),    &
&                   rmdnns, rmdnps, rmdnbns, rmdnbps, rmdnhes, rmdnhs,     &
&                   coh, cohb)

      IF ( n == 1  .or.  n == 3 ) THEN
        ymdnps(k)    = rmdnps
        ymdnns(k)    = rmdnns
      ELSE
        ymdnps(k)    = rmdnbps
        ymdnns(k)    = rmdnbns
      END IF ! n == 1  .or.  n == 3

      ymdnhs(k)      = rmdnhes
      ymdnas(k)      = rmdnhs

      CALL scatical_ed( rho(j), t(j), unu(j,k), aesv(j,7,ij_ray,ik_ray),   &
&                   aesv(j,8,ij_ray,ik_ray), xhe, aesv(j,9,ij_ray,ik_ray), &
&                   aesv(j,10,ij_ray,ik_ray), aesv(j,11,ij_ray,ik_ray),    &
&                   rmdnns_ed, rmdnps_ed, dum1, dum2, dum3 )

      ymdnps_ed(k)   = rmdnps_ed
      ymdnns_ed(k)   = rmdnns_ed

!-----------------------------------------------------------------------
!  Weak magnetism corrections to neutral current scattering on nucleons
!-----------------------------------------------------------------------

      CALL nc_weak_mag( unu(j,k), xi_p_wm, xi_n_wm, xib_p_wm, xib_n_wm )

      IF ( n == 1  .or.  n == 3 ) THEN
        wk_mag_p(k)  = xi_p_wm
        wk_mag_n(k)  = xi_n_wm
      ELSE
        wk_mag_p(k)  = xib_p_wm
        wk_mag_n(k)  = xib_n_wm
      END IF ! n == 1  .or.  n == 3

!-----------------------------------------------------------------------
!  n-type neutrino-electron scattering inverse mean free paths (cm^{-1}).
!-----------------------------------------------------------------------

      CALL sctednes( n, j, ij_ray,ik_ray, k, rho(j), t(j), ye(j), a0w,  &
&           b0w, c0w, a1w, b1w, c1w, rmdnes0, rmdnes1, rmdnes, arnes,   &
&           brnes, crnes, aenes, benes, cenes, arscte, brscte, crscte,  &
&           arsctu, brsctu, crsctu, arscti, brscti, crscti, arsctd,     &
&           brsctd, crsctd, aesctu, besctu, cesctu, aescti, bescti,     &
&           cescti, aesctd, besctd, cesctd )

      ymdnes(k)      = rmdnes
      dpsi0nesin     = arscte * psi0(j,k,n) + brscte * psi1(j,k,n) + crscte
      dpsi0nesnet    = a0w    * psi0(j,k,n) + b0w    * psi1(j,k,n) + c0w

!-----------------------------------------------------------------------
!  Rate of nes per n-neutrino of energy group k (/sec) Fraction of
!   n-neutrino energy transferred to electrons by nes per scattering
!   event
!-----------------------------------------------------------------------

      rnes(k)        = cvel * ( arnes * psi0(j,k,n) )/( psi0(j,k,n) + epsilon )
      efinal         = aenes
      einitial       = arnes * unu(j,k)
      qnes(k)        = - ( efinal - einitial )/( einitial + epsilon )

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by neutrino-electron "in" scattering
!   per unit volume per unit energy, and net rate of n-neutrino
!   production by neutrino-electron scattering (1.e-30 /cm3*sec*mev).
!-----------------------------------------------------------------------

      rpnesin(k)     = nstate * cvel * dpsi0nesin  * w2
      rpnesnet(k)    = nstate * cvel * dpsi0nesnet * w2

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering inverse mean free paths (cm^{-1}).
!-----------------------------------------------------------------------

      CALL sctednncs( n, j, k, rho(j), t(j), ye(j), a0w, b0w, c0w, a1w, b1w, c1w, &
&          rmdnncs0, rmdnncs1, rmdnncs, arnncs, brnncs, crnncs, aenncs, benncs, cenncs, &
&          arscte, brscte, crscte, arsctu, brsctu, crsctu, arscti, brscti, crscti, &
&          arsctd, brsctd, crsctd, aesctu, besctu, cesctu, aescti, bescti, cescti, &
&          aesctd, besctd, cesctd)

      ymdnncs(k)     = rmdnncs
      dpsi0nncsin    = arscte * psi0(j,k,n) + brscte * psi1(j,k,n) + crscte
      dpsi0nncsnet   = a0w    * psi0(j,k,n) + b0w    * psi1(j,k,n) + c0w   

!-----------------------------------------------------------------------
!  Rate of nncs per n-neutrino of energy group k (s^{-1}) Fraction of
!   n-neutrino energy transferred to electrons by nncs per scattering
!   event.
!-----------------------------------------------------------------------

      rnncs(k)       = cvel * ( arnncs * psi0(j,k,n) )/( psi0(j,k,n) + epsilon )
      efinal         = aenncs
      einitial       = arnncs * unu(j,k)
      qnncs(k)       = - ( efinal - einitial )/( einitial + epsilon )

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by neutrino-nucleus inelastic "in"
!   scattering per unit volume per unit energy, and net rate of 
!  n-neutrino production by neutrino-nucleus inelastic scattering
!   per unit volume per unit energy (1.e-30 /cm3*sec*mev).
!-----------------------------------------------------------------------

      rpnncsin(k)    = nstate * cvel * dpsi0nncsin  * w2
      rpnncsnet(k)   = nstate * cvel * dpsi0nncsnet * w2

!-----------------------------------------------------------------------
!  Neutrino-nucleon elastic scattering inverse mean free paths (cm^{-1}).
!-----------------------------------------------------------------------

      CALL sctedns( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), a0w, &
&           b0w, c0w, a1w, b1w, c1w, rmdnn0, rmdnn1, rmdnn, arnes, brnes, &
&           crnes, aenes, benes, cenes, arscte, brscte, crscte, arsctu, &
&           brsctu, crsctu, arscti, brscti, crscti, arsctd, brsctd, crsctd, &
&           aesctu, besctu, cesctu, aescti, bescti, cescti, aesctd, besctd, &
&           cesctd )

      ymdnn(k)       = rmdnn

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by neutrino-nucleon elastic "in"
!   scattering per unit volume per unit energy, and net rate of
!   n-neutrino production by neutrino-nucleon elastic scattering
!   per unit volume per unit energy. (cm^{-3} s^{-1} MeV^{-1}).
!-----------------------------------------------------------------------

      dpsi0nnsin     = arscte * psi0(j,k,n) + brscte * psi1(j,k,n) + crscte
      dpsi0nnsnet    = a0w    * psi0(j,k,n) + b0w    * psi1(j,k,n) + c0w   
      rpnnsin(k)     = nstate * cvel * dpsi0nnsin   * w2
      rpnnsnet(k)    = nstate * cvel * dpsi0nnsnet  * w2

!-----------------------------------------------------------------------
!  Rate of neutrino-nucleon elastic scattering per n-neutrino of energy
!   group k (/sec).
!  Fraction of n-neutrino energy transferred to nucleons by 
!  neutrino-nucleon elastic scattering per scattering event.
!-----------------------------------------------------------------------

      rns(k)         = cvel * ( arnes * psi0(j,k,n) )/( psi0(j,k,n) + epsilon )
      efinal         = aenes
      einitial       = arnes * unu(j,k)
      qns(k)         = - ( efinal - einitial )/( einitial + epsilon )

!-----------------------------------------------------------------------
!  Neutrino-nucleon inelastic  scatteringinverse mean free paths (cm^{-1}).
!-----------------------------------------------------------------------

      CALL sctednns( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), a0w, b0w, &
&           c0w, a1w, b1w, c1w, rmdnnn0, rmdnnn1, rmdnnn, arnes, brnes, crnes, &
&           aenes, benes, cenes, arscte, brscte, crscte, arsctu, brsctu, crsctu, &
&           arscti, brscti, crscti, arsctd, brsctd, crsctd, aesctu, besctu,  &
&           cesctu, aescti, bescti, cescti,aesctd, besctd, cesctd )

      ymdnnn(k)      = rmdnnn

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by neutrino-nucleon inelastic "in"
!   scattering per unit volume per unit energy, and net rate of 
!   n-neutrino production by neutrino-nucleon inelastic scattering
!   per unit volume per unit energy (cm^{-3} s^{-1} MeV^{-1}).
!-----------------------------------------------------------------------

      dpsi0nnnsin    = arscte * psi0(j,k,n) + brscte * psi1(j,k,n) + crscte
      dpsi0nnnsnet   = a0w    * psi0(j,k,n) + b0w    * psi1(j,k,n) + c0w   
      rpnnnsin(k)    = nstate * cvel * dpsi0nnnsin  * w2
      rpnnnsnet(k)   = nstate * cvel * dpsi0nnnsnet * w2

!-----------------------------------------------------------------------
!  Rate of neutrino-nucleon inelastic scattering per n-neutrino of
!   energy group k (/sec).
!  Fraction of n-neutrino energy transferred to nucleons by
!   neutrino-nucleon inelastic scattering per scattering event.
!-----------------------------------------------------------------------

      rnns(k)      = cvel * ( arnes * psi0(j,k,n) )/( psi0(j,k,n) + epsilon )
      efinal       = aenes
      einitial     = arnes * unu(j,k)
      qnns(k)      = - ( efinal - einitial )/( einitial + epsilon )

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scatteringinverse mean free paths (cm^{-1}).
!-----------------------------------------------------------------------

      CALL sctednAs( n, j, ij_ray, ik_ray, k, rho(j), t(j), ye(j), a0w, b0w, &
&           c0w, a1w, b1w, c1w, rmdnnA0, rmdnnA1, rmdnnA, arnes, brnes, crnes, &
&           aenes, benes, cenes, arscte, brscte, crscte, arsctu, brsctu, crsctu, &
&           arscti, brscti, crscti, arsctd, brsctd, crsctd, aesctu, besctu, &
&           cesctu, aescti, bescti, cescti, aesctd, besctd, cesctd )

      ymdnAn(k)      = rmdnnA

!-----------------------------------------------------------------------
!  Rate of n-neutrino production by neutrino-nucleus inelastic "in"
!   scattering per unit volume per unit energy, and net rate of
!   n-neutrino production by neutrino-nucleus inelastic scattering per
!   unit volume per unit energy (cm^{-3} s^{-1} MeV^{-1}).
!-----------------------------------------------------------------------

      dpsi0nnAsin    = arscte * psi0(j,k,n) + brscte * psi1(j,k,n) + crscte
      dpsi0nnAsnet   = a0w    * psi0(j,k,n) + b0w    * psi1(j,k,n) + c0w   
      rpnnAsin(k)    = nstate * cvel * dpsi0nnAsin  * w2
      rpnnAsnet(k)   = nstate * cvel * dpsi0nnAsnet * w2

!-----------------------------------------------------------------------
!  Rate of neutrino-nucleus inelastic scattering per n-neutrino of 
!   energy group k (/sec).
!  Fraction of n-neutrino energy transferred to nuclei by
!   neutrino-nucleus inelastic scattering per scattering event.
!-----------------------------------------------------------------------

      rnAs(k)      = cvel * ( arnes * psi0(j,k,n) )/( psi0(j,k,n) + epsilon )
      efinal       = aenes
      einitial     = arnes * unu(j,k)
      qnAs(k)      = - ( efinal - einitial )/( einitial + epsilon )

!-----------------------------------------------------------------------
!  Ion-ion correlation correction to isoenergetic neutrino-nucleus
!   scattering.
!-----------------------------------------------------------------------

     CALL scatiicr( rho(j), t(j), unu(j,k), aesv(j,9,ij_ray,ik_ray),     &
&          aesv(j,10,ij_ray,ik_ray), aesv(j,11,ij_ray,ik_ray), ciicra(k) )

!-----------------------------------------------------------------------
!  Pair production inverse mean free paths.
!-----------------------------------------------------------------------

      CALL paired( n, j, ij_ray, ik_ray, k, r(j), rho(j), t(j), ye(j), a0w, &
&           b0w, c0w, a1w, b1w, c1w, rmdnts0, rmdnts1, rmdnts, artpe, brtpe, &
&           crtpe, artae, brtae, crtae )

      ymdnts(k)      = rmdnts

!-----------------------------------------------------------------------
!  Net rate of n-neutrino production by electron-positron pair
!   annihilation per unit volume per unit energy (/cm3*sec*mev).
!-----------------------------------------------------------------------

      dpsi0pain      = artpe * psi0(j,k,n) + brtpe * psi1(j,k,n) + crtpe
      dpsi0panet     = a0w   * psi0(j,k,n) + b0w   * psi1(j,k,n) + c0w  
      rpain(k)       = nstate * cvel * dpsi0pain  * w2
      rpanet(k)      = nstate * cvel * dpsi0panet * w2

!-----------------------------------------------------------------------
!  Bremsstrahlung pair production inverse mean free paths.
!-----------------------------------------------------------------------

      CALL bremed( n, j, ij_ray, ik_ray, k, r(j), rho(j), t(j), ye(j), a0w,  &
&           b0w, c0w, a1w, b1w, c1w, rmdnbs0, rmdnbs1, rmdnbs, artpe, brtpe, &
&           crtpe, artae, brtae, crtae )

      ymdnbs(k)      = rmdnbs

!-----------------------------------------------------------------------
!  Net rate of n-neutrino production by bremsstrahlungn pair
!   annihilation per unit volume per unit energy (/cm3*sec*mev).
!-----------------------------------------------------------------------

      dpsi0bain      = artpe * psi0(j,k,n) + brtpe * psi1(j,k,n) + crtpe
      dpsi0banet     = a0w   * psi0(j,k,n) + b0w   * psi1(j,k,n) + c0w  
      rbain(k)       = nstate * cvel * dpsi0bain  * w2
      rbanet(k)      = nstate * cvel * dpsi0banet * w2

!-----------------------------------------------------------------------
!  Total rate of n-neutrino production due to all local sources.
!-----------------------------------------------------------------------

      rlocal(k)    = reanpnet(k) + reancnet(k) + reahncnet(k) + rpnesnet(k) &
&                  + rpnncsnet(k) + rpnnsnet(k) + rpnnnsnet(k) + rpanet(k)  &
&                  + rbanet(k) + dndt_v(j,k,n,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Total absorption, emission ond nonisoenergetic scattering optical
!   depth across a mass zone.
!-----------------------------------------------------------------------

      tauajk(j,k)  = ( taueanp(k) + taueanc(k) + taueahnc(k) + ymdnes(k) &
&                  + ymdnncs(k) + ymdnts(k) + ymdnbs(k) + ymdnnn(k)      &
&                  + ymdnAn(k) ) * drj

!-----------------------------------------------------------------------
!  Total isoenergetic optical depth across a mass zone.
!-----------------------------------------------------------------------

      tausjk(j,k)  = ( ymdnhs(k) + ymdnas(k) ) * drj

!-----------------------------------------------------------------------
!  Total transport optical depth across a mass zone and the total
!   transport optical depth from the surface.
!-----------------------------------------------------------------------

      taujk(j,k)   = ( taueanp(k) + taueanc(k) + taueahnc(k) + ymdnhs(k)  &
&                  + ymdnas(k) + ymdnes(k) + ymdnncs(k) + ymdnts(k) + ymdnbs(k) &
&                  + ymdnn(k) + ymdnnn(k) + ymdnAn(k) ) * drj
      tautrns(k)   = tautrns(k) + taujk(j,k)

!-----------------------------------------------------------------------
!  Total neutrino outbeam optical depth from the surface.
!-----------------------------------------------------------------------

      tauopt(k)    = tauopt(k) + ( taueanp(k) + taueanc(k) + taueahnc(k) &
&                  + rmdnhes0 + rmdnhs0 + rmdnes0 + rmdnncs0 + rmdnts0   &
&                  + rmdnbs0 + rmdnn0 + rmdnnn0 ) * drj

!-----------------------------------------------------------------------
!  e-neutrino drift velocity.
!-----------------------------------------------------------------------

      vdrift(k)    = cvel * psi1(j,k,n)/( psi0(j,k,n) + epsilon )
      

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                   ||||| End energy loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    END DO energy

!-----------------------------------------------------------------------
!
!                         \\\\\ EDIT /////
!
!-----------------------------------------------------------------------

    IF ( jedit == idxedn(n)  .or.  j == jr_min ) THEN

      WRITE (nprint,201) j,( unu(j,k),     k = kmin(l),kmax(l) )
      WRITE (nprint,203) ( dunu(j,k),      k = kmin(l),kmax(l) )

      IF ( neden(1) == niedn(1) )   WRITE (nprint,205) ( psi0eq(k),      k = kmin(l),kmax(l) )
      IF ( neden(2) == niedn(2) )   WRITE (nprint,207) ( psi0(j,k,n),    k = kmin(l),kmax(l) )
      IF ( neden(3) == niedn(3) )   WRITE (nprint,209) ( psi0l(k),       k = kmin(l),kmax(l) )

      IF ( j == jr_max  .and. neden(4) == niedn(4) ) WRITE (nprint,211) ( psi0(jp1,k,n),  k = kmin(l),kmax(l) )

      IF ( neden(5) == niedn(5) )   WRITE (nprint,213) ( psi1(j,k,n),    k = kmin(l),kmax(l) )
      IF ( neden(6) == niedn(6) )   WRITE (nprint,215) ( psi1l(k),       k = kmin(l),kmax(l) )
      IF ( neden(7) == niedn(7) )   WRITE (nprint,217) ( flxf(j,k,n),    k = kmin(l),kmax(l) )
      IF ( neden(7) == niedn(7) )   WRITE (nprint,219) ( Ef(j,k,n),      k = kmin(l),kmax(l) )
      IF ( neden(8) == niedn(8) )   WRITE (nprint,221) ( psiteq(k),      k = kmin(l),kmax(l) )
      IF ( neden(9) == niedn(9) )   WRITE (nprint,223) ( reanpin(k),     k = kmin(l),kmax(l) )
      IF ( neden(9) == niedn(9) )   WRITE (nprint,225) ( reanpnet(k),    k = kmin(l),kmax(l) )
      IF ( neden(10) == niedn(10) ) WRITE (nprint,227) ( reancin(k),     k = kmin(l),kmax(l) )
      IF ( neden(10) == niedn(10) ) WRITE (nprint,229) ( reancnet(k),    k = kmin(l),kmax(l) )
      IF ( neden(11) == niedn(11) ) WRITE (nprint,231) ( reahncin(k),    k = kmin(l),kmax(l) )
      IF ( neden(11) == niedn(11) ) WRITE (nprint,233) ( reahncnet(k),   k = kmin(l),kmax(l) )
      IF ( neden(12) == niedn(12) ) WRITE (nprint,235) ( rpnesin(k),     k = kmin(l),kmax(l) )
      IF ( neden(12) == niedn(12) ) WRITE (nprint,237) ( rpnesnet(k),    k = kmin(l),kmax(l) )
      IF ( neden(13) == niedn(13) ) WRITE (nprint,239) ( rpnncsin(k),    k = kmin(l),kmax(l) )
      IF ( neden(13) == niedn(13) ) WRITE (nprint,241) ( rpnncsnet(k),   k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,243) ( rpain(k),       k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,245) ( rpanet(k),      k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,247) ( rbain(k),       k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,249) ( rbanet(k),      k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,251) ( rpnnsin(k),     k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,253) ( rpnnsnet(k),    k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,255) ( rpnnnsin(k),    k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,257) ( rpnnnsnet(k),   k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,259) ( dndt_v(j,k,n,ij_ray,ik_ray), k = kmin(l),kmax(l) )
      IF ( neden(14) == niedn(14) ) WRITE (nprint,261) ( rlocal(k),      k = kmin(l),kmax(l) )
      IF ( neden(15) == niedn(15) ) WRITE (nprint,263) ( rtr(k),         k = kmin(l),kmax(l) )
      IF ( neden(16) == niedn(16) ) WRITE (nprint,265) ( edotnp(k),      k = kmin(l),kmax(l) )
      IF ( neden(17) == niedn(17) ) WRITE (nprint,267) ( rnes(k),        k = kmin(l),kmax(l) )
      IF ( neden(18) == niedn(18) ) WRITE (nprint,269) ( qnes(k),        k = kmin(l),kmax(l) )
      IF ( neden(19) == niedn(19) ) WRITE (nprint,271) ( rnncs(k),       k = kmin(l),kmax(l) )
      IF ( neden(20) == niedn(20) ) WRITE (nprint,273) ( qnncs(k),       k = kmin(l),kmax(l) )
      IF ( neden(19) == niedn(19) ) WRITE (nprint,275) ( rns(k),         k = kmin(l),kmax(l) )
      IF ( neden(19) == niedn(19) ) WRITE (nprint,277) ( qns(k),         k = kmin(l),kmax(l) )
      IF ( neden(20) == niedn(20) ) WRITE (nprint,279) ( rnns(k),        k = kmin(l),kmax(l) )
      IF ( neden(20) == niedn(20) ) WRITE (nprint,281) ( qnns(k),        k = kmin(l),kmax(l) )
      IF ( neden(20) == niedn(20) ) WRITE (nprint,283) ( rnAs(k),        k = kmin(l),kmax(l) )
      IF ( neden(20) == niedn(20) ) WRITE (nprint,285) ( qnAs(k),        k = kmin(l),kmax(l) )
      IF ( neden(21) == niedn(21) ) WRITE (nprint,287) ( ymefrnp(k),     k = kmin(l),kmax(l) )
      IF ( neden(22) == niedn(22) ) WRITE (nprint,289) ( ymenc(k),       k = kmin(l),kmax(l) )
      IF ( neden(23) == niedn(23) ) WRITE (nprint,291) ( ymehnc(k),      k = kmin(l),kmax(l) )
      IF ( neden(24) == niedn(24) ) WRITE (nprint,293) ( taueanp(k),     k = kmin(l),kmax(l) )
      IF ( neden(25) == niedn(25) ) WRITE (nprint,295) ( taueanc(k),     k = kmin(l),kmax(l) )
      IF ( neden(26) == niedn(26) ) WRITE (nprint,297) ( taueahnc(k),    k = kmin(l),kmax(l) )
      IF ( neden(27) == niedn(27) ) WRITE (nprint,299) ( ymdnps(k),      k = kmin(l),kmax(l) )
      IF ( neden(27) == niedn(27) ) WRITE (nprint,301) ( ymdnps_ed(k),   k = kmin(l),kmax(l) )
      IF ( neden(28) == niedn(28) ) WRITE (nprint,303) ( ymdnns(k),      k = kmin(l),kmax(l) )
      IF ( neden(28) == niedn(28) ) WRITE (nprint,305) ( ymdnns_ed(k),   k = kmin(l),kmax(l) )
      IF ( neden(29) == niedn(29) ) WRITE (nprint,307) ( ymdnhs(k),      k = kmin(l),kmax(l) )
      IF ( neden(30) == niedn(30) ) WRITE (nprint,309) ( ymdnas(k),      k = kmin(l),kmax(l) )
      IF ( neden(31) == niedn(31) ) WRITE (nprint,311) ( ymdnes(k),      k = kmin(l),kmax(l) )
      IF ( neden(32) == niedn(32) ) WRITE (nprint,313) ( ymdnncs(k),     k = kmin(l),kmax(l) )
      IF ( neden(33) == niedn(33) ) WRITE (nprint,315) ( ymdnts(k),      k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,317) ( ymdnbs(k),      k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,319) ( ymdnn(k),       k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,321) ( ymdnnn(k),      k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,323) ( ymdnAn(k),      k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,325) ( wk_mag_p(k),    k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,327) ( wk_mag_n(k),    k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,329) ( tautrns(k),     k = kmin(l),kmax(l) )
      IF ( neden(34) == niedn(34) ) WRITE (nprint,331) ( tauopt(k),      k = kmin(l),kmax(l) )
      IF ( neden(35) == niedn(35) ) WRITE (nprint,333) ( ciicra(k),      k = kmin(l),kmax(l) )
      IF ( neden(36) == niedn(36) ) WRITE (nprint,335) ( fluxn(k),       k = kmin(l),kmax(l) )
      IF ( neden(37) == niedn(37) ) WRITE (nprint,237) ( dnuradjn(k),    k = kmin(l),kmax(l) )
      IF ( neden(38) == niedn(38) ) WRITE (nprint,339) ( vdrIFt(k),      k = kmin(l),kmax(l) )
      IF ( neden(39) == niedn(39) ) WRITE (nprint,341) ( fluxlm(j,k,n),  k = kmin(l),kmax(l) )
      IF ( neden(39) == niedn(39) ) WRITE (nprint,343) ( tmfp_j(j,k,n),  k = kmin(l),kmax(l) )
      IF ( neden(39) == niedn(39) ) WRITE (nprint,345) ( dc(j,k,n),      k = kmin(l),kmax(l) )
      IF ( neden(40) == niedn(40) ) WRITE (nprint,347) ( strsnujn(k),    k = kmin(l),kmax(l) )

    END IF ! jedit = idxedn(n)  or  j = jr_min

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!          ||||| End loop for sets of 10 energy groups |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  END DO energy_ct

  IF ( jedit == idxedn(n) ) jedit = 0
  jedit            = jedit + 1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!              ||||| End radial loop |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END DO radial

!-----------------------------------------------------------------------
!
!                  \\\\\ EMISSION SUEFACES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!         ||||| Begin loop for sets of 10 energy groups |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO l = 1,icnt

  IF (l == 1) THEN
    WRITE (nprint,1)
    WRITE (nprint,3) head
    WRITE (nprint,5)
    CALL date_and_time_print( n_time )
    WRITE (nprint,101) type
    WRITE (nprint,103) under
    WRITE (nprint,401)
  END IF ! l = 1

  IF ( l == 2 ) THEN
    WRITE (nprint,105)
    WRITE (nprint,403)
  END IF ! l = 2

  IF (l == 3) THEN
    WRITE (nprint,105)
    WRITE (nprint,405)
  END IF ! l = 3

  WRITE (nprint,421) ( unu(j,k) ,            k = kmin(l),kmax(l) )
  WRITE (nprint,423) ( j_sphere(k,n,ij_ray,ik_ray),  k = kmin(l),kmax(l) )
  WRITE (nprint,425) ( r_sphere(k,n,ij_ray,ik_ray),  k = kmin(l),kmax(l) )
  WRITE (nprint,427) ( d_sphere(k,n,ij_ray,ik_ray),  k = kmin(l),kmax(l) )
  WRITE (nprint,429) ( t_sphere(k,n,ij_ray,ik_ray),  k = kmin(l),kmax(l) )
  WRITE (nprint,431) ( m_sphere(k,n,ij_ray,ik_ray),  k = kmin(l),kmax(l) )

!-----------------------------------------------------------------------
!  Total energy flux in e-type neutrinos
!-----------------------------------------------------------------------

  tfluxc           = frpi * r(jr_max)**2

  DO k = kmin(l),kmax(l)
    tflux(k)       = tfluxc * psi1(jr_max,k,n) * fcoef * unu(j,k)**3 * dunu(j,k) * ergfoe
    uenulf(k)      = unukrad(k,n,ij_ray,ik_ray) * ergfoe
  END DO

  WRITE (nprint,433)
  WRITE (nprint,435) ( tflux(k),    k = kmin(l),kmax(l) )
  WRITE (nprint,437)
  WRITE (nprint,435) ( tflux(k)/( dunu(j,k) ), k = kmin(l),kmax(l) )
  WRITE (nprint,439)
  WRITE (nprint,435) ( uenulf(k),   k = kmin(l),kmax(l) )
  WRITE (nprint,441)
  WRITE (nprint,435) ( uenulf(k)/( dunu(j,k) ), k = kmin(l),kmax(l) )

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!         ||||| End loop for sets of 10 energy groups |||||
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

END DO

!-----------------------------------------------------------------------
!  Number and energy transfers of n-type neutrinos
!-----------------------------------------------------------------------

wenulr               = unulum(n)                 * ergfoe
wenul                = unurad  (n,ij_ray,ik_ray) * ergfoe
wenut                = unucr(n)                  * ergfoe
wenuta               = unucrea (n,ij_ray,ik_ray) * ergfoe
wenuts               = unucrnis(n,ij_ray,ik_ray) * ergfoe
wenutp               = unucrpa (n,ij_ray,ik_ray) * ergfoe
wenutb               = unucrba (n,ij_ray,ik_ray) * ergfoe
wenutnn              = unucrnns(n,ij_ray,ik_ray) * ergfoe
wenutt               = unucrt  (n,ij_ray,ik_ray) * ergfoe
wenutv               = unucrv  (n,ij_ray,ik_ray) * ergfoe
wnnurad              = nnurad  (n,ij_ray,ik_ray)
wnnucr               = nnucr   (n,ij_ray,ik_ray)

WRITE (nprint,105)                                                      
WRITE (nprint,501) nutype(n),wenulr                                              
WRITE (nprint,503) nutype(n),wenul                                               
WRITE (nprint,505) nutype(n),wenut                                               
WRITE (nprint,507) nutype(n),wenuta                                              
WRITE (nprint,509) nutype(n),wenuts                                              
WRITE (nprint,511) nutype(n),wenutp
WRITE (nprint,513) nutype(n),wenutb
WRITE (nprint,515) nutype(n),wenutnn
WRITE (nprint,517) nutype(n),wenutt
WRITE (nprint,519) nutype(n),wenutv
WRITE (nprint,521) nutype(n),wnnurad 
WRITE (nprint,523) nutype(n),wnnucr 

!-----------------------------------------------------------------------
!  Total lepton number
!-----------------------------------------------------------------------

elecn                = zero

DO j = 2,jr_max
  elecn              = elecn + dmrst(j) * ye(j)/rmu
END DO

totlpn               = nnurad(1,ij_ray,ik_ray) + nnucr(1,ij_ray,ik_ray) &
&                    + elecn - nnurad(2,ij_ray,ik_ray) - nnucr(2,ij_ray,ik_ray)

WRITE (nprint,601) elecn
WRITE (nprint,603) totlpn

!-----------------------------------------------------------------------
!  Planck averaged emission surfaces
!-----------------------------------------------------------------------

taut                 = zero
itaut                = .false.
jj                   = jr_min

DO j = jr_max,jr_min,-1

  tautp              = taut
  ymdavt             = zero
  flux               = zero

  DO k = 1,nnugpmx
    ymdavt           = ymdavt + taujk(j,k) * cd3 * ecoefp(j,k) * half * ( psi1(j-1,k,n) + psi1(j,k,n) )
    flux             = flux   +              cd3 * ecoefp(j,k) * half * ( psi1(j-1,k,n) + psi1(j,k,n) )
  END DO

  ymdave             = ymdavt/( flux + epsilon )
  taut               = taut + ymdave

  IF ( taut >= tthird ) THEN
    itaut            = .true.
    jj               = j
    EXIT
  END IF ! taut ge one

END DO

IF ( itaut ) THEN

  IF ( jj == jr_min ) THEN
    remsf            = r(jr_min)
    demsf            = rho(jr_min)
    gemst            = rstmss(jr_min)
    temsf            = t(jr_min)
  ELSE
    remsf            = r(jj)      + ( r(jj-1)      - r(jj)      ) * ( tthird - tautp )/( taut - tautp + epsilon )
    demsf            = rho(jj)    + ( rho(jj-1)    - rho(jj)    ) * ( tthird - tautp )/( taut - tautp + epsilon )
    gemst            = rstmss(jj) + ( rstmss(jj-1) - rstmss(jj) ) * ( tthird - tautp )/( taut - tautp + epsilon )
    temsf            = t(jj)      + ( t(jj-1)      - t(jj)      ) * ( tthird - tautp )/( taut - tautp + epsilon )
  END IF ! jj eq jr_min

ELSE
  remsf             = zero
  demsf             = rho(jr_min)
  temsf             = t(jr_min)
  gemst             = zero

END IF ! itaut

gemsf               = gemst/msolar

WRITE (nprint,105)
WRITE (nprint,445) nutype(n),remsf,demsf,temsf,gemsf

!-----------------------------------------------------------------------
!  Planck averaged thermalization surfaces
!-----------------------------------------------------------------------

ymdavat              = zero
ymdavst              = zero
flux                 = zero

DO k = 1,nnugpmx
  ymdavat            = ymdavat + tauajk(jj,k) * cd3 * ecoefp(jj,k) * half * ( DABS(psi1(jj-1,k,n)) &
&                    + DABS(psi1(jj,k,n)) )
  ymdavst            = ymdavst + tausjk(jj,k) * cd3 * ecoefp(jj,k) * half * ( DABS(psi1(jj-1,k,n)) &
&                    + DABS(psi1(jj,k,n)) )
  flux               = flux + cd3 * ecoefp(jj,k) * half * ( psi1(jj-1,k,n) + psi1(jj,k,n) )
END DO

ymdava               = ymdavat/( dabs(flux) + epsilon )
ymdavs               = ymdavst/( dabs(flux) + epsilon )
tautherm             = DSQRT( (ymdava + ymdavs )/( 3.d0 * ymdava + epsilon ) )

taut                 = zero
itaut                = .false.

DO j = jr_max,jr_min,-1
  tautp              = taut
  ymdavt             = zero
  flux               = zero

  DO k = 1,nnugpmx
    ymdavt           = ymdavt + taujk(j,k)*cd3*ecoefp(j,k) * half * ( DABS(psi1(j-1,k,n)) + DABS(psi1(j,k,n)) )
    flux             = flux + cd3 * ecoefp(j,k) * half * ( DABS(psi1(j-1,k,n)) + DABS(psi1(j,k,n)) )
  END DO

  ymdave             = ymdavt/( dabs(flux) + epsilon )
  taut               = taut + ymdave

  IF ( taut >= tautherm ) THEN
    itaut            = .true.
    jj               = j
    EXIT
  END IF ! taut ge tautherm

END DO

IF ( itaut ) THEN

  IF ( jj == jr_min ) THEN
    remsf            = r(jr_min)
    demsf            = rho(jr_min)
    gemst            = rstmss(jr_min)
    temsf            = t(jr_min)
  ELSE
    remsf            = r(jj)      + ( r(jj-1)      - r(jj)      ) * ( tautherm - tautp )/( taut - tautp + epsilon )
    demsf            = rho(jj)    + ( rho(jj-1)    - rho(jj)    ) * ( tautherm - tautp )/( taut - tautp + epsilon )
    gemst            = rstmss(jj) + ( rstmss(jj-1) - rstmss(jj) ) * ( tautherm - tautp )/( taut - tautp + epsilon )
    temsf            = t(jj)      + ( t(jj-1)      - t(jj)      ) * ( tautherm - tautp )/( taut - tautp + epsilon )
  END IF ! jj = jr_min

ELSE

  remsf             = zero
  demsf             = rho(jr_min)
  temsf             = t(jr_min)
  gemst             = zero

END IF

gemsf               = gemst/msolar

WRITE (nprint,447) nutype(n),remsf,demsf,temsf,gemsf

IF ( prnttest ) THEN
  DO i = 1,40
    IF ( neden(i) >= niedn(i) ) neden(i) = 0
  END DO
ELSE
  DO i = 1,40
    neden(i)        = nedens(i)
  END DO
END IF ! prnttest

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (psi0l, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0l     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi1l, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi1l     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psiteq, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psiteq    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi10, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi10     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (tautrns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tautrns   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (tauopt, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tauopt    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (taueanp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taueanp   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (taueanc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taueanc   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (taueahnc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taueahnc  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxn     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (tauajk, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tauajk    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (taujk, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'taujk     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (tausjk, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tausjk    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dnuradjn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dnuradjn  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (strsnujn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'strsnujn  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (edotnp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'edotnp    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (psi0eq, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi0eq    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (reanpin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reanpin   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (reanpnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reanpnet  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (reancin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reancin   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (reancnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reancnet  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (reahncin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reahncin  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (reahncnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'reahncnet '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rlocal, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rlocal    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rtr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rtr       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (tflux, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'tflux     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (vdrift, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vdrift    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ciicra, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ciicra    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rpnesin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnesin   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rpnesnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnesnet  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (qnes, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnes      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rnes, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnes      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rpnncsin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnncsin  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rpnncsnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnncsnet '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (qnncs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnncs     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rnncs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnncs     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rpnnsin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnsin   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rpnnsnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnsnet  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (qns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qns       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rns       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rpnnnsin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnnsin  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rpnnnsnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnnsnet '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (qnns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnns      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rnns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnns      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rpnnAsin, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnAsin  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rpnnAsnet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpnnAsnet '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (qnAs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'qnAs      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rnAs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rnAs      '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (wk_mag_p, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wk_mag_p  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (wk_mag_n, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wk_mag_n  '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rpain, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpain     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rpanet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rpanet    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rbain, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rbain     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (rbanet, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rbanet    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ymefrnp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymefrnp   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymenc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymenc     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymehnc, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymehnc    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnps, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnps    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnns, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnns    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnhs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnhs    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnas, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnas    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnes, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnes    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnncs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnncs   '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnn     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnnn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnnn    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnAn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnAn    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnts, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnts    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnbs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnbs    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnps_ed, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnps_ed '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ymdnns_ed, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ymdnns_ed '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ecoefp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ecoefp    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (uenulf, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uenulf    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE editn
