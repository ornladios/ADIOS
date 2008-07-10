SUBROUTINE transport_write_adios( ndump, nnu )
!-----------------------------------------------------------------------
!
!    File:         transport_write
!    Module:       transport_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/6/00
!
!    Purpose:
!      To dump the transport keys and control parameters.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!   ndump       : unit number to dump the data
!
!    Output arguments:
!        none
!
!    Include files:
!      cycle_module, eos_snc_x_module, it_tol_module, nu_dist_module,
!      nu_energy_grid_module, prb_cntl_module, radial_ray_module,
!      t_cntrl_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle, intnu_trns, ncynu_trns, intnur_trns
USE eos_snc_x_module, ONLY : rhopnu
USE it_tol_module, ONLY : iternu, tolnut, tolnuye, tolnupsi, tolpsimin, &
& a_prec
USE nu_dist_module, ONLY : unumn, unumx, stwt
USE nu_energy_grid_module, ONLY : nnugp, unui, dunui, unubi, nnugpmx
USE prb_cntl_module, ONLY : inutrn, idiff, ireltrns, iyenu, itnu, jnumin, &
& jnumax, iaefnp, rhoaefnp, iaence, edmpe, iaenca, edmpa, iaencnu, roaencnu, &
& iaenct, roaenct, i_aeps, icc_wm, iscat, in, ip, ihe, iheavy, iicor, &
& inc_wm, ietann, nes, rhonesmn, rhonesmx, nncs, rhonncsmn, rhonncsmx, &
& ipair, rhopairemn, rhopairemx, rhopairtmn, rhopairtmx, &
& ibrem, rhobrememn, rhobrememx, rhobremtmn, rhobremtmx, &
& isctn, rhosctnemn, rhosctnemx, rhosctntmn, rhosctntmx, &
& isctnn, rhosctnnemn, rhosctnnemx, rhosctnntmn, rhosctnntmx, &
& isctnA, rhosctnAemn, rhosctnAemx, rhosctnAtmn, rhosctnAtmx
USE radial_ray_module, ONLY : imin, imax, nse_c, nprint
USE t_cntrl_module, ONLY : rdtmax, psimin, psipmin, tcntrl

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER*8, INTENT(in)              :: ndump           ! unit number to write restart file
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! do index
INTEGER                          :: k               ! radial zone index
INTEGER                          :: n               ! neutrino flavor index

!-----------------------------------------------------------------------
!       Formats
!-----------------------------------------------------------------------

    1 FORMAT ('!                  \\\\\                ///// ')
    2 FORMAT ('!                        TRANSPORT KEYS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!-----------------------------------------------------------------------
!
!           \\\\\ TRANSPORT PARAMETERS AND SWITCHES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transport switches
!-----------------------------------------------------------------------

   21 FORMAT ('!      inutrn: transport switch.')
   22 FORMAT ('!      idiff : transport boundary conditions toggle..')
   23 FORMAT ('!      ireltrns : Newtonian - GR transport toggle.')
   24 FORMAT ('!      iyenu : transport electron fraction change switch.')
   25 FORMAT ('!      itnu(n) : transport temperature change switch.')
   26 FORMAT ('inutrn',14x,i10,42x,'inutrn')
   27 FORMAT ('idiff ',14x,i10,42x,'idiff')
   28 FORMAT ('irel  ',14x,i10,42x,'ireltrns')
   29 FORMAT ('iyenu ',14x,i10,42x,'iyenu')
   30 FORMAT ('itnu  ',4x,2i10,42x,'itnu')

!-----------------------------------------------------------------------
!  Transport radial zones
!-----------------------------------------------------------------------

   41 FORMAT ('!      jnumin : inner radial zone for which neutrino transport is computed.')
   42 FORMAT ('!      jnumax : outer radial zone for which neutrino transport is computed.')
   43 FORMAT ('jnu   ',14x,i10,42x,'jnumin')
   44 FORMAT ('jnu   ',14x,i10,42x,'jnumax')

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------

   51 FORMAT ('!      nnugp(n) : the number of energy zones for neutrinos of type n.')
   52 FORMAT ('!      unumn : the minimum value of the zone-centered neutrino energy (MeV).')
   53 FORMAT ('!      unumx : the maximum value of the zone-centered neutrino energy (MeV).')
   54 FORMAT ('!      stwt(n) : the statistical weight of neutrinos of type n.')
   55 FORMAT ('unu   ',4x,2i10,42x,'nnugp')
   56 FORMAT ('unu   ',29x,1pe15.8,22x,'unumn')
   57 FORMAT ('unu   ',29x,1pe15.8,22x,'unumx')
   58 FORMAT ('unu   ',14x,i10,5x,1pe15.8,22x,'unu')
   59 FORMAT ('unu   ',14x,i10,5x,1pe15.8,22x,'unub')
   60 FORMAT ('stwt  ',14x,i10,5x,1pe15.8,22x,'stwt')

!-----------------------------------------------------------------------
!  Emission and absorption
!-----------------------------------------------------------------------

   71 FORMAT ('!      iaefnp : emission and absorption  on free nucleons switch.')
   72 FORMAT ('!      rhoaefnp :  density above which iaefnp is turned off.')
   73 FORMAT ('!      iaence : emission and absorption of nu_e s on nuclei, FFN prescription switch.')
   74 FORMAT ('!      edmpe : daughter-parent mean excited state energy difference (MeV).')
   75 FORMAT ('!      iaenca : emission and absorption of nu_ebars s on nuclei, FFN prescription switch.')
   76 FORMAT ('!      edmpa : daughter-parent mean excited state energy difference (MeV).')
   77 FORMAT ('!      iaencnu : emission and absorption of nu-e s and nu_ebars s on nuclei, Haxtons prescription switch.')
   78 FORMAT ('!      roaencnu : density above which iaencnu is turned off.')
   79 FORMAT ('!      iaenct : emission and absorption of e-neutrinos on nuclei (Hix et al., 2003).')
   80 FORMAT ('!      roaenct : density above which Hix et al. emission and absorption of e-neutrinos on nuclei is turned off.')
   81 FORMAT ('!      i_aeps : inclusion of recoil, thermal motions, nucleon blocking factors switch.')
   82 FORMAT ('!      icc_wm: weak magnetism correction for charged current interactions switch.')
   83 FORMAT ('iaefnp',14x,i10,5x,1pe15.8,22x,'iaefnp')
   84 FORMAT ('aence ',14x,i10,5x,1pe15.8,22x,'aence ')
   85 FORMAT ('aenca ',14x,i10,5x,1pe15.8,22x,'aenca ')
   86 FORMAT ('aencnu',14x,i10,5x,1pe15.8,22x,'aencnu')
   87 FORMAT ('aenct ',14x,i10,5x,1pe15.8,22x,'aenct ')
   88 FORMAT ('i_aeps',14x,i10,42x,'i_aeps')
   89 FORMAT ('icc_wm',14x,i10,42x,'icc_wm')

!-----------------------------------------------------------------------
!  Scattering switch
!-----------------------------------------------------------------------

   91 FORMAT ('!      iscat : general scattering switch.')
   92 FORMAT ('iscat ',14x,i10,42x,'iscat')

!-----------------------------------------------------------------------
!  Isoenergetic scattering
!-----------------------------------------------------------------------

  101 FORMAT ('!      in : neutrino-neutron scattering switch.')
  102 FORMAT ('!      ip : neutrino-proton scattering switch.')
  103 FORMAT ('!      ihe : neutrino-helium isoenergetic scattering switch.')
  104 FORMAT ('!      iheavy : neutrino-heavy nucleus isoenergetic scattering switch.')
  105 FORMAT ('!      iicor : ion-ion correlation switch.')
  106 FORMAT ('!      inc_wm : weak magnetism correction for neutral current interactions switch.')
  107 FORMAT ('!      ietann : final state blocking switch.')
  108 FORMAT ('in    ',14x,i10,42x,'in')
  109 FORMAT ('ip    ',14x,i10,42x,'ip')
  110 FORMAT ('ihe   ',14x,i10,42x,'ihe')
  111 FORMAT ('iheavy',14x,i10,42x,'iheavy')
  112 FORMAT ('iicor ',14x,i10,42x,'iicor')
  113 FORMAT ('inc_wm',14x,i10,42x,'inc_wm')
  114 FORMAT ('ietann',14x,i10,42x,'ietann')

!-----------------------------------------------------------------------
!  Neutrino-electron scattering
!-----------------------------------------------------------------------

  121 FORMAT ('!      nes : neutrino-electron scattering switch.')
  122 FORMAT ('!      rhonesmn : density below which n-neutrino-electron scattering is turned off.')
  123 FORMAT ('!      rhonesmx : density above which n-neutrino-electron scattering is turned off.')
  124 FORMAT ('nes   ',4x,i10,2(5x,1pe15.8),12x,'nes')

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering (Haxton rates)
!-----------------------------------------------------------------------

  131 FORMAT ('!      nncs : neutrino-nucleus inelastic scattering switch. (Haxton rates)')
  132 FORMAT ('!      rhonncsmn : density below which n-neutrino-nucleus scattering is turned off.')
  133 FORMAT ('!      rhonncsmx : density above which n-neutrino-nucleus scattering is turned off.')
  134 FORMAT ('sncnu ',4x,i10,2(5x,1pe15.8),12x,'sncnu')

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation
!-----------------------------------------------------------------------

  141 FORMAT ('!      ipair : electron-positron pair annihilation switch.')
  142 FORMAT ('!      rhopairemn : density below which nu_e - nu_ebar ipair is turned off.')
  143 FORMAT ('!      rhopairemx : density above which nu_e - nu_ebar ipair is turned off.')
  144 FORMAT ('!      rhopairtmn : density below which nu_x - nu_xbar ipair is turned off.')
  145 FORMAT ('!      rhopairtmx : density above which nu_x - nu_xbar ipair is turned off.')
  146 FORMAT ('ipair ',4x,i10,4(1pe10.2),12x,'ipair')

!-----------------------------------------------------------------------
!  Nucleon-nucleon pair bremsstrahlung
!-----------------------------------------------------------------------

  151 FORMAT ('!      ibrem : bremsstrahlung pair annihilation switch.')
  152 FORMAT ('!      rhobrememn : density below which nu_e - nu_ebar ibrem is turned off.')
  153 FORMAT ('!      rhobrememx : density above which nu_e - nu_ebar ibrem is turned off.')
  154 FORMAT ('!      rhobremtmn : density below which nu_x - nu_xbar ibrem is turned off.')
  155 FORMAT ('!      rhobremtmx : density above which nu_x - nu_xbar ibrem is turned off.')
  156 FORMAT ('ibrem ',4x,i10,4(1pe10.2),12x,'ibrem')

!-----------------------------------------------------------------------
!  Nucleon-nucleon elastic scattering
!-----------------------------------------------------------------------

  161 FORMAT ('!      isctn : neutrino-nucleon elastic scattering switch.')
  162 FORMAT ('!      rhosctnemn : density below which nu_e - nu_ebar isctn is turned off.')
  163 FORMAT ('!      rhosctnemx : density above which nu_e - nu_ebar isctn is turned off.')
  164 FORMAT ('!      rhosctntmn : density below which nu_x - nu_xbar isctn is turned off.')
  165 FORMAT ('!      rhosctntmx : density above which nu_x - nu_xbar isctn is turned off.')
  166 FORMAT ('isctn ',4x,i10,4(1pe10.2),12x,'isctn ')

!-----------------------------------------------------------------------
!  Nucleon-nucleon inelastic scattering
!-----------------------------------------------------------------------

  171 FORMAT ('!      isctnn : neutrino-nucleon inelastic scattering switch.')
  172 FORMAT ('!      rhosctnnemn : density below which nu_e - nu_ebar isctnn is turned off.')
  173 FORMAT ('!      rhosctnnemx : density above which nu_e - nu_ebar isctnn is turned off.')
  174 FORMAT ('!      rhosctnntmn : density below which nu_x - nu_xbar isctnn is turned off.')
  175 FORMAT ('!      rhosctnntmx : density above which nu_x - nu_xbar isctnn is turned off.')
  176 FORMAT ('isctnn',4x,i10,4(1pe10.2),12x,'isctnn')

!-----------------------------------------------------------------------
!  Nucleon-nucleus inelastic scattering
!-----------------------------------------------------------------------

  181 FORMAT ('!      isctnA : neutrino-nucleus inelastic scattering switch.')
  182 FORMAT ('!      rhosctnAemn : density below which nu_e - nu_ebar isctnA is turned off.')
  183 FORMAT ('!      rhosctnAemx : density above which nu_e - nu_ebar isctnA is turned off.')
  184 FORMAT ('!      rhosctnAtmn : density below which nu_x - nu_xbar isctnA is turned off.')
  185 FORMAT ('!      rhosctnAtmx : density above which nu_x - nu_xbar isctnA is turned off.')
  186 FORMAT ('isctnA',4x,i10,4(1pe10.2),12x,'isctnA')

!-----------------------------------------------------------------------
!  Neutrino pressure
!-----------------------------------------------------------------------

  201 FORMAT ('!      rhopnu : density above which neutrino pressure is computed isotropically.')
  202 FORMAT ('rhopnu',29x,1pe15.8,22x,'rhopnu')

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------

  211 FORMAT ('!      tcntrl(10+n) : n-neutrino net temperature change time step criterion.')
  212 FORMAT ('!      tcntrl(15+n) : n-neutrino net electron fraction change time step criterion.')
  213 FORMAT ('!      tcntrl(20+n) : n-neutrino zero-moment change time step criterion due to absorption, emission, and transport.')
  214 FORMAT ('!      tcntrl(49): maximum increase in the neutrino transport time step.')
  215 FORMAT ('!      tcntrl(50): time step is set equal to tcntrl(50) if jdt(50) = -1.')
  216 FORMAT ('tcntrl',14x,i10,5x,1pe15.8,22x,'tcntrl_trans')

!-----------------------------------------------------------------------
!  Hydro subcycling
!-----------------------------------------------------------------------

  221 FORMAT ('!      intnu_trns : hydro subcycling parameter for sources and transport.')
  222 FORMAT ('!      ncynu_trns : the number of hydro updates since the last intnu_trns update.')
  223 FORMAT ('!      intnur_trns : value of intnu_trns for the preceding update.')
  224 FORMAT ('intnu ',14x,i10,42x,'intnu_trns')
  225 FORMAT ('intnu ',14x,i10,42x,'ncynu_trns')
  226 FORMAT ('intnu ',14x,i10,42x,'intnur_trns')

!-----------------------------------------------------------------------
!  Time step tolarances
!-----------------------------------------------------------------------

  231 FORMAT ('!      rdtmax : neutrino transport is implemented if delt + rdtmax*dtnph > dtnphn.')
  232 FORMAT ('!      psimin : dpsisp(n) = max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psimin(n) ) ).')
  233 FORMAT ('!      psipmin : dpsisp(n) = max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psipmin(n) ) ).')
  234 FORMAT ('tcntrl',29x,1pe15.8,22x,'rdtmax')
  235 FORMAT ('psitol',14x,i10,5x,1pe15.8,22x,'psimin')
  236 FORMAT ('psitol',14x,i10,5x,1pe15.8,22x,'psipmin')

!-----------------------------------------------------------------------
!  Neutrino transport tolerances
!-----------------------------------------------------------------------

  241 FORMAT ('!      iternu : maximum number of iterations TO SOLVE neutrino transport variables.')
  242 FORMAT ('!      tolnut: Temperature convergence parameter for neutrino transport.')
  243 FORMAT ('!      tolnuye : Electron fraction convergence parameter for neutrino transport.')
  244 FORMAT ('!      tolnupsi : psi0 convergence parameters for neutrino transport.')
  245 FORMAT ('!      tolpsimin : abs( dpsi0/max( psi0, tolpsimin ) < tolnupsi.')
  246 FORMAT ('!      a_prec: recompute neutrino rates if dabs( agr(j) - agr_prev(j) )/agr(j) > a_prec.')
  247 FORMAT ('iternu',14x,i10,42x,'iternu')
  248 FORMAT ('tolnu ',29x,1pe15.8,22x,'tolnut')
  249 FORMAT ('tolnu ',29x,1pe15.8,22x,'tolnuye')
  250 FORMAT ('tolnu ',29x,1pe15.8,22x,'tolnupsi')
  251 FORMAT ('tolnu ',29x,1pe15.8,22x,'tolpsimin')
  252 FORMAT ('a_prec',29x,1pe15.8,22x,'a_prec')

!-----------------------------------------------------------------------
!  Document the dump
!-----------------------------------------------------------------------

 1001 FORMAT (' ***Transport keys dump written at cycle  ',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!
!           \\\\\ TRANSPORT PARAMETERS AND SWITCHES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Header
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!-----------------------------------------------------------------------
!  Transport switches
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,22)
WRITE (ndump,23)
WRITE (ndump,24)
WRITE (ndump,25)
WRITE (ndump,15)
WRITE (ndump,26) inutrn
WRITE (ndump,27) idiff
WRITE (ndump,28) ireltrns
WRITE (ndump,29) iyenu
WRITE (ndump,30) (n,itnu(n),n = 1,nnu)

!-----------------------------------------------------------------------
!  Transport radial zones
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,41)
WRITE (ndump,42)
WRITE (ndump,15)
WRITE (ndump,43) jnumin
WRITE (ndump,44) jnumax

!-----------------------------------------------------------------------
!  Neutrino group energies
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,51)
WRITE (ndump,52)
WRITE (ndump,53)
WRITE (ndump,54)
WRITE (ndump,15)
WRITE (ndump,55) (n,nnugp(n),n = 1,nnu)
WRITE (ndump,56) unumn
WRITE (ndump,57) unumx
WRITE (ndump,58) (k,unui(k),k = 1,nnugpmx)
WRITE (ndump,59) (k,unubi(k),k = 1,nnugpmx+1)
WRITE (ndump,60) (n,stwt(n),n = 1,nnu)

!-----------------------------------------------------------------------
!  Emission and absorption
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,71)
WRITE (ndump,72)
WRITE (ndump,73)
WRITE (ndump,74)
WRITE (ndump,75)
WRITE (ndump,76)
WRITE (ndump,77)
WRITE (ndump,78)
WRITE (ndump,79)
WRITE (ndump,80)
WRITE (ndump,81)
WRITE (ndump,82)
WRITE (ndump,15)
WRITE (ndump,83) iaefnp,rhoaefnp
WRITE (ndump,84) iaence,edmpe
WRITE (ndump,85) iaenca,edmpa
WRITE (ndump,86) iaencnu,roaencnu
WRITE (ndump,87) iaenct,roaenct
WRITE (ndump,88) i_aeps
WRITE (ndump,89) icc_wm

!-----------------------------------------------------------------------
!  Scattering switch
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,91)
WRITE (ndump,15)
WRITE (ndump,92) iscat

!-----------------------------------------------------------------------
!  Isoenergetic scattering
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,101)
WRITE (ndump,102)
WRITE (ndump,103)
WRITE (ndump,104)
WRITE (ndump,105)
WRITE (ndump,106)
WRITE (ndump,107)
WRITE (ndump,15)
WRITE (ndump,108) in
WRITE (ndump,109) ip
WRITE (ndump,110) ihe
WRITE (ndump,111) iheavy
WRITE (ndump,112) iicor
WRITE (ndump,113) inc_wm
WRITE (ndump,114) ietann

!-----------------------------------------------------------------------
!  Neutrino-electron scattering
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,121)
WRITE (ndump,122)
WRITE (ndump,123)
WRITE (ndump,15)
WRITE (ndump,124) nes,rhonesmn,rhonesmx

!-----------------------------------------------------------------------
!  Neutrino-nucleus inelastic scattering (Haxton rates)
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,131)
WRITE (ndump,132)
WRITE (ndump,133)
WRITE (ndump,15)
WRITE (ndump,134) nncs,rhonncsmn,rhonncsmx

!-----------------------------------------------------------------------
!  Electron-positron pair annihilation
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,141)
WRITE (ndump,142)
WRITE (ndump,143)
WRITE (ndump,144)
WRITE (ndump,145)
WRITE (ndump,15)
WRITE (ndump,146) ipair,rhopairemn,rhopairemx,rhopairtmn,rhopairtmx

!-----------------------------------------------------------------------
!  Nucleon-nucleon pair bremsstrahlung
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,151)
WRITE (ndump,152)
WRITE (ndump,153)
WRITE (ndump,154)
WRITE (ndump,155)
WRITE (ndump,15)
WRITE (ndump,156) ibrem,rhobrememn,rhobrememx,rhobremtmn,rhobremtmx

!-----------------------------------------------------------------------
!  Nucleon-nucleon elastic scattering
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,161)
WRITE (ndump,162)
WRITE (ndump,163)
WRITE (ndump,164)
WRITE (ndump,165)
WRITE (ndump,15)
WRITE (ndump,166) isctn,rhosctnemn,rhosctnemx,rhosctntmn,rhosctntmx

!-----------------------------------------------------------------------
!  Nucleon-nucleon inelastic scattering
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,171)
WRITE (ndump,172)
WRITE (ndump,173)
WRITE (ndump,174)
WRITE (ndump,175)
WRITE (ndump,15)
WRITE (ndump,176) isctnn,rhosctnnemn,rhosctnnemx,rhosctnntmn,rhosctnntmx

!-----------------------------------------------------------------------
!  Nucleon-nucleus inelastic scattering
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,181)
WRITE (ndump,182)
WRITE (ndump,183)
WRITE (ndump,184)
WRITE (ndump,185)
WRITE (ndump,15)
WRITE (ndump,186) isctnA,rhosctnAemn,rhosctnAemx,rhosctnAtmn,rhosctnAtmx

!-----------------------------------------------------------------------
!  Neutrino pressure
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,201)
WRITE (ndump,15)
WRITE (ndump,202) rhopnu

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,211)
WRITE (ndump,212)
WRITE (ndump,213)
WRITE (ndump,214)
WRITE (ndump,215)
WRITE (ndump,15)
WRITE (ndump,216) (i,tcntrl(i),i = 11,24)
WRITE (ndump,216) (i,tcntrl(i),i = 49,50)

!-----------------------------------------------------------------------
!  Hydro subcycling
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,221)
WRITE (ndump,222)
WRITE (ndump,223)
WRITE (ndump,15)
WRITE (ndump,224) intnu_trns
WRITE (ndump,225) ncynu_trns
WRITE (ndump,226) intnur_trns

!-----------------------------------------------------------------------
!  Time step tolarances
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,231)
WRITE (ndump,232)
WRITE (ndump,233)
WRITE (ndump,15)
WRITE (ndump,234) rdtmax
WRITE (ndump,235) (n,psimin(n),n = 1,nnu)
WRITE (ndump,236) (n,psipmin(n),n = 1,nnu)

!-----------------------------------------------------------------------
!  Neutrino transport tolerances
!-----------------------------------------------------------------------

WRITE (ndump,13)
WRITE (ndump,241)
WRITE (ndump,242)
WRITE (ndump,243)
WRITE (ndump,244)
WRITE (ndump,245)
WRITE (ndump,246)
WRITE (ndump,15)
WRITE (ndump,247) iternu
WRITE (ndump,248) tolnut
WRITE (ndump,249) tolnuye
WRITE (ndump,250) tolnupsi
WRITE (ndump,251) tolpsimin
WRITE (ndump,252) a_prec

!-----------------------------------------------------------------------
!  Document the dump
!-----------------------------------------------------------------------

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE transport_write_adios
