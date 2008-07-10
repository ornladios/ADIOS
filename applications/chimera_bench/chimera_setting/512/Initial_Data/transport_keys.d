!-----------------------------------------------------------------------
!      inutrn: Transport switch.
!
!          inutrn = 0: neutrino transport bypassed.
!          inutrn = 1: neutrino transport included.
!-----------------------------------------------------------------------

inutrn                       1                                          inutrn

!-----------------------------------------------------------------------
!      idiff: Transport boundary conditions toggle.
!
!          idiff = 0: dc(j,k,n) computed normally at each iteration.
!          idiff = 1: dc(j,k,n) is the average of the current dc and the preceding timestep dc.
!          idiff = 2: dc(j,k,n) = 0 for j = jnumax, dc(j,k,n) computed normally at each iteration for 
!           j ne jnumax.
!          idiff = 3: dc(j,k,n) = 0 for j = jnumax, dc(j,k,n) computed normally only at the first
!           iteration for j ne jnumax.
!          idiff = 4: dc(j,k,n) = 0 for all j.
!          idiff = 5: dc(j,k,n) computed normally only at the first iteration.
!-----------------------------------------------------------------------

idiff                        5                                          idiff

!-----------------------------------------------------------------------
!      ireltrns: Newtonian - GR transport toggle.
!
!          ireltrns = 0: Newtonian neutrino transport.
!          ireltrns = 1: general relativistic neutrino transport.
!-----------------------------------------------------------------------

irel                         1                                          ireltrns

!-----------------------------------------------------------------------
!      iyenu: transport electron fraction change switch.
!
!          iyenu = 0: electron fraction change due to n-type neutrinos bypassed.
!          iyenu = 1: electron fraction change due to n-type neutrinos included.
!-----------------------------------------------------------------------

iyenu                        1                                          iyenu

!-----------------------------------------------------------------------
!      itnu(n): transport temperature change switch.
!
!          itnu(n) = 0: temperature change due to n-type neutrinos bypassed.
!          itnu(n) = 1: temperature change due to n-type neutrinos included.
!-----------------------------------------------------------------------

itnu               1         1                                          itnu
itnu               2         1                                          itnu
itnu               3         1                                          itnu
itnu               4         1                                          itnu

!-----------------------------------------------------------------------
!      jnumin: inner radial zone for which neutrino transport is computed.
!      jnumin < 2  ==> jnumin = 2
!-----------------------------------------------------------------------

jnu                          2                                          jnumin

!-----------------------------------------------------------------------
!      jnumax: outer radial zone for which neutrino transport is computed.
!      jnumax > jmax  ==> jnumax = jmax
!-----------------------------------------------------------------------

jnu                    9000000                                          jnumax

!-----------------------------------------------------------------------
!         Neutrino distribution parameters
!
!      nnugp(n): the number of energy zones for neutrinos of type n.
!
!      unumn: the minimum value of the zone-centered neutrino energy (MeV). 
!
!      unumx: the maximum value of the zone-centered neutrino energy (MeV).
!
!      stwt(n) is the statistical weight of neutrinos of type n, i.e., the number of distinct multiplicities of
!       neutrinos subsumed as neutrinos of type n. (Typically, stwt(1) = 1., stwt(2) = 1., and stwt(3) = 4.)
!-----------------------------------------------------------------------

unu                1        20                                          nnugp
unu                2        20                                          nnugp
unu                3        20                                          nnugp
unu                4        20                                          nnugp
unu                                 4.00000000e+00                      unumn
unu                                 4.00000000e+02                      unumx

stwt                         1      1.00000000e+00                      stwt
stwt                         2      1.00000000e+00                      stwt
stwt                         3      2.00000000e+00                      stwt
stwt                         4      2.00000000e+00                      stwt

!-----------------------------------------------------------------------
!      iaefnp: emission and absorption of e-neutrinos and e-antineutrinos on free neutrons and
!       protons switch.
!
!          iaefnp = 0: emission and absorption of e-neutrinos and e-antineutrinos on free neutrons
!           and protons turned off.
!          iaefnp = 1: emission and absorption of e-neutrinos and e-antineutrinos on free neutrons
!           and protons turned on.
!
!      rhoaefnp: density above which emission and absorption of e-neutrinos and e-antineutrinos on 
!       free neutrons and protons is turned off.
!-----------------------------------------------------------------------

iaefnp                       1      1.00000000e+15                      iaefnp

!-----------------------------------------------------------------------
!      iaence: emission and absorption of e-neutrinos on nuclei as given by the FFN prescription 
!       switch
!
!          iaence = 0: emission and absorption of e-neutrinos on nuclei as given by the FFN 
!           prescription turned off.
!          iaence = 1: emission and absorption of e-neutrinos on nuclei as given by the FFN 
!           prescription turned on.
!
!      edmpe: difference between the mean excited state energy of the daughter nucleus and that of
!       the parent nucleus (MeV).
!-----------------------------------------------------------------------

aence                        0      3.00000000e+00                      aence

!-----------------------------------------------------------------------
!      iaenca: emission and absorption of e-antineutrinos on nuclei as given by the FFN prescription 
!       switch
!
!          iaenca = 0: emission and absorption of e-antineutrinos on nuclei as given by the FFN
!           prescription turned off.
!          iaenca = 1: emission and absorption of e-antineutrinos on nuclei as given by the FFN
!           prescription turned on.
!
!      edmpa: difference between the mean excited state energy of the daughter nucleus and that of 
!       the parent nucleus (MeV).
!-----------------------------------------------------------------------

aenca                        0      3.00000000e+00                      aenca

!-----------------------------------------------------------------------
!      iaencnu: emission and absorption of e-neutrinos and e-antineutrinos on nuclei as given by
!        Haxton's prescription
!
!          iaencnu = 0: emission and absorption of e-neutrinos and e-antineutrinos on nuclei as
!           given by Haxton's prescription turned off.
!          iaencnu = 1: emission and absorption of e-neutrinos and e-antineutrinos on nuclei as
!           given by Haxton's prescription turned on.
!
!      roaencnu: density above which emission and absorption of e-neutrinos and e-antineutrinos on
!       nuclei as given by Haxton is turned off.
!-----------------------------------------------------------------------

aencnu                       0      1.00000000e+15                      aencnu

!-----------------------------------------------------------------------
!      iaenct: emission and absorption of e-neutrinos on nuclei as given by
!       the NSE-folded tabular data for Hix et al (2003)                
!
!          iaenct  = 0: emission and absorption of e-neutrinos on nuclei as
!           given by the tabular data is turned off.                    
!          iaenct  = 1: emission and absorption of e-neutrinos on nuclei as
!           given by the tabular data is turned on.                     
!
!      roaenct: density above which emission and absorption of e-neutrinos on
!       nuclei as given by the is turned off.
!-----------------------------------------------------------------------

iaenct                       1      1.00000000e+13                      iaenct

!-----------------------------------------------------------------------
!      i_aeps : emission and aborption of e-neutrinos and e-antineutrinos
!       on free neutrons and protons phase space switch.
!
!        i_aeps = 0: emission and absorption of e-neutrinos and
!         e-antineutrinos on free neutrons and protons are computed
!         assuming no recoil, or thermal motions, and approximate
!         nucleon phase space blocking.
!        iaefnp = 1: emission and absorption of e-neutrinos and
!         e-antineutrinos on free neutrons and protons are computed
!         taking into account recoil, thermal motions, and nucleon
!         phase space blocking.
!-----------------------------------------------------------------------

i_aeps                       0                                          i_aeps

!-----------------------------------------------------------------------
!      icc_wm: charged current weak magnetism switch.
!
!          icc_wm = 0: weak magnetism correction for charged current
!           interactions turned off.
!          icc_wm = 1: weak magnetism correction for charged current
!          interactions turned on.
!-----------------------------------------------------------------------

icc_wm                       1                                          icc_wm

!-----------------------------------------------------------------------
!      iscat: general scattering switch.
!
!          iscat = 0: all scattering processes turned off.
!          iscat = 1: any scattering process is allowed to proceed if its key is turned on.
!-----------------------------------------------------------------------

iscat                        1                                          iscat

!-----------------------------------------------------------------------
!      in: neutrino-neutron scattering switch
!
!          in = 0: neutrino-neutron isoenergetic scattering turned off.
!          in = 1: neutrino-neutron isoenergetic scattering turned on.
!-----------------------------------------------------------------------

in                           1                                          in

!-----------------------------------------------------------------------
!      ip: neutrino-proton scattering switch
!
!          ip = 0: neutrino-proton isoenergetic scattering turned off.
!          ip = 1: neutrino-proton isoenergetic scattering turned on.
!-----------------------------------------------------------------------

ip                           1                                          ip

!-----------------------------------------------------------------------
!      ihe: neutrino-helium isoenergetic scattering switch
!
!          ihe = 0: neutrino-helium isoenergetic scattering turned off.
!          ihe = 1: neutrino-helium isoenergetic scattering turned on.
!-----------------------------------------------------------------------

ihe                          1                                          ihe

!-----------------------------------------------------------------------
!      iheavy: neutrino-heavy nucleus isoenergetic scattering switch
!
!          iheavy = 0: neutrino-heavy nucleus isoenergetic scattering turned off.
!          iheavy = 1: neutrino-heavy nucleus isoenergetic scattering turned on.
!-----------------------------------------------------------------------

iheavy                       1                                          iheavy

!-----------------------------------------------------------------------
!      iicor: ion-ion correlation switch.
!
!          iicor = 0: neutrino-heavy nucleus isoenergetic scattering correction due to ion-ion
!           correlation effects turned off.
!          iicor = 1: neutrino-heavy nucleus isoenergetic scattering correction due to ion-ion
!           correlation effects turned on.
!-----------------------------------------------------------------------

iicor                        1                                          iicor

!-----------------------------------------------------------------------
!      inc_wm: neutral current weak magnetism switch.
!
!          inc_wm = 0: weak magnetism correction for neutral current
!           interactions turned off.
!          inc_wm = 1: weak magnetism correction for neutral current
!          interactions turned on.
!-----------------------------------------------------------------------

inc_wm                       1                                          inc_wm

!-----------------------------------------------------------------------
!      ietann: final state blocking switch
!
!      ietann = 0: final state blocking due to neutron or proton final state occupancy in
!       neutrino-neutron and neutrino-proton isoenergetic scattering turned off.
!      ietann = 1: final state blocking due to neutron or proton final state occupancy in
!       neutrino-neutron and neutrino-proton isoenergetic scattering turned on.
!-----------------------------------------------------------------------

ietann                       1                                          ietann

!-----------------------------------------------------------------------
!      nes: neutrino-electron scattering switch.
!
!          nes = 0: neutrino-electron scattering turned off.
!          nes = 1: neutrino-electron scattering turned on.
!
!      rhonesmn: density below which n-neutrino-electron scattering is turned off.
!
!      rhonesmx: density above which n-neutrino-electron scattering is turned off.
!-----------------------------------------------------------------------

nes                1      1.00000000e+07      1.00000000E+15            nes

!-----------------------------------------------------------------------
!      nncs: neutrino-nucleus inelastic scattering switch. (Haxton rates)
!
!          nncs = 0: neutrino-nucleus inelastic scattering turned off.
!          nncs = 1: neutrino-nucleus inelastic scattering turned on.
!
!      rhonncsmn: density below which n-neutrino-electron scattering is turned off.
!
!      rhonncsmx: density above which n-neutrino-electron scattering is turned off.
!-----------------------------------------------------------------------

sncnu              0      1.00000000e+04      1.20000000E+13            sncnu

!-----------------------------------------------------------------------
!      ipair: electron-positron pair annihilation switch.
!
!          ipair = 0: electron-positron pair annihilation process turned off.
!          ipair = 1: electron-positron pair annihilation process turned on.
!
!      rhopairemn: density below which e-neutrino - e-antineutrino pair production by
!       electron-positron pair annihilation is turned off.
!
!      rhopairemx: density above which e-neutrino - e-antineutrino pair production by
!       electron-positron pair annihilation is turned off.
!
!      rhopairtmn: density below which x-neutrino - x-antineutrino pair production by
!       electron-positron pair annihilation is turned off..
!
!      rhopairtmx: density above which x-neutrino - x-antineutrino pair production by
!       electron-positron pair annihilation is turned off.
!-----------------------------------------------------------------------

ipair              1  1.00e+07  1.00e+15  1.00e+07  1.00e+15            ipair

!-----------------------------------------------------------------------
!      ibrem: bremsstrahlung pair annihilation switch.
!
!          ibrem = 0: bremsstrahlung pair annihilation process turned off.
!          ibrem = 1: bremsstrahlung pair annihilation process turned on.
!
!      rhobrememn: density below which e-neutrino - e-antineutrino pair production by
!       bremsstrahlung pair annihilation is turned off.
!
!      rhobrememx: density above which e-neutrino - e-antineutrino pair production by
!       bremsstrahlung pair annihilation is turned off.
!
!      rhobremtmn: density below which x-neutrino - x-antineutrino pair production by
!       bremsstrahlung pair annihilation is turned off..
!
!      rhobremtmx: density above which x-neutrino - x-antineutrino pair production by
!       bremsstrahlung pair annihilation is turned off.
!-----------------------------------------------------------------------

ibrem              1  1.00e+07  3.00e+11  1.00e+07  1.00e+15            ibrem

!-----------------------------------------------------------------------
!      isctn: neutrino-nucleon inelastic scattering switch.
!
!          isctn = 0: neutrino-nucleon inelastic scattering turned off.
!          isctn = 1: neutrino-nucleon inelastic scattering turned on.
!
!      rhosctnemn: density below which e-neutrino-neutrino-nucleon, e-antineutrino-neutrino-nucleon 
!       elastic scattering is turned off.
!
!      rhosctnemx: density above which e-neutrino-neutrino-nucleon, e-antineutrino-neutrino-nucleon
!       elastic scattering is turned off.
!
!      rhosctntmn: density below which x-neutrino-neutrino-nucleon, x-antineutrino-neutrino-nucleon 
!       elastic scattering is turned off.
!
!      rhosctntmx: density above which x-neutrino-neutrino-nucleon, x-antineutrino-neutrino-nucleon
!       elastic scattering is turned off.
!-----------------------------------------------------------------------

isctn              1  1.00e+07  1.00e+15  1.00e+07  1.00e+15            isctn

!-----------------------------------------------------------------------
!      isctnn: neutrino-nucleon inelastic scattering switch.
!
!          isctnn = 0: neutrino-nucleon inelastic scattering turned off.
!          isctnn = 1: neutrino-nucleon inelastic scattering turned on.
!
!      rhosctnnemn: density below which e-neutrino-neutrino-nucleon, e-antineutrino-neutrino-nucleon 
!       inelastic scattering is turned off.
!
!      rhosctnnemx: density above which e-neutrino-neutrino-nucleon, e-antineutrino-neutrino-nucleon
!       inelastic scattering is turned off.
!
!      rhosctnntmn: density below which x-neutrino-neutrino-nucleon, x-antineutrino-neutrino-nucleon 
!       inelastic scattering is turned off.
!
!      rhosctnntmx: density above which x-neutrino-neutrino-nucleon, x-antineutrino-neutrino-nucleon
!       scattering is turned off.
!-----------------------------------------------------------------------

isctnn             0  1.00e+07  1.00e+15  1.00e+07  1.00e+15            isctnn

!-----------------------------------------------------------------------
!      rhopnu: density above which neutrino pressure is computed isotropically as a equilibrated
!       relativistic gas.
!-----------------------------------------------------------------------

rhopnu                              1.00000000e+16                      rhopnu

!-----------------------------------------------------------------------
!         Time step control criteria
!
!      tcntrl(i): numerical criterion for the ith time step control.
!
!          tcntrl(10+n): n-neutrino net temperature change time step criterion,
!           i.e., the maximum permitted
!               abs( tf(j) - ti(j) )/ti(j),
!           due to n-neutrinos. (n=1,2,3,4).
!
!          tcntrl(15+n): n-neutrino net electron fraction change time step
!           criterion, i.e., the maximum permitted
!               abs( yef(j,n,i_ray) - yei(j,n,i_ray) )/yei(j) ),
!           where dye(j,1,i_ray) is the change of ye in zone j due to all neutrinos.
!
!          tcntrl(20+n): n-neutrino zero-moment change time step criterion due
!           to absorption, emission, scattering, production and transport of
!           n-neutrinos, i.e., the maximum permitted
!               abs( psi0f(j,k,n) - psi0i(j,k,n) )/max( psi0f(j,k,n), tolpsimin )
!
!          tcntrl(49): maximum increase in the 'neutrino transport' time step,
!           i.e., the maximum allowed value of dtnphn/dtnmhn.
!
!          tcntrl(50): time step is set equal to tcntrl(50) if jdt(50) = -1.
!
!      dt(i): time step as given by time step criterion i.
!
!      jdt(i): radial zone responsible, when appropriate, for the value of dt(i).
!
!          jdt(50) = -1: time step used in the calculation (all time step controls are bypassed).
!          jdt(50) ne -1: maximum possible time step.
!-----------------------------------------------------------------------

tcntrl                      11      1.00000000E-02                      tcntrl_trans
tcntrl                      12      1.00000000E-02                      tcntrl_trans
tcntrl                      13      1.00000000E-02                      tcntrl_trans
tcntrl                      14      1.00000000E-02                      tcntrl_trans
tcntrl                      15      1.00000000E-02                      tcntrl_trans
tcntrl                      16      1.00000000E-02                      tcntrl_trans
tcntrl                      17      1.00000000E-02                      tcntrl_trans
tcntrl                      18      1.00000000E-02                      tcntrl_trans
tcntrl                      19      1.00000000E-02                      tcntrl_trans
tcntrl                      20      1.00000000E-02                      tcntrl_trans
tcntrl                      21      1.00000000E-02                      tcntrl_trans
tcntrl                      22      1.00000000E-02                      tcntrl_trans
tcntrl                      23      1.00000000E-02                      tcntrl_trans
tcntrl                      24      1.00000000E-02                      tcntrl_trans
tcntrl                      49      1.20000000E+00                      tcntrl_trans
tcntrl                      50      1.00000000E+02                      tcntrl_trans

!-----------------------------------------------------------------------
!         Hydro subcycling.
!
!      intnur_trns : hydro subcycling parameters
!
!          > 0: intnu_i is the number of 'hydro' cycles per neutrino process i cycle at the current
!           time.
!         < 0: separate timesteps dtnph and dtnphn_i are computed for 'hydro' and neutrino process i,
!           respectively.
!
!         dtnphn_i < dtnph: dtnph is set equal to dtnphn_i and neutrino process i is implemented.
!
!         dtnpnn_i > dtnph: neutrino process i is bypassed until the accumulated time, delt_i, since
!          the last update for neutrino process i plus the current 'hydro' time step could exceed
!          dtnphn_i in the next cycle, i.e., until delt_i + dtnph*rdtmax > dtnphn_i. When the above
!          inequality is satisfied, neutrino process i is called with dtnphn_i = delt_i.
!
!      intnur_trns : hydro subcycling parameters for the preceding time step. i.e.,
!       the number of hydro updates per neutrino process i update at the last neutrino transport
!       update.
!
!      ncynu_trns : the number of hydro updates since the last neutrino process i update.
!-----------------------------------------------------------------------

intnu                       10                                          intnu_trns
intnu                        0                                          ncynu_trns
intnu                        1                                          intnur_trns

!-----------------------------------------------------------------------
!         Time step tolarances.
!
!      rdtmax: dtnphn > dtnph: neutrino transport is bypassed until the accumulated time since the
!       last neutrino transport update could exceed dtnphn by a specified amount in the subsequent
!       time step. Thus, if delt is the accumulated time from the last implementation of neutrino
!       transport to the beginning of the current time cycle, neutrino transport is implemented in
!       the current time cycle with dtnphn set equal to the the accumulated time if
!
!          delt + rdtmax*dtnph > dtnphn .
!
!      To ensure that the accumulated time will not exceed the predicted nutrino transport time step
!       dtnphn, rdtmax must be given by
!
!          rdtmax = 1 + tcntrl(6) . 
!
!      Setting rdtmax to a negative number in the data file will instruct the code to reset rdtmax
!       to the above value.
!
!      dpsaet(n), psimin(n): parameters used in determining the psi0 change time step due to
!       absorption, emission and transport.
!
!          dpsaet(n) = max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psimin(n) ) ) .
!
!      dpsisp(n), psipmin(n): parameters used in determining the psi0 change time step due to
!       inelastic scattering and pair production.
!
!      dpsisp(n) = max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psipmin(n) ) ) .
!-----------------------------------------------------------------------

tcntrl                              1.50000000E+00                      rdtmax
psitol                       1      1.00000000e-01                      psimin
psitol                       2      1.00000000e-01                      psimin
psitol                       3      1.00000000e-01                      psimin
psitol                       4      1.00000000e-01                      psimin
psitol                       1      1.00000000e-01                      psipmin
psitol                       2      1.00000000e-01                      psipmin
psitol                       3      1.00000000e-01                      psipmin
psitol                       4      1.00000000e-01                      psipmin

!-----------------------------------------------------------------------
!         Neutrino transport tolerances.
!
!      iternu: maximum number of iterations to attempt in order to obtain convergence of the neutrino 
!       transport variables (i.e., t, ye, and the psi0's). If iternu = 1, the variables are assumed
!       to have converged after the first iteration attempt.                                          c
!
!      nutdiv: (Obsolete)
!
!      itfail: iteration failure switch. (Obsolete)
!
!          itfail: 0 - iteration failure does not stop calculation.
!          itfail: 1 - iteration failure stops calculation.
!
!      tolnut: Temperature convergence parameter for neutrino transport. The criterion for 
!       temperature convergence is that
!
!          abs(dt/t) < tolnut .
!
!      tolnuye: Electron fraction convergence parameter for neutrino transport. The criterion for
!       electron fraction convergence is that
!
!               abs( dye/ye ) < tolnuye .
!
!      tolnupsi, tolpsimin: zero-monent neutrino occupation distribution convergence parameters for
!       neutrino transport. The criteria for neutrino occupation convergence is that
!
!               abs( dpsi0/max( psi0, tolpsimin ) < tolnupsi .
!
!      a_prec: recompute neutrino rates if
!
!               dabs( agr(j) - agr_prev(j) )/agr(j) > a_prec (Obsolete)
!
!      dyealrm:  (Obsolete)
!      dtalrm:   (Obsolete)
!      dpsialrm: (Obsolete)
!-----------------------------------------------------------------------

iternu                      30                                          iternu
iternu                       1                                          nutdiv
iternu                       0                                          itfail

tolnu                               1.00000000e-06                      tolnut
tolnu                               1.00000000e-06                      tolnuye
tolnu                               1.00000000e-06                      tolnupsi
tolnu                               1.00000000e+00                      tolpsimin

a_prec                              1.00000000e-02                      a_prec

nualrm                              1.00000000e+01                      dyealrm
nualrm                              1.00000000e+01                      dtalrm
nualrm                              1.00000000e+01                      dpsialrm

!-----------------------------------------------------------------------
!         Diagnostic control parameters. (Obsolete)
!-----------------------------------------------------------------------

diagnu                 9000000                                          nktw
diagnu                       1                                          ktt
diagnu                      18                                          kimin
diagnu                      19                                          kimax
diagnu                       1                                          kwmin
diagnu                      20                                          kwmax
diagnu                      18                                          ktmin
diagnu                      19                                          ktmax
