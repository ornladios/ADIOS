
!-----------------------------------------------------------------------
!          ihydro = 0: hydrodynamics bypassed.
!          ihydro = 1: hydrodynamics included.
!-----------------------------------------------------------------------

ihydro                       1                                          ihydro

!-----------------------------------------------------------------------
!          irelhy = 0: Newtonian hydrodynamics.
!          irelhy = 1: PN hydrodynamics.
!          irelhy = 2: general relativistic hydrodynamics.
!-----------------------------------------------------------------------

irel                         1                                          irelhy

!-----------------------------------------------------------------------
!          ilapsehy = 0: Lapse = 1 in the hydrodynamics if irelhy = 1.
!          ilapsehy = 1: Lapse computed in the hydrodynamics if irelhy = 1.
!-----------------------------------------------------------------------

irel                         1                                          ilapsehy

!-----------------------------------------------------------------------
!          degen = degenerate - nondegenerate criterion
!-----------------------------------------------------------------------

degen                               2.00000000E+01                      degen


!-----------------------------------------------------------------------
!         Times
!
!      t_bounce: the time of core bounce.
!-----------------------------------------------------------------------

sptime                              0.00000000E+00                      t_bounce

!-----------------------------------------------------------------------
!         Individual time steps. (Not applicable for EVH1)
!
!      idtj,jtv,jtvshk,ncyshk,ncyrho,rhojtv: parameters that can be used to circumvent the severe
!       time step restriction due to the Courant condition applied to zones compressed to very high
!       densities.
!
!      If idtj = 0 or jtv = 0, each radial zone has a common 'hydro' time step, dtnph, which is the
!       minimum of all 'hydro' time step criteria for all zones. (This is the usual time stepping
!       procedure.)
!
!      If ncycle ge ncyshk, then idtj is set to 1, and radial zones jmin to jtv have individual
!       'hydro' time steps, dtj(j). dtj(j) is the minimum of the minimum 'hydro' time step criteria
!       of zone j and that of zone j+1; i.e.,
!
!          dtj(j) le dtj(j+1);
!
!       radial zones jtv+1 to jm have a common 'hydro' time step, dtnph, which is the minimum of the 
!       'hydro' time step restrictions applied to these zones. The quantity jtv is determined by
!       computing the location of the shock (if one is present), and using the equation
!
!          jtv = jshock - jtvshk
!
!      If ncycle ge ncyrho, then idtj is set to 1, and radial zones jmin to jtv have individual 
!       'hydro' time steps, dtj(j). dtj(j) is the minimum of the minimum 'hydro' time step criteria
!       of zone j and that of zone j+1; i.e.,
!
!          dtj(j) le dtj(j+1);
!
!       radial zones jtv+1 to jm have a common 'hydro' time step, dtnph, which is the minimum of the 
!       'hydro' time step restrictions applied to these zones. The quantity jtv is the maximum j such 
!       that rho(j) > rhojtv.
!-----------------------------------------------------------------------

vtime                        0                                          idtj
vtime                        0                                          jtv
vtime                        6                                          jtvshk
vtime                  9000000                                          ncyshk
vtime                  9000000                                          ncyrho
vtime                               1.00000000E+13                      rhojtv

!-----------------------------------------------------------------------
!         Time step control criteria
!
!      tcntrl(i): numerical criterion for the ith time step control.
!
!          tcntrl(1)     : Courant condition - tcntrl(1) is the fraction
!           of a zone width a sonic disturbance is allowed to propagate
!           in the radial direction in one time step.
!
!          tcntrl(2)     : density time step criterion,i.e., the maximum
!           permitted abs( d(rho)/rho ) during an x-sweep.
!
!          tcntrl(3)     : hydro temperature change time step criterion,
!           i.e., the maximum permitted abs( tl(j,i_ray) - ti(j,i_ray)/t(j,i_ray) )
!           i.e., the temperature change of radial zone j due to the hydro
!           x-sweeps.
!
!          tcntrl(4)     :  Courant condition - tcntrl(4) is the fraction
!           of a zone width a sonic disturbance is allowed to propagate
!           in the angular direction in one time step.
!
!          tcntrl(5)     : density time step criterion,i.e., the maximum
!           permitted abs( d(rho)/rho ) during a y-sweep.
!
!          tcntrl(6)     : hydro temperature change time step criterion, i.e.,
!           the maximum permitted abs( tl(j,j_ray) - ti(j,j_ray)/ti(j),j_ray )
!           i.e., the temperature change of angular zone j due to the hydro
!           y-sweeps.
!
!          tcntrl(7)     : Courant condition - tcntrl(7) is the fraction
!           of a zone width a sonic disturbance is allowed to propagate
!           in the azimuthal direction in one time step.
!
!          tcntrl(8)     : density time step criterion,i.e., the maximum
!           permitted abs( d(rho)/rho ) during a z-sweep.
!
!          tcntrl(9)     : hydro temperature change time step criterion, i.e.,
!           the maximum permitted abs( tl(j,j_ray) - ti(j,j_ray)/ti(j),j_ray )
!           i.e., the temperature change of azimuthal zone k due to the hydro
!           z-sweeps.
!
!          tcntrl(10)    : maximum increase in the 'hydro' time step, i.e.,
!           the maximum allowed value of dtnph/dtnmh.
!
!          tcntrl(31)    : convective time step control - tcntrl(7) is the
!           fraction of a zone width a convective blob is allowed to propagate
!           in one time step.
!
!          jdt(50) = -1  : time step used in the calculation (all time step
!           controls are bypassed).
!          jdt(50) ne -1 : maximum possible time step.
!
!      dtst1: density change time step control is bypassed if rho(j) < dtst1.
!
!      dtst2: (obsolete)
!
!      dtst3: (obsolete)
!
!      ttst1: temperature change time step control is bypassed if t(j) < ttst1.
!-----------------------------------------------------------------------

tcntrl                       1      5.00000000E-01                      tcntrl_hydro
tcntrl                       2      3.00000000E-02                      tcntrl_hydro
tcntrl                       3      1.00000000E-02                      tcntrl_hydro
tcntrl                       4      5.00000000E-01                      tcntrl_hydro
tcntrl                       5      3.00000000E-02                      tcntrl_hydro
tcntrl                       6      1.00000000E-02                      tcntrl_hydro
tcntrl                       7      5.00000000E-01                      tcntrl_hydro
tcntrl                       8      3.00000000E-02                      tcntrl_hydro
tcntrl                       9      1.00000000E-02                      tcntrl_hydro
tcntrl                      10      1.20000000E+00                      tcntrl_hydro
tcntrl                      31      1.00000000E-01                      tcntrl_hydro
dcntrl                              1.00000000e+06                      dtst1
dcntrl                              1.00000000e-10                      dtst2
dcntrl                              1.00000000e+17                      dtst3
dcntrl                              1.00000000e+03                      ttst1

!-----------------------------------------------------------------------
!        Inner velocity boundary conditions.
!
!      iubcjmn: inner velocity boundary codition switch.
!
!          iubcjmn = 0: inner velocity (j=1) not imposed.
!          iubcjmn = 1: inner velocity (j=1) imposed.
!
!      ubcjmn: the value of the inner velocity.
!-----------------------------------------------------------------------

ubcjm                        0      0.00000000e+00                      ubcjmn

!-----------------------------------------------------------------------
!        Outer velocity boundary conditions.
!
!      iubcjmx: outer velocity boundary codition switch.
!
!          iubcjmx = 0: outer velocity (j=jmax) not imposed.
!          iubcjmx = 1: outer velocity (j=jmax) imposed.
!
!      ubcjmx: the value of the outer velocity.
!-----------------------------------------------------------------------

ubcjm                        0      0.00000000e+00                      ubcjmx

!-----------------------------------------------------------------------
!        Intermediate ( jmin < j < jmax) velocity boundary conditions.
!
!      iubcjmn: intermediate zone velocity boundary codition switch.
!
!          iuset = 0: no velocities imposed for 2 < j < jmax.
!          iuset = 1: velocities imposed for 2 < j < jmax.
!
!          uset(j): the value of the velocity imposed for zone j.
!-----------------------------------------------------------------------

iuset                        0                                          iuset

!-----------------------------------------------------------------------
!      r_u: zero velocity cutoff
!
!          r_u >= 0.: velocity computed normally
!          r_u <  0.: velocity set to zero for r(j) > abs(r_u)
!-----------------------------------------------------------------------

r_u                                 1.00000000E+20                      r_u

!-----------------------------------------------------------------------
!      ipbnd: outer pressure boundary codition switch.
!
!          ipbnd = 0: outer pressure boundary condition not imposed.
!          ipbnd = 1: outer pressure boundary condition imposed.
!
!      pbound: the value of the outer pressure boundary condition.
!-----------------------------------------------------------------------

ipbnd                        1      0.00000000e+00                      ipbnd

!-----------------------------------------------------------------------
!         Mixing length convection controls
!
!      ilcnvct: mixing length convection switch.
!
!          ilcnvct = 0: no convection.
!          ilcnvct = 1: convection computed.
!
!      iemix: energy mixing switch.
!
!          iemix = 0: no energy mixing.
!          iemix = 1: energy mixing computed.
!
!      iyemix: ye mixing switch.
!
!          iyemix = 0: no ye mixing.
!          iyemix = 1: ye mixing computed.
!
!      ipsimix: psi0 mixing switch.
!
!          ipsimix = 0: no psi0 mixing.
!          ipsimix = 1: psi0 mixing computed.
!
!      ipcnvct: turbulent pressure switch.
!
!          ipcnvct = 0: no convective, turbulent pressure.
!          ipcnvct = 1: convective, turbulent pressure computed.
!
!      iscnvct: convective energy dissipation switch.
!
!          iscnvct = 0: no convective energy dissipation.
!          iscnvct = 1: convective energy dissipation computed.
!
!      alphlmix:  multiple of pressure scale height used to compute lmix.
!
!      adrftcv: neutrino advection parameter -
!
!          if vdriftk < adrftcv*ulcnvct
!
!       where vdriftk is the drift velocity of neutrinos of energy zone k, ulcnvct is the convective
!       velocity, and adrftcv is an arbitrary parameter of order unity, neutrinos are advected with
!       the matter.
!
!          if vdriftk > adrftcv*ulcnvct
!
!       neutrinos are not advected with the matter. In this case psi0aa(k) is set equal to psi0(jf,k,n)
!       so the the neutrinos of energy k delivered to the zone are equal to those removed, resulting
!       in no net advection.
!-----------------------------------------------------------------------

conv                         0                                          ilcnvct
conv                         0                                          iemix
conv                         0                                          iyemix
conv                         0                                          ipsimix
conv                         0                                          ipcnvct
conv                         0                                          iscnvct
conv                                1.00000000E+00                      adrftcv
conv                                5.00000000E+00                      alphlmix

!-----------------------------------------------------------------------
!      ipq: pseudoviscosity switch
!
!          ipq = 0: pseudoviscosity (pq(j)) computed normally.
!          ipq = 1: pq(j)=0 unless ua(jp) gt 0 for some jp, i.e., pq(j)=0 during infall.
!          ipq = 2: pq(j)=0.
!
!      q0(j): pseudoviscous pressure multiplyer for radial zone j.
!-----------------------------------------------------------------------

shock                        0      1.00000000e+00                      shock
