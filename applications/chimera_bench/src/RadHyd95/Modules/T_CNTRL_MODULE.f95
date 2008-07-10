!-----------------------------------------------------------------------
!    Module:       t_cntrl_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!-----------------------------------------------------------------------

MODULE t_cntrl_module

USE kind_module

SAVE

!-----------------------------------------------------------------------
!  Times
!-----------------------------------------------------------------------
!  time     : the elapsed time since the initiation of the calculation
!   (i.e., since ncycle = 0)
!
!  t_start  : the time at the initiation of the calculation (typically
!   t_start = 0.0).
!
!  t_bounce : the time of core bounce.
!
!  t_stop   : the elapsed time after which calculation is terminated.
!
!  tb_stop  : the elapsed time from bounce after which calculation is
!   terminated.
!-----------------------------------------------------------------------

REAL(KIND=double)                                  :: time
REAL(KIND=double)                                  :: t_start
REAL(KIND=double)                                  :: t_bounce
REAL(KIND=double)                                  :: t_stop
REAL(KIND=double)                                  :: tb_stop


!-----------------------------------------------------------------------
!  Time steps
!-----------------------------------------------------------------------
!
!  dtnph : the coordinate 'hydro' time step, set by hydro and nuclear
!   reactions (if the material is not in nse), between time cycle m and
!   time cycle m + 1.
!
!  dtnmh : the coordinate 'hydro' time step between time cycle m - 1 and
!   time cycle m.
!
!  dtj(j) : current coordinate 'hydro' time step of radial zone j.
!
!  dtjr(j) : preceding coordinate 'hydro' time step of radial zone j.
!
!  dt_sub : the coordinate time step for processes being subcycled.
!
!  dtau_nuadvct(j) : current proper neutrino advection time step of radial
!   zone j - 1/2.
!
!  dtnph_trans : the time step for absorption, emission, and scattering,
!   set by the neutrino occupation distributions and variables affected
!   by the neutrino-matter interactions, between time cycle m and time cycle
!   m + 1.
!
!  dtnmhm_trans : the coordinate time step for absorption, emission, and
!   scattering between time cycle m - 1 and time cycle m.
!
!  delt_trans : elapsed time from last implementation of source and transport
!
!  dtau_nutrns(j) : current proper time step of neutrino absorption, emission,
!   and transport of radial zone j - 1/2.
!
!  dtauj_nutrns(j) : current proper time step of neutrino source and transport
!   of radial zone j.
!
!  dtime_hydro : maximum time step set by hydro criterion
!
!  dtime_nuadvct : maximum time step set by neutrino energy advection
!   criterion
!
!  dtime_nutrns : maximum time step set by neutrino source  and transport
!   criterion
!-----------------------------------------------------------------------

REAL(KIND=double)                                  :: dtnph
REAL(KIND=double)                                  :: dtnmh
REAL(KIND=double)                                  :: dt_sub
REAL(KIND=double)                                  :: dtnph_trans
REAL(KIND=double)                                  :: dtnmhn_trans
REAL(KIND=double)                                  :: delt_trans
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: dtj
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: dtjr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: dtau_nuadvct
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: dtau_nutrns
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: dtauj_nutrns
REAL(KIND=double)                                  :: dtime_hydro
REAL(KIND=double)                                  :: dtime_nuadvct
REAL(KIND=double)                                  :: dtime_nutrns

!........c x dt.........................................................

!      cdt(j) : c * dtau_nutrns(j).
!      cdtinv(j) : 1/cdt(j).

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: cdt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: cdtinv


!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------
!  tcntrl(i) : numerical criterion for the ith time step control.
!
!     tcntrl(1)    : Courant condition - tcntrl(1) is the fraction of a
!      zone width a sonic disturbance is allowed to propagate in the
!      radial direction in one time step.
!
!     tcntrl(2)    : density time step criterion,i.e., the maximum permitted
!      abs( d(rho)/rho ) during an x-sweep.
!
!     tcntrl(3)    : hydro temperature change time step criterion, i.e.,
!      the maximum permitted abs( tl(j,i_ray) - ti(j,i_ray)/t(j,i_ray) )
!      i.e., the temperature change of radial zone j due to the hydro
!      x-sweeps.
!
!     tcntrl(4)    : Courant condition - tcntrl(4) is the fraction of a
!      zone width a sonic disturbance is allowed to propagate in the
!      angular direction in one time step.
!
!     tcntrl(5)    : density time step criterion,i.e., the maximum permitted
!      abs( d(rho)/rho ) during a y-sweep.
!
!     tcntrl(6)    : hydro temperature change time step criterion, i.e.,
!      the maximum permitted abs( tl(j,j_ray) - ti(j,j_ray)/ti(j),j_ray )
!      i.e., the temperature change of angular zone j due to the hydro
!      y-sweeps.
!
!     tcntrl(7)    : Courant condition - tcntrl(7) is the fraction of a
!      zone width a sonic disturbance is allowed to propagate in the
!      azimuthal direction in one time step.
!
!     tcntrl(8)    : density time step criterion,i.e., the maximum permitted
!      abs( d(rho)/rho ) during a z-sweep.
!
!     tcntrl(9)    : hydro temperature change time step criterion, i.e.,
!      the maximum permitted abs( tl(j,j_ray) - ti(j,j_ray)/ti(j),j_ray )
!      i.e., the temperature change of azimuthal zone j due to the hydro
!      z-sweeps.
!
!     tcntrl(10)   : maximum increase in the 'hydro' time step, i.e., the
!      maximum allowed value of dtnph/dtnmh.
!
!     tcntrl(11) : neutrino net temperature change time step criterion,
!     i.e., the maximum permitted abs( dtmpnn(j,1,i_ray)/t(j) ), where
!     dtmpnn(j,1,i_ray) is the temperature change of radial zone j due to
!     all neutrinos.
!
!     tcntrl(16) : n-neutrino net electron fraction change time step
!      criterion, i.e., the maximum permitted abs( dye(j,1,i_ray)/ye(j) ), where
!      dye(j,1,i_ray) is the change of ye in zone j due to all neutrinos.
!
!     tcntrl(20+n) : n-neutrino zero-moment change time step criterion due
!      to absorption, emission, scattering, production and transport of
!      n-neutrinos.
!
!     tcntrl(31)   : convective time step control - tcntrl(31) is the fraction
!      of a zone width a convective blob is allowed to propagate in one
!      time step.
!
!     tcntrl(33)    : nuclear burn composition change time step criterion,
!      i.e., the maximum permitted abs( dyn(j,i,i_ray)/yn(j,i) ), where
!      dyn(j,i,i_ray) is the abundance change of specie i in radial zone j.
!
!     tcntrl(34)    : convective time step control - tcntrl(9) is the fraction
!      of a zone width a convective blob is allowed to propagate in one
!      time step.
!
!     tcntrl(49)   : maximum increase in the 'neutrino transport' time step,
!      i.e., the maximum allowed value of dtnphn/dtnmhn.
!
!     tcntrl(50)   : time step is set equal to tcntrl(50) if jdt(50) = -1.
!
!  dt(i)     : time step as given by time step criterion i.
!
!  jrdt(i)   : radial zone responsible, when appropriate, for the value of dt(i).
!
!   jrdt(50) = -1: time step used in the calculation (all time step controls
!                         are bypassed).
!   jrdt(50) /= -1: maximum possible time step.
!
!  jadt(i)   : angular zone responsible, when appropriate, for the value
!   of dt(i).
!
!  jzdt(i)   : azimuthal zone responsible, when appropriate, for the value
!   of dt(i).
!-----------------------------------------------------------------------

INTEGER, DIMENSION(50)                             :: jrdt
INTEGER, DIMENSION(50)                             :: jadt
INTEGER, DIMENSION(50)                             :: jzdt

REAL(KIND=double), DIMENSION(50)                   :: tcntrl
REAL(KIND=double), DIMENSION(50)                   :: dt

!-----------------------------------------------------------------------
!  Global hydro time step option
!-----------------------------------------------------------------------
!  t_step_xyz    : hydro time-step option
!
!     t_step_xyz = ye : hydro time step set to minimum of xyz hydro
!                        no yx hydro sybcycling
!     t_step_xyz = no : hydro time step set to minimum of x:
!                        yz hydro sybcycled (sub_cy_yz = ye)
!                        yz hydro t-step set to min( yz t-step, x t-step )
!                         (sub_cy_yz = no)
!-----------------------------------------------------------------------

CHARACTER (len=2)                                  :: t_step_xyz

!-----------------------------------------------------------------------
!  Time step tolarances
!-----------------------------------------------------------------------
!  rdtmax     : dtnphn > dtnph: neutrino transport is bypassed until the
!   accumulated time since the last neutrino transport update could exceed
!   dtnphn by a specified amount in the subsequent time step. Thus, if delt
!   is the accumulated time from the last implementation of neutrino
!   transport to the beginning of the current time cycle, neutrino transport
!   is implemented in the current time cycle with dtnphn set equal to the
!   the accumulated time if
!
!    delt + rdtmax*dtnph > dtnphn .
!
!   To ensure that the accumulated time will not exceed the predicted nutrino
!    transport time step dtnphn, rdtmax must be given by
!
!   rdtmax = 1 + tcntrl(10) . 
!
!   Setting rdtmax to a negative number in the data file will instruct the
!    code to reset rdtmax to the above value.
!
!  dpsisp(n), psipmin(n): parameters used in determining the psi0 change
!   time step due to  inelastic scattering and pair production.
!
!  dpsisp(n)  : max( abs( dpsi0(j,k,n) )/( psi0(j,k,n) + psipmin(n) ) ) .
!
!  dtst1      : density change time step control is bypassed if rho(j) < dtst1.
!
!  dtst2      :
!
!  dtst3      :
!
!  ttst1      : temperature change time step control is bypassed if t(j) < ttst1.
!-----------------------------------------------------------------------

REAL(KIND=double)                                  :: rdtmax
REAL(KIND=double)                                  :: dtst1
REAL(KIND=double)                                  :: dtst2
REAL(KIND=double)                                  :: dtst3
REAL(KIND=double)                                  :: ttst1
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: psimin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: psipmin
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: dpsisp

!-----------------------------------------------------------------------
!  Individual time steps
!-----------------------------------------------------------------------
!  idtj,jtv,jtvshk,ncyshk,ncyrho,rhojtv: parameters that can be used to
!   circumvent the severe time step restriction due to the Courant condition
!   applied to zones compressed to very high densities.
!
!     If idtj = 0 or jtv = 0, each radial zone has a common 'hydro' time
!      step, dtnph, which is the minimum of all 'hydro' time step criteria
!      for all zones. (This is the usual time stepping procedure.)
!
!     If ncycle ge ncyshk, then idtj is set to 1, and radial zones jmin to
!      jtv have individual 'hydro' time steps, dtj(j). dtj(j) is the minimum
!      of the minimum 'hydro' time step criteria of zone j and that of zone
!      j+1; i.e.,
!
!          dtj(j) le dtj(j+1);
!
!      radial zones jtv+1 to jr_max have a common 'hydro' time step, dtnph,
!      which is the minimum of the 'hydro' time step restrictions applied to
!      these zones. The quantity jtv is determined by computing the location
!      of the shock (if one is present), and using the equation
!
!          jtv = jshock - jtvshk
!
!      If ncycle ge ncyrho, then idtj is set to 1, and radial zones jmin to
!       jtv have individual'hydro' time steps, dtj(j). dtj(j) is the minimum
!       of the minimum 'hydro' time step criteria of zone j and that of zone
!       j+1; i.e.,
!
!          dtj(j) le dtj(j+1);
!
!      radial zones jtv+1 to jr_max have a common 'hydro' time step, dtnph,
!       which is the minimum of the 'hydro' time step restrictions applied to
!       these zones. The quantity jtv is the maximum j such that rho(j) > rhojtv.
!
!      dth(j): minimum 'hydro' time step of radial zone j. When idtj ne 0 and
!       jtv ne 0, dtj(j) is given by
!
!          dtj(j) = min( dth(j) , dtj(j+1) ) .
!-----------------------------------------------------------------------

INTEGER                                            :: idtj
INTEGER                                            :: jtv
INTEGER                                            :: jtvshk
INTEGER                                            :: ncyrho
INTEGER                                            :: ncyshk

REAL(KIND=double)                                  :: rhojtv
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)       :: dth


END module t_cntrl_module
