SUBROUTINE read_pack_hydro_keys( nreadp, nprint, iskip, nx, i_hydro_data, &
& d_hydro_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_hydro_keys
!    Module:       read_pack_hydro_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/03/03
!
!    Purpose:
!      To unpack the hydro key arrays and restore the values
!       to the appropriate variables in the appropriate modules.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp       : unit number from which to read
!  nprint       : unit number from which to print
!  iskip        : echo data read flag
!  nx           : radial array dimension
!  nrst         : cycle number at start or restart
!
!    Output arguments:
!  i_hydro_data : integer array of hydro keys
!  d_hydro_data : real*8 array of hydro keys
!
!    Include files:
!  kind_module, numerical_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp       ! unit number to read from
INTEGER, INTENT(in)              :: nprint       ! unit number to print to
INTEGER, INTENT(in)              :: iskip        ! echo data read flag
INTEGER, INTENT(in)              :: nx           ! radial array dimension
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(30)                :: i_hydro_data  ! integer array of transport keys

REAL(KIND=double), INTENT(out), DIMENSION((30+nx)) :: d_hydro_data  ! 64 bit real array of transport keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line
CHARACTER (len=10)               :: var_name

INTEGER                          :: i             ! radial do index
INTEGER                          :: j             ! radial zone index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

INTEGER                          :: int           ! integer data variable to read in an interger datum

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!  Hydrodynamics controls
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
!  degen  : Degeneracy - nondegeneracy criterion.
!-----------------------------------------------------------------------

INTEGER                                        :: ihydro
INTEGER                                        :: irelhy
INTEGER                                        :: ilapsehy

REAL(KIND=double)                              :: degen

!-----------------------------------------------------------------------
!  Outer boundary
!-----------------------------------------------------------------------
!  ipbnd   : outer pressure boundary codition switch.
!
!     ipbnd = 0 : outer pressure boundary condition not imposed.
!     ipbnd = 1 : outer pressure boundary condition imposed.
!
!  pbound  : the value of the outer pressure boundary condition.
!
!  iubcjmx : outer velocity boundary codition switch.
!
!     iubcjmx = 0 : outer velocity (j=jr_max) not imposed.
!     iubcjmx = 1 : outer velocity (j=jr_max) imposed.
!
!  ubcjmx  : the value of the outer velocity.
!-----------------------------------------------------------------------

INTEGER                                        :: ipbnd
INTEGER                                        :: iubcjmx

REAL(KIND=double)                              :: pbound
REAL(KIND=double)                              :: ubcjmx

!-----------------------------------------------------------------------
!  Inner boundary
!-----------------------------------------------------------------------
!  iubcjmn : inner velocity boundary codition switch.
!
!     iubcjmn = 0 : inner velocity (j=1) not imposed.
!     iubcjmn = 1 : inner velocity (j=1) imposed.
!
!  ubcjmn  : the value of the inner velocity.
!-----------------------------------------------------------------------

INTEGER                                        :: iubcjmn
REAL(KIND=double)                              :: ubcjmn


!-----------------------------------------------------------------------
!  Intermediate zones
!-----------------------------------------------------------------------
!  iubcjmn : intermediate zone velocity boundary codition switch.
!
!     iuset = 0 : no velocities imposed for 1 < j < jr_max.
!     iuset = 1 : velocities imposed for 1 < j < jr_max.
!
!  uset(j) : the value of the velocity imposed for zone j.
!-----------------------------------------------------------------------

INTEGER                                        :: iuset
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: uset

!-----------------------------------------------------------------------
!  Zero velocity cutoff
!-----------------------------------------------------------------------
!  r_u : zero velocity cutoff
!
!     r_u >= 0. : velocity computed normally
!     r_u <  0. : velocity set to zero for r(j) > abs(r_u)
!-----------------------------------------------------------------------

REAL(KIND=double)                              :: r_u

!-----------------------------------------------------------------------
!
!                          \\\\\ TIMES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Individual time steps
!-----------------------------------------------------------------------
!  idtj,jtv,jtvshk,ncyshk,ncyrho,rhojtv : parameters that can be used to circumvent the severe
!   time step restriction due to the Courant condition applied to zones compressed to very high
!   densities.
!
!  If idtj = 0 or jtv = 0, each radial zone has a common 'hydro' time step, dtnph, which is the
!   minimum of all 'hydro' time step criteria for all zones. (This is the usual time stepping
!   procedure.)
!
!  If ncycle ge ncyshk, then idtj is set to 1, and radial zones jr_min to jtv have individual
!   'hydro' time steps, dtj(j). dtj(j) is the minimum of the minimum 'hydro' time step criteria
!   of zone j and that of zone j+1; i.e.,
!
!     dtj(j) <= dtj(j+1);
!
!  Rradial zones jtv+1 to jr_max have a common 'hydro' time step, dtnph, which is the minimum of the 
!   'hydro' time step restrictions applied to these zones. The quantity jtv is determined by
!   computing the location of the shock (if one is present), and using the equation
!
!     jtv = jshock - jtvshk
!
!  If ncycle ge ncyrho, then idtj is set to 1, and radial zones jr_min to jtv have individual 
!   'hydro' time steps, dtj(j). dtj(j) is the minimum of the minimum 'hydro' time step criteria
!   of zone j and that of zone j+1; i.e.,
!
!     dtj(j) <= dtj(j+1);
!
!  Radial zones jtv+1 to jr_max have a common 'hydro' time step, dtnph, which is the minimum of the 
!   'hydro' time step restrictions applied to these zones. The quantity jtv is the maximum j such 
!   that rho(j) > rhojtv.
!
!  dth(j) : minimum 'hydro' time step of radial zone j. When idtj ne 0 and jtv ne 0, dtj(j) is
!   given by
!
!     dtj(j) = min( dth(j) , dtj(j+1) ) .
!-----------------------------------------------------------------------

INTEGER                                        :: idtj
INTEGER                                        :: jtv
INTEGER                                        :: jtvshk
INTEGER                                        :: ncyrho
INTEGER                                        :: ncyshk

REAL(KIND=double)                              :: rhojtv

!-----------------------------------------------------------------------
!  Time step control criteria
!-----------------------------------------------------------------------
!     tcntrl(i): numerical criterion for the ith time step control.
!
!  tcntrl(1) : Courant condition - tcntrl(1) is the fraction of a zone
!   width a sonic disturbance is allowed to propagate in one time step.
!
!  tcntrl(2) : density time step criterion,i.e., the maximum permitted
!   abs( d(rho)/rho ).
!
!  tcntrl(3) : hydro temperature change time step criterion, i.e., the
!   maximum permitted abs( dtmpmn(j,1,i_ray)/t(j) ), where dtmpmn(j,1,i_ray)
!   is the hydro temperature change of radial zone j.
!
!  tcntrl(9) : convective time step control - tcntrl(9) is the fraction
!   of a zone width a convective blob is allowed to propagate in one time step.
!
!  tcntrl(10) : maximum increase in the 'hydro' time step, i.e., the
!   maximum allowed value of dtnph/dtnmh.
!
!  dtst1 : density change time step control is bypassed if rho(j) < dtst1.
!
!  dtst2 :
!
!  dtst3 :
!
!  ttst1 : temperature change time step control is bypassed if t(j) < ttst1.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(50)               :: tcntrl

REAL(KIND=double)                              :: dtst1
REAL(KIND=double)                              :: dtst2
REAL(KIND=double)                              :: dtst3
REAL(KIND=double)                              :: ttst1

!-----------------------------------------------------------------------
!  Pseudoviscous shock parameters
!-----------------------------------------------------------------------
!  ipq     : pseudoviscosity switch
!
!     ipq = 0 : pseudoviscosity (pq_x(j)) computed normally.
!     ipq = 1 : pq_x(j)=0 unless ua(jp) gt 0 for some jp, i.e., pq_x(j)=0 during infall.
!     ipq = 2 : pq_x(j)=0.
!
!  q0_x(1) : pseudoviscous pressure multiplyer.
!-----------------------------------------------------------------------

INTEGER                                        :: ipq

REAL(KIND=double)                              :: q0_x

!-----------------------------------------------------------------------
!  Convection controls
!-----------------------------------------------------------------------
!  ilcnvct  : mixing length convection switch.
!
!     ilcnvct = 0 : no convection.
!     ilcnvct = 1 : convection computed.
!
!  iemix    : energy mixing switch.
!
!     iemix = 0: no energy mixing.
!     iemix = 1: energy mixing computed.
!
!  iyemix   : ye mixing switch.
!
!     iyemix = 0: no ye mixing.
!     iyemix = 1: ye mixing computed.
!
!  ipsimix  : psi0 mixing switch.
!
!     ipsimix = 0: no psi0 mixing.
!     ipsimix = 1: psi0 mixing computed.
!
!  ipcnvct  : turbulent pressure switch.
!
!     ipcnvct = 0: no convective, turbulent pressure.
!     ipcnvct = 1: convective, turbulent pressure computed.
!
!  iscnvct  : convective energy dissipation switch.
!
!     iscnvct = 0: no convective energy dissipation.
!     iscnvct = 1: convective energy dissipation computed.
!
!  alphlmix :  multiple of pressure scale height used to compute lmix.
!-----------------------------------------------------------------------

INTEGER                                         :: ilcnvct
INTEGER                                         :: iemix
INTEGER                                         :: iyemix
INTEGER                                         :: ipsimix
INTEGER                                         :: ipcnvct
INTEGER                                         :: iscnvct

!-----------------------------------------------------------------------
!  Neutrino advection parameter
!-----------------------------------------------------------------------
!  adrftcv: neutrino advection parameter -
!
!     if vdriftk < adrftcv*ulcnvct
!
!  where vdriftk is the drift velocity of neutrinos of energy zone k, ulcnvct is the convective
!   velocity, and adrftcv is an arbitrary parameter of order unity, neutrinos are advected with
!   the matter.
!
!     if vdriftk > adrftcv*ulcnvct
!
!  neutrinos are not advected with the matter. In this case psi0aa(k) is set equal to psi0(jf,k,n)
!   so the the neutrinos of energy k delivered to the zone are equal to those removed, resulting
!   in no net advection.
!-----------------------------------------------------------------------

REAL(KIND=double)                               :: adrftcv
REAL(KIND=double)                               :: alphlmix

!-----------------------------------------------------------------------
!  Explosion simulation control
!-----------------------------------------------------------------------
!  i_bomb       : 0 - no bomb
!                 1 - bomb
!
!  e_bomb       : energy added to generate explosion (ergs).
!
!  bomb_time    : time during which energy is added.
!
!  t_start_bomb : time at which energy addition begins.
!
!  jexpl_min    : inner zone for energy addition.
!
!  jexpl_max    : outer zone for energy addition.
!-----------------------------------------------------------------------

INTEGER                                         :: i_bomb
INTEGER                                         :: jexpl_min
INTEGER                                         :: jexpl_max

REAL(KIND=double)                               :: e_bomb
REAL(KIND=double)                               :: bomb_time
REAL(KIND=double)                               :: t_start_bomb

  101 FORMAT (a128)
  103 FORMAT (a6)
  105 FORMAT (1x,a128)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)
  131 FORMAT (20x,i10,5x,e15.8)
  133 FORMAT (1x,a6,14x,i10,5x,es15.8,22x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)
 1001 FORMAT (' Allocation problem for array ',a10,' in read_pack_hydro_keys')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_pack_hydro_keys')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (uset(nx), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uset      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

ihydro               = 1
irelhy               = 0
ilapsehy             = 0
iubcjmn              = 0
iubcjmx              = 0
iuset                = 0
ipbnd                = 0
idtj                 = 0
jtv                  = 0
jtvshk               = 0
ncyshk               = 0
ncyrho               = 0
ipq                  = 1
ilcnvct              = 0
iemix                = 0
iyemix               = 0
ipsimix              = 0
ipcnvct              = 0
iscnvct              = 0
i_bomb               = 0
jexpl_min            = 0
jexpl_max            = 0

ubcjmn               = zero
ubcjmx               = zero
r_u                  = zero
pbound               = zero
rhojtv               = zero
tcntrl(1)            = 0.5d0
tcntrl(2)            = 1.d-02
tcntrl(3)            = 3.d-02
tcntrl(4)            = 0.5d0
tcntrl(5)            = 1.d-02
tcntrl(6)            = 3.d-02
tcntrl(7)            = 0.5d0
tcntrl(8)            = 1.d-02
tcntrl(9)            = 3.d-02
tcntrl(10)           = 1.2d0
tcntrl(31)           = 1.d-01
dtst1                = 1.d+06
dtst2                = 1.d-10
dtst3                = 1.d+17
ttst1                = 1.0d-3
q0_x                 = 1.d0
adrftcv              = zero
alphlmix             = zero
e_bomb               = zero
bomb_time            = zero
t_start_bomb         = zero
degen                = 2.d+01

uset                 = zero

i_hydro_data         = 0
d_hydro_data         = zero

!-----------------------------------------------------------------------
!
!                   \\\\\ READ HYDRO KEYS /////
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
!              \\\\\ PROBLEM CONTROL PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ihydro
!-----------------------------------------------------------------------

  IF ( type == 'ihydro' ) THEN
    READ (line ,111) ihydro
    IF ( iskip == 0 ) WRITE (nprint,113) type,ihydro,name
    CYCLE
  END IF ! type = 'ihydro'

!-----------------------------------------------------------------------
!  irel
!-----------------------------------------------------------------------

  IF ( type == 'irel  ' ) THEN
    READ (line ,111) int

    IF ( name == 'irelhy  ' ) THEN
      irelhy         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,irelhy,name
      CYCLE
    END IF ! name = 'irelhy  '

    IF ( name == 'ilapsehy' ) THEN
      ilapsehy       = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ilapsehy,name
      CYCLE
    END IF ! name = 'ilapsehy'

  END IF ! type = 'irel  '

!-----------------------------------------------------------------------
!  degen
!-----------------------------------------------------------------------

  IF ( type == 'degen ' ) THEN
    READ (line ,141) degen
    IF ( iskip == 0 ) WRITE (nprint,143) type,degen,name
    CYCLE
  END IF ! type = 'degen'

!-----------------------------------------------------------------------
!
!                 \\\\\ BOUNDARY CONDITIONS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ubcjm
!-----------------------------------------------------------------------

  IF ( type == 'ubcjm ' ) THEN
    READ (line ,131) int,rl

    IF ( name == 'ubcjmn  ' ) THEN
      iubcjmn        = int
      ubcjmn         = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,iubcjmn,ubcjmn,name
      CYCLE
    END IF ! name = 'ubcjmn  '

    IF ( name == 'ubcjmx  ' ) THEN
      iubcjmx        = int
      ubcjmx         = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,iubcjmx,ubcjmx,name
      CYCLE
    END IF ! name = 'ubcjmx  '

  END IF ! type = 'ubcjm '

!-----------------------------------------------------------------------
!  iuset
!-----------------------------------------------------------------------

  IF ( type == 'iuset ' ) THEN

    IF ( name == 'iuset   ') THEN
      READ (line ,111) iuset
      IF ( iskip == 0 ) WRITE (nprint,113) type,iuset,name
      CYCLE
    END IF ! name = 'iuset   '

    IF ( name == 'uset    ' ) THEN
      READ (line,131) j,rl
      uset(j)        = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,j,uset(j),name
      CYCLE
    END IF ! name = 'uset    '

  END IF ! type = 'iuset '

!-----------------------------------------------------------------------
!  r_u
!-----------------------------------------------------------------------

  IF ( type == 'r_u   ' ) THEN
    READ (line ,141) r_u
    IF ( iskip == 0 ) WRITE (nprint,143) type,r_u,name
    CYCLE
  END IF ! type = 'r_u   '

!-----------------------------------------------------------------------
!  ipbnd
!-----------------------------------------------------------------------

  IF ( type == 'ipbnd ' ) THEN
    READ (line ,131) int,rl
    ipbnd            = int
    pbound           = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,ipbnd,pbound,name
    CYCLE
  END IF ! type = 'ipbnd '

!-----------------------------------------------------------------------
!
!              \\\\\ TIME AND TIME STEP CONTROLS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  vtime
!-----------------------------------------------------------------------

  IF ( type == 'vtime ' ) THEN

    IF ( name == 'idtj    ') THEN
      READ (line ,111) int
      idtj           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,idtj,name
      CYCLE
    END IF ! name = 'idtj    '

    IF ( name == 'jtv     ') THEN
      READ (line ,111) int
      jtv            = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jtv,name
      CYCLE
    END IF ! name  'jtv     '

    IF ( name == 'jtvshk  ') THEN
      READ (line ,111) int
      jtvshk         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jtvshk,name
      CYCLE
    END IF ! name = 'jtvshk  '

    IF ( name == 'ncyshk  ') THEN
      READ (line ,111) int
      ncyshk         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ncyshk,name
      CYCLE
    END IF ! name = 'ncyshk  '

    IF ( name == 'ncyrho  ') THEN
      READ (line ,111) int
      ncyrho         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ncyrho,name
      CYCLE
    END IF ! name = 'ncyrho  '

    IF ( name == 'rhojtv  ') THEN
      READ (line ,141) rl
      rhojtv         = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,rhojtv,name
      CYCLE
    END IF ! name = 'rhojtv  '

  END IF ! type = 'vtime '

!-----------------------------------------------------------------------
!  tcntrl
!-----------------------------------------------------------------------

  IF ( type == 'tcntrl' ) THEN

    IF ( name == 'tcntrl_hydro' ) THEN
      READ (line ,131) n,rl
      tcntrl(n)     = rl
      IF ( iskip == 0 ) WRITE (nprint,133) type,n,tcntrl(n),name
      CYCLE
    END IF ! name = 'tcntrl_hydro'

  END IF ! type = 'tcntrl'

!-----------------------------------------------------------------------
!  dcntrl
!-----------------------------------------------------------------------

  IF ( type == 'dcntrl' ) THEN
    READ (line ,141) rl

    IF ( name == 'dtst1   ' ) THEN
      dtst1          = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtst1,name
      CYCLE
    END IF ! name = 'dtst1   '

    IF ( name == 'dtst2   ' ) THEN
     dtst2          = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtst2,name
      CYCLE
    END IF ! name = 'dtst2   '

    IF ( name == 'dtst3   ' ) THEN
      dtst3          = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtst3,name
      CYCLE
    END IF ! name = 'dtst3   '

    IF ( name == 'ttst1   ' ) THEN
      ttst1          = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,ttst1,name
      CYCLE
    END IF ! name = 'ttst1   '

  END IF ! type = 'dcntrl'

!-----------------------------------------------------------------------
!
!                     \\\\\ SHOCK ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  shock
!-----------------------------------------------------------------------

  IF ( type == 'shock ' ) THEN
    READ (line ,131) int,rl
    ipq              = int
    q0_x             = rl
    IF ( iskip == 0 ) WRITE (nprint,133) type,ipq,q0_x,name
    CYCLE
  END IF ! type = 'shock '

!-----------------------------------------------------------------------
!
!                \\\\\ CONVECTION PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  conv
!-----------------------------------------------------------------------

  IF ( type == 'conv  ' ) THEN
    IF ( name == 'ilcnvct') THEN
      READ (line ,111) ilcnvct
      IF ( iskip == 0 ) WRITE (nprint,113) type,ilcnvct,name
      CYCLE
    END IF ! name = 'ilcnvct'
    IF ( name == 'iemix') THEN
      READ (line ,111) iemix
      IF ( iskip == 0 ) WRITE (nprint,113) type,iemix,name
      CYCLE
    END IF ! name = 'iemix'
    IF ( name == 'iyemix') THEN
      READ (line ,111) iyemix
      IF ( iskip == 0 ) WRITE (nprint,113) type,iyemix,name
      CYCLE
    END IF ! name = 'iyemix'
    IF ( name == 'ipsimix') THEN
      READ (line ,111) ipsimix
      IF ( iskip == 0 ) WRITE (nprint,113) type,ipsimix,name
      CYCLE
    END IF ! name = 'ipsimix'
    IF ( name == 'ipcnvct') THEN
      READ (line ,111) ipcnvct
      IF ( iskip == 0 ) WRITE (nprint,113) type,ipcnvct,name
      CYCLE
    END IF ! name = 'ipcnvct'
    IF ( name == 'iscnvct') THEN
      READ (line ,111) iscnvct
      IF ( iskip == 0 ) WRITE (nprint,113) type,iscnvct,name
      CYCLE
    END IF ! name = 'iscnvct'
    IF ( name == 'adrftcv') THEN
      READ (line ,141) adrftcv
      IF ( iskip == 0 ) WRITE (nprint,143) type,adrftcv,name
      CYCLE
    END IF ! name = 'adrftcv'
    IF ( name == 'alphlmix') THEN
      READ (line ,141) alphlmix
      IF ( iskip == 0 ) WRITE (nprint,143) type,alphlmix,name
      CYCLE
    END IF ! name = 'alphlmix'
  END IF ! type = 'conv  '

!-----------------------------------------------------------------------
!
!                \\\\\ BOMB CONTROL PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  bomb
!-----------------------------------------------------------------------

  IF ( type == 'bomb  ' ) THEN

    IF ( name == 'i_bomb') THEN
      READ (line ,111) i_bomb
      IF ( iskip == 0 ) WRITE (nprint,113) type,i_bomb,name
      CYCLE
    END IF ! name = 'i_bomb'

    IF ( name == 'e_bomb') THEN
      READ (line ,141) e_bomb
      IF ( iskip == 0 ) WRITE (nprint,143) type,e_bomb,name
      CYCLE
    END IF ! name = 'e_bomb'

    IF ( name == 'bomb_time') THEN
      READ (line ,141) bomb_time
      IF ( iskip == 0 ) WRITE (nprint,143) type,bomb_time,name
      CYCLE
    END IF ! name = 'bomb_time'

    IF ( name == 't_start_bomb') THEN
      READ (line ,141) t_start_bomb
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_start_bomb,name
      CYCLE
    END IF ! name = 't_start_bomb'

    IF ( name == 'jexpl_min') THEN
      READ (line ,111) jexpl_min
      IF ( iskip == 0 ) WRITE (nprint,113) type,jexpl_min,name
      CYCLE
    END IF ! name = 'jexpl_min'

    IF ( name == 'jexpl_max') THEN
      READ (line ,111) jexpl_max
      IF ( iskip == 0 ) WRITE (nprint,113) type,jexpl_max,name
      CYCLE
    END IF ! name = 'jexpl_max'
    
  END IF ! type == 'bomb  '

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
!                   \\\\\ PACK HYDRO KEYS /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE

i_hydro_data(1)                   = ihydro
i_hydro_data(2)                   = irelhy
i_hydro_data(3)                   = ilapsehy
i_hydro_data(4)                   = iubcjmn
i_hydro_data(5)                   = iubcjmx
i_hydro_data(6)                   = iuset
i_hydro_data(7)                   = ipbnd
i_hydro_data(8)                   = idtj
i_hydro_data(9)                   = jtv
i_hydro_data(10)                  = jtvshk
i_hydro_data(11)                  = ncyshk
i_hydro_data(12)                  = ncyrho
i_hydro_data(13)                  = ipq
i_hydro_data(14)                  = ilcnvct
i_hydro_data(15)                  = iemix
i_hydro_data(16)                  = iyemix
i_hydro_data(17)                  = ipsimix
i_hydro_data(18)                  = ipcnvct
i_hydro_data(19)                  = iscnvct
i_hydro_data(20)                  = i_bomb
i_hydro_data(21)                  = jexpl_min
i_hydro_data(22)                  = jexpl_max

d_hydro_data(1)                   = ubcjmn
d_hydro_data(2)                   = ubcjmx
d_hydro_data(3)                   = r_u
d_hydro_data(4)                   = pbound
d_hydro_data(6)                   = rhojtv
d_hydro_data(7)                   = tcntrl(1)
d_hydro_data(8)                   = tcntrl(2)
d_hydro_data(9)                   = tcntrl(3)
d_hydro_data(10)                  = tcntrl(4)
d_hydro_data(11)                  = tcntrl(5)
d_hydro_data(12)                  = tcntrl(6)
d_hydro_data(13)                  = tcntrl(7)
d_hydro_data(14)                  = tcntrl(8)
d_hydro_data(15)                  = tcntrl(9)
d_hydro_data(16)                  = tcntrl(10)
d_hydro_data(17)                  = dtst1
d_hydro_data(18)                  = dtst2
d_hydro_data(19)                  = dtst3
d_hydro_data(20)                  = ttst1
d_hydro_data(21)                  = q0_x
d_hydro_data(22)                  = adrftcv
d_hydro_data(23)                  = alphlmix
d_hydro_data(24)                  = e_bomb
d_hydro_data(25)                  = bomb_time
d_hydro_data(26)                  = t_start_bomb
d_hydro_data(27)                  = degen
d_hydro_data(28)                  = tcntrl(31)

DO i = 1,nx
  d_hydro_data(30+i)              = uset(i)
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (uset, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'uset      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE read_pack_hydro_keys
