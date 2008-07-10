SUBROUTINE read_pack_radhyd_keys( nreadp, nprint, iskip, c_radhyd_data, &
& i_radhyd_data, d_radhyd_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_radhyd_keys
!    Module:       read_pack_radhyd_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/02/04
!
!    Purpose:
!      To read in the radhyd keys defining the dimensions, geometry,
!       and principal parameters of the probkem, and pack them in a
!       character array and an integer array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nreadp        : unit number to read from
!  nprint        : unit number to print to
!  iskip         : echo data read flag
!  nrst          : cycle number at start or restart
!
!    Output arguments:
!  c_radhyd_data : character array of radhyd keys
!  i_radhyd_data : integer array of radhyd keys
!  d_radhyd_data : 64 bit real array of radhyd keys
!
!    Include files:
!      kind_module, numerical_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nreadp        ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

CHARACTER (len=2), INTENT(out), DIMENSION(20) :: c_radhyd_data ! character array of radhyd keys

INTEGER, INTENT(out), DIMENSION(50)           :: i_radhyd_data ! integer array of radhyd keys

REAL(KIND=double), INTENT(out), DIMENSION(50) :: d_radhyd_data ! 64 bit real array of radhyd keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line

INTEGER                          :: int           ! integer data variable to read in an interger datum

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!
!                    ----- CYCLE PARAMETERS -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ncycle : cycle number
!  ncymax : maximum cycle number
!-----------------------------------------------------------------------

INTEGER                                              :: ncycle
INTEGER                                              :: ncymax

!-----------------------------------------------------------------------
!
!               ----- TIME AND TIME STEP CONTROLS -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Times
!-----------------------------------------------------------------------
!  time     : the elapsed time since the initiation of the calculation
!   (i.e., since ncycle = 0)
!
!  t_start  : the time at the initiation of the calculation (typically
!   t_start = 0.0.
!
!  t_bounce : the time of core bounce.
!
!  t_stop   : the elapsed time after which calculation is terminated.
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: time
REAL(KIND=double)                                    :: t_start
REAL(KIND=double)                                    :: t_bounce
REAL(KIND=double)                                    :: t_stop
REAL(KIND=double)                                    :: tb_stop

!-----------------------------------------------------------------------
!  Time steps
!-----------------------------------------------------------------------
!  dtnph : the coordinate 'hydro' time step, set by hydro and nuclear
!   reactions (if the material is not in nse), between time cycle m and
!   time cycle m + 1.
!
!  dtnmh : the coordinate 'hydro' time step between time cycle m - 1
!   and time cycle m.
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: dtnph
REAL(KIND=double)                                    :: dtnmh

!-----------------------------------------------------------------------
!
!            ----- DIMENSIONS AND GEOMETRY OF PROBLEM -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  ndim    : number of geometric dimensions
!  ngeomx  : x-geometry flag
!  ngeomy  : y-geometry flag
!  ngeomz  : z-geometry flag
!  nleftx  : lower x-boundary condition flag
!  nlefty  : lower y-boundary condition flag
!  nleftz  : lower z-boundary condition flag
!  nrightx : upper x-boundary condition flag
!  nrighty : upper y-boundary condition flag
!  nrightz : upper z-boundary condition flag
!-----------------------------------------------------------------------

INTEGER                                              :: ndim
INTEGER                                              :: ngeomx
INTEGER                                              :: ngeomy
INTEGER                                              :: ngeomz
INTEGER                                              :: nleftx
INTEGER                                              :: nlefty
INTEGER                                              :: nleftz
INTEGER                                              :: nrightx
INTEGER                                              :: nrighty
INTEGER                                              :: nrightz

!-----------------------------------------------------------------------
!
!                 ----- LAGRANGIAN OR EULERIAN -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  lagr = ye : Lagrangian hydrodynamics.
!  lagr = no : Eulerian hydrodynamics.
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: lagr

!-----------------------------------------------------------------------
!
!                       ----- MOVING GRID -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  m_grid = ye : Moving radial grid if lagr = no.
!  m_grid = no : Eulerian if lagr = no.
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: m_grid

!-----------------------------------------------------------------------
!
!                      ----- REGRID OPTION -----
!
!-----------------------------------------------------------------------
!  regrid            : regrid grid switch
!
!     m_grid = ye : regrid every int_pre_b or int_post_b cycles
!     m_grid = no : regrid option off
!
!  int_pre_b  : number of cycles between successive regrids before bounce
!  int_post_b : number of cycles between successive regrids after bounce
!  grid_frac  : fraction of a grid width grid is allowed to move per time step
!  rho_regrid : regrid up to density rho_regrid or 5 zones behind shock,
!                whichever is larger [g cm^{-3}]
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: regrid

INTEGER                                              :: int_pre_b
INTEGER                                              :: int_post_b

REAL(KIND=double)                                    :: grid_frac
REAL(KIND=double)                                    :: rho_regrid

!-----------------------------------------------------------------------
!
!                ----- INPOSED ROTATION OPTION ----- 
!
!-----------------------------------------------------------------------
!  rot  : imposed rotation switch
!
!      rot = ye : impart to the initial model a rotation with constant
!       angular velocity on cylinders according to the rotation law
!
!                                _      _ -1
!                               |      2 |
!                               |     r  |
!          Omega(r) = Omega_{0} | 1 + -- |
!                               |      2 |
!                               |_    A _|
!
!       where Omega(r) is the angular velocity, r the distance from the
!       rotation axis. A and beta (the ratio of the rotational energy to
!       the gravitational binding energy are free parameters that determine
!       the  rotational energy of the model and the distribution of angular
!       momentum. Omega_{0} is iterated until beta achieves the value
!       selected.
!
!     rot = no : No rotation imposed on the initial model.
!     m_grid = no : Eulerian if lagr = no.
!
!  A    : differential rotation parameter [km]
!  beta : ratio of the rotational energy to the gravitational binding
!   energy
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: rot

REAL(KIND=double)                                    :: A
REAL(KIND=double)                                    :: beta

!-----------------------------------------------------------------------
!
!                       ----- GRID WIGGLE -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  y_shft = ye : zone wiggle on
!  y_shft = no : zone wiggle off
!  dy_shift    : fraction of a zone width to wiggle grid
!  ncy_shift   : cycle number at which to commence zone wiggling
!  tb_dy_shift : time after bounce at which to stop zone wiggling
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: y_shft

INTEGER                                              :: ncy_shift

REAL(KIND=double)                                    :: dy_shift
REAL(KIND=double)                                    :: tb_dy_shift

!-----------------------------------------------------------------------
!
!                     ----- SHOCK SMOOTHING -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  v_diff : fraction of quantity difference to add to flux ingrid
!   aligned shocks for eliminating odd-even effect
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: v_diff

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

CHARACTER (len=2)                                    :: v_trans_0

!-----------------------------------------------------------------------
!
!               ----- YZ HYDRO SUBCYCLING OPTION -----
!
!-----------------------------------------------------------------------
!  sub_cy_yz     : yz-subcycling switch
!
!     sub_cy_xyz = ye : subcycle yz hydrodynamics relative to x hydrodynamics.
!     sub_cy_xyz = no : yz subvycle option off
!
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: sub_cy_yz

!-----------------------------------------------------------------------
!
!              ----- GLOBAL HYDRO TIME STEP OPTION -----
!
!-----------------------------------------------------------------------
!  t_step_xyz    : global hydro time-step switch
!
!     t_step_xyz = ye : hydro time step set to minimum of xyz hydro
!                        no yx hydro sybcycling
!     t_step_xyz = no : hydro time step set to minimum of x hydro; then
!                        yz hydro sybcycled (sub_cy_yz = ye)
!                        yz hydro t-step set to min( yz t-step, x t-step )
!                         (sub_cy_yz = no)
!
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: t_step_xyz

!-----------------------------------------------------------------------
!
!                  \\\\\ GRAVITATION OPTIONS /////
!
!-----------------------------------------------------------------------
!  i_grav        : gravitation key
!
!     i_grav = 1 : spherical symmetric gravity
!     i_grav = 2 : nonspherical components computed from by a Newtonian
!                   Poisson solver
!
!-----------------------------------------------------------------------

INTEGER                                              :: i_grav

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO EQUILIBRATION OPTION /////
!
!-----------------------------------------------------------------------
!  nu_equil        : neutrino equilibration switch
!
!     nu_equil = ye : neutrinos equilibrated with matter above density
!                      rho_equilibrate
!     nu_equil = no : neutrino equilibration step omitted
!
!  rho_equilibrate : density above whicn neutrinos equilibrated with
!                     mattew
!
!  t_equilibrate   : time after bounce at which neutrino equilibration
!                     turned on if off
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: nu_equil

REAL(KIND=double)                                    :: rho_equilibrate
REAL(KIND=double)                                    :: t_equilibrate

!-----------------------------------------------------------------------
!
!              \\\\\ GALILEAN TRANSFORMATION OPTION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  G_trns : Galilean transformation switch
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: G_trns

!-----------------------------------------------------------------------
!
!                       ----- REGRIDDING -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  rezn = ye : regrid on startup or restart.
!  rezn = no : no regriding on startup or restart.
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: rezn

!-----------------------------------------------------------------------
!  Lagrangian rezoning controls
!-----------------------------------------------------------------------
!  n_lgrgrid : 1 - generate a smooth Lagrangian grid
!              2 - generate three separate grids from m_1 separated m_2 and m_3
!
!  m_1       : mass of first zone
!
!  m_2       : mass separating grid 1 from grid 2
!
!  m_3       : mass separating grid 2 from grid 3
!
!  n1zoom    : number of zones between m_1 and m_2
!
!  n2zoom    : number of zones between m_2 and m_3
!
!  n3zoom    : number of zones between m_3 and the outer edge
!
!  zoome1    : mass ratio for rezoning the first n1 zones
!
!  zoome2    : mass ratio for rezoning the remaining zones
!
!  zoome3    : mass ratio for rezoning the n3zoom zones
!
!-----------------------------------------------------------------------
!         Eulerian rezoning controls
!-----------------------------------------------------------------------
!  n_eulgrid : 1 - generate a smooth Eulerian grid
!                   2 - generate three separate grids from r_1 separated r_2 and r_3
!
!  r_1       : radius of first zone
!
!  r_2       : radius separating grid 1 from grid 2
!
!  r_3       : radius separating grid 2 from grid 3
!
!  n1zoom    : number of zones between r_1 and r_2
!
!  n2zoom    : number of zones between r_2 and r_3
!
!  n3zoom    : number of zones between r_3 and the outer edge
!
!  zoome1    : radius ratio for rezoning the first n1 zones
!
!  zoome2    : radius ratio for rezoning the remaining zones
!
!  zoome3    : radius ratio for rezoning the n3zoom zones
!-----------------------------------------------------------------------

INTEGER                                              :: n_lgrgrid
INTEGER                                              :: n_eulgrid
INTEGER                                              :: n1zoom
INTEGER                                              :: n2zoom
INTEGER                                              :: n3zoom

REAL(KIND=double)                                    :: r_1
REAL(KIND=double)                                    :: r_2
REAL(KIND=double)                                    :: r_3
REAL(KIND=double)                                    :: m_1
REAL(KIND=double)                                    :: m_2
REAL(KIND=double)                                    :: m_3
REAL(KIND=double)                                    :: zoome1
REAL(KIND=double)                                    :: zoome2
REAL(KIND=double)                                    :: zoome3
REAL(KIND=double)                                    :: t_bounce_lagr_chg
REAL(KIND=double)                                    :: t_bounce_mgrd_chg

!-----------------------------------------------------------------------
!
!                       ----- GRID EXTENTS -----
!
!  imin : minimum unshifted x-array index
!
!  imax : maximum unshifted x-array index
!
!  jmin : minimum unshifted y-array index
!
!  jmax : maximum unshifted y-array index
!
!  kmin : minimum unshifted z-array index
!
!  kmax : maximum unshifted z-array index
!-----------------------------------------------------------------------

INTEGER                                              :: imin
INTEGER                                              :: imax
INTEGER                                              :: jmin
INTEGER                                              :: jmax
INTEGER                                              :: kmin
INTEGER                                              :: kmax

!-----------------------------------------------------------------------
!
!                    ----- NEUTRON STAR VELOCITY -----
!
!  vel_ns : velocity of the neutron star (cm s^{-2})
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: vel_ns

!-----------------------------------------------------------------------
!
!                    ----- PROBLEM DIMENSIONS -----
!
!  xmin : minimum value of x-coordinate
!
!  xmax : maximum value of x-coordinate
!
!  ymin : minimum value of y-coordinate
!
!  ymax : maximum value of y-coordinate
!
!  zmin : minimum value of z-coordinate
!
!  zmax : maximum value of z-coordinate
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: xmin
REAL(KIND=double)                                    :: xmax
REAL(KIND=double)                                    :: ymin
REAL(KIND=double)                                    :: ymax
REAL(KIND=double)                                    :: zmin
REAL(KIND=double)                                    :: zmax

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (a128)
  105 FORMAT (1x,a128)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)
  115 FORMAT (28x,2a)
  117 FORMAT (1x,a6,22x,a2,42x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE RADHYD DATA /////
!
!-----------------------------------------------------------------------

lagr                            = 'no'
m_grid                          = 'no'
regrid                          = 'no'
rot                             = 'no'
y_shft                          = 'no'
rezn                            = 'no'
v_trans_0                       = 'no'
sub_cy_yz                       = 'no'
t_step_xyz                      = 'no'
nu_equil                        = 'no'
G_trns                          = 'no'
i_grav                          = 1

ncycle                          = 0
ncymax                          = 9000000

ndim                            = -1
ngeomx                          = -1
ngeomy                          = -1
ngeomz                          = -1
nleftx                          = -1
nrightx                         = -1
nlefty                          = -1
nrighty                         = -1
nleftz                          = -1
nrightz                         = -1

imin                            = 1
imax                            = 1
jmin                            = 1
jmax                            = 1
kmin                            = 1
kmax                            = 1

vel_ns                          = zero

xmin                            = zero
xmax                            = zero
ymin                            = zero
ymax                            = 1.d0
zmin                            = zero
zmax                            = 2.d0

n_lgrgrid                       = 0
n_eulgrid                       = 0
n1zoom                          = 0
n2zoom                          = 0
n3zoom                          = 0


time                            = zero
t_start                         = zero
t_bounce                        = zero
t_stop                          = zero
tb_stop                         = zero

dtnph                           = 1.d+20
dtnmh                           = 1.d+20

t_bounce_lagr_chg               = 1.d-3
t_bounce_mgrd_chg               = 1.d-3

A                               = 1.d+20
beta                            = zero

int_pre_b                       = 25
int_post_b                      = 25
grid_frac                       = 1.d-01
rho_regrid                      = 1.d+06

ncy_shift                       = 100
dy_shift                        = 1.d-01
tb_dy_shift                     = 5.d-3

v_diff                          = 0.075d+0

rho_equilibrate                 = 1.d+14
t_equilibrate                   = 1.d+20

r_1                             = zero
r_2                             = zero
r_3                             = zero
m_1                             = zero
m_2                             = zero
m_3                             = zero
zoome1                          = zero
zoome2                          = zero
zoome3                          = zero


!-----------------------------------------------------------------------
!
!                   \\\\\ READ RADHYD DATA /////
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
!                    ----- CYCLE PARAMETERS -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  cycle
!-----------------------------------------------------------------------

  IF ( type == 'cycle ' ) THEN
    READ (line ,111) int

    IF ( name == 'ncycle  ') THEN
      ncycle     = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ncycle,name
      CYCLE
    END IF ! name = 'ncycle  '

    IF ( name == 'ncymax  ') THEN
      ncymax     = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ncymax,name
      CYCLE
    END IF ! name = 'ncymax  '

  END IF ! type = 'cycle '

!-----------------------------------------------------------------------
!
!               ----- TIME AND TIME STEP CONTROLS -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  sptime
!-----------------------------------------------------------------------

  IF ( type == 'sptime' ) THEN
    READ (line ,141) rl

    IF ( name == 'time    ') THEN
      time           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,time,name
      CYCLE
    END IF ! name = 'time    '

    IF ( name == 't_start ') THEN
      t_start         = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_start,name
      CYCLE
    END IF ! name = 't_start  '

    IF ( name == 't_bounce') THEN
      t_bounce       = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_bounce,name
      CYCLE
    END IF ! name = 't_bounce'

    IF ( name == 't_stop') THEN
      t_stop         = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_stop,name
      CYCLE
    END IF ! name = 't_stop'

    IF ( name == 'tb_stop') THEN
      tb_stop        = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,tb_stop,name
      CYCLE
    END IF ! name = 'tb_stop'

    IF ( name == 'dtnph   ') THEN
      dtnph          = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtnph,name
      CYCLE
    END IF ! name = 'dtnph   '

    IF ( name == 'dtnmh   ') THEN
      dtnmh          = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtnmh,name
      CYCLE
    END IF ! name = 'dtnmh   '

  END IF ! type = 'sptime'

!-----------------------------------------------------------------------
!
!            ----- DIMENSIONS AND GEOMETRY OF PROBLEM -----
!
!-----------------------------------------------------------------------

  IF ( type == 'geom  ' ) THEN

    IF ( name == 'ndim  ') THEN
      READ (line ,111) ndim
      IF ( iskip == 0 ) WRITE (nprint,113) type,ndim,name
      CYCLE
    END IF ! name = 'ndim  '

    IF ( name == 'ngeomx') THEN
      READ (line ,111) ngeomx
      IF ( iskip == 0 ) WRITE (nprint,113) type,ngeomx,name
      CYCLE
    END IF ! name = 'ngeomx'

    IF ( name == 'ngeomy') THEN
      READ (line ,111) ngeomy
      IF ( iskip == 0 ) WRITE (nprint,113) type,ngeomy,name
      CYCLE
    END IF ! name = 'ngeomy'

    IF ( name == 'ngeomz') THEN
      READ (line ,111) ngeomz
      IF ( iskip == 0 ) WRITE (nprint,113) type,ngeomz,name
      CYCLE
    END IF ! name = 'ngeomz'

    IF ( name == 'nleftx') THEN
      READ (line ,111) nleftx
      IF ( iskip == 0 ) WRITE (nprint,113) type,nleftx,name
      CYCLE
    END IF ! name = 'nleftx'

    IF ( name == 'nrightx') THEN
      READ (line ,111) nrightx
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrightx,name
      CYCLE
    END IF ! name = 'nrightx'

    IF ( name == 'nlefty') THEN
      READ (line ,111) nlefty
      IF ( iskip == 0 ) WRITE (nprint,113) type,nlefty,name
      CYCLE
    END IF ! name = 'nlefty'

    IF ( name == 'nrighty') THEN
      READ (line ,111) nrighty
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrighty,name
      CYCLE
    END IF ! name = 'nrighty'

    IF ( name == 'nleftz') THEN
      READ (line ,111) nleftz
      IF ( iskip == 0 ) WRITE (nprint,113) type,nleftz,name
      CYCLE
    END IF ! name = 'nleftz'

    IF ( name == 'nrightz') THEN
      READ (line ,111) nrightz
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrightz,name
      CYCLE
    END IF ! name = 'nrightz'

  END IF ! type = 'geom  '

!-----------------------------------------------------------------------
!
!                 ----- LAGRANGIAN OR EULERIAN -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  lagr
!-----------------------------------------------------------------------

  IF ( type == 'lagr  ' ) THEN

    IF ( name == 'lagr  ') THEN
      READ (line ,115) lagr 
      IF ( iskip == 0 ) WRITE (nprint,117) type,lagr,name
      CYCLE
    END IF ! name = 'lagr  '

    IF ( name == 'tb_lagr') THEN
      READ (line ,141) t_bounce_lagr_chg
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_bounce_lagr_chg,name
      CYCLE
    END IF ! name = 'tb_lagr'

  END IF ! type = 'lagr  '

    IF ( name == 'ndim  ') THEN
      READ (line ,111) ndim
      IF ( iskip == 0 ) WRITE (nprint,113) type,ndim,name
      CYCLE
    END IF ! name = 'ndim  '

!-----------------------------------------------------------------------
!
!                       ----- MOVING GRID -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  m_grid
!-----------------------------------------------------------------------

  IF ( type == 'm_grid' ) THEN

    IF ( name == 'm_grid  ') THEN
      READ (line ,115) m_grid 
      IF ( iskip == 0 ) WRITE (nprint,117) type,m_grid,name
      CYCLE
    END IF ! name = 'm_grid'

    IF ( name == 'tb_mgrd') THEN
      READ (line ,141) t_bounce_mgrd_chg
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_bounce_mgrd_chg,name
      CYCLE
    END IF ! name = 'tb_mgrd'

  END IF ! type = 'm_grid'

!-----------------------------------------------------------------------
!
!                       ----- REGRID OPTION -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  regrid
!-----------------------------------------------------------------------

  IF ( type == 'regrid' ) THEN

    IF ( name == 'regrid  ') THEN
      READ (line ,115) regrid 
      IF ( iskip == 0 ) WRITE (nprint,117) type,regrid,name
      CYCLE
    END IF ! name = 'regrid'

    IF ( name == 'grid_frac') THEN
      READ (line ,141) grid_frac
      IF ( iskip == 0 ) WRITE (nprint,143) type,grid_frac,name
      CYCLE
    END IF ! name = 'grid_frac'

    IF ( name == 'int_post_b') THEN
      READ (line ,111) int_post_b
      IF ( iskip == 0 ) WRITE (nprint,113) type,int_post_b,name
      CYCLE
    END IF ! name = 'int_post_b'

    IF ( name == 'int_pre_b') THEN
      READ (line ,111) int_pre_b
      IF ( iskip == 0 ) WRITE (nprint,113) type,int_pre_b,name
      CYCLE
    END IF ! name = 'dt_regrid'

    IF ( name == 'rho_regrid') THEN
      READ (line ,141) rho_regrid
      IF ( iskip == 0 ) WRITE (nprint,143) type,rho_regrid,name
      CYCLE
    END IF ! name = 'rho_regrid'

  END IF ! type = 'm_grid'

!-----------------------------------------------------------------------
!
!                    ----- INPOSED ROTATION ----- 
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  rot
!-----------------------------------------------------------------------

  IF ( type == 'rot   ' ) THEN

    IF ( name == 'rot   ') THEN
      READ (line ,115) rot 
      IF ( iskip == 0 ) WRITE (nprint,117) type,rot,name
      CYCLE
    END IF ! name = 'rot   '

    IF ( name == 'A     ') THEN
      READ (line ,141) A
      IF ( iskip == 0 ) WRITE (nprint,143) type,A,name
      CYCLE
    END IF ! name = 'A     '

    IF ( name == 'beta  ') THEN
      READ (line ,141) beta
      IF ( iskip == 0 ) WRITE (nprint,143) type,beta,name
      CYCLE
    END IF ! name = 'beta  '

  END IF ! type = 'rot   '

!-----------------------------------------------------------------------
!
!                       ----- GRID WIGGLE -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  y_shft
!-----------------------------------------------------------------------

  IF ( type == 'y_shft' ) THEN

    IF ( name == 'y_shft  ') THEN
      READ (line ,115) y_shft 
      IF ( iskip == 0 ) WRITE (nprint,117) type,y_shft,name
      CYCLE
    END IF ! name = 'y_shft'

    IF ( name == 'dy_shift') THEN
      READ (line ,141) dy_shift
      IF ( iskip == 0 ) WRITE (nprint,143) type,dy_shift,name
      CYCLE
    END IF ! name = 'dy_shift'

    IF ( name == 'ncy_shift') THEN
      READ (line ,111) ncy_shift 
      IF ( iskip == 0 ) WRITE (nprint,113) type,ncy_shift,name
      CYCLE
    END IF ! name = 'ncy_shift'

    IF ( name == 'tb_dy_shift') THEN
      READ (line ,141) tb_dy_shift
      IF ( iskip == 0 ) WRITE (nprint,143) type,tb_dy_shift,name
      CYCLE
    END IF ! name = 'tb_dy_shift'

  END IF ! type = 'y_shft'

!-----------------------------------------------------------------------
!
!                     ----- SHOCK SMOOTHING -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  v_diff
!-----------------------------------------------------------------------

  IF ( type == 'v_diff' ) THEN
    READ (line ,141) v_diff
    IF ( iskip == 0 ) WRITE (nprint,143) type,v_diff,name
    CYCLE
  END IF ! type = 'v_diff'

!-----------------------------------------------------------------------
!
!             ----- TRANSVERSE VELOCITIES ABOVE SHOCK -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  v_tran
!-----------------------------------------------------------------------

  IF ( type == 'v_tran' ) THEN
    READ (line ,115) v_trans_0
    IF ( iskip == 0 ) WRITE (nprint,117) type,v_trans_0,name
    CYCLE
  END IF ! type == 'v_tran'

!-----------------------------------------------------------------------
!
!                ----- YZ HYDRO SUBCYCLING OPTION -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  sub_cy
!-----------------------------------------------------------------------

  IF ( type == 'sub_cy' ) THEN
    READ (line ,115) sub_cy_yz
    IF ( iskip == 0 ) WRITE (nprint,117) type,sub_cy_yz,name
    CYCLE
  END IF ! type == 'sub_cy'

!-----------------------------------------------------------------------
!
!               ----- GLOBAL HYDRO TIME STEP OPTION -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  g_tstp
!-----------------------------------------------------------------------

  IF ( type == 'g_tstp' ) THEN
    READ (line ,115) t_step_xyz
    IF ( iskip == 0 ) WRITE (nprint,117) type,t_step_xyz,name
    CYCLE
  END IF ! type == 'g_tstp'

!-----------------------------------------------------------------------
!
!                      ----- GRAVITATION KEY -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  i_grav
!-----------------------------------------------------------------------

  IF ( type == 'i_grav' ) THEN
    READ (line ,111) i_grav
    IF ( iskip == 0 ) WRITE (nprint,113) type,i_grav,name
    CYCLE
  END IF ! type == 'i_grav'

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO EQUILIBRATION OPTION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  nequil
!-----------------------------------------------------------------------

  IF ( type == 'nequil' ) THEN

    IF ( name == 'nu_equil') THEN
      READ (line ,115) nu_equil 
      IF ( iskip == 0 ) WRITE (nprint,117) type,nu_equil,name
      CYCLE
    END IF ! name = 'nequil'

    IF ( name == 'rho_equilibrate') THEN
      READ (line ,141) rho_equilibrate
      IF ( iskip == 0 ) WRITE (nprint,143) type,rho_equilibrate,name
      CYCLE
    END IF ! name = 'rho_equilibrate'

    IF ( name == 't_equilibrate') THEN
      READ (line ,141) t_equilibrate
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_equilibrate,name
      CYCLE
    END IF ! name = 't_equilibrate'

  END IF ! type = 'nequil'

!-----------------------------------------------------------------------
!
!              \\\\\ GALILEAN TRANSFORMATION OPTION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  G_trns
!-----------------------------------------------------------------------

  IF ( type == 'G_trns' ) THEN
    READ (line ,115) G_trns
    IF ( iskip == 0 ) WRITE (nprint,117) type,G_trns,name
    CYCLE
  END IF ! type == 'G_trns'

!-----------------------------------------------------------------------
!
!                       ----- REGRIDDING -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  rezn
!-----------------------------------------------------------------------

  IF ( type == 'rezn  ' ) THEN

    IF ( name == 'rezn  ') THEN
      READ (line ,115) rezn
      IF ( iskip == 0 ) WRITE (nprint,117) type,rezn,name
      CYCLE
    END IF ! name = 'rezn  '

    IF ( name == 'n_lgrgrid') THEN
      READ (line ,111) n_eulgrid
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_lgrgrid,name
      CYCLE
    END IF ! name = 'n_lgrgrid'

    IF ( name == 'n_eulgrid') THEN
      READ (line ,111) n_eulgrid
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_eulgrid,name
      CYCLE
    END IF ! name = 'n_eulgrid'

    IF ( name == 'n1zoom ') THEN
      READ (line ,111) n1zoom
      IF ( iskip == 0 ) WRITE (nprint,113) type,n1zoom,name
      CYCLE
    END IF ! name = 'n1zoom '

    IF ( name == 'n2zoom ') THEN
      READ (line ,111) n2zoom
      IF ( iskip == 0 ) WRITE (nprint,113) type,n2zoom,name
      CYCLE
    END IF ! name = 'n2zoom '

    IF ( name == 'n3zoom ') THEN
      READ (line ,111) n3zoom
      IF ( iskip == 0 ) WRITE (nprint,113) type,n3zoom,name
      CYCLE
    END IF ! name = 'n3zoom '

    IF ( name == 'r_1   ') THEN
      READ (line ,141) r_1
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_1,name
      CYCLE
    END IF ! name = 'r_1   '

    IF ( name == 'r_2   ') THEN
      READ (line ,141) r_2
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_2,name
      CYCLE
    END IF ! name = 'r_2   '

    IF ( name == 'r_3   ') THEN
      READ (line ,141) r_3
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_3,name
      CYCLE
    END IF ! name = 'r_3   '

    IF ( name == 'm_1   ') THEN
      READ (line ,141) m_1
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_1,name
      CYCLE
    END IF ! name = 'm_1   '

    IF ( name == 'm_2   ') THEN
      READ (line ,141) m_2
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_2,name
      CYCLE
    END IF ! name = 'm_2   '

    IF ( name == 'm_3   ') THEN
      READ (line ,141) m_3
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_3,name
      CYCLE
    END IF ! name = 'm_3   '

    IF ( name == 'zoome1') THEN
      READ (line ,141) zoome1
      IF ( iskip == 0 ) WRITE (nprint,143) type,zoome1,name
      CYCLE
    END IF ! name = 'zoome1'

    IF ( name == 'zoome2') THEN
      READ (line ,141) zoome2
      IF ( iskip == 0 ) WRITE (nprint,143) type,zoome2,name
      CYCLE
    END IF ! name = 'zoome2'

    IF ( name == 'zoome3') THEN
      READ (line ,141) zoome3
      IF ( iskip == 0 ) WRITE (nprint,143) type,zoome3,name
      CYCLE
    END IF ! name = 'zoome3'

  END IF ! type = 'rezn  '

!-----------------------------------------------------------------------
!
!                       ----- GRID EXTENTS -----
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Minimum and maximum indices of unshifted arrays
!-----------------------------------------------------------------------

  IF ( type == 'evh1zn' ) THEN
    READ (line ,111) int

    IF ( name == 'imin    ' ) THEN
      imin           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,imin,name
      CYCLE
    END IF ! name = 'imin    '

    IF ( name == 'imax    ' ) THEN
      imax           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,imax,name
      CYCLE
    END IF ! name = 'imax    '

    IF ( name == 'jmin    ' ) THEN
      jmin           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jmin,name
      CYCLE
    END IF ! name = 'jmin    '


    IF ( name == 'jmax    ' ) THEN
      jmax           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,jmax,name
      CYCLE
    END IF ! name = 'jmax    '

    IF ( name == 'kmin    ' ) THEN
      kmin           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,kmin,name
      CYCLE
    END IF ! name = 'kmin    '

    IF ( name == 'kmax    ' ) THEN
      kmax           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,kmax,name
      CYCLE
    END IF ! name = 'kmax    '

  END IF ! type = 'evh1zn'

!-----------------------------------------------------------------------
!  Neutron star velocity
!-----------------------------------------------------------------------

  IF ( type == 'vel_ns' ) THEN
    READ (line ,141) vel_ns
    IF ( iskip == 0 ) WRITE (nprint,143) type,vel_ns,name
    CYCLE
  END IF ! type == 'vel_ns'

!-----------------------------------------------------------------------
!  Minimum and maximum coordinate dimensions
!-----------------------------------------------------------------------

  IF ( type == 'pb_dim' ) THEN
    READ (line ,141) rl

    IF ( name == 'xmin    ' ) THEN
      xmin           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,xmin,name
      CYCLE
    END IF ! name = 'xmin    '

    IF ( name == 'xmax    ' ) THEN
      xmax           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,xmax,name
      CYCLE
    END IF ! name = 'xmax    '

    IF ( name == 'ymin    ' ) THEN
      ymin           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,ymin,name
      CYCLE
    END IF ! name = 'ymin    '

    IF ( name == 'ymax    ' ) THEN
      ymax           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,ymax,name
      CYCLE
    END IF ! name = 'ymax    '

    IF ( name == 'zmin    ' ) THEN
      zmin           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,zmin,name
      CYCLE
    END IF ! name = 'zmin    '

    IF ( name == 'zmax    ' ) THEN
      zmax           = rl
      IF ( iskip == 0 ) WRITE (nprint,143) type,zmax,name
      CYCLE
    END IF ! name = 'zmax    '

  END IF ! type = 'pb_dim'

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
!                   \\\\\ PACK RADHYD DATA /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE
 
c_radhyd_data(1)                = lagr
c_radhyd_data(2)                = rezn
c_radhyd_data(3)                = m_grid
c_radhyd_data(4)                = regrid
c_radhyd_data(5)                = y_shft
c_radhyd_data(6)                = v_trans_0
c_radhyd_data(7)                = sub_cy_yz
c_radhyd_data(8)                = t_step_xyz
c_radhyd_data(9)                = nu_equil
c_radhyd_data(10)               = G_trns
c_radhyd_data(11)               = rot

i_radhyd_data(1)                = ncycle
i_radhyd_data(2)                = ncymax
i_radhyd_data(3)                = ndim
i_radhyd_data(4)                = ngeomx
i_radhyd_data(5)                = ngeomy
i_radhyd_data(6)                = ngeomz
i_radhyd_data(7)                = nleftx
i_radhyd_data(8)                = nrightx
i_radhyd_data(9)                = nlefty
i_radhyd_data(10)               = nrighty
i_radhyd_data(11)               = nleftz
i_radhyd_data(12)               = nrightz
i_radhyd_data(13)               = imin
i_radhyd_data(14)               = imax
i_radhyd_data(15)               = jmin
i_radhyd_data(16)               = jmax
i_radhyd_data(17)               = kmin
i_radhyd_data(18)               = kmax
i_radhyd_data(19)               = n_lgrgrid
i_radhyd_data(20)               = n_eulgrid
i_radhyd_data(21)               = n1zoom
i_radhyd_data(22)               = n2zoom
i_radhyd_data(23)               = n3zoom
i_radhyd_data(24)               = ncy_shift
i_radhyd_data(25)               = i_grav
i_radhyd_data(26)               = int_pre_b
i_radhyd_data(27)               = int_post_b

d_radhyd_data(1)                = time
d_radhyd_data(2)                = t_start
d_radhyd_data(3)                = t_bounce
d_radhyd_data(4)                = t_stop
d_radhyd_data(5)                = tb_stop
d_radhyd_data(6)                = dtnph
d_radhyd_data(7)                = dtnmh
d_radhyd_data(8)                = r_1
d_radhyd_data(9)                = r_2
d_radhyd_data(10)               = r_3
d_radhyd_data(11)               = m_1
d_radhyd_data(12)               = m_2
d_radhyd_data(13)               = m_3
d_radhyd_data(14)               = zoome1
d_radhyd_data(15)               = zoome2
d_radhyd_data(16)               = zoome3
d_radhyd_data(17)               = t_bounce_mgrd_chg
d_radhyd_data(18)               = t_bounce_lagr_chg
d_radhyd_data(19)               = grid_frac
d_radhyd_data(20)               = rho_regrid
d_radhyd_data(21)               = dy_shift
d_radhyd_data(22)               = tb_dy_shift
d_radhyd_data(23)               = v_diff
d_radhyd_data(24)               = xmin
d_radhyd_data(25)               = xmax
d_radhyd_data(26)               = ymin
d_radhyd_data(27)               = ymax
d_radhyd_data(28)               = zmin
d_radhyd_data(29)               = zmax
d_radhyd_data(30)               = rho_equilibrate
d_radhyd_data(31)               = t_equilibrate
d_radhyd_data(32)               = vel_ns
d_radhyd_data(33)               = A
d_radhyd_data(34)               = beta

RETURN
END SUBROUTINE read_pack_radhyd_keys
