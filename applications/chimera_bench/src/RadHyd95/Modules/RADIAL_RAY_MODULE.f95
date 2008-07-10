MODULE radial_ray_module
!-----------------------------------------------------------------------
!
!    File:         radial_ray_module
!    Module:       radial_ray_module
!    Type:         Module
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/01/04
!
!    Purpose:
!      Contains the radhyd variables on a given processor.
!
!    Include files:
!      kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
SAVE

!-----------------------------------------------------------------------
!
!                  \\\\\ COMPUTATIONAL CYCLES /////
!
!-----------------------------------------------------------------------
!  ncycle : the cycle number of the calculation.
!
!  ncymax : probelm termination crierion. Calculation terminated when
!   ncycle > ncymax.
!-----------------------------------------------------------------------

INTEGER                                              :: ncycle
INTEGER                                              :: ncymax

!-----------------------------------------------------------------------
!
!                       \\\\\ TIMES /////
!
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

REAL(KIND=double)                                    :: time
REAL(KIND=double)                                    :: t_start
REAL(KIND=double)                                    :: t_bounce
REAL(KIND=double)                                    :: t_stop
REAL(KIND=double)                                    :: tb_stop

!-----------------------------------------------------------------------
!
!           \\\\\ TIME STEPS AND TIME STEP CONSTRAINTS /////
!
!-----------------------------------------------------------------------
!  jdt            : radial zone for which dt(j,ij_ray,ik_ray) us a minimum
!
!  j_radial_dt    : x (radial ray) responsible for each dt_process
! 
!  j_angular_dt   : y (angular) ray responsible for each dt_process
!
!  j_azimuthal_dt : z (azimuhal) ray responsible for each dt_process
!
!  dtnph          : the coordinate 'hydro' time step, set by hydro and
!   nuclear reactions (if the material is not in nse), between time
!   cycle m and time cycle m + 1.
!
!  dtnmh          : the coordinate 'hydro' time step between time cycle
!   m - 1 and time cycle m.
!
!  dtnph_trans  : the time step for absorption, emission, and
!   scattering, set by the neutrino occupation distributions and
!   variables affected by the neutrino-matter interactions, between
!   time cycle m and time cycle m + 1.
!
!  dt_process     : minimum time step for process i
!
!  dt             : minimum time step for process i along radial ray
!   ij_ray,ik_ray
!
!  dtime_trans    : minimum time step source and transport from radial
!   ray ij_ray,ik_ray
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: jdt
INTEGER, ALLOCATABLE, DIMENSION(:)                   :: j_radial_dt
INTEGER, ALLOCATABLE, DIMENSION(:)                   :: j_angular_dt
INTEGER, ALLOCATABLE, DIMENSION(:)                   :: j_azimuthal_dt

REAL(KIND=double)                                    :: dtnph
REAL(KIND=double)                                    :: dtnmh
REAL(KIND=double)                                    :: dtnph_trans
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dt_process
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: dt
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)       :: dtime_trans

!-----------------------------------------------------------------------
!
!              \\\\\ GEOMETRY AND BOUNDARY FLAGS /////
!
!-----------------------------------------------------------------------
!  ndim    : number of geometric dimensions
!
!  ngeomx  : x-geometry flag
!
!  ngeomy  : y-geometry flag
!
!  ngeomz  : z-geometry flag
!
!  nleftx  : lower x-boundary condition flag
!
!  nlefty  : lower y-boundary condition flag
!
!  nleftz  : lower z-boundary condition flag
!
!  nrightx : upper x-boundary condition flag
!
!  nrighty : upper y-boundary condition flag
!
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
!         \\\\\ MINIMUM AND MAXUMUM SPATIAL ARRAY INDICES /////
!
!-----------------------------------------------------------------------
!  imin : minimum x-array index
!
!  imax : maximum x-array index
!
!  jmin : minimum y-array index
!
!  jmax : maximum y-array index
!
!  kmin : minimum z-array index
!
!  kmax : maximum z-array index
!-----------------------------------------------------------------------

INTEGER                                              :: imin
INTEGER                                              :: imax
INTEGER                                              :: jmin
INTEGER                                              :: jmax
INTEGER                                              :: kmin
INTEGER                                              :: kmax

!-----------------------------------------------------------------------
!
!                \\\\\ NEUTRON STAR PARAMETERS /////
!
!-----------------------------------------------------------------------
!  d_vel_ns : change in the neutron star velocity in a time step
!
!  vel_ns   : velocity of the neutron star
!
!  mass_ns  : mass of the neutron star
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: d_x_vel_ns
REAL(KIND=double)                                    :: d_y_vel_ns
REAL(KIND=double)                                    :: d_z_vel_ns
REAL(KIND=double)                                    :: d_vel_ns
REAL(KIND=double)                                    :: x_vel_ns
REAL(KIND=double)                                    :: y_vel_ns
REAL(KIND=double)                                    :: z_vel_ns
REAL(KIND=double)                                    :: vel_ns
REAL(KIND=double)                                    :: mass_ns

!-----------------------------------------------------------------------
!
!                    \\\\\ PROBLEM DIMENSIONS /////
!
!-----------------------------------------------------------------------
!  xmin   : minimum value of x-coordinate
!
!  xmax   : maximum value of x-coordinate
!
!  ymin   : minimum value of y-coordinate
!
!  ymax   : maximum value of y-coordinate
!
!  zmin   : minimum value of z-coordinate
!
!  zmax   : maximum value of z-coordinate
!
!  xmin_i : minimum value of x-coordinate before modified in rezone
!
!  xmax_i : maximum value of x-coordinate before modified in rezone
!
!  ymin_i : minimum value of y-coordinate before modified in rezone
!
!  ymax_i : maximum value of y-coordinate before modified in rezone
!
!  zmin_i : minimum value of z-coordinate before modified in rezone
!
!  zmax_i : maximum value of z-coordinate before modified in rezone
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: xmin
REAL(KIND=double)                                    :: xmax
REAL(KIND=double)                                    :: ymin
REAL(KIND=double)                                    :: ymax
REAL(KIND=double)                                    :: zmin
REAL(KIND=double)                                    :: zmax

REAL(KIND=double)                                    :: xmin_i
REAL(KIND=double)                                    :: xmax_i
REAL(KIND=double)                                    :: ymin_i
REAL(KIND=double)                                    :: ymax_i
REAL(KIND=double)                                    :: zmin_i
REAL(KIND=double)                                    :: zmax_i

!-----------------------------------------------------------------------
!
!                  \\\\\ LAGRANGIAN SWITCH /////
!
!-----------------------------------------------------------------------
!  lagr              : Lagrangian switch
!
!     lagr = ye : Lagrangian hydrodynamics.
!     lagr = no : Hydrodynamics via moving grid or Eulerian
!
!  t_bounce_lagr_chg : time after bounce for lagr = 'no'
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: lagr

REAL(KIND=double)                                    :: t_bounce_lagr_chg

!-----------------------------------------------------------------------
!
!                  \\\\\ MOVING GRID OPTION /////
!
!-----------------------------------------------------------------------
!  m_grid            : moving grid switch
!
!     m_grid = ye : Moving radial grid if lagr = no.
!     m_grid = no : Eulerian if lagr = no.
!
!  t_bounce_mgrd_chg : time after bounce for m_grid = 'no'
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: m_grid

REAL(KIND=double)                                    :: t_bounce_mgrd_chg

!-----------------------------------------------------------------------
!
!                      \\\\\ REGRID OPTION /////
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
!  rho_regrid : regrid up to density rho_regrid or 2 zones behind shock,
!                whichever is larger [g cm^{-3}]
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: regrid

INTEGER                                              :: int_pre_b
INTEGER                                              :: int_post_b

REAL(KIND=double)                                    :: grid_frac
REAL(KIND=double)                                    :: rho_regrid

!-----------------------------------------------------------------------
!
!                \\\\\ INPOSED ROTATION OPTION /////
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
!          energy
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: rot

REAL(KIND=double)                                    :: A
REAL(KIND=double)                                    :: beta

!-----------------------------------------------------------------------
!
!                 \\\\\ GRID WIGGLING OPTION /////
!
!-----------------------------------------------------------------------
!  y_shft      : grid wiggling toggle
!
!     y_shft = ye : grid wiggling turned on
!     y_shft = no : grid wiggling turned off
!
!  ncy_shift   : cycle number to commence grid wigling
!
!  dy_shift    : fraction of zone width to wiggle zone
!
!  tb_dy_shift : time after bounce tp stop zone wiggling
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: y_shft

INTEGER                                              :: ncy_shift

REAL(KIND=double)                                    :: dy_shift
REAL(KIND=double)                                    :: tb_dy_shift

!-----------------------------------------------------------------------
!
!                     \\\\\ SHOCK SMOOTHING /////
!
!-----------------------------------------------------------------------
!  v_diff : fraction of quantity difference to add to flux in grid-
!   aligned shocks for eliminating odd-even effect
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: v_diff

!-----------------------------------------------------------------------
!
!             \\\\\ ZERO TRANSVERSE VELOCITY OPTION /////
!
!-----------------------------------------------------------------------
!  v_trans_0     : zero transverse velocity above shock switch
!
!     v_trans_0 = ye : transverse velocities 0 above shock
!     v_trans_0 = no : transverse velocities 0 above shock computed
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: v_trans_0

!-----------------------------------------------------------------------
!
!               \\\\\ YZ HYDRO SUBCYCLING OPTION /////
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
!              \\\\\ GLOBAL HYDRO TIME STEP OPTION /////
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
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: t_step_xyz

!-----------------------------------------------------------------------
!
!                  \\\\\ GRAVITATION OPTIONS /////
!
!-----------------------------------------------------------------------
!  i_grav     : gravitation key
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
!              \\\\\ GALILEAN TRANSFORMATION SWITCH /////
!
!-----------------------------------------------------------------------
!  G_trns            : Gelilean transformation switch
!
!     G_trns = ye : Perform a Galilean transformation after each time
!                    step to keep the neutron star at rest in the grid
!                    center
!     G_trns = no : Omit the Galilean transformation
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: G_trns

!-----------------------------------------------------------------------
!
!                  \\\\\ REGRIDDING OPTION /////
!
!-----------------------------------------------------------------------
!  rezn      : regridding toggle
!
!     rezn = ye : regrid on startup or restart.
!     rezn = no : no regriding on startup or restart.
!-----------------------------------------------------------------------

CHARACTER (len=2)                                    :: rezn

!-----------------------------------------------------------------------
!
!             \\\\\ LAGRANGIAN REZONING CONTROLS /////
!
!-----------------------------------------------------------------------
!  n_lgrgrid : Lagrangian gridding flag
!
!     n_lgrgrid : 1 - generate a smooth Lagrangian grid
!     n_lgrgrid : 2 - generate three separate grids from m_1 to m_2,
!      m_2 to m_3, and m_3 to the outer edge
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
!-----------------------------------------------------------------------
!
!               \\\\\ EULERIAN REZONING CONTROLS /////
!
!-----------------------------------------------------------------------
!  n_eulgrid : Eulerian gridding flag
!         
!     n_eulgrid : 1 - generate a smooth Eulerian grid
!     n_eulgrid : 2 - generate three separate grids from r_1 to r_2,
!      r_2 to r_3, and r_3 to the outer edge
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
!
!
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

!-----------------------------------------------------------------------
!
!                   \\\\\ BOUNDARY VALUES /////
!
!-----------------------------------------------------------------------
!  u_bcl    : right boundary x-velocity
!
!  v_bcl    : right boundary y-velocity
!
!  w_bcl    : right boundary z-velocity
!
!  r_bcl    : right boundary density
!
!  p_bcl    : right boundary pressure
!
!  ei_bcl   : right boundary internal energy
!
!  ye_bcl   : right boundary electron fraction
!
!  temp_bcl : right boundary temperature
!
!  gc_bcl   : right boundary gamma
!
!  ge_bcl   : right boundary gamma
!
!  psi0_bcl : right boundary psi0
!
!  comp_bcl : right boundary composition
!
!  u_bcr    : right boundary x-velocity
!
!  v_bcr    : right boundary y-velocity
!
!  w_bcr    : right boundary z-velocity
!
!  r_bcr    : right boundary density
!
!  p_bcr    : right boundary pressure
!
!  ei_bcr   : right boundary internal energy
!
!  ye_bcr   : right boundary electron fraction
!
!  temp_bcr : right boundary temperature
!
!  gc_bcr   : right boundary gamma
!
!  ge_bcr   : right boundary gamma
!
!  psi0_bcr : right boundary psi0
!
!  comp_bcr : right boundary composition
!-----------------------------------------------------------------------

REAL(KIND=double)                                    :: u_bcl
REAL(KIND=double)                                    :: v_bcl
REAL(KIND=double)                                    :: w_bcl
REAL(KIND=double)                                    :: r_bcl
REAL(KIND=double)                                    :: p_bcl
REAL(KIND=double)                                    :: ei_bcl
REAL(KIND=double)                                    :: ye_bcl
REAL(KIND=double)                                    :: temp_bcl
REAL(KIND=double)                                    :: gc_bcl
REAL(KIND=double)                                    :: ge_bcl
REAL(KIND=double)                                    :: psi0_bcl
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: comp_bcl

REAL(KIND=double)                                    :: u_bcr
REAL(KIND=double)                                    :: v_bcr
REAL(KIND=double)                                    :: w_bcr
REAL(KIND=double)                                    :: r_bcr
REAL(KIND=double)                                    :: p_bcr
REAL(KIND=double)                                    :: ei_bcr
REAL(KIND=double)                                    :: ye_bcr
REAL(KIND=double)                                    :: temp_bcr
REAL(KIND=double)                                    :: gc_bcr
REAL(KIND=double)                                    :: ge_bcr
REAL(KIND=double)                                    :: psi0_bcr
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: comp_bcr

!-----------------------------------------------------------------------
!
!            \\\\\ SOLID ANGLES SUBTENDED BY RAYS /////
!
!-----------------------------------------------------------------------
!  d_omega : solid angles subtended by radial rays
!
!  omega   : total solid angle subtended by all the radial rays
!
!  cos_theta : cos( theta )
!
!  sin_theta : sin( theta )
!
!  cos_phi   : cos ( phi )
!
!  sin_phi   : sin( phi )
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:)       :: d_omega
REAL(KIND=double)                                    :: omega
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: cos_theta
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: sin_theta
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: cos_phi
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: sin_phi

!-----------------------------------------------------------------------
!
!         \\\\\ COORDINATES - VALUES AT CYCLE BEGINNING /////
!
!-----------------------------------------------------------------------
!  x_ei  : x grid zone face locations at beginning of cycle
!
!  y_ei  : y grid zone face locations at beginning of cycle
!
!  z_ei  : z grid zone face locations at beginning of cycle
!
!  dx_ci : x_ei(i+1) - x_ei(i)
!
!  dy_ci : y_ei(i+1) - y_ei(i)
!
!  dz_ci : z_ei(i+1) - z_ei(i)
!
!  x_ci  : x grid zone midpoint locations at beginning of cycle
!
!  y_ci  : y grid zone midpoint locations at beginning of cycle
!
!  z_ci  : z grid zone midpoint locations at beginning of cycle
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: x_ei
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: y_ei
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: z_ei
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dx_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dy_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dz_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: x_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: y_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: z_ci

!-----------------------------------------------------------------------
!
!       \\\\\ COORDINATES - VALUES AFTER LAGRANGIAN UPDATE /////
!
!-----------------------------------------------------------------------
!  x_el  : x grid zone face locations after Lagrangian update
!
!  y_el  : y grid zone face locations after Lagrangian update
!
!  z_el  : z grid zone face locations after Lagrangian update
!
!  dx_cl : x_el(i+1) - x_el(i)
!
!  dy_cl : y_el(i+1) - y_el(i)
!
!  dz_cl : z_el(i+1) - z_el(i)
!
!  x_cl  : x grid zone midpoint locations after Lagrangian update
!
!  y_cl  : y grid zone midpoint locations after Lagrangian update
!
!  z_cl  : z grid zone midpoint locations after Lagrangian update
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: x_el
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: y_el
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: z_el
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: dx_cl
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: dy_cl
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: dz_cl
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: x_cl
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: y_cl
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: z_cl

!-----------------------------------------------------------------------
!
!            \\\\\ COORDINATES - VALUES AT CYCLE END /////
!
!-----------------------------------------------------------------------
!  x_ef  : x grid zone face locations at cycle end [cm]
!
!  y_ef  : y grid zone face locations at cycle end
!
!  z_ef  : z grid zone face locations at cycle end
!
!  dx_cf : x_ef(i+1) - x_ef(i) [cm]
!
!  dy_cf : y_ef(i+1) - y_ef(i)
!
!  dz_cf : z_ef(i+1) - z_ef(i)
!
!  x_cf  : x grid zone midpoint locations at cycle end [cm]
!
!  y_cf  : y grid zone midpoint locations at cycle end
!
!  z_cf  : z grid zone midpoint locations at cycle end
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: x_ef
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: y_ef
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: z_ef
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dx_cf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dy_cf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: dz_cf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: x_cf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: y_cf
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: z_cf

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  rho_c    : density: zone average [g cm^{-3}]
!
!  t_c      : temperature: zone average [K]
!
!  ye_c     : electron fraction: zone average
!
!  ei_c     : internal energy: zone average (ergs g^{-1})
!
!  e_v_c    : total energy minus kinetic energy [ergs g^{-1}]
!
!  p_c      : pressure: zone average [ergs cm^{-3}]
!
!  gc_c     : 1st adiabatic index
!
!  ge_c     : 1 + p/(ei*rho)
!
!  u_c      : zone centered average velocity x direction [cm s^{-1}]
!
!  u_e      : x grid edge velocity for moving grid option [cm s^{-1}]
!
!  v_c      : zone centered average velocity y direction [cm s^{-1}]
!
!  w_c      : zone centered average velocity z direction [cm s^{-1}]
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: rho_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: t_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: ye_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: ei_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: e_v_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: p_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: gc_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: ge_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: u_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: u_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: v_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: w_c

!-----------------------------------------------------------------------
!
!        \\\\\ STATE VARIABLES - AFTER LAGRANGIAN UPDATE /////
!
!-----------------------------------------------------------------------
!  rho_l    : density: zone average [g cm^{-3}]
!
!  u_l      : zone centered average velocity x direction [cm s^{-1}]
!
!  v_l      : zone centered average velocity y direction [cm s^{-1}]
!
!  w_l      : zone centered average velocity z direction [cm s^{-1}]
!-----------------------------------------------------------------------


REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: u_l
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: v_l
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: w_l
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: rho_l

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - ANGULAR AVERAGES /////
!
!-----------------------------------------------------------------------
!  rhobar(j)           : angular averaged density [g cm^{-3}]
!
!  p_bar (j)           : angular averaged matter pressure [ergs cm^{-3}]
!
!  e_bar (j)           : angular averaged matter internal energy [ergs]
!   cm^{-3})
!
!  vx_bar(j)           : angular averaged radial velocity [cm s^{-1}]
!
!  v2_bar(j)           : angular averaged velocity square [cm^{2} s^{-2}]
!
!  e_nu_c_bar(j)       : angular averaged neutrino energy density
!   [ergs cm^{-3}].
!
!  f_nu_e_bar(j)       : angular averaged neutrino energy flux
!   [ergs cm^{-2} s^{-1}].
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: rhobar
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: p_bar
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: e_bar
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: vx_bar
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: v2_bar
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: e_nu_c_bar
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: f_nu_e_bar

!-----------------------------------------------------------------------
!
!            \\\\\ STATE VARIABLES - INITIAL VALUES /////
!
!-----------------------------------------------------------------------
!  rho_ci : density: zone average (g cm^{-3}) at beginning of cycle
!
!  t_ci   : temperature: zone average (K) at beginning of cycle
!
!  ye_ci  : electron fraction: zone average at beginning of cycle
!
!  ei_ci  : internal energy: zone average (ergs g^{-1}) at beginning of
!   cycle
!
!  u_ci   : zone centered average velocity x direction (cm s^{-1}) at
!   beginning of cycle
!
!  v_ci   : zone centered average velocity y direction (cm s^{-1}) at
!   beginning of cycle
!
!  w_ci   : zone centered average velocity z direction (cm s^{-1}) at
!   beginning of cycle
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: rho_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: t_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: ye_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: ei_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: u_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: v_ci
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: w_ci

!-----------------------------------------------------------------------
!
!                     \\\\\ GR VARIABLES /////
!
!-----------------------------------------------------------------------
!  agr_e   : zone-edged value of the lapse function after gravitational
!   update
!
!  agr_c   : zone-centered value of the lapse function after
!   gravitational update
!
!  agr_e_r : zone-edged value of the lapse function prior to
!   gravitational update
!
!  agr_e_r : zone-centered value of the lapse function  prior to
!   gravitational update
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: agr_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: agr_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: agr_e_r
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: agr_c_r

!-----------------------------------------------------------------------
!
!              \\\\\ SHOCK STABILIZING VARIABLES /////
!
!-----------------------------------------------------------------------
!  flat_x   : variables aligned radial rays indicating the presence of a
!   radial shock
!
!  flat_y_x : variables aligned along radial rays indicating the presence
!   of an angular shock
!
!  flat_z_x : variables aligned along radial rays indicating the presence
!   of an azimuthal shock
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: flat_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: flat_y_x
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: flat_z_x

!-----------------------------------------------------------------------
!
!         \\\\\ GRAVITATION ACCELERATIONS AND POTENTIALS /////
!
!-----------------------------------------------------------------------
!  grav_x_c   : zone-centered x-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_y_c   : zone-centered y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_z_c   : zone-centered z-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_pot_c : zone-centered gravitational potential (energy required
!   to move a unit mass from x to infinity (erg g^{-1})
!
!  grav_pot_c_i : zone-centered gravitational potential (energy required
!   to move a unit mass from x to infinity (erg g^{-1}) saved for
!   subtraction from grav_pot_c computed at the end of the remap to
!   calculate dPhi/dt
!
!  gtot_pot_c : zone-centered gravitational potential (energy required
!   to move a unit mass from x to infinity with no mass overlying it
!   (erg g^{-1})
!
!  grav_x_e   : zone-edged x-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_y_e   : zone-edged y-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_z_e   : zone-edged z-component of gravitational acceleration
!   (cm s^{-2} g^{-1})
!
!  grav_pot_e : zone-edged gravitational potential (energy required
!   to move a unit mass from x to infinity (erg g^{-1})
!
!  gtot_pot_e : zone-edged gravitational potential (energy required
!   to move a unit mass from x to infinity with no mass overlying it
!   (erg g^{-1})
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_x_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_y_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_z_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_pot_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_pot_c_i
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: gtot_pot_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_x_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_y_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_z_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: grav_pot_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: gtot_pot_e

!-----------------------------------------------------------------------
!
!          \\\\\ COMPOSITION VARIABLES - CURRENT VALUES /////
!
!-----------------------------------------------------------------------
!  nuc_number   : number of nuclear species (not counting representative
!   heavy nucleus).
!
!  xn_c         : mass fraction of the ith nucleus.
!
!  uburn_c      : cumulative energy generated in zone j by nuclear
!   reactions (ergs/gm).
!
!  be_nuc_rep_c : binding energy of the representative heavy nucleus
!   (MeV).
!
!  a_nuc_rep_c  : mass number of the representative heavy nucleus.
!
!  z_nuc_rep_c  : charge number of the representative heavy nucleus.
!
!  a_nuc_c      : mass number of the nth nuclear species.
!
!  z_nuc_c      : charge number of the nth nuclear species.
!
!  be_nuc_c     : binding energy of the nth nuclear species (MeV).
!
!  e_bind_zn_c  : total nuclear binding energy initially in each mass
!                  shell (ergs)
!
!  fluxbe_c     : total nuclear binding energy transferred during remap
!                  (ergs)
!
!  eb_c         : nuclear binding energy (ergs g^{-1})
!-----------------------------------------------------------------------

INTEGER                                              :: nuc_number

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: xn_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: uburn_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: be_nuc_rep_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: a_nuc_rep_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: z_nuc_rep_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: a_nuc_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: z_nuc_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)         :: be_nuc_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: e_bind_zn_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: fluxbe_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: eb_c

!-----------------------------------------------------------------------
!
!                      \\\\\ NSE FLAG /////
!
!-----------------------------------------------------------------------
!  nse_c : nuclear statistical equilibrium flag
!
!     nse(i,ij_ray,ik_ray) = 0 : material not in nuclear statistical
!      equilibrium; nuclear reaction network must be turned on to evolve
!      the matter composition.
!     nse(i,ij_ray,ik_ray) = 1 : material in nuclear statistical
!      equilibrium; nuclear reaction network turned off.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nse_c

!-----------------------------------------------------------------------
!
!                  \\\\\ EOS AND OPACITY GRID /////
!
!-----------------------------------------------------------------------
!  dgrid(i), tgrid(i), ygrid(i) : log(rho), log(t), ye space is overlain
!   with a uniform grid of 'dgrid(i)' divisions per unit change in
!   log(rho), 'tgrid(i)' divisions per unit change in log(t), and
!   'ygrid(i)' divisions per 0.5 change in ye. Equation of state, nuclear
!   reaction rate, and neutrino interaction variables at each radial zone
!   are interpolated from values at the nearest corners on this grid.
!
!  rhoes(k) : The variables dgrid, tgrid, and ygrid are each 3 element
!   arrays permitting different partitionings of log(rho), log(t), ye
!   space in different density regimes delimited by rhoes. These different
!   regimes are
!
!     regime 1 :             rho < rhoes(1)
!     regime 2 :        rhoes(1) < rho < rhoes(2)
!     regime 3 :             rhoes(2) < rho
!
!  idty(j,ij_ray,ik_ray) : the density regime (i.e., 1, 2, or 3) of radial zone j
!   as given by the above inequalities.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: idty

REAL(KIND=double), DIMENSION(3)                      :: dgrid
REAL(KIND=double), DIMENSION(3)                      :: tgrid
REAL(KIND=double), DIMENSION(3)                      :: ygrid
REAL(KIND=double), DIMENSION(2)                      :: rhoes

!-----------------------------------------------------------------------
!
!        \\\\\ INTERPOLATED EQUATION OF STATE VARIABLES /////
!
!-----------------------------------------------------------------------
!  aesv(j,i,ij_ray,ik_ray)  : equation of state dependent variable i at
!   radial zone j.
!
!     i = 1   : pressure
!     i = 2   : energy
!     i = 3   : entropy
!     i = 4   : neutron chemical potential
!     i = 5   : proton chemical potential
!     i = 6   : electron chemical potential
!     i = 7   : free neutron mass fraction
!     i = 8   : free proton mass fraction
!     i = 9   : heavy nucleus mass fraction
!     i = 10  : heavy nucleus mass number
!     i = 11  : heavy nucleus charge number
!     i = 12  : gamma1
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: aesv_c

!-----------------------------------------------------------------------
!
!                   \\\\\ NEUTRINO ARRAYS /////
!
!-----------------------------------------------------------------------
!  psi0_c(j,k,n,ij_ray,ik_ray) : the zero moment of the occupation
!   distribution for neutrinos at the midpoint of radial zone j, of energy
!   zone k, and of type n.
!
!  psi1_e(j,k,n,ij_ray,ik_ray) : the first moment of the occupation
!   distribution for neutrinos at the outer boundary radial zone j, of
!   energy zone k, and of type n.
!
!  nu_str_c(j,j_ray)   : x-component of the zone-centered neutrino stress
!   (dynes g^{-1})
!
!  nu_str_e(j,j_ray)   : x-component of the zone-edge neutrino stress
!   (dynes g^{-1})
!
!  rhs1_c(j,k,n,j_ray) : the right-hand side of the first moment of the
!   transport equation divided by psi1(j,k,n)/psi0(j,k,n)
!
!  dc_e(j,k,n,j_ray)   : the diffusion coefficient for neutrinos of
!   energy zone k, type n, at radial zone j and at time m + 1.
!
!  e_nu_c(j,j_ray)     : neutrino energy density (ergs cm^{-3}).
!
!  f_nu_e(j,j_ray)     : neutrino energy flux (ergs cm^{-2} s^{-1}).
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: psi0_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: psi1_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: nu_str_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: nu_str_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: rhs1_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: dc_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: e_nu_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:)     :: f_nu_e

!-----------------------------------------------------------------------
!
!                   \\\\\ NEUTRINO ENERGIES /////
!
!-----------------------------------------------------------------------
!  unu_c(j,k,ij_ray,ik_ray)   : the value of the zone-centered neutrino
!   energy of energy zone k at radial zone center j at time m (MeV).
!
!  unub_c(j,k,ij_ray,ik_ray)  : the value of the inner zone-edge neutrino
!   energy of energy zone k at radial zone center j at time m (MeV).
!
!  dunu_c(j,k,ij_ray,ik_ray)) : the width of energy zone k at radial zone
!   center j (MeV),  dunu(j,k) = unub(j,k+1) - unub(j,k) at time m
!
!  unue_e(j,k,ij_ray,ik_ray)  : the value of the zone-centered neutrino
!   energy of energy zone k at radial zone edge j at time m (MeV).
!
!  unube_e(j,k,ij_ray,ik_ray) : the value of the inner zone-edge neutrino
!   energy of energy zone k at radial zone edge j at time m (MeV).
!
!  dunue_e(j,k,ij_ray,ik_ray) : the width of energy zone k at radial zone
!   edge j (MeV), dunu(j,k) = unub(j,k+1) - unub(j,k)
!-----------------------------------------------------------------------

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: unu_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: unub_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: dunu_c
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: unue_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: unube_e
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:,:,:)   :: dunue_e


!-----------------------------------------------------------------------
!
!            \\\\\ EDIT PRINT UNIT NUMBER AND COUNTERS /////
!
!-----------------------------------------------------------------------
!  nprint : a parameter used in the call to an edit subroutine denoting
!   the unit number of the file to which the print file is to be sent.
!   It is usually redefined before a call but nominally carries the
!   value denoting the unit number of the file 'superdump' to which any
!   diagnostics are written during the course of a run.
!
!  nedc  : number of cycles since the last configuration edit
!
!  nede  : number of cycles since the last energy conservation edit
!
!  nedmi :
!
!  nedma : number of cycles since the last mass average edit
!
!  nedh  : number of cycles since the last hydrodynamic edit
!
!  nedps : number of cycles since the last pressure - stress edit
!
!  nedu  : number of cycles since the last energy edit
!
!  nedy  : number of cycles since the last composition edit
!
!  nedsc : number of cycles since the last entropy - chemical potential
!   edit
!
!  nedn  : number of cycles since the last differential neutrino edit
!
!  nedng : number of cycles since the last integral edit
!-----------------------------------------------------------------------

INTEGER                                              :: nprint

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedc
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nede
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedmi
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedma
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedh
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedps
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedu
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedy
INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedsc

INTEGER, ALLOCATABLE, DIMENSION(:,:,:)               :: nedn
INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:)             :: nedng

END MODULE radial_ray_module
