SUBROUTINE evh1_x_lagr_inout( imin, imax, nx, ny, nz, ij_ray, ik_ray,         &
& ij_ray_dim, ik_ray_dim, nez, nnu, x_ei, dx_ci, x_ci, x_el, dx_cl, x_cl,     &
& rhop, rho_l, tp, yep, eip, e_vp, up, u_e, vp, wp, u_l, v_l, w_l, nu_strcp,  &
& nu_strep, dtimep, dt, jdt, rhobarp, flat_x, grav_x_c, grav_pot_c, grav_x_e, &
& grav_pot_e, e_nu_c_bar, f_nu_e_bar, agr_c, agr_e, e_nu_c )
!-----------------------------------------------------------------------
!
!    File:         evh1_x_lagr_inout
!    Module:       evh1_x_lagr_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive the EVH1 arrays from radial_ray_module, execute the x-array
!       Lagrangian hydrodynamics, and return to radial_ray_module the updated
!       variables.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  imin         : lower unshifted x-array index
!  imax         : upper unshifted x-array index
!  nx           : x-array extent
!  ny           : y-array extent
!  nz           : z-array extent
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping with y
!  ik_ray_dim   : number of z-zones on a processor before swapping with z
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  x_ei         : initial x grid zone faces
!  dx_ci        : initial x_ei(i+1) - x_ei(i)
!  x_ci         : initial x grid zone midpoints
!  rhop         : initial densities (g cm^{-3})
!  tp           : initial temperatures (K)
!  yep          : initial electron fractions
!  eip          : initial internal energies (ergs g^{-1})
!  e_vp         : total energy minus kinetic energy (ergs g^{-1})
!  up           : initial velocity x-components [cm s^{-1}]
!  vp           : initial velocity y-components [cm s^{-1}]
!  wp           : initial velocity z-components [cm s^{-1}]
!  dtimep       : time step (s)
!  rhobarp      : mean density as a function of radius (g cm^{-3})
!  psi0         : zero moment of the neutrino occupation probability
!  nu_strcp     : zone-centered neutrino stress (dynes g^{-1})
!  nu_strep     : zone-edge neutrino stress (dynes g^{-1})
!  flat_x       : variable indicating the presence of radial shocks
!  grav_x_c     : zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_pot_c   : zone-centered gravitational potential (erg g^{-1})
!  grav_x_e     : zone-edged x-component of gravitational acceleration (cm s^{-2} g^{-1})
!  grav_pot_e   : zone-edged gravitational potential (erg g^{-1})
!  e_nu_c_bar   : angular averaged neutrino energy density (ergs cm^{-3})
!  f_nu_e_bar   : angular averaged neutrino energy flux (ergs cm^{-2} s^{-1})
!  agr_c        : zone-centered lapse function
!  agr_e        : zone-edged lapse function
!  e_nu_c       : neutrino energy density [ergs cm^{-3}]
!
!    Output arguments:
!  x_el         : updated x grid zone faces
!  dx_cl        : updated x_el(i+1,ij_ray,ik_ray) - x_el(i,ij_ray,ik_ray)
!  x_cl         : updated x grid zone midpoints
!  rhop         : updated densities (g cm^{-3})
!  rho_l        : updated densities (g cm^{-3})
!  tp           : updated temperatures (K)
!  yep          : updated electron fractions
!  eip          : updated internal energies (ergs g^{-1})
!  up           : updated centered velocity x-components [cm s^{-1}]
!  u_l          : velocity x-components after Lagrangian update [cm s^{-1}]
!  v_l          : velocity y-components after Lagrangian update [cm s^{-1}]
!  w_l          : velocity z-components after Lagrangian update [cm s^{-1}]
!  u_e          : updated edged velocity x-components [cm s^{-1}]
!  flat_x       : variable indicating the presence of radial shocks
!  dt           : minimum hydro time step restrictions (s)
!  jdt          : zones causing minimum time step restrictions dt
!
!    Subprograms called:
!  etotal                  : computes the total energies
!  sweepx                  : performs the x-Lagrangian hydro step
!  hydro_x_sweep_time_step : computes the x-sweep hydro time step restrictions
!      
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_x_module, evh1_sweep, evh1_sweep,
!  incrmnt_module, mdl_cnfg_module, mgfld_remap_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, half, third

USE e_advct_module, ONLY : rhomin_y_eadvect
USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY: aesv
USE evh1_global, ONLY: dt_g=>dt
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, u, v, w, ei, ye, temp, xa0, &
& dx0, xa, dx, nu_strc, nu_stre, rhobar, entrop, flat_s, g_force_c, &
& g_pot_c, g_force_e, g_pot_e, e_nu_cs=>e_nu_c, f_nu_e, lapse_c, lapse_e, &
& u_edge, egrav, e_v, p_nu, e_nu
USE incrmnt_module, ONLY: dtmpmn
USE mdl_cnfg_module, ONLY: rho_m=>rho,t_m=>t,ye_m=>ye
USE mgfld_remap_module, ONLY : r0i, t0i, r0l, dxl, xal, t0l, ei0l, ye0l, u0l


IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                       :: imin        ! lower unshifted x-array index
INTEGER, INTENT(in)                                       :: imax        ! upper unshifted x-array index
INTEGER, INTENT(in)                                       :: nx          ! x-array extent
INTEGER, INTENT(in)                                       :: ny          ! y-array extent
INTEGER, INTENT(in)                                       :: nz          ! z-array extent

INTEGER, INTENT(in)                                       :: ij_ray      ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                                       :: ik_ray      ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                                       :: ij_ray_dim  ! number of y-zones on a processor before swapping with y
INTEGER, INTENT(in)                                       :: ik_ray_dim  ! number of z-zones on a processor before swapping with z

INTEGER, INTENT(in)                                       :: nez         ! neutrino energy array extent
INTEGER, INTENT(in)                                       :: nnu         ! neutrino flavor array extent

REAL(KIND=double), INTENT(in)                             :: dtimep      ! time step

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)            :: x_ei        ! x-coordinate zone edge
REAL(KIND=double), INTENT(in), DIMENSION(nx)              :: dx_ci       ! x-coordinate zone thickness
REAL(KIND=double), INTENT(in), DIMENSION(nx)              :: x_ci        ! x-coordinate zone center

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: vp          ! y-velocity of zone [cm]

REAL(KIND=double), INTENT(in), DIMENSION(nx)                          :: rhobarp     ! mean density as a function of radius (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: nu_strcp    ! neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: nu_strep    ! neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: grav_x_c    ! zone-centered x-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: grav_pot_c  ! zone-centered gravitational potential (erg g^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)  :: grav_x_e    ! zone-edged x-component of gravitational acceleration (cm s^{-2} g^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)  :: grav_pot_e  ! zone-edged gravitational potential (erg g^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                          :: e_nu_c_bar  ! angular averaged neutrino energy density (ergs cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                        :: f_nu_e_bar  ! angular averaged neutrino energy flux (ergs cm^{-2} s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: agr_c       ! zone-centered lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)  :: agr_e       ! zone-edged lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: e_nu_c      ! neutrino energy density [ergs cm^{-3}]

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)             :: jdt         ! zone causing dt

REAL(KIND=double), INTENT(out), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: x_el        ! x-coordinate zone edge
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: dx_cl       ! x-coordinate zone thickness
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: x_cl        ! x-coordinate zone center
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: u_e         ! updated edged velocity x-components [cm s^{-1}]
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: u_l         ! velocity x-components after Lagrangian update [cm s^{-1}]
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: v_l         ! velocity y-components after Lagrangian update [cm s^{-1}]
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: w_l         ! velocity z-components after Lagrangian update [cm s^{-1}]
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: rho_l       ! density after Lagrangian update (g cm^{-3}) 
REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)   :: dt          ! minimum allowed time step
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: flat_x      ! variable indicating the presence of radial shocks

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: rhop        ! density (g cm^{-3})
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: tp          ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: yep         ! electron fraction
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: eip         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: e_vp        ! internal energy (ergs/g)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: up          ! radial velocity of zone [cm]
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: wp          ! z-velocity of zone [cm]

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                                        :: l_write = .false.

INTEGER                                                   :: ntot        ! parabolic array dimension
INTEGER                                                   :: jr_min      ! minimum radial zone index
INTEGER                                                   :: jr_max      ! maximum radial zone index
INTEGER                                                   :: i_nu        ! the maximum index for which rhobarp > rhomin_y_eadvect
INTEGER                                                   :: n_nu        ! the maximum padded index for which rhobarp > rhomin_y_eadvect

REAL(KIND=double)                                         :: dtime_hydro ! minimum hydro time_step for ij_ray,ik_ray

REAL(KIND=double), DIMENSION(nx)                          :: drl         ! zone width after Lagrangian step [cm]
REAL(KIND=double), DIMENSION(nx)                          :: rhol        ! density after Lagrangian step (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)                          :: rhoi        ! density before Lagrangian step (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)                          :: tl          ! temperature after Lagrangian step (K)
REAL(KIND=double), DIMENSION(nx)                          :: ti          ! temperature before Lagrangian step (K)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!      \\\\\ PACK VARIABLES IN 1-D ARRAYS WITH GHOST ZONES /////
!
!  |||||  zone-edged  : x_edge(i)   = x_edge(j)     = x_edge(n+6)  |||||
!  ||||| zone_center : x_center(i) = x_center(j+1) = x_center(n+6) |||||
!
!-----------------------------------------------------------------------

sweep                   = 'x'

!-----------------------------------------------------------------------
!  Add six ghost zones to 1D sweep arrays
!-----------------------------------------------------------------------

nmin                    = imin + 6
nmax                    = imax + 6
ntot                    = imax + 12
jr_min                  = imin + 1
jr_max                  = imax + 1

nu_strc                 = 0.d0
nu_stre                 = 0.d0
g_force_c               = zero
g_force_e               = zero

!-----------------------------------------------------------------------
!  Put state variables into 1D arrays, padding with 6 ghost zones.
!-----------------------------------------------------------------------

r        (nmin:nmax)    = rhop      (imin:imax,ij_ray,ik_ray)
temp     (nmin:nmax)    = tp        (imin:imax,ij_ray,ik_ray)
ye       (nmin:nmax)    = yep       (imin:imax,ij_ray,ik_ray)
u        (nmin:nmax)    = up        (imin:imax,ij_ray,ik_ray)
v        (nmin:nmax)    = vp        (imin:imax,ij_ray,ik_ray)
w        (nmin:nmax)    = wp        (imin:imax,ij_ray,ik_ray)
ei       (nmin:nmax)    = eip       (imin:imax,ij_ray,ik_ray)
rhobar   (nmin:nmax)    = rhobarp   (imin:imax)
entrop   (nmin:nmax)    = aesv      (jr_min:jr_max,3,ij_ray,ik_ray)
nu_strc  (nmin:nmax)    = nu_strcp  (imin:imax,ij_ray,ik_ray)
nu_stre  (nmin:nmax)    = nu_strep  (imin:imax,ij_ray,ik_ray)
g_force_c(nmin:nmax)    = grav_x_c  (imin:imax,ij_ray,ik_ray)
g_pot_c  (nmin:nmax)    = grav_pot_c(imin:imax,ij_ray,ik_ray)
g_force_e(nmin:nmax)    = grav_x_e  (imin:imax,ij_ray,ik_ray)
g_pot_e  (nmin:nmax)    = grav_pot_e(imin:imax,ij_ray,ik_ray)
egrav    (nmin:nmax)    = grav_pot_c(imin:imax,ij_ray,ik_ray)
e_nu_cs  (nmin:nmax)    = e_nu_c_bar(imin:imax)
f_nu_e   (nmin:nmax+1)  = f_nu_e_bar(imin:imax+1)
lapse_c  (nmin:nmax)    = agr_c     (imin:imax,ij_ray,ik_ray)
lapse_e  (nmin:nmax+1)  = agr_e     (imin:imax+1,ij_ray,ik_ray)

rho_m(jr_min:jr_max)    = rhop      (imin:imax,ij_ray,ik_ray)
t_m  (jr_min:jr_max)    = tp        (imin:imax,ij_ray,ik_ray)
ye_m (jr_min:jr_max)    = yep       (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Store neutrino pressure and energy if rhobarp >= rhomin_y_eadvect
!-----------------------------------------------------------------------

IF ( rhobarp(1) <= rhomin_y_eadvect ) THEN
  i_nu                   = 1
ELSE ! rhobarp(1) > rhomin_y_eadvect
  i_nu                   = MINLOC( rhobarp, DIM = 1, MASK = rhobarp > rhomin_y_eadvect )
END IF ! rhobarp(1) <= rhomin_y_eadvect
p_nu                     = zero
e_nu                     = zero
n_nu                     = i_nu + 6
iF ( i_nu > 1 ) THEN
  p_nu(nmin:n_nu)        = third * e_nu_c(imin:i_nu,ij_ray,ik_ray)
  e_nu(nmin:n_nu)        = e_nu_c(imin:i_nu,ij_ray,ik_ray)/r(nmin:n_nu)
END IF ! i_nu > 1

!-----------------------------------------------------------------------
!  Load Grid Coordinates
!-----------------------------------------------------------------------

xa0(nmin:nmax+1)        = x_ei (imin:imax+1)
dx0(nmin:nmax)          = dx_ci(imin:imax)
xa (nmin:nmax+1)        = x_ei (imin:imax+1)
dx (nmin:nmax)          = dx_ci(imin:imax)
 
!-----------------------------------------------------------------------
!  Store initial values
!-----------------------------------------------------------------------

r0i (nmin:nmax)         = rhop (imin:imax,ij_ray,ik_ray)
t0i (nmin:nmax)         = tp   (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

dt_g                    = dtimep

!-----------------------------------------------------------------------
!  Compute energies
!-----------------------------------------------------------------------

CALL etotal( .false. )

!-----------------------------------------------------------------------
!
!                   \\\\\ LAGRANGIAN HYDRO /////
!
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling sweepx, ij_ray, ik_ray=',2i4 )") &
& ij_ray, ik_ray
CALL sweepx( ij_ray, ik_ray )
IF ( l_write )  WRITE (nlog,"(' Returning from sweepx, ij_ray, ik_ray=',2i4 )") &
& ij_ray, ik_ray

!-----------------------------------------------------------------------
!
!                \\\\\ RETURN UPDATED VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Store variables after Lagrangian step
!-----------------------------------------------------------------------

r0l (nmin:nmax)         = r   (nmin:nmax)
dxl (nmin:nmax)         = dx  (nmin:nmax)
xal (nmin:nmax)         = xa  (nmin:nmax)
t0l (nmin:nmax)         = temp(nmin:nmax)
ei0l(nmin:nmax)         = ei  (nmin:nmax)
ye0l(nmin:nmax)         = ye  (nmin:nmax)
u0l (nmin:nmax)         = u   (nmin:nmax)

!-----------------------------------------------------------------------
!  Store temperature increment
!-----------------------------------------------------------------------

dtmpmn(jr_min:jr_max,1,ij_ray,ik_ray) = t0l(nmin:nmax) - t0i(nmin:nmax)

!-----------------------------------------------------------------------
!  Put variables advanced by Lagrangian hydrodynamics back in radhyd
!   arrays
!-----------------------------------------------------------------------

rhop  (imin:imax,ij_ray,ik_ray) = r     (nmin:nmax)
rho_l (imin:imax,ij_ray,ik_ray) = r     (nmin:nmax)
tp    (imin:imax,ij_ray,ik_ray) = temp  (nmin:nmax)
yep   (imin:imax,ij_ray,ik_ray) = ye    (nmin:nmax)
up    (imin:imax,ij_ray,ik_ray) = u     (nmin:nmax)
u_l   (imin:imax,ij_ray,ik_ray) = u     (nmin:nmax)
v_l   (imin:imax,ij_ray,ik_ray) = vp    (imin:imax,ij_ray,ik_ray)
w_l   (imin:imax,ij_ray,ik_ray) = wp    (imin:imax,ij_ray,ik_ray)
eip   (imin:imax,ij_ray,ik_ray) = ei    (nmin:nmax)
e_vp  (imin:imax,ij_ray,ik_ray) = e_v   (nmin:nmax)
x_el  (imin:imax,ij_ray,ik_ray) = xa    (nmin:nmax)
dx_cl (imin:imax,ij_ray,ik_ray) = dx    (nmin:nmax)
flat_x(imin:imax,ij_ray,ik_ray) = flat_s(nmin:nmax)

x_el(imax+1,ij_ray,ik_ray)      = xa(nmax+1)

u_e(imin:imax+1,ij_ray,ik_ray)  = u_edge(nmin:nmax+1)
x_cl(imin:imax,ij_ray,ik_ray)   = half * ( x_el(imin:imax,ij_ray,ik_ray) &
&                               + x_el(imin+1:imax+1,ij_ray,ik_ray) )

!-----------------------------------------------------------------------
!  Laod arrays for time_step determination
!-----------------------------------------------------------------------

drl (jr_min:jr_max)     = x_el(imin+1:imax+1,ij_ray,ik_ray)             &
&                       - x_el(imin:imax,ij_ray,ik_ray)
rhol(jr_min:jr_max)     = r0l(nmin:nmax)
rhoi(jr_min:jr_max)     = r0i(nmin:nmax)
tl  (jr_min:jr_max)     = t0l(nmin:nmax)
ti  (jr_min:jr_max)     = t0i(nmin:nmax)

!-----------------------------------------------------------------------
!  Hydro time step
!-----------------------------------------------------------------------

CALL hydro_x_sweep_time_step( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, drl, rhol, rhoi, tl, ti, nx, dt, jdt, dtime_hydro) 

RETURN
END SUBROUTINE evh1_x_lagr_inout
