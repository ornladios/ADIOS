SUBROUTINE evh1_y_lagr_inout( jmin, jmax, nx, ny, nz, ji_ray, jk_ray, &
& j_ray_dim, ik_ray_dim, i_radial, nez, nnu, x_cf, y_ei, dy_ci, y_ci, &
& y_el, dy_cl, y_cl, dyphys_c, rho_y, t_y, ye_y, ei_y, u_y, v_y, w_y, &
& nu_str_cy, nu_str_ey, dtime, dt_y, jdt_y, rhobarp, flat_y, agr_y,   &
& grav_y_cy, e_nu_y  )
!-----------------------------------------------------------------------
!
!    File:         evh1_y_lagr_inout
!    Module:       evh1_y_lagr_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To receive the EVH1 arrays from radial_ray_module and angular_ray_module,
!       execute the y-array Lagrangian hydrodynamics, and return to radial_ray_module
!       angular_ray_module the updated variables.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  jmin         : minimum y-array index
!  jmax         : maximum y-array index
!  nx           : x-array extent
!  ny           : y-array extent
!  nz           : z-array extent
!  ji_ray       : x (radial) index of a specific y (angular) ray
!  jk_ray       : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim    : the number of radial zones on a processor after swapping with y
!  ik_ray_dim   : the number of z-zones on a processor before swapping with z
!  i_radial     : the unshifted radial zone (angular ray) corresponding to ji_ray,jk_ray
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  x_cf         : final x grid zone midpoints
!  y_ei         : initial y grid zone faces
!  dy_ci        : initial y_ei(i+1) - y_ei(i)
!  y_ci         : initial y grid zone midpoints
!  dyphys_c     : physical coordinate difference of y_ei(i_1) - uei(i)
!  rho_y        : initial densities (g cm^{-3})
!  t_y          : initial temperatures (K)
!  ye_y         : initial electron fractions
!  ei_y         : initial internal energies (ergs g^{-1})
!  u_y          : initial velocity x-components (cm s^{-1})
!  v_y          : initial velocity y-components (cm s^{-1})
!  w_y          : initial velocity z-components (cm s^{-1})
!  nu_str_cy    : y-component of zone-cenered neutrino stress (dynes g^{-1})
!  nu_str_ey    : y-component of zone-edge neutrino stress (dynes g^{-1})
!  dtime        : time step (s)
!  rhobarp      : mean density as a function of radius (g cm^{-3})
!  agr_y        : lapse function
!  grav_y_cy    : zone-centered y-component of gravitational acceleration 
!  e_nu_y       : neutrino energy density [ergs cm^{-3}]
!
!    Output arguments:
!
!  y_el         : updated y grid zone faces
!  dy_cl        : updated y_ei(i+1) - y_ei(i)
!  y_cl         : updated y grid zone midpoints
!  rho_y        : updated densities (g cm^{-3})
!  t_y          : updated temperatures (K)
!  ye_y         : updated electron fractions
!  ei_y         : updated internal energies (ergs g^{-1})
!  u_y          : updated velocity x-components (cm s^{-1})
!  flat_y       : variables indicating the presence of angular shocks
!  dt           : minimum hydro time step restrictions (s)
!  jdt          : zones causing minimum time step restrictions dt
!
!    Subprograms called:
!  sweepy                        : performs the y-Lagrangian hydro step
!  hydro_y_sweep_state_time_step : computes the y-sweep hydro time step restrictions due
!                                   to density and temperature change criteria
!      
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, edit_module, eos_snc_y_module, evh1_sweep,
!  evh1_sweep, mgfld_remap_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, half, third

USE e_advct_module, ONLY : rhomin_y_eadvect
USE edit_module, ONLY : nlog
USE eos_snc_y_module, ONLY: aesv
USE evh1_global, ONLY: dt_g=>dt, i_radial_g=>i_radial
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, u, v, w, ei, ye, temp, xa0, &
& dx0, xa, dx, radius, nu_strc, nu_stre, entrop, rhobar, flat_s, lapse_c, &
& lapse_e, g_force_c, g_force_e, p_nu, e_nu
USE mgfld_remap_module, ONLY : r0i, t0i, r0l, dxl, xal, t0l, ei0l, ye0l, u0l
USE parallel_module, ONLY : myid

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                         :: jmin       ! inner physical x-zone
INTEGER, INTENT(in)                                         :: jmax       ! outer physical x-zone
INTEGER, INTENT(in)                                         :: ny         ! y-array extent

INTEGER, INTENT(in)                                         :: ji_ray     ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)                                         :: jk_ray     ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)                                         :: j_ray_dim  ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)                                         :: ik_ray_dim ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)                                         :: i_radial   ! the unshifted radial zone corresponding to ji_ray,jk_ray

INTEGER, INTENT(in)                                         :: nx         ! x-array extent
INTEGER, INTENT(in)                                         :: nz         ! z-array extent
INTEGER, INTENT(in)                                         :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)                                         :: nnu        ! neutrino flavor array extent

REAL(KIND=double), INTENT(in)                               :: dtime      ! time step

REAL(KIND=double), INTENT(in), DIMENSION(nx)                :: x_cf       ! x-coordinate zone center

REAL(KIND=double), INTENT(in), DIMENSION(ny+1)              :: y_ei       ! y-coordinate zone edge
REAL(KIND=double), INTENT(in), DIMENSION(ny)                :: dy_ci      ! y-coordinate zone thickness
REAL(KIND=double), INTENT(in), DIMENSION(ny)                :: y_ci       ! y-coordinate zone center

REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)      :: u_y        ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)      :: w_y        ! z-velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)      :: nu_str_cy  ! zone-centered neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(ny+1,j_ray_dim,ik_ray_dim)    :: nu_str_ey  ! zone-edged neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(nx)                :: rhobarp    ! mean density as a function of radius (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)      :: agr_y     ! lapse function
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)      :: grav_y_cy ! zone-centered y-component of gravitational acceleration
REAL(KIND=double), INTENT(in), DIMENSION(ny,j_ray_dim,ik_ray_dim)      :: e_nu_y    ! neutrino energy density [ergs cm^{-3}]

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(3,j_ray_dim,ik_ray_dim)                :: jdt_y      ! zone causing dt

REAL(KIND=double), INTENT(out), DIMENSION(ny+1,j_ray_dim,ik_ray_dim)   :: y_el       ! y-coordinate zone edge after Lagrangian step
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)     :: dy_cl      ! y-coordinate zone thickness after Lagrangian step
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)     :: y_cl       ! y-coordinate zone center after Lagrangian step
REAL(KIND=double), INTENT(out), DIMENSION(3,j_ray_dim,ik_ray_dim)      :: dt_y       ! minimum allowed time step
REAL(KIND=double), INTENT(out), DIMENSION(ny,j_ray_dim,ik_ray_dim)     :: flat_y     ! variables indicating the presence of angular shocks

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)   :: dyphys_c   ! physical coordinate difference of y_ei(i_1) - uei(i)

REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)   :: rho_y      ! density (cm^{-3})
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)   :: t_y        ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)   :: ye_y       ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)   :: v_y        ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)   :: ei_y       ! internal energy (ergs/g)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                                          :: l_write = .false.

REAL(KIND=double), DIMENSION(ny)                            :: rhol       ! density after Lagrangian step
REAL(KIND=double), DIMENSION(ny)                            :: rhoi       ! density before Lagrangian step
REAL(KIND=double), DIMENSION(ny)                            :: tl         ! temperature after Lagrangian step
REAL(KIND=double), DIMENSION(ny)                            :: ti         ! temperature before Lagrangian step

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!      \\\\\ PACK VARIABLES IN 1-D ARRAYS WITH GHOST ZONES /////
!
!-----------------------------------------------------------------------

sweep                    = 'y'

!-----------------------------------------------------------------------
!  Add six ghost zones to 1D sweep arrays
!-----------------------------------------------------------------------

nmin                     = jmin + 6
nmax                     = jmax + 6

nu_strc                  = zero
nu_stre                  = zero
g_force_c                = zero
g_force_e                = zero

!-----------------------------------------------------------------------
!  Put state variables into 1D arrays, padding with 6 ghost zones.
!-----------------------------------------------------------------------

r        (nmin:nmax)     = rho_y    (jmin:jmax,  ji_ray,jk_ray)
temp     (nmin:nmax)     = t_y      (jmin:jmax,  ji_ray,jk_ray)
ye       (nmin:nmax)     = ye_y     (jmin:jmax,  ji_ray,jk_ray)
u        (nmin:nmax)     = v_y      (jmin:jmax,  ji_ray,jk_ray)
v        (nmin:nmax)     = w_y      (jmin:jmax,  ji_ray,jk_ray)
w        (nmin:nmax)     = u_y      (jmin:jmax,  ji_ray,jk_ray)
ei       (nmin:nmax)     = ei_y     (jmin:jmax,  ji_ray,jk_ray)
entrop   (nmin:nmax)     = aesv     (jmin:jmax,3,ji_ray,jk_ray)
nu_strc  (nmin:nmax)     = nu_str_cy(jmin:jmax,  ji_ray,jk_ray)
nu_stre  (nmin:nmax+1)   = nu_str_ey(jmin:jmax+1,ji_ray,jk_ray)
lapse_c  (nmin:nmax)     = agr_y    (jmin:jmax,  ji_ray,jk_ray)
lapse_e  (nmin:nmax)     = agr_y    (jmin:jmax,  ji_ray,jk_ray)
g_force_c(nmin:nmax)     = grav_y_cy(jmin:jmax,  ji_ray,jk_ray)
g_force_e(nmin+1:nmax)   = half * ( grav_y_cy(jmin:jmax-1,ji_ray,jk_ray) &
&                        + grav_y_cy(jmin+1:jmax,ji_ray,jk_ray) )

lapse_e(nmax+1)          = lapse_e(nmax)

!-----------------------------------------------------------------------
!  Store neutrino pressure and energy if rhobarp >= rhomin_y_eadvect
!-----------------------------------------------------------------------

p_nu                     = zero
e_nu                     = zero
IF ( rhobarp(i_radial) >= rhomin_y_eadvect ) THEN
  p_nu(nmin:nmax)        = third * e_nu_y(jmin:jmax,ji_ray,jk_ray)
  e_nu(nmin:nmax)        = e_nu_y(jmin:jmax,ji_ray,jk_ray)/r(nmin:nmax)
END IF ! rhobarp( i_radial ) >= rhomin_y_eadvect

!-----------------------------------------------------------------------
!  Load Grid Coordinates
!-----------------------------------------------------------------------

xa0(nmin:nmax+1)         = y_ei (jmin:jmax+1)
dx0(nmin:nmax)           = dy_ci(jmin:jmax)
xa (nmin:nmax+1)         = y_ei (jmin:jmax+1)
dx (nmin:nmax)           = dy_ci(jmin:jmax)

!-----------------------------------------------------------------------
!  Store initial values
!-----------------------------------------------------------------------

r0i (nmin:nmax)          = rho_y(jmin:jmax,ji_ray,jk_ray)
t0i (nmin:nmax)          = t_y  (jmin:jmax,ji_ray,jk_ray)

radius                   = x_cf(i_radial)

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

dt_g                     = dtime
i_radial_g               = i_radial

!-----------------------------------------------------------------------
!  Load rhobar
!-----------------------------------------------------------------------

rhobar(7:nx+6)           = rhobarp(1:nx)

!-----------------------------------------------------------------------
!  Compute energies
!-----------------------------------------------------------------------

!CALL etotal(.false.)

!-----------------------------------------------------------------------
!
!                   \\\\\ LAGRANGIAN HYDRO /////
!
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling sweepy, ji_ray, jk_ray=',2i4 )") &
& ji_ray, jk_ray
CALL sweepy( ji_ray, jk_ray )
IF ( l_write )  WRITE (nlog,"(' Returning from sweepy, ji_ray, jk_ray=', &
& 2i4 )") &
& ji_ray, jk_ray

!-----------------------------------------------------------------------
!
!                \\\\\ RETURN UPDATED VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Store variables after Lagrangian step
!-----------------------------------------------------------------------

r0l (nmin:nmax)          = r   (nmin:nmax)
dxl (nmin:nmax)          = dx  (nmin:nmax)
xal (nmin:nmax)          = xa  (nmin:nmax)
t0l (nmin:nmax)          = temp(nmin:nmax)
ei0l(nmin:nmax)          = ei  (nmin:nmax)
ye0l(nmin:nmax)          = ye  (nmin:nmax)
u0l (nmin:nmax)          = u   (nmin:nmax)

!-----------------------------------------------------------------------
!  Put variables advanced by Lagrangian hydrodynamics back in radhyd
!   arrays
!-----------------------------------------------------------------------

rho_y (jmin:jmax,  ji_ray,jk_ray) = r     (nmin:nmax)
t_y   (jmin:jmax,  ji_ray,jk_ray) = temp  (nmin:nmax)
ye_y  (jmin:jmax,  ji_ray,jk_ray) = ye    (nmin:nmax)
v_y   (jmin:jmax,  ji_ray,jk_ray) = u     (nmin:nmax)
ei_y  (jmin:jmax,  ji_ray,jk_ray) = ei    (nmin:nmax)
y_el  (jmin:jmax+1,ji_ray,jk_ray) = xa    (nmin:nmax+1)
dy_cl (jmin:jmax,  ji_ray,jk_ray) = dx    (nmin:nmax)
flat_y(jmin:jmax,  ji_ray,jk_ray) = flat_s(nmin:nmax)

y_cl(jmin:jmax,ji_ray,jk_ray)     = half * ( y_el(jmin  :jmax  ,ji_ray,jk_ray) &
&                                 +          y_el(jmin+1:jmax+1,ji_ray,jk_ray) )

!-----------------------------------------------------------------------
!  Laod arrays for time_step determination
!-----------------------------------------------------------------------

rhol(jmin:jmax)          = r0l(nmin:nmax)
rhoi(jmin:jmax)          = r0i(nmin:nmax)
tl  (jmin:jmax)          = t0l(nmin:nmax)
ti  (jmin:jmax)          = t0i(nmin:nmax)

!-----------------------------------------------------------------------
!  Hydro time step
!-----------------------------------------------------------------------

CALL hydro_y_sweep_state_time_step( jmin, jmax, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, rhol, rhoi, tl, ti, ny, dt_y, jdt_y ) 

RETURN
END SUBROUTINE evh1_y_lagr_inout
