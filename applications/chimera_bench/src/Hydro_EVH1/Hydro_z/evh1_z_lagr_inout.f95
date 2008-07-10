SUBROUTINE evh1_z_lagr_inout( kmin, kmax, nx, ny, nz, ki_ray, kj_ray, &
& ij_ray_dim, k_ray_dim, i_radial, nez, nnu, x_cf, z_ei, dz_ci, z_ci, &
& z_el, dz_cl, z_cl, dzphys_c, rho_z, t_z, ye_z, ei_z, u_z, v_z, w_z, &
& nu_str_cz, nu_str_ez, dtime, dt_z, jdt_z, rhobarp, flat_z, agr_z,   &
& grav_z_cz, cos_theta, sin_theta, e_nu_z )
!-----------------------------------------------------------------------
!
!    File:         evh1_z_lagr_inout
!    Module:       evh1_z_lagr_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/08/07
!
!    Purpose:
!      To receive the EVH1 arrays from radial_ray_module and azimuthal_ray_module,
!       execute the z-array Lagrangian hydrodynamics, and return to radial_ray_module
!       azimuthal_ray_module the updated variables.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  kmin         : minimum z-array index
!  kmax         : maximum z-array index
!  nx           : x-array extent
!  ny           : y-array extent
!  nz           : z-array extent
!  ki_ray       : x (radial) index of a specific z (azimuthal) ray
!  kj_ray       : y (angular) index of a specific z (azimuthal) ray
!  ij_ray_dim   : the number of radial zones on a processor before swapping with y
!  k_ray_dim    : the number of z-zones on a processor after swapping with z
!  i_radial     : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  x_cf         : final x grid zone midpoints
!  z_ei         : initial z grid zone faces
!  dz_ci        : initial z_ei(i+1) - z_eii)
!  z_ci         : initial z grid zone midpoints
!  dzphys_c     : physical coordinate difference of z_ei(i+1) - z_ei(i)
!  rho_z        : initial densities (g cm^{-3})
!  t_z          : initial temperatures (K)
!  ye_z         : initial electron fractions
!  ei_z         : initial internal energies (ergs g^{-1})
!  u_z          : initial velocity x-components (cm s^{-1})
!  v_z          : initial velocity y-components (cm s^{-1})
!  w_z          : initial velocity z-components (cm s^{-1})
!  nu_str_cz    : z-component of zone-cenered neutrino stress (dynes g^{-1})
!  nu_str_ez    : z-component of zone-edge neutrino stress (dynes g^{-1})
!  dtime        : time step (s)
!  rhobarp      : mean density as a function of radius (g cm^{-3})
!  cos_theta    : cos(theta(j))
!  sin_theta    : sin(theta(j))
!  agr_z        : lapse function
!  grav_zy_cz   : zone-centered z-component of gravitational acceleration 
!  e_nu_y       : neutrino energy density [ergs cm^{-3}]
!
!    Output arguments:
!
!  z_el         : updated z grid zone faces
!  dz_cl        : updated z_ei(i+1) - y_ei(i)
!  z_cl         : updated z grid zone midpoints
!  rho_z        : updated densities (g cm^{-3})
!  t_z          : updated temperatures (K)
!  ye_z         : updated electron fractions
!  ei_z         : updated internal energies (ergs g^{-1})
!  u_z          : updated velocity x-components (cm s^{-1})
!  flat_z       : variables indicating the presence of angular shocks
!  dt           : minimum hydro time step restrictions (s)
!  jdt          : zones causing minimum time step restrictions dt
!
!    Subprograms called:
!  sweepz                        : performs the z-Lagrangian hydro step
!  hydro_z_sweep_state_time_step : computes the z-sweep hydro time step restrictions due
!                                   to density and temperature change criteria
!      
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_y_module, evh1_sweep, evh1_sweep,
!  mgfld_remap_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, half, third

USE e_advct_module, ONLY : rhomin_z_eadvect
USE edit_module, ONLY : nlog
USE eos_snc_z_module, ONLY: aesv
USE evh1_global, ONLY: dt_g=>dt, i_radial_g=>i_radial
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, u, v, w, ei, ye, temp, xa0, &
& dx0, xa, dx, radius, nu_strc, nu_stre, entrop, rhobar, flat_s, lapse_c, &
& lapse_e, g_force_c, g_force_e, ctheta, stheta, p_nu, e_nu
USE mgfld_remap_module, ONLY : r0i, t0i, r0l, dxl, xal, t0l, ei0l, ye0l, u0l
USE parallel_module, ONLY : myid_y

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                                         :: kmin       ! minimum z-array index
INTEGER, INTENT(in)                                         :: kmax       ! maximum z-array index
INTEGER, INTENT(in)                                         :: nz         ! z-array extent

INTEGER, INTENT(in)                                         :: ki_ray     ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                                         :: kj_ray     ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                                         :: ij_ray_dim ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)                                         :: k_ray_dim  ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)                                         :: i_radial   ! the unshifted radial zone corresponding to ki_ray, kj_ray

INTEGER, INTENT(in)                                         :: nx         ! x-array extent
INTEGER, INTENT(in)                                         :: ny         ! y-array extent
INTEGER, INTENT(in)                                         :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)                                         :: nnu        ! neutrino flavor array extent

REAL(KIND=double), INTENT(in)                               :: dtime      ! time step

REAL(KIND=double), INTENT(in), DIMENSION(nx)                :: x_cf       ! x-coordinate zone center

REAL(KIND=double), INTENT(in), DIMENSION(nz+1)              :: z_ei       ! z-coordinate zone edge
REAL(KIND=double), INTENT(in), DIMENSION(nz)                :: dz_ci      ! z-coordinate zone thickness
REAL(KIND=double), INTENT(in), DIMENSION(nz)                :: z_ci       ! z-coordinate zone center

REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)      :: u_z        ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)      :: v_z        ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)      :: nu_str_cz  ! zone-centered neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(nz+1,ij_ray_dim,k_ray_dim)    :: nu_str_ez  ! zone-edged neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(nx)                :: rhobarp    ! mean density as a function of radius (g cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)      :: agr_z     ! lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)      :: grav_z_cz ! zone-centered z-component of gravitational acceleration
REAL(KIND=double), INTENT(in), DIMENSION(ny)                :: cos_theta  ! cos(theta(j))
REAL(KIND=double), INTENT(in), DIMENSION(ny)                :: sin_theta  ! sin(theta(j))
REAL(KIND=double), INTENT(in), DIMENSION(nz,ij_ray_dim,k_ray_dim)      :: e_nu_z    ! neutrino energy density [ergs cm^{-3}]

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(3,ij_ray_dim,k_ray_dim)                :: jdt_z      ! zone causing dt

REAL(KIND=double), INTENT(out), DIMENSION(nz+1,ij_ray_dim,k_ray_dim)   :: z_el       ! z-coordinate zone edge after Lagrangian step
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)     :: dz_cl      ! z-coordinate zone thickness after Lagrangian step
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)     :: z_cl       ! z-coordinate zone center after Lagrangian step
REAL(KIND=double), INTENT(out), DIMENSION(3,ij_ray_dim,k_ray_dim)      :: dt_z       ! minimum allowed time step
REAL(KIND=double), INTENT(out), DIMENSION(nz,ij_ray_dim,k_ray_dim)     :: flat_z     ! variables indicating the presence of angular shocks

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)   :: dzphys_c   ! zphys_e(i+1) - zphys_e(i)

REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)   :: rho_z      ! density (cm^{-3})
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)   :: t_z        ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)   :: ye_z       ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)      :: w_z        ! z-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)   :: ei_z       ! internal energy (ergs/g)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                                          :: l_write = .false.

REAL(KIND=double), DIMENSION(nz)                            :: rhol       ! density after Lagrangian step
REAL(KIND=double), DIMENSION(nz)                            :: rhoi       ! density before Lagrangian step
REAL(KIND=double), DIMENSION(nz)                            :: tl         ! temperature after Lagrangian step
REAL(KIND=double), DIMENSION(nz)                            :: ti         ! temperature before Lagrangian step

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!      \\\\\ PACK VARIABLES IN 1-D ARRAYS WITH GHOST ZONES /////
!
!-----------------------------------------------------------------------

sweep                    = 'z'

!-----------------------------------------------------------------------
!  Add six ghost zones to 1D sweep arrays
!-----------------------------------------------------------------------

nmin                     = kmin + 6
nmax                     = kmax + 6

nu_strc                  = zero
nu_stre                  = zero
g_force_c                = zero
g_force_e                = zero

!-----------------------------------------------------------------------
!  Put state variables into 1D arrays, padding with 6 ghost zones.
!-----------------------------------------------------------------------

r      (nmin:nmax)       = rho_z    (kmin:kmax,  kj_ray,ki_ray)
temp   (nmin:nmax)       = t_z      (kmin:kmax,  kj_ray,ki_ray)
ye     (nmin:nmax)       = ye_z     (kmin:kmax,  kj_ray,ki_ray)
u      (nmin:nmax)       = w_z      (kmin:kmax,  kj_ray,ki_ray)
v      (nmin:nmax)       = u_z      (kmin:kmax,  kj_ray,ki_ray)
w      (nmin:nmax)       = v_z      (kmin:kmax,  kj_ray,ki_ray)
ei     (nmin:nmax)       = ei_z     (kmin:kmax,  kj_ray,ki_ray)
entrop (nmin:nmax)       = aesv     (kmin:kmax,3,kj_ray,ki_ray)
nu_strc(nmin:nmax)       = nu_str_cz(kmin:kmax,  kj_ray,ki_ray)
nu_stre(nmin:nmax+1)     = nu_str_ez(kmin:kmax+1,kj_ray,ki_ray)
lapse_c  (nmin:nmax)     = agr_z    (kmin:kmax,  kj_ray,ki_ray)
lapse_e  (nmin:nmax)     = agr_z    (kmin:kmax,  kj_ray,ki_ray)
g_force_c(nmin:nmax)     = grav_z_cz(kmin:kmax,  kj_ray,ki_ray)
g_force_e(nmin+1:nmax)   = half * ( grav_z_cz(kmin:kmax-1,kj_ray,ki_ray) &
&                        + grav_z_cz(kmin+1:kmax,kj_ray,ki_ray) )

lapse_e(nmax+1)          = lapse_e(nmax)

!-----------------------------------------------------------------------
!  Store neutrino pressure and energy if rhobarp >= rhomin_y_eadvect
!-----------------------------------------------------------------------

p_nu                     = zero
e_nu                     = zero
IF ( rhobarp(i_radial) >= rhomin_z_eadvect ) THEN
  p_nu(nmin:nmax)        = third * e_nu_z(kmin:kmax,kj_ray,ki_ray)
  e_nu(nmin:nmax)        = e_nu_z(kmin:kmax,kj_ray,ki_ray)/r(nmin:nmax)
END IF ! rhobarp( i_radial ) >= rhomin_y_eadvect

!-----------------------------------------------------------------------
!  Load Grid Coordinates
!-----------------------------------------------------------------------

xa0(nmin:nmax+1)         = z_ei (kmin:kmax+1)
dx0(nmin:nmax)           = dz_ci(kmin:kmax)
xa (nmin:nmax+1)         = z_ei (kmin:kmax+1)
dx (nmin:nmax)           = dz_ci(kmin:kmax)

!-----------------------------------------------------------------------
!  Store initial values
!-----------------------------------------------------------------------

r0i (nmin:nmax)          = rho_z(kmin:kmax,kj_ray,ki_ray)
t0i (nmin:nmax)          = t_z  (kmin:kmax,kj_ray,ki_ray)

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
!  Cosine and sine of theta
!-----------------------------------------------------------------------

ctheta                   = cos_theta(kj_ray+myid_y*ij_ray_dim)
stheta                   = sin_theta(kj_ray+myid_y*ij_ray_dim)

!-----------------------------------------------------------------------
!  Radius (distance to z-axis)
!-----------------------------------------------------------------------

radius                   = x_cf(i_radial) * stheta

!-----------------------------------------------------------------------
!
!                   \\\\\ LAGRANGIAN HYDRO /////
!
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling sweepy, ki_ray, kj_ray=',2i4 )") &
& ki_ray, kj_ray
CALL sweepz( ki_ray, kj_ray )
IF ( l_write )  WRITE (nlog,"(' Returning from sweepy, ki_ray, kj_ray=', &
& 2i4 )") &
& ki_ray, kj_ray

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

rho_z (kmin:kmax,  kj_ray,ki_ray) = r     (nmin:nmax)
t_z   (kmin:kmax,  kj_ray,ki_ray) = temp  (nmin:nmax)
ye_z  (kmin:kmax,  kj_ray,ki_ray) = ye    (nmin:nmax)
w_z   (kmin:kmax,  kj_ray,ki_ray) = u     (nmin:nmax)
ei_z  (kmin:kmax,  kj_ray,ki_ray) = ei    (nmin:nmax)
z_el  (kmin:kmax+1,kj_ray,ki_ray) = xa    (nmin:nmax+1)
dz_cl (kmin:kmax,  kj_ray,ki_ray) = dx    (nmin:nmax)
flat_z(kmin:kmax,  kj_ray,ki_ray) = flat_s(nmin:nmax)

z_cl(kmin:kmax,kj_ray,ki_ray)     = half * ( z_el(kmin  :kmax  ,kj_ray,ki_ray) &
&                                 +          z_el(kmin+1:kmax+1,kj_ray,ki_ray) )

!-----------------------------------------------------------------------
!  Laod arrays for time_step determination
!-----------------------------------------------------------------------

rhol(kmin:kmax)          = r0l(nmin:nmax)
rhoi(kmin:kmax)          = r0i(nmin:nmax)
tl  (kmin:kmax)          = t0l(nmin:nmax)
ti  (kmin:kmax)          = t0i(nmin:nmax)

!-----------------------------------------------------------------------
!  Hydro time step
!-----------------------------------------------------------------------

CALL hydro_z_sweep_state_time_step( kmin, kmax, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, rhol, rhoi, tl, ti, nz, dt_z, jdt_z ) 

RETURN
END SUBROUTINE evh1_z_lagr_inout
