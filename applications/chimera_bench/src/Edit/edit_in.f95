SUBROUTINE edit_in( imin, imax, jmin, jmax, kmin, kmax, nx, ny, nz,   &
& ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, nez, nnu, nnc, rhorp, rhop, &
& trp, tp, yerp, yep, rhobarp, rp, up, vp, wp, psi0p, psi1p, dtime,   &
& time_elapsed, i_editp, cycle_number, xnp, be_nuc_repp, a_nuc_repp,  &
& z_nuc_repp, nu_strep, nu_strcp, agr_e, agr_c, grav_x_e, grav_x_c,   &
& grav_y_c, grav_pot_c, mass_ns, vel_ns, nsep, nedc_in, nedmi_in,     &
& nedma_in, nedh_in, nedps_in, nedu_in, nedy_in, nedsc_in, nedn_in,   &
& nedng_in, d_omega, first )
!-----------------------------------------------------------------------
!
!    File:         edit_in
!    Module:       edit_in
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To load variables for the mgfld edits.
!
!    Subprograms called:
!
!  eqstz_x      : updates equation of state quantities
!  gammaz_x     : updates equation of state gammas
!  flux         : computes the neutrino flux
!  time_bounce  : sets time of t_bounce
!  edit_exec    : performs the edit
!
!    Input arguments:
!
!  imin         : minimum x-array index for the edit
!  imax         : maximum x-array index for the edit
!  jmin         : minimum y-array index for the edit
!  jmax         : maximum y-array index for the edit
!  kmin         : minimum z-array index for the edit
!  kmax         : maximum z-array index for the edit
!  nx           : x_array extent
!  ny           : y_array extent
!  nz           : z_array extent
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  nnc          : neutrino abundance array extent
!  rhop         : density (cm^{-3})
!  tp           : temperature (MeV)
!  yep          : electron fraction
!  rhobarp      : mean density (cm^{-3})
!  rp           : radial zone radii (cm)
!  up           : radial velocity of zone (cm)
!  vp           : angular velocity of zone (cm)
!  wp           : azimuthal velocity of zone (cm)
!  psi0p        : zeroth angular moment of the NDS
!  dtime_h      : hydro time step
!  dtime        : time step used at current cycle
!  time_elapsed : elapsed time
!  i_editp       : edit parameter
!                 O - edits are executed on the basis of the counter values
!                 1 - brief model edit
!                 2 - full model edit
!                 3 - full model edit plus differential neutrino edit
!                 4 - full model edit plus differential plus integrated neutrino edit
!  cycle_number : cycle number
!  xnp          : initial mass fractions
!  be_nuc_repp  : binding energy of mean heavy nucleus
!  a_nuc_repp   : mass number of mean heavy nucleus
!  z_nuc_repp   : charge number of mean heavy nucleus
!  nu_strep     : total zone-edged neutrino stress (dynes g^{-1})
!  nu_strcp     : total zone-centered neutrino stress (dynes g^{-1})
!  nesp         : nuclear statistical equilibrium flag
!  agr_e        : unshifted zone-edged lapse function
!  agr_c        : unshifted zone-centered lapse function
!  grav_x_e     : unshifted zone-edged gravtational x-acceleration
!  grav_x_c     : unshifted zone-centered gravtational x-acceleration
!  grav_y_c     : unshifted zone-centered gravtational y-acceleration
!  grav_pot_c   : unshifted zone-centered gravitational potential
!  mass_ns      : mass of the neutron star
!  vel_ns       : velocity of the neutron star
!  nedc_in      : the number of cycles since the last implementation of subsection i of subroutine editc
!  nedmi_in     : the number of cycles since the last implementation of subsection i of subroutine editmi
!  nedma_in     : the number of cycles since the last implementation of subsection i of subroutine editma
!  nedh_in      : the number of cycles since the last implementation of subsection i of subroutine edith
!  nedps_in     : the number of cycles since the last implementation of subsection i of subroutine editps
!  nedu_in      : the number of cycles since the last implementation of subsection i of subroutine editu
!  nedy_in      : the number of cycles since the last implementation of subsection i of subroutine edity
!  nedsc_in     : the number of cycles since the last implementation of subsection i of subroutine editsc
!  nedn_in      : the number of cycles since the last implementation of subroutine editn for neutrinos
!                  of type n.
!  nedng_in     : the number of cycles since the last implementation of subsection i of subroutine 
!                  editng for neutrinos of type n.
!  d_omega      : solid angles subtended by radial rays
!  first        : initial call flag
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, array_module, numerical_module, physcnst_module
!  cycle_module, edit_module, eos_snc_x_module, mdl_cnfg_module, 
!  nucbrn_module, nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : kmev_inv, msolar

USE cycle_module, ONLY : ncycle
USE edit_module, ONLY : nedc, nedmi, nedma, nedh, nedps, nedu, nedy, nedsc, &
& nedn, nedng, gstrss, gstrss_cx, gstrss_cy, g_pot, mass_ns_e=>mass_ns, &
& vel_ns_e=>vel_ns, d_omega_e=>d_omega, rhobar
USE eos_snc_x_module, ONLY : xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, nse
USE mdl_cnfg_module, ONLY : rhor, rho, tr, t, yer, ye, dr, r, u, v, w, dmrst, &
& rstmss, agr, agrh, jr_min, jr_max
USE nucbrn_module, ONLY : xn_n=>xn, be_nuc_rep_n=>be_nuc_rep, a_nuc_rep_n=>a_nuc_rep, &
& z_nuc_rep_n=>z_nuc_rep, fescrn, fascrn, uburn_n=>uburn, nse_n=>nse
USE nu_dist_module, ONLY : psi0, psi1, vol, nu_str_ex, nu_str_cx
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY: dtime_hydro, dtnph, time

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

LOGICAL, INTENT(in)                               :: first         ! initial model print flag

INTEGER, INTENT(in)                               :: imin          ! minimum x-array index for the edit
INTEGER, INTENT(in)                               :: imax          ! maximum x-array index for the edit
INTEGER, INTENT(in)                               :: jmin          ! minimum y-array index for the edit
INTEGER, INTENT(in)                               :: jmax          ! maximum y-array index for the edit
INTEGER, INTENT(in)                               :: kmin          ! minimum z-array index for the edit
INTEGER, INTENT(in)                               :: kmax          ! maximum z-array index for the edit
INTEGER, INTENT(in)                               :: nx            ! x-array extent
INTEGER, INTENT(in)                               :: ny            ! y-array extent
INTEGER, INTENT(in)                               :: nz            ! z-array extent

INTEGER, INTENT(in)                               :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)                               :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)                               :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                               :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                               :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)                               :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)                               :: nnc           ! composition array extent

INTEGER, INTENT(in)                               :: i_editp       ! edit parameter
INTEGER, INTENT(in)                               :: cycle_number  ! cycle number

INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedc_in       ! edit counter for subroutine editc
INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedmi_in      ! edit counter for subroutine editmi
INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedma_in      ! edit counter for subroutine editma
INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedh_in       ! edit counter for subroutine edith
INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedps_in      ! edit counter for subroutine editps
INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedu_in       ! edit counter for subroutine editu
INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedy_in       ! edit counter for subroutine edity
INTEGER, INTENT(in), DIMENSION(20,ij_ray_dim,ik_ray_dim)      :: nedsc_in      ! edit counter for subroutine editsc

INTEGER, INTENT(in), DIMENSION(nnu,ij_ray_dim,ik_ray_dim)     :: nedn_in       ! edit counter for subroutine editn
INTEGER, INTENT(in), DIMENSION(100,nnu,ij_ray_dim,ik_ray_dim) :: nedng_in      ! edit counter for subroutine editng

INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: nsep          ! nuclear sttimintical equilibrium flag

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rhorp        ! density (cm^{-3}) before hydro step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rhop         ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: trp          ! temperature (MeV) before hydro step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: tp           ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: yerp         ! electron fraction before hydro step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: yep          ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx)                               :: rhobarp      ! mean density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                             :: rp           ! radial zone radii (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: up           ! radial velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: vp           ! angular velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: wp           ! azimuthal velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: nu_strep     ! neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: nu_strcp     ! neutrino stress
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)       :: agr_e        ! unshifted zone-edged lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: agr_c        ! unshifted zone-entered lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)       :: grav_x_e     ! unshifted zone-edged gravtational x-acceleration
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: grav_x_c     ! unshifted zone-centered gravtational x-acceleration
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: grav_y_c     ! unshifted zone-centered gravtational y-acceleration
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: grav_pot_c   ! unshifted zone-centered gravitational potential
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0p        ! zeroth angular moment of the NDS
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1p        ! first angular moment of the NDS

REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim)     :: xnp          ! mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: be_nuc_repp  ! nuclear binding energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: a_nuc_repp   ! nuclear mass number
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: z_nuc_repp   ! nuclear charge number

REAL(KIND=double), INTENT(in)                     :: mass_ns      ! mass of the neutron star (g)
REAL(KIND=double), INTENT(in)                     :: vel_ns       ! velocity of the neutron star (cm s^{-1})

REAL(KIND=double), INTENT(in)                     :: dtime        ! time step used at current cycle
REAL(KIND=double), INTENT(in)                     :: time_elapsed ! elapsed time

REAL(KIND=double), INTENT(in), DIMENSION(ny,nz)   :: d_omega      ! solid angles subtended by radial rays

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                           :: jr_maxp       ! jr_max + 1
INTEGER                                           :: j             ! radial zone index
INTEGER                                           :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!            \\\\\ LOAD VARIABLES INTO MGFLD ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer integer scalar parameters
!-----------------------------------------------------------------------

jr_min                      = imin + 1
jr_maxp                     = jr_max + 1
jr_max                      = imax + 1
jr_maxp                     = jr_max + 1
ncycle                      = cycle_number

!-----------------------------------------------------------------------
!  Transfer edit counters
!-----------------------------------------------------------------------

nedc (1:20)                 = nedc_in (1:20,ij_ray,ik_ray)
nedmi(1:20)                 = nedmi_in(1:20,ij_ray,ik_ray)
nedma(1:20)                 = nedma_in(1:20,ij_ray,ik_ray)
nedh (1:20)                 = nedh_in (1:20,ij_ray,ik_ray)
nedps(1:20)                 = nedps_in(1:20,ij_ray,ik_ray)
nedu (1:20)                 = nedu_in (1:20,ij_ray,ik_ray)
nedy (1:20)                 = nedy_in (1:20,ij_ray,ik_ray)
nedsc(1:20)                 = nedsc_in(1:20,ij_ray,ik_ray)

nedn(:)                     = nedn_in(:,    ij_ray,ik_ray)
nedng(1:40,:)               = nedng_in(1:40,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to index shifted mgfld
!   arrays
!-----------------------------------------------------------------------

rhor  (jr_min:jr_max)       = rhorp     (imin:imax,ij_ray,ik_ray)
rho   (jr_min:jr_max)       = rhop      (imin:imax,ij_ray,ik_ray)
tr    (jr_min:jr_max)       = trp       (imin:imax,ij_ray,ik_ray)
t     (jr_min:jr_max)       = tp        (imin:imax,ij_ray,ik_ray)
yer   (jr_min:jr_max)       = yerp      (imin:imax,ij_ray,ik_ray)
ye    (jr_min:jr_max)       = yep       (imin:imax,ij_ray,ik_ray)
u     (jr_min:jr_max)       = up        (imin:imax,ij_ray,ik_ray)
v     (jr_min:jr_max)       = vp        (imin:imax,ij_ray,ik_ray)
w     (jr_min:jr_max)       = wp        (imin:imax,ij_ray,ik_ray)
rhobar(jr_min:jr_max)       = rhobarp   (imin:imax)
agrh  (jr_min:jr_max)       = agr_c     (imin:imax,ij_ray,ik_ray)
g_pot (jr_min:jr_max)       = grav_pot_c(imin:imax,ij_ray,ik_ray)
gstrss_cx(jr_min:jr_max)    = grav_x_c  (imin:imax,ij_ray,ik_ray)
gstrss_cy(jr_min:jr_max)    = grav_y_c  (imin:imax,ij_ray,ik_ray)

nu_str_cx(jr_min:jr_max)    = nu_strcp(imin:imax,ij_ray,ik_ray)
psi0(jr_min:jr_max,:,:)     = psi0p(imin:imax,:,:,ij_ray,ik_ray)
psi1(imin:imax+1,:,:)       = psi1p(imin:imax+1,:,:,ij_ray,ik_ray)

rho(jr_maxp)                = rho(jr_max)
t  (jr_maxp)                = t  (jr_max)
ye (jr_maxp)                = ye (jr_max)

!-----------------------------------------------------------------------
!  Transfer zone-edgeed independent variables to mgfld arrays
!-----------------------------------------------------------------------

r        (imin:imax+1)      = rp      (imin:imax+1)
agr      (imin:imax+1)      = agr_e   (imin:imax+1,ij_ray,ik_ray)
gstrss   (imin:imax+1)      = grav_x_e(imin:imax+1,ij_ray,ik_ray)
nu_str_ex(imin:imax+1)      = nu_strep(imin:imax+1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer zone-centered composition variables
!-----------------------------------------------------------------------

nse(jr_min:jr_max,ij_ray,ik_ray) = nsep  (imin:imax,ij_ray,ik_ray)
be_nuc_rep(jr_min:jr_max)   = be_nuc_repp(imin:imax,ij_ray,ik_ray)
a_nuc_rep (jr_min:jr_max)   = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep (jr_min:jr_max)   = z_nuc_repp (imin:imax,ij_ray,ik_ray)

xn(jr_min:jr_max,:)         = xnp(imin:imax,:,ij_ray,ik_ray)

nse_n   (jr_min:jr_max)     = nsep   (imin:imax,ij_ray,ik_ray)
be_nuc_rep_n(jr_min:jr_max) = be_nuc_repp(imin:imax,ij_ray,ik_ray)
a_nuc_rep_n (jr_min:jr_max) = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep_n (jr_min:jr_max) = z_nuc_repp (imin:imax,ij_ray,ik_ray)

xn_n(jr_min:jr_max,:)       = xnp(imin:imax,:,ij_ray,ik_ray)
 
!-----------------------------------------------------------------------
!  Derived zone-centered and zone-edged variables
!-----------------------------------------------------------------------

!........radial zone thickness

dr(jr_min:jr_max)           = r(jr_min:jr_max) - r(jr_min-1:jr_max-1)

!........Newtonian rest masses

IF ( r(1) == zero ) THEN
  dmrst(1)                  = zero
ELSE
  dmrst(1)                  = frpith * r(1)**3 * rhobar(2)
END IF

rstmss(1)                   = dmrst(1)

DO j    = jr_min,jr_max
  vol(j)                    = frpith * ( r(j) * ( r(j) + r(j-1) ) + r(j-1) * r(j-1) ) * dr(j)
  dmrst(j)                  = vol(j) * rhobar(j)
  rstmss(j)                 = rstmss(j-1) + dmrst(j)
END DO

!-----------------------------------------------------------------------
!  Neutron star parameters
!-----------------------------------------------------------------------

mass_ns_e                   = mass_ns/msolar
vel_ns_e                    = vel_ns/1.d+05

!-----------------------------------------------------------------------
!  Solid angles
!-----------------------------------------------------------------------

d_omega_e                   = d_omega

!-----------------------------------------------------------------------
!  Set timestep
!-----------------------------------------------------------------------

dtnph                       = dtime
time                        = time_elapsed

!-----------------------------------------------------------------------
!  Update EOS variables
!-----------------------------------------------------------------------

CALL eqstz_x( jr_min, jr_maxp, rho, t, ye, ij_ray, ik_ray )
CALL gammaz_x( jr_min, jr_max, rho, t, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Compute the neutrino fluxes
!-----------------------------------------------------------------------

DO n    = 1,nnu
  CALL flux( jr_min, jr_max, n )
END DO

!-----------------------------------------------------------------------
!  Set time of t_bounce
!-----------------------------------------------------------------------

CALL time_bounce( rhobar, nx )

!-----------------------------------------------------------------------
!
!                  \\\\\ EXECUTE THE EDITS /////
!
!-----------------------------------------------------------------------

CALL edit_exec( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& jmin, jmax, kmin, kmax, i_editp, first, nx, ny, nz, nez, nnu  )

RETURN
END SUBROUTINE edit_in
