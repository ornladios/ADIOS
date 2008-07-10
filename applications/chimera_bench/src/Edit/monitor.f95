SUBROUTINE monitor( imin, imax, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& nx, nez, nnu, nnc, rhorp, rhop, trp, tp, yerp, yep, rhobarp, rp, up, &
& psi0p, psi1p, dtime, time_elapsed, cycle_number, gtot_pot_c, xnp, &
& be_nuc_repp, a_nuc_repp, z_nuc_repp, nsep, u_ge_tot, u_ie_tot, u_ke_tot, &
& u_ne_tot, u_ke_x_tot, u_ke_y_tot, u_ke_z_tot, e_radp, radtot, u_tot, &
& elecn, elec_radp, nnucr_enu, nnucr_enubar, nnurad_enu, nnurad_enubar, &
& totlpn )
!-----------------------------------------------------------------------
!
!    File:         monitor
!    Module:       monitor
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To calculate quantities that should be conserved during the
!          course of a simulation.
!
!    Subprograms called:
!  eqstz_x      : interpolates quantities in a local EOS table
!  gammaz_x     : computes the EOS gammas
!  flux         : computes neutrino fluxes
!  cnfg_energy  : computes the total energy and lepton number of the model
!
!    Input arguments:
!
!  imin         : inner physical x-zone center
!  imax         : outer physical x-zone center
!  ij_ray       : j-index of a radial ray
!  ik_ray       : k-index of a radial ray
!  ij_ray_dim   : number of y-zones on a processor before swapping
!  ik_ray_dim   : number of z-zones on a processor before swapping
!  nx           : x_array extent
!  nez          : neutrino energy array extent
!  nnu          : neutrino flavor array extent
!  nnc          : neutrino abundance array extent
!  rhop         : density (cm^{-3})
!  tp           : temperature (MeV)
!  yep          : electron fraction
!  rhobarp      : mean density (cm^{-3})
!  rp           : radial zone radii (cm)
!  up           : radial velocity of zone (cm)
!  psi0p        : zeroth angular moment of the NDS
!  dtime_h      : hydro time step
!  dtime        : time step used at current cycle
!  time_elapsed : elapsed time
!  cycle_number : cycle number
!  gtot_pot_c   : unshifted zone-centered gravitational potential energy
!  xnp          : initial mass fractions
!  be_nuc_repp  : binding energy of mean heavy nucleus
!  a_nuc_repp   : mass number of mean heavy nucleus
!  z_nuc_repp   : charge number of mean heavy nucleus
!  nesp         : nuclear statistical equilibrium flag
!
!    Output arguments:
!
!  u_ge_tot     : total NT gravitational potential energy of i_ray (ergs)
!  u_ie_tot     : total internal energy of i_ray (ergs)
!  u_ke_tot     : total NT kinetic energy of i_ray (ergs)
!  u_ne_tot     : total NT neutrino energy of i_ray (ergs)
!  u_ke_x_tot   : total NT x-kinetic energy of i_ray (ergs)
!  u_ke_y_tot   : total NT y-kinetic energy of i_ray (ergs)
!  u_ke_z_tot   : total NT z-kinetic energy of i_ray (ergs)
!  e_radp       : the cumulative material energy entering i_ray grid (ergs)
!  radtot       : energy radiated by neutrinos in i_ray (ergs)
!  u_tot        : total energy of i_ray (ergs)
!  elecn        : total electron number of i_ray
!  elec_rad     : the net number of electrons advected into the ij_ray,ik_ray grid
!  nnucr        : the net number of n-type neutrinos currently residing in ij_ray,ik_ray grid
!  nnurad       : the cumulative number of n-type neutrinos emitted by the ij_ray,ik_ray grid
!  totlpn       : total number of leptons emitted or residing in i_ray
!  u_ge_tot     : total NT gravitational potential energy
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  cycle_module, eos_snc_x_module, edit_module, mdl_cnfg_module,
!  nucbrn_module, nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : kmev_inv

USE cycle_module, ONLY : ncycle
USE eos_snc_x_module, ONLY : xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, nse
USE edit_module, ONLY : g_pot
USE mdl_cnfg_module, ONLY : rhor, rho, tr, t, yer, ye, dr, r, u, dmrst, rstmss, &
& jr_min, jr_max
USE nucbrn_module, ONLY : xn_n=>xn, be_nuc_rep_n=>be_nuc_rep, a_nuc_rep_n=>a_nuc_rep, &
& z_nuc_rep_n=>z_nuc_rep, fescrn, fascrn, uburn_n=>uburn, nse_n=>nse
USE nu_dist_module, ONLY : psi0, psi1, vol, e_rad, elec_rad, nnucr, nnurad
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY: dtime_hydro, dtnph, time

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                             :: imin          ! inner physical x-zone center
INTEGER, INTENT(in)                             :: imax          ! outer physical x-zone center

INTEGER, INTENT(in)                             :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)                             :: ik_ray        ! k-index of a radial ray
INTEGER, INTENT(in)                             :: ij_ray_dim    ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                             :: ik_ray_dim    ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                             :: nx            ! x-array extent
INTEGER, INTENT(in)                             :: nez           ! neutrino energy array extent
INTEGER, INTENT(in)                             :: nnu           ! neutrino flavor array extent
INTEGER, INTENT(in)                             :: nnc           ! composition array extent

INTEGER, INTENT(in)                             :: cycle_number  ! cycle number

INTEGER, INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim) :: nsep ! nuclear sttistical equilibrium flag

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rhorp        ! density (cm^{-3}) before hydro step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: rhop         ! current density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: trp          ! temperature (MeV) before hydro step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: tp           ! current temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: yerp         ! electron fraction before hydro step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: yep          ! current electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx)                               :: rhobarp      ! mean density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                               :: rp           ! radial zone radii (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: up           ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: gtot_pot_c   ! unshifted zone-centered gravitational potential energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0p        ! zeroth angular moment of the NDS
REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1p        ! first angular moment of the NDS

REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim)     :: xnp          ! mass fractions
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: be_nuc_repp  ! nuclear binding energy
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: a_nuc_repp   ! nuclear mass number
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)         :: z_nuc_repp   ! nuclear charge number

REAL(KIND=double), INTENT(in)                   :: dtime         ! time step used at current cycle
REAL(KIND=double), INTENT(in)                   :: time_elapsed  ! elapsed time

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)                  :: u_ge_tot      ! total NT gravitational potential energy
REAL(KIND=double), INTENT(out)                  :: u_ie_tot      ! total internal energy
REAL(KIND=double), INTENT(out)                  :: u_ke_tot      ! total NT kinetic energy
REAL(KIND=double), INTENT(out)                  :: u_ne_tot      ! total NT neutrino energy
REAL(KIND=double), INTENT(out)                  :: u_ke_x_tot    ! total NT x-kinetic energy
REAL(KIND=double), INTENT(out)                  :: u_ke_y_tot    ! total NT y-kinetic energy
REAL(KIND=double), INTENT(out)                  :: u_ke_z_tot    ! total NT z-kinetic energy
REAL(KIND=double), INTENT(out)                  :: e_radp        ! the cumulative material energy entering i_ray grid (ergs)
REAL(KIND=double), INTENT(out)                  :: radtot        ! total energy radiated by neutrinos
REAL(KIND=double), INTENT(out)                  :: u_tot         ! total energy
REAL(KIND=double), INTENT(out)                  :: elecn         ! total electron number
REAL(KIND=double), INTENT(out)                  :: elec_radp     ! the net number of electrons advected into the ij_ray,ik_ray grid
REAL(KIND=double), INTENT(out)                  :: nnucr_enu     ! the net number of e-type neutrinos currently residing in ij_ray,ik_ray grid
REAL(KIND=double), INTENT(out)                  :: nnucr_enubar  ! the net number of e-type antineutrinos currently residing in ij_ray,ik_ray grid
REAL(KIND=double), INTENT(out)                  :: nnurad_enu    ! the cumulative number of e-type neutrinos emitted by the ij_ray,ik_ray grid
REAL(KIND=double), INTENT(out)                  :: nnurad_enubar ! the cumulative number of e-type antineutrinos emitted by the ij_ray,ik_ray grid
REAL(KIND=double), INTENT(out)                  :: totlpn        ! total number of leptons emitted or residing in the core

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                         :: jr_maxp       ! jr_max + 1
INTEGER                                         :: j             ! radial zone index
INTEGER                                         :: n             ! neutrino flavor index

REAL(KIND=double), DIMENSION(nx)                :: rhobar        ! mean density (MGFLD indexed( (cm^{-3})

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

jr_min                    = imin + 1
jr_maxp                   = jr_max + 1
jr_max                    = imax + 1
ncycle                    = cycle_number

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to mgfld arrays
!-----------------------------------------------------------------------

rhor  (jr_min:jr_max)     = rhorp     (imin:imax,ij_ray,ik_ray)
rho   (jr_min:jr_max)     = rhop      (imin:imax,ij_ray,ik_ray)
tr    (jr_min:jr_max)     = trp       (imin:imax,ij_ray,ik_ray)
t     (jr_min:jr_max)     = tp        (imin:imax,ij_ray,ik_ray)
yer   (jr_min:jr_max)     = yerp      (imin:imax,ij_ray,ik_ray)
ye    (jr_min:jr_max)     = yep       (imin:imax,ij_ray,ik_ray)
u     (jr_min:jr_max)     = up        (imin:imax,ij_ray,ik_ray)
g_pot (jr_min:jr_max)     = gtot_pot_c(imin:imax,ij_ray,ik_ray)
rhobar(jr_min:jr_max)     = rhobarp(imin:imax)

psi0(jr_min:jr_max,:,:)   = psi0p(imin:imax,  :,:,ij_ray,ik_ray)
psi1(imin:imax+1,  :,:)   = psi1p(imin:imax+1,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer zone-edgeed independent variables to mgfld arrays
!-----------------------------------------------------------------------

r(imin:imax+1)            = rp(imin:imax+1)

!-----------------------------------------------------------------------
!  Transfer zone-edgeed composition variables
!-----------------------------------------------------------------------

nse(jr_min:jr_max,ij_ray,ik_ray) = nsep(imin:imax,ij_ray,ik_ray)
be_nuc_rep(jr_min:jr_max) = be_nuc_repp(imin:imax,ij_ray,ik_ray)
a_nuc_rep (jr_min:jr_max) = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep (jr_min:jr_max) = z_nuc_repp (imin:imax,ij_ray,ik_ray)

xn(jr_min:jr_max,:)       = xnp(imin:imax,:,ij_ray,ik_ray)

nse_n       (jr_min:jr_max) = nsep       (imin:imax,ij_ray,ik_ray)
be_nuc_rep_n(jr_min:jr_max) = be_nuc_repp(imin:imax,ij_ray,ik_ray)
a_nuc_rep_n (jr_min:jr_max) = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep_n (jr_min:jr_max) = z_nuc_repp (imin:imax,ij_ray,ik_ray)

xn_n(jr_min:jr_max,:)     = xnp(imin:imax,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Derived zone-centered and zone-edged variables
!-----------------------------------------------------------------------

!........radial zone thickness

dr(jr_min:jr_max)         = r(jr_min:jr_max) - r(jr_min-1:jr_max-1)

!........Newtonian rest masses

IF ( r(1) == zero ) THEN
  dmrst(1)                = zero
ELSE
  dmrst(1)                = frpith * r(1)**3 * rho(2)
END IF

rstmss(1)                 = dmrst(1)

vol   (jr_min:jr_max)     = frpith * ( r(jr_min:jr_max) * ( r(jr_min:jr_max) + r(jr_min-1:jr_max-1) ) &
&                         + r(jr_min-1:jr_max-1) * r(jr_min-1:jr_max-1) ) * dr(jr_min:jr_max)
dmrst (jr_min:jr_max)     = vol(jr_min:jr_max) * rho(jr_min:jr_max)
DO j = jr_min,jr_max
  rstmss(j)               = rstmss(j-1) + dmrst(j)
END DO

!-----------------------------------------------------------------------
!  Set timestep
!-----------------------------------------------------------------------

dtnph                     = dtime
time                      = time_elapsed

!-----------------------------------------------------------------------
!  Update EOS variables
!-----------------------------------------------------------------------

CALL eqstz_x( jr_min, jr_maxp, rho, t, ye, ij_ray, ik_ray )
CALL gammaz_x( jr_min, jr_max, rho, t, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Compute the neutrino fluxes
!-----------------------------------------------------------------------

DO n = 1,nnu
  CALL flux( jr_min, jr_max, n )
END DO

!-----------------------------------------------------------------------
!
!          \\\\\ MONITOR ENERGY AND LEPTON CONSERVATION /////
!
!-----------------------------------------------------------------------

CALL cnfg_energy( ij_ray, ik_ray, u_ge_tot, u_ie_tot, u_ke_tot, u_ne_tot, &
& u_ke_x_tot, u_ke_y_tot, u_ke_z_tot, radtot, elecn, totlpn, nx, nnu )

u_tot                     = u_ge_tot + u_ie_tot + u_ke_tot + u_ne_tot &
&                         + e_rad(ij_ray,ik_ray) + radtot
e_radp                    = e_rad   (ij_ray,ik_ray)
elec_radp                 = elec_rad(ij_ray,ik_ray)
nnucr_enu                 = nnucr (1,ij_ray,ik_ray)
nnucr_enubar              = nnucr (2,ij_ray,ik_ray)
nnurad_enu                = nnurad(1,ij_ray,ik_ray)
nnurad_enubar             = nnurad(2,ij_ray,ik_ray)

RETURN
END SUBROUTINE monitor
