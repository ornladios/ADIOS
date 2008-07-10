SUBROUTINE setup_initial_trans( imin, imax, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nx, nez, nnu, rhop, tp, yep, rp, drp, up, agr_e, agr_c, &
& psi0p, psi1p, rhobar )
!-----------------------------------------------------------------------
!
!    File:         setup_initial_trans
!    Module:       setup_initial_trans
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/04/03
!
!    Purpose:
!        To compute, from the initial problem configuration, the neutrino
!         absorption, emission, scattering, pair production, and bremsstrahlung
!         pair production function tables
!
!    Subprograms called:
!  abemrate      : computes absorption and emission rates
!  abemset       : computes absorption and emission kernels on cube corners
!  agr_nu_cal    : computes neutrino lapse functions
!  bremrate      : computes bremsstrahlung rates
!  bremset       : computes bremsstrahlung kernels on cube corners
!  diffc         : computes neutrino diffusion coefficients
!  eddington     : computes neutrino Eddington factors
!  enu_cal       : computes functions of the neutrino energy zones
!  gamgr_nu_cal  : computes neutrino relativistic parameters
!  gennur        : initializes Haxtons rates
!  mfp_cal       : computes neutrino inverse mean free paths
!  nu_number     : computes neutrino numbers and energies
!  nu_sphere     : computes neutrinospheres
!  nu_stress_x   : computes neutrino stresses at radial zone interfaces
!  nu_U          : computes neutrino pressures and energy densities
!  pairrate      : computes electron-positron pair annihilation rates
!  pairset       : computes electron-positron pair annihilation kernels on cube corners
!  pblmst2       : option to change problem after neutrino quantities are calculated
!  pre_trans     : computes quantities needed for neutrino transport
!  scataset      : computes neutrino-nucleus inelastic scattering kernels on cube corners
!  scateset      : computes neutrino-electron elastic scattering kernels on cube corners
!  scatiset      : computes neutrino-nucleon isoenergetic scattering kernels on cube corners
!  scatnnset     : computes neutrino-nucleon inelastic scattering kernels on cube corners
!  scatnAset     : computes neutrino-nucleus inelastic scattering kernels on cube corners
!  scatnset      : computes neutrino-nucleon elastic scattering kernels on cube corners
!  sctarate      : computes neutrino-nucleus inelastic scattering rates
!  scterate      : computes neutrino-electron elastic scattering rates
!  sctirate      : computes neutrino-nucleon isoenergetic scattering rates
!  sctnnrate     : computes neutrino-nucleon inelastic scattering rates
!  sctnrate      : computes neutrino-nucleon elastic scattering rates
!  sctnArate     : computes neutrino-electron inelastic scattering rates
!    
!    Input arguments:
!  imin          : inner x-zone array index
!  imax          : outer x-zone array index
!  ij_ray        : j-index of a radial ray
!  ik_ray        : k-index of a radial ray
!  ij_ray_dim    : number of y-zones on a processor before swapping
!  ik_ray_dim    : number of z-zones on a processor before swapping
!  nx            : x_array extent
!  nez           : neutrino energy array extent
!  nnu           : neutrino flavor array extent
!  nnc           : neutrino abundance array extent
!  rhop          : density (cm^{-3})
!  tp            : temperature (MeV)
!  yep           : electron fraction
!  rp            : radial zone radii (cm)
!  up            : radial velocity of zone (cm)
!  agr_e         : unsifted zone-edged lapse function
!  agr_c         : unsifted zone-centered lapse function
!  psi0p         : zero angular moments of the neutrino occupation number
!  rhobar        : angularly averaged density [g cm^{-3}]
!
!    Output arguments:
!  psi1p         : first angular moments of the neutrino occupation number
!
!    Include files:
!  kind_module, numerical_module
!  cycle_module, mdl_cnfg_module, nu_dist_module, edit_module,
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module, prb_cntl_module
!     
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith

USE cycle_module, ONLY : ncynu_trns, nutrans_trns
USE edit_module, ONLY : nlog
USE mdl_cnfg_module, ONLY : jr_min, jr_max, r, rstmss, rho, t, ye, u, dr, &
& dmrst, agr, agrh, agra, agrah
USE nu_dist_module, ONLY : j_sphere, r_sphere, d_sphere, t_sphere, m_sphere, &
& psi0, psi1, rhs1, dc
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx
USE prb_cntl_module, ONLY : jnumin, jnumax, iaefnp, iaencnu, nncs, nes, &
& in, ip, ihe, iheavy, ipair, ibrem, isctn, isctnn, isctnA

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin            ! inner x_array index
INTEGER, INTENT(in)              :: imax            ! outer x_array index

INTEGER, INTENT(in)              :: ij_ray          ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray          ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim      ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim      ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)              :: nx              ! x-array extent
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: tp          ! temperature (MeV)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: yep         ! electron fraction
REAL(KIND=double), INTENT(in), DIMENSION(nx+1)                         :: rp          ! radial zone radii (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                           :: drp         ! radial zone thickness (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: up          ! radial velocity of zone (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                           :: rhobar      ! angularly averaged density [g cm^{-3}]

REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)   :: agr_e       ! unsifted zone-edged lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: agr_c       ! unsifted zone-centered lapse function

REAL(KIND=double), INTENT(in), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0p   ! zero moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1p   ! first moment of the neutrino occupation probability

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: rhop        ! density (cm^{-3})

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j               ! radial zone index
INTEGER                          :: jr_maxp         ! jr_max + 1
INTEGER                          :: jnumaxp         ! maximum neutrino radial zone index
INTEGER                          :: n               ! neutrino flavor index

REAL(KIND=double), DIMENSION(nx) :: nu_strs_c       ! zone-centered x-component of neutrino stress
REAL(KIND=double), DIMENSION(nx) :: nu_strs_e       ! zone-edged x-component of neutrino stress

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' Variable initialization complete in setup_initial_trans')
  103 FORMAT (' Variable transfer complete in setup_initial_trans')
  111 FORMAT (' Neutrino gammas set in setup_initial_trans')
  113 FORMAT (' Neutrino lapse functions set in setup_initial_trans')
  115 FORMAT (' Relativistic neutrino group energies computed in setup_initial_trans')
  117 FORMAT (' Transport quantities computed in setup_initial_trans')
  119 FORMAT (' pblmst2 called in setup_initial_trans')
  121 FORMAT (' Haxton rates cmputed in setup_initial_trans')
  123 FORMAT (' Tabular ec rates read in setup_initial_trans')

  131 FORMAT (' Emission and absorption rates are being set in setup_initial_trans')
  133 FORMAT (' Emission and absorption rates have been deselected by options')
  135 FORMAT (' Emission and absorption rates have been set in setup_initial_trans')
  141 FORMAT (' Haxton rates are being set in setup_initial_trans')
  143 FORMAT (' Haxton rates have been deselected by options')
  145 FORMAT (' Haxton rates have been set in setup_initial_trans')
  151 FORMAT (' Neutrino-electron scattering rates are being set in setup_initial_trans')
  153 FORMAT (' Neutrino-electron scattering rates have been deselected by options')
  155 FORMAT (' Neutrino-electron scattering rates have been set in setup_initial_trans')
  161 FORMAT (' Neutrino-nucleon and nucleus isoenergetic scattering rates are being set in setup_initial_trans')
  163 FORMAT (' Neutrino-nucleon and nucleus isoenergetic scattering rates have been deselected by options')
  165 FORMAT (' Neutrino-nucleon and nucleus isoenergetic rates have been set in setup_initial_trans')
  171 FORMAT (' Electron-positron pair annihilation rates are being set in setup_initial_trans')
  173 FORMAT (' Electron-positron pair annihilation rates have been deselected by options')
  175 FORMAT (' Electron-positron pair annihilation rates have been set in setup_initial_trans')
  181 FORMAT (' Nucleon-nucleon bremsstrahlung rates are being set in setup_initial_trans')
  183 FORMAT (' Nucleon-nucleon bremsstrahlung rates have been deselected by options')
  185 FORMAT (' Nucleon-nucleon bremsstrahlung rates have been set in setup_initial_trans')
  191 FORMAT (' Neutrino-nucleon elastic scattering rates are being set in setup_initial_trans')
  193 FORMAT (' Neutrino-nucleon elastic scattering rates have been deselected by options')
  195 FORMAT (' Neutrino-nucleon elastic rates have been set in setup_initial_trans')
  201 FORMAT (' Neutrino-nucleon inelastic scattering rates are being set in setup_initial_trans')
  203 FORMAT (' Neutrino-nucleon inelastic scattering rates have been deselected by options')
  205 FORMAT (' Neutrino-nucleon inelastic rates have been set in setup_initial_trans')
  211 FORMAT (' Neutrino-nucleus inelastic scattering rates are being set in setup_initial_trans')
  213 FORMAT (' Neutrino-nucleus inelastic scattering rates have been deselected by options')
  215 FORMAT (' Neutrino-nucleus inelastic rates have been set in setup_initial_trans')

  301 FORMAT (' Emission and absorption rates interpolated in setup_initial_trans')
  303 FORMAT (' Isoenergetic scattering rates interpolated in setup_initial_trans')
  305 FORMAT (' Haxtons neutrino-nucleus inelastic scattering rates interpolated in setup_initial_trans')
  307 FORMAT (' Neutrino-electron elastic scattering rates interpolated in setup_initial_trans')
  309 FORMAT (' Electron-positron pair annihilation rates interpolated in setup_initial_trans')
  311 FORMAT (' Nucleon-nucleon bremsstrahlung rates interpolated in setup_initial_trans')
  313 FORMAT (' Neutrino-nucleon elastic scattering rates interpolated in setup_initial_trans')
  315 FORMAT (' Neutrino-nucleon inelastic scattering rates interpolated in setup_initial_trans')
  317 FORMAT (' Neutrino-nucleus inelastic scattering rates interpolated in setup_initial_trans')

  401 FORMAT (' Neutrino mean free paths computed in setup_initial_trans')
  403 FORMAT (' Neutrinospheres computed in setup_initial_trans')
  405 FORMAT (' Neutrinos diffusion coefficients computed in setup_initial_trans')
  407 FORMAT (' Neutrino-number data computed in setup_initial_trans')
  409 FORMAT (' Neutrino stresses computed in setup_initial_trans')
  411 FORMAT (' Neutrino Eddington factors computed in setup_initial_trans')
  413 FORMAT (' Neutrino energy densities computed in setup_initial_trans')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize shifted radial array index boundaries
!-----------------------------------------------------------------------

jr_min                    = imin + 1
jr_max                    = imax + 1
jr_maxp                   = jr_max + 1
IF ( jnumax > jr_max ) jnumax = jr_max
jnumaxp                   = jnumax + 1

IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,101)


!-----------------------------------------------------------------------
!
!                  \\\\\ TRANSFER VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to shifted mgfld arrays
!-----------------------------------------------------------------------

rho  (imin+1:imax+1)      = rhop (imin:imax,ij_ray,ik_ray)
t    (imin+1:imax+1)      = tp   (imin:imax,ij_ray,ik_ray)
ye   (imin+1:imax+1)      = yep  (imin:imax,ij_ray,ik_ray)
u    (imin+1:imax+1)      = up   (imin:imax,ij_ray,ik_ray)
agr  (imin+1:imax+1)      = agr_e(imin:imax,ij_ray,ik_ray)
agrh (imin+1:imax+1)      = agr_c(imin:imax,ij_ray,ik_ray)
agra (imin+1:imax+1)      = agr_e(imin:imax,ij_ray,ik_ray)
agrah(imin+1:imax+1)      = agr_c(imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Transfer zone-edgeed independent variables to shifted mgfld arrays
!-----------------------------------------------------------------------

r(imin:imax+1)            = rp(imin:imax+1)
r(imax-imin+3)            = r(imax-imin+2) + ( r(imax-imin+2) - r(imax-imin+1) )

!-----------------------------------------------------------------------
!  Derived zone-centered and zone-edged variables
!-----------------------------------------------------------------------

!........radial zone thickness

dr(imin+1:imax+1)         = drp(imin:imax)

!........Newtonian rest masses

IF ( r(1) == zero ) THEN
  dmrst(1)                = zero
ELSE
  dmrst(1)                = frpith * r(1)**3 * rho(2)
END IF

rstmss(1)                 = dmrst(1)

DO j = jr_min,jr_max
  dmrst(j)                = frpith * ( r(j) * ( r(j) + r(j-1) ) + r(j-1) * r(j-1) ) * dr(j) * rho(j)
  rstmss(j)               = rstmss(j-1) + dmrst(j)
END DO

IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,103)

!-----------------------------------------------------------------------
!
!                   \\\\\ NEUTRINO VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Load neutrino occupation distribution
!-----------------------------------------------------------------------

psi0(jr_min:jr_max,:,:)   = psi0p(imin:imax,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Neutrino relativistic variables
!-----------------------------------------------------------------------

CALL gamgr_nu_cal( jr_min, jr_max )
IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,111)

CALL agr_nu_cal( jr_min, jr_max )
IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,113)

IF ( nnugpmx > 0 )THEN

!-----------------------------------------------------------------------
!  GR neutrino energies
!-----------------------------------------------------------------------

  CALL enu_cal( jnumin, jnumaxp )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,115)

!-----------------------------------------------------------------------
!  Quantities needed for neutrino transport
!-----------------------------------------------------------------------

  CALL pre_trans( jnumin, jnumax, rho, r, nx, nnu )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,117)

!-----------------------------------------------------------------------
!  Modify problem given the neutrino energies.
!  The modifications are contained in subroutine pblmst2. If no
!   modifications, subroutine pblmst2 is a dummy subroutine.
!-----------------------------------------------------------------------

  CALL pblmst2
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,119)

!-----------------------------------------------------------------------
!  Read in and regrid Wick's neutrino interaction rates
!-----------------------------------------------------------------------

  CALL gennur
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,121)

!-----------------------------------------------------------------------
!  Read in and regrid tabular electron capture rates
!-----------------------------------------------------------------------

  CALL read_ec_table
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,123)

!-----------------------------------------------------------------------
!  Construct neutrino interaction rate tables
!-----------------------------------------------------------------------

!........Emission and absorption

  jnumaxp               = jnumax + 1
  IF ( iaefnp /= 0  .or.  iaencnu /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,131)
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,133)
  END IF
  DO j = jnumin,jnumaxp
    CALL abemset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) ) 
  END DO
  IF ( iaefnp /= 0  .or.  iaencnu /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,135)
  END IF

!........Neutrino-nucleus inelastic scattering (Haxton)

  IF ( nncs /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,141)
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,143)
  END IF
  DO j = jnumin,jnumaxp
    CALL scataset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( nncs /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,145)
  END IF

!........Neutrino-electron scattering

  IF ( nes /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,151)
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,153)
  END IF
  DO j = jnumin,jnumaxp
    CALL scateset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( nes /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,155)
  END IF

!........Neutrino-nucleon and neutrino-nucleus isoenergetic scattering

  IF ( in /= 0  .or.  ip /= 0  .or.  ihe /= 0  .or.  iheavy /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,161)
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,163)
  END IF
  DO j = jnumin,jnumaxp
    CALL scatiset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( in /= 0  .or.  ip /= 0  .or.  ihe /= 0  .or.  iheavy /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,165)
  END IF

!........Electron-positron pair annihilation

  IF ( ipair /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,171)
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,173)
  END IF
  DO j = jnumin,jnumaxp
    CALL pairset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( ipair /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,175)
  END IF

!........Nucleon-nucleon bremsstrahlung

  IF ( ibrem /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,181)
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,183)
  END IF
  DO j = jnumin,jnumaxp
    CALL bremset  ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( ibrem /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,185)
  END IF

!........Neutrino-nucleon elastic scattering

  IF ( isctn /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim )  WRITE (nlog,191)
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,193) 
  END IF
  DO j = jnumin,jnumaxp
    CALL scatnset ( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( isctn /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,195)
  END IF

!........Neutrino-nucleon inelastic scattering

  IF ( isctnn /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,201) 
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,203)
  END IF
  DO j = jnumin,jnumaxp
    CALL scatnnset( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( isctnn /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,205)
  END IF

!........Neutrino-nucleus inelastic scattering (Fuller & Meyer)

  IF ( isctnA /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,211) 
  ELSE
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,213)
  END IF
  DO j = jnumin,jnumaxp
    CALL scatnAset( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  END DO
  IF ( isctnA /= 0 ) THEN
    IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,215)
  END IF

!-----------------------------------------------------------------------
!  Interpolate neutrino interaction rates
!-----------------------------------------------------------------------

  CALL abemrate( jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,301)

  DO n = 1,nnu
    CALL sctirate( jr_min, jr_max, ij_ray, ik_ray, n, rho, t, ye, nx )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,303)

  DO n = 1,nnu
    CALL sctarate( n, jr_min, jr_max, ij_ray, ik_ray )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,305)

  DO n = 1,nnu
    CALL scterate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,307)

  DO n = 1,nnu
    CALL pairrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, r, nx)
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,309)

  DO n = 1,nnu
    CALL bremrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, r, nx )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,311)

  DO n = 1,nnu
    CALL sctnrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,313)

  DO n = 1,nnu
    CALL sctnnrate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,315)

  DO n = 1,nnu
    CALL sctnArate( n, jr_min, jr_max, ij_ray, ik_ray, rho, t, ye, nx )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,317)

!-----------------------------------------------------------------------
!  Compute functions of neutrino interaction rates
!-----------------------------------------------------------------------
 
  CALL mfp_cal( jr_min, jr_max, ij_ray, ik_ray )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,401)

  CALL nu_sphere( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim,ik_ray_dim, &
&  r, rho, t, rstmss, nx, nez, nnu, j_sphere, r_sphere, d_sphere, t_sphere, &
&  m_sphere )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,403)

  CALL psi1_cal( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
&  rho, t, ye, r, rstmss, u, psi0, psi1, nx, nez, nnu, 1 )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,405)

  DO n = 1,nnu
    CALL nu_number( jr_min, jr_max, n, ij_ray, ik_ray, nx, nez, nnu, r, u, &
&    psi0, psi1 )
  END DO
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,407)

  CALL nu_stress_x( jr_min, jr_max, ij_ray, ik_ray, nx, nez, nnu, rho, &
&  rhobar, r, rhs1, dc, psi0, nu_strs_c, nu_strs_e )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,409)

  CALL eddington( jr_min, jr_max )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,411)

  CALL nu_U( jr_min, jr_max, rho, nx, nnu )
  IF ( ij_ray == ij_ray_dim  .and.  ik_ray == ik_ray_dim ) WRITE (nlog,413)

!-----------------------------------------------------------------------
!  Initialize transport implementation keys and transfer variables back
!   to radial_ray_module
!-----------------------------------------------------------------------

  nutrans_trns            = .true.
  ncynu_trns              = 0

  psi1p(imin:imax+1,:,:,ij_ray, ik_ray) = psi1(imin:imax+1,:,:)

END IF ! nnugpmx > 0

RETURN
END SUBROUTINE setup_initial_trans
