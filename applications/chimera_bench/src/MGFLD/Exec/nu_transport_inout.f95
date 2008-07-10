SUBROUTINE nu_transport_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim,   &
& ik_ray_dim, nez, nnu, nnc, nprintp, rhop, tp, yep, rhobarp, rp, up, agr_e, &
& agr_c, psi0p, psi1p, rhs1_c, dc_e, dt, jdt, dtnph_trans, dtime_trans, xn_c )
!-----------------------------------------------------------------------
!
!    File:         nu_transport_inout
!    Module:       nu_transport_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/02/04
!
!    Purpose:
!        To load variables for the mgfld transport step.
!
!    Input arguments:
!
!  imin               : unshifted inner x-array index
!  imax               : unshifted outer x-array index
!  nx                 : x-array extent
!  ij_ray             : j-index of a radial ray
!  ik_ray             : k-index of a radial ray
!  ij_ray_dim         : number of y-zones on a processor before swapping
!  ik_ray_dim         : number of z-zones on a processor before swapping
!  nez                : neutrino energy array extent
!  nnu                : neutrino flavor array extent
!  nnc                : composition array extent
!  nprintp            : unit number to print diagnostics
!  rhop               : density after hydro advance (cm^{-3})
!  tp                 : initial temperatures (K)
!  yep                : initial electron fractions
!  rhobarp            : mean density (cm^{-3})
!  rp                 : x grid zone faces (cm)
!  up                 : velocity x-components (cm s^{-1})
!  agr_e              : unshifted zone-edged lapse function
!  agr_c              : unshifted zone-centered lapse function
!  psi0p              : initial zero angular moments of the neutrino occupation number
!  psi1p              : initial first angular moments of the neutrino occupation number
!  dtnph_trans        : source and transport time step
!  xn_c               : composition mass fractions
!
!    Output arguments:
!
!  tp                 : updated temperatures (K)
!  yep                : updated electron fractions
!  psi0p              : updated zero angular moments of the neutrino occupation number
!  psi1p              : updated first angular moments of the neutrino occupation number
!  rhs1_c             : the right-hand side of transport equation, first moment
!  dc_e               : the diffusion coefficient for neutrinos
!  dt                 : minimum transport time step restrictions for criteria i, radial ray ij_ray, ik_ray (s)
!  jdt                : radial zone causing minimum time step for criteria i, radial ray ij_ray, ik_ray
!  dtime_trans        : new source and transport time step given by radial ray ij_ray, ik_ray
!
!    Input-Output arguments:
!
!
!    Subprograms called:
!
!  eqstt_x            : compues the entropy before the transport step
!  mgfld_reset        : updates neutrino interaction rates
!  mgfld_transport    : performs the neutrino transport step
!  tgvndsye_x         : computes the temperature given rho, s, and ye
!  time_step_nu_trans : calculates the new transport time step
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : msolar

USE mdl_cnfg_module, ONLY : rho, ye, t, dr, r, u, dmrst, rstmss, agr, agrh, &
& jr_min, jr_max
USE nu_dist_module, ONLY : psi0, psi1, stress_x, dc, rhs1, xn
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY : dtnph_trans_t=>dtnph_trans

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin             ! inner x-array index
INTEGER, INTENT(in)              :: imax             ! outer x-array index
INTEGER, INTENT(in)              :: nx               ! x-array extent

INTEGER, INTENT(in)              :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray           ! k-index of a radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent
INTEGER, INTENT(in)              :: nnc              ! composition array extent

INTEGER, INTENT(in)              :: nprintp          ! unit number to print diagnostics

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: rhop            ! density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                         :: rhobarp         ! mean density (cm^{-3})
REAL(KIND=double), INTENT(in), DIMENSION(nx)                         :: rp              ! x-ccordinates (cm)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: up              ! zone-centered x-velocity velocity (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim) :: agr_e           ! unshifted zone-edged lapse function
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: agr_c           ! unshifted zone-entered lapse function
REAL(KIND=double), INTENT(in)                            :: dtnph_trans     ! source and transport time step

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)            :: jdt             ! zone causing dt

REAL(KIND=double)                                                    :: dtime_trans     ! new source and transport time step
REAL(KIND=double), INTENT(out), DIMENSION(50,ij_ray_dim,ik_ray_dim)  :: dt              ! minimum allowed time step
REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim)   :: rhs1_c ! right-hand side of transport equation, first moment
REAL(KIND=double), INTENT(out), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim)   :: dc_e   ! diffusion coefficient for neutrinos

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)          :: tp    ! temperature (K)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)          :: yep   ! electron fraction
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim)  :: psi0p ! zero moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim)  :: psi1p ! first moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim)      :: xn_c  ! composition mass fractions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j               ! shifted radial zone index
INTEGER                          :: ivar = 3        ! EOS entropy index

REAL(KIND=double), DIMENSION(nx) :: s               ! entropy
REAL(KIND=double)                :: dsdd            ! d(entropy)/d(rho)
REAL(KIND=double)                :: dsdt            ! d(entropy)/d(t)
REAL(KIND=double)                :: dsdy            ! d(entropy)/d(ye)
REAL(KIND=double)                :: t_guess         ! temperature guess for conserving entropy
REAL(KIND=double), DIMENSION(nx) :: rhobar          ! mean density (MGFLD indexed( (cm^{-3})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Add shift to minimum and maximum x-array extent
!-----------------------------------------------------------------------

jr_min                        = imin + 1
jr_max                        = imax + 1

!-----------------------------------------------------------------------
!
!               \\\\\ SET UP FOR NEURINO TRANSPORT ///
!
!        Load variables received from radial_ray_module into
!         mdl_cnfg_module and nu_dist_module. MGFLD varialbles
!         are alligned from jr_min = 2 to jr_max.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer arrays to shifted arrays for MGFLD
!-----------------------------------------------------------------------

rho   (jr_min:jr_max)         = rhop (imin:imax,ij_ray,ik_ray)
t     (jr_min:jr_max)         = tp   (imin:imax,ij_ray,ik_ray)
ye    (jr_min:jr_max)         = yep  (imin:imax,ij_ray,ik_ray)
u     (jr_min:jr_max)         = up   (imin:imax,ij_ray,ik_ray)
agrh  (jr_min:jr_max)         = agr_c(imin:imax,ij_ray,ik_ray)

rhobar(jr_min:jr_max)         = rhobarp(imin:imax)

r  (imin:imax+1)              = rp   (imin:imax+1)
agr(imin:imax+1)              = agr_e(imin:imax+1,ij_ray,ik_ray)

psi0(jr_min:jr_max,:,:)       = psi0p(imin:imax,  :,:,ij_ray,ik_ray)
psi1(imin:imax+1  ,:,:)       = psi1p(imin:imax+1,:,:,ij_ray,ik_ray)

xn(jr_min:jr_max,:)           = xn_c(imin:imax,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Derived zone-centered and zone-edged variables
!-----------------------------------------------------------------------

!........radial zone thickness

dr(jr_min:jr_max)             = r(jr_min:jr_max) - r(jr_min-1:jr_max-1)

!........radial ghost zone

r(jr_max+1)                   = r(jr_max) + dr(jr_max)

!........Newtonian rest masses

IF ( r(1) == zero ) THEN
  dmrst(1)                    = zero
ELSE
  dmrst(1)                    = frpith * r(1)**3 * rho(2)
END IF

rstmss(1)                     = dmrst(1)

DO j = jr_min,jr_max
  dmrst(j)                    = frpith * ( r(j) * ( r(j) + r(j-1) ) + r(j-1) * r(j-1) ) &
&                             * dr(j) * rho(j)
  rstmss(j)                   = rstmss(j-1) + dmrst(j)
END DO
dmrst(jr_max+1)               = dmrst(jr_max)

!........neutrino source and transport time step

dtnph_trans_t                 = dtnph_trans

!-----------------------------------------------------------------------
!  Get entropies before the transport step
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  CALL eqstt_x( ivar, j, ij_ray, ik_ray, rho(j), t(j), ye(j), s(j), dsdd, &
&  dsdt, dsdy )
END DO

!-----------------------------------------------------------------------
!
!             \\\\\ NEUTRINO TRANSPORT STEP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Reset the neutrino rate tables
!-----------------------------------------------------------------------

CALL mgfld_reset( ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Execute the neutrino transport
!-----------------------------------------------------------------------

CALL mgfld_transport( ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, nx, nez, &
& nnu, jdt, dtnph_trans )

!-----------------------------------------------------------------------
!  Adjust temperature
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  IF ( rho(j) > 1.d+16 ) THEN
    t_guess                   = t(j)
    CALL tgvndsye_x( j, ij_ray, ik_ray, rho(j), s(j), ye(j), t_guess, t(j) )
  END IF
END DO

!-----------------------------------------------------------------------
!
!                \\\\\  TRANSPORT TIME STEP /////
!
!-----------------------------------------------------------------------

CALL time_step_nu_trans( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, dt, jdt, dtime_trans, nnu )

!-----------------------------------------------------------------------
!
!             \\\\\ RETURN UPDAATED VARIABLES /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer shifted arrays to unshifted arrays for return
!-----------------------------------------------------------------------

tp    (imin:imax,ij_ray,ik_ray)       = t (jr_min:jr_max)
yep   (imin:imax,ij_ray,ik_ray)       = ye(jr_min:jr_max)

psi0p (imin:imax,  :,:,ij_ray,ik_ray) = psi0(jr_min:jr_max,:,:)
rhs1_c(imin:imax,  :,:,ij_ray,ik_ray) = rhs1(jr_min:jr_max,:,:)

psi1p (imin:imax+1,:,:,ij_ray,ik_ray) = psi1(imin:imax+1,  :,:)
dc_e  (imin:imax+1,:,:,ij_ray,ik_ray) = dc  (imin:imax+1,  :,:)

xn_c  (imin:imax,:,ij_ray,ik_ray)     = xn(jr_min:jr_max,:)

RETURN
END SUBROUTINE nu_transport_inout
