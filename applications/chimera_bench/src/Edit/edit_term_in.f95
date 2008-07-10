SUBROUTINE edit_term_in( imin, imax, idim, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nx, nez, nnu, nnc, rhorp, rhop, trp, tp, yerp, yep, rp,  &
& up, psi0p, psi1p, dtime, time_elapsed, cycle_number, xnp, be_nuc_repp, &
& a_nuc_repp, z_nuc_repp, nsep )
!-----------------------------------------------------------------------
!
!    File:         edit_term_in
!    Module:       edit_term_in
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
!  esrgnz_comp_x : updates EOS on cube corners if nse    = 0 (to incorporate latest xn)
!  eqstz_x       : updates EOS variables
!  gammaz_x      : updates EOS gammas
!  edit_term     : performs the final edit
!
!    Input arguments:
!
!  imin          : inner physical x-zone center
!  imax          : outer physical x-zone center
!  idim          : logical dimension of x-zone
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
!  psi0p         : zeroth angular moment of the NDS
!  dtime_h       : hydro time step
!  dtime         : time step used at current cycle
!  time_elapsed  : elapsed time
!  cycle_number  : cycle number
!  xnp           : initial mass fractions
!  be_nuc_repp   : binding energy of mean heavy nucleus
!  a_nuc_repp    : mass number of mean heavy nucleus
!  z_nuc_repp    : charge number of mean heavy nucleus
!  nesp          : nuclear statistical equilibrium flag
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  cycle_module, eos_snc_x_module, mdl_cnfg_module, nucbrn_module,
!  nu_dist_module, nu_energy_grid_module, t_cntrl_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, frpith
USE physcnst_module, ONLY : kmev_inv

USE cycle_module, ONLY : ncycle
USE eos_snc_x_module, ONLY : xn, be_nuc_rep, a_nuc_rep, z_nuc_rep, nse
USE mdl_cnfg_module, ONLY : rhor, rho, tr, t, yer, ye, dr, r, u, dmrst, rstmss, &
& jr_min_m   =>jr_min, jr_max_m   =>jr_max
USE nucbrn_module, ONLY : xn_n   =>xn, be_nuc_rep_n   =>be_nuc_rep, a_nuc_rep_n   =>a_nuc_rep, &
&  z_nuc_rep_n   =>z_nuc_rep, fescrn, fascrn, uburn_n   =>uburn, nse_n   =>nse
USE nu_dist_module, ONLY : psi0, psi1
USE nu_energy_grid_module, ONLY : nnugp
USE t_cntrl_module, ONLY: dtime_hydro, dtnph, time

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                   :: imin             ! inner physical x-zone center
INTEGER, INTENT(in)                   :: imax             ! outer physical x-zone center
INTEGER, INTENT(in)                   :: idim             ! logical dimension of x-zone

INTEGER, INTENT(in)                   :: ij_ray           ! j-index of a radial ray
INTEGER, INTENT(in)                   :: ik_ray           ! k-index of a radial ray
INTEGER, INTENT(in)                   :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)                   :: ik_ray_dim       ! number of z-zones on a processor before swapping
INTEGER, INTENT(in)                   :: nx               ! x-array extent
INTEGER, INTENT(in)                   :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)                   :: nnu              ! neutrino flavor array extent
INTEGER, INTENT(in)                   :: nnc              ! composition array extent

INTEGER, INTENT(in)                   :: cycle_number     ! cycle number

INTEGER, INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)  :: nsep   ! nuclear sttistical equilibrium flag

REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: rhorp        ! density (cm^{-3}) before hydro step
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: rhop         ! density (cm^{-3})
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: trp          ! temperature (MeV) before hydro step
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: tp           ! temperature (MeV)
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: yerp         ! electron fraction before hydro step
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: yep          ! electron fraction
REAL(KIND   =double), INTENT(in), DIMENSION(idim)                               :: rp           ! radial zone radii (cm)
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: up           ! radial velocity of zone (cm)
REAL(KIND   =double), INTENT(in), DIMENSION(idim,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0p        ! zeroth angular moment of the NDS
REAL(KIND   =double), INTENT(in), DIMENSION(idim,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1p        ! first angular moment of the NDS

REAL(KIND   =double), INTENT(in), DIMENSION(idim,nnc,ij_ray_dim,ik_ray_dim)     :: xnp          ! mass fractions
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: be_nuc_repp  ! density (cm^{-3})
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: a_nuc_repp   ! temperature (MeV)
REAL(KIND   =double), INTENT(in), DIMENSION(idim,ij_ray_dim,ik_ray_dim)         :: z_nuc_repp   ! electron fraction

REAL(KIND   =double), INTENT(in)         :: dtime         ! time step used at current cycle
REAL(KIND   =double), INTENT(in)         :: time_elapsed  ! elapsed time

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: jr_min          ! minimum radial inddex
INTEGER                          :: jr_max          ! maximum radial inddex
INTEGER                          :: jr_maxp         ! jr_max + 1
INTEGER                          :: i               ! do index
INTEGER                          :: j               ! radial zone index
INTEGER                          :: k               ! neutrino energy index
INTEGER                          :: n               ! neutrino flavor index
INTEGER                          :: nc              ! compoisition index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Transfer integer scalar parameters
!-----------------------------------------------------------------------

jr_min                    = imin + 1
jr_min_m                  = jr_min
jr_max                    = imax + 1
jr_max_m                  = jr_max
jr_maxp                   = jr_max + 1
ncycle                    = cycle_number

!-----------------------------------------------------------------------
!  Transfer zone-centered independent variables to mgfld arrays
!-----------------------------------------------------------------------

DO i    = imin,imax
  j                       = i-imin+2
  rhor(j)                 = rhorp(i,ij_ray,ik_ray)
  rho (j)                 = rhop (i,ij_ray,ik_ray)
  tr  (j)                 = trp  (i,ij_ray,ik_ray)
  t   (j)                 = tp   (i,ij_ray,ik_ray)
  yer (j)                 = yerp (i,ij_ray,ik_ray)
  ye  (j)                 = yep  (i,ij_ray,ik_ray)
  u   (j)                 = up   (i,ij_ray,ik_ray)
END DO

DO n    = 1,nnu
  IF ( nnugp(n) == 0 ) CYCLE
  DO k    = 1,nnugp(n)
    DO i    = imin,imax
      psi0(i-imin+2,k,n)  = psi0p(i,k,n,ij_ray,ik_ray)
      psi1(i-imin+1,k,n)  = psi1p(i,k,n,ij_ray,ik_ray)
    END DO
    psi1(imax-imin+2,k,n) = psi1p(imax+1,k,n,ij_ray,ik_ray)
  END DO
END DO

!-----------------------------------------------------------------------
!  Transfer zone-edgeed independent variables to mgfld arrays
!-----------------------------------------------------------------------

DO i    = imin,imax+1
  r(i-imin+1)             = rp(i)
END DO

!-----------------------------------------------------------------------
!  Transfer zone-edgeed composition variables
!-----------------------------------------------------------------------

DO i    = imin,imax
  j                       = i-imin+2
  nse(j,ij_ray,ik_ray)    = nsep   (i,ij_ray,ik_ray)
  be_nuc_rep(j)           = be_nuc_repp(i,ij_ray,ik_ray)
  a_nuc_rep (j)           = a_nuc_repp (i,ij_ray,ik_ray)
  z_nuc_rep (j)           = z_nuc_repp (i,ij_ray,ik_ray)
END DO

DO nc    = 1,nnc
  DO i    = imin,imax+1
    xn(i-imin+2,nc)         = xnp(i,nc,ij_ray,ik_ray)
  END DO
END DO

DO i    = imin,imax
  j                       = i-imin+2
  nse_n   (j)             = nsep   (i,ij_ray,ik_ray)
  be_nuc_rep_n(j)         = be_nuc_repp(i,ij_ray,ik_ray)
  a_nuc_rep_n (j)         = a_nuc_repp(i,ij_ray,ik_ray)
  z_nuc_rep_n (j)         = z_nuc_repp(i,ij_ray,ik_ray)
END DO

DO nc    = 1,nnc
  DO i    = imin,imax
    xn_n(i-imin+2,nc)     = xnp(i,nc,ij_ray,ik_ray)
  END DO
END DO

!-----------------------------------------------------------------------
!  Derived zone-centered and zone-edged variables
!-----------------------------------------------------------------------

!........radial zone thickness

DO j    = jr_min,jr_max
  dr(j)                   = r(j) - r(j-1)
END DO

!........Newtonian rest masses

IF ( r(1) == zero ) THEN
  dmrst(1)                = zero
ELSE
  dmrst(1)                = frpith * r(1)**3 * rho(2)
END IF

rstmss(1)                 = dmrst(1)

DO j    = jr_min,jr_max
  dmrst(j)                = frpith * ( r(j) * ( r(j) + r(j-1) ) + r(j-1) * r(j-1) ) * dr(j) * rho(j)
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

CALL esrgnz_comp_x( jr_min, jr_maxp, rho, t, ye, ij_ray, ik_ray )
CALL eqstz_x( jr_min, jr_maxp, rho, t, ye, ij_ray, ik_ray )
CALL gammaz_x( jr_min, jr_max, rho, t, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Execute the Edit
!-----------------------------------------------------------------------

CALL edit_term( jr_min, jr_max, ij_ray, ik_ray, ij_ray_dim, ik_ray_dim, &
& nx, nnu )

RETURN
END SUBROUTINE edit_term_in
