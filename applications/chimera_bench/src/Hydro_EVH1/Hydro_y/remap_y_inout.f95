SUBROUTINE remap_y_inout( nx, ny, jmin, jmax, ji_ray, jk_ray, j_ray_dim, &
& ik_ray_dim, i_radial, nez, nnu, ls, le, ldim, y_el, dy_cl, y_cl, y_ef, &
& dy_cf, y_cf, rhop, tp, yep, eip, up, vp, wp, psi0p, xnp, a_nuc_repp,   &
& z_nuc_repp, be_nuc_repp, flat_x_y, time, t_bounce, tb_dy_shift, rhobar )
!-----------------------------------------------------------------------
!
!    File:         remap_y_inout
!    Module:       remap_y_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive the EVH1 arrays from radhyd_angular_module, execute the
!       y-array remap, and return to radhyd_angular_module the updated
!       variables
!
!    Input arguments:
!  nx             : x-array extent
!  ny             : y-array extent
!  jmin           : lower y-array index
!  jmax           : upper y-array index
!  ji_ray         : x (radial) index of a specific y (angular) ray
!  jk_ray         : z (azimuthal) index of a specific y (angular) ray
!  j_ray_dim      : the number of radial zones on a processor after swapping with y
!  ik_ray_dim     : the number of z-zones on a processor before swapping with z
!  i_radial       : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!  y_el           : y grid zone left interfaces after Lagrangian update
!  dy_cl          : y_el(j+1) - y_el(j)
!  y_cl           : y grid zone midpoints after Lagrangian update
!  y_ef           : final y grid zone left interfaces
!  dy_cf          : y_ef(j+1) - y_ef(j)
!  y_cf           : final y grid zone midpoints
!  rhop           : densities (g cm^{-3}) after Lagrangian update
!  tp             : temperatures (K) after Lagrangian update
!  yep            : electron fractions after Lagrangian update
!  eip            : internal energies (ergs g^{-1}) after Lagrangian update
!  up             : velocity x-components (cm s^{-1}) after Lagrangian update
!  vp             : velocity y-components (cm s^{-1}) after Lagrangian update
!  wp             : velocity z-components (cm s^{-1}) after Lagrangian update
!  psi0p          : zero moment of the neutrino distribution function after Lagrangian update
!  xnp            : mass fractions of nuclei not in NSE after Lagrangian update
!  a_nuc_repp     : mass number of auxiliar nucleus after Lagrangian update
!  z_nuc_repp     : charge number of auxiliar nucleus after Lagrangian update
!  be_nuc_repp    : binding energy of auxiliar nucleus after Lagrangian update
!  flat_x_y       : variables indicating the presence of radial shocks
!  time           : elapsed time
!  t_bounce       : time from bounce
!  tb_dy_shift    : time from bounce to trun off grid wiggle and/or psi0 diffusion
!  rhobar         : mean density at a given radius
!
!    Output arguments:
!  rhop           : updated densities (g cm^{-3})
!  tp             : updated temperatures (K)
!  yep            : updated electron fractions
!  eip            : updated internal energies (ergs g^{-1})
!  up             : updated velocity x-components (cm s^{-1})
!  vp             : updated velocity y-components (cm s^{-1})
!  wp             : updated velocity z-components (cm s^{-1})
!  psi0p          : updated zero moment of the neutrino distribution function
!  xnp            : updated mass fractions of nuclei not in NSE
!  a_nuc_repp     : updated mass number of auxiliar nucleus
!  z_nuc_repp     : updated charge number of auxiliar nucleus
!  be_nuc_repp    : updated binding energy of auxiliar nucleus
!
!    Subprograms called:
!  sweepbc               : inserts boundary values inro ghost zones
!  volume                : computes zone volumes
!  paraset               : computes the PPM coefficients for the new grid
!  pre_remap_psi         : sets up boundary conditions for the psi0 remap
!  remap_psi_y           : remaps the psi0's in the y-direction
!  shock_smooth_y        : marks zones lying along a shock for additional diffusion
!  pre_remap_comp        : sets up boundary conditions for the composition remap
!  remap_comp_y          : remaps the composition and the binding energies
!  remap_y               : remaps the hydro variables
!  tgvndeye_sweep_comp_y : computes the temperature given rho, e, and ye
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, edit_module eos_snc_y_module, evh1_global,
!  evh1_sweep, evh1_zone, mgfld_remap_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero

USE e_advct_module, ONLY : rhomin_y_eadvect
USE edit_module, ONLY : nlog
USE eos_snc_y_module, ONLY : xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, nse, aesv
USE evh1_global, ONLY : ngeomy, nlefty, nrighty
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, temp, u, v, w, ei, ye, xa0, &
& dx0, xa, dx, entrop
USE evh1_zone, ONLY : zparay
USE mgfld_remap_module, ONLY : psi0_re, r_r=>r, temp_r=>temp, ye_r=>ye, &
& xa_r=>xa, dx_r=>dx, xa0_r=>xa0, dx0_r=>dx0, comp
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx               ! x-array extent

INTEGER, INTENT(in)              :: jmin             ! minimum y-array index
INTEGER, INTENT(in)              :: jmax             ! maximum y-array index
INTEGER, INTENT(in)              :: ny               ! y-array extent

INTEGER, INTENT(in)              :: ji_ray           ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray           ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: j_ray_dim        ! number of radial zones on a processor after swapping with y
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of radial zones on a processor before swapping with z
INTEGER, INTENT(in)              :: i_radial         ! the unshifted radial zone corresponding to ki_ray, kj_ray

INTEGER, INTENT(in)              :: ls               ! minimum composition index
INTEGER, INTENT(in)              :: le               ! maximum composition index
INTEGER, INTENT(in)              :: ldim             ! composition array extent

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

REAL(KIND=double), INTENT(in),    DIMENSION(ny+1)                    :: y_ef        ! Eulerian x-coordinate zone edge
REAL(KIND=double), INTENT(in),    DIMENSION(ny)                      :: dy_cf       ! Eulerian x-coordinate zone thickness
REAL(KIND=double), INTENT(in),    DIMENSION(ny)                      :: y_cf        ! Eulerian x-coordinate zone center
REAL(KIND=double), INTENT(in),    DIMENSION(ny,j_ray_dim,ik_ray_dim) :: flat_x_y    ! variables indicating the presence of radial shocks
REAL(KIND=double), INTENT(in)    :: time             ! elapsed time
REAL(KIND=double), INTENT(in)    :: t_bounce         ! time from bounce
REAL(KIND=double), INTENT(in)    :: tb_dy_shift      ! time from bounce to trun off grid wiggle and/or psi0 diffusion
REAL(KIND=double), INTENT(in),    DIMENSION(nx)                      :: rhobar      ! mean density at a given radius

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(ny+1,j_ray_dim,ik_ray_dim+1)     :: y_el        ! y-coordinate zone edge after Lagrangian step
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: dy_cl       ! y-coordinate zone thickness after Lagrangian step
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: y_cl        ! y-coordinate zone center after Lagrangian step

REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: rhop        ! density (cm^{-3})
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: tp          ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: yep         ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: eip         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: up          ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: vp          ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: wp          ! z-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(ny,nez,nnu,j_ray_dim,ik_ray_dim) :: psi0p       ! zero moment of the neutrino occupation probability

REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: a_nuc_repp  ! nuclear mass number of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: z_nuc_repp  ! nuclear charge number of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(ny,j_ray_dim,ik_ray_dim)         :: be_nuc_repp ! binding energy of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(ny,ldim,j_ray_dim,ik_ray_dim)    :: xnp         ! mass fractions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                  :: l_write = .false.

INTEGER                             :: ntot           ! nmax + 6
INTEGER, DIMENSION(ny)              :: j_shock        ! zones marked for added y-diffusion

REAL(KIND=double), DIMENSION(ny)    :: rho_m          ! density (g cm^{-3})
REAL(KIND=double), DIMENSION(ny)    :: t_m            ! temperature (MeV)
REAL(KIND=double), DIMENSION(ny)    :: ye_m           ! electron fraction
REAL(KIND=double), DIMENSION(ny+12) :: rho_i          ! density after the Lagrangian step (g cm^{-3})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rho_i                       = zero

!-----------------------------------------------------------------------
!
!                   \\\\\ SETUP FOR Y-REMAP /////
!
!        Load variables received from radial_ray_module into
!         evh1_sweep and mgfld_remap_module
!
!        Subroutine remap uses the variables stored in evh1_sweep.
!        Subroutine remap_comp_y uses the variables stored in nucbrn_module
!         and mgfld_remap_module.
!        Subroutine remap_psi uses the variables stored in mgfld_remap_module
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute boundary indices of padded zones
!-----------------------------------------------------------------------

nmin                         = jmin + 6
nmax                         = jmax + 6
ntot                         = nmax + 6

!-----------------------------------------------------------------------
!  Load final coordinates
!-----------------------------------------------------------------------

xa0(nmin:nmax+1)             = y_ef(jmin:jmax+1)
dx0(nmin:nmax)               = dy_cf(jmin:jmax)
xa0_r(nmin:nmax+1)           = y_ef(jmin:jmax+1)
dx0_r(nmin:nmax)             = dy_cf(jmin:jmax)

!-----------------------------------------------------------------------
!  Load coordinates resulting from the Lagrangian step
!-----------------------------------------------------------------------

xa  (nmin:nmax+1)            = y_el (jmin:jmax+1,ji_ray,jk_ray)
dx  (nmin:nmax)              = dy_cl(jmin:jmax  ,ji_ray,jk_ray)
xa_r(nmin:nmax+1)            = y_el (jmin:jmax+1,ji_ray,jk_ray)
dx_r(nmin:nmax)              = dy_cl(jmin:jmax  ,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Put state variables into 1D arrays, padding with 6 ghost zones
!-----------------------------------------------------------------------

rho_m (jmin:jmax)            = rhop(jmin:jmax,ji_ray,jk_ray)
r     (nmin:nmax)            = rhop(jmin:jmax,ji_ray,jk_ray)
r_r   (nmin:nmax)            = rhop(jmin:jmax,ji_ray,jk_ray)
rho_i (nmin:nmax)            = rhop(jmin:jmax,ji_ray,jk_ray)
t_m   (jmin:jmax)            = tp  (jmin:jmax,ji_ray,jk_ray)
temp  (nmin:nmax)            = tp  (jmin:jmax,ji_ray,jk_ray)
temp_r(nmin:nmax)            = tp  (jmin:jmax,ji_ray,jk_ray)
ye_m  (jmin:jmax)            = yep (jmin:jmax,ji_ray,jk_ray)
ye    (nmin:nmax)            = yep (jmin:jmax,ji_ray,jk_ray)
ye_r  (nmin:nmax)            = yep (jmin:jmax,ji_ray,jk_ray)
u     (nmin:nmax)            = vp  (jmin:jmax,ji_ray,jk_ray)
v     (nmin:nmax)            = wp  (jmin:jmax,ji_ray,jk_ray)
w     (nmin:nmax)            = up  (jmin:jmax,ji_ray,jk_ray)
ei    (nmin:nmax)            = eip (jmin:jmax,ji_ray,jk_ray)
entrop(nmin:nmax)            = aesv(jmin:jmax,3,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Load neutrino distribution function, padding with 6 ghost zones
!-----------------------------------------------------------------------

psi0_re(nmin:nmax,:,:)       = psi0p(jmin:jmax,:,:,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!  Load abundances
!-----------------------------------------------------------------------

a_nuc_rep (jmin:jmax)        = a_nuc_repp (jmin:jmax,ji_ray,jk_ray)
z_nuc_rep (jmin:jmax)        = z_nuc_repp (jmin:jmax,ji_ray,jk_ray)
be_nuc_rep(jmin:jmax)        = be_nuc_repp(jmin:jmax,ji_ray,jk_ray)

xn(jmin:jmax,ls:le)          = xnp(jmin:jmax,ls:le,ji_ray,jk_ray)

!-----------------------------------------------------------------------
!
!                      \\\\\ Y-REMAP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute boundary values
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling sweepbc, nlefty, nrighty, nmin, &
&nmax, ji_ray, jk_ray=',6i4 )") nlefty, nrighty, nmin, nmax, ji_ray, jk_ray

CALL sweepbc( nlefty, nrighty, nmin, nmax, ji_ray, jk_ray )

!-----------------------------------------------------------------------
!  Compute volume elements
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling volume, ngeomy=',i4 )") ngeomy

CALL volume ( ngeomy )

!-----------------------------------------------------------------------
!  Compute piecewise parabolic coefficents for the updated
!   Lagrangian grid
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling paraset, ntot, nmin-2, nmax+2, ngeomy=', &
& 4i4 )") ntot, nmin-2, nmax+2, ngeomy

CALL paraset( ntot, zparay, dx, xa, nmin-3, nmax+3, ngeomy )

!-----------------------------------------------------------------------
!  Prepare neutrino distribution for remap if neutrinos have been
!   advected in energy, that is, if
!
!      rhobar(i_radial) >= rhomin_y_eadvect
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' rhobar(i_radial), rhomin_y_eadvect',2es11.3)") &
& rhobar(i_radial), rhomin_y_eadvect

IF ( rhobar(i_radial) >= rhomin_y_eadvect ) THEN

  IF ( l_write ) WRITE (nlog,"(' Calling pre_remap_psi, nlefty, nrighty, &
&  jmin, jmax, nnu=',5i4 )") nlefty, nrighty, jmin, jmax, nnu

  CALL pre_remap_psi( nlefty, nrighty, jmin, jmax, nnu )

!-----------------------------------------------------------------------
!  Remap neutrinos
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling remap_psi_y, jmin, jmax, ji_ray,   &
&  jk_ray, ny, nez, nnu, time, t_bounce, tb_dy_shift=',7i4,3es11.3 )") jmin, &
&  jmax, ji_ray, jk_ray, ny, nez, nnu, time, t_bounce, tb_dy_shift

  CALL remap_psi_y( ngeomy, jmin, jmax, ji_ray, jk_ray, ny, nez, nnu,        &
&  j_shock, time, t_bounce, tb_dy_shift )

END IF ! rhobar(i_radial) >= rhomin_y_eadvect

!-----------------------------------------------------------------------
!  Mark zones for shock stabilixation
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling shock_smooth_y, jmin, jmax, ji_ray,   &
&jk_ray, ny, j_ray_dim, ik_ray_dim=',7i4 )") jmin, jmax, ji_ray, jk_ray, ny, &
& j_ray_dim, ik_ray_dim

CALL shock_smooth_y( jmin, jmax, ji_ray, jk_ray, ny, j_ray_dim, ik_ray_dim,  &
& flat_x_y, rhop, j_shock )

!-----------------------------------------------------------------------
!  Prepare composition variables for remap
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling pre_remap_comp_y, nlefty, nrighty, &
&jmin, jmax, ji_ray, jk_ray, ldim=',7i4 )") nlefty, nrighty, jmin, jmax,  &
& ji_ray, jk_ray, ldim

CALL pre_remap_comp_y( nlefty, nrighty, jmin, jmax, rho_m, t_m, ye_m,     &
& ji_ray, jk_ray, ldim )

!-----------------------------------------------------------------------
!  Remap composition
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling remap_comp_y, ngeomy, jmin, jmax, &
&ji_ray, jk_ray, ny, ldim =',6i4 )") ngeomy, jmin, jmax, ji_ray, jk_ray, &
& ny, ldim

CALL remap_comp_y( ngeomy, jmin, jmax, ji_ray, jk_ray, ny, ldim )

!-----------------------------------------------------------------------
!  Remap variables
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling remap_y, ngeomy, jmin, jmax, ji_ray,  &
&jk_ray, ny, nez, nnu, ldim =',9i4 )") ngeomy, jmin, jmax, ji_ray, jk_ray, &
& ny, nez, nnu, ldim

CALL remap_y( ngeomy, jmin, jmax, ji_ray, jk_ray, ny, nez, nnu, ldim, j_shock )

!-----------------------------------------------------------------------
!  Update the temperatures taking into account possible composition
!   changes for matter not in NSE
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling tgvndeye_sweep_comp_y, nmin, nmax, &
&ji_ray, jk_ray =',4i4 )") nmin, nmax, ji_ray, jk_ray

CALL tgvndeye_sweep_comp_y( nmin, nmax, ji_ray, jk_ray, r, rho_i )

!-----------------------------------------------------------------------
!
!             \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return updated state variables
!-----------------------------------------------------------------------

rhop(jmin:jmax,ji_ray,jk_ray) = r   (nmin:nmax)
tp  (jmin:jmax,ji_ray,jk_ray) = temp(nmin:nmax)
yep (jmin:jmax,ji_ray,jk_ray) = ye  (nmin:nmax)
vp  (jmin:jmax,ji_ray,jk_ray) = u   (nmin:nmax)
wp  (jmin:jmax,ji_ray,jk_ray) = v   (nmin:nmax)
up  (jmin:jmax,ji_ray,jk_ray) = w   (nmin:nmax)
eip (jmin:jmax,ji_ray,jk_ray) = ei  (nmin:nmax)

!-----------------------------------------------------------------------
!  Return neutrino distribution
!-----------------------------------------------------------------------

psi0p(jmin:jmax,:,:,ji_ray,jk_ray)   = psi0_re(nmin:nmax,:,:)

!-----------------------------------------------------------------------
!  Return abundances
!-----------------------------------------------------------------------

a_nuc_repp (jmin:jmax,ji_ray,jk_ray) = a_nuc_rep (jmin:jmax)
z_nuc_repp (jmin:jmax,ji_ray,jk_ray) = z_nuc_rep (jmin:jmax)
be_nuc_repp(jmin:jmax,ji_ray,jk_ray) = be_nuc_rep(jmin:jmax)

xnp(jmin:jmax,ls:le,ji_ray,jk_ray)   = xn(jmin:jmax,ls:le)

RETURN
END SUBROUTINE remap_y_inout
