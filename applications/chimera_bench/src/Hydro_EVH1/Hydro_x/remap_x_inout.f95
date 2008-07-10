SUBROUTINE remap_x_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nez, nnu, ls, le, ldim, x_el, dx_cl, x_cl, x_ef, dx_cf, &
& x_cf, rhop, tp, yep, eip, up, vp, wp, psi0p, psi1p, xnp, a_nuc_repp, &
& z_nuc_repp, be_nuc_repp, agr_c, e_bind_zn_c, eb_c, fluxbe_c )
!-----------------------------------------------------------------------
!
!    File:         remap_x_inout
!    Module:       remap_x_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         9/24/04
!
!    Purpose:
!      To receive the EVH1 arrays from radial_ray_module, execute the
!       x-array remap, and return to radhyd_angular_module the updated
!       variables
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  imin                  : lower x-array index
!  imax                  : upper x-array index
!  nx                    : x-array dimension
!  ij_ray                : index denoting the j-index of a specific radial ray
!  ik_ray                : index denoting the k-index of a specific radial ray
!  ij_ray_dim            : number of y-zones on a processor before swapping with y
!  ik_ray_dim            : number of z-zones on a processor before swapping with z
!  nez                   : neutrino energy array extent
!  nnu                   : neutrino flavor array extent
!  x_el                  : x grid zone left interfaces after Lagrangian update
!  dx_cl                 : x_el(i+1,ij_ray, ik_ray) - x_el(i,ij_ray, ik_ray)
!  x_cl                  : x grid zone midpoints after Lagrangian update
!  x_ef                  : final x grid zone left interfaces
!  dx_cf                 : final x_ef(i+1) - x_ef(i)
!  x_cf                  : final x grid zone midpoints
!  rhop                  : densities (g cm^{-3}) after Lagrangian update
!  tp                    : temperatures (K) after Lagrangian update
!  yep                   : electron fractions after Lagrangian update
!  eip                   : internal energies (ergs g^{-1}) after Lagrangian update
!  up                    : velocity x-components (cm s^{-1}) after Lagrangian update
!  vp                    : velocity y-components (cm s^{-1}) after Lagrangian update
!  wp                    : velocity z-components (cm s^{-1}) after Lagrangian update
!  psi0p                 : zero moment of the neutrino distribution function after Lagrangian update
!  psi1p                 : first moment of the neutrino distribution function after Lagrangian update
!  xnp                   : mass fractions of nuclei not in NSE after Lagrangian update
!  a_nuc_repp            : mass number of auxiliar nucleus after Lagrangian update
!  z_nuc_repp            : charge number of auxiliar nucleus after Lagrangian update
!  be_nuc_repp           : binding energy of auxiliar nucleus after Lagrangian update
!  agr_c                 : zone-centered lapse function
!
!    Output arguments:
!  rhop                  : updated densities (g cm^{-3})
!  tp                    : updated temperatures (K)
!  yep                   : updated electron fractions
!  eip                   : updated internal energies (ergs g^{-1})
!  up                    : updated velocity x-components (cm s^{-1})
!  vp                    : updated velocity y-components (cm s^{-1})
!  wp                    : updated velocity z-components (cm s^{-1})
!  psi0p                 : updated zero moment of the neutrino distribution function
!  psi1p                 : updated first moment of the neutrino distribution function
!  xnp                   : updated mass fractions of nuclei not in NSE
!  a_nuc_repp            : updated mass number of auxiliar nucleus
!  z_nuc_repp            : updated charge number of auxiliar nucleus
!  be_nuc_repp           : updated binding energy of auxiliar nucleus
!  e_bind_zn_c           : total nuclear binding energy initially in each mass shell (ergs)
!  eb_c                  : nuclear binding energy (ergs g^{-1})
!  fluxbe_c              : total nuclear binding energy transferred during remap (ergs)
!
!    Subprograms called:
!  sweepbc               : inserts boundary values into ghost zones
!  volume                : compute zone volumes
!  paraset               : computes the PPM coefficients for the new grid
!  pre_remap_psi         : sets up boundary conditions for the psi0 remap
!  remap_psi_x           : remaps the psi0's along the x-aray
!  pre_remap_comp_x      : sets up boundary conditions for the composition remap
!  remap_comp_x          : remaps the composition and the binding energies
!  remap_x               : remaps the hydro variables
!  tgvndeye_sweep_comp_x : computes the temperature given rho, e, and ye
!      
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, evh1_global, evh1_sweep,
!  evh1_zone, mgfld_remap_module, nucbrn_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, half
USE physcnst_module, ONLY : cvel

USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY: aesv, xn_e=>xn, a_nuc_rep_e=>a_nuc_rep, &
& z_nuc_rep_e=>z_nuc_rep, be_nuc_rep_e=>be_nuc_rep
USE evh1_global, ONLY : ngeomx, nleftx, nrightx, lagrangian
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, temp, u, v, w, ei, ye, xa0, &
& dx0, xa, dx, entrop, lapse_c
USE evh1_zone, ONLY : imax_z=>imax, zparax
USE mgfld_remap_module, ONLY : psi0_re, r_r=>r, temp_r=>temp, ye_r=>ye, &
& xa_r=>xa, dx_r=>dx, xa0_r=>xa0, dx0_r=>dx0, fluxbe, e_bind_zn0, eb
USE nucbrn_module, ONLY : xn, a_nuc_rep, z_nuc_rep, be_nuc_rep
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin             ! minimum x-array index
INTEGER, INTENT(in)              :: imax             ! maximum x-array index
INTEGER, INTENT(in)              :: nx               ! x-array extent

INTEGER, INTENT(in)              :: ij_ray           ! iindex denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray           ! iindex denoting the k-index of a specific radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping

INTEGER, INTENT(in)              :: ls               ! minimum composition index
INTEGER, INTENT(in)              :: le               ! maximum composition index
INTEGER, INTENT(in)              :: ldim             ! composition array extent

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)             :: x_ef           ! Eulerian x-coordinate zone edge
REAL(KIND=double), INTENT(in), DIMENSION(nx)               :: dx_cf          ! Eulerian x-coordinate zone thickness
REAL(KIND=double), INTENT(in), DIMENSION(nx)               :: x_cf           ! Eulerian x-coordinate zone center

REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)   :: x_el           ! x-coordinate zone edge after Lagrangian step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: dx_cl          ! x-coordinate zone thickness after Lagrangian step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: x_cl           ! x-coordinate zone center after Lagrangian step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)     :: agr_c          ! zone-centered lapse function

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: e_bind_zn_c    ! total nuclear binding energy initially in each mass shell (ergs)
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: eb_c           ! nuclear binding energy (ergs g^{-1})
REAL(KIND=double), INTENT(out), DIMENSION(nx,ij_ray_dim,ik_ray_dim)    :: fluxbe_c       ! total nuclear binding energy transferred during remap (ergs)

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: rhop           ! density (cm^{-3})
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: tp             ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: yep            ! electron fraction
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: eip            ! internal energy (ergs/g)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: up             ! x (radial) velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: vp             ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: wp             ! z-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi0p   ! zero moment of the neutrino occupation probability
REAL(KIND=double), INTENT(inout), DIMENSION(nx,nez,nnu,ij_ray_dim,ik_ray_dim) :: psi1p   ! first moment of the neutrino occupation probability

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: a_nuc_repp     ! nuclear mass number of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: z_nuc_repp     ! nuclear charge number of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)  :: be_nuc_repp    ! binding energy of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ldim,ij_ray_dim,ik_ray_dim) :: xnp        ! mass fractions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                  :: l_write = .false.

INTEGER                             :: i_nnse       ! lowest non-nse radial zone (1,imax)
INTEGER                             :: ntot         ! nmax + 6
INTEGER                             :: jr_min       ! minimum shifted radial index
INTEGER                             :: jr_max       ! maximum shifted radial index

REAL(KIND=double), DIMENSION(nx)    :: rho_m        ! density (g cm^{-3})
REAL(KIND=double), DIMENSION(nx)    :: t_m          ! temperature (MeV)
REAL(KIND=double), DIMENSION(nx)    :: ye_m         ! electron fraction
REAL(KIND=double), DIMENSION(nx+12) :: rho_i        ! density after the Lagrangian step (g cm^{-3})

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if lagrangian = 'ye'
!-----------------------------------------------------------------------

IF ( lagrangian ) RETURN

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rho_i                       = zero

!-----------------------------------------------------------------------
!
!                   \\\\\ SETUP FOR X-REMAP /////
!
!        Load variables received from radial_ray_module into
!         evh1_sweep and mgfld_remap_module
!
!        Subroutine remap uses the variables stored in evh1_sweep.
!        Subroutine remap_comp_x uses the variables stored in nucbrn_module
!         and mgfld_remap_module.
!        Subroutine remap_psi uses the variables stored in mgfld_remap_module
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute boundary indices of padded zones
!-----------------------------------------------------------------------

nmin                        = imin + 6
nmax                        = imax + 6
ntot                        = nmax + 6
imax_z                      = imax
jr_min                      = imin + 1
jr_max                      = imax + 1

!-----------------------------------------------------------------------
!  Load final coordinates
!-----------------------------------------------------------------------

xa0  (nmin:nmax+1)          = x_ef(imin:imax+1)
dx0  (nmin:nmax)            = dx_cf(imin:imax)
xa0_r(nmin:nmax+1)          = x_ef(imin:imax+1)
dx0_r(nmin:nmax)            = dx_cf(imin:imax)

!-----------------------------------------------------------------------
!  Load coordinates resulting from the Lagrangian step
!-----------------------------------------------------------------------

xa  (nmin:nmax+1)           = x_el (imin:imax+1,ij_ray,ik_ray)
dx  (nmin:nmax)             = dx_cl(imin:imax,ij_ray,ik_ray)
xa_r(nmin:nmax+1)           = x_el (imin:imax+1,ij_ray,ik_ray)
dx_r(nmin:nmax)             = dx_cl(imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Put state variables into 1D arrays, padding with 6 ghost zones
!-----------------------------------------------------------------------

rho_m (jr_min:jr_max)       = rhop   (imin:imax,ij_ray,ik_ray)
r     (nmin:nmax)           = rhop   (imin:imax,ij_ray,ik_ray)
r_r   (nmin:nmax)           = rhop   (imin:imax,ij_ray,ik_ray)
rho_i (nmin:nmax)           = rhop   (imin:imax,ij_ray,ik_ray)
t_m   (jr_min:jr_max)       = tp     (imin:imax,ij_ray,ik_ray)
temp  (nmin:nmax)           = tp     (imin:imax,ij_ray,ik_ray)
temp_r(nmin:nmax)           = tp     (imin:imax,ij_ray,ik_ray)
ye_m  (jr_min:jr_max)       = yep    (imin:imax,ij_ray,ik_ray)
ye    (nmin:nmax)           = yep    (imin:imax,ij_ray,ik_ray)
ye_r  (nmin:nmax)           = yep    (imin:imax,ij_ray,ik_ray)
u     (nmin:nmax)           = up     (imin:imax,ij_ray,ik_ray)
v     (nmin:nmax)           = vp     (imin:imax,ij_ray,ik_ray)
w     (nmin:nmax)           = wp     (imin:imax,ij_ray,ik_ray)
ei    (nmin:nmax)           = eip    (imin:imax,ij_ray,ik_ray)
entrop(nmin:nmax)           = aesv   (imin+1:imax+1,3,ij_ray,ik_ray)
lapse_c(nmin:nmax)          = agr_c  (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Load abundances
!-----------------------------------------------------------------------

a_nuc_rep   (jr_min:jr_max) = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep   (jr_min:jr_max) = z_nuc_repp (imin:imax,ij_ray,ik_ray)
be_nuc_rep  (jr_min:jr_max) = be_nuc_repp(imin:imax,ij_ray,ik_ray)
a_nuc_rep_e (jr_min:jr_max) = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep_e (jr_min:jr_max) = z_nuc_repp (imin:imax,ij_ray,ik_ray)
be_nuc_rep_e(jr_min:jr_max) = be_nuc_repp(imin:imax,ij_ray,ik_ray)

xn  (jr_min:jr_max,ls:le)   = xnp(imin:imax,ls:le,ij_ray,ik_ray)
xn_e(jr_min:jr_max,ls:le)   = xnp(imin:imax,ls:le,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!
!                      \\\\\ X-REMAP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute boundary values
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling sweepbc, nleftx, nrightx, nmin, &
&nmax, ij_ray, ik_ray=',6i4 )") nleftx, nrightx, nmin, nmax, ij_ray, ik_ray

CALL sweepbc( nleftx, nrightx, nmin, nmax, ij_ray, ik_ray )

!-----------------------------------------------------------------------
!  Compute volume elements
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling volume, ngeomx=',i4 )") ngeomx

CALL volume( ngeomx )

!-----------------------------------------------------------------------
!  Compute piecewise parabolic coefficents for the updated
!   Lagrangian grid
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling paraset, ntot, nmin-2, nmax+2, ngeomx=', &
& 4i4 )") ntot, nmin-2, nmax+2, ngeomx

CALL paraset( ntot, zparax, dx, xa, nmin-3, nmax+3, ngeomx )

!-----------------------------------------------------------------------
!  Prepare composition variables for remap
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling pre_remap_psi, nleftx, nrightx,    &
&imin, imax, ij_ray, ik_ray, i_nnse, ldim=',8i4 )") nleftx, nrightx, imin, &
& imax, ij_ray, ik_ray, i_nnse, ldim

CALL pre_remap_comp_x( nleftx, nrightx, imin, imax, rho_m, t_m, ye_m, &
& ij_ray, ik_ray, i_nnse, ldim )

!-----------------------------------------------------------------------
!  Remap composition
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling remap_comp_x, ngeomx, i_nnse,       &
&imin, imax, ij_ray, ik_ray, nx, ldim=',8i4 )") ngeomx, i_nnse, imin, imax, &
& ij_ray, ik_ray, nx, ldim

CALL remap_comp_x( ngeomx, i_nnse, imin, imax, ij_ray, ik_ray, nx, ldim )

!-----------------------------------------------------------------------
!  Remap variables
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling remap_x, ngeomx, imin, imax, ij_ray, &
& ik_ray, nx, nez, nnu, ldim=',8i4 )") ngeomx, imin, imax, ij_ray, ik_ray,  &
& nx, nez, nnu, ldim

CALL remap_x( ngeomx, imin, imax, ij_ray, ik_ray, nx, nez, nnu, ldim )

!-----------------------------------------------------------------------
!  Load neutrino distribution function, padding with 6 ghost zones
!  Transform to the rest frame using remaped velocities
!-----------------------------------------------------------------------

psi0_re(imin+6:imax+6,:,:)  = psi0p(imin:imax,:,:,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Prepare neutrino distribution for remap
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling pre_remap_psi, nleftx, nrightx, &
& imin, imax, nnu',5i4 )") nleftx, nrightx, imin, imax, nnu

CALL pre_remap_psi( nleftx, nrightx, imin, imax, nnu )

!-----------------------------------------------------------------------
!  Remap neutrinos
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling remap_psi_x, ngeomx, imin, imax,     &
& ij_ray, ik_ray, nx, nez, nnu=',8i4 )") ngeomx, imin, imax, ij_ray, ik_ray, &
& nx, nez, nnu

CALL remap_psi_x( ngeomx, imin, imax, ij_ray, ik_ray, nx, nez, nnu )

!-----------------------------------------------------------------------
!  Update the temperatures taking into account possible composition
!   changes for matter not in NSE
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling tgvndeye_sweep_comp_x, nmin, nmax, &
&ij_ray, ik_ray =',4i4 )") nmin, nmax, ij_ray, ik_ray

CALL tgvndeye_sweep_comp_x( nmin, nmax, ij_ray, ik_ray, r, rho_i )

!-----------------------------------------------------------------------
!
!             \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return updated state variables
!-----------------------------------------------------------------------

rhop(imin:imax,ij_ray, ik_ray) = r   (nmin:nmax)
tp  (imin:imax,ij_ray, ik_ray) = temp(nmin:nmax)
yep (imin:imax,ij_ray, ik_ray) = ye  (nmin:nmax)
up  (imin:imax,ij_ray, ik_ray) = u   (nmin:nmax)
vp  (imin:imax,ij_ray, ik_ray) = v   (nmin:nmax)
wp  (imin:imax,ij_ray, ik_ray) = w   (nmin:nmax)
eip (imin:imax,ij_ray, ik_ray) = ei  (nmin:nmax)

!-----------------------------------------------------------------------
!  Return neutrino distribution
!-----------------------------------------------------------------------

psi0p(imin:imax,:,:,ij_ray,ik_ray) = psi0_re(nmin:nmax,:,:)

!-----------------------------------------------------------------------
!  Return abundances
!-----------------------------------------------------------------------

a_nuc_repp (imin:imax,ij_ray,ik_ray) = a_nuc_rep (jr_min:jr_max)
z_nuc_repp (imin:imax,ij_ray,ik_ray) = z_nuc_rep (jr_min:jr_max)
be_nuc_repp(imin:imax,ij_ray,ik_ray) = be_nuc_rep(jr_min:jr_max)

xnp(imin:imax,ls:le,ij_ray,ik_ray)   = xn(jr_min:jr_max,ls:le)
!-----------------------------------------------------------------------
!  Return binding energies and fluxes for use in energy remap
!-----------------------------------------------------------------------

e_bind_zn_c(imin:imax,ij_ray,ik_ray) = e_bind_zn0(nmin:nmax)
eb_c       (imin:imax,ij_ray,ik_ray) = eb        (nmin:nmax)
fluxbe_c   (imin:imax,ij_ray,ik_ray) = fluxbe    (nmin:nmax)


RETURN
END SUBROUTINE remap_x_inout
