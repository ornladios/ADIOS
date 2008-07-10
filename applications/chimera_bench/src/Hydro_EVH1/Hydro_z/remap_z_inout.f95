SUBROUTINE remap_z_inout( nx, nz, kmin, kmax, ki_ray, kj_ray, ij_ray_dim, &
& k_ray_dim, i_radial, nez, nnu, ls, le, ldim, z_el, dz_cl, z_cl, z_ef,   &
& dz_cf, z_cf, rhop, tp, yep, eip, up, vp, wp, psi0p, xnp, a_nuc_repp,    &
& z_nuc_repp, be_nuc_repp, flat_x_z, time, t_bounce, tb_dy_shift, rhobar )
!-----------------------------------------------------------------------
!
!    File:         remap_z_inout
!    Module:       remap_z_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/12/07
!
!    Purpose:
!      To port angular_ray_module variables into, and updated variabbles
!       out of, the z-array remap modules via subroutine remap_z_inout.
!
!    Input arguments:
!  nx             : x-array extent
!  nz             : z-array extent
!  kmin           : lower y-array index
!  kmax           : upper y-array index
!  ki_ray         : x (radial) index of a specific z (azimuthal) ray
!  kj_ray         : y (azimuthal) index of a specific z (azimuthal) ray
!  ij_ray_dim     : the number of radial zones on a processor before swapping with y
!  k_ray_dim      : the number of z-zones on a processor after swapping with z
!  i_radial       : the unshifted radial zone (angular ray) corresponding to ki_ray, kj_ray
!  nez            : neutrino energy array extent
!  nnu            : neutrino flavor array extent
!  z_el           : z grid zone left interfaces after Lagrangian update
!  dz_cl          : z_el(j+1) - z_el(j)
!  z_cl           : z grid zone midpoints after Lagrangian update
!  z_ef           : final z grid zone left interfaces
!  dz_cf          : z_ef(j+1) - z_ef(j)
!  z_cf           : final z grid zone midpoints
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
!  flat_x_z       : variables indicating the presence of radial shocks
!  time           : elapsed time
!  t_bounce       : time from bounce
!  tb_dy_shift    : time from bounce to turn off grid wiggle and/or psi0 diffusion
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
!  sweepbc               : inserts boundary values into ghost zones
!  volume                : compute zone volumes
!  paraset               : computes the PPM coefficients for the new grid
!  pre_remap_psi         : sets up boundary conditions for the psi0 remap
!  remap_psi_z           : remaps the psi0's along the z-aray
!  shock_smooth_z        : marks zones lying along a shock for additional diffusion
!  pre_remap_comp_z      : sets up boundary conditions for the composition remap
!  remap_comp_z          : remaps the composition and the binding energies
!  remap_z               : remaps the hydro variables
!  tgvndeye_sweep_comp_z : computes the temperature given rho, e, and ye
!      
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, edit_module eos_snc_z_module, evh1_global,
!  evh1_sweep, evh1_zone, mgfld_remap_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero

USE e_advct_module, ONLY : rhomin_z_eadvect
USE edit_module, ONLY : nlog
USE eos_snc_z_module, ONLY : xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, nse, aesv
USE evh1_global, ONLY : ngeomz, nleftz, nrightz
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, temp, u, v, w, ei, ye, xa0, &
& dx0, xa, dx, entrop
USE evh1_zone, ONLY : zparay
USE mgfld_remap_module, ONLY : psi0_re, r_r=>r, temp_r=>temp, ye_r=>ye, &
& xa_r=>xa, dx_r=>dx, xa0_r=>xa0, dx0_r=>dx0, comp
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nx               ! x-array extent

INTEGER, INTENT(in)              :: kmin             ! minimum y-array index
INTEGER, INTENT(in)              :: kmax             ! maximum y-array index
INTEGER, INTENT(in)              :: nz               ! z-array extent

INTEGER, INTENT(in)              :: ki_ray           ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray           ! y (azimuthal) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of radial zones on a processor before swapping with y
INTEGER, INTENT(in)              :: k_ray_dim        ! number of radial zones on a processor after swapping with z
INTEGER, INTENT(in)              :: i_radial         ! the unshifted radial zone corresponding to ki_ray, kj_ray

INTEGER, INTENT(in)              :: ls               ! minimum composition index
INTEGER, INTENT(in)              :: le               ! maximum composition index
INTEGER, INTENT(in)              :: ldim             ! composition array extent

INTEGER, INTENT(in)              :: nez              ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu              ! neutrino energy flavor extent

REAL(KIND=double), INTENT(in),    DIMENSION(nz+1)                            :: z_ef        ! Eulerian z-coordinate zone edge
REAL(KIND=double), INTENT(in),    DIMENSION(nz)                              :: dz_cf       ! Eulerian z-coordinate zone thickness
REAL(KIND=double), INTENT(in),    DIMENSION(nz)                              :: z_cf        ! Eulerian z-coordinate zone center
REAL(KIND=double), INTENT(in),    DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: flat_x_z    ! variables indicating the presence of radial shocks
REAL(KIND=double), INTENT(in)    :: time             ! elapsed time
REAL(KIND=double), INTENT(in)    :: t_bounce         ! time from bounce
REAL(KIND=double), INTENT(in)    :: tb_dy_shift      ! time from bounce to trun off grid wiggle and/or psi0 diffusion
REAL(KIND=double), INTENT(in),    DIMENSION(nx)                              :: rhobar      ! mean density at a given radius

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nz+1,ij_ray_dim,k_ray_dim)       :: z_el        ! z-coordinate zone edge after Lagrangian step
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: dz_cl       ! z-coordinate zone thickness after Lagrangian step
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: z_cl        ! z-coordinate zone center after Lagrangian step

REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: rhop        ! density (cm^{-3})
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: tp          ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: yep         ! temperature (MeV)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: eip         ! internal energy (ergs/g)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: up          ! radial velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: vp          ! y-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: wp          ! z-velocity of zone (cm)
REAL(KIND=double), INTENT(inout), DIMENSION(nz,nez,nnu,ij_ray_dim,k_ray_dim) :: psi0p       ! zero moment of the neutrino occupation probability

REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: a_nuc_repp  ! nuclear mass number of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: z_nuc_repp  ! nuclear charge number of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ij_ray_dim,k_ray_dim)         :: be_nuc_repp ! binding energy of auxiliary nucleus
REAL(KIND=double), INTENT(inout), DIMENSION(nz,ldim,ij_ray_dim,k_ray_dim)    :: xnp         ! mass fractions

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                  :: l_write = .false.

INTEGER                             :: ntot           ! nmax + 6
INTEGER, DIMENSION(nz)              :: k_shock        ! zones marked for added y-diffusion

REAL(KIND=double), DIMENSION(nz)    :: rho_m          ! density (g cm^{-3})
REAL(KIND=double), DIMENSION(nz)    :: t_m            ! temperature (MeV)
REAL(KIND=double), DIMENSION(nz)    :: ye_m           ! electron fraction
REAL(KIND=double), DIMENSION(nz+12) :: rho_i          ! density after the Lagrangian step (g cm^{-3})

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
!        Subroutine remap_comp_z uses the variables stored in nucbrn_module
!         and mgfld_remap_module.
!        Subroutine remap_psi uses the variables stored in mgfld_remap_module
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute boundary indices of padded zones
!-----------------------------------------------------------------------

nmin                         = kmin + 6
nmax                         = kmax + 6
ntot                         = nmax + 6

!-----------------------------------------------------------------------
!  Load final coordinates
!-----------------------------------------------------------------------

xa0(nmin:nmax+1)             = z_ef(kmin:kmax+1)
dx0(nmin:nmax)               = dz_cf(kmin:kmax)
xa0_r(nmin:nmax+1)           = z_ef(kmin:kmax+1)
dx0_r(nmin:nmax)             = dz_cf(kmin:kmax)

!-----------------------------------------------------------------------
!  Load coordinates resulting from the Lagrangian step
!-----------------------------------------------------------------------

xa  (nmin:nmax+1)            = z_el (kmin:kmax+1,kj_ray,ki_ray)
dx  (nmin:nmax)              = dz_cl(kmin:kmax  ,kj_ray,ki_ray)
xa_r(nmin:nmax+1)            = z_el (kmin:kmax+1,kj_ray,ki_ray)
dx_r(nmin:nmax)              = dz_cl(kmin:kmax  ,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Put state variables into 1D arrays, padding with 6 ghost zones
!-----------------------------------------------------------------------

rho_m (kmin:kmax)            = rhop(kmin:kmax,kj_ray,ki_ray)
r     (nmin:nmax)            = rhop(kmin:kmax,kj_ray,ki_ray)
r_r   (nmin:nmax)            = rhop(kmin:kmax,kj_ray,ki_ray)
rho_i (nmin:nmax)            = rhop(kmin:kmax,kj_ray,ki_ray)
t_m   (kmin:kmax)            = tp  (kmin:kmax,kj_ray,ki_ray)
temp  (nmin:nmax)            = tp  (kmin:kmax,kj_ray,ki_ray)
temp_r(nmin:nmax)            = tp  (kmin:kmax,kj_ray,ki_ray)
ye_m  (kmin:kmax)            = yep (kmin:kmax,kj_ray,ki_ray)
ye    (nmin:nmax)            = yep (kmin:kmax,kj_ray,ki_ray)
ye_r  (nmin:nmax)            = yep (kmin:kmax,kj_ray,ki_ray)
u     (nmin:nmax)            = wp  (kmin:kmax,kj_ray,ki_ray)
v     (nmin:nmax)            = up  (kmin:kmax,kj_ray,ki_ray)
w     (nmin:nmax)            = vp  (kmin:kmax,kj_ray,ki_ray)
ei    (nmin:nmax)            = eip (kmin:kmax,kj_ray,ki_ray)
entrop(nmin:nmax)            = aesv(kmin:kmax,3,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Load neutrino distribution function, padding with 6 ghost zones
!-----------------------------------------------------------------------

psi0_re(nmin:nmax,:,:)       = psi0p(kmin:kmax,:,:,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!  Load abundances
!-----------------------------------------------------------------------

a_nuc_rep (kmin:kmax)        = a_nuc_repp (kmin:kmax,kj_ray,ki_ray)
z_nuc_rep (kmin:kmax)        = z_nuc_repp (kmin:kmax,kj_ray,ki_ray)
be_nuc_rep(kmin:kmax)        = be_nuc_repp(kmin:kmax,kj_ray,ki_ray)

xn(kmin:kmax,ls:le)          = xnp(kmin:kmax,ls:le,kj_ray,ki_ray)

!-----------------------------------------------------------------------
!
!                      \\\\\ Y-REMAP /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute boundary values
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling sweepbc, nleftz, nrightz, nmin, &
&nmax, ki_ray, kj_ray=',6i4 )") nleftz, nrightz, nmin, nmax, ki_ray, kj_ray

CALL sweepbc( nleftz, nrightz, nmin, nmax, ki_ray, kj_ray )

!-----------------------------------------------------------------------
!  Compute volume elements
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling volume, ngeomz=',i4 )") ngeomz

CALL volume ( ngeomz )

!-----------------------------------------------------------------------
!  Compute piecewise parabolic coefficents for the updated
!   Lagrangian grid
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling paraset, ntot, nmin-2, nmax+2, ngeomz=', &
& 4i4 )") ntot, nmin-2, nmax+2, ngeomz

CALL paraset( ntot, zparay, dx, xa, nmin-3, nmax+3, ngeomz )

!-----------------------------------------------------------------------
!  Prepare neutrino distribution for remap if neutrinos have been
!   advected in energy, that is, if
!
!      rhobar(i_radial) >= rhomin_z_eadvect
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' rhobar(i_radial), rhomin_z_eadvect',2es11.3)") &
& rhobar(i_radial), rhomin_z_eadvect

IF ( rhobar(i_radial) >= rhomin_z_eadvect ) THEN

  IF ( l_write ) WRITE (nlog,"(' Calling pre_remap_psi, nleftz, nrightz, &
&  kmin, kmax, nnu=',5i4 )") nleftz, nrightz, kmin, kmax, nnu

  CALL pre_remap_psi( nleftz, nrightz, kmin, kmax, nnu )

!-----------------------------------------------------------------------
!  Remap neutrinos
!-----------------------------------------------------------------------

  IF ( l_write )  WRITE (nlog,"(' Calling remap_psi_z, kmin, kmax, ki_ray, &
&  kj_ray, nz, nez, nnu, time, t_bounce, tb_dy_shift=',7i4,3es11.3 )") kmin, &
&  kmax, ki_ray, kj_ray, nz, nez, nnu, time, t_bounce, tb_dy_shift

  CALL remap_psi_z( ngeomz, kmin, kmax, ki_ray, kj_ray, nz, nez, nnu,        &
&  k_shock, time, t_bounce, tb_dy_shift )

END IF ! rhobar(i_radial) >= rhomin_z_eadvect

!-----------------------------------------------------------------------
!  Mark zones for shock stabilixation
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling shock_smooth_z, kmin, kmax, ki_ray,   &
&kj_ray, nz, ij_ray_dim, k_ray_dim=',7i4 )") kmin, kmax, ki_ray, kj_ray, nz, &
& ij_ray_dim, k_ray_dim


CALL shock_smooth_z( kmin, kmax, ki_ray, kj_ray, nz, ij_ray_dim, k_ray_dim,   &
& flat_x_z, rhop, k_shock )

!-----------------------------------------------------------------------
!  Prepare composition variables for remap
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling pre_remap_comp_z, nleftz, nrightz, &
&kmin, kmax, ki_ray, kj_ray, ldim=',7i4 )") nleftz, nrightz, kmin, kmax,  &
& ki_ray, kj_ray, ldim

CALL pre_remap_comp_z( nleftz, nrightz, kmin, kmax, rho_m, t_m, ye_m,     &
& ki_ray, kj_ray, ldim )

!-----------------------------------------------------------------------
!  Remap composition
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling remap_comp_z, ngeomz, kmin, kmax, &
&ki_ray, kj_ray, nz, ldim =',7i4 )") ngeomz, kmin, kmax, ki_ray, kj_ray, &
& nz, ldim

CALL remap_comp_z( ngeomz, kmin, kmax, ki_ray, kj_ray, nz, ldim )

!-----------------------------------------------------------------------
!  Remap variables
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling remap_z, ngeomz, kmin, kmax, ki_ray,  &
&kj_ray, nz, nez, nnu, ldim =',9i4 )") ngeomz, kmin, kmax, ki_ray, kj_ray, &
& nz, nez, nnu, ldim

CALL remap_z( ngeomz, kmin, kmax, ki_ray, kj_ray, nz, nez, nnu, ldim, &
& k_shock )

!-----------------------------------------------------------------------
!  Update the temperatures taking into account possible composition
!   changes for matter not in NSE
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling tgvndeye_sweep_comp_z, nmin, nmax, &
&ki_ray, kj_ray =',4i4 )") nmin, nmax, ki_ray, kj_ray

CALL tgvndeye_sweep_comp_z( nmin, nmax, ki_ray, kj_ray, r, rho_i )

!-----------------------------------------------------------------------
!
!             \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return updated state variables
!-----------------------------------------------------------------------

rhop(kmin:kmax,kj_ray,ki_ray)        = r   (nmin:nmax)
tp  (kmin:kmax,kj_ray,ki_ray)        = temp(nmin:nmax)
yep (kmin:kmax,kj_ray,ki_ray)        = ye  (nmin:nmax)
wp  (kmin:kmax,kj_ray,ki_ray)        = u   (nmin:nmax)
up  (kmin:kmax,kj_ray,ki_ray)        = v   (nmin:nmax)
vp  (kmin:kmax,kj_ray,ki_ray)        = w   (nmin:nmax)
eip (kmin:kmax,kj_ray,ki_ray)        = ei  (nmin:nmax)

!-----------------------------------------------------------------------
!  Return neutrino distribution
!-----------------------------------------------------------------------

psi0p(kmin:kmax,:,:,kj_ray,ki_ray)   = psi0_re(nmin:nmax,:,:)

!-----------------------------------------------------------------------
!  Return abundances
!-----------------------------------------------------------------------

a_nuc_repp (kmin:kmax,kj_ray,ki_ray) = a_nuc_rep (kmin:kmax)
z_nuc_repp (kmin:kmax,kj_ray,ki_ray) = z_nuc_rep (kmin:kmax)
be_nuc_repp(kmin:kmax,kj_ray,ki_ray) = be_nuc_rep(kmin:kmax)

xnp(kmin:kmax,ls:le,kj_ray,ki_ray)   = xn(kmin:kmax,ls:le)

RETURN
END SUBROUTINE remap_z_inout
