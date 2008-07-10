SUBROUTINE remap_x_e_inout( imin, imax, nx, ij_ray, ik_ray, ij_ray_dim, &
& ik_ray_dim, nnc, x_ci, x_el, dx_cl, x_cl, x_ef, dx_cf, x_cf, rho_l, &
& rhop, tp, yep, eip, e_vp, u_l, v_l, w_l, u_c, v_c, w_c, p_c, xnp, &
& a_nuc_repp, z_nuc_repp, be_nuc_repp, e_bind_zn_c, eb_c, fluxbe_c, &
& grav_pot_c, grav_pot_c_i )
!-----------------------------------------------------------------------
!
!    File:         remap_x_e_inout
!    Module:       remap_x_e_inout
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/27/07
!
!    Purpose:
!      To receive the EVH1 arrays from radial_ray_module, execute the
!       x-array remapping of the total energy, and return to 
!       radhyd_angular_module the updated variables
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
!  x_ci                  : x grid zone midpoints at the beginning of the time step
!  x_el                  : x grid zone left interfaces after Lagrangian update
!  dx_cl                 : x_el(i+1,ij_ray,ik_ray) - x_el(i,ij_ray,ik_ray)
!  x_cl                  : x grid zone midpoints after Lagrangian update
!  x_ef                  : final x grid zone left interfaces
!  dx_cf                 : x_ef(i+1) - x_ef(i)
!  x_cf                  : final x grid zone midpoints
!  rho_i                 : densities (g cm^{-3}) after Lagrangian update
!  rhop                  : densities (g cm^{-3}) after remap
!  yep                   : electron fractions after remap
!  e_vp                  : partial energies (ergs g^{-1}) after Lagrangian update
!  u_l                   : velocity x-components (cm s^{-1}) after Lagrangian update
!  v_l                   : velocity y-components (cm s^{-1}) after Lagrangian update
!  w_l                   : velocity z-components (cm s^{-1}) after Lagrangian update
!  u_c                   : velocity x-components (cm s^{-1}) after rampping
!  v_c                   : velocity y-components (cm s^{-1}) after rampping
!  w_c                   : velocity z-components (cm s^{-1}) after rampping
!  xnp                   : mass fractions of nuclei not in NSE after Lagrangian update
!  a_nuc_repp            : mass number of auxiliar nucleus after Lagrangian update
!  z_nuc_repp            : charge number of auxiliar nucleus after Lagrangian update
!  be_nuc_repp           : binding energy of auxiliar nucleus after Lagrangian update
!  e_bind_zn_c           : total nuclear binding energy initially in each mass shell (ergs)
!  eb_c                  : nuclear binding energy (ergs g^{-1})
!  fluxbe_c              : total nuclear binding energy transferred during remap (ergs)
!  grav_pot_c            : final zone-centered zone-centered gravitational potential (erg g^{-1})
!  grav_pot_c_i          : initial zone-centered zone-centered gravitational potential (erg g^{-1})
!
!    Output arguments:
!  tp                    : updated temperatures (K)
!  eip                   : updated internal energies (ergs g^{-1})
!  p_c                   : updated pressuress (ergs cm^{-3})
!
!    Subprograms called:
!  sweepbc               : inserts boundary values into ghost zones
!  volume                : compute zone volumes
!  paraset               : computes the PPM coefficients for the new grid
!  remap_x_e             : remaps the total energy
!  tgvndeye_sweep_comp_x : computes the temperature given rho, e, and ye
!      
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, evh1_global, evh1_sweep,
!  evh1_zone, mgfld_remap_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero, half
USE physcnst_module, ONLY : cvel

USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY: aesv, xn, a_nuc_rep, z_nuc_rep, be_nuc_rep
USE evh1_global, ONLY : ngeomx, nleftx, nrightx, lagrangian
USE evh1_sweep, ONLY: sweep, nmin, nmax, r, temp, u, v, w, ei, e_v, ye, &
& p, xa0, dx0, xa, dx, egrav, degrav
USE evh1_zone, ONLY : imax_z=>imax, zparax
USE mgfld_remap_module, ONLY : fluxbe, e_bind_zn0, eb

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: imin             ! minimum x-array index
INTEGER, INTENT(in)              :: imax             ! maximum x-array index
INTEGER, INTENT(in)              :: nx               ! x-array extent
INTEGER, INTENT(in)              :: nnc              ! composition array extent

INTEGER, INTENT(in)              :: ij_ray           ! iindex denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray           ! iindex denoting the k-index of a specific radial ray
INTEGER, INTENT(in)              :: ij_ray_dim       ! number of y-zones on a processor before swapping
INTEGER, INTENT(in)              :: ik_ray_dim       ! number of z-zones on a processor before swapping

REAL(KIND=double), INTENT(in), DIMENSION(nx)                :: x_ci           ! x grid zone midpoints at the beginning of the time step

REAL(KIND=double), INTENT(in), DIMENSION(nx+1)              :: x_ef           ! final x grid zone left interfaces
REAL(KIND=double), INTENT(in), DIMENSION(nx)                :: dx_cf          ! x_ef(i+1) - x_ef(i)
REAL(KIND=double), INTENT(in), DIMENSION(nx)                :: x_cf           ! final x grid zone midpoints

REAL(KIND=double), INTENT(in), DIMENSION(nx+1,ij_ray_dim,ik_ray_dim)    :: x_el           ! x-coordinate zone edge after Lagrangian step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: dx_cl          ! x-coordinate zone thickness after Lagrangian step
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: x_cl           ! x-coordinate zone center after Lagrangian step

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: rho_l          ! densities (cm^{-3}) after Lagrangian update
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: rhop           ! densities (cm^{-3}) after remap
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: yep            ! electron fractions after remap
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: e_vp           ! partial energies after Lagrangian update

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: e_bind_zn_c    ! total nuclear binding energy initially in each mass shell (ergs)
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: eb_c           ! nuclear binding energy (ergs g^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: fluxbe_c       ! total nuclear binding energy transferred during remap (ergs)

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: grav_pot_c     ! initial zone-centered zone-centered gravitational potential measured from r = 0 (erg g^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: grav_pot_c_i   ! final zone-centered zone-centered gravitational potential measured from r = 0 (erg g^{-1})

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: u_l            ! radial velocity of zone after Lagrangian update (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: v_l            ! y-velocity of zone after Lagrangian update (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: w_l            ! z-velocity of zone after Lagrangian update (cm s^{-1})

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: u_c            ! radial velocity of zone after remapping (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: v_c            ! y-velocity of zone after remapping (cm s^{-1})
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: w_c            ! z-velocity of zone after remapping (cm s^{-1})

REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: a_nuc_repp     ! nuclear mass number of auxiliary nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: z_nuc_repp     ! nuclear charge number of auxiliary nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,ij_ray_dim,ik_ray_dim)      :: be_nuc_repp    ! binding energy of auxiliary nucleus
REAL(KIND=double), INTENT(in), DIMENSION(nx,nnc,ij_ray_dim,ik_ray_dim)  :: xnp            ! mass fractions

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: eip            ! internal energy (ergs g^{-1})
REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: p_c            ! pressure (ergs cm^{-3})

!-----------------------------------------------------------------------
!        Input-Output variables
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nx,ij_ray_dim,ik_ray_dim)   :: tp             ! temperature (MeV)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL, PARAMETER                  :: l_write = .false.

INTEGER                             :: ntot         ! nmax + 6
INTEGER                             :: jr_min       ! minimum shifted radial index
INTEGER                             :: jr_max       ! maximum shifted radial index

REAL(KIND=double), DIMENSION(nx+12) :: xi           ! padded x grid zone midpoints at the beginning of the time step
REAL(KIND=double), DIMENSION(nx+12) :: xf           ! padded x grid zone midpoints of the final grid
REAL(KIND=double), DIMENSION(nx+12) :: rho_i        ! density after the Lagrangian step (g cm^{-3})
REAL(KIND=double), DIMENSION(nx+12) :: egrav_i      ! padded gravitational potential before Lgrangian step
REAL(KIND=double), DIMENSION(nx+12) :: egrav_if     ! padded gravitational potential before Lgrangian step interpolated to the final grid

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
!                   \\\\\ SETUP FOR X_E-REMAP /////
!
!  Load variables received from radial_ray_module into
!   evh1_sweep and mgfld_remap_module
!
!  Subroutine remap uses the variables stored in evh1_sweep.
!  Subroutine remap_comp_x uses the variables stored in nucbrn_module
!   and mgfld_remap_module.
!  Subroutine remap_psi uses the variables stored in mgfld_remap_module
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

!-----------------------------------------------------------------------
!  Load coordinates resulting from the Lagrangian step into 1D arrays,
!   padding with 6 ghost zones
!-----------------------------------------------------------------------

xa  (nmin:nmax+1)           = x_el (imin:imax+1,ij_ray,ik_ray)
dx  (nmin:nmax)             = dx_cl(imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Put state variables into 1D arrays, padding with 6 ghost zones
!-----------------------------------------------------------------------

r      (nmin:nmax)          = rho_l  (imin:imax,ij_ray,ik_ray)
rho_i  (nmin:nmax)          = rho_l  (imin:imax,ij_ray,ik_ray)
u      (nmin:nmax)          = u_l    (imin:imax,ij_ray,ik_ray)
v      (nmin:nmax)          = v_l    (imin:imax,ij_ray,ik_ray)
w      (nmin:nmax)          = w_l    (imin:imax,ij_ray,ik_ray)
e_v    (nmin:nmax)          = e_vp   (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Load abundances
!-----------------------------------------------------------------------

a_nuc_rep   (jr_min:jr_max) = a_nuc_repp (imin:imax,ij_ray,ik_ray)
z_nuc_rep   (jr_min:jr_max) = z_nuc_repp (imin:imax,ij_ray,ik_ray)
be_nuc_rep  (jr_min:jr_max) = be_nuc_repp(imin:imax,ij_ray,ik_ray)

xn  (jr_min:jr_max,1:nnc)   = xnp(imin:imax,1:nnc,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Load binding energies and fluxes into 1D arrays, padding with 6 ghost
!   zones
!-----------------------------------------------------------------------

e_bind_zn0(nmin:nmax)       = e_bind_zn_c(imin:imax,ij_ray,ik_ray)
eb        (nmin:nmax)       = eb_c       (imin:imax,ij_ray,ik_ray)
fluxbe    (nmin:nmax)       = fluxbe_c   (imin:imax,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!  Load gravitational potential into 1D arrays, padding with 6 ghost
!   zones
!-----------------------------------------------------------------------

egrav    (nmin:nmax)        = grav_pot_c(imin:imax,ij_ray,ik_ray)


!-----------------------------------------------------------------------
!
!                   \\\\\ TOTAL ENERGY REMAP /////
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

IF ( l_write )  WRITE (nlog,"(' Calling volume, ngeomx',i4 )") ngeomx

CALL volume( ngeomx )

!-----------------------------------------------------------------------
!  Compute piecewise parabolic coefficents for the updated
!   Lagrangian grid
!-----------------------------------------------------------------------

IF ( l_write )  WRITE (nlog,"(' Calling paraset, ntot, nmin-2, nmax+2, ngeomx=', &
& 4i4 )") ntot, nmin-2, nmax+2, ngeomx

CALL paraset( ntot, zparax, dx, xa, nmin-3, nmax+3, ngeomx )

!-----------------------------------------------------------------------
!  Remap the total energy
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling remap_x_e, ngeomx, imin, imax, &
& ij_ray, ik_ray, nx=',6i4 )") ngeomx, imin, imax, ij_ray, ik_ray, nx

CALL remap_x_e( ngeomx, imin, imax, ij_ray, ik_ray, nx )

!-----------------------------------------------------------------------
!  Put remapped state variables into 1D arrays, padding with 6 ghost
!   zones
!-----------------------------------------------------------------------

r      (nmin:nmax)          = rhop   (imin:imax,ij_ray,ik_ray)
temp   (nmin:nmax)          = tp     (imin:imax,ij_ray,ik_ray)
ye     (nmin:nmax)          = yep    (imin:imax,ij_ray,ik_ray)
u      (nmin:nmax)          = u_c    (imin:imax,ij_ray,ik_ray)
v      (nmin:nmax)          = v_c    (imin:imax,ij_ray,ik_ray)
w      (nmin:nmax)          = w_c    (imin:imax,ij_ray,ik_ray)
!-----------------------------------------------------------------------
!  Compute egrav_i, the initial gravitational potential interpolated to
!   the final grid, and compute degrav, the change in the gravitational
!   potential at constant radius
!-----------------------------------------------------------------------

xi(nmin:nmax)               = x_ci(imin:imax)
xf(nmin:nmax)               = x_cf(imin:imax)
egrav_i(nmin:nmax)          = grav_pot_c_i(imin:imax,ij_ray,ik_ray)

CALL phi_intrp( nmin, nmax, ntot, xi, xf, egrav_i, egrav_if )

degrav (nmin:nmax)          = egrav(nmin:nmax) - egrav_if(nmin:nmax)

!-----------------------------------------------------------------------
!  Extract the internal energy from the total energy
!-----------------------------------------------------------------------

CALL e_decompose_g( nmin, nmax )

!-----------------------------------------------------------------------
!  Update the temperatures taking into account possible composition
!   changes for matter not in NSE
!-----------------------------------------------------------------------

IF ( l_write ) WRITE (nlog,"(' Calling tgvndeye_sweep_x, nmin, nmax, &
&ij_ray, ik_ray =',4i4 )") nmin, nmax, ij_ray, ik_ray

CALL tgvndeye_sweep_x_rho( nmin, nmax, ij_ray, ik_ray, r, rho_i )

!-----------------------------------------------------------------------
!
!             \\\\\ RETURN UPDAATED VARIABLES ////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return updated state variables
!-----------------------------------------------------------------------

tp  (imin:imax,ij_ray,ik_ray) = temp(nmin:nmax)
eip (imin:imax,ij_ray,ik_ray) = ei  (nmin:nmax)
p_c (imin:imax,ij_ray,ik_ray) = p   (nmin:nmax)

RETURN
END SUBROUTINE remap_x_e_inout
