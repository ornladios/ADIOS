SUBROUTINE remap_comp_x( ngeom, i_nnse, imin, imax, ij_ray, ik_ray, nx, &
& nnc )
!-----------------------------------------------------------------------
!
!    File:         remap_comp_x
!    Module:       remap_comp_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/16/03
!
!    Purpose:
!      To remap the composition and binding energy along radial rays
!       using piecewise parabolic functions. No flattening is used on remap
!       (pass dummy array: dum to parabola.f).
!
!    Input arguments:
!  ngeom       : geometry index
!  i_nnse      : i-index of first zone not in nse
!  imin        : minimum x-array index
!  imax        : maximim x-array index
!  ij_ray      : index denoting the j-index of a specific radial ray
!  ik_ray      : index denoting the k-index of a specific radial ray
!  nx          : x-array extent
!  nnc         : composition array extent
!
!    Output arguments:
!        none
!
!    Subprograms called:
!
!  volume_zone : Computes volumes in the Lagrangian and Eulerian grid
!  eos_nnse_e  : Computes the mean nuclear binding energy for zones not in nse
!  parabola    : Computes piecewise parabolic fits for the densit and composition
!  eos_nnse_e  : Computes the mean binding energy for the mass being trransferred
!  deflash     : Deflashes material in material being transferred across a nse to non-nse interface
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, evh1_global, evh1_sweep,
!  evh1_zone, mgfld_remap_module, mdl_cnfg_module, nucbrn_module
!
!-----------------------------------------------------------------------
!
!    Method:
!
!    nse = 0 advected into nse = 0:
!  Find the overlapping mass fractions and advect
!
!    nse = 0 advected into nse = 1:
!  Find the overlapping mass fractions and advect out of nse = 0; flash
!   this material and advect into nse = 1
!
!    nse = 1 advected into nse = 0:
!  Find the overlapping mass fractions, deflash the material and advect
!   into nse = 0 region
!
!    nse = 1 advected into nse = 1
!  Only rho, e, and ye are advected, and this is accomplished in remap_y
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, epsilon
USE physcnst_module, ONLY: pi

USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY: nse, xn_e=>xn, a_nuc_rep_e=>a_nuc_rep, &
& z_nuc_rep_e=>z_nuc_rep, be_nuc_rep_e=>be_nuc_rep
USE evh1_global, ONLY: smallr
USE evh1_sweep, ONLY: radius, ye
USE evh1_zone, ONLY : zparax
USE mgfld_remap_module, ONLY : comp, r, temp, ye_r=>ye, xa, dx, xa0, dx0, &
& fluxbe, e_bind_zn0, eb, fluxye_comp
USE nucbrn_module, ONLY : xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, nuc_number, &
& a_nuc, z_nuc, be_nuc

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                            :: ngeom        ! geometry index
INTEGER, INTENT(in)                            :: i_nnse       ! index of first non-nse zone (1,imax)
INTEGER, INTENT(in)                            :: imin         ! minimum x-array index
INTEGER, INTENT(in)                            :: imax         ! maximim x-array index
INTEGER, INTENT(in)                            :: ij_ray       ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                            :: ik_ray       ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                            :: nx           ! x-array extent
INTEGER, INTENT(in)                            :: nnc          ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=12)                             :: var_name

LOGICAL                                        :: first = .true.

INTEGER                                        :: istat        ! allocation status flag
INTEGER                                        :: nc           ! composition index
INTEGER                                        :: n            ! padded aray index
INTEGER                                        :: k            ! k-1 padded aray index
INTEGER                                        :: j            ! MGFLD radial zone index
INTEGER                                        :: i            ! unpadded aray index
INTEGER                                        :: n_nucp1      ! nuc_number + 1
INTEGER                                        :: nmin         ! minimum padded index
INTEGER                                        :: nmax         ! maximum padded index
INTEGER                                        :: ntot         ! number of zones (real plus ghost)
INTEGER                                        :: nminc        ! minimum padded zone index for composition
INTEGER                                        :: nmincm1      ! nminc - 1

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: rl           ! density at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: r6           ! density parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dr           ! density slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: cmpl         ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: cmp6         ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dcmp         ! izero neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: compl        ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: comp6        ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: dcomp        ! izero neutrino moment slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dvol         ! volume after Lagr step
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dvol0        ! volume after Eul remap

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dm           ! mass after lagrangian step
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dm0          ! mass after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: delta        ! volume of overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: re           ! density after Eul remap

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: cmp          ! 1-D storage array for the composition array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: flux_nc_rep  ! proportional to heavy nuclei number advected
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: flux_ba_rep  ! proportional to heavy nuclei baryon number advected
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: flux_z_rep   ! proportional to heavy nuclei charge advected
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: flux_b_rep   ! proportional to binding energy advected
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: n_rep        ! proportional to initial number of rep heavy nuclei/zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: ba_rep       ! proportional to initial number of rep heavy nuclei barions/zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: n0_rep       ! proportional to the final number of rep heavy nuclei/zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: z0_rep       ! proportional to the initial total charge of rep heavy nuclei/zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: z_rep        ! proportional to the total charge of rep heavy nuclei/zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: b0_rep       ! proportional to theinitial total binding energy of rep heavy nuclei/zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: b_rep        ! proportional to the total binding energy of rep heavy nuclei/zone

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: fluxr        ! mass contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: fluxcmp      ! composition contained in overlapping subshells

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: xn_t         ! mass fractions contained in overlapping subshells
REAL(KIND=double)                              :: a_nuc_rep_t  ! mass number of specie nnc contained in overlapping subshells
REAL(KIND=double)                              :: z_nuc_rep_t  ! charge number of specie nnc contained in overlapping subshells
REAL(KIND=double)                              :: be_nuc_rep_t ! binding energy of specie nnc contained in overlapping subshells

REAL(KIND=double)                              :: fluxcmp_t    ! total mass advected computed from mass fractions
REAL(KIND=double)                              :: phi          ! ratio of mass advected to mass fractions advected

REAL(KIND=double)                              :: deltx        ! grid displacement during Lagrangian step
REAL(KIND=double)                              :: fractn       ! half the zone Lagrangian and Eulerian frid overlap
REAL(KIND=double)                              :: fractn2      ! 1 - 4/3*fractn

REAL(KIND=double)                              :: flux_z       ! total proton number advected

REAL(KIND=double)                              :: e_ph         ! photon energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_elec       ! electron energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_drip       ! drip energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_hvy        ! nuclear excited state energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_bind       ! nuclear binding energy (ergs g^{-1})
REAL(KIND=double)                              :: e_no_bind    ! internal energy - binding energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_total      ! internal energy (ergs g^{-1}) (not used here

REAL(KIND=double)                              :: third        ! 1/3
REAL(KIND=double)                              :: twothd       ! 2/3
REAL(KIND=double)                              :: fourthd      ! 4/3
REAL(KIND=double)                              :: fopi         ! 4*pi

REAL(KIND=double)                              :: x_sum        ! mass fraction sum
REAL(KIND=double)                              :: norm         ! normalization factor
REAL(KIND=double)                              :: ye_tot       ! total electron fraction

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in remap_comp_x')
 2001 FORMAT (' Deallocation problem for array ',a10,' in remap_comp_x')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rl          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r6          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dr(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr          '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (cmpl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmpl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cmp6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dcmp(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcmp        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (compl(nx+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'compl       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (comp6(nx+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp6       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dcomp(nx+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcomp       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dvol(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dm(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm          '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (delta(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (re(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 're          '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (cmp(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp         '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flux_nc_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_nc_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flux_ba_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_ba_rep '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flux_z_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_z_rep  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flux_b_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_b_rep  '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n_rep       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ba_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ba_rep      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n0_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n0_rep      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z0_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z0_rep      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_rep       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (b0_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b0_rep      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (b_rep(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b_rep       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (fluxr(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxr       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxcmp(nx+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxcmp     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (xn_t(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn_t        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

rl                         = zero
r6                         = zero
dr                         = zero

cmpl                       = zero
cmp6                       = zero
dcmp                       = zero
compl                      = zero
comp6                      = zero
dcomp                      = zero

dvol                       = zero
dvol0                      = zero

dm                         = zero
dm0                        = zero
delta                      = zero
re                         = zero

cmp                        = zero
flux_nc_rep                = zero
flux_ba_rep                = zero
flux_z_rep                 = zero
flux_b_rep                 = zero
n0_rep                     = zero
z_rep                      = zero
b_rep                      = zero
n_rep                      = zero

fluxr                      = zero
fluxcmp                    = zero
fluxbe                     = zero
fluxye_comp                = zero

xn_t                       = zero


IF ( first ) THEN
  first                    = .false.
  third                    = 1.d0/3.d0
  twothd                   = 2.d0 * third
  fourthd                  = 4.d0 * third
  fopi                     = 4.d0 * pi
  n_nucp1                  = nuc_number + 1
END IF

nmin                       = imin + 6
nmax                       = imax + 6
ntot                       = imax + 12

!-----------------------------------------------------------------------
!  Return if i_nnse = 0 (i.e., of all the zones are in NSE)
!-----------------------------------------------------------------------

IF ( i_nnse == 0 ) RETURN

!-----------------------------------------------------------------------
!  Specify nminc, the minimum padded index for matter in non-NSE
!-----------------------------------------------------------------------

nminc                      = i_nnse + 6
nmincm1                    = nminc - 1
 
!-----------------------------------------------------------------------
!  Calculate volumes before and after remap
!-----------------------------------------------------------------------

CALL volume_zone( ngeom, imin, imax, xa, dx, xa0, dx0, dvol, dvol0 )

!-----------------------------------------------------------------------
!  Compute the total binding energy for zones not in nse and the
!   outermost zone in nse. The composition of the outermost zone in nse
!   is determined in subroutine pre_remap_comp_x by deflashing it.
!  This is done in order to perform a consistent advection of
!   composition and binding energy.
!-----------------------------------------------------------------------

DO n = nmax, nmin, -1
  j                        = n - 5

  IF ( n >= nmincm1 ) THEN

    xn_t(1:n_nucp1)        = xn(j,1:n_nucp1)
    CALL eos_nnse_e( j, r(n), temp(n), ye(n), xn_t, nnc, a_nuc_rep(j), z_nuc_rep(j), &
&    be_nuc_rep(j), e_ph, e_elec, e_drip, e_hvy, e_bind, e_no_bind, e_total )

  END IF ! n >= nmincm1

!-----------------------------------------------------------------------
!  Set eb(n) = eb(nmincm1) ( = e_bind) for n < nminc
!-----------------------------------------------------------------------

  e_bind_zn0(n)            = e_bind * r(n) * dvol(n)
  eb(n)                    = e_bind

END DO ! nmax, nmin, -1

!-----------------------------------------------------------------------
!  Load boundary values of binding energy
!-----------------------------------------------------------------------

eb(nmin-4:nmin-1)          = eb(nmin)
eb(nmax+1:nmax+4)          = eb(nmax)
e_bind_zn0(nmin-4:nmin-1)  = e_bind_zn0(nmin)
e_bind_zn0(nmax+1:nmax+4)  = e_bind_zn0(nmax)

!-----------------------------------------------------------------------
!
!  Generate interpolation functions for the composition.
!
!
!        a(x) = a    + x[ da  = a   (1 - x) ]
!                L,j        j    6,j
!
!        x = ( xi - xi     )/dxi       xi     < xi < xi
!                     j-1/2     j        j-1/2         j+1/2
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cubic spline the density
!-----------------------------------------------------------------------

CALL parabola( nmin-3, nmax+3, ntot, zparax, r, dr, r6, rl, dm, 0, 0,   &
& ngeom )

!-----------------------------------------------------------------------
!  Cubic spline the composition
!   cmp from i_nnse-1 to nminc-1 have been padded with cmp(nminc-1)
!   cmp from nmax+1 to ntot = imax+12 have been padded with cmp(nmax)
!-----------------------------------------------------------------------

DO nc = 1,n_nucp1

  cmp(i_nnse-4:ntot)       = comp(i_nnse-4:ntot,nc)

  CALL parabola( nminc-4, nmax+3, imax+12, zparax, cmp, dcmp, cmp6, cmpl, &
&  dm, 0, 0, ngeom )

  dcomp(i_nnse-1:ntot,nc)  = dcmp(i_nnse-1:ntot)
  comp6(i_nnse-1:ntot,nc)  = cmp6(i_nnse-1:ntot)
  compl(i_nnse-1:ntot,nc)  = cmpl(i_nnse-1:ntot)

END DO ! nc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
!  Calculate the volume of the overlapping subshells (delta)
!-----------------------------------------------------------------------

IF ( ngeom == 0 ) THEN

  DO n = nmin, nmax+1
    delta(n)               = xa(n) - xa0(n)
  END DO

ELSE IF ( ngeom == 1 ) THEN

  DO n = nmin, nmax+1
    delta(n)               = xa(n) - xa0(n)
    delta(n)               = delta(n) * ( xa0(n) + .5d0 * delta(n) )
  END DO

ELSE IF ( ngeom == 2) THEN ! delta(n) = third * ( xa**3 - xa0**3 )

  DO n = nmin, nmax+1
    delta(n)               = xa(n) - xa0(n)
    delta(n)               = delta(n) * ( xa0(n) * ( xa0(n) + delta(n) )  + third * delta(n) * delta(n) )
  END DO

ELSE IF( ngeom == 3 ) THEN

  DO n = nmin, nmax+1
    delta(n)               = ( xa(n) - xa0(n) )  * radius
  END DO

ELSE IF( ngeom == 4 ) THEN

  DO n = nmin, nmax+1
    delta(n)               = ( cos(xa0(n)) - cos(xa(n)) ) * radius
  END DO

ELSE IF( ngeom == 5 ) THEN

  DO n = nmin, nmax+1
    delta(n)               = ( xa(n) - xa0(n) ) * radius 
  END DO

END IF ! ngeom

!-----------------------------------------------------------------------
!  Calculate the total mass (fluxr), and mass of each specie (fluxcmp)
!   in the subshell created by the overlap of the Lagrangian and
!   Eulerian grids.
!  If the zone face has moved to the left (deltx < 0), use the integral
!   from the left side of zone n (fluxrr).  If the zone face has moved
!   to the right (deltx > 0), use the integral from the right side of
!   zone k=n-1 (fluxrl).
!-----------------------------------------------------------------------

temp(nmax+1)               = temp(nmax)
ye(nmax+1)                 = ye(nmax)

DO n = nmin, nmax + 1
  deltx                    = xa(n) - xa0(n)

  IF ( deltx >= 0.0d0 ) THEN

!-----------------------------------------------------------------------
!
!        ...deltx > 0.0...
!
!
!        a   (x) = al    + z(da    + a6 (1-z)),    z = (x - x      )/dx
!         n-1        n-1       n-1     n-1                    n-1      n-1
!                                                                Lagr     Lagr
!
!             x
!              n
!           /   Lagr
!     1    |
!   _____  |          a   (x) dx = al    + da    - del z [ da    - a6    ( 1 -  4/3 del z ) ]
!   del x  |           n-1           n-1     n-1             n-1     n-1
!          /
!            x        - del x
!             n
!              Lagr
!
!........Material is remaped outwards
!
!                  |   n-1 --> n   |
!                  |     j --> j+1 |
!                  |               |
!                  |  fluxr(n) > 0 |
!         n-1      |       --------|-->         n
!      (j = n-6)   |               |       (j+1 = n-5)
!                  n               n
!                   Eulr            Lagr
!
!                  <---- del x ---->
!
!       n-1 material in overlapping shell added to n
!
!       n-1 = j is material on left; n = j+1 is the material on right
!
!-----------------------------------------------------------------------

    k                      = n - 1
    j                      = k - 5

    fractn                 = 0.5d0 * deltx/dx(k)
    fractn2                = 1.d0 - fourthd * fractn

!-----------------------------------------------------------------------
!  fluxr is the overlap mass
!-----------------------------------------------------------------------

    fluxr(n)               = delta(n) * ( rl(k) + dr(k) - fractn * ( dr(k) - fractn2 * r6(k) ) )

!-----------------------------------------------------------------------
!
!              \\\\\   nse(j) = 0 --> mse(j+1) = 0   /////
!              |||||   nse(j) = 0 --> mse(j+1) = 1   |||||
!              /////   nse(j) = 1 --> mse(j+1) = 0   \\\\\
!
!  If nse(j = n-6,i_ray) = 0 or 1, composition profile in zone j (zone n-1)
!   is integrated over causal domain to be transferred out of zone j.
!   If  nse(j+1) = 0 this composition is transferred into zone j+1 (zone n).
!   If nse(j+1) = 1, rho, e, and ye are transferred into zone j+1 in
!   subroutine remap_x.
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!  fluxcmp(n,nc) is the overlap mass of composition comp(k,nc) that will
!   be transferred from zone n-1 to zone n
!  fluxcmp_t is the overlap mass given by the sum of the composition 
!   masses before rescaling
!------------------------------------------------------------------------

    IF ( nse(j,ij_ray,ik_ray) /= 1  .or.  nse(j+1,ij_ray,ik_ray) /= 1 ) THEN

      fluxcmp(n,1:n_nucp1) = ( compl(k,1:n_nucp1) + dcomp(k,1:n_nucp1) &
&                          - fractn * ( dcomp(k,1:n_nucp1) - fractn2 * comp6(k,1:n_nucp1) ) ) * fluxr(n)
      fluxcmp_t            = SUM( fluxcmp(n,:) )

!-----------------------------------------------------------------------
!  Rescale factor phi modifyimg mass fractions to make them comsistent
!   with total mass advected if advection is from non-nse to non-nse
!-----------------------------------------------------------------------

      phi                  = fluxr(n)/( fluxcmp_t + epsilon )

!-----------------------------------------------------------------------
!  Rescale the overlap composition masses so that sum of mass fractions
!   is unity
!-----------------------------------------------------------------------

      fluxcmp(n,1:n_nucp1) = phi * fluxcmp(n,1:n_nucp1)

!-----------------------------------------------------------------------
!  Compute the total binding energy of the mass being transferred
!-----------------------------------------------------------------------

      xn_t(1:n_nucp1)      = fluxcmp(n,1:n_nucp1)/( fluxr(n) + epsilon )
      a_nuc_rep_t          = a_nuc_rep(j)
      z_nuc_rep_t          = z_nuc_rep(j)
      be_nuc_rep_t         = be_nuc_rep(j)

      CALL eos_nnse_e( j, r(n), temp(n), ye(n), xn_t, nnc, a_nuc_rep_t ,z_nuc_rep_t, &
&      be_nuc_rep_t, e_ph, e_elec, e_drip, e_hvy, e_bind, e_no_bind, e_total )

      fluxbe(n)            = e_bind * fluxr(n)

!-----------------------------------------------------------------------
!  Compute the total electron number being transferred
!-----------------------------------------------------------------------

      flux_z               = SUM( z_nuc(1:nuc_number) * fluxcmp(n,1:nuc_number)/a_nuc(1:nuc_number) )
      fluxye_comp(n)       = flux_z + z_nuc_rep_t * fluxcmp(n,n_nucp1)/a_nuc_rep_t

!------------------------------------------------------------------------
!  For representative heavy nucleus, mass number is advected so that
!   nuclei number, N , is computed by
!                   A
!
!               M X(A)
!         N  = _______
!          A    m  A
!                B
!
!        This requires that
!
!         M X(A)      M1 X(A1)     M2 X(A2)
!         ______   =  ________  +  ________
!           A            A1           A2
!
!        Charge number is advected to so that
!
!          M X(A)        M1 X(A1)       M2 X(A2)
!        Z ______   = Z1 ________  + Z2 ________
!            A              A1             A2
!
!        Likewise, binding energy is advected to so that
!
!           M X(A)         M1 X(A1)        M2 X(A2)
!        BE ______   = BE1 ________  + BE2 ________
!             A               A1              A2
!
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!  If a_nuc_rep(j) < 1, compl(k,n_nucp1) should be zero, but will be
!   interpolated as not quite zero. Set fluxcmp(n,n_nucp1) to zero to
!   avoid unphysical advection of a_nuc_rep, z_nuc_rep, and be_nuc_rep
!-----------------------------------------------------------------------

      IF ( a_nuc_rep(j) < 1.d0 ) fluxcmp(n,n_nucp1) = zero
      flux_nc_rep(n)       = fluxcmp(n,n_nucp1)/( a_nuc_rep(j) + epsilon )
      flux_ba_rep(n)       = fluxcmp(n,n_nucp1)
      flux_z_rep(n)        = z_nuc_rep(j)  * flux_nc_rep(n)
      flux_b_rep(n)        = be_nuc_rep(j) * flux_nc_rep(n)
  
    END IF ! nse(j,i_ray) /= 1  .or.  nse(j+1,i_ray) /= 1

  ELSE ! deltx < 0.0

!-----------------------------------------------------------------------
!
!        ...deltx < 0.0...
!
!        a (x) = al  + z(da  + a6 (1-z))    z = (x - x      )/dx
!         n        n       n     n                    n         n
!                                                      Lagr      Lagr
!
!             x      + del x
!              n
!           /   Lagr
!     1    |
!   _____  |          a (x) dx = al  + del z [ da  + a6  ( 1 -  4/3 del z ) ]
!   del x  |           n           n             n     n
!          |  
!          /
!            x
!             n
!              Lagr
!
!
!........Material is remaped inwards
!
!                  |    n-1 <-- n    |
!                  |    j-1 <-- j    |
!                  |                 |
!                  |  fluxr(n) <  0  |
!        n-1    <--|---------        |            n
!   (j-1 = n-6)    |                 |         (j = n-5)
!                  n                 n
!                   Lagr              Eulr
!
!
!                  <----- del x ----->
!
!         n material in overlapping shell added to n-1
!
!       n-1 = j-1 is material on left; n = j is the material on right
!
!-----------------------------------------------------------------------

    j                      = n - 5
  
    fractn                 = 0.5d0 * deltx/dx(n)
    fractn2                = 1.d0 + fourthd * fractn

!-----------------------------------------------------------------------
!  fluxr is the overlap mass
!-----------------------------------------------------------------------

    fluxr(n)               = delta(n) * ( rl(n) - fractn * ( dr(n) + fractn2 * r6(n) ) )

!-----------------------------------------------------------------------
!
!              \\\\\   nse(j) = 0 --> mse(j-1) = 0   /////
!              |||||   nse(j) = 0 --> mse(j-1) = 1   |||||
!              /////   nse(j) = 1 --> mse(j-1) = 0   \\\\\
!
!  If nse(j = n-5,i_ray) = 0 or 1, composition profile in zone j (zone n) 
!   is integrated over causal domain to be transferred out of zone j.
!  If  nse(j-1) = 0 the composition from j is transferred into zone j-1.
!  If nse(j-1) = 1, rho, e, and ye are transferred into zone j-1 in
!   subroutine remap_x.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  fluxcmp(n,nc) is the negative of the overlap mass of composition
!   comp(k,nc); it will be transferred from zone n to zone n-1
!  fluxcmp_t is the negative of the overlap mass of composition before
!   rescaling
!-----------------------------------------------------------------------

    IF ( nse(j,ij_ray,ik_ray) /= 1  .or.  nse(j-1,ij_ray,ik_ray) /= 1 ) THEN
    
      fluxcmp(n,1:n_nucp1) = ( compl(n,1:n_nucp1) &
&                          - fractn * ( dcomp(n,1:n_nucp1) + fractn2 * comp6(n,1:n_nucp1) ) ) * fluxr(n)
      fluxcmp_t            = SUM( fluxcmp(n,:) )

!-----------------------------------------------------------------------
!  Rescale factor phi modifyimg mass fractions to make them comsistent
!   with total mass advected if advection is from non-nse to non-nse
!-----------------------------------------------------------------------

      phi                  = fluxr(n)/( fluxcmp_t + epsilon )

!-----------------------------------------------------------------------
!  Rescale the overlap composition masses so that sum of mass fractions
!   is unity
!-----------------------------------------------------------------------

      fluxcmp(n,1:n_nucp1) = phi * fluxcmp(n,1:n_nucp1)

!-----------------------------------------------------------------------
!  Compute the total binding energy of the mass being transferred
!-----------------------------------------------------------------------

      xn_t(1:n_nucp1)      = fluxcmp(n,1:n_nucp1)/( fluxr(n) + epsilon )
      a_nuc_rep_t          = a_nuc_rep(j)
      z_nuc_rep_t          = z_nuc_rep(j)
      be_nuc_rep_t         = be_nuc_rep(j)

      CALL eos_nnse_e( j, r(n), temp(n), ye(n), xn_t, nnc, a_nuc_rep_t, &
&      z_nuc_rep_t, be_nuc_rep_t, e_ph, e_elec, e_drip, e_hvy, e_bind, &
&      e_no_bind, e_total )

      fluxbe(n)            = e_bind * fluxr(n)

!-----------------------------------------------------------------------
!  Compute the total electron fraction being transferred
!-----------------------------------------------------------------------

      flux_z               = SUM( z_nuc(1:nuc_number) * fluxcmp(n,1:nuc_number)/a_nuc(1:nuc_number) )
      fluxye_comp(n)       = flux_z + z_nuc_rep_t * fluxcmp(n,n_nucp1)/a_nuc_rep_t

!------------------------------------------------------------------------
!  For representative heavy nucleus, mass number is advected so that
!   nuclei number, N , is computed by
!                   A
!
!               M X(A)
!         N  = _______
!          A    m  A
!                B
!
!        This requires that
!
!         M X(A)      M1 X(A1)     M2 X(A2)
!         ______   =  ________  +  ________
!           A            A1           A2
!
!        Charge number is advected to so that
!
!          M X(A)        M1 X(A1)       M2 X(A2)
!        Z ______   = Z1 ________  + Z2 ________
!            A              A1             A2
!
!        Likewise, binding energy is advected to so that
!
!           M X(A)         M1 X(A1)        M2 X(A2)
!        BE ______   = BE1 ________  + BE2 ________
!             A               A1              A2
!
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!  If a_nuc_rep(j) < 1, compl(k,n_nucp1) should be zero, but will be
!   interpoated as not quite zero. Set fluxcmp(n,n_nucp1) to zero to
!   avoid unphysical advection of a_nuc_rep, z_nuc_rep, and be_nuc_rep
!-----------------------------------------------------------------------

      IF ( a_nuc_rep(j) < 1.d0 ) fluxcmp(n,n_nucp1) = zero
      flux_nc_rep(n)       = fluxcmp(n,n_nucp1)/( a_nuc_rep(j) + epsilon )
      flux_ba_rep(n)       = fluxcmp(n,n_nucp1)
      flux_z_rep(n)        = z_nuc_rep(j)  * flux_nc_rep(n)
      flux_b_rep(n)        = be_nuc_rep(j) * flux_nc_rep(n)

    END IF ! nse(j,i_ray) == 0  .or.  nse(j-1,i_ray) == 0

  END IF ! deltx >= 0.0

END DO ! n = nmin, nmax + 1

!-----------------------------------------------------------------------
!
!                   \\\\\ DO THE ADVECTION /////
!
!-----------------------------------------------------------------------

DO n = nmin, nmax

!-----------------------------------------------------------------------
!  Determine the masses and densities in the remapped grid.
!-----------------------------------------------------------------------

  dm(n)                    = r(n) * dvol(n)
  dm0(n)                   = dm(n) + fluxr(n) - fluxr(n+1)
  re(n)                    = dm0(n)/dvol0(n)
  re(n)                    = DMAX1( smallr, re(n) )
  dm0(n)                   = 1.0d0/( dm0(n) )

  IF ( n >= nminc ) THEN
  
    j                      = n - 5

!-----------------------------------------------------------------------
!  Determine the mass number, charge number, and binding energy of the
!   auxiliary heavy nucleus in the remapped grid.
!-----------------------------------------------------------------------

    n_rep(n)               = comp(n,n_nucp1) * dm(n)/( a_nuc_rep(j) + epsilon )
    n0_rep(n)              = n_rep(n)                 + flux_nc_rep(n) - flux_nc_rep(n+1)
    ba_rep(n)              = comp(n,n_nucp1) * dm(n)  + flux_ba_rep(n) - flux_ba_rep(n+1)
    z_rep(n)               = z_nuc_rep(j)  * n_rep(n) + flux_z_rep(n)  - flux_z_rep(n+1)
    b_rep(n)               = be_nuc_rep(j) * n_rep(n) + flux_b_rep(n)  - flux_b_rep(n+1)
    IF ( comp(n,n_nucp1) > epsilon ) THEN
      a_nuc_rep  (j)       = ba_rep(n)/( n0_rep(n) + epsilon )
      a_nuc_rep_e(j)       = a_nuc_rep(j)
      z_nuc_rep (j)        = z_rep(n)       /( n0_rep(n) + epsilon )
      z_nuc_rep_e(j)       = z_nuc_rep (j)
      be_nuc_rep(j)        = b_rep(n)       /( n0_rep(n) + epsilon )
      be_nuc_rep_e(j)      = be_nuc_rep(j)
    ELSE
      a_nuc_rep_e(j)       = a_nuc_rep(j)
      z_nuc_rep_e(j)       = z_nuc_rep(j)
      be_nuc_rep_e(j)      = be_nuc_rep(j)
    END IF      

!-----------------------------------------------------------------------
!  Determine the composition mass fractions in the remapped grid.
!-----------------------------------------------------------------------

    comp(n,1:n_nucp1)      = ( comp(n,1:n_nucp1) * dm(n) + fluxcmp(n,1:n_nucp1) - fluxcmp(n+1,1:n_nucp1) ) * dm0(n)

!-----------------------------------------------------------------------
!  Normalize composition mass fractions
!-----------------------------------------------------------------------

    x_sum                  = SUM( comp(n,1:n_nucp1) )
    norm                   = 1.d0/x_sum
    comp(n,1:n_nucp1)      = norm * comp(n,1:n_nucp1)

  END IF !  n >= nminc

END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
!
!                \\\\\ RESTORE COMPOSITION TO XN /////
!
!-----------------------------------------------------------------------

DO i = i_nnse,imax
  n = i + 6
  xn  (i+1,:)              = comp(n,:)
  xn_e(i+1,:)              = comp(n,:)
END DO ! i = i_nnse,imax

DO n = nminc,nmax
  j                        = n - 5
  ye_tot                   = SUM( xn(j,1:nuc_number) * z_nuc(1:nuc_number)/a_nuc(1:nuc_number) )
  ye_tot                   = ye_tot + xn(j,n_nucp1) * z_nuc_rep(j)/a_nuc_rep(j)
  ye(n)                    = ye_tot
  ye_r(n)                  = ye(n)
END DO ! n = nminc,nmax

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------


DEALLOCATE (rl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rl          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (r6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r6          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr          '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (cmpl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmpl        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (cmp6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dcmp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcmp        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (compl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'compl       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (comp6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp6       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dcomp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcomp       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dvol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dm, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm          '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dm0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (delta, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (re, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 're          '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (cmp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp         '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (flux_nc_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_nc_rep '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (flux_ba_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_ba_rep '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (flux_z_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_z_rep  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (flux_b_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_b_rep  '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (n_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n_rep      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ba_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ba_rep      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (n0_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n0_rep       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (z0_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z0_rep      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (z_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_rep       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (b0_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b0_rep      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (b_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b_rep       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (fluxr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxr       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxcmp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxcmp     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (xn_t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn_t        '; WRITE (nlog,2001) var_name; END IF


RETURN
END SUBROUTINE remap_comp_x















