SUBROUTINE remap_comp_y( ngeom, jmin, jmax, ji_ray, jk_ray, ny, nnc )
!-----------------------------------------------------------------------
!
!    File:         remap_comp_y
!    Module:       remap_comp_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/03/05
!
!    Purpose:
!      To remap the composition and binding energy along angular rays
!       using piecewise parabolic functions. No flattening is used on remap
!       (pass dummy array: dum to parabola.f).
!
!    Input arguments:
!  ngeom       : geometry index
!  jmin        : minimum y-array index
!  jmax        : maximim y-array index
!  ji_ray      : x (radial) index of a specific y (angular) ray
!  jk_ray      : z (azimuthal) index of a specific y (angular) ray
!  ny          : y_array extent
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
!  edit_module, eos_snc_y_module, evh1_global, evh1_sweep,
!  evh1_zone, mgfld_remap_module
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
USE numerical_module, ONLY : zero, epsilon, third, frpith
USE physcnst_module, ONLY: pi

USE edit_module, ONLY : nlog
USE eos_snc_y_module, ONLY: nse, xn, a_nuc_rep, z_nuc_rep, be_nuc_rep, &
& nuc_number, a_nuc, z_nuc, be_nuc
USE evh1_global, ONLY: smallr
USE evh1_sweep, ONLY: radius, ye
USE evh1_zone, ONLY : zparay
USE mgfld_remap_module, ONLY : comp, r, temp, ye_r=>ye, xa, dx, xa0, dx0, &
& fluxbe, e_bind_zn0, eb, fluxye_comp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                            :: ngeom        ! geometry index
INTEGER, INTENT(in)                            :: jmin         ! minimum y-array index
INTEGER, INTENT(in)                            :: jmax         ! maximim y-array index
INTEGER, INTENT(in)                            :: ji_ray       ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)                            :: jk_ray       ! z (azimuthal) index of a specific y (angular) ray
INTEGER, INTENT(in)                            :: ny           ! y-array extent
INTEGER, INTENT(in)                            :: nnc          ! composition array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                             :: var_name

LOGICAL                                        :: first = .true.

INTEGER                                        :: istat        ! allocation status flag
INTEGER                                        :: nc           ! composition index
INTEGER                                        :: n            ! padded aray index
INTEGER                                        :: k            ! k-1 padded aray index
INTEGER                                        :: j            ! y-zone index
INTEGER                                        :: i            ! unpadded aray index
INTEGER                                        :: n_nucp1      ! nuc_number + 1
INTEGER                                        :: nmin         ! minimum padded index
INTEGER                                        :: nmax         ! maximum padded index
INTEGER                                        :: ntot         ! number of zones (real plus ghost)
INTEGER                                        :: nse_zn       ! 1 if all zones or in nse; otherwise 0

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
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dum          ! dummy array

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: cmp          ! 1-D storage array for the composition array
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: flux_n_rep   ! proportional to heavy nuclei number advected
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: flux_z_rep   ! proportional to heavy nuclei charge advected
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: flux_b_rep   ! proportional to binding energy advected
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: n0_rep       ! proportional to initial number of rep heavy nuclei/zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: n_rep        ! proportional to number of rep heavy nuclei/zone
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

REAL(KIND=double)                              :: flux_z       ! total charge number advected
REAL(KIND=double)                              :: flux_a       ! total mass number advected

REAL(KIND=double)                              :: e_ph         ! photon energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_elec       ! electron energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_drip       ! drip energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_hvy        ! nuclear excited state energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_bind       ! nuclear binding energy (ergs g^{-1})
REAL(KIND=double)                              :: e_no_bind    ! internal energy - binding energy (ergs g^{-1}) (not used here)
REAL(KIND=double)                              :: e_total      ! internal energy (ergs g^{-1}) (not used here

REAL(KIND=double)                              :: x_sum        ! mass fraction sum
REAL(KIND=double)                              :: norm         ! normalization factor
REAL(KIND=double)                              :: ye_tot       ! total electron fraction

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in remap_comp_y')
 2001 FORMAT (' Deallocation problem for array ',a10,' in remap_comp_y')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Determine nse_zn, and return if nse_zn = 1, i.e., all zone in NSE
!-----------------------------------------------------------------------

nse_zn                 = 1
DO j = jmin,jmax
  IF ( nse(j,ji_ray,jk_ray) == 0 ) THEN
    nse_zn             = 0
    EXIT
  END IF
END DO

IF ( nse_zn == 1 ) RETURN

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (rl(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r6(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dr(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (cmpl(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmpl      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (cmp6(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp6      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dcmp(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcmp      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (compl(ny+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'compl     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (comp6(ny+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp6     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dcomp(ny+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcomp     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dvol(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dvol0(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dm(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm0(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (delta(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dum(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dum       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (cmp(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flux_n_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_n_rep'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flux_z_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_z_rep'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (flux_b_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_b_rep'; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n0_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n0_rep    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (n_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n_rep     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z0_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z0_rep    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (z_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_rep     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (b0_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b0_rep    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (b_rep(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b_rep     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (fluxr(ny+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxr     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxcmp(ny+12,nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxcmp   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (xn_t(nnc), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn_t      '; WRITE (nlog,1001) var_name; END IF

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
dum                        = zero

cmp                        = zero
flux_n_rep                 = zero
flux_z_rep                 = zero
flux_b_rep                 = zero
n_rep                      = zero
z_rep                      = zero
b_rep                      = zero

fluxr                      = zero
fluxcmp                    = zero
fluxbe                     = zero
fluxye_comp                = zero

IF ( first ) THEN
  first                    = .false.
  n_nucp1                  = nuc_number + 1
END IF

nmin                       = jmin + 6
nmax                       = jmax + 6
ntot                       = jmax + 12

!-----------------------------------------------------------------------
!  Calculate volumes before and after remap
!-----------------------------------------------------------------------

CALL volume_zone( ngeom, jmin, jmax, xa, dx, xa0, dx0, dvol, dvol0 )

!-----------------------------------------------------------------------
!  Compute the total binding energy for zones not in nse.
!  This is done in order to perform a consistent advection of
!   composition and binding energy.
!-----------------------------------------------------------------------

DO n = nmin, nmax
  j                        = n - 6
 
  xn_t(1:n_nucp1)          = xn(j,1:n_nucp1)

  CALL eos_nnse_e( j, r(n), temp(n), ye(n), xn_t, nnc, a_nuc_rep(j), z_nuc_rep(j), &
&  be_nuc_rep(j), e_ph, e_elec, e_drip, e_hvy, e_bind, e_no_bind, e_total )

  e_bind_zn0(n)            = e_bind * r(n) * dvol(n)
  eb(n)                    = e_bind

END DO ! n = nmin, nmax

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

CALL parabola( nmin-3, nmax+3, ntot, zparay, r, dr, r6, rl, dum, 0, 0,  &
& ngeom )

!-----------------------------------------------------------------------
!  Cubic spline the composition
!-----------------------------------------------------------------------

DO nc = 1,n_nucp1

  cmp(jmin:ntot)           = comp(jmin:ntot,nc)

  CALL parabola( nmin-3, nmax+3, ntot, zparay, cmp, dcmp, cmp6, cmpl,   &
&  dum, 0, 0, ngeom )
  dcomp(jmin:ntot,nc)      = dcmp(jmin:ntot)
  comp6(jmin:ntot,nc)      = cmp6(jmin:ntot)
  compl(jmin:ntot,nc)      = cmpl(jmin:ntot)

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
!  If the zone face has moved to the left (deltx < 0), USE the integral
!   from the left side of zone n (fluxrr).  If the zone face has moved
!   to the right (deltx > 0), USE the integral from the right side of
!   zone k=n-1 (fluxrl).
!-----------------------------------------------------------------------

DO n = nmin+1, nmax
  deltx                    = xa(n) - xa0(n)

  IF ( deltx >= 0.0 ) THEN

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
!      (j = n-7)   |               |       (j+1 = n-6)
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
    j                      = k - 6

    fractn                 = 0.5d0 * deltx/dx(k)
    fractn2                = 1.d0 - frpith * fractn

!-----------------------------------------------------------------------
!  fluxr is the overlap mass
!-----------------------------------------------------------------------

    fluxr(n)               = delta(n) * ( rl(k) + dr(k) - fractn * ( dr(k) - fractn2 * r6(k) ) )

!-----------------------------------------------------------------------
!
!        \\\\\ nse(j) = 0 --> mse(j+1) = 0 or nse(j+1) = 1 /////
!        /////        nse(j) = 1 --> mse(j+1) = 0          \\\\\
!
!  If nse(j = n-7,ji_ray,jk_ray) = 0 or 1, composition profile in zone j
!   (zone n-1) is integrated over causal domain to be transferred out of
!   zone j.
!   If  nse(j+1) = 0 this composition is transferred into zone j+1 (zone n).
!   If nse(j+1) = 1, rho, e, and ye are transferred into zone j+1 in
!   subroutine remap_y.
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!  fluxcmp(n,nc) is the overlap mass of composition comp(k,nc) that will
!   be transferred from zone n-1 to zone n
!  fluxcmp_t is the overlap mass of composition before rescaling
!------------------------------------------------------------------------

    IF ( nse(j,ji_ray,jk_ray) /= 1  .or.  nse(j+1,ji_ray,jk_ray) /= 1 ) THEN

      fluxcmp(n,1:n_nucp1) = ( compl(k,1:n_nucp1) + dcomp(k,1:n_nucp1) &
&                          - fractn * ( dcomp(k,1:n_nucp1) - fractn2 * comp6(k,1:n_nucp1) ) ) * fluxr(n)
      fluxcmp_t            = SUM( fluxcmp(n,:) )

!-----------------------------------------------------------------------
!  Rescale factor phi modifyimg mass fractions to make them comsistent
!   with total mass advected if advection is from non-nse to non-nse
!-----------------------------------------------------------------------

      phi                  = fluxr(n)/( fluxcmp_t + epsilon )

!-----------------------------------------------------------------------
!  Rescale mass fractions
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

      flux_z               = SUM( z_nuc(1:nuc_number) * xn_t(1:nuc_number) )
      flux_a               = SUM( a_nuc(1:nuc_number) * xn_t(1:nuc_number) )
      flux_z               = flux_z + z_nuc_rep_t * xn_t(n_nucp1)
      flux_a               = flux_a + a_nuc_rep_t * xn_t(n_nucp1)
      fluxye_comp(n)       = fluxr(n) * flux_z/( flux_a + 1.d-20 )
 
    END IF ! nse(j,ji_ray,jk_ray) /= 1  .or.  nse(j+1,ji_ray,jk_ray) /= 1

!------------------------------------------------------------------------
!        Mass number is advected so that nuclei number, N , is cmputed by
!                                                        A
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
    flux_n_rep(n)          = fluxcmp(n,n_nucp1)/( a_nuc_rep(j) + epsilon )
    flux_z_rep(n)          = z_nuc_rep(j)  * flux_n_rep(n)
    flux_b_rep(n)          = be_nuc_rep(j) * flux_n_rep(n)

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
!   (j-1 = n-7)    |                 |         (j = n-6)
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

    j                      = n - 6
 
    fractn                 = 0.5d0 * deltx/dx(n)
    fractn2                = 1.d0 + frpith * fractn

!-----------------------------------------------------------------------
!  fluxr is the overlap mass
!-----------------------------------------------------------------------

    fluxr(n)               = delta(n) * ( rl(n) - fractn * ( dr(n) + fractn2 * r6(n) ) )

!-----------------------------------------------------------------------
!
!        \\\\\ nse(j) = 0 --> mse(j-1) = 0 or nse(j-1) = 1 /////
!        /////        nse(j) = 1 --> mse(j-1) = 0          \\\\\
!
!  If nse(j = n-6,ji_ray,jk_ray) = 0 or 1, composition profile in zone j
!   (zone n) is integrated over causal domain to be transferred out of
!   zone j.
!  If  nse(j-1) = 0 the composition from j is transferred into zone j-1.
!  If nse(j-1) = 1, rho, e, and ye are transferred into zone j-1 in
!   subroutine remap_y.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  fluxcmp(n,nc) is the negative of the overlap mass of composition
!   comp(k,nc); it will be transferred from zone n to zone n-1
!  fluxcmp_t is the negative of the overlap mass of composition before
!   rescaling
!-----------------------------------------------------------------------

    IF ( nse(j,ji_ray,jk_ray) /= 1  .or.  nse(j-1,ji_ray,jk_ray) /= 1 ) THEN
    
      fluxcmp(n,1:n_nucp1) = ( compl(n,1:n_nucp1) &
&                          - fractn * ( dcomp(n,1:n_nucp1) + fractn2 * comp6(n,1:n_nucp1) ) ) * fluxr(n)
      fluxcmp_t            = SUM( fluxcmp(n,:) )

!-----------------------------------------------------------------------
!  Rescale factor phi modifyimg mass fractions to make them comsistent
!   with total mass advected if advection is from non-nse to non-nse
!-----------------------------------------------------------------------

      phi                  = fluxr(n)/( fluxcmp_t + epsilon )

!-----------------------------------------------------------------------
!  Rescale mass fractions
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

      flux_z               = SUM( z_nuc(1:nuc_number) * xn_t(1:nuc_number) )
      flux_a               = SUM( a_nuc(1:nuc_number) * xn_t(1:nuc_number) )
      flux_z               = flux_z + z_nuc_rep_t * xn_t(n_nucp1)
      flux_a               = flux_a + a_nuc_rep_t * xn_t(n_nucp1)
      fluxye_comp(n)       = fluxr(n) * flux_z/( flux_a + 1.d-20 )

    END IF ! nse(j-1,ji_ray,jk_ray) == 0

!------------------------------------------------------------------------
!  Mass number is advected so that nuclei number, N , is computed by
!                                                  A
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
    flux_n_rep(n)          = fluxcmp(n,n_nucp1)/( a_nuc_rep(j) + epsilon )
    flux_z_rep(n)          = z_nuc_rep(j)  * flux_n_rep(n)
    flux_b_rep(n)          = be_nuc_rep(j) * flux_n_rep(n)

  END IF ! deltx >= 0.0

END DO ! n = nmin+1, nmax

!-----------------------------------------------------------------------
!
!                   \\\\\ DO THE ADVECTION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Advect mass and composition mass fractions by moving the subshell
!   quantities into the appropriate Eulerian zone. 
!-----------------------------------------------------------------------

DO n = nmin, nmax

  dm(n)                    = r(n) * dvol(n)
  dm0(n)                   = ( dm(n) + fluxr(n) - fluxr(n+1) )
  dm0(n)                   = 1.0d0/( dm0(n) )
  
  j                        = n - 6
 
  n0_rep(n)                = comp(n,n_nucp1) * dm(n)/( a_nuc_rep(j) + epsilon )
  n_rep(n)                 = n0_rep(n)                 + flux_n_rep(n) - flux_n_rep(n+1)
  z_rep(n)                 = z_nuc_rep(j)  * n0_rep(n) + flux_z_rep(n) - flux_z_rep(n+1)
  b_rep(n)                 = be_nuc_rep(j) * n0_rep(n) + flux_b_rep(n) - flux_b_rep(n+1)

  comp(n,1:n_nucp1)        = ( comp(n,1:n_nucp1) * dm(n) + fluxcmp(n,1:n_nucp1) - fluxcmp(n+1,1:n_nucp1) ) * dm0(n)

  IF ( comp(n,n_nucp1) > epsilon ) THEN
    a_nuc_rep  (j)         = comp(n,n_nucp1)/( n_rep(n) * dm0(n) + epsilon )
    z_nuc_rep (j)          = z_rep(n)       /( n_rep(n) + epsilon )
    be_nuc_rep(j)          = b_rep(n)       /( n_rep(n) + epsilon )
  END IF

!-----------------------------------------------------------------------
!  Normalize composition mass fractions
!-----------------------------------------------------------------------

  x_sum                    = SUM( comp(n,1:n_nucp1) )
  norm                     = 1.d0/x_sum
  comp(n,1:n_nucp1)        = norm * comp(n,1:n_nucp1)

END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
!
!                \\\\\ RESTORE COMPOSITION TO XN /////
!
!-----------------------------------------------------------------------

DO nc = 1,n_nucp1
  DO i = jmin,jmax
    n = i + 6
    xn(i,nc)               = DMAX1( comp(n,nc), zero )
  END DO
END DO

DO n = nmin, nmax
  j                        = n - 6
  ye_tot                   = SUM( xn(j,1:nuc_number) * z_nuc(1:nuc_number)/( a_nuc(1:nuc_number) + epsilon ) )
  ye_tot                   = ye_tot + xn(j,n_nucp1) * z_nuc_rep(j)/( a_nuc_rep(j) + epsilon )
  ye(n)                    = ye_tot
  ye_r(n)                  = ye(n)
END DO ! n = nminc,nmax

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (rl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rl        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (r6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (cmpl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmpl      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (cmp6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp6      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dcmp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcmp      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (compl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'compl     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (comp6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'comp6     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dcomp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dcomp     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dvol, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dvol0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dvol0     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dm, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dm0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (delta, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dum, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dum       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (cmp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'cmp       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (flux_n_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_n_rep'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (flux_z_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_z_rep'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (flux_b_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_b_rep'; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (n0_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n0_rep    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (n_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'n_rep     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (z0_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z0_rep    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (z_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'z_rep     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (b0_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b0_rep    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (b_rep, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'b_rep     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (fluxr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxr     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxcmp, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxcmp   '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (xn_t, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'xn_t      '; WRITE (nlog,2001) var_name; END IF


RETURN
END SUBROUTINE remap_comp_y















