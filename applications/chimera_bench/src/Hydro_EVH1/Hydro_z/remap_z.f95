SUBROUTINE remap_z( ngeom, kmin, kmax, ki_ray, kj_ray, nz, nez, nnu, &
& nnc, j_shock )
!-----------------------------------------------------------------------
!
!    File:         remap_z
!    Module:       remap_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/23/07
!
!    Purpose:
!      To remap the mass, momentum, and energy and other state variables
!       along a radial ray from the updated lagrangian grid to the fixed
!       Eulerian grid, using piecewise parabolic functions.
!
!    Input arguments:
!
!  ngeom       : geometry index
!  kmin        : minimum z-array index
!  kmax        : maximim z-array index
!  ki_ray      : x (radial) index of a specific z (azimuthal) ray
!  kj_ray      : y (angular) index of a specific z (azimuthal) ray
!  nz          : z-array extent
!  nez         : neutrino energy array extent
!  nnu         : neutrino flavor array extent
!  nnc         : composition array extent
!  j_shock     : zones marked for additional diffusion
!
!    Output arguments:
!        none
!
!    Subprograms called:
!
!  e_compose   : computes the specific total energy
!  parabola    : computes piecewise parabolic fits to the profiles of variables to be remapped
!  sweepbc     : fills ghost zones with appropriate boundary conditions
!  paraset     : computes parabolic coefficients for the Eulerian grid
!  e_decompose : extracts the internal energy after the remap of the total energy
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, eos_snc_z_module, evh1_global, evh1_sweep,
!  evh1_zone, mgfld_remap_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY: zero
USE physcnst_module, ONLY: pi, rmu, ergfoe

USE edit_module, ONLY : nlog
USE eos_snc_z_module, ONLY: nse, xn, nuc_number, a_nuc, z_nuc, a_nuc_rep, z_nuc_rep
USE evh1_global, ONLY: smallr, nleftz, nrightz, v_diff
USE evh1_sweep, ONLY: xa, dx, r, u, v, w, ei, e, ye, xa0, dx0, &
& radius, dvol, dvol0, se=>entrop
USE evh1_zone, ONLY : zparaz
USE mgfld_remap_module, ONLY : fluxbe, e_bind_zn0, eb, fluxye_comp

USE parallel_module, ONLY : myid

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                          :: ngeom      ! geometry index
INTEGER, INTENT(in)                          :: kmin       ! minimum z-array index
INTEGER, INTENT(in)                          :: kmax       ! maximim z-array index
INTEGER, INTENT(in)                          :: ki_ray     ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                          :: kj_ray     ! y (angular) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)                          :: nz         ! z-array extent
INTEGER, INTENT(in)                          :: nez        ! neutrino energy array extent
INTEGER, INTENT(in)                          :: nnu        ! neutrino energy flavor extent
INTEGER, INTENT(in)                          :: nnc        ! composition array extent
INTEGER, INTENT(in), DIMENSION(nz)           :: j_shock    ! zones marked for added z-diffusion

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                           :: var_name

LOGICAL                                      :: first = .true.

INTEGER                                      :: istat      ! allocation status flag
INTEGER                                      :: nmin       ! minimum padded index
INTEGER                                      :: nmax       ! maximum padded index
INTEGER                                      :: ntot       ! total number of padded indices
INTEGER                                      :: n          ! padded aray index
INTEGER                                      :: k          ! k-1 padded aray index
INTEGER                                      :: j          ! angular zone index
INTEGER                                      :: nc         ! abundance index
INTEGER                                      :: n_nucp1    ! nuc_number + 1
INTEGER, PARAMETER                           :: n_ei = 0   ! number of inner zones to remap internal energy

REAL(KIND=double)                            :: ur         ! x-velocity on right-hand side of zone
REAL(KIND=double)                            :: vr         ! y-velocity on right-hand side of zone
REAL(KIND=double)                            :: wr         ! z-velocity on right-hand side of zone
REAL(KIND=double)                            :: eir        ! inernal energy on right-hand side of zone
REAL(KIND=double)                            :: er         ! total energy on right-hand side of zone

REAL(KIND=double)                            :: deltx      ! grid displacement during Lagrangian step
REAL(KIND=double)                            :: fractn     ! half the zone Lagrangian and Eulerian frid overlap
REAL(KIND=double)                            :: fractn2    ! 1 - 4/3*fractn

REAL(KIND=double)                            :: third      ! 1/3
REAL(KIND=double)                            :: twothd     ! 2/3
REAL(KIND=double)                            :: fourthd    ! 4/3
REAL(KIND=double)                            :: fopi       ! 4*pi

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ul         ! x-velocity at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: u6         ! x-velocity parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: du         ! x-velocity slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: vl         ! y-velocity at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: v6         ! y-velocity parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dv         ! y-velocity slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: wl         ! z-velocity at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: w6         ! z-velocity parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dw         ! z-velocity slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: rl         ! density at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: r6         ! density parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dr         ! density slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: el         ! energy at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: e6         ! energy parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: de         ! energy slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ei_b       ! internal energy minus the nuclear binding energy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: eil        ! internal energy at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ei6        ! internal energy parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dei        ! internal energy slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: yel        ! electron fraction at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ye6        ! electron fraction energy parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dye        ! electron fraction energy slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: sl         ! entropy at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: s6         ! entropy energy parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ds         ! entropy energy slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dm         ! mass after lagrangian step
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dm0        ! mass after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: delta      ! volume of overlapping subshells

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxr      ! mass contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxu      ! x-momentum contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxv      ! y-momentum contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxw      ! z-momentum contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxe      ! energy contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxei     ! internal energy contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxye     ! electrons contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxs      ! entropy contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: var1_i     ! initial value of a variable to numerically diffuse
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: var2_i     ! initial value of a variable to numerically diffuse
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: var3_i     ! initial value of a variable to numerically diffuse

REAL(KIND=double)                            :: z_tot
REAL(KIND=double)                            :: a_tot

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: es         ! remapped interal energy

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in remap_z')
 2001 FORMAT (' Deallocation problem for array ',a10,' in remap_z')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (ul(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ul        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (du(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'du        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (vl(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dv(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dv        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (wl(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dw(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dw        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rl(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dr(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (el(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'el        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (ei_b(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_b      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (eil(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eil       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ei6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei6       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dei(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dei       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (yel(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yel       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ye6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye6       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dye(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (sl(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (s6(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 's6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ds(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ds        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dm(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm0(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (delta(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (fluxr(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxr     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxu(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxu     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxv(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxv     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxw(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxw     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxe(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxe     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxei(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxei    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxye(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxye    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxs(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxs     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (var1_i(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'var1_i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (var2_i(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'var2_i    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (var3_i(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'var3_i    '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (es(nz+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'es        '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

ul               = zero
u6               = zero
du               = zero

vl               = zero
v6               = zero
dv               = zero

wl               = zero
w6               = zero
dr               = zero

rl               = zero
r6               = zero
ul               = zero

el               = zero
e6               = zero
de               = zero

ei_b             = zero
eil              = zero
ei6              = zero
dei              = zero

yel              = zero
ye6              = zero
dye              = zero

sl               = zero
s6               = zero
ds               = zero

dm               = zero
dm0              = zero
delta            = zero

fluxr            = zero
fluxu            = zero
fluxv            = zero
fluxw            = zero
fluxe            = zero
fluxei           = zero
fluxye           = zero
fluxs            = zero

es               = zero

IF ( first ) THEN
  first          = .false.
  third          = 1.d0/3.d0
  twothd         = 2.d0 * third
  fourthd        = 4.d0 * third
  fopi           = 4.d0 * pi
  n_nucp1        = nuc_number + 1
END IF

nmin             = kmin + 6
nmax             = kmax + 6
ntot             = kmax + 12

!-----------------------------------------------------------------------
!  Subtract the nuclear binding energy from the internal energy if the
!   matter is not in nse.
!  If the matter is in nse, subtract the binding energy of the first
!   non-nse zone to ensure continuity at the nse - non-nse boundary.
!-----------------------------------------------------------------------

CALL e_compose( nmin, nmax )

DO n = nmax+3, nmin-3, -1
  ei_b(n)        = ei(n) - eb(n)
  e(n)           = e(n)  - eb(n)
END DO ! n = nmax+3, nmin-3, -1

!-----------------------------------------------------------------------
!
!  Generate interpolation functions, saving da, al for constructing 
!   left and right total energy states. (dmu is passed as dummy array)
!
!
!        a(x) = a    + x[ da  = a   (1 - x) ]
!                L,j        j    6,j
!
!        x = ( xi - xi     )/dxi       xi     < xi < xi
!                     j-1/2     j        j-1/2         j+1/2
!
!-----------------------------------------------------------------------

CALL parabola( nmin-3, nmax+3, ntot, zparaz, r,    dr,  r6,  rl,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparaz, u,    du,  u6,  ul,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparaz, v,    dv,  v6,  vl,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparaz, w,    dw,  w6,  wl,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparaz, ei_b, dei, ei6, eil, dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparaz, ye,   dye, ye6, yel, dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparaz, se,   ds,  s6,  sl,  dm, 0, &
& 0, ngeom )

!-----------------------------------------------------------------------
!  Use the profiles for density, pressure, and velocities to calculate
!   consistent values of the left and right values of total energy
!-----------------------------------------------------------------------

DO n = nmin-1, nmax+1
  ur             = ul(n)  + du(n)
  vr             = vl(n)  + dv(n)
  wr             = wl(n)  + dw(n)
  eir            = eil(n) + dei(n)
  el(n)          = eil(n) + 0.5d0 * ( ul(n)**2 + vl(n)**2 + wl(n)**2 )
  er             = eir    + 0.5d0 * ( ur**2    + vr**2    + wr**2    )
  de(n)          = er - el(n)
END DO ! n = nmin-1, nmax+1

!-----------------------------------------------------------------------
!  Use the constructed left and right states to get parabolas for total
!   energy
!-----------------------------------------------------------------------

CALL parabola( nmin-3, nmax+3, ntot, zparaz, e, de, e6, el, dm, 0, 1,   &
& ngeom )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!-----------------------------------------------------------------------
!  Calculate the volume of the overlapping subshells (delta).
!-----------------------------------------------------------------------

IF ( ngeom == 0 ) THEN

  DO n = nmin, nmax+1
    delta(n)     = xa(n) - xa0(n)
  END DO

ELSE IF ( ngeom == 1 ) THEN

  DO n = nmin, nmax+1
    delta(n)     = xa(n) - xa0(n)
    delta(n)     = delta(n) * ( xa0(n) + .5d0 * delta(n) )
  END DO

ELSE IF ( ngeom == 2) THEN ! delta(n) = third * ( xa**3 - xa0**3 )

  DO n = nmin, nmax+1
    delta(n)     = xa(n) - xa0(n)
    delta(n)     = delta(n) * ( xa0(n) * ( xa0(n) + delta(n))  + third * delta(n)**2 )
  END DO

ELSE IF( ngeom == 3 ) THEN

  DO n = nmin, nmax+1
    delta(n)     = ( xa(n) - xa0(n) )  * radius
  END DO

ELSE IF( ngeom == 4 ) THEN

  DO n = nmin, nmax+1
    delta(n)     = ( cos(xa0(n)) - cos(xa(n)) ) * radius
  END DO

ELSE IF( ngeom == 5 ) THEN

  DO n = nmin, nmax+1
    delta(n)     = ( xa(n) - xa0(n) ) * radius 
  END DO

END IF ! ngeom

!-----------------------------------------------------------------------
!  Calculate the total mass (fluxr), momentum (fluxu), and energy (fluxe)
!   in the subshell created by the overlap of the Lagrangian and Eulerian
!   grids.
!  If the zone face has moved to the left (deltx < 0), use the integral
!   from the left side of zone n (fluxrr).  If the zone face has moved
!   to the right (deltx > 0), USE the integral from the right side of
!   zone k=n-1 (fluxrl).
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Unflashed material is advected into flashed material without
!   modification as there is no change in ye or the internal energy
!-----------------------------------------------------------------------

DO n = nmin+1, nmax
  deltx          = xa(n) - xa0(n)

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
!                  |               |
!                  |  fluxr(n) > 0 |
!         n-1      |       --------|-->         n
!        (j-7)     |               |          (j-6)
!                  n               n
!                   Eulr            Lagr
!
!                  <---- del x ---->
!
!       n-1 material in overlapping shell added to n
!
!-----------------------------------------------------------------------

    k            = n - 1
    j            = n - 7
    fractn       = 0.5d0 * deltx/dx(k)
    fractn2      = 1.d0 - fourthd * fractn

    fluxr(n)     = delta(n) * ( rl(k) + dr(k) - fractn * ( dr(k) - fractn2 * r6(k) ) )

    fluxu(n)     = ( ul(k)  + du(k)  -fractn * ( du(k)  -fractn2 * u6(k)  ) ) * fluxr(n)
    fluxv(n)     = ( vl(k)  + dv(k)  -fractn * ( dv(k)  -fractn2 * v6(k)  ) ) * fluxr(n)
    fluxw(n)     = ( wl(k)  + dw(k)  -fractn * ( dw(k)  -fractn2 * w6(k)  ) ) * fluxr(n)
    fluxe(n)     = ( el(k)  + de(k)  -fractn * ( de(k)  -fractn2 * e6(k)  ) ) * fluxr(n)
    fluxei(n)    = ( eil(k) + dei(k) -fractn * ( dei(k) -fractn2 * ei6(k) ) ) * fluxr(n)
    fluxye(n)    = ( yel(k) + dye(k) -fractn * ( dye(k) -fractn2 * ye6(k) ) ) * fluxr(n)
    fluxs(n)     = ( sl(k)  + ds(k)  -fractn * ( ds(k)  -fractn2 * s6(k)  ) ) * fluxr(n)

    IF ( nse(j,kj_ray,ki_ray) == 1 ) fluxbe(n) = eb(n) * fluxr(n)
    IF ( nse(j,kj_ray,ki_ray) == 0 ) fluxye(n) = fluxye_comp(n)

  ELSE ! deltx < 0.0d0

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
!                  |                 |
!                  |  fluxr(n) <  0  |
!        n-1    <--|---------        |            n
!       (j-7)      |                 |          (j-6)
!                  n                 n
!                   Lagr              Eulr
!
!
!                  <----- del x ----->
!
!         n material in overlapping shell added to n-1
!
!-----------------------------------------------------------------------

    j            = n - 6
    fractn       = 0.5d0 * deltx/dx(n)
    fractn2      = 1.d0 + fourthd * fractn

    fluxr(n)     = delta(n) * ( rl(n) - fractn * ( dr(n) + fractn2 * r6(n) ) )

    fluxu(n)     = ( ul(n)  - fractn * ( du(n)  + fractn2 * u6(n)  ) ) * fluxr(n)
    fluxv(n)     = ( vl(n)  - fractn * ( dv(n)  + fractn2 * v6(n)  ) ) * fluxr(n)
    fluxw(n)     = ( wl(n)  - fractn * ( dw(n)  + fractn2 * w6(n)  ) ) * fluxr(n)
    fluxe(n)     = ( el(n)  - fractn * ( de(n)  + fractn2 * e6(n)  ) ) * fluxr(n)
    fluxei(n)    = ( eil(n) - fractn * ( dei(n) + fractn2 * ei6(n) ) ) * fluxr(n)
    fluxye(n)    = ( yel(n) - fractn * ( dye(n) + fractn2 * ye6(n) ) ) * fluxr(n)
    fluxs(n)     = ( sl(n)  - fractn * ( ds(n)  + fractn2 * s6(n)  ) ) * fluxr(n)

    IF ( nse(j,kj_ray,ki_ray) == 1 ) fluxbe(n) = eb(n) * fluxr(n)
    IF ( nse(j,kj_ray,ki_ray) == 0 ) fluxye(n) = fluxye_comp(n)

  END IF ! deltx >= 0.0d0
END DO ! n = nmin+1, nmax

!-----------------------------------------------------------------------
!  Add diffusion to a shock aligned parallel to the z-axis
!-----------------------------------------------------------------------

DO n = nmin, nmax
  dm(n)          = r(n) * dvol(n)
END DO ! n = nmin, nmax

DO n = nmin, nmax-1
  j              = n - 6
  IF ( j_shock(j) == 1 ) THEN
    fluxr (n+1)  = fluxr (n+1) - v_diff * ( r (n+1) - r (n) ) * DMIN1( dvol(n), dvol(n+1) )
    fluxu (n+1)  = fluxu (n+1) - v_diff * ( u (n+1) - u (n) ) * DMIN1( dm(n), dm(n+1) )
    fluxv (n+1)  = fluxv (n+1) - v_diff * ( v (n+1) - v (n) ) * DMIN1( dm(n), dm(n+1) )
    fluxw (n+1)  = fluxw (n+1) - v_diff * ( w (n+1) - w (n) ) * DMIN1( dm(n), dm(n+1) )
    fluxye(n+1)  = fluxye(n+1) - v_diff * ( ye(n+1) - ye(n) ) * DMIN1( dm(n), dm(n+1) )
    fluxe (n+1)  = fluxe (n+1) - v_diff * ( e (n+1) + e_bind_zn0(n+1)/dm(n+1) &
&                - e (n) - e_bind_zn0(n)/dm(n) ) * DMIN1( dm(n),  dm(n+1) )
  END IF ! j_shock(j) == 1
END DO ! n = nmin, nmax-1

!-----------------------------------------------------------------------
!  Advect mass, momentum, and energy by moving the subshell quantities 
!   into the appropriate Eulerian zone. 
!-----------------------------------------------------------------------

DO n = nmin, nmax

  j              = n - 6

  dm(n)          = r(n) * dvol(n)
  dm0(n)         = ( dm(n) + fluxr(n) - fluxr(n+1) )
  r(n)           = dm0(n)/dvol0(n)
  r(n)           = DMAX1( smallr, r(n) )
  dm0(n)         = 1.0d0/( r(n) * dvol0(n) )

  u(n)           = ( u(n)  * dm(n) + fluxu(n)  - fluxu(n+1)  ) * dm0(n)
  v(n)           = ( v(n)  * dm(n) + fluxv(n)  - fluxv(n+1)  ) * dm0(n)
  w(n)           = ( w(n)  * dm(n) + fluxw(n)  - fluxw(n+1)  ) * dm0(n)
  es(n)          = ( ei(n) * dm(n) + fluxei(n) - fluxei(n+1) ) * dm0(n)
  se(n)          = ( se(n) * dm(n) + fluxs(n)  - fluxs(n+1)  ) * dm0(n)
  e(n)           = ( e(n)  * dm(n) + e_bind_zn0(n) + fluxe(n) + fluxbe(n) &
&                - fluxe(n+1) - fluxbe(n+1) ) * dm0(n)

  IF ( nse(j,kj_ray,ki_ray) == 1 ) THEN
    ye(n)        = ( ye(n) * dm(n) + fluxye(n) - fluxye(n+1) ) * dm0(n)
  END IF

END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
!  Grid change requires updated parabolic coeff 
!-----------------------------------------------------------------------

CALL sweepbc( nleftz, nrightz, nmin, nmax, ki_ray, kj_ray )
CALL paraset( ntot, zparaz, dx0, xa0, nmin-3, nmax+3, ngeom )
CALL e_decompose( nmin, nmax )

IF ( n_ei > 0 ) THEN
  DO n = nmin,nmin+n_ei
    ei(n)          = es(n)
  END DO
END IF

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (ul, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ul        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (u6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (du, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'du        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (vl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vl        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (v6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dv, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dv        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (wl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wl        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (w6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dw        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (rl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rl        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (r6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (el, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'el        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (e6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (de, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (ei_b, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei_b      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (eil, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'eil       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ei6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ei6       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dei, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dei       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (yel, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'yel       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ye6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ye6       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dye       '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (sl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'sl        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (s6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 's6        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ds, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ds        '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (dm, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dm0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (delta, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (fluxr, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxr     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxu, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxu     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxv, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxv     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxw, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxw     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxe, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxe     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxei, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxei    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxye, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxye    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxs, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxs     '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (var1_i, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'var1_i    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (var2_i, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'var2_i    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (var3_i, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'var3_i    '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (es, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'es        '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE remap_z
