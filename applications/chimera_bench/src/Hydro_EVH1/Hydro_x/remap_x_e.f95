SUBROUTINE remap_x_e( ngeom, imin, imax, ij_ray, ik_ray, nx )
!-----------------------------------------------------------------------
!
!    File:         remap_x_e
!    Module:       remap_x_e
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/28/07
!
!    Purpose:
!      To remap the total energy (internal, kinetic, and gravitational) 
!       from the updated lagrangian grid to the fixed final grid, using
!       piecewise parabolic functions.
!
!    Input arguments:
!
!  ngeom       : geometry index
!  imin        : minimum x-array index
!  imax        : maximim x-array index
!  ij_ray      : index denoting the j-index of a specific radial ray
!  ik_ray      : index denoting the k-index of a specific radial ray
!  nx          : x-array dimension
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

!    Include files:
!  kind_module, numerical_module, physcnst_module
!  edit_module, eos_snc_x_module, evh1_global, evh1_sweep,
!  evh1_zone, mdl_cnfg_module, mgfld_remap_module, nucbrn_module,
!  nu_dist_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: zero, frpi, epsilon
USE physcnst_module, ONLY: pi

USE edit_module, ONLY : nlog
USE eos_snc_x_module, ONLY: nse
USE evh1_global, ONLY: smallr, nleftx, nrightx
USE evh1_sweep, ONLY: xa, dx, r, u, v, w, ei, e_v, e, xa0, dx0, &
& radius, dvol, dvol0
USE evh1_zone, ONLY : zparax
USE mgfld_remap_module, ONLY : fluxbe, e_bind_zn0, eb

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                          :: ngeom      ! geometry index
INTEGER, INTENT(in)                          :: imin       ! minimum x-array index
INTEGER, INTENT(in)                          :: imax       ! maximim x-array index
INTEGER, INTENT(in)                          :: ij_ray     ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)                          :: ik_ray     ! index denoting the k-index of a specific radial ray
INTEGER, INTENT(in)                          :: nx         ! x-array extent

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
INTEGER                                      :: j          ! mgfld radial zone index
INTEGER, PARAMETER                           :: n_ei = 0   ! number of inner zones to remap internal energy

REAL(KIND=double)                            :: ur         ! x-velocity on right-hand side of zone
REAL(KIND=double)                            :: vr         ! y-velocity on right-hand side of zone
REAL(KIND=double)                            :: wr         ! z-velocity on right-hand side of zone
REAL(KIND=double)                            :: evr        ! inernal energy on right-hand side of zone
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

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: e_v_b      ! partial energy minus the nuclear binding energy
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: evl        ! partial energy at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: ev6        ! partial energy parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dev        ! partial energy slope

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dm         ! mass after lagrangian step
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: dm0        ! mass after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: delta      ! volume of overlapping subshells

REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxr      ! mass contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxu      ! x-momentum contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxv      ! y-momentum contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxw      ! z-momentum contained in overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:) :: fluxe      ! energy contained in overlapping subshells

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in remap_x')
 2001 FORMAT (' Deallocation problem for array ',a10,' in remap_x')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (ul(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ul        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (u6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'u6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (du(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'du        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (vl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'vl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (v6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'v6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dv(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dv        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (wl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'wl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (w6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'w6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dw(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dw        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (rl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'rl        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (r6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'r6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dr(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dr        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (el(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'el        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (e6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e6        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (de(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'de        '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (e_v_b(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_v_b     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (evl(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'evl       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (ev6(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ev6       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dev(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dev       '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (dm(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm0(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (delta(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (fluxr(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxr     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxu(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxu     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxv(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxv     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxw(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxw     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxe(nx+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxe     '; WRITE (nlog,1001) var_name; END IF

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

e_v_b            = zero
evl              = zero
ev6              = zero
dev              = zero

dm               = zero
dm0              = zero
delta            = zero

fluxr            = zero
fluxu            = zero
fluxv            = zero
fluxw            = zero
fluxe            = zero

IF ( first ) THEN
  first          = .false.
  third          = 1.d0/3.d0
  twothd         = 2.d0 * third
  fourthd        = 4.d0 * third
  fopi           = 4.d0 * pi
END IF

nmin             = imin + 6
nmax             = imax + 6
ntot             = imax + 12

!-----------------------------------------------------------------------
!  Assemble the total energy (including gravitational) from the partial
!   energy and the kinetic energy, and store in e(n)
!-----------------------------------------------------------------------

CALL e_compose_v( nmin, nmax )

!-----------------------------------------------------------------------
!  Subtract the nuclear binding energy from the internal energy if the
!   matter is not in nse.
!  If the matter is in nse, subtract the binding energy of the first
!   non-nse zone to ensure continuity at the nse - non-nse boundary.
!-----------------------------------------------------------------------

DO n = nmax+3, nmin-3, -1
  e_v_b(n)       = e_v(n) - eb(n)
  e(n)           = e(n)   - eb(n)
END DO

!-----------------------------------------------------------------------
!
! Generate interpolation functions, saving da, al for constructing 
! left and right total energy states. (dmu is passed as dummy array)
!
!
!        a(x) = a    + x[ da  = a   (1 - x) ]
!                L,j        j    6,j
!
!        x = ( xi - xi     )/dxi       xi     < xi < xi
!                     j-1/2     j        j-1/2         j+1/2
!
!-----------------------------------------------------------------------

CALL parabola( nmin-3, nmax+3, ntot, zparax, r,     dr,  r6,  rl,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparax, u,     du,  u6,  ul,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparax, v,     dv,  v6,  vl,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparax, w,     dw,  w6,  wl,  dm, 0, &
& 0, ngeom )
CALL parabola( nmin-3, nmax+3, ntot, zparax, e_v_b, dev, ev6, evl, dm, 0, &
& 0, ngeom )

!-----------------------------------------------------------------------
!  Use the profiles for density, pressure, and velocities to calculate
!   consistent values of the left and right values of total energy
!-----------------------------------------------------------------------

DO n = nmin-1, nmax+1
  ur             = ul(n)     + du(n)
  vr             = vl(n)     + dv(n)
  wr             = wl(n)     + dw(n)
  evr            = evl(n)    + dev(n)
  el(n)          = evl(n) + 0.5d0 * ( ul(n)**2 + vl(n)**2 + wl(n)**2 )
  er             = evr    + 0.5d0 * ( ur**2    + vr**2    + wr**2    )
  de(n)          = er - el(n)
END DO

!-----------------------------------------------------------------------
!  Use the constructed left and right states to get parabolas for total
!   energy
!-----------------------------------------------------------------------

CALL parabola( nmin-3, nmax+3, ntot, zparax, e, de, e6, el, dm, 0, 1,   &
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
  delta(n)       = ( xa(n) - xa0(n) ) * radius 
END DO

END IF

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

DO n = nmin, nmax + 1
  deltx          = xa(n) - xa0(n)

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
!                  |               |
!                  |  fluxr(n) > 0 |
!         n-1      |       --------|-->         n
!        (j-6)     |               |          (j-5)
!                  n               n
!                   Eulr            Lagr
!
!                  <---- del x ---->
!
!       n-1 material in overlapping shell added to n
!
!-----------------------------------------------------------------------

    k            = n - 1
    j            = n - 6
    fractn       = 0.5d0 * deltx/dx(k)
    fractn2      = 1.d0 - fourthd * fractn

    fluxr(n)     = delta(n) * ( rl(k) + dr(k) - fractn * ( dr(k) - fractn2 * r6(k) ) )

    fluxu(n)     = ( ul(k)  + du(k)  -fractn * ( du(k)  -fractn2 * u6(k)  ) ) * fluxr(n)
    fluxv(n)     = ( vl(k)  + dv(k)  -fractn * ( dv(k)  -fractn2 * v6(k)  ) ) * fluxr(n)
    fluxw(n)     = ( wl(k)  + dw(k)  -fractn * ( dw(k)  -fractn2 * w6(k)  ) ) * fluxr(n)
    fluxe(n)     = ( el(k)  + de(k)  -fractn * ( de(k)  -fractn2 * e6(k)  ) ) * fluxr(n)

    IF ( nse(j,ij_ray,ik_ray) == 1 ) fluxbe(n) = eb(n) * fluxr(n)

  ELSE

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
!       (j-6)      |                 |          (j-5)
!                  n                 n
!                   Lagr              Eulr
!
!
!                  <----- del x ----->
!
!         n material in overlapping shell added to n-1
!
!-----------------------------------------------------------------------

    j            = n - 5
    fractn       = 0.5d0 * deltx/dx(n)
    fractn2      = 1.d0 + fourthd * fractn

    fluxr(n)     = delta(n) * ( rl(n) - fractn * ( dr(n) + fractn2 * r6(n) ) )

    fluxu(n)     = ( ul(n)  - fractn * ( du(n)  + fractn2 * u6(n)  ) ) * fluxr(n)
    fluxv(n)     = ( vl(n)  - fractn * ( dv(n)  + fractn2 * v6(n)  ) ) * fluxr(n)
    fluxw(n)     = ( wl(n)  - fractn * ( dw(n)  + fractn2 * w6(n)  ) ) * fluxr(n)
    fluxe(n)     = ( el(n)  - fractn * ( de(n)  + fractn2 * e6(n)  ) ) * fluxr(n)

    IF ( nse(j,ij_ray,ik_ray) == 1 ) fluxbe(n) = eb(n) * fluxr(n)

  END IF
END DO

!-----------------------------------------------------------------------
!  Advect mass, momentum, and energy by moving the subshell quantities 
!   into the appropriate Eulerian zone. 
!-----------------------------------------------------------------------

DO n = nmin, nmax

  j              = n - 5

  dm(n)          = r(n) * dvol(n)
  dm0(n)         = ( dm(n) + fluxr(n) - fluxr(n+1) )
  r(n)           = dm0(n)/dvol0(n)
  r(n)           = DMAX1( smallr, r(n) )
  dm0(n)         = 1.0d0/( r(n) * dvol0(n) )

  u(n)           = ( u(n)  * dm(n) + fluxu(n)  - fluxu(n+1)  ) * dm0(n)
  v(n)           = ( v(n)  * dm(n) + fluxv(n)  - fluxv(n+1)  ) * dm0(n)
  w(n)           = ( w(n)  * dm(n) + fluxw(n)  - fluxw(n+1)  ) * dm0(n)
  e(n)           = ( e(n)  * dm(n) + e_bind_zn0(n) + fluxe(n) + fluxbe(n)  - fluxe(n+1) - fluxbe(n+1) ) * dm0(n)

END DO !  n = nmin, nmax

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

DEALLOCATE (e_v_b, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'e_v_b      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (evl, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'evl       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (ev6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'ev6       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dev, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dev       '; WRITE (nlog,2001) var_name; END IF

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

RETURN
END SUBROUTINE remap_x_e
