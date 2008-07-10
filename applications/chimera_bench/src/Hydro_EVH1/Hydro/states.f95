SUBROUTINE states( ngeom, pl, ul, rl, gel, gcl, p6, u6, r6, ge6, gc6, &
& dp, du, dr, dge, dgc, plft, ulft, rlft, gelft, gclft, prgh, urgh, rrgh, &
& gergh, gcrgh, ij_ray, ik_ray )

!-----------------------------------------------------------------------
! This subroutine takes the values of rho, u, and P at the left hand
! side of the zone, the change accross the zone, and the parabolic 
! coefficients, p6, u6, and rho6, and computes the left and right states
! (integrated over the charachteristics) for each variable for input 
! to the Riemann solver.
!-----------------------------------------------------------------------
!
! Input variables are:
!  pl     : pressure at left edge of zone
!  dp     : pressure slope
!  p6     : pressure parabola coeffecient
!  ul     : velocity at left edge of zone
!  du     : velocity slope
!  u6     : velocity parabola coeffecient
!  rl     : density at left edge of zone
!  dr     : density slope
!  r6     : density parabola coeffecient
!  ij_ray : y-index of an x ray, x-index of an y ray, x index of a z array
!  ik_ray : z-index of an x ray, z-index of an y ray, y index of a z array
!
! Output variables are:
!  plft   : average pressure in the + wave state (left of interface)
!  prgh   : average pressure in the - wave state (right of interface)
!  ulft   : average velocity in the + wave state (left of interface)
!  urgh   : average velocity in the - wave state (right of interface)
!  rlft   : average density in the  + wave state (left of interface)
!  rrgh   : average density in the  - wave state (right of interface)
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE array_module, ONLY: max_12

USE edit_module, ONLY : nlog
USE evh1_global, ONLY : svel, smallp, smallr, dt
USE evh1_sweep, ONLY : nmin, nmax, radius, gc, dx, xa, p, r, u, v, w, &
& l_shock, l_rho, sweep

IMPLICIT none
SAVE

INTEGER, INTENT(in)              :: ij_ray   ! y-index of an x ray, x-index of an y ray, x index of a z array
INTEGER, INTENT(in)              :: ik_ray   ! z-index of an x ray, z-index of an y ray, y index of a z array

INTEGER                          :: n        ! zone index for padded array
INTEGER                          :: k        ! shifted zone index for padded array
INTEGER, INTENT(in)              :: ngeom    ! geometry index

REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: pl       ! pressure at left edge of zone
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: p6       ! pressure parabola coeffecient
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: dp       ! pressure slope
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: ul       ! velocity at left edge of zone
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: u6       ! velocity parabola coeffecient
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: du       ! velocity slope
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: rl       ! density at left edge of zone
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: r6       ! density parabola coeffecient
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: dr       ! density slope
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: gcl      ! gamma_1 at left edge of zone
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: gc6      ! gamma_1 parabola coeffecient
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: dgc      ! gamma_1 slope
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: gel      ! gamma_p at left edge of zone
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: ge6      ! gamma_p parabola coeffecient
REAL(KIND=double), INTENT(in), DIMENSION(max_12)  :: dge      ! gamma_p slope

REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: plft     ! average pressure in the + wave state (left of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: prgh     ! average pressure in the - wave state (right of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: ulft     ! average velocity in the + wave state (left of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: urgh     ! average velocity in the - wave state (right of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: rlft     ! average density in the + wave state (left of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: rrgh     ! average density in the - wave state (right of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: gclft    ! average gamma_1 in the + wave state (left of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: gcrgh    ! average gamma_1 in the - wave state (right of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: gelft    ! average gamma_p in the + wave state (left of interface)
REAL(KIND=double), INTENT(out), DIMENSION(max_12) :: gergh    ! average gamma_p in the - wave state (right of interface)

REAL(KIND=double), DIMENSION(max_12)              :: grav     ! gravitational force
REAL(KIND=double), DIMENSION(max_12)              :: fict     ! fictitious force
REAL(KIND=double), DIMENSION(max_12)              :: nuf      ! neutrino stress
REAL(KIND=double), DIMENSION(max_12)              :: Crhodt   ! fraction of a zone width traversed  by a wave
REAL(KIND=double), DIMENSION(max_12)              :: Cdtdx    ! fraction of a zone width traversed  by a wave
REAL(KIND=double), DIMENSION(max_12)              :: fCdtdx   ! 1 - 4/3 * Cdtdx

REAL(KIND=double)                                 :: C_sound  ! 4/3
REAL(KIND=double)                                 :: fothd    ! 4/3
REAL(KIND=double)                                 :: hdt      ! 0.5 * dt

fothd       = 4.d0/3.d0
hdt         = 0.5d0 * dt

!--------------------------------------------------------------------------
!
! Calculate the domain of dependence along space coordinate,
! c*dt, by multiplying the Eulerian sound speed in 
! each zone by dt.  Divide by two to save flops later.
!
! Cdtdx is half the distance a sound wave travels in time dt
! 
!--------------------------------------------------------------------------

DO n = nmin-4, nmax+4
  C_sound   = DSQRT( gc(n) * p(n)/r(n) )
  Cdtdx (n) = C_sound/( dx(n) * radius )
  Crhodt(n) = r(n) * C_sound * hdt
  svel      = DMAX1( svel, Cdtdx(n) )
  Cdtdx (n) = Cdtdx(n) * hdt
  fCdtdx(n) = 1.d0 - fothd * Cdtdx(n)
END DO

!--------------------------------------------------------------------------
! Obtain averages of rho, u, and P over the domain (+/-)Cdt
!                    lft is the + wave on the left  side of the boundary
!                    rgh is the - wave on the right side of the boundary
!
! Include gravitational, ficticious, and neutrino forces
!--------------------------------------------------------------------------

!
!..........................................................................
!           !.                 .!.                    !
!           ! .               . ! .                   !
!           !  .             .  !  .                  !
!           !   .           .   !   .                 !
!           !    .         .    !    .                !
!           !     .       .     !     .               !
!           !      .     .      !      .              !
!...........|...................|.....................|....................
!           n         n        n+1        n+1         n+2
!                               k          k          k+1
!
!--------------------------------------------------------------------------

CALL forces_edge( ij_ray, ik_ray, xa, u, v, w, r, grav, fict, nuf )

DO n = nmin-4, nmax+4

  k = n + 1
  plft(k)   = pl(n)  + dp(n)  - Cdtdx(n) * ( dp(n)  - fCdtdx(n) * p6(n)  )
  ulft(k)   = ul(n)  + du(n)  - Cdtdx(n) * ( du(n)  - fCdtdx(n) * u6(n)  )
  rlft(k)   = rl(n)  + dr(n)  - Cdtdx(n) * ( dr(n)  - fCdtdx(n) * r6(n)  )
  gelft(k)  = gel(n) + dge(n) - Cdtdx(n) * ( dge(n) - fCdtdx(n) * ge6(n) )
  gclft(k)  = gcl(n) + dgc(n) - Cdtdx(n) * ( dgc(n) - fCdtdx(n) * gc6(n) )
  plft(k)   = DMAX1( smallp, plft(k) )
  rlft(k)   = DMAX1( smallr, rlft(k) )
  plft(k)   = plft(k) + Crhodt(n) * ( grav(k) + fict(k) + nuf(k) )
  plft(k)   = DMAX1( smallp, plft(k) )

  prgh(n)   = pl(n)  + Cdtdx(n) * ( dp(n)  + fCdtdx(n) * p6(n)  )
  urgh(n)   = ul(n)  + Cdtdx(n) * ( du(n)  + fCdtdx(n) * u6(n)  )
  rrgh(n)   = rl(n)  + Cdtdx(n) * ( dr(n)  + fCdtdx(n) * r6(n)  )
  gergh(n)  = gel(n) + Cdtdx(n) * ( dge(n) + fCdtdx(n) * ge6(n) )
  gcrgh(n)  = gcl(n) + Cdtdx(n) * ( dgc(n) + fCdtdx(n) * gc6(n) )
  prgh(n)   = DMAX1( smallp, prgh(n) )
  rrgh(n)   = DMAX1( smallr, rrgh(n) )
  prgh(n)   = prgh(n) - Crhodt(n) * ( grav(n) + fict(n) + nuf(n) )
  prgh(n)   = DMAX1( smallp, prgh(n) )

END DO

RETURN
END
