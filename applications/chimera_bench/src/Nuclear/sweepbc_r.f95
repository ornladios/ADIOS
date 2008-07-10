SUBROUTINE sweepbc_r( nleft, nright, nmin, nmax )
!-----------------------------------------------------------------------
!
!    File:         sweepbc_r
!    Module:       sweepbc_r
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/02/04
!
!    Purpose:
!      To impose boundary conditions at the ends of the 1-D swweep arrays,
!       for use with composition remapping.
!
!    Subprograms called:
!  coord_bc         : computes the ghost coordinates
!
!    Input arguments:
!  nleft            : left-hand boundary flags
!  nright           : right-hand boundary flags
!  nmin             : minimum paddded array index
!  nmin             : maximum paddded array index
!
!    Output arguments:
!      none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  evh1_bound, evh1_global, evh1_sweep, evh1_zone,
!  mgfld_remap_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, half
USE physcnst_module, ONLY: pi, G => g

USE evh1_bound, ONLY : r_bcl, r_bcr
USE evh1_global, ONLY : smallp, ngeomx
USE evh1_sweep, ONLY : p, ge, gc, dvol, u, xa0, dx0
USE evh1_zone, ONLY : imax
USE mgfld_remap_module, ONLY : r, xa, dx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: nleft    ! left boundary condition key
INTEGER, INTENT(in)                    :: nright   ! right boundary condition key
INTEGER, INTENT(in)                    :: nmin     ! minimum paddded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum paddded array index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: n        ! zone index

REAL(KIND=double)                      :: massi    ! enclosed mass
REAL(KIND=double)                      :: mass_co  ! central mass
REAL(KIND=double)                      :: gcr      ! 1/gc(nmax
REAL(KIND=double)                      :: dpdr     ! pressure gradient
REAL(KIND=double)                      :: dp       ! change in pressure
REAL(KIND=double)                      :: prat     ! pressure ratio
REAL(KIND=double)                      :: rrat     ! density ratio

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : fixed (eg, u_bcl,p_bcr,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))
!-----------------------------------------------------------------------

!........Find ghost coordinates.........................................

CALL coord_bc( nleft, nright, nmin, nmax, xa, dx, xa0, dx0, imax+12 )

!........Load left (inner) ghosts.......................................

IF( nleft == 0 ) THEN          ! symmetric accross left (inner) edge

  DO n = 1, 6
    r (nmin-n)       = r (nmin+n-1)
  END DO

ELSE IF ( nleft == 1 ) THEN    ! Zero Gradient 

  DO n = 1, 6
    r (nmin-n)       = r (nmin)
  END DO

ELSE IF ( nleft == 2 ) THEN   ! Externally Fixed

  DO n = 1, 6
    r (nmin-n)       = r_bcl
  END DO

ELSE IF ( nleft == 3 ) THEN   ! Periodic

  DO n = 1, 6
    r (nmin-n)       = r (nmax+1-n)
  END DO

END IF ! nleft

!........Load right (outer) ghosts......................................

IF ( nright == 0 ) THEN          ! symmetric accross right (outer) edge

  DO n = 1, 6
    r (nmax+n)       = r (nmax+1-n)
  END DO

ELSE IF ( nright == 1 ) THEN     ! Zero Gradient

  DO n = 1, 6
    r (nmax+n)       = r (nmax)
  END DO

ELSE IF ( nright == 2 ) THEN     ! Externally Fixed

  DO n = 1, 6
    r (nmax+n)       = r_bcr
  END DO

ELSE IF ( nright == 3 ) THEN     ! Externally Fixed

  DO n = 1, 6
    r (nmax+n)       = r (nmin+n-1)
  END DO

ELSE IF ( nright == 4 ) THEN     ! Balance pressure with gravity  

!........Compute volume elements

  CALL volume ( ngeomx ) 

!........Calculate mass interior

  massi              = mass_co + 4.0d0 * pi * sum( r(nmin:nmax) * dvol(nmin:nmax) )
  gcr                = 1.0d0/gc(nmax)
  dpdr               = ( p(nmax) - p(nmax-1) )/( xa(nmax) - xa(nmax-1 ) + 0.5d0 * ( dx(nmax) - dx(nmax-1) ) )

  DO n = 1, 6

!........Gamma constant

    ge(nmax+n)       = ge(nmax)
    gc(nmax+n)       = gc(nmax)

!........Approximate the hydrostatic dp, using the previous density

    dp               = -G * massi * dx(nmax+n) * r(nmax+n-1)/xa(nmax+n)**2
    p(nmax+n)        = p(nmax+n-1) + dp
    p(nmax+n)        = DMAX1( smallp, p(nmax+n) )

!  Calculate P assuming constant dP/dr
!         p (nmax+n) = p (nmax)+ dpdr*
!    &      (xa(nmax+n)-xa(nmax)+.5*(dx(nmax+n)-dx(nmax)))
!	  p (nmax+n) = max(smallp,p(nmax+n))

!........Set rho to keep P/rho**gamma constant

    r(nmax+n)        = r(nmax+n-1) * (p(nmax+n)/p(nmax+n-1))**gcr

  END DO

ELSE IF ( nright == 5 ) THEN  ! Scaled Down P and rho

  prat               = p(nmax)/p(nmax-1)
  rrat               = r(nmax)/r(nmax-1)
  gcr                = 1.0d0/gc(nmax)                  

  DO n=1,6

!........Constant gamma

    ge(nmax+n)       = ge(nmax)
    gc(nmax+n)       = gc(nmax)

!........Scale down the pressure by the ratio of the last real zones

    p(nmax+n)        = p(nmax+n-1) * prat
    p(nmax+n)        = DMAX1( smallp, p(nmax+n) )

!........Set rho to keep P/rho**gamma constant

    r(nmax+n)        = r(nmax+n-1) * (p(nmax+n)/p(nmax+n-1))**gcr

  END DO ! nright

ELSE IF ( nright == 6 ) THEN     ! Externally Fixed

  IF ( u(nmax) < zero  .and.  r(nmax) < r_bcr ) THEN
    r_bcr            = r(nmax)
  ELSE IF ( u(nmax) >= zero ) THEN
    r_bcr            = r(nmax)
  END IF

  DO n = 1, 6
    r (nmax+n)       = r_bcr
  END DO

END IF

RETURN
END SUBROUTINE sweepbc_r

