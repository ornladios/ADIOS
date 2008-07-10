SUBROUTINE parabola_nu( nmin, nmax, nnz, para, a, deltaa, a6, al, flat, &
& iflat, iskip, ngeom )
!-----------------------------------------------------------------------
! Colella and Woodward, JCompPhys 54, 174-201 (1984) eq 1.5-1.8,1.10
!
! parabola calculates the parabolas themselves. call paraset first
! for a given grid-spacing to set the constants, which can be reused
! each time parabola is called.
!
! flatening coefficients are calculated externally in flaten. 
!
! nmin/nmax are indicies over which the parabolae are calculated
!-----------------------------------------------------------------------

USE kind_module, ONLY: double

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                  :: nnz      ! array dimension
INTEGER, INTENT(in)                  :: nmin     ! minimum index over which the parabolae are calculated
INTEGER, INTENT(in)                  :: nmax     ! minimum index over which the parabolae are calculated
INTEGER, INTENT(in)                  :: iflat    ! flaten interpolation flag
INTEGER, INTENT(in)                  :: iskip    ! skip monotonicity flag
INTEGER, INTENT(in)                  :: ngeom    ! geometry parameter

REAL(KIND=double), INTENT(in), DIMENSION(10,nnz)  :: para     ! parabolic coefficients
REAL(KIND=double), INTENT(in), DIMENSION(nnz) :: a        ! working array
REAL(KIND=double), INTENT(in), DIMENSION(nnz) :: flat     ! flaten array

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(inout), DIMENSION(nnz) :: deltaa ! right - left
REAL(KIND=double), INTENT(inout), DIMENSION(nnz) :: al     ! left-hand value
REAL(KIND=double), INTENT(inout), DIMENSION(nnz) :: a6     ! parabolic coefficient

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                              :: n        ! padded zone index

REAL, PARAMETER                      :: twthrd = 2.0d0/3.0d0

REAL(KIND=double), DIMENSION(nnz)    :: ar       ! right-hand value
REAL(KIND=double), DIMENSION(nnz)    :: da       ! left-hand value - right-hand value
REAL(KIND=double), DIMENSION(nnz)    :: diffa    ! working array
REAL(KIND=double), DIMENSION(nnz)    :: scrch1   ! working array
REAL(KIND=double), DIMENSION(nnz)    :: scrch2   ! working array
REAL(KIND=double), DIMENSION(nnz)    :: almon    ! working array
REAL(KIND=double), DIMENSION(nnz)    :: armon    ! working array

REAL(KIND=double)                    :: onemfl   ! 

!-----------------------------------------------------------------------
!  Skip to monotonicity when interpolating total Energy for remap
!   since deltaa is already set up for us.
!-----------------------------------------------------------------------

IF ( iskip == 1 ) THEN
  DO n = nmin, nmax
    ar(n)            = deltaa(n) + al(n)
  END DO

ELSE

  DO n = nmin-2, nmax+1
    diffa(n)         = a(n+1) * para(6,n+1) - a(n) * para(6,n)
  END DO

!........Equation 1.7 of C&W............................................

!-----------------------------------------------------------------------
!     da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))
!     zero out da(n) IF a(n) is a local max/min
!-----------------------------------------------------------------------

  DO n = nmin-1, nmax+1
    IF ( diffa(n-1) * diffa(n) < 0.0d0 ) THEN
      da(n)          = 0.0
    ELSE
      da(n)          = para(4,n) * diffa(n) + para(5,n) * diffa(n-1)
      da(n)          = DSIGN( DMIN1( DABS(da(n)), 2.d0 * DABS(diffa(n-1) ), 2.d0 * DABS(diffa(n)) ), da(n) )
    END IF
  END DO

!........Equation 1.6 of C&W............................................

!-----------------------------------------------------------------------
!     a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * dma(j+1) + C3 * dma(j)
! MONOT: Limit ar(n) to the range defined by a(n) and a(n+1)
!-----------------------------------------------------------------------

  DO n = nmin-1, nmax
    ar(n)            = (a(n) + para(1,n) * diffa(n) + para(2,n) * da(n+1) + para(3,n) * da(n)) * para(10,n+1)
    ar(n)            = MAX( ar(n), MIN( a(n), a(n+1) ) )
    ar(n)            = MIN( ar(n), MAX( a(n), a(n+1) ) )
    al(n+1)          = ar(n)
  END DO
  
!-----------------------------------------------------------------------
! Reset zone-face value at the r=0 boundary  
!-----------------------------------------------------------------------

  IF ( ngeom == 1 ) THEN
    ar(nmin-1)       = 2.0d0 * a(nmin) - al(nmin+1)
    IF ( a(nmin-1) * a(nmin) < 0.0d0 ) ar(nmin-1) = 0.0d0  
    al(nmin)         = ar(nmin-1)
  ELSE IF ( ngeom==2 ) THEN
    ar(nmin-1)       = 2.5d0 * a(nmin) - 1.5d0 * al(nmin+1)
    if ( a(nmin-1) * a(nmin) < 0.0d0 ) ar(nmin-1) = 0.0d0  
    al(nmin)         = ar(nmin-1)
  ELSE IF ( ngeom==4 ) then
    ar(nmin-1)       = ( 6.0*a(nmin) + al(nmin+1 ) )/7.0d0
    if ( a(nmin-1) * a(nmin) < 0.0d0) ar(nmin-1) = 0.0d0
    al(nmin)         = ar(nmin-1)
  END IF

!-----------------------------------------------------------------------
! eqn. 4.1 - flaten interpolation in zones with a shock ( flat(n)->1. )
!-----------------------------------------------------------------------

  IF( iflat == 1 ) THEN
    DO n = nmin, nmax
      onemfl          = 1.0d0 - flat(n)
      ar(n)           = flat(n) * a(n) + onemfl * ar(n)
      al(n)           = flat(n) * a(n) + onemfl * al(n)
    END DO
  END IF

END IF ! iskip == 1

!-----------------------------------------------------------------------
! MONOTONICITY constraints:
!
! compute delta_a, a_6
! MONOT: IF a is a local max/min, flaten zone structure ar,al -> a.
! MONOT: compute monotonized values using eq. 1.10 of C&W
!        IF parabola exceeds al/ar, reset ar/al so that slope -> 0.
! Recalculate delta_a and a_6
!-----------------------------------------------------------------------

IF ( ngeom == 1 ) THEN

  DO n = nmin, nmax
    deltaa(n)         = ar(n) - al(n)
    a6(n)             = 6.0d0 * (a(n) - 0.5d0 * ( al(n) * ( 1.0d0 - para(7,n)) &
&                     + ar(n) * ( 1.0d0 + para(7,n) ) ) )
  END DO ! n = nmin, nmax

ELSE IF ( ngeom == 2 ) THEN

  DO n = nmin, nmax
    deltaa(n)         = ar(n) - al(n)
    a6(n)             = 6.0d0 * ( a(n) - 0.5d0 * ( al(n) * ( 1.0 - para(7,n) ) &
&                     + ar(n) * ( 1.0d0 + para(7,n) ) ) ) * para(8,n)
  END DO ! n = nmin, nmax

ELSE IF ( ngeom == 4 ) THEN

  DO n = nmin, nmax
    deltaa(n)         = ar(n) - al(n)
    a6(n)             = ( a(n) * ( para(8,n) - para(7,n) ) - al(n) * para(8,n) &
&                     + ar(n) * para(7,n))/para(9,n)
  END DO ! n = nmin, nmax

ELSE

  DO n = nmin, nmax
    deltaa(n)         = ar(n) - al(n)
    a6(n)             = 6.0d0 * ( a(n) - 0.5d0 * ( al(n) + ar(n) ) )
  END DO ! n = nmin, nmax

END IF ! ngeom == 1

!-----------------------------------------------------------------------
! MONOT: IF a is a local max/min, flaten zone structure ar,al -> a.
!-----------------------------------------------------------------------

DO n = nmin, nmax
  scrch1(n)           = ( ar(n) - a(n) ) * ( a(n) - al(n) ) 
  IF ( scrch1(n) < 0.0d0 ) THEN
    ar(n)             = a(n)
    al(n)             = a(n)
  END IF ! scrch1(n) < 0.0d0
END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
! MONOT: compute monotonized values using eq. 1.10 of C&W
!-----------------------------------------------------------------------

IF ( ngeom == 1 ) THEN

  DO n = nmin, nmax
    almon(n)          = ( 2.d0 * a(n) - ( 2d0 * twthrd + para(7,n)) * ar(n) )/( twthrd - para(7,n) )
    armon(n)          = ( 2.d0 * a(n) - ( 2d0 * twthrd - para(7,n)) * al(n) )/( twthrd + para(7,n) )
  END DO ! n = nmin, nmax

ELSE IF ( ngeom == 2 ) THEN

  DO n = nmin, nmax
    almon(n)          = ( ar(n) * ( 1.d0 + .5d0 * para(7,n) * para(8,n)) - para(7,n) * a(n) ) &
&                     / ( 1.0d0 - 0.5d0 * para(7,n) * (2.0d0 - para(8,n) ) )
    armon(n)          = ( al(n) * ( 1.d0 +.5d0 * para(7,n) * (2.0d0 - para(8,n) ) ) - para(7,n) * a(n) ) &
&                     /( 1.0d0 - 0.5d0 * para(7,n) * para(8,n) )
  END DO ! n = nmin, nmax

ELSE IF ( ngeom == 4 ) THEN

  DO n = nmin, nmax
    almon(n)          = ( ar(n) * ( para(7,n) - para(9,n) ) + a(n) * ( para(8,n) - para(7,n) ) ) &
&                     /( para(8,n) - para(9,n) )
    armon(n)          = ( al(n) * ( para(8,n) - para(9,n) ) - a(n) * ( para(8,n) - para(7,n) ) ) &
&                     /( para(7,n) - para(9,n) )
  END DO ! n = nmin, nmax

ELSE

  DO n = nmin, nmax
    almon(n)          = 3.0d0 * a(n) - 2.0d0 * ar(n)
    armon(n)          = 3.0d0 * a(n) - 2.0d0 * al(n)
  END DO ! n = nmin, nmax

END IF ! ngeom == 1


!-----------------------------------------------------------------------
! If parabola exceeds al/ar, reset ar/al so that slope -> 0.
!-----------------------------------------------------------------------

DO n = nmin, nmax
  scrch1(n)           = deltaa(n) * deltaa(n)
  scrch2(n)           = deltaa(n) * a6(n)
  IF ( scrch1(n) < +scrch2(n) ) al(n) = almon(n)       
  IF ( scrch1(n) < -scrch2(n) ) ar(n) = armon(n)
END DO ! n = nmin, nmax

!-----------------------------------------------------------------------
! Recalculate delta_a and a_6
!-----------------------------------------------------------------------

IF ( ngeom == 1 ) THEN

  DO n = nmin, nmax
    deltaa(n)         = ar(n) - al(n)
    a6(n)             = 6.0d0 * ( a(n) - 0.5d0 * ( al(n) * (1.0d0 - para(7,n) ) &
&                     + ar(n) * ( 1.0d0 + para(7,n) ) ) )
  END DO ! n = nmin, nmax

ELSE IF ( ngeom == 2 ) THEN

  DO n = nmin, nmax    
    deltaa(n)         = ar(n) - al(n)
    a6(n)             = 6.0d0 * ( a(n) - 0.5d0 * ( al(n) * ( 1.0d0 - para(7,n)) &
&                     + ar(n) * ( 1.0d0 + para(7,n) ) ) ) * para(8,n)
  END DO ! n = nmin, nmax

ELSE IF (ngeom==4) THEN

  DO n = nmin, nmax    
    deltaa(n)         = ar(n) - al(n)
    a6(n)             = ( a(n) * ( para(8,n) - para(7,n) ) - al(n) * para(8,n) + ar(n) *para(7,n) )/para(9,n)
  END DO ! n = nmin, nmax

ELSE

  DO n = nmin, nmax
    deltaa(n)        = ar(n) - al(n)
    a6(n)            = 6.0d0 * ( a(n) - 0.5d0 * ( al(n) + ar(n) ) )
  END DO ! n = nmin, nmax

END IF ! ngeom == 1

RETURN
END SUBROUTINE parabola_nu
