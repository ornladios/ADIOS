SUBROUTINE volume( ngeom )

USE kind_module

USE evh1_sweep, ONLY : nmin, nmax, dx, dx0, xa, xa0, dvol, dvol0, radius
     
IMPLICIT none
SAVE

INTEGER                                :: l        ! zone index
INTEGER                                :: ngeom    ! geometry flag

REAL(KIND=double), PARAMETER           :: third = 1.d0/3.d0

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Calculate dvolume and average area based on geometry of sweep
!-----------------------------------------------------------------------

IF ( ngeom == 0 ) THEN

  radius             = 1.0d0

  DO l = nmin-3, nmax+4
    dvol (l)         = dx (l)
    dvol0(l)         = dx0(l)
  END DO

ELSE IF ( ngeom == 1 ) THEN

  radius             = 1.0d0

  DO l = nmin-3, nmax+4
    dvol (l)         = dx (l) * ( xa (l) + 0.5d0 * dx (l) )
    dvol0(l)         = dx0(l) * ( xa0(l) + 0.5d0 * dx0(l) )
  END DO

ELSE IF( ngeom == 2 ) THEN

  radius             = 1.0d0

  DO l = nmin-3, nmax+4
    dvol (l)         = dx(l)  * ( xa(l)  * ( xa(l)  + dx(l) )  + dx(l)  * dx(l)  * third ) 
    dvol0(l)         = dx0(l) * ( xa0(l) * ( xa0(l) + dx0(l) ) + dx0(l) * dx0(l) * third ) 
  END DO

ELSE IF ( ngeom == 3 ) THEN

  DO l = nmin-3, nmax+4
    dvol (l)         = dx (l) * radius
    dvol0(l)         = dx0(l) * radius
  END DO

ELSE IF ( ngeom == 4 ) THEN

  DO l = nmin-3, nmax+4
    dvol (l)         = ( cos(xa (l)) - cos(xa (l+1)) ) * radius
    dvol0(l)         = ( cos(xa0(l)) - cos(xa0(l+1)) ) * radius
  END DO

ELSE IF ( ngeom == 5 ) THEN

  DO l = nmin-3, nmax+4
    dvol (l)         = dx (l) * radius
    dvol0(l)         = dx0(l) * radius
  END DO

END IF              

RETURN
END SUBROUTINE volume
