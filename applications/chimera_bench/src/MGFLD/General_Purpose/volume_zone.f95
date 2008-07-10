SUBROUTINE volume_zone( ngeom, imin, imax, xa, dx, xa0, dx0, dvol, dvol0 )
!-----------------------------------------------------------------------
!
!    File:         volume_zone
!    Module:       volume_zone
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/03/05
!
!    Purpose:
!      To compute the zone volumes in the xa and xa0 grid.
!
!    Input arguments:
!  ngeom       : problem geometry flag
!  imin        : minimum x-array index
!  imax        : maximim x-array index
!  xa          : current value of the coordinate
!  dx          ; current value of the zone thickness
!  xa0         : initial value of the coordinate
!  dx0         : initial value of the zone thickness
!
!    Output arguments:
!  dvol        : current volume of mass shells
!  dvol0       : initial volume of mass shells
!
!    Subprograms called:
!      none
!
!    Include files:
!  kind_module
!  evh1_sweep
!
!-----------------------------------------------------------------------

USE kind_module

USE evh1_sweep, ONLY : radius
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------
     
INTEGER, INTENT(in)                                :: ngeom    ! geometry flag
INTEGER, INTENT(in)                                :: imin     ! lower index of array
INTEGER, INTENT(in)                                :: imax     ! upper index of array

REAL(KIND=double), DIMENSION(imax+12), INTENT(in)  :: xa       ! padded coordinate after Lagr update
REAL(KIND=double), DIMENSION(imax+12), INTENT(in)  :: dx       ! padded zone thickness after Lagr update
REAL(KIND=double), DIMENSION(imax+12), INTENT(in)  :: xa0      ! padded final coordinate
REAL(KIND=double), DIMENSION(imax+12), INTENT(in)  :: dx0      ! padded final zone thickness

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(imax+12), INTENT(out) :: dvol     ! volume after Lagr step
REAL(KIND=double), DIMENSION(imax+12), INTENT(out) :: dvol0    ! final volume

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: nmin     ! minimum padded array index
INTEGER                                :: nmax     ! maximum padded array index
INTEGER                                :: l        ! zone index

REAL(KIND=double), PARAMETER           :: third = 1.d0/3.d0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

nmin                 = imin + 6
nmax                 = imax + 6

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
    dvol (l)         = third * ( (xa (l) + dx (l) )**3 - xa (l)**3 )
    dvol0(l)         = third * ( (xa0(l) + dx0(l) )**3 - xa0(l)**3 )
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
END SUBROUTINE volume_zone
