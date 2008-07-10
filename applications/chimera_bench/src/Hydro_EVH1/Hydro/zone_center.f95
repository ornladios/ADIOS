SUBROUTINE zone_center( ngeom, nmn, nmx, xedge, xcen )

USE kind_module, ONLY: double
USE array_module, ONLY: max_12
     
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                               :: ngeom    ! geometry index
INTEGER, INTENT(in)                               :: nmn      ! minimum zone index to obtain zone center average
INTEGER, INTENT(in)                               :: nmx      ! maximum zone index to obtain zone center average

REAL(KIND=double), DIMENSION(max_12), INTENT(in)  :: xedge    ! position of zone edge

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(max_12), INTENT(out) :: xcen     ! position of zone center

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                           :: n        ! padded zone index

REAL(KIND=double), PARAMETER                      :: third = 1.d0/3.d0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Based on the geometry, calculate volume center..................

SELECT CASE (ngeom)

  CASE(0,4:)

    DO n = nmn,nmx 
      xcen(n)       = 0.5d0 * ( xedge(n+1) + xedge(n) )
    END DO
    
  CASE(1) ! Cylindrical geometry, watch negative r

    DO n = nmn,nmx 
      xcen(n)       = 0.5d0 * ( xedge(n+1)**2 + xedge(n)**2 )
      xcen(n)       = DSIGN( DSQRT(DABS(xcen(n))), xcen(n) )
    END DO

  CASE(2) ! Spherical geometry, watch negative r

    DO n = nmn,nmx 
      xcen(n)       = 0.5d0 * ( xedge(n+1)**3 + xedge(n)**3 )
      xcen(n)       = DSIGN( DABS(xcen(n))**third, xcen(n) )
    END DO

END SELECT

RETURN
END SUBROUTINE zone_center
