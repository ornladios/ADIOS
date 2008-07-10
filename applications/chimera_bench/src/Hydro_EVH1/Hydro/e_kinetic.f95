SUBROUTINE e_kinetic( nmin, nmax )
!=======================================================================  
!  Calculate the kinetic energy, ekin.
!=======================================================================

USE numerical_module, ONLY : half

USE evh1_sweep, ONLY: u, v, w, ekin
     
IMPLICIT none

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: nmin     ! minimum padded array index
INTEGER, INTENT(in)                    :: nmax     ! maximum padded array index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                :: n        ! padded zone index

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO n = nmin,nmax
  ekin (n)       = half * ( u(n) * u(n) + v(n) * v(n) + w(n) * w(n) )
END DO
     
RETURN
END SUBROUTINE e_kinetic
