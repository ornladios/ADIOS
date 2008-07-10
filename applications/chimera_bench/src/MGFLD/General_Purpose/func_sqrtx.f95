FUNCTION func_sqrtx(x)
!-----------------------------------------------------------------------
!
!    File:         func_sqrtx
!    Module:       func_sqrtx
!    Type:         Subrogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/02
!
!    Purpose:
!      Test function for subroutine int_gauslegdr. 
!
!    Subprograms called:
!      none
!
!    Mudules used: 
!  kind_module
!-----------------------------------------------------------------------

USE kind_module, ONLY : double

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: x            ! independent variable

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

REAL(KIND=double)              :: func_sqrtx   ! x**1/2

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

func_sqrtx     = DSQRT(x)

END FUNCTION func_sqrtx
