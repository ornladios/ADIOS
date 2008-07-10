Function func_x2(x)
!-----------------------------------------------------------------------
!
!    File:         func_x2
!    Module:       func_x2
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

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: x            ! independent variable

!-----------------------------------------------------------------------
!        Function variable.
!-----------------------------------------------------------------------

REAL(KIND=double)              :: func_x2      ! x**2

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

func_x2            = x * x

END FUNCTION func_x2
