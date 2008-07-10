FUNCTION f10(x)
!-----------------------------------------------------------------------
!
!    File:         f10
!    Module:       f10
!    Type:         Function
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/28/02
!
!    Purpose:
!      To exponentiate a given value without over or underflow
!
!----------------------------------------------------------------------c

USE kind_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: x               ! function argument

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double)              :: f10             ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER   :: expmax = 300.d0 ! maximum absolute value of function argument

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

f10             = 10.d0 ** ( DMIN1( expmax, DMAX1( -expmax, x ) ) )

END FUNCTION f10