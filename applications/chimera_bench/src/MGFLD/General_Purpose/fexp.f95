FUNCTION fexp(x)
!-----------------------------------------------------------------------
!
!    File:         fexp
!    Module:       fexp
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

USE kind_module, ONLY : double

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: x               ! function argument

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double)              :: fexp            ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER   :: expmax = 300.d0 ! maximum absolute value of function argument

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

fexp            = DEXP( DMIN1( expmax, DMAX1( -expmax, x ) ) )

END FUNCTION fexp