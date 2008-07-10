FUNCTION ff( x ,y )
!-----------------------------------------------------------------------
!
!    File:         ff
!    Module:       ff
!    Type:         Function
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/28/02
!
!    Purpose:
!      To multiply two fermi functions
!
!----------------------------------------------------------------------c

USE kind_module, ONLY : double
USE numerical_module, ONLY : one

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(IN)  :: x               ! function argument
REAL(KIND=double), INTENT(IN)  :: y               ! function argument

!-----------------------------------------------------------------------
!        Output variables
!-----------------------------------------------------------------------

REAL(KIND=double)              :: ff              ! declare function

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)              :: fexp            ! exponential

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

ff              = one/( ( one + fexp(x) ) * ( one + fexp(y) ) )

END FUNCTION ff