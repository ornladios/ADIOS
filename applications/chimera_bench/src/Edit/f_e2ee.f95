FUNCTION f_e2ee( eta )
!-----------------------------------------------------------------------
!
!    File:         f_e2ee
!    Module:       f_e2ee
!    Type:         Subrogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/13/00
!
!    Purpose:
!      To calculate F2*F4/F3**2.
!
!    Subprograms called:
!      F_eta
!
!    Input arguments:
!      none
!
!    Output arguments:
!      none
!
!    Variables that must be passed through common:
!      none
!
!    Include files:
!      none
!
!-----------------------------------------------------------------------

USE kind_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: eta          ! chemical potential / kT

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                        :: n            ! Fermi-function index

REAL(KIND=double)              :: f2           ! Fermi integral of order 2
REAL(KIND=double)              :: f3           ! Fermi integral of order 3
REAL(KIND=double)              :: f4           ! Fermi integral of order 4
REAL(KIND=double)              :: f_e2ee       ! function value

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

n           = 2
CALL F_eta( n, eta, f2 )
n           = 3
CALL F_eta( n, eta, f3 )
n           = 4
CALL F_eta( n, eta, f4 )
f_e2ee      = f2 * f4/f3**2

RETURN
END FUNCTION f_e2ee
