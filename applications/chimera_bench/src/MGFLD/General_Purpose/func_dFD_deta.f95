FUNCTION func_dFD_deta(x)
!-----------------------------------------------------------------------
!
!    File:         func_dFD_deta
!    Module:       func_dFD_deta
!    Type:         Subrogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/02
!
!    Purpose:
!      Fermi-Dirac function.
!
!                                                      n
!                      d func_FD       exp( x - eta ) x
!        func_FDdeta = ---------  =  ---------------------
!                        d eta                           2
!                                    (exp( x - eta ) + 1)
!
!    Subprograms called:
!      none
!
!    Mudules used: 
!  kind_module, numerical_module
!  FD_module
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : one

USE FD_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: x             ! independent variable

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)              :: farg          ! exp(x - eta)
REAL(KIND=double)              :: func_dFD_deta ! x**n_fermi*farg/( 1 + farg )**2
REAL(KIND=double), EXTERNAL    :: fexp          ! exp with upper and lower bounds

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  farg            = fexp( x - eta_FD )
  func_dFD_deta   = x**n_FD * farg/( one + farg )**2

END FUNCTION func_dFD_deta
