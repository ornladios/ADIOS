FUNCTION func_FD(x)
!-----------------------------------------------------------------------
!
!    File:         func_FD
!    Module:       func_FD
!    Type:         Subrogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/16/02
!
!    Purpose:
!      Fermi-Dirac function.
!
!                           n
!                         x
!        func_FD = ------------------
!                  exp( x - eta ) + 1
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
!        Local variable.
!-----------------------------------------------------------------------

REAL(KIND=double)              :: func_FD      ! x**n_fermi/( 1 + exp(x - eta) )

!-----------------------------------------------------------------------
!        Function variable.
!-----------------------------------------------------------------------

REAL(KIND=double), EXTERNAL    :: fexp         ! exp with upper and lower bounds

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

func_FD            = x**n_FD/( one + fexp(x - eta_FD) )

END FUNCTION func_FD
