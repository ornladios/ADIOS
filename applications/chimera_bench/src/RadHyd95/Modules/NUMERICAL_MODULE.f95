!-----------------------------------------------------------------------
!    Module:       numerical_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!
!    Common numerical constants
!-----------------------------------------------------------------------

module numerical_module

USE kind_module
USE physcnst_module
SAVE


REAL(KIND=double), PARAMETER  :: zero      = 0.0d0
REAL(KIND=double), PARAMETER  :: half      = 0.5d0
REAL(KIND=double), PARAMETER  :: one       = 1.0d0
REAL(KIND=double), PARAMETER  :: epsilon   = 1.0d-100

!-----------------------------------------------------------------------
!  Derived constants
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER  :: pi2       = pi * pi
REAL(KIND=double), PARAMETER  :: twpi      = 2.d0 * pi
REAL(KIND=double), PARAMETER  :: frpi      = 4.d0 * pi
REAL(KIND=double), PARAMETER  :: third     = 1.d0/3.d0
REAL(KIND=double), PARAMETER  :: frpith    = third * frpi
REAL(KIND=double), PARAMETER  :: sxtnpi    = 16.d0 * pi

!-----------------------------------------------------------------------
!  Coefficients for neutrino number and energy densities, and fluxes
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER  ::  ncoef    = 4.d+00 * pi/( h * cvel )**3
REAL(KIND=double), PARAMETER  ::  ecoef    = 4.d+00 * pi * ergmev/( h * cvel )**3
REAL(KIND=double), PARAMETER  ::  csqinv   = 1.d0/cvel**2

END module numerical_module
