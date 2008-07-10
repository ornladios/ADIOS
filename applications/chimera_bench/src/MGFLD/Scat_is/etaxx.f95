SUBROUTINE etaxx(rho,t,xn,xp,etann,etapp)
!-----------------------------------------------------------------------
!
!    File:         etaxx
!    Module:       etaxx
!    Type:         subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/15/93
!
!    Purpose:
!      To compute etann and etapp, the free neutron and free proton
!       number densities corrected for free neutron and free proton
!       blocking factors, respectively, as a function of rho, t,
!       xn, and xp. These quantities are interpolated between
!       expressions valid for nondegenerate and degenerate nucleons.
!
!    Input arguments:
!
!   rho       : matter density (g/cm**3)
!   t         : matter temperature (K)
!   xn        : mass fraction of free neutrons
!   xp        : mass fraction of free protons
!
!    Output arguments:
!
!   etann     : number density if free neutrons (ietann = 0); number
!                density of free neutrons corrected for the free
!                neutron blocking factor (ietann = 1).
!   etapp     : number density if free neutrons (ietapp = 0); number
!                density of free neutrons corrected for the free
!                neutron blocking factor (ietapp = 1).
!
!    Input arguments (common):
!
!   ietann    : 0, (1) degeneracy corrections to the free neutron and
!                free proton densities for isoenegetic scattering
!                omitted (computed).
!
!    Subprograms called:
!        none
!
!    Include files:
!  kind_module, numerical_module, physcnst_module
!  prb_cntl_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero
USE physcnst_module, ONLY : rmu, hbarc, mn, mp, pi, kmev

USE prb_cntl_module, ONLY : ietann

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(IN)    :: rho           ! density (g/cm3)
REAL(KIND=double), INTENT(IN)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(IN)    :: xn            ! neutron mass fraction
REAL(KIND=double), INTENT(IN)    :: xp            ! proton mass fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT)   :: etann         ! neutron number corrected for blocking
REAL(KIND=double), INTENT(OUT)   :: etapp         ! proton number corrected for blocking

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)                :: nn            ! neutron number uncorrected for blocking (cm^{-3})
REAL(KIND=double)                :: np            ! proton number uncorrected for blocking (cm^{-3})
REAL(KIND=double)                :: dn            ! neutron number uncorrected for blocking (fm^{-3})
REAL(KIND=double)                :: dp            ! proton number uncorrected for blocking (fm^{-3})
REAL(KIND=double)                :: efn           ! degenerate expression
REAL(KIND=double)                :: efp           ! degenerate expression
REAL(KIND=double)                :: etanndgnt     ! nondegenerate expression
REAL(KIND=double)                :: etappdgnt     ! nondegenerate expression

REAL(KIND=double), PARAMETER     :: tthird = 2.d0/3.d0

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  nn, np
!-----------------------------------------------------------------------

nn                 = xn * rho/rmu
np                 = xp * rho/rmu

!-----------------------------------------------------------------------
!  etann, etanp (zeroth approximation)
!-----------------------------------------------------------------------

etann              = nn
etapp              = np

!-----------------------------------------------------------------------
!  Return if ietann = 0, nn <= 0, or np <= 0
!-----------------------------------------------------------------------

IF ( ietann == 0  .or.  nn <= zero  .or.  np <= zero ) RETURN

!-----------------------------------------------------------------------
!  etann, etanp (analytic approximation)
!-----------------------------------------------------------------------

dn                 = nn * 1.d-39
dp                 = np * 1.d-39
efn                = ( hbarc**2/( 2.d+00 * mn ) ) * ( 3.d+00 * pi**2 * dn )**tthird
efp                = ( hbarc**2/( 2.d+00 * mp ) ) * ( 3.d+00 * pi**2 * dp )**tthird
etanndgnt          = 1.5d+00 * ( kmev * t/efn )
etappdgnt          = 1.5d+00 * ( kmev * t/efp )
etann              = nn * etanndgnt/DSQRT( 1.d+00 + etanndgnt**2 )
etapp              = np * etappdgnt/DSQRT( 1.d+00 + etappdgnt**2 )

RETURN
END SUBROUTINE etaxx
