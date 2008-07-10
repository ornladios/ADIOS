SUBROUTINE cmp_cal( brydns, tmev, y_nucleon, cmp_nucleon )
!-----------------------------------------------------------------------
!
!    File:         cmp_cal
!    Module:       cmp_cal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         4/15/00
!
!    Purpose:
!      To compute the chemical potential of free neutrons or
!       protons.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  brydns      : number density of nucleons (free or bound) (/fm3)
!  tmev        : temperature (MeV)
!  y_nucleon   : mass fraction of free neutrons or protons!  !
!
!    Output arguments:
!  cmp_nucleon : free nucleon chenical potential
!
!    Input arguments (common):
!        none
!
!    Output arguments (common)
!        none
!
!-----------------------------------------------------------------------

USE kind_module

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)    :: brydns        ! number density of baryons (fm^{-3})
REAL(KIND=double), INTENT(in)    :: tmev          ! temperature (MeV)
REAL(KIND=double), INTENT(in)    :: y_nucleon     ! ratio of free neutrons or protons to baryons

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: cmp_nucleon   ! nucleon chemical potential (MeV)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double), PARAMETER     :: c0  = 2.3798d-4
REAL(KIND=double), PARAMETER     :: g_n = 2.0d+0
REAL(KIND=double)                :: therm         ! temporary variable

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

therm              = c0 * tmev * dsqrt(tmev)/brydns
cmp_nucleon        = tmev * dlog( y_nucleon/( therm * g_n ) )

RETURN
END SUBROUTINE cmp_cal