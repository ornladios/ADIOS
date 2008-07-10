SUBROUTINE bremcal( enu, enubar, rho, t, xn, xp, s_a )
!-----------------------------------------------------------------------
!
!    File:         bremcal
!    Module:       bremcal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/23/02
!
!    Purpose:
!      To compute the NN Bremsstrahlung by the analytic fitting
!       formula given by Hannestad & Raffelt 1998, Apj, 507, 339.
!
!    Subprograms called:
!  g_brem, s_brem
!
!    Input arguments:
!  e_sum   : enu + enubar (MeV)
!  rho     : matter density (g/cm3)
!  tmev    : matter temperature (MeV)
!  xn      : free neutron mass fraction
!  xp      : free proton mass fraction
!
!    Output arguments:
!  s_p     : differential neutrino pair production kernal
!
!    Include files: 
!  kind_module, numerical_module, physcnst.cmn
!
!----------------------------------------------------------------------c

USE kind_module
USE numerical_module, ONLY : half
USE physcnst_module, ONLY : cvel, kmev, pi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(in)  :: enu                  ! incoming neutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: enubar               ! incoming antineutrino energy (MeV)
REAL(KIND=double), INTENT(in)  :: rho                  ! density (g/cm3)
REAL(KIND=double), INTENT(in)  :: t                    ! temperature (K)
REAL(KIND=double), INTENT(in)  :: xn                   ! free neutron mass fraction
REAL(KIND=double), INTENT(in)  :: xp                   ! free proton mass fraction

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out) :: s_a                  ! differential absorption kernal

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                        :: first = .true.

REAL(KIND=double), PARAMETER   :: tthird = 2./3.

REAL(KIND=double)              :: e_sum                ! enu + enubar (MeV)
REAL(KIND=double)              :: C_AG_FnB             ! Raffelt constant
REAL(KIND=double)              :: conv1                ! 1/(hc)**3 in natural units
REAL(KIND=double)              :: conv2                ! integrates ( 3 - costh ) over solid angle
REAL(KIND=double)              :: coef                 ! C_AG_FnB * conv1 * conv2
REAL(KIND=double)              :: rho_14               ! rho/1e+14
REAL(KIND=double)              :: t_10                 ! t/1e+10
REAL(KIND=double)              :: x                    ! e_sum/tmev
REAL(KIND=double)              :: tmev                 ! kt (MeV)
REAL(KIND=double)              :: eta_star             ! degeneracy parameter
REAL(KIND=double)              :: gamma                ! spin fluctuation rate
REAL(KIND=double)              :: y                    ! pion mass parameter
REAL(KIND=double)              :: sb                   ! dimensionless fitting parameter
REAL(KIND=double)              :: gb                   ! dimensionless fitting parameter

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Initialize....................................................

IF ( first ) THEN
  C_AG_FnB         = 3.79d+4/cvel
  conv1            = 1.d0 / ( 2.d0 * pi )**3
  conv2            = 3 * 4 * pi
  coef             = conv1 * conv2 * C_AG_FnB
  first            = .false.
END IF

rho_14             = rho/1.d+14
t_10               = t/1.d+10
tmev               = kmev * t
e_sum              = enu + enubar

!........Compute eta_star, the neutron effective degeneracy parameter..

eta_star           = 3.04d0 * rho_14**tthird / t_10

!........Compute gamma, the spin fluctuation rate......................

gamma              = 8.6d0 * rho_14 / sqrt(t_10)

!........Compute y, the pion mass parameter.

y                  = 1.94d0 / t_10

!........Compute x, the dimensionless neutrino energy sum..............

x                  = e_sum/tmev

!........Compute the dimensionless fitting parameter, sb...............

CALL s_brem(x,y,eta_star,sb)

!........Compute the dimensionless fitting parameter, gb...............

CALL g_brem(y,eta_star,gb)

!........Compute the differential absorption kernal, s_a...............

s_a                = gamma/( x**2 + ( half * gamma * gb )**2 )         &
&                  * sb/tmev * coef * rho_14 * ( xn + xp )


RETURN
END SUBROUTINE bremcal
