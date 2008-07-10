SUBROUTINE snu( j, n, s )
!-----------------------------------------------------------------------
!
!    File:         snu
!    Module:       snu
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/14/00
!
!    Purpose:
!      To compute the neutrino entropy
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Output arguments:
!        none
!
!    Include files:
!  numerical_module, physcnst_module,
!  mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE numerical_module
USE physcnst_module, ONLY : rmu, pi, ergmev, cvel, h

USE mdl_cnfg_module, ONLY : rho
USE nu_dist_module, ONLY : unu, dunu, psi0
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: n             ! neutrino flavor index

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: s             ! neutrino entropy

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

LOGICAL                          :: first = .true.

INTEGER                          :: k             ! neutrino energy index

REAL(KIND=double)                :: coefs         ! coef/rho
REAL(KIND=double)                :: coef          ! coefficient
REAL(KIND=double)                :: psi0m         ! psi0

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

IF ( first ) THEN
  coef             = 4.d0 * rmu * pi * ergmev**3/( cvel * h )**3
  first            = .false.
END IF

!-----------------------------------------------------------------------
!  Calculate theneutrino entropy
!-----------------------------------------------------------------------

s                  = zero
IF ( nnugp(n) == 0 ) RETURN

coefs              = coef/rho(j)
DO k = 1,nnugp(n)
  psi0m            = DMAX1( 1.d-100, DMIN1( 0.9999999d+00, psi0(j,k,n) ) )
  s                = s + ( psi0m * DLOG(psi0m) + ( one - psi0m )        &
&                  * dlog( one - psi0m ) ) * unu(j,k)**2 * dunu(j,k)
END DO

s                  = - coefs * s

RETURN
END SUBROUTINE snu
