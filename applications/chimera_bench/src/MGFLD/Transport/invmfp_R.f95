SUBROUTINE invmfp_R( j, n, invmfp_R_j )
!-----------------------------------------------------------------------
!
!    File:         invmfp_R
!    Module:       invmfp_R
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         5/17/03
!
!    Purpose:
!      To compute the Rosseland inverse mean free path for radial zone j
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  j           : radial zone index
!  n           : neutrino type index
!
!    Output arguments:
!  mfp_R_j     : Rosseland mean inverse mean free path
!
!    Modules used:
!        kind_module, numerical_module, physcnst_module
!        mdl_cnfg_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero
USE physcnst_module, ONLY : kmev

USE mdl_cnfg_module, ONLY : t
USE nu_dist_module, ONLY : unu, dunu, rhs1
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(IN)            :: j                    ! radial zone index
INTEGER, INTENT(IN)            :: n                    ! neutrino type index

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(OUT) :: invmfp_R_j           ! Rosseland mean inverse mean free path 

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                        :: k                    ! neutrino energy index
REAL(KIND=double)              :: denom                ! denominator
REAL(KIND=double)              :: tmev                 ! temperature (MeV)
REAL(KIND=double)              :: int_num              ! integral in numerator
REAL(KIND=double)              :: int_denom            ! integral in denominator
REAL(KIND=double)              :: fexp                 ! exponential function

EXTERNAL fexp

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set Rosseland mean inverse mean free path to zero if nnugp = 0
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) THEN
  invmfp_R_j  = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Compute Rosseland mean inverse mean free path
!-----------------------------------------------------------------------

invmfp_R_j    = zero
int_num       = zero
int_denom     = zero
tmev          = kmev * t(j)

DO k = 1,nnugp(n)

  denom       = ( fexp( unu(j,k)/tmev ) + 1 )**2
  int_num     = int_num - rhs1(j,k,n) * unu(j,k)**4 * dunu(j,k) * denom
  int_denom   = int_denom + unu(j,k)**4 * dunu(j,k) * denom

END DO ! k

invmfp_R_j    = int_num/int_denom

RETURN
END SUBROUTINE invmfp_R
