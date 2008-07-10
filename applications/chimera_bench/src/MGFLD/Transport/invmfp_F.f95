SUBROUTINE invmfp_F( j, n, invmfp_F_j )
!-----------------------------------------------------------------------
!
!
!    File:         invmfp_R
!    Module:       invmfp_R
!    Type:         Subroutine
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/21/03
!
!    Purpose:
!      To compute the flux inverse mean free path for radial zone j
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
!  invmfp_F_j  : Flux mean inverse mean free path
!
!    Modules used:
!  kind_module, numerical_module
!  edit_module, nu_dist_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, epsilon

USE edit_module, ONLY : nlog
USE nu_dist_module, ONLY : unu, dunu, rhs1, psi0
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)            :: j                    ! radial zone index
INTEGER, INTENT(in)            :: n                    ! neutrino type index

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out) :: invmfp_F_j           ! flux averaged mean inverse mean free path 

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

CHARACTER (len=10)             :: var_name

INTEGER                        :: istat                ! allocation status

REAL(KIND=double),ALLOCATABLE  :: flux_k(:)            ! flux of eneregy zone k
REAL(KIND=double)              :: int_num              ! integral in numerator
REAL(KIND=double)              :: int_denom            ! integral in denominator

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' Allocation problem for array ',a10,' in invmfp_F')
 2001 FORMAT (' Deallocation problem for array ',a10,' in invmfp_F')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set flux mean inverse mean free path to zero if nnugp = 0.
!-----------------------------------------------------------------------

IF ( nnugp(n) == 0 ) THEN
  invmfp_F_j       = zero
  RETURN
END IF

!-----------------------------------------------------------------------
!  Allocate flux_k
!-----------------------------------------------------------------------

ALLOCATE (flux_k(nnugp(n)), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_k    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

flux_k             = zero

!-----------------------------------------------------------------------
!  Compute the flux mean inverse mean free path
!-----------------------------------------------------------------------

invmfp_F_j         = zero
int_num            = zero
int_denom          = zero

flux_k(1:nnugp(n)) = psi0(j,1:nnugp(n),n) * unu(j,1:nnugp(n))**4 * dunu(j,1:nnugp(n))
int_num            = - SUM(rhs1(j,1:nnugp(n),n) * flux_k)
int_denom          = SUM(flux_k,1)

invmfp_F_j         = int_num/( int_denom + epsilon )

!-----------------------------------------------------------------------
!  Deallocate flux_k
!-----------------------------------------------------------------------

DEALLOCATE (flux_k, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'flux_k    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE invmfp_F
