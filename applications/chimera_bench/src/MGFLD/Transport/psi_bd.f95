SUBROUTINE psi_bd( jr_max, k, n, r, nx, psrtio )
!-----------------------------------------------------------------------
!
!    File:         psi_bd
!    Module:       psi_bd
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         7/22/00
!
!    Purpose:
!      To apply the surface boundary conditions to compute the
!       ratio of psi0(jr_max+1,k,n) to psi0(jr_max,k,n).
!
!        The surface boundary conditions are given by
!
!             psi1J  =  0.5 * ( psi0Jp + psi0Jm )
!             psi1J  = -dcJ * ( psi0Jp - psi0Jm) /raJ
!
!     --->    psi0Jp = psrtio * psi0Jm
!
!  where
!
!             psrtio = ( dcJ/raJ - 0.5 )/( dcJ/raJ + 0.5 )
!
!    Subprograms called:
!        none
!
!    Input arguments:
!
!  jr_max  : outer radial zone for computing the diffusion coefficient.
!  k       : energy zone index
!  n       : the neutrino flavor index
!  r       : radius array (cm)
!  nx      : x-array extent
!
!    Output arguments:
!
!  psrtio  : ratio of psi0(jr_max+1,k,n) to psi0(jr_max,k,n)
!
!    Include files:
!  kind_module, numerical_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : half

USE nu_dist_module, ONLY : dc

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_max        ! maximum radial zone index
INTEGER, INTENT(in)              :: k             ! neutrino energy index
INTEGER, INTENT(in)              :: n             ! neutrino flavor index
INTEGER, INTENT(in)              :: nx            ! radial array extent

REAL(KIND=double), INTENT(in), DIMENSION(nx) :: r ! radius (cm)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: psrtio        ! psi0Jp = psrtio*psi0Jm + cpsi0Jp

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

REAL(KIND=double)                :: dc_r          ! dc/delta r

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute psrtio
!-----------------------------------------------------------------------

dc_r               = dc(jr_max,k,n)/( half * ( r(jr_max+1) - r(jr_max-1) ) )
psrtio             = ( dc_r - half )/( dc_r + half )

RETURN
END SUBROUTINE psi_bd
