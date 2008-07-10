SUBROUTINE rebal( jr_min, jr_max , n )
!-----------------------------------------------------------------------
!
!    File:         rebal
!    Module:       rebal
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         2/25/02
!
!    Purpose:
!      To prevent overfilling of neutrino occupation numbers.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  n           : the neutrino flavor index
!  jr_min      : inner zone for which calculation agr is to be
!                 made
!  jr_max      : outer zone for which calculation agr is to be
!                 made
!
!    Input arguments (common):
!  psi0(j,k,n) : zero angular moment of the neutrino occupation
!                 number
!  ncoefa(j,k) : 4.*pi/((h*!)**3)*w**2*dw
!  ecoefa(j,k) : 4.*pi/((h*!)**3)*w**3*dw
!
!    Output arguments:
!        none
!
!    Output arguments (common):
!  psi0(j,k,n) : zero angular moment of the neutrino occupation
!                 number
!
!    Include files:
!  kind_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY : zero, one

USE e_advct_module, ONLY : psi0, ncoefaa, ecoefaa
USE nu_energy_grid_module, ONLY : nnugp

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)               :: jr_min         ! minimum radial zone index
INTEGER, INTENT(in)               :: jr_max         ! maximum radial zone index
INTEGER, INTENT(in)               :: n              ! neutrino flavor index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                           :: j              ! radial zone index
INTEGER                           :: k              ! neutrino energy index

REAL(KIND=double)                 :: dpsk           ! amount by which psi0(j,k,n) exceeds 1
REAL(KIND=double)                 :: dpskp1m        ! change in psi0(j,k+1,n) if extra neutrinos are moved to k+1
REAL(KIND=double)                 :: dpskp1p        ! change in psi0(j,k+1,n) from psi0(j,k+2,n) to conserve energy
REAL(KIND=double)                 :: dpskp2         ! amout by which psi0(j,k+2,n) must be reduced to conserve energy
REAL(KIND=double)                 :: psi0k          ! psi0(j,k,n)
REAL(KIND=double)                 :: psi0kp1        ! psi0(j,k+1,n)
REAL(KIND=double)                 :: psi0kp2        ! psi0(j,k+2,n)

REAL(KIND=double)                 :: C1_dpsk        ! coefficient of dpsk in computing dpskp1m
REAL(KIND=double)                 :: C2_dpskp1m     ! coefficient of dpskp1m in computing dpskp2
REAL(KIND=double)                 :: C2_dpsk        ! coefficient of dpsk in computing dpskp2

REAL(KIND=double), PARAMETER      :: psi0min = 0.d0 ! minimum value of psi0(j,k+2,n)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Return if nnugp(n) <= 2
!-----------------------------------------------------------------------

IF ( nnugp(n) <= 2 ) RETURN


!-----------------------------------------------------------------------
!
!                   \\\\\ REPBALANCE PSI0 /////
!
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  DO k = 1,nnugp(n) - 2

!-----------------------------------------------------------------------
!  If psi0a > 1.0 and k le nnugp(n) - 2, set psi0a = 1.0 and shift
!   remaining n-neutrinos to the adjacent energy zone of the higher
!   energy so that the total n-neutrino number is conserved.
!   Then add some neutrinos to the adjacent energy zone of the higher
!   energy and remove some from the 2nd higher energy so that both
!   neutrino number and energy are conserved.
!-----------------------------------------------------------------------

    psi0k             = psi0(j,k,n)
    psi0kp1           = psi0(j,k+1,n)
    psi0kp2           = psi0(j,k+2,n)

    IF ( psi0k > one ) THEN

!-----------------------------------------------------------------------
!  Compute the number of neutrinos needed to be removed from the 2nd
!   higher energy zone.
!-----------------------------------------------------------------------

      dpsk            = psi0k - one
      C1_dpsk         = ncoefaa(j,k)/ncoefaa(j,k+1)
      C2_dpskp1m      = ecoefaa(j,k+1)/( ecoefaa(j,k+2) &
&                     - ecoefaa(j,k+1) * ncoefaa(j,k+2)/ncoefaa(j,k+1) )
      C2_dpsk         = - ecoefaa(j,k)/( ecoefaa(j,k+2) &
&                     - ecoefaa(j,k+1) * ncoefaa(j,k+2)/ncoefaa(j,k+1) )
      dpskp2          = ( C1_dpsk * C2_dpskp1m + C2_dpsk ) * dpsk

!-----------------------------------------------------------------------
!  If removing neutrinos from the 2nd higher energy zone would render
!   the distribution there negative, pick dpsk so that the distribution
!   there goes to zero.
!-----------------------------------------------------------------------
      
      IF ( psi0kp2 - dpskp2 < zero ) THEN
        dpskp2        = psi0kp2
        dpsk          = dpskp2/( C1_dpsk * C2_dpskp1m + C2_dpsk )
      END IF

!-----------------------------------------------------------------------
!  Now use the dpsk to compute everything else needed.
!-----------------------------------------------------------------------

      dpskp1m         = C1_dpsk * dpsk
      dpskp1p         = ncoefaa(j,k+2) * dpskp2/ncoefaa(j,k+1)
      psi0k           = psi0k - dpsk
      psi0kp1         = psi0kp1 + dpskp1m + dpskp1p
      psi0kp2         = psi0kp2 - dpskp2
      psi0(j,k,n)     = psi0k
      psi0(j,k+1,n)   = psi0kp1
      psi0(j,k+2,n)   = DMAX1( psi0kp2, psi0min )

    END IF

  END DO
END DO

RETURN
END SUBROUTINE rebal
