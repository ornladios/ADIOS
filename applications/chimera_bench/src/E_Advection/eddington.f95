SUBROUTINE eddington( jr_min, jr_max )
!-----------------------------------------------------------------------
!
!    File:         eddington
!    Module:       eddington
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         11/02/01
!
!    Purpose:
!      To compute the flux and Eddington factors.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  jr_min        : inner zone for which calculation w is to be made
!  jr_max        : outer zone for which calculation w is to be made
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  psi0(j,k,n) : zeroth angular moment of the neutrino distribution
!                function
!  psi1(j,k,n) : first angular moment of the neutrino distribution
!                function
!  nnugp(n)    : number of neutrino energy groups for n-neutrinos
!
!
!    Output arguments (common):
!  flxf(j,k,n) : flux factor as a function of j,k,n
!  Ef(j,k,n)   : Eddington factor as a function of j,k,n
!
!    Include files:
!  kind_module, array_module, numerical_module
!  e_advct_module, nu_energy_grid_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nnu
USE numerical_module, ONLY : half, one, epsilon

USE e_advct_module, ONLY : psi0, psi1, flxf, Ef
USE nu_energy_grid_module, ONLY : nnugp, nnugpmx

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: jr_min          ! minimum radial zone index
INTEGER, INTENT(in)              :: jr_max          ! maximum radial zone index

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: j             ! radial zone index
INTEGER                          :: k             ! neutrino energy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: kmax          ! nnugp(n)

REAL(KIND=double)                :: fft           ! flux factor (psi^(1)/psi^(0))

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Compute the flux and Eddington factors
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  DO k = 1,nnugpmx
    DO n = 1,nnu
      IF ( nnugp(n) == 0 ) CYCLE
      fft            = DABS( psi1(j,k,n) )/( half * ( psi0(j,k,n) + psi0(j+1,k,n) ) + epsilon )
      fft            = DMAX1( DMIN1( fft, one ), -one )
      flxf(j,k,n)    = fft
      Ef(j,k,n)      = ( one + 2.d0 * fft**2 )/3.d0
    END DO !  n = 1,nnu
  END DO ! k = 1,nnugp(n)
END DO ! j = jr_min,jr_max

!-----------------------------------------------------------------------
!  Boundary values
!-----------------------------------------------------------------------

DO j = jr_min,jr_max
  DO n = 1,nnu
    IF ( nnugp(n) == 0 ) CYCLE
    kmax             = nnugp(n)
    flxf(j,kmax+1,n) = flxf(j,kmax,n)
    Ef(j,kmax+1,n)   = Ef(j,kmax,n)
  END DO !  n = 1,nnu
END DO ! j = jr_min,jr_max

RETURN
END SUBROUTINE eddington
