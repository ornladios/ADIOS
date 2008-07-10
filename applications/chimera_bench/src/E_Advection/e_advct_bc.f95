SUBROUTINE e_advct_bc
!-----------------------------------------------------------------------
!
!    File:         e_advct_bc
!    Module:       e_advct_bc
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/27/01
!
!    Purpose:
!      To impose the boundary conditions at the ends of the array
!       nu(n)
!
!    Subprograms called:
!        none
!
!    Input arguments:
!        none
!
!    Input arguments (common)
!  nmin        : number of first real zone (=1+6)
!  nmax        : number of last  real zone (=nnugp(n)+6)
!  xk(k)       : energy grid zone edge locations
!  dk(k)       : xk(k+1) - xk(k)
!  xk_0(k)     : energy grid zone edge locations prior to lagrangian
!                 update
!  dk_0(k)     : dk(k) prior to lagrangian update
!  xk_1(k)     : energy grid zone edge locations at time m+1
!  dk_1(k)     : dk(k) time at m+1
!  psi(k)      : neutrino occupation number
!
!    Output arguments:
!        none
!
!    Output arguments (common):
!  xk(k)       : ghost zones filled with bc values
!  dk(k)       : ghost zones filled with bc values
!  xk_0(k)     : ghost zones filled with bc values
!  dk_0(k)     : ghost zones filled with bc values
!  xk_1(k)     : ghost zones filled with bc values
!  dk_1(k)     : ghost zones filled with bc values
!  psi(k)      : ghost zones filled with bc values
!
!    Include files:
!  e_advct_module
!
!-----------------------------------------------------------------------

USE e_advct_module, ONLY : nmin, nmax, xk_0, dk_0, xk_1, dk_1, xk, dk, &
& vk, psi

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                           :: n             ! zone padding index
INTEGER                           :: nmax1n        ! nmax+1-n
INTEGER                           :: nminn1        ! nmin+n-1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Impose lower boundary conditions
!-----------------------------------------------------------------------

DO n = 1,6
  nminn1            = MIN( nmin + n - 1, nmax )
  dk  (nmin-n)      = dk  (nminn1)
  xk  (nmin-n)      = xk  (nminn1) - dk  (nmin-n)
  dk_0(nmin-n)      = dk_0(nminn1)
  xk_0(nmin-n)      = xk_0(nminn1) - dk_0(nmin-n)
  dk_1(nmin-n)      = dk_1(nminn1)
  xk_1(nmin-n)      = xk_1(nminn1) - dk_1(nmin-n)
  psi (nmin-n)      = psi (nminn1)
END DO

!-----------------------------------------------------------------------
!  Impose upper boundary conditions
!-----------------------------------------------------------------------

DO n = 1,6
  nmax1n            = MAX( nmax + 1 - n, nmin )
  dk  (nmax+n)      = dk  (nmax1n)
  xk  (nmax+n)      = xk  (nmax1n) + dk (nmax+n-1)
  dk_0(nmax+n)      = dk_0(nmax1n)
  xk_0(nmax+n)      = xk_0(nmax1n) + dk_0(nmax+n-1)
  dk_1(nmax+n)      = dk_1(nmax1n)
  xk_1(nmax+n)      = xk_1(nmax1n) + dk_1(nmax+n-1)
  psi (nmax+n)      = psi (nmax1n)
END DO

RETURN
END SUBROUTINE e_advct_bc
