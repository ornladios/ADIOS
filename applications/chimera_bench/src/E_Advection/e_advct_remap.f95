SUBROUTINE e_advct_remap( ngeom )
!-----------------------------------------------------------------------
!
!    File:         e_advct_remap
!    Module:       e_advct_remap
!    Type:         Subprogram
!    Author:       S. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/28/01
!
!    Purpose:
!      To remap the neutrino number from the updated lagrangian grid
!       to the fixed Eulerian grid, using piecewise parabolic functions.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ngeom       : geometry flag
!
!    Output arguments:
!        none
!
!    Input arguments (common):
!  nmin        : number of first real zone (=1+6)
!  nmax        : number of last  real zone (=nnugp(n)+6)
!  psi(k)      : neutrino occupation after lagrangian update
!  xk(k)       : energy grid zone edge locations  after lagrangian update
!  dk(k)       : xk(k+1) - xk(k) after lagrangian update
!  xk_1(k)     : energy grid zone edge locations at time m+1
!  dvolk(l)    : energy space volume  after lagrangian update
!  dvolk_1(l)  : energy space volume at time m+1
!
!    Output arguments (common):
!  psi(k)      : neutrino occupation at time m+1
!
!    Include files:
!  kind_module, array_module, numerical_module
!  e_advct_module, edit_module
!
!-----------------------------------------------------------------------

USE kind_module
USE array_module, ONLY : nez
USE numerical_module, ONLY : zero, half, third, one

USE e_advct_module, ONLY : nmin, nmax, dk, xk, xk_1, psi, dvolk, dvolk_1
USE edit_module, ONLY : nlog

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables
!-----------------------------------------------------------------------

INTEGER                                :: ngeom    ! geometry parameter

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=10)                     :: var_name

INTEGER                                :: k        ! neutrino energy index
INTEGER                                :: n        ! padded zone index

INTEGER                                :: istat         ! allocation status

REAL(KIND=double), PARAMETER           :: fourthd = 4.d0/3.d0

REAL(KIND=double)                      :: deltx    ! change in position
REAL(KIND=double)                      :: fractn   ! half the zone Lagrangian and Eulerian frid overlap
REAL(KIND=double)                      :: fractn2  ! 1 - 4/3*fractn

REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: psil     ! zero neutrino moment at left edge of zone
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: psi6     ! zero neutrino moment parabola coeffecient
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dpsi     ! zero neutrino moment slope
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dm       ! number of neutrinos after Lagrangian step
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: dm0      ! number of neutrinos after remap
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: delta    ! volume of overlapping subshells
REAL(KIND=double), ALLOCATABLE, DIMENSION(:)   :: fluxpsi  ! psi0 flux at the zone boundary

REAL(KIND=double), ALLOCATABLE, DIMENSION(:,:) :: para     ! array of parabolic interpolation coefficients

 1001 FORMAT (' Allocation problem for array ',a10,' in e_advct_remap')
 2001 FORMAT (' Deallocation problem for array ',a10,' in e_advct_remap')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (psil(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (psi6(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dpsi(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (dm0(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0       '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (delta(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (fluxpsi(nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi   '; WRITE (nlog,1001) var_name; END IF

ALLOCATE (para(10,nez+12), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para      '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

psil                 = zero
psi6                 = zero
dpsi                 = zero
dm                   = zero
dm0                  = zero
delta                = zero
fluxpsi              = zero

para                 = zero

!-----------------------------------------------------------------------
!  Generate interpolation functions. (dmu is passed as dummy array.)
!-----------------------------------------------------------------------

CALL paraset_nu( nmax+6, para, dk, xk, nmin-4, nmax+4, ngeom )
CALL parabola_nu( nmin-4, nmax+4, nmax+6, para, psi, dpsi, psi6, psil,  &
& dm, 0, 0, ngeom )

!-----------------------------------------------------------------------
!  Calculate the volume of the overlapping subshells (delta)
!-----------------------------------------------------------------------

DO n = nmin, nmax+1
  delta(n)           = xk(n) - xk_1(n)
  delta(n)           = delta(n) * ( xk_1(n) * (xk_1(n) + delta(n) ) + delta(n)**2 * third )
END DO

!-----------------------------------------------------------------------
!  Calculate the total number of neutrinos (fluxpsi) in the subshell
!   created by the overlap of the Lagrangian and Eulerican grids.
!   If the zone face has moved to the left (deltx > 0), use the 
!   integral from the left side of zone n (fluxpsir). If the zone face
!   has moved to the right (deltx < 0), use the integral from the
!   right side of zone k=n-1 (fluxrl).
!-----------------------------------------------------------------------

DO n = nmin, nmax + 1
  deltx              = xk(n) - xk_1(n)
  IF( deltx >= zero ) THEN
    k                = n - 1
    fractn           = half * deltx/dk(k)
    fractn2          = one - fourthd * fractn
    fluxpsi(n)       = delta(n) * ( psil(k) + dpsi(k) - fractn * ( dpsi(k) - fractn2 * psi6(k) ) )
  ELSE
    fractn           = half * deltx/dk(n)
    fractn2          = one + fourthd * fractn
    fluxpsi(n)       = delta(n) * ( psil(n) - fractn * ( dpsi(n) + fractn2 * psi6(n) ) )
  END IF
END DO

!-----------------------------------------------------------------------
!  Advect neutrinos by moving the subshell quantities into the
!   appropriate Eulerian zone.
!-----------------------------------------------------------------------

DO n = nmin, nmax
  dm (n)             = psi(n) * dvolk(n)
  dm0(n)             = ( dm(n) + fluxpsi(n) - fluxpsi(n+1) )
  psi(n)             = dm0(n)/dvolk_1(n)
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (psil, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psil      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (psi6, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'psi6      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dpsi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dpsi      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dm, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm        '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (dm0, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'dm0       '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (delta, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'delta     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (fluxpsi, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'fluxpsi   '; WRITE (nlog,2001) var_name; END IF

DEALLOCATE (para, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'para      '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE e_advct_remap
