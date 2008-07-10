SUBROUTINE linear_fit( x, y, ndim, n_min, n_max, sig, nwt, a, b )
!-----------------------------------------------------------------------
!
!    File:         linear_fit
!    Module:       linear_fit
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         12/25/07
!
!    Purpose:
!      Given a set of ndata points x(i), y(i), with standard deviations
!       sig(i), fit them to a straight line y = a + bx by minimizing
!       chi^{2}. If nwt = 0 on input, then the standard deviations are
!       assumed to be unavailable.
!
!    Input arguments:
!  x           : set of values of the independent variable
!  y           : set of values of the dependent variable
!  ndim        : extent of arrays x and y
!  n_min       : minimum index of data to be fit
!  n_max       : maximum indec of data to be fit
!  sig         : set of standard deviations
!  nwt         : standard deviation flag
!
!    Output arguments:
!  a           ; parameter a in the straight line fit y = a + bx
!  b           ; parameter b in the straight line fit y = a + bx
!
!    Subprograms called:
!      none
!
!    Include files:
!  kind_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY : zero, one
    
IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------
     
INTEGER, INTENT(in)                              :: ndim    ! extent of arrays x and y
INTEGER, INTENT(in)                              :: n_min   ! minimum index of data to be fit
INTEGER, INTENT(in)                              :: n_max   ! maximum index of data to be fit
INTEGER, INTENT(in)                              :: nwt     ! lower index of array

REAL(KIND=double), DIMENSION(ndim), INTENT(in)   :: x       ! set of values of the independent variable
REAL(KIND=double), DIMENSION(ndim), INTENT(in)   :: y       ! set of values of the dependent variable
REAL(KIND=double), DIMENSION(ndim)               :: sig     ! set of standard deviations

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)                   :: a       ! parameter a in the straight line fit y = a + bx
REAL(KIND=double), INTENT(out)                   :: b       ! parameter b in the straight line fit y = a + bx

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                                          :: i       ! do index

REAL(KIND=double)                                :: sx      ! sum of x(i) * wt(i)
REAL(KIND=double)                                :: sy      ! sum of y(i) * wt(i)
REAL(KIND=double)                                :: ss      ! weight sum
REAL(KIND=double), DIMENSION(ndim)               :: wt      ! weights
REAL(KIND=double)                                :: sxoss   ! sx/ss
REAL(KIND=double)                                :: t       ! ( x - sxoss )/sig(i)
REAL(KIND=double)                                :: st2     ! sum of t * t

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

sx                 = zero
sy                 = zero
st2                = zero
b                  = zero

!-----------------------------------------------------------------------
!  Weights
!-----------------------------------------------------------------------

IF ( nwt == 0 ) THEN
  wt               = one
  sig              = one
ELSE
  wt(n_min:n_max)  = one/( sig(n_min:n_max) * sig(n_min:n_max) )
END IF

!-----------------------------------------------------------------------
!  Accumulate sums
!-----------------------------------------------------------------------

ss                 = SUM( wt(n_min:n_max) )
sx                 = SUM( x(n_min:n_max) * wt(n_min:n_max) )
sy                 = SUM( y(n_min:n_max) * wt(n_min:n_max) )

!-----------------------------------------------------------------------
!  Solve for a and b
!-----------------------------------------------------------------------

sxoss              = sx/ss

DO i = n_min, n_max
  t                = ( x(i) - sxoss )/sig(i)
  st2              = st2 + t * t
  b                = b + t * y(i)/sig(i)
END DO

b                  = b/st2
a                  = ( sy - sx * b )/ss

RETURN
END SUBROUTINE linear_fit
