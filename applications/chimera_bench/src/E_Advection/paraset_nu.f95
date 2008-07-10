SUBROUTINE paraset_nu( nnz, para, dx, xa, nmin, nmax, ngeom )
!-----------------------------------------------------------------------
! Colella and Woodward, JCompPhys 54, 174-201 (1984) eq 1.6, 1.7
!
! paraset sets up constants which are re-used each time we want to
! interpolate a parabola on a quantity. First pull out constants
! A, B, and !, and THEN compute all the equations in terms of those
! quantities. 
!
! the quantities calculated here are stored in a array para,
! to be read by parabola()
!
! nmin/nmax are index range for which one will calculate parabolae
!-----------------------------------------------------------------------------------

USE kind_module, ONLY: double

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)                    :: nnz      ! array dimension
INTEGER, INTENT(in)                    :: nmin     ! minimum index over which the parabolae are calculated
INTEGER, INTENT(in)                    :: nmax     ! minimum index over which the parabolae are calculated
INTEGER, INTENT(in)                    :: ngeom    ! geometry parameter

REAL(KIND=double), INTENT(in), DIMENSION(nnz)    :: dx    ! zone width of zone
REAL(KIND=double), INTENT(in), DIMENSION(nnz+1)  :: xa    ! position of zone interfaces

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out), DIMENSION(10,nnz) :: para ! parabolic coefficients

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER                                :: n        ! padded zone index

REAL(KIND=double)                      :: y        ! working scalar
REAL(KIND=double)                      :: denomy   ! working scalar

REAL(KIND=double), DIMENSION(nnz)      :: a        ! working array, dX_j +  dX_j+1
REAL(KIND=double), DIMENSION(nnz)      :: b        ! working array, 2dX_j +  dX_j+1
REAL(KIND=double), DIMENSION(nnz)      :: c        ! working array, dX_j + 2dX_j+1
REAL(KIND=double), DIMENSION(nnz)      :: d        ! working array

REAL(KIND=double), DIMENSION(nnz)      :: ai       ! working array
REAL(KIND=double), DIMENSION(nnz)      :: bi       ! working array
REAL(KIND=double), DIMENSION(nnz)      :: ci       ! working array

!-----------------------------------------------------------------------
!        A =  dX_j +  dX_j+1
!        B = 2dX_j +  dX_j+1
!        C =  dX_j + 2dX_j+1
!        ai, bi, and ci are inverse quantities
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

DO n = nmin-2, nmax+1
  a (n) = dx(n) + dx(n+1)
  ai(n) = 1.0/a(n)
  b (n) = a(n) + dx(n)
  bi(n) = 1.0/b(n)
  c (n) = a(n) + dx(n+1)
  ci(n) = 1.0/c(n)
END DO

!                                        constants for equation 1.6
!     a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * da(j+1) + C3 * da(j)

DO n = nmin-1, nmax
  d(n)      = 1. / (a(n-1) + a(n+1))
  para(1,n) = dx(n) * ai(n) + 2. * dx(n+1) * dx(n) * d(n) * ai(n) * ( a(n-1) * bi(n) - a(n+1) * ci(n) )
  para(2,n) = - d(n) * dx(n)   * a(n-1) * bi(n)
  para(3,n) =   d(n) * dx(n+1) * a(n+1) * ci(n)
END DO

!                                        constants for equation 1.7
!     da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))

DO n = nmin-1, nmax+1
  d(n) = dx(n) / ( a(n-1) + dx(n+1) )
  para(4,n) = d(n) * b(n-1) * ai(n)
  para(5,n) = d(n) * c(n)   * ai(n-1)
END DO

!******************************************************************************
! Calculate geometry factors based on current ngeom:
!  (Only implement IF problem involves coordinate singularies)
!           para(6,n) = dvol/dx
!           para(7,n) = a6 factor
!           para(8,n) = a6 factor
!           para(9,n) = a6 factor
!           para(10,n) = 1 / x**alpha
!******************************************************************************

IF(ngeom==0) THEN
  DO n = nmin, nmax+2
    para(6,n) = 1.0
    para(10,n) = 1.0
  END DO
ELSE IF(ngeom==1) THEN
  DO n = nmin, nmax+2
    para(6,n) = xa(n) + 0.5*dx(n)
    para(7,n) = 1. / (6.*xa(n)/dx(n) + 3.)    ! = F in B&L
    IF(xa(n) == 0.0) THEN
      para(10,n) = 1. 
    ELSE
      para(10,n) = 1. / xa(n)
    END IF
  END DO
  para(6,nmin-1) = -para(6,nmin)
  para(6,nmin-2) = -para(6,nmin+1)
ELSE IF(ngeom==2) THEN
  DO n = nmin, nmax+2
    para(6,n) = xa(n) * (xa(n) + dx(n)) + dx(n)**2/3.
    y = xa(n) / dx(n)
    para(8,n) = (3.0*y*(y+1.0) + 1.0) / (3.0*y*(y+1.0) + 0.9) ! = G^{-1} in B&L
    para(7,n) = (y+0.5)/(3.0*y*(y+1.0) + 1.0)      ! = H in B&L
    IF(xa(n) == 0.0) THEN
       para(10,n) = 1.0
    ELSE
       para(10,n) = 1.0 / xa(n)**2
    END IF
  END DO
  para(6,nmin-1) = para(6,nmin)
  para(6,nmin-2) = para(6,nmin+1)
ELSE IF(ngeom==3) THEN
  DO n = nmin, nmax+2
    para(6,n) = 1.0
    para(10,n) = 1.0
  END DO
ELSE IF(ngeom==4) THEN
  DO n = nmin, nmax+2
    para(6,n) = (cos(xa(n)) - cos(xa(n+1)))/dx(n)
    para(7,n) = cos(xa(n+1)) - (sin(xa(n+1))-sin(xa(n)))/dx(n)
    para(8,n) = cos(xa(n))   - (sin(xa(n+1))-sin(xa(n)))/dx(n)
    para(9,n) = 2.0*para(6,n)/dx(n) - (sin(xa(n+1))+sin(xa(n)))/dx(n)
    IF( xa(n) == 0.0 ) THEN
       para(10,n) = 1.0d0
    ELSE
       para(10,n) = 1.0d0 / sin(xa(n))
    END IF
  END DO
ELSE IF(ngeom==5) THEN
  DO n = nmin, nmax+2
    para(6,n) = 1.0
    para(10,n) = 1.0
  END DO
END IF

RETURN
END SUBROUTINE paraset_nu
