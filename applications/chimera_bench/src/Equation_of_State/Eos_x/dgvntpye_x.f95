SUBROUTINE dgvntpye_x( j, ij_ray, ik_ray, t, p, ye, rho_prev, rho )
!-----------------------------------------------------------------------
!
!    File:         dgvntpye_x
!    Module:       dgvntpye_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/08/03
!
!    Purpose:
!      To compute the temperature given rho, p, and ye. Iteration
!       is by means of the bisection method if an initial guess
!       of the temperature is unavailable or is not within fraction
!       of the previous value, otherwise by Newton-Rhapson.
!
!    Variables that must be passed through common:
!  nse(j,ij_ray,ik_ray) : nuclear statistical equilibrium flag for radial zone j.
!
!    Subprograms called:
!  esrgn_x      : regenerates local EOS table if necessary
!  eqstt_x      : interpolates EOS values in local table
!
!    Input arguments:
!  j            : radial zone index.
!  ij_ray       : index denoting the j-index of a specific radial ray
!  ik_ray       : index denoting the k-index of a specific radial ray
!  rho_prev     : guess of the density (g/cm**3).
!  t            : temperature (K)
!  p            : pressure (dynes/cm2)
!  ye           : electron fraction.
!
!    Output arguments:
!  rho          : density (g/cm**3).
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_x_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: epsilon, half

USE edit_module, ONLY: nread, nprint
USE eos_snc_x_module, ONLY: nse

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray        ! index denoting the k-index of a specific radial ray

REAL(KIND=double), INTENT(in)    :: t             ! temperature (K)
REAL(KIND=double), INTENT(in)    :: p             ! pressure (dynes/cm2)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction
REAL(KIND=double), INTENT(in)    :: rho_prev      ! guess of the density (g/cm3)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: rho           ! density (g/cm3)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER               :: itmax = 40    ! maximum number of iterations
INTEGER                          :: it            ! iteration index
INTEGER                          :: ivarp = 1     ! EOS pressure index

REAL(KIND=double), PARAMETER     :: tol = 1.d-5   ! convergence criterion
REAL(KIND=double)                :: rho_min       ! minimum density for bisection
REAL(KIND=double)                :: rho_max       ! maximum density for bisection
REAL(KIND=double)                :: rho_test      ! density guess for iteration
REAL(KIND=double)                :: p_test        ! pressure computed using rho_test
REAL(KIND=double)                :: dpdd          ! derivative of the pressure wrt density
REAL(KIND=double)                :: dpdt          ! derivative of the pressure wrt temperature
REAL(KIND=double)                :: dpdy          ! derivative of the pressure wrt electron fraction
REAL(KIND=double)                :: drho          ! increment in density

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' NR iteration for p will not converge in subroutine dgvntpye, will try Bisection')
 1003 FORMAT (' p=',1pe14.7,' p_test=',1pe14.7,' rho_test=',1pe14.7,' rho_prev=',1pe14.7)
 2001 FORMAT (' Bisection iteration for p will not converge in subroutine dgvntpye')
 2003 FORMAT (' p=',1pe14.7,' p_test=',1pe14.7,' rho_test=',1pe14.7,' rho_prev=',1pe14.7)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

nse(j,ij_ray,ik_ray) = 1

!-----------------------------------------------------------------------
!  Newton-Rhapson iteration
!-----------------------------------------------------------------------

rho_test             = rho_prev

DO it = 1,itmax

  CALL esrgn_x( j, ij_ray, ik_ray, rho_test, t, ye )
  CALL eqstt_x( ivarp, j, ij_ray, ik_ray, rho_test, t, ye, p_test, dpdd, &
& dpdt, dpdy )

  IF ( DABS( p - p_test ) <= tol * p ) THEN
    rho              = rho_test
    RETURN
  END IF

  drho               = ( p - p_test )/( dpdd + epsilon )
  rho_test           = rho_test + drho

END DO

WRITE (nprint,1001)
WRITE (nprint,1003) p,p_test,rho_test,rho_prev

!-----------------------------------------------------------------------
!  Bisection iteration
!-----------------------------------------------------------------------

rho_test             = rho_prev
rho_min              = half * rho_prev
rho_max              = 2.d0 * rho_prev

DO it = 1,itmax

  rho_test           = half * ( rho_min + rho_max )

  CALL esrgn_x( j, ij_ray, ik_ray, rho_test, t, ye )
  CALL eqstt_x( ivarp, j, ij_ray, ik_ray, rho_test, t, ye, p_test, dpdd, dpdt, dpdy )

  IF ( DABS( p - p_test ) <= tol * p ) THEN
    rho              = rho_test
    RETURN
  END IF

  IF ( p_test <= p ) THEN
    rho_min          = rho_test
  ELSE
    rho_max          = rho_test
  END IF

END DO

WRITE (nprint,2001)
WRITE (nprint,2003) p,p_test,rho_test,rho_prev

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE dgvntpye_x
