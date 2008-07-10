SUBROUTINE tdgvnspye_x( j, ij_ray, ik_ray, s, p, ye, rho_prev, t_prev, &
& rho, t )
!-----------------------------------------------------------------------
!
!    File:         tdgvnspye_x
!    Module:       tdgvnspye_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         6/08/03
!
!    Purpose:
!      To compute the temperature and rho given s, p, and ye. Iteration
!       is by means of Newton-Rhapson
!
!    Variables that must be passed through common:
!   nse(j,ij_ray, ik_ray): nuclear statistical equilibrium flag for radial zone j.
!
!    Subprograms called:
!  esrgn_x   : regenerates the local EOS table if necessary
!  eqstt_x   : interpolates quantities in the local EOS table
!
!    Input arguments:
!  j         : radial zone index
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!  s         : entropy
!  p         : pressure (dynes/cm2)
!  ye        : electron fraction
!  rho_prev  : input density (g/cm**3).
!  t_prev    : input temperature (K)
!
!    Output arguments:
!        rho       : density (g/cm**3).
!  t         : temperature (K)
!
!    Include files:
!  kind_module, array_module, numerical_module
!  eos_snc_x_module, edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: epsilon, half

USE eos_snc_x_module, ONLY: nse
USE edit_module, ONLY: nread, nprint
USE parallel_module, ONLY : myid

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)    :: s             ! entropy
REAL(KIND=double), INTENT(in)    :: p             ! pressure (dynes/cm2)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction
REAL(KIND=double), INTENT(in)    :: rho_prev      ! guess of the density (g/cm3)
REAL(KIND=double), INTENT(in)    :: t_prev        ! guess of the temperature (K)

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: t             ! temperature (K)
REAL(KIND=double), INTENT(out)   :: rho           ! density (g/cm3)

!-----------------------------------------------------------------------
!        Local variables.
!-----------------------------------------------------------------------

INTEGER, PARAMETER               :: itmax = 40    ! maximum number of iterations
INTEGER, PARAMETER               :: itrgn = 20    ! maximum number of iterations in which esrgn_x is called
INTEGER                          :: it            ! iteration index
INTEGER                          :: ivarp = 1     ! EOS pressure index
INTEGER                          :: ivars = 3     ! EOS entropy index

REAL(KIND=double), PARAMETER     :: tol = 1.d-5   ! convergence criterion
REAL(KIND=double)                :: t_test        ! temperature guess for iteration
REAL(KIND=double)                :: rho_test      ! energy computed using t_test
REAL(KIND=double)                :: p_test        ! iterated pressure
REAL(KIND=double)                :: s_test        ! iterated entropy
REAL(KIND=double)                :: dpdd          ! derivative of the pressure wrt density
REAL(KIND=double)                :: dpdt          ! derivative of the pressure wrt temperature
REAL(KIND=double)                :: dpdy          ! derivative of the pressure wrt electron fraction
REAL(KIND=double)                :: dsdd          ! derivative of the entropy wrt density
REAL(KIND=double)                :: dsdt          ! derivative of the entropy wrt the temperature
REAL(KIND=double)                :: dsdy          ! derivative of the entropy wrt the electron fraction
REAL(KIND=double)                :: dp            ! increment in pressure
REAL(KIND=double)                :: ds            ! increment in entropy

REAL(KIND=double), PARAMETER     :: two = 2.d0    ! convergence criterion
REAL(KIND=double)                :: det           ! determinant
REAL(KIND=double)                :: drho          ! increment in density
REAL(KIND=double)                :: dt            ! increment in temperature

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' p and s will not converge in subroutine tdgvnspye_x')
 1003 FORMAT (' j=',i4,' ij_ray=',i4,' ik_ray=',i4,' myid=',i4,' p=',1pe14.7,       &
 & ' p_test=',1pe14.7, ' s=',1pe14.7,' s_test=',1pe14.7)
 1005 FORMAT (' rho_test=',1pe14.7,' drho=',1pe14.7,' t_test=',1pe14.7, &
 & ' dt=',1pe14.7)

!-----------------------------------------------------------------------
!  Newton-Rhapson iteration
!-----------------------------------------------------------------------

nse(j,ij_ray,ik_ray) = 1
rho_test             = rho_prev
t_test               = t_prev

DO it = 1,itmax

!  IF ( it <= itrgn ) CALL esrgn_x(j,rho_test,t_test,ye)
  CALL eqstt_x( ivarp, j, ij_ray, ik_ray, rho_test, t_test, ye, p_test, &
&  dpdd, dpdt, dpdy )
  CALL eqstt_x( ivars, j, ij_ray, ik_ray, rho_test, t_test, ye, s_test, &
&  dsdd, dsdt, dsdy )
  
  dp                 = p - p_test
  ds                 = s - s_test

  IF ( DABS(dp) <= tol * p  .and.  DABS(ds) <= tol * s ) THEN
    rho              = rho_test
    t                = t_test
    RETURN
  END IF

  det                = dpdd * dsdt - dsdd * dpdt
  drho               = ( dsdt * dp - dpdt * ds )/det
  dt                 = ( dpdd * ds - dsdd * dp )/det

  rho_test           = rho_test + drho
  t_test             = t_test + dt
  rho_test           = DMAX1( half * rho_prev, DMIN1( two * rho_prev, rho_test ) )
  t_test             = DMAX1( half * t_prev  , DMIN1( two * t_prev  , t_test   ) )

END DO

WRITE (nprint,1001)
WRITE (nprint,1003) j, ij_ray, ik_ray, myid, p, p_test, s, s_test
WRITE (nprint,1005) rho_test, drho, t_test, dt

STOP
END SUBROUTINE tdgvnspye_x
