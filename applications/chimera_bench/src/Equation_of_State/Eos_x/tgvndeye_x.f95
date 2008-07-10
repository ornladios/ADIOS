SUBROUTINE tgvndeye_x( j, ij_ray, ik_ray, rho, e, ye, t_prev, t )
!-----------------------------------------------------------------------
!
!    File:         tgvndeye_x
!    Module:       tgvndeye_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/29/03
!
!    Purpose:
!      To compute the temperature given rho, s, and ye. Iteration
!       is by means of the bisection method if an initial guess
!       of the temperature is unavailable or is not within fraction
!       of the previous value, otherwise by Newton-Rhapson.
!
!    Variables that must be passed through common:
!        none
!
!    Subprograms called:
!  esrgn_x_t : Tests whether the temperature requires a new local cube
!               to be evaluated, and does so if required
!  eqstt_x   : interpolates between cube corners to a given state point
!
!    Input arguments:
!  j         : radial zone index.
!  ij_ray    : j-index of a radial ray
!  ik_ray    : k-index of a radial ray
!  rho       : density (g/cm**3).
!  e         : inernal energy (ergs/g)
!  ye        : electron fraction.
!  t_prev    : temperature guess.
!
!    Output arguments:
!  t         : temperature (K)
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_x_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: epsilon, half
USE array_module, ONLY : n_proc_y, ij_ray_dim, ik_ray_dim

USE edit_module, ONLY: nread, nprint, nlog
USE eos_snc_x_module, ONLY : nse
USE parallel_module, ONLY : myid

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

REAL(KIND=double), INTENT(in)    :: rho           ! ensity (g/cm**3)
REAL(KIND=double), INTENT(in)    :: e             ! inernal energy (ergs/g)
REAL(KIND=double), INTENT(in)    :: ye            ! electron fraction
REAL(KIND=double), INTENT(in)    :: t_prev        ! guess of the temperature

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: t             ! temperature (K)

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER, PARAMETER               :: itmax = 50    ! maximum number of iterations
INTEGER                          :: it            ! iteration index
INTEGER                          :: ivar = 2      ! EOS energy index
INTEGER                          :: j_ray         ! polar index of the radial ray
INTEGER                          :: k_ray         ! azimuthal index of the radial ray

REAL(KIND=double), PARAMETER     :: tol = 1.d-8   ! convergence criterion
REAL(KIND=double), PARAMETER     :: t_min = 1.d+07! minimum temperature for bisection
REAL(KIND=double), PARAMETER     :: t_max = 5.d+11! maximum temperature for bisection
REAL(KIND=double)                :: tmin          ! minimum temperature during bisection iteration
REAL(KIND=double)                :: tmax          ! maximum temperature during bisection iteration
REAL(KIND=double)                :: t_test        ! temperature guess for iteration
REAL(KIND=double)                :: e_test        ! energy computed using t_test
REAL(KIND=double)                :: dedd          ! derivative of the energy wrt density
REAL(KIND=double)                :: dedt          ! derivative of the temperature wrt density
REAL(KIND=double)                :: dedy          ! derivative of the electron fraction wrt density
REAL(KIND=double)                :: dt            ! increment in temperature

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' NR iteration for e will not converge in subroutine tgvndeye_x, will try Bisection')
 1003 FORMAT (' j=',i4,' j_ray=',i4,' k_ray=',i4,' myid=',i4,' rho=',es11.3, &
 & ' t_prev=',es11.3,' ye=',es11.3,' e=',es14.7,' e_test=',es14.7,' t_test=',es14.7)
 1005 FORMAT (' nes(',i4,')=',i4,' nes(',i4,')=',i4,' nes(',i4,')=',i4)
 2001 FORMAT (' it=',i4,' rho=',es12.4,' t_test=',es12.4,' ye=',es12.4, &
& ' e=',es12.4,' e_test=',es12.4,' tmin=',es12.4,' tmax=',es12.4)
 2003 FORMAT (' Bisection iteration for e will not converge in subroutine tgvndeye_x')
 2005 FORMAT (' j=',i4,' j_ray=',i4,' k_ray=',i4,' myid=',i4,' e=',es14.7, &
 & ' e_test=',es18.11,' tmin=',es18.11,' tmax=',es18.11)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!               \\\\\ NEWTON-RHAPSON ITERATION /////
!
!-----------------------------------------------------------------------

t_test             = t_prev

DO it = 1,itmax

!-----------------------------------------------------------------------
!  esrgn_x not called to avoid discontinuity in t
!-----------------------------------------------------------------------

  CALL eqstt_x( ivar, j, ij_ray, ik_ray, rho, t_test, ye, e_test, dedd, &
&  dedt, dedy )

  IF ( DABS( e - e_test ) <= tol * DABS(e)  .and.  it /= 1 ) THEN
    t              = t_test
    RETURN
  END IF

  dt               = ( e - e_test )/( dedt + epsilon )
  t_test           = DMIN1( DMAX1( t_test + dt, 0.9d0 * t_test ), 1.1d0 * t_test )

END DO

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_ray              = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
k_ray              = ( myid/n_proc_y ) * ik_ray_dim + ik_ray

!-----------------------------------------------------------------------
!  Edit diagnostics
!-----------------------------------------------------------------------

WRITE (nprint,1001)
WRITE (nprint,1003) j, j_ray, k_ray, myid, rho, t_prev, ye, e, e_test, t_test
WRITE (nprint,1005) j-1, nse(j-1,ij_ray,ik_ray), j, nse(j,ij_ray,ik_ray), &
& j+1, nse(j+1,ij_ray,ik_ray)
WRITE (nlog,1001)
WRITE (nlog,1003) j, j_ray, k_ray, myid, rho, t_prev, ye, e, e_test, t_test
WRITE (nlog,1005) j-1, nse(j-1,ij_ray,ik_ray), j, nse(j,ij_ray,ik_ray), &
& j+1, nse(j+1,ij_ray,ik_ray)

!-----------------------------------------------------------------------
!
!                  \\\\\ BISECTION ITERATION /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Set iteration boundaries
!-----------------------------------------------------------------------

tmin               = t_min
tmax               = t_max

!-----------------------------------------------------------------------
!  Iterate
!-----------------------------------------------------------------------

DO it = 1,itmax

  t_test           = half * ( tmin + tmax )

  CALL esrgn_x_t( j, ij_ray, ik_ray, rho, t_test, ye )
  CALL eqstt_x( ivar, j, ij_ray, ik_ray, rho, t_test, ye, e_test, dedd, &
&  dedt, dedy )

  WRITE (nlog,2001) it, rho, t_test, ye, e, e_test, tmin, tmax

  IF ( DABS( e - e_test ) <= tol * DABS(e) ) THEN
    t              = t_test
    RETURN
  END IF

  IF ( e_test <= e ) THEN
    tmin           = t_test
  ELSE
    tmax           = t_test
  END IF

END DO

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_ray              = MOD( myid, n_proc_y ) * ij_ray_dim + ij_ray
k_ray              = ( myid/n_proc_y ) * ik_ray_dim + ik_ray

!-----------------------------------------------------------------------
!  Edit diagnostics
!-----------------------------------------------------------------------

WRITE (nprint,2003)
WRITE (nprint,2005) j, j_ray, k_ray, myid, e, e_test, tmin, tmax
WRITE (nlog,2003)
WRITE (nlog,2005) j, j_ray, k_ray, myid, e, e_test, tmin, tmax
STOP

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE tgvndeye_x
