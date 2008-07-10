SUBROUTINE tgvndsye_x( j, ij_ray, ik_ray, rho, s, ye, t_prev, t )
!-----------------------------------------------------------------------
!                                                                      !
!    File:         tgvndsye_x
!    Module:       tgvndsye_x
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
!  esrgn_x  : test whether a new EOS-cube is needed
!  eqstt_x  : interpolates between cube corners to a given state point
!
!    Input arguments:
!  j        : radial zone index.
!  ij_ray   : j-index of a radial ray
!  ik_ray   : k-index of a radial ray
!  rho      : density (g/cm**3).
!  s        : entropy
!  ye       : electron fraction.
!  t_prev   : temperature guess.
!
!    Output arguments:
!  t        : temperature (K)
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: epsilon, half
USE array_module, ONLY : n_proc_y, ij_ray_dim, ik_ray_dim

USE edit_module, ONLY: nread, nprint, nlog
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
REAL(KIND=double), INTENT(in)    :: s             ! entropy
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
INTEGER                          :: ivar = 3      ! EOS entropy index
INTEGER                          :: j_ray         ! polar index of the radial ray
INTEGER                          :: k_ray         ! azimuthal index of the radial ray

REAL(KIND=double), PARAMETER     :: tol = 1.d-8   ! convergence criterion
REAL(KIND=double), PARAMETER     :: t_min = 1.d+07! minimum temperature for bisection
REAL(KIND=double), PARAMETER     :: t_max = 5.d+11! maximum temperature for bisection
REAL(KIND=double)                :: tmin          ! minimum temperature during bisection iteration
REAL(KIND=double)                :: tmax          ! maximum temperature during bisection iteration
REAL(KIND=double)                :: t_test        ! temperature guess for iteration
REAL(KIND=double)                :: s_test        ! energy computed using t_test
REAL(KIND=double)                :: dsdd          ! derivative of the entropy wrt density
REAL(KIND=double)                :: dsdt          ! derivative of the entropy wrt temperature
REAL(KIND=double)                :: dsdy          ! derivative of the entropy wrt electron fraction
REAL(KIND=double)                :: dt            ! increment in temperature

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

 1001 FORMAT (' NR iteration for s will not converge in subroutine tgvndsye_x, &
 &will try Bisection')
 1003 FORMAT (' j=',i4,' j_ray=',i4,' k_ray=',i4,' myid=',i4,' rho=',es11.3, &
 & ' t_prev=',es11.3,' ye=',es11.3,' s=',es14.7,' s_test=',es14.7,' t_test=',es14.7)
 2001 FORMAT (' it=',i4,' rho=',es12.4,' t_test=',es12.4,' ye=',es12.4, &
& ' s=',es12.4,' s_test=',es12.4,' tmin=',es12.4,' tmax=',es12.4)
 2003 FORMAT (' Bisection iteration for s will not converge in subroutine tgvndsye_x')
 2005 FORMAT (' j=',i4,' j_ray=',i4,' k_ray=',i4,' myid=',i4,' s=',es18.11, &
 & ' s_test=',es18.11,' tmin=',es18.11,' tmax=',es18.11)

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

  CALL eqstt_x( ivar, j, ij_ray, ik_ray, rho, t_test, ye, s_test, dsdd, &
&  dsdt, dsdy )
  IF ( DABS( s - s_test ) <= tol * DABS(s) ) THEN
    t              = t_test
    RETURN
  END IF

  dt               = ( s - s_test )/( dsdt + epsilon )
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
WRITE (nprint,1003) j, j_ray, k_ray, myid, rho, t_prev, ye, s, s_test, t_test
WRITE (nlog,1001)
WRITE (nlog,1003) j, j_ray, k_ray, myid, rho, t_prev, ye, s, s_test, t_test

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

  CALL esrgn_x( j, ij_ray, ik_ray, rho, t_test, ye )
  CALL eqstt_x( ivar, j, ij_ray, ik_ray, rho, t_test, ye, s_test, dsdd, &
&  dsdt, dsdy )

  WRITE (nlog,2001) it, rho, t_test, ye, s, s_test, tmin, tmax

  IF ( DABS( s - s_test ) <= tol * DABS(s) ) THEN
    t              = t_test
    RETURN
  END IF

  IF ( s_test <= s ) THEN
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
WRITE (nprint,2005) j, j_ray, k_ray, myid, s, s_test, tmin, tmax
WRITE (nlog,2003)
WRITE (nlog,2005) j, j_ray, k_ray, myid, s, s_test, tmin, tmax
STOP

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE tgvndsye_x
