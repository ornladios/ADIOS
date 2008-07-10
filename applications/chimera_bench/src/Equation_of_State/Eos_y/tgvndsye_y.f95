SUBROUTINE tgvndsye_y( j, ji_ray, jk_ray, rho, s, ye, t_prev, t )
!-----------------------------------------------------------------------
!
!    File:         tgvndsye_y
!    Module:       tgvndsye_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         3/19/05
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
!  esrgn_y_t : recomputes, if the temperature warrents, the EOS table entries
!  eqstt_y   : interpolates between cube corners to a given state point
!
!    Input arguments:
!  j         : y (angular) zone index.
!  ji_ray    : x (radial) index of a specific y (angular) ray
!  jk_ray    : z (azimuthal) index of a specific y (angular) ray
!  rho       : density (g cm^{-3}).
!  s         : inernal energy (ergs g^{-1})
!  ye        : electron fraction.
!  t_prev    : temperature guess.
!
!    Output arguments:
!  t         : temperature (K)
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_y_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY: n_proc_y, ij_ray_dim, j_ray_dim, ik_ray_dim, ny
USE numerical_module, ONLY: epsilon, half

USE edit_module, ONLY: nprint, nlog
USE eos_snc_y_module, ONLY : nse
USE parallel_module, ONLY : myid, myid_y, myid_z

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! y (angular) zone index
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific y (angular) ray

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
INTEGER                          :: ivar = 3      ! EOS energy index
INTEGER                          :: j_angular     ! angular ray index
INTEGER                          :: k_angular     ! azimuthal ray index
INTEGER                          :: j_ray_bndl    ! polar index of the radial ray bundle
INTEGER                          :: k_ray_bndl    ! azimuthal index of the radial ray bundle
INTEGER                          :: jm1           ! mimimum of j-1 and 1
INTEGER                          :: jp1           ! maximum of j+1 and ny

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

 1001 FORMAT (' NR iteration for s will not converge in subroutine tgvndsye_y, &
 &will try Bisection')
  1003 FORMAT (' j=',i4,' j_ray_bndl=',i4,' k_ray_bndl=',i4,' j_angular=',i4, &
& ' k_angular=',i4,' myid=',i4,' rho=',es11.3,' t_prev=',es11.3, ' ye=',es11.3, &
& ' s=',es20.11,' s_test=',es20.11,' t_test=',es20.11)
 1005 FORMAT (' nes(',i4,')=',i4,' nes(',i4,')=',i4,' nes(',i4,')=',i4)
 2001 FORMAT (' it=',i4,' rho=',es15.8,' t_test=',es15.8,' ye=',es15.8, &
& ' s=',es15.8,' s_test=',es15.8,' tmin=',es15.8,' tmax=',es15.8)
 2003 FORMAT (' Bisection iteration for e will not converge in subroutine tgvndsye_y')
 2005 FORMAT (' j=',i4,' j_ray_bndl=',i4,' k_ray_bndl=',i4,' j_angular=',i4, &
& ' k_angular=',i4,' s=',es20.11, ' s_test=',es20.11,' tmin=',es20.11, &
& ' tmax=',es20.11)

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
!  esrgn_y not called to avoid discontinuity in t
!-----------------------------------------------------------------------

  CALL eqstt_y( ivar, j, ji_ray, jk_ray, rho, t_test, ye, s_test, dsdd, &
&  dsdt, dsdy )
  IF ( DABS( s - s_test ) <= tol * DABS(s)  .and.  it /= 1 ) THEN
    t              = t_test
    RETURN
  END IF

  dt               = ( s - s_test )/( dsdt + epsilon )
  t_test           = DMIN1( DMAX1( t_test + dt, 0.9d0 * t_test ), 1.1d0 * t_test )

END DO

!-----------------------------------------------------------------------
!  Ray coordinates
!-----------------------------------------------------------------------

j_angular          = myid_y * j_ray_dim + ji_ray
k_angular          = myid_z * ik_ray_dim + jk_ray
j_ray_bndl         = MOD( myid, n_proc_y ) * ij_ray_dim + 1
k_ray_bndl         = ( myid/n_proc_y ) * ik_ray_dim + 1

!-----------------------------------------------------------------------
!  Edit diagnostics
!-----------------------------------------------------------------------

jm1                = MAX( j - 1, 1 )
jp1                = MIN( j + 1, ny )

WRITE (nprint,1001)
WRITE (nprint,1003) j, j_ray_bndl, k_ray_bndl, j_angular, k_angular, myid, &
& rho, t_prev, ye, s, s_test, t_test
WRITE (nprint,1005) jm1, nse(jm1,ji_ray,jk_ray), j, nse(j,ji_ray,jk_ray), &
& jp1, nse(jp1,ji_ray,jk_ray)
WRITE (nlog,1001)
WRITE (nlog,1003) j, j_ray_bndl, k_ray_bndl, j_angular, k_angular, myid, &
& rho, t_prev, ye, s, s_test, t_test
WRITE (nlog,1005) jm1, nse(jm1,ji_ray,jk_ray), j, nse(j,ji_ray,jk_ray), &
& jp1, nse(jp1,ji_ray,jk_ray)

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

  CALL esrgn_y_t( j, ji_ray, jk_ray, rho, t_test, ye )
  CALL eqstt_y( ivar, j, ji_ray, jk_ray, rho, t_test, ye, s_test, dsdd, &
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

j_angular          = myid_y * j_ray_dim + ji_ray
k_angular          = myid_z * ik_ray_dim + jk_ray
j_ray_bndl         = MOD( myid, n_proc_y ) * ij_ray_dim + 1
k_ray_bndl         = ( myid/n_proc_y ) * ik_ray_dim + 1

!-----------------------------------------------------------------------
!  Edit diagnostics
!-----------------------------------------------------------------------

WRITE (nprint,2003)
WRITE (nprint,2005) j, j_ray_bndl, k_ray_bndl, j_angular, k_angular, s, &
& s_test, tmin, tmax
WRITE (nlog,2003)
WRITE (nlog,2005) j, j_ray_bndl, k_ray_bndl, j_angular, k_angular, s, &
& s_test, tmin, tmax
STOP

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE tgvndsye_y
