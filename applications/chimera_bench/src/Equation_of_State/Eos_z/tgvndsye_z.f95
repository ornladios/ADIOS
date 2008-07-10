SUBROUTINE tgvndsye_z( k, ki_ray, kj_ray, rho, s, ye, t_prev, t )
!-----------------------------------------------------------------------
!
!    File:         tgvndsye_z
!    Module:       tgvndsye_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
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
!  esrgn_z_t : recomputes, if the temperature warrents, the EOS table entries
!  eqstt_z   : interpolates between cube corners to a given state point
!
!    Input arguments:
!  k         : z (azimuthal) zone index
!  ki_ray    : x (radial) index of a specific z (azimuthal) ray
!  kj_ray    : y (angular) index of a specific z (azimuthal) rayy
!  rho       : density (g/cm**3).
!  s         : inernal energy (ergs/g)
!  ye        : electron fraction.
!  t_prev    : temperature guess.
!
!    Output arguments:
!  t         : temperature (K)
!
!    Include files:
!  kind_module, array_module, numerical_module
!  edit_module, eos_snc_z_module, parallel_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE array_module, ONLY: n_proc_y, ij_ray_dim, ik_ray_dim, k_ray_dim, nz
USE numerical_module, ONLY: epsilon, half

USE edit_module, ONLY: nprint, nlog
USE eos_snc_z_module, ONLY : nse
USE parallel_module, ONLY : myid, myid_y, myid_z

IMPLICIT NONE
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: k             ! z (azimuthal) zone index
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray

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
INTEGER                          :: km1           ! mimimum of k-1 and 1
INTEGER                          :: kp1           ! maximum of k+1 and nz

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

 1001 FORMAT (' NR iteration for s will not converge in subroutine tgvndsye_z, &
 &will try Bisection')
  1003 FORMAT (' j=',i4,' j_ray_bndl=',i4,' k_ray_bndl=',i4,' j_angular=',i4, &
& ' k_angular=',i4,' myid=',i4,' rho=',es11.3,' t_prev=',es11.3, ' ye=',es11.3, &
& ' s=',es20.11,' s_test=',es20.11,' t_test=',es20.11)
 1005 FORMAT (' nes(',i4,')=',i4,' nes(',i4,')=',i4,' nes(',i4,')=',i4)
 2001 FORMAT (' it=',i4,' rho=',es15.8,' t_test=',es15.8,' ye=',es15.8, &
& ' s=',es15.8,' s_test=',es15.8,' tmin=',es15.8,' tmax=',es15.8)
 2003 FORMAT (' Bisection iteration for e will not converge in subroutine tgvndsye_z')
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
!  esrgn_z not called to avoid discontinuity in t
!-----------------------------------------------------------------------

  CALL eqstt_z( ivar, k, ki_ray, kj_ray, rho, t_test, ye, s_test, dsdd, dsdt, dsdy )
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

j_angular          = myid_y * ij_ray_dim + ki_ray
k_angular          = myid_z * k_ray_dim + kj_ray
j_ray_bndl         = MOD( myid, n_proc_y ) * ij_ray_dim + 1
k_ray_bndl         = ( myid/n_proc_y ) * ik_ray_dim + 1

!-----------------------------------------------------------------------
!  Edit diagnostics
!-----------------------------------------------------------------------

km1                = MAX( k - 1, 1 )
kp1                = MIN( k + 1, nz )

WRITE (nprint,1001)
WRITE (nprint,1003) k, j_ray_bndl, k_ray_bndl, j_angular, k_angular, myid, &
& rho, t_prev, ye, s, s_test, t_test
WRITE (nprint,1005) km1, nse(km1,ki_ray,kj_ray), k, nse(k,ki_ray,kj_ray), &
& kp1, nse(kp1,ki_ray,kj_ray)
WRITE (nlog,1001)
WRITE (nlog,1003) k, j_ray_bndl, k_ray_bndl, j_angular, k_angular, myid, &
& rho, t_prev, ye, s, s_test, t_test
WRITE (nlog,1005) km1, nse(km1,ki_ray,kj_ray), k, nse(k,ki_ray,kj_ray), &
& kp1, nse(kp1,ki_ray,kj_ray)

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

  CALL esrgn_z_t( k, ki_ray, kj_ray, rho, t_test, ye )
  CALL eqstt_z( ivar, k, ki_ray, kj_ray, rho, t_test, ye, s_test, dsdd, dsdt, dsdy )

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

j_angular          = myid_y * ij_ray_dim + ki_ray
k_angular          = myid_z * k_ray_dim + kj_ray
j_ray_bndl         = MOD( myid, n_proc_y ) * ij_ray_dim + 1
k_ray_bndl         = ( myid/n_proc_y ) * ik_ray_dim + 1

!-----------------------------------------------------------------------
!  Edit diagnostics
!-----------------------------------------------------------------------

WRITE (nprint,2001)
WRITE (nprint,2003) k, j_ray_bndl, k_ray_bndl, j_angular, k_angular, myid, &
& s, s_test, tmin, tmax
WRITE (nlog,2001)
WRITE (nlog,2003) k, j_ray_bndl, k_ray_bndl, j_angular, k_angular, myid, &
& s, s_test, tmin, tmax
STOP

!-----------------------------------------------------------------------
!  Done
!-----------------------------------------------------------------------

RETURN
END SUBROUTINE tgvndsye_z
