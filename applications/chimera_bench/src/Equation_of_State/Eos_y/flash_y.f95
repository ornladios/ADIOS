SUBROUTINE flash_y( j, ji_ray, jk_ray )
!-----------------------------------------------------------------------
!
!    File:         flash_y
!    Module:       flash_y
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         19/06/03
!
!    Purpose:
!      To flash a mass zone to nse. The new temperature
!       of the mass zone is computed such that the internal
!       energy of the zone after being switched to nse is the
!       same as before. 
!
!    Subprograms called:
!  esrgn_y : regenerates the local EOS table if necessary
!  eqstt_y : interpolates quantities in the local EOS table
!
!    Input arguments:
!  j       : y (angular) zone index of zone to be flashed
!  ji_ray  : x (radial) index of a specific y (angular) ray
!  jk_ray  : z (azimuthal) index of a specific y (angular) ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_y_module, mdl_cnfg_y_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: half

USE edit_module, ONLY: nprint, nlog
USE eos_snc_y_module, ONLY: nse
USE mdl_cnfg_y_module, ONLY: rho, t, ye

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j             ! y (angular) zone index of zone to be flashed
INTEGER, INTENT(in)              :: ji_ray        ! x (radial) index of a specific y (angular) ray
INTEGER, INTENT(in)              :: jk_ray        ! z (azimuthal) index of a specific y (angular) ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: it           ! iteration index
INTEGER                          :: itmax = 30   ! maximum number of iterations

REAL(KIND=double), PARAMETER     :: t_minimum = 1.d+9  ! minimum value of temperature for bisection iteration
REAL(KIND=double), PARAMETER     :: t_maximum = 1.d+11 ! maximum value of temperature for bisection iteration

REAL(KIND=double)                :: t_non_nse    ! value of temperature before flashing
REAL(KIND=double)                :: t_min        ! lower bound of temperature during bisection iteration
REAL(KIND=double)                :: t_max        ! upper bound of temperature during bisection iteration
REAL(KIND=double)                :: t_nse        ! trial temperature during bisection iteration
REAL(KIND=double)                :: tol = 1.d-5  ! dummy variable

REAL(KIND=double)                :: u_non_nse    ! value of internal energy before flashing
REAL(KIND=double)                :: u_nse        ! value of internal energy in NSE
REAL(KIND=double)                :: u1           ! dummy variable
REAL(KIND=double)                :: u2           ! dummy variable
REAL(KIND=double)                :: u3           ! dummy variable

  101 format (' The temperature will not converge in subroutine flash_y')
  103 format (' t(',i3,')=',1pe10.3,' t_min=',1pe10.3,' t_max=',1pe10.3,' u_non_nse=',1pe10.3,' u_nse=',1pe10.3)


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Iterate to determine nse t(j) by equating nes energy to non-nse
!   energy
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Store t_non_nse and u_non_nse
!-----------------------------------------------------------------------

t_non_nse           = t(j)

CALL esrgn_y( j, ji_ray, jk_ray, rho(j), t(j), ye(j) )
CALL eqstt_y( 2, j, ji_ray, jk_ray, rho(j), t(j), ye(j), u_non_nse, u1, u2, u3 )

!-----------------------------------------------------------------------
!  Set up bisection interval
!-----------------------------------------------------------------------

t_min               = t_minimum
t_max               = t_maximum
nse(j,ji_ray,jk_ray) = 1

DO it = 1,itmax

!-----------------------------------------------------------------------
!  Bisect interval
!-----------------------------------------------------------------------

  t_nse             = half * ( t_min + t_max )
  CALL esrgn_y( j, ji_ray, jk_ray, rho(j), t_nse, ye(j) )
  CALL eqstt_y( 2, j, ji_ray, jk_ray,rho(j), t_nse,ye(j), u_nse, u1, u2, u3 )

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

  IF ( ABS( u_nse - u_non_nse ) <= tol * DABS( u_non_nse ) ) EXIT

!-----------------------------------------------------------------------
!  Reduce bisection interval if not converged
!-----------------------------------------------------------------------

  IF ( u_nse <= u_non_nse ) THEN
    t_min           = t_nse
  ELSE
    t_max           = t_nse
  END IF ! u_nse < u_non_nse
  
  IF ( it == itmax ) THEN
    WRITE (nprint,101)
    WRITE (nprint,103) j, t(j), t_min, t_max, u_non_nse ,u_nse
    WRITE (nlog,101)
    WRITE (nlog,103) j, t(j), t_min, t_max, u_non_nse, u_nse
  END IF

END DO

!-----------------------------------------------------------------------
!  Store nse temperature
!-----------------------------------------------------------------------

t(j)               = t_nse

RETURN
END SUBROUTINE flash_y