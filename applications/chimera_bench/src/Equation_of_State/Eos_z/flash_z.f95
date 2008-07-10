SUBROUTINE flash_z( k, ki_ray, kj_ray )
!-----------------------------------------------------------------------
!
!    File:         flash_z
!    Module:       flash_z
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/04/07
!
!    Purpose:
!      To flash a mass zone to nse. The new temperature
!       of the mass zone is computed such that the internal
!       energy of the zone after being switched to nse is the
!       same as before. 
!
!    Subprograms called:
!  esrgn_z : regenerates the local EOS table if necessary
!  eqstt_z : interpolates quantities in the local EOS table
!
!    Input arguments:
!  k       : z (azimuthal) zone index of zone to be flashed
!  ki_ray  : x (radial) index of a specific z (azimuthal) ray
!  kj_ray  : y (angular) index of a specific z (azimuthal) ray
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_z_module, mdl_cnfg_z_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: half

USE edit_module, ONLY: nprint, nlog
USE eos_snc_z_module, ONLY: nse
USE mdl_cnfg_z_module, ONLY: rho, t, ye

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: k             ! z (azimuthal) zone index of zone to be flashed
INTEGER, INTENT(in)              :: ki_ray        ! x (radial) index of a specific z (azimuthal) ray
INTEGER, INTENT(in)              :: kj_ray        ! y (angular) index of a specific z (azimuthal) ray

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

  101 format (' The temperature will not converge in subroutine flash_z')
  103 format (' t(',i3,')=',1pe10.3,' t_min=',1pe10.3,' t_max=',1pe10.3,' u_non_nse=',1pe10.3,' u_nse=',1pe10.3)


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Iterate to determine nse t(k) by equating nes energy to non-nse
!   energy
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Store t_non_nse and u_non_nse
!-----------------------------------------------------------------------

t_non_nse           = t(k)

CALL esrgn_z( k, ki_ray, kj_ray, rho(k), t(k), ye(k) )
CALL eqstt_z( 2, k, ki_ray, kj_ray, rho(k), t(k), ye(k), u_non_nse, u1, u2, u3 )

!-----------------------------------------------------------------------
!  Set up bisection interval
!-----------------------------------------------------------------------

t_min               = t_minimum
t_max               = t_maximum
nse(k,kj_ray,ki_ray) = 1

DO it = 1,itmax

!-----------------------------------------------------------------------
!  Bisect interval
!-----------------------------------------------------------------------

  t_nse             = half * ( t_min + t_max )
  CALL esrgn_z( k, ki_ray, kj_ray, rho(k), t_nse, ye(k) )
  CALL eqstt_z( 2, k, ki_ray, kj_ray, rho(k), t_nse,ye(k), u_nse, u1, u2, u3 )

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
    WRITE (nprint,103) k, t(k), t_min, t_max, u_non_nse ,u_nse
    WRITE (nlog,101)
    WRITE (nlog,103) k, t(k), t_min, t_max, u_non_nse, u_nse
  END IF

END DO

!-----------------------------------------------------------------------
!  Store nse temperature
!-----------------------------------------------------------------------

t(k)               = t_nse

RETURN
END SUBROUTINE flash_z