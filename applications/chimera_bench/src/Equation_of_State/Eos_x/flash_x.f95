SUBROUTINE flash_x( j, ij_ray, ik_ray, v_const )
!-----------------------------------------------------------------------
!
!    File:         flash_x
!    Module:       flash_x
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         19/06/03
!
!    Purpose:
!      To flash_x a mass zone to nse. The new temperature
!       of the mass zone is computed such that the internal
!       energy of the zone after being switched to nse is the
!       same as before. 
!
!    Subprograms called:
!  esrgn_x : regenerates the local EOS table if necessary
!  eqstt_x : interpolates quantities in the local EOS table
!
!    Input arguments:
!  j       : radial zone index of zone to be flashed
!  ij_ray  : j-index of a radial ray
!  ik_ray  : k-index of a radial ray
!  v_const : parameter (pressure or energy) that is kept constant during flash
!
!    Output arguments:
!        none
!
!    Include files:
!  kind_module, numerical_module
!  edit_module, eos_snc_x_module, mdl_cnfg_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY : double
USE numerical_module, ONLY: half

USE edit_module, ONLY: nprint
USE eos_snc_x_module, ONLY: nse
USE mdl_cnfg_module, ONLY: rho, t, ye, tr

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

CHARACTER (len=1), INTENT(in)    :: v_const       ! parameter that is kept constant during flash

INTEGER, INTENT(in)              :: j             ! radial zone index
INTEGER, INTENT(in)              :: ij_ray        ! j-index of a radial ray
INTEGER, INTENT(in)              :: ik_ray        ! k-index of a radial ray

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: it           ! iteration index
INTEGER                          :: itmax = 40   ! maximum number of iterations

REAL(KIND=double), PARAMETER     :: t_minimum = 1.d+9  ! minimum value of temperature for bisection iteration
REAL(KIND=double), PARAMETER     :: t_maximum = 1.d+11 ! maximum value of temperature for bisection iteration

REAL(KIND=double)                :: t_non_nse    ! value of temperature before flashing
REAL(KIND=double)                :: t_min        ! lower bound of temperature during bisection iteration
REAL(KIND=double)                :: t_max        ! upper bound of temperature during bisection iteration
REAL(KIND=double)                :: t_nse        ! trial temperature during bisection iteration
REAL(KIND=double)                :: tol = 1.d-8  ! dummy variable

REAL(KIND=double)                :: u_non_nse    ! value of internal energy before flashing
REAL(KIND=double)                :: u_nse        ! value of internal energy in NSE
REAL(KIND=double)                :: u1           ! dummy variable
REAL(KIND=double)                :: u2           ! dummy variable
REAL(KIND=double)                :: u3           ! dummy variable

REAL(KIND=double)                :: p_non_nse    ! value of pressure before flashing
REAL(KIND=double)                :: p_nse        ! value of pressure in NSE
REAL(KIND=double)                :: p1           ! dummy variable
REAL(KIND=double)                :: p2           ! dummy variable
REAL(KIND=double)                :: p3           ! dummy variable

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 FORMAT (' The temperature will not converge in subroutine flash_x')
  103 FORMAT (' t(',i3,')=',1pe10.3,' t_min=',1pe10.3,' t_max=',1pe10.3,' u_non_nse=',1pe10.3,' u_nse=',1pe10.3)
  105 FORMAT (' t(',i3,')=',1pe10.3,' t_min=',1pe10.3,' t_max=',1pe10.3,' p_non_nse=',1pe10.3,' p_nse=',1pe10.3)
  107 FORMAT (' v_const is neither e nor p in subroutine flash_x')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!       \\\\\ ITERATE TO DETERONE NSE T ASSUMING CONSTANT E /////
!
!-----------------------------------------------------------------------

IF ( v_const == 'e' ) THEN

!-----------------------------------------------------------------------
!  Store t_non_nse and u_non_nse
!-----------------------------------------------------------------------

  t_non_nse           = t(j)
  CALL esrgn_x( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), u_non_nse, &
&  u1, u2, u3 )

!-----------------------------------------------------------------------
!  Set up bisection interval
!-----------------------------------------------------------------------

  t_min               = t_minimum
  t_max               = t_maximum
  nse(j,ij_ray,ik_ray) = 1

  DO it = 1,itmax

!-----------------------------------------------------------------------
!  Bisect interval
!-----------------------------------------------------------------------

    t_nse             = half * ( t_min + t_max )
    CALL esrgn_x( j, ij_ray, ik_ray, rho(j), t_nse, ye(j) )
    CALL eqstt_x( 2, j, ij_ray, ik_ray, rho(j), t_nse,ye(j), u_nse, &
&    u1, u2, u3 )

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
      WRITE (nprint,103) j, t(j), t_min, t_max, u_non_nse, u_nse
    END IF

  END DO

!-----------------------------------------------------------------------
!  Store nse temperature
!-----------------------------------------------------------------------

  t(j)               = t_nse
  tr(j)              = t_nse

ELSE IF ( v_const == 'p' ) THEN

!-----------------------------------------------------------------------
!
!       \\\\\ ITERATE TO DETERONE NSE T ASSUMING CONSTANT P /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Store t_non_nse and p_non_nse
!-----------------------------------------------------------------------

  t_non_nse           = t(j)
  CALL esrgn_x( j, ij_ray, ik_ray, rho(j), t(j), ye(j) )
  CALL eqstt_x( 1, j, ij_ray, ik_ray, rho(j), t(j), ye(j), p_non_nse, &
&  p1, p2, p3 )

!-----------------------------------------------------------------------
!  Set up bisection interval
!-----------------------------------------------------------------------

  t_min               = t_minimum
  t_max               = t_maximum
  nse(j,ij_ray,ik_ray)        = 1

  DO it = 1,itmax

!-----------------------------------------------------------------------
!  Bisect interval
!-----------------------------------------------------------------------

    t_nse             = half * ( t_min + t_max )
    CALL esrgn_x( j, ij_ray, ik_ray, rho(j), t_nse, ye(j) )
    CALL eqstt_x( 1, j, ij_ray, ik_ray, rho(j), t_nse,ye(j), p_nse, &
&    p1, p2, p3 )

!-----------------------------------------------------------------------
!  Test for convergence
!-----------------------------------------------------------------------

    IF ( ABS( p_nse - p_non_nse ) <= tol * DABS( p_non_nse ) ) EXIT

!-----------------------------------------------------------------------
!  Reduce bisection interval if not converged
!-----------------------------------------------------------------------

    IF ( p_nse <= p_non_nse ) THEN
      t_min           = t_nse
    ELSE
      t_max           = t_nse
    END IF ! u_nse < u_non_nse
  
    IF ( it == itmax ) THEN
      WRITE (nprint,101)
      WRITE (nprint,103) j, t(j), t_min, t_max, p_non_nse, p_nse
    END IF

  END DO

!-----------------------------------------------------------------------
!  Store nse temperature
!-----------------------------------------------------------------------

  t(j)               = t_nse
  tr(j)              = t_nse

ELSE

  WRITE ( nprint,107 )
  STOP

END IF ! v_const == 'e'

RETURN
END SUBROUTINE flash_x