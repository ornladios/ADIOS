SUBROUTINE flash_grid( j, ij_ray, ik_ray, t_nse )
!-----------------------------------------------------------------------
!
!    File:         flash_grid
!    Module:       flash_grid
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         8/08/04
!
!    Purpose:
!      To flash a mass zone to nse. The new temperature
!       of the mass zone is computed such that the internal
!       energy of the zone after being switched to nse is the
!       same as before. 
!
!    Subprograms called:
!  eqstta_x
!
!    Input arguments:
!  j      : radial zone index of zone to be flashed
!  ij_ray : index denoting the j-index of a specific radial ray
!  ik_ray : index denoting the k-index of a specific radial ray
!
!    Output arguments:
!  t_nse  : temperature of material in nse
!
!    Include files:
!  numerical_module, physcnst_module,
!  edit_module, eos_snc_x_module, mdl_cnfg_module
!
!-----------------------------------------------------------------------

USE kind_module
USE numerical_module, ONLY: half

USE edit_module, ONLY: nprint
USE eos_snc_x_module, ONLY: nse
USE mdl_cnfg_module, ONLY: rho, t, ye

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: j            ! radial zone index
INTEGER, INTENT(in)              :: ij_ray       ! index denoting the j-index of a specific radial ray
INTEGER, INTENT(in)              :: ik_ray       ! index denoting the k-index of a specific radial ray

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

REAL(KIND=double), INTENT(out)   :: t_nse        ! temperature of material in nse

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
REAL(KIND=double)                :: tol = 1.d-5  ! dummy variable

REAL(KIND=double)                :: u_non_nse    ! value of internal energy before flashing
REAL(KIND=double)                :: u_nse        ! value of internal energy in NSE
REAL(KIND=double)                :: u1           ! dummy variable
REAL(KIND=double)                :: u2           ! dummy variable
REAL(KIND=double)                :: u3           ! dummy variable

!-----------------------------------------------------------------------
!        Formats
!-----------------------------------------------------------------------

  101 format (' The temperature will not converge in subroutine flash_grid')
  103 format (' t(',i3,')=',1pe10.3,' t_min=',1pe10.3,' t_max=',1pe10.3,' u_non_nse=',1pe10.3,' u_nse=',1pe10.3)


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!........Iterate to determine nse t(j)..................................
!.......................................................................

!........Store t_non_nse and u_non_nse

t_non_nse           = t(j)

CALL eqstta_x( 2, j, ij_ray, ik_ray, rho(j), t(j), ye(j), u_non_nse, u1, u2, u3 )

!........Set up bisection interval

t_min               = t_minimum
t_max               = t_maximum
nse(j,ij_ray,ik_ray) = 1

DO it = 1,itmax

!........Bisect interval

  t_nse             = half * ( t_min + t_max )
  CALL eqstta_x( 2, j, ij_ray, ik_ray, rho(j), t_nse, ye(j), u_nse, u1, u2, u3 )

!........Test for convergence

  IF ( ABS( u_nse - u_non_nse ) <= tol * DABS( u_non_nse ) ) EXIT

!........Reduce bisection interval if not converged

  IF ( u_nse <= u_non_nse ) THEN
    t_min           = t_nse
  ELSE
    t_max           = t_nse
  END IF ! u_nse < u_non_nse
  
  IF ( it == itmax ) THEN
    WRITE (nprint,101)
    WRITE (nprint,103) j,t(j),t_min,t_max,u_non_nse,u_nse
  END IF

END DO

RETURN
END SUBROUTINE flash_grid